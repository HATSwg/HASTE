!-------------------------------------------------------------------------------
!   Copyright (C) 2017  Whitman T. Dailey
!   
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License version 3 as 
!   published by the Free Software Foundation.
!   
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
Module Find_Trajectory
    
    Implicit None
    Private
    Public :: Next_Event_Trajectory
    Public :: Simple_Trajectory
    Public :: Prev_Event_Trajectory
    
Contains

Subroutine Next_Event_Trajectory(sat, Gravity, r1, t1, s1cm, u_vec, outbound_search, Found, r2, tof, v1cm, v2sat, vS2)
    Use Kinds, Only: dp
    Use Satellite_Motion, Only: Satellite_Position_Type
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Converged
    Use Astro_Utilities, Only: Lambert => Lambert_Gooding
    Use Astro_Utilities, Only: Hits_Center
    Use Global, Only: n_life
    Implicit None
    Type(Satellite_Position_Type), Intent(In) :: sat
    Logical, Intent(In) :: Gravity
    Real(dp), Intent(In):: r1(1:3)       ! [km]   Location of scatter
    Real(dp), Intent(In) :: t1           ! [s]    time of scatter
    Real(dp), Intent(In) :: s1cm         ! [km/s] speed of neutron in cm system after scatter
    Real(dp), Intent(In) :: u_vec(1:3)   ! [km/s] velocity of center of mass of atom and neutron    
    Logical, Intent(In) :: outbound_search ! Whether to search for solutions on outbout or inbound trajectories
    Logical, Intent(Out) :: Found        ! Whether solution is found (neutron can reach satellite)
    Real(dp), Intent(Out) :: r2(1:3)     ! [km]   neutron/target position at time of arrival
    Real(dp), Intent(Out) :: tof         ! [s] time of flight from r1 to r2 on transfer trajectory
    Real(dp), Intent(Out) :: v1cm(1:3)   ! [km/s] v1 in CM system
    Real(dp), Intent(Out) :: v2sat(1:3)  ! [km/s]  neutron velocity at arrival in rest frame of satellite
    Real(dp), Intent(Out) :: vS2(1:3)    ! [km/s] satellite velocity at time of arrival
    Real(dp) :: tof_old          ! [s]    copy of tof used to test for convergence of tof
    Real(dp) :: v1(1:3)          ! [km/s] velocity of neutron leaving scatter (inertial frame)
    Real(dp) :: vScm(1:3)        ! [km/s] vS in CM frame
    Real(dp) :: v2(1:3)          ! [km/s] velocity of neutron at satellite rendezvous (inertial frame)
    Real(dp) :: ds_1,ds_2,ds,ds_old ![km/s] error in s1cm at points for root-solving
    Real(dp) :: tof_1,tof_2         ![s] min and max vaules for time of flight
    Real(dp) :: s1cm_1,s1cm_2       ![km/s] speed in CM frame for tof_1 and tof_2
    Real(dp) :: u_speed             ![km/s] speed of CM frame
    Real(dp) :: m   !relaxation factor for Illinois fix to false position
    Real(dp) :: x
    Integer :: i    !counter for Illinois fixes to false position
    Integer :: i_lim

    Found = .FALSE.
    If (Gravity) Then
        u_speed = Vector_Length(u_vec)
        r2 = sat%R(t1)
        If (outbound_search) Then
            !estimate minimum TOF: starting distance to target divided by max possible closing speed
            tof_1 = Vector_Length(r2 - r1) / (s1cm + u_speed + sat%vp)  !this is less than or equal to the true minimum
            !check if this minimum TOF bounds true TOF
            Call Lambert(r1,sat%R(t1+tof_1),tof_1,v1,v2)
            s1cm_1 = Vector_Length(v1 - u_vec)  !this s1cm_i needs to be higher than s1cm
            If (s1cm_1 .LT. s1cm) Return  !there is not a trajectory that joins the two points with the available neutron speed
            !estimate maximum TOF: depends on transfer orbit type...
            Call Find_max_TOF(tof_2)
            !check if this maximum TOF bounds true TOF
            Call Lambert(r1,sat%R(t1+tof_2),tof_2,v1,v2)
            s1cm_2 = Vector_Length(v1 - u_vec)  !this s1cm_i needs to be lower than s1cm
            If (s1cm_2 .GT. s1cm) Return  !there is not a trajectory that joins the two points with the available neutron speed
        Else !search for an inbound rendezvous, only valid when elliptical trajectories are possible
            !estimate minimum TOF: same criteria as the maximum TOF for the outbound case
            Call Find_max_TOF(tof_1)
            Call Lambert(r1,sat%R(t1+tof_1),tof_1,v1,v2)
            s1cm_1 = Vector_Length(v1 - u_vec)
            !Since the placement and number of roots can vary, a practical approach is to check for a sign change in the difference
            !between speed required and speed available for an intercept on a reasonable interval; if a sign change is indicated,
            !then search for that root.
            !A conservative upper limit to time of flight is arbitrarily 100 neutron lifetimes (the probability of a neutron NOT
            !decaying during a flight this long is approximately 3.72E-44... sufficently close to zero to ignore trajectories with 
            !TOFs longer than this)
            !Test at logarithmically spaced intervals until a sign change is found...
            x = Log10(tof_1) + ( Log10(100._dp*n_life) - Log10(tof_1) )
            Do i = 1,100
                tof_2 = 10._dp**( x * Real(i,dp)*0.01_dp )
                Call Lambert(r1,sat%R(t1+tof_2),tof_2,v1,v2)
                s1cm_2 = Vector_Length(v1 - u_vec)
                If ( (s1cm-s1cm_1)*(s1cm-s1cm_2) .LT. 0._dp ) Then !sign change detected, a root exists on this interval
                    Exit
                ElseIf (i .EQ. 100) Then !a sign change was not found on any interval, no root to search for
                    Return
                End If
                !advance tof_1 to check for the next interval on the next iteration
                tof_1 = tof_2
                s1cm_1 = s1cm_2
            End Do
        End If
        !min and max TOF are established and bracket a TOF for an achievable transfer trajectory
        !Use False-Position w/ Illinois modification to find TOF of transfer
        ds_1 = s1cm_1 - s1cm
        ds_2 = s1cm_2 - s1cm
        !False-Position starter
        Call New_TOF_and_DS(tof_1,tof_2,ds_1,ds_2,tof,ds)
        If (ds*ds_1 .GE. 0._dp) Then !ds and ds_1 have same sign
            tof_1 = tof
            ds_1 = ds
        Else  !ds and ds_2 have same sign
            tof_2 = tof
            ds_2 = ds
        End If
        !Proceed with False-Position w/ Illinois
        !HACK Iteration limit prevents stalling of the bracketing method at expense of missing some rendezvous
        Do i_lim = 1,1000
            tof_old = tof
            ds_old = ds
            !False-Position for next TOF
            Call New_TOF_and_DS(tof_1,tof_2,ds_1,ds_2,tof,ds)
            !compare to desired post-scatter speed and convergence of TOF interval
            !N2H Study convergence criteria for speed and tof, current criteria is probably overkill
            !w/in 1 cm/s of desired neutron speed (a 1 micro-eV neutron goes about 14 m/s..) AND w/in 1 microsecond of actual TOF
            If ( (Abs(ds).LT.1.E-6_dp .AND. Converged(tof_1,tof_2,rTol=1.E-12_dp,aTol=1.E-6_dp)) & 
               & .OR. ds.EQ.0._dp .OR. tof_1.EQ.tof_2 ) Then  !OR iteration landed on answer or closed interval
                !check that trajectory does not intersect Earth
                If (Hits_Center(r1,r2,v1,v2)) Return
                !If we get this far, a valid trajectory has been found!
                Found = .TRUE.
                v2Sat = v2 - vS2
                Return
            End If
            !check for poor behavior of false position
            If (ds*ds_old .GE. 0._dp) Then !ds and ds_old has same sign, False-Position is behaving poorly, apply Illinois
                m = 0.5_dp
                If (ds*ds_1 .GE. 0._dp) Then !function value ds_2 will be modified
                    Do i = 1,5
                        !move left bracket
                        tof_1 = tof
                        ds_1 = ds
                        !False-Position-Illinois for next TOF
                        Call New_TOF_and_DS(tof_1,tof_2,ds_1,ds_2*m,tof,ds)
                        !check for sign change in ds
                        If (ds*ds_1 .LT. 0._dp) Then  !Illinois sucessfully moved BOTH brackets w/in 5 attempts
                            !move right bracket
                            tof_2 = tof
                            ds_2 = ds
                            Exit
                        End If
                        m = m * 0.5_dp
                    End Do
                Else !function value ds_1 will be modified
                    Do i = 1,5
                        !move right bracket
                        tof_2 = tof
                        ds_2 = ds
                        !False-Position-Illinois for next TOF
                        Call New_TOF_and_DS(tof_1,tof_2,ds_1*m,ds_2,tof,ds)
                        !check for sign change in ds
                        If (ds*ds_2 .LT. 0._dp) Then  !Illinois sucessfully moved BOTH brackets w/in 5 attempts
                            !move left bracket
                            tof_1 = tof
                            ds_1 = ds
                            Exit
                        End If
                        m = m * 0.5_dp
                    End Do
                End If
            Else !False-Position is behaving well, replace appropriate bracket
                If (ds*ds_1 .GE. 0._dp) Then !ds and ds_1 have same sign
                    tof_1 = tof
                    ds_1 = ds
                Else  !ds and ds_2 have same sign
                    tof_2 = tof
                    ds_2 = ds
                End If
            End If
        End Do
    Else  !straight (no gravity) trajectory
        If (sat%is_conic) Then
            u_speed = Vector_Length(u_vec)
            r2 = sat%R(t1)
            !estimate minimum TOF: starting distance to target divided by max possible closing speed
            tof_1 = Vector_Length(r2 - r1) / (s1cm + u_speed + sat%vp)  !this is less than or equal to the true minimum
            !check if this minimum TOF bounds true TOF
            v1 = (r2 - r1) / tof_1
            s1cm_1 = Vector_Length(v1 - u_vec)  !this s1cm_i needs to be higher than s1cm
            If (s1cm_1 .LT. s1cm) Return  !there is not a trajectory that joins the two points with the available neutron speed
            !estimate maximum TOF: use same criteria when gravity is on, checked against minimum possible closing speed w/o gravity
            Call Find_max_TOF(tof_2)
            !check if straight flight time at min closing speed is a better bound
            If ( Abs(s1cm - (u_speed + sat%vp)) .GT. 0._dp ) Then
                tof_2 = Min( tof_2 , Vector_Length(r2 - r1) / Abs(s1cm - (u_speed + sat%vp)) )
            End If
            !check if this maximum TOF bounds true TOF
            v1 = (r2 - r1) / tof_2
            s1cm_2 = Vector_Length(v1 - u_vec)  !this s1cm_i needs to be lower than s1cm
            If (s1cm_2 .GT. s1cm) Return  !there is not a trajectory that joins the two points with the available neutron speed
            !min and max TOF are established and bracket a TOF for an achievable transfer trajectory
            !Use False-Position w/ Illinois modification to find TOF of transfer
            ds_1 = s1cm_1 - s1cm
            ds_2 = s1cm_2 - s1cm
            !False-Position starter
            Call New_TOF_and_DS_straight(tof_1,tof_2,ds_1,ds_2,tof,ds)
            If (ds*ds_1 .GE. 0._dp) Then !ds and ds_1 have same sign
                tof_1 = tof
                ds_1 = ds
            Else  !ds and ds_2 have same sign
                tof_2 = tof
                ds_2 = ds
            End If
            !Proceed with False-Position w/ Illinois
            !HACK Iteration limit prevents stalling of the bracketing method at expense of missing some rendezvous
            Do i_lim = 1,1000
                tof_old = tof
                ds_old = ds
                !False-Position for next TOF
                Call New_TOF_and_DS_straight(tof_1,tof_2,ds_1,ds_2,tof,ds)
                !compare to desired post-scatter speed and convergence of TOF interval
                !N2H Study convergence criteria for speed and tof, current criteria is probably overkill
                !w/in 1 cm/s of desired neutron speed (1 micro-eV neutron goes about 14 m/s..) AND w/in 1 microsecond of actual TOF
                If ( (Abs(ds).LT.1.E-6_dp .AND. Converged(tof_1,tof_2,rTol=1.E-12_dp,aTol=1.E-6_dp)) & 
                   & .OR. ds.EQ.0._dp .OR. tof_1.EQ.tof_2 ) Then  !OR iteration landed on answer or closed interval
                    Found = .TRUE.
                    vScm = vS2 - u_vec
                    v2Sat = v1cm - vScm
                    Return
                End If
                !check for poor behavior of false position
                If (ds*ds_old .GE. 0._dp) Then !ds and ds_old has same sign, False-Position is behaving poorly, apply Illinois
                    m = 0.5_dp
                    If (ds*ds_1 .GE. 0._dp) Then !function value ds_2 will be modified
                        Do i = 1,5
                            !move left bracket
                            tof_1 = tof
                            ds_1 = ds
                            !False-Position-Illinois for next TOF
                            Call New_TOF_and_DS_straight(tof_1,tof_2,ds_1,ds_2*m,tof,ds)
                            !check for sign change in ds
                            If (ds*ds_1 .LT. 0._dp) Then  !Illinois sucessfully moved BOTH brackets w/in 5 attempts
                                !move right bracket
                                tof_2 = tof
                                ds_2 = ds
                                Exit
                            End If
                            m = m * 0.5_dp
                        End Do
                    Else !function value ds_1 will be modified
                        Do i = 1,5
                            !move right bracket
                            tof_2 = tof
                            ds_2 = ds
                            !False-Position-Illinois for next TOF
                            Call New_TOF_and_DS_straight(tof_1,tof_2,ds_1*m,ds_2,tof,ds)
                            !check for sign change in ds
                            If (ds*ds_2 .LT. 0._dp) Then  !Illinois sucessfully moved BOTH brackets w/in 5 attempts
                                !move left bracket
                                tof_1 = tof
                                ds_1 = ds
                                Exit
                            End If
                            m = m * 0.5_dp
                        End Do
                    End If
                Else !False-Position is behaving well, replace appropriate bracket
                    If (ds*ds_1 .GE. 0._dp) Then !ds and ds_1 have same sign
                        tof_1 = tof
                        ds_1 = ds
                    Else  !ds and ds_2 have same sign
                        tof_2 = tof
                        ds_2 = ds
                    End If
                End If
            End Do
        Else
            !For stationary and linear detector motion, result from Simple_Trajectory is a valid trajectory
            Call sat%R_and_V(t1,r2,vS2)
            vScm = vS2 - u_vec
            Call Simple_Trajectory(r1, r2, s1cm, vScm, v1cm, tof, Found)
        End If
        If (Found) v2Sat = v1cm - vScm
    End If
Contains
    !Contains construct gives these subroutines access to parent routine variables and USE statements
    !NOTE:  These contained subroutines INTENTIONALLY CAUSE SIDE-EFFECTS in the calling routine OTHER THAN the In and Out variables
    !NOTE:  Removing these subroutines from the contained construct will change the functionality of the containing routine
    Subroutine Find_max_TOF(max_TOF)
        Use Astro_Utilities, Only: SME
        Use Astro_Utilities, Only: Parabolic_TOF
        Use Astro_Utilities, Only: Lambert_minV
        Implicit None
        Real(dp), Intent(Out) :: max_TOF
        
        If (SME(Vector_Length(r1),s1cm-u_speed) .GE. 0._dp) Then
        !neutron must be on a parabolic or hyperbolic trajectory, max TOF occurs where transfer SME is zero
            !find transfer orbit TOF where SME=0 (parabolic transfer)
            !initial guess is parabolic flight time from r1 to sat at time of scatter
            max_TOF = Parabolic_TOF(r1,r2)
            If (.NOT. sat%is_stationary) Then  !need to iterate to refine
                Do
                    tof = max_TOF
                    r2 = sat%R(t1+max_TOF)  !target position after tof_2
                    max_TOF = Parabolic_TOF(r1,r2)
                    If (Converged(max_TOF,tof,rTol=1.E-9_dp,aTol=1.E-3_dp)) Exit  !w/in 1 millisecond of TOF for transfer with SME=0
                End Do
            End If
        Else  !transfer orbits include elliptical 
        !HACK Elliptical transfers and large u_speed introduce the complication of multiple roots, currently only one TOF is found
        !UNDONE The elliptical cases in this routine also do not account for complication of the available velocity function:
        !UNDONE Multiple roots are ignored (either not found, or could potentially prevent convergence)
            !find transfer orbit TOF where minV TOF is same as Lambert TOF
            !initial guess is minimum energy flight time from r1 to sat at time of scatter
            Call Lambert_minV(r1,r2,tof = max_TOF)
            If (.NOT. sat%is_stationary) Then  !need to iterate to refine
                Do
                    tof = max_TOF
                    r2 = sat%R(t1+max_TOF)  !target position after tof_max
                    Call Lambert_minV(r1,r2,tof = max_TOF)
                    If (Converged(max_TOF,tof,rTol=1.E-9_dp,aTol=1.E-3_dp)) Then
                    !w/in 1 ms of TOF for transfer where minV TOF is same as Lambert TOF
                        Exit
                    End If
                End Do
            End If
        End If
    End Subroutine Find_max_TOF
    Subroutine New_TOF_and_DS(tof1,tof2,ds1,ds2,tof3,ds3)
        Implicit None
        Real(dp), Intent(In) :: tof1,tof2
        Real(dp), Intent(In) :: ds1,ds2
        Real(dp), Intent(Out) :: tof3,ds3

        If (ds1 .NE. ds2) Then  !False-Position for next TOF
            tof3 = (tof1*ds2 - tof2*ds1) / (ds2 - ds1)
        Else !ds1 and ds2 too close, bisect for TOF
            tof3 = 0.5_dp * (tof1 + tof2)
        End If
        !update target position
        Call sat%R_and_V(t1+tof3,r2,vS2)
        !compute Lambert velocities
        Call Lambert(r1,r2,tof3,v1,v2)
        !convert Lambert-V into neutron speed in CM frame
        v1cm = v1 - u_vec
        ds3 = Vector_Length(v1cm) - s1cm
    End Subroutine New_TOF_and_DS
    Subroutine New_TOF_and_DS_straight(tof1,tof2,ds1,ds2,tof3,ds3)
        Implicit None
        Real(dp), Intent(In) :: tof1,tof2
        Real(dp), Intent(In) :: ds1,ds2
        Real(dp), Intent(Out) :: tof3,ds3

        If (ds1 .NE. ds2) Then  !False-Position for next TOF
            tof3 = (tof1*ds2 - tof2*ds1) / (ds2 - ds1)
        Else !ds1 and ds2 too close, bisect for TOF
            tof3 = 0.5_dp * (tof1 + tof2)
        End If
        !update target position
        Call sat%R_and_V(t1+tof3,r2,vS2)
        !compute velocity to intercept
        v1 = (r2 - r1) / tof3
        !convert intercept-V into neutron speed in CM frame
        v1cm = v1 - u_vec
        ds3 = Vector_Length(v1cm) - s1cm
    End Subroutine New_TOF_and_DS_straight
    !NOTE:  The preceeding contained subroutines INTENTIONALLY CAUSE SIDE-EFFECTS in the calling routine
    !NOTE:  Removing these subroutines from the contained construct will change the functionality of the containing routine
End Subroutine Next_Event_Trajectory
    
Subroutine Simple_Trajectory(r1, r2, s1cm, vScm, v1cm, dt, Found)
    Use Kinds, Only: dp
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Cross_Product
    Implicit None
    Real(dp), Intent(In) :: r1(1:3)      ! [km]   Location of scatter
    Real(dp), Intent(In) :: r2(1:3)      ! [km]   satellite position at time of scatter (actual or adjusted)
    Real(dp), Intent(In) :: s1cm         ! [km/s] speed of neutron in cm system after scatter
    Real(dp), Intent(In) :: vScm(1:3)    ! [km/s] velocity of satellite in cm frame
    Real(dp), Intent(Out) :: v1cm(1:3)   ! [km/s] velocity of neutron leaving scatter in cm frame
    Real(dp), Intent(Out) :: dt          ! [s]    time of flight
    Logical, Intent(Out) :: Found
    ! A, B, C are basis vectors for Cartesian coordinates in CM system
    Real(dp) :: A(1:3)           ! points from scatter point to satellite position at t1
    Real(dp) :: B(1:3)           ! basis vector orthogonal to A and vScm
    Real(dp) :: C(1:3)           ! basis vector in plane of A and vScm
    Real(dp) :: s_Closing        ! [km/s] speed of neutron at arrival in satellite frame
    Real(dp) :: vScm_A, vScm_C   ! components of vScm (No B component by construction of A, B, C)
    Real(dp) :: v1cm_A, v1cm_C   ! components of v1cm (No B component by construction of A, B, C)

    Found = .FALSE.
    A = Unit_Vector(r2-r1)
    If (Any(vScm .NE. 0._dp)) Then
        ! Criteria for neutron to meet satellite are: v1cm_A > vScm and v1cm_A = vScm_A
        B = Unit_Vector(Cross_Product(vScm, A))
        C = Cross_Product(A, B)
        vScm_A = Dot_Product(vScm, A)
        vScm_C = Dot_Product(vScm, C)
        v1cm_C = vScm_C     ! Criterion for neutron to meet satellite
        If (s1cm .LT. v1cm_C) Return
        v1cm_A = Sqrt(s1cm**2 - v1cm_C**2)  ! No B component by construction of A, B, and C
        v1cm = v1cm_A * A + v1cm_C * C      ! needed for computing scatter cosines mu0cm and mu0
        s_Closing = v1cm_A - vScm_A
        If (s_Closing .LE. 0._dp) Return
        dt = Vector_Length(r2-r1) / s_Closing
        Found = .TRUE.
    Else !special case for apparently stationary target
        v1cm = s1cm * A
        dt = Vector_Length(r2-r1) / s1cm
        Found = .TRUE.
    End If
End Subroutine Simple_Trajectory

Subroutine Prev_Event_Trajectory(sat, Gravity, t2, v2sat, Found, r1, v1, tof)
    !Given an arrival velocity (speed & direction) to a detector in the 
    !detector frame of reference, this routine determines the properties of a 
    !neutron originating in the source region to produce the detection event
    !UNDONE This routine finds tof to the surface of the central body (adequate for initial intent of this procedure)
    !UNDONE This routine does not find the loaction of emission in a scattering medium (it puts the emission on a sphere with r=Rc)
    Use Kinds, Only: dp
    Use Satellite_Motion, Only: Satellite_Position_Type
    Use Astro_Utilities, Only: Radius_of_Apoapsis
    Use Astro_Utilities, Only: Radius_of_Periapsis
    Use Astro_Utilities, Only: Time_Since_Periapsis
    Use Astro_Utilities, Only: Period
    Use Astro_Utilities, Only: SME
    Use Astro_Utilities, Only: Time_to_R
    Use Astro_Utilities, Only: Kepler_Gooding
    Use Global, Only: Rc => R_center
    Use Global, Only: SOIc => SOI_center
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Unit_Vector
    Implicit None
    Type(Satellite_Position_Type), Intent(In) :: sat
    Logical, Intent(In) :: Gravity
    Real(dp), Intent(In) :: t2  ![s]  time of intercept
    Real(dp), Intent(In) :: v2sat(1:3)  ! [km/s]  neutron velocity at arrival in rest frame of satellite
    Logical, Intent(Out) :: Found        ! Whether solution is found (neutron can reach satellite)
    Real(dp), Intent(Out):: r1(1:3)       ! [km]   Location of emission resulting in intercept
    Real(dp), Intent(Out):: v1(1:3)       ! [km]   Emission velocity (in intertial frame) resulting in intercept from r1
    Real(dp), Intent(Out) :: tof          ! [s]    time of flight from emission at r1 with velocity v1 (time BEFORE t2)
    Real(dp) :: r2(1:3)          ! [km]   location of intercept
    Real(dp) :: v2(1:3)          ! [km/s] velocity of neutron at satellite (inertial frame)
    Real(dp) :: r2mag, zeta, r_ca
    Real(dp) :: t_min,t_max,half_P

    Found = .FALSE.
    !Get detector position and velocity at time of intercept
    Call sat%R_and_V(t2,r2,v2)
    !Shift detected neutron into intertial frame
    v2 = v2sat + v2
    If (Gravity) Then
        If (Dot_Product(r2,v2) .LE. 0._dp) Then !the neutron is on the INBOUND or RETURN portion of the trajectory
            If (SME(r2,v2).GE.0._dp) Then
                !The neutron required for this intercept comes from the INBOUND portion of an escape trajectory...
                !The intercept is not possible since the neutron must originate from the surface BELOW
                Return
            Else  !this neutron is on the return portion of a closed trajectory
                If (Radius_of_Apoapsis(r2,v2) .GT. SOIc) Then
                    !this trajectory departs the gravitational SOI of the attracting body
                    !the two-body approximation is no good, and the tof is stupid long
                    Return
                End If
                !the search needs to be along the 'long' way around the orbit
                half_P = 0.5_dp*Period(r2,v2)
                t_min = half_P + Time_Since_Periapsis(r2,v2)
                t_max = t_min + half_P
            End If
        Else  !the neutron in on the OUTBOUND or DEPARTING portion of the trajectory
            t_min = 0._dp
            t_max = Time_Since_Periapsis(r2,v2)
        End If
        If (Radius_of_Periapsis(r2,v2) .LE. Rc) Then  !there is a possible trajectory from the surface at Rc
            found = .TRUE.
            !Reflect the velocity vector to change orbit direction so we're searching "forward"
            v2 = -v2
            !intercept ocurrs between t2 (by definition above the surface) and time of periapsis (by construction below the surface)
            tof = Time_to_R(r2,v2,Rc,t_min,t_max)
            Call Kepler_Gooding(r2,v2,tof,r1,v1)
            v1 = -v1  !emission velocity for intercept is in the opposite direction on this orbit
        Else  !trajectory does not intersect the central body or the scattering medium
            Return
        End If
    Else  !straight (no gravity) trajectory
    !UNDONE Straight (no gravity) trajectory calculations not necessary for initial implementation

        r2mag = Vector_Length(r2)
        zeta = Dot_Product(r2/r2mag,Unit_Vector(v2))
        If (zeta .LT. 0._dp) Then !only downward trajectories can hit the surface
            r_ca = r2mag * Sqrt(1._dp - zeta**2)
            If (r_ca .LT. Rc) Then !the trajectory intersects the surface
                !Found = .TRUE.
                   !UNDONE
                   !UNDONE
                   !UNDONE
            End If
        End If
    End If
End Subroutine Prev_Event_Trajectory

End Module Find_Trajectory
