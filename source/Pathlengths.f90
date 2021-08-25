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
Module Pathlengths

    Use Kinds, Only: dp
    Implicit None
    Private
    Public :: S_and_L_to_edge
    Public :: Next_Event_L_to_edge
    Public :: L_to_S
    Public :: EPL
    Public :: EPL_upward
    Public :: R_close_approach
    Public :: deltaZ
#   if TEPL
     Public :: zeta_downward
     Public :: S_upward
#   endif

    Interface Next_Event_L_to_edge
        Module Procedure NE_L_to_edge_straight
        Module Procedure NE_L_to_edge_orbit
    End Interface Next_Event_L_to_edge

    Interface L_to_S
        Module Procedure L_to_S_exact
        !Module Procedure L_to_S_fast
    End Interface L_to_S

    Interface zeta_downward
        Module Procedure zeta_downward_r1r_ca
        Module Procedure zeta_downward_r0r1zeta0
    End Interface zeta_downward

    Interface EPL_upward
        Module Procedure EPL_straight_upward_layers
        Module Procedure EPL_straight_upward_full
        Module Procedure EPL_orbit_upward_layers
        Module Procedure EPL_orbit_upward_full
    End Interface EPL_upward

Contains

Subroutine S_and_L_to_edge(atm,r0,z0,zeta0,s,L,nb,bb,Lb,db)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Global, Only: Rc => R_center
    Implicit None
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: r0,z0
    Real(dp), Intent(In) :: zeta0
    Real(dp), Intent(Out) :: s
    Real(dp), Intent(Out) :: L
    Integer, Intent(Out) :: nb
    Integer, Intent(Out), Allocatable :: bb(:)
    Real(dp), Intent(Out), Allocatable :: Lb(:)
    Logical, Intent(Out), Allocatable :: db(:)
    Real(dp) :: r_ca
    Real(dp) :: zeta1
    Integer :: nb1,nb2
    Integer, Allocatable :: bb1(:),bb2(:)
    Real(dp), Allocatable :: Lb1(:),Lb2(:)
    Integer :: i,j
    Logical, Parameter :: highfi = .FALSE.

    If (zeta0 .LT. 0._dp) Then  !direction is downward
        r_ca = R_close_approach(r0,zeta0)
        If (r_ca .LT. atm%R_bot) Then  !downward to earth
            !convert to an upward ray and compute
            zeta1 = -zeta_downward(atm%R_bot,r_ca)
            s = S_upward( atm%R_bot,r0-atm%R_bot,zeta1 )
            Call EPL_upward(atm,atm%R_bot,zeta1,atm%Z_bot,nb,bb1,Lb1,z0)
            Allocate(bb(1:nb))
            Allocate(Lb(1:nb))
            !reorder back to downward ray
            j = 1
            Do i = nb,1,-1
                bb(j) = bb1(i)
                Lb(j) = Lb1(i)
                j = j + 1
            End Do
            Allocate(db(1:nb))
            db = .TRUE.
        Else  !downward to point of closest approach and then upward to top of atmosphere
            !add downward and upward segments
            s = -zeta0 * r0 + S_upward(r_ca,atm%R_top-r_ca,0._dp)
            Call EPL_upward(atm,r_ca,0._dp,r_ca-Rc,nb1,bb1,Lb1,z0)
            Call EPL_upward(atm,r_ca,0._dp,r_ca-Rc,nb2,bb2,Lb2)
            nb = nb1 + nb2
            Allocate(bb(1:nb))
            Allocate(Lb(1:nb))
            !reorder back to downward ray
            j = 1
            Do i = nb1,1,-1
                bb(j) = bb1(i)
                Lb(j) = Lb1(i)
                j = j + 1
            End Do
            Lb(nb1+1:nb) = Lb2
            bb(nb1+1:nb) = bb2
            Allocate(db(1:nb))
            db(1:nb1) = .TRUE.
            db(nb1+1:nb) = .FALSE.
        End If
    Else  !direction is upward to top of atmosphere
        s = S_upward(r0,atm%R_top-r0,zeta0)
        Call EPL_upward(atm,r0,zeta0,z0,nb,bb,Lb)
        Allocate(db(1:nb))
        db = .FALSE.
    End If
    L = Sum(Lb)
End Subroutine S_and_L_to_edge

Function EPL(atm,r0,z0,zeta0,s) Result(L)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Global, Only: Rc => R_center
    Implicit None
    Real(dp) :: L
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: r0
    Real(dp), Intent(In) :: z0
    Real(dp), Intent(In) :: zeta0
    Real(dp), Intent(In) :: s
    Real(dp) :: L1,L2
    Real(dp) :: s_ca,r_ca,z_ca
    Real(dp) :: z1,r1,zeta1

    If (zeta0 .LT. 0._dp) Then  !direction is downward
        s_ca = -zeta0 * r0
        r_ca = R_close_approach(r0,zeta0)
        If (s .GT. s_ca) Then !downward to point of closest approach and upward to z1
            z_ca = r_ca - Rc
            !add downward and upward segments
            Call EPL_upward(atm,r_ca,0._dp,z_ca,L1,z0)
            z1 = z_ca + deltaZ(r_ca,0._dp,s-s_ca)
            Call EPL_upward(atm,r_ca,0._dp,z_ca,L2,z1)
            L = L1 + L2
        Else  !just a downward segment
            z1 = z0 + deltaZ(r0,zeta0,s)
            r1 = Rc + z1
            zeta1 = -zeta_downward(r1,r_ca)
            Call EPL_upward(atm,r1,zeta1,z1,L,z0)
        End If
    Else  !direction is upward
        z1 = z0 + deltaZ(r0,zeta0,s)
        Call EPL_upward(atm,r0,zeta0,z0,L,z1)
    End If
End Function EPL

Function NE_L_to_edge_straight(atm,r0,z0,zeta0,L) Result(path_to_ground)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Global, Only: Rc => R_center
    Implicit None
    Logical :: path_to_ground
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: r0,z0
    Real(dp), Intent(In) :: zeta0
    Real(dp), Intent(Out) :: L
    Real(dp) :: r_ca
    Real(dp) :: L1,L2
    Logical, Parameter :: highfi = .FALSE.

    path_to_ground = .FALSE.
    If (zeta0 .LT. 0._dp) Then  !direction is downward
        r_ca = R_close_approach(r0,zeta0)
        If (r_ca .LT. atm%R_bot) Then  !downward to earth
            path_to_ground = .TRUE.  !set flag and return
            L = Huge(L)
        Else  !downward to point of closest approach and then upward to top of atmosphere
            !add downward and upward segments
            Call EPL_straight_upward_full(atm,r_ca,0._dp,r_ca-Rc,L1,z0)
            Call EPL_straight_upward_full(atm,r_ca,0._dp,r_ca-Rc,L2)
            L = L1 + L2
        End If
    Else  !direction is upward to top of atmosphere
        Call EPL_straight_upward_full(atm,r0,zeta0,z0,L)
    End If
End Function NE_L_to_edge_straight

Function NE_L_to_edge_orbit(atm,big_r,Sn,z0,zeta0,L) Result(path_to_ground)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Global, Only: Rc => R_center
    Use Global, Only: mu => grav_param
    Implicit None
    Logical :: path_to_ground
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: big_r
    Real(dp), Intent(In) :: Sn
    Real(dp), Intent(In) :: z0
    Real(dp), Intent(In) :: zeta0
    Real(dp), Intent(Out) :: L
    Real(dp) :: h,xi
    Real(dp) :: p,e,rp
    Real(dp) :: L1,L2
    Logical, Parameter :: highfi = .FALSE.

    path_to_ground = .FALSE.  !for curved path implementation, LOS has already been checked, return default of FALSE
    h = big_r * Sn * Sqrt(1._dp - zeta0**2)
    xi = 0.5_dp * Sn**2 - mu / big_r
    p = h*h / mu
    e = Sqrt(1._dp + 2._dp * xi * p / mu)
    rp = p / (1._dp + e)
    If (zeta0 .LT. 0._dp) Then  !direction is downward
        !downward to pariapsis and then upward to top of atmosphere
        !add downward and upward segments
        Call EPL_orbit_upward_full(atm,rp-Rc,0._dp,h,xi,p,e,rp,L1,z0)
        Call EPL_orbit_upward_full(atm,rp-Rc,0._dp,h,xi,p,e,rp,L2)
        L = L1 + L2
    Else  !direction is upward to top of atmosphere
        Call EPL_orbit_upward_full(atm,z0,zeta0,h,xi,p,e,rp,L)
    End If
End Function NE_L_to_edge_orbit

Function R_close_approach(big_r,zeta) Result(r_ca)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: r_ca
    Real(dp), Intent(In) :: big_r,zeta

    r_ca = big_r * Sqrt(1._dp - zeta**2)
End Function R_close_approach

Function S_upward(big_R,dZ,zeta) Result(s)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: s
    Real(dp), Intent(In) :: big_R  !R_center plus base altitude
    Real(dp), Intent(In) :: dZ  !change in altitude from big_R after distance s, must be positive
    Real(dp), Intent(In) :: zeta  !zenith cosine at R

    s = dZ * (2._dp * big_R + dZ) / ( zeta * big_R + Sqrt( (zeta * big_R)**2 + dZ * (2._dp * big_R + dZ) ) )
End Function S_upward

Function zeta_upward(r0,dZ,zeta0) Result(zeta1)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: zeta1
    Real(dp), Intent(In) :: r0  !R_center plus base altitude
    Real(dp), Intent(In) :: dZ  !change in altitude from big_R after distance s, must be positive
    Real(dp), Intent(In) :: zeta0  !zenith cosine at r0

    zeta1 = Sqrt( (r0*zeta0)**2 + 2._dp*r0*dZ + dZ**2 ) / (r0 + dZ)
End Function zeta_upward

Function zeta_orbit_upward(rp,vp,p,e,r1) Result(zeta)
    Use Kinds, Only: dp
    Use Global, Only: mu => grav_param
    Implicit None
    Real(dp) :: zeta
    Real(dp), Intent(In) :: rp,vp,p,e
    Real(dp), Intent(In) :: r1
    Real(dp) :: v1

    v1 =  Sqrt(mu * (2._dp / r1 - (1._dp - e**2) / p))
    zeta = 1._dp - (rp*vp / (r1*v1))**2
End Function zeta_orbit_upward

Function zeta_downward_r1r_ca(r1,r_ca) Result(zeta1)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: zeta1
    Real(dp), Intent(In) :: r1  !radius at which new zeta is wanted
    Real(dp), Intent(In) :: r_ca  !radius of closest approach along the path
    Real(dp) :: x

    !TODO Numerical conditioning can be improved
    x = 1._dp - (r_ca / r1)**2
    If (Abs(x) .LT. 10._dp*Spacing(1._dp)) Then
        zeta1 = 0._dp
    Else
        zeta1 = -Sqrt(x)
    End If
End Function zeta_downward_r1r_ca

Function deltaZ(big_r,zeta,s)
    Use Kinds, Only: dp
    Use Utilities, Only: Smaller_Quadratic_Root
    Implicit None
    Real(dp) :: deltaZ
    Real(dp), Intent(In) :: big_r
    Real(dp), Intent(In) :: zeta
    Real(dp), Intent(In) :: s
        
    deltaZ = Smaller_Quadratic_root(big_r,s*(2._dp*big_r*zeta + s))
End Function

Function zeta_downward_r0r1zeta0(r0,r1,zeta0) Result(zeta1)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: zeta1
    Real(dp), Intent(In) :: r0  !R_center plus base altitude
    Real(dp), Intent(In) :: r1  !radius at which new zeta is wanted
    Real(dp), Intent(In) :: zeta0  !zenith cosine at r0

    zeta1 = zeta_downward_r1r_ca(r1,R_close_approach(r0,zeta0))
End Function zeta_downward_r0r1zeta0

Subroutine EPL_straight_upward_full(atm,r0,zeta0,z0,L,z1_in)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Implicit None
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: r0
    Real(dp), Intent(In) :: zeta0
    Real(dp), Intent(In) :: z0
    Real(dp), Intent(Out) :: L
    Real(dp), Intent(In), Optional :: z1_in
    Integer :: nb
    Integer, Allocatable :: bb(:)
    Real(dp), Allocatable :: Lb(:)

    If (Present(z1_in)) Then
        Call EPL_straight_upward_layers(atm,r0,zeta0,z0,nb,bb,Lb,z1_in)
    Else
        Call EPL_straight_upward_layers(atm,r0,zeta0,z0,nb,bb,Lb)
    End If
    L = Sum(Lb)
End Subroutine EPL_straight_upward_full

Subroutine EPL_orbit_upward_full(atm,z0,zeta0,h,xi,p,e,rp,L,z1_in)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Implicit None
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: z0,zeta0
    Real(dp), Intent(In) :: h,xi,p,e,rp
    Real(dp), Intent(Out) :: L
    Real(dp), Intent(In), Optional :: z1_in
    Integer :: nb
    Integer, Allocatable :: bb(:)
    Real(dp), Allocatable :: Lb(:)

    If (Present(z1_in)) Then
        Call EPL_orbit_upward_layers(atm,z0,zeta0,h,xi,p,e,rp,nb,bb,Lb,z1_in)
    Else
        Call EPL_orbit_upward_layers(atm,z0,zeta0,h,xi,p,e,rp,nb,bb,Lb)
    End If
    L = Sum(Lb)
End Subroutine EPL_orbit_upward_full

# if CHECK_L
Subroutine Check_EPL(L,z0,z1,zeta0,b,atm)
    !UNSTANDARD: ABORT (GFORT) is an extension
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Global, Only: Rc => R_center
    Use Quadratures, Only: Romberg_Quad
    Use Utilities, Only: Prec
    Implicit None
    Real(dp), Intent(In) :: L
    Real(dp), Intent(In) :: z0,z1,zeta0
    Integer, Intent(In) :: b
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp) :: r0,dZ,Smax
    Real(dp) :: L0
    Real(dp) :: Pgoal
    Real(dp) :: zeta00
    Real(dp), Parameter :: two_thirds = 2._dp / 3._dp

    r0 = Rc + z0
    dZ = z1 - z0
    Smax = dZ * (2._dp * r0 + dZ) / ( zeta0 * r0 + Sqrt( (zeta0 * r0)**2 + dZ * (2._dp * r0 + dZ) ) )
    L0 = Romberg_Quad(EPL_Integrand,0._dp,Smax,aTol=0._dp,rTol=1.E-9_dp)
    Pgoal = two_thirds * Real(atm%EPL_Prec,dp)
    If (Prec(L0,L) .GT. Pgoal) Return !computed EPLs agree to at least two thirds of the desired precision
    Write(*,*)
    Write(*,'(A,ES24.16)')           'ERROR:  PATHLENGTHS failed integral check:  L = ',L
    Write(*,'(A,ES24.16)')           '                                           L0 = ',L0
    Write(*,'(A,ES24.16)')           '                                            p = ',Prec(L0,L)
    Write(*,'(A,ES24.16,A,ES24.16)') '        z0 = ',z0,         '  z1 = ',z1
    Write(*,'(A,ES24.16,A,ES24.16)') '       Zbs = ',atm%zb(b-1),'  ---  ',atm%zb(b)
    Write(*,'(A,ES24.16)')           '     zeta0 = ',zeta0
    dZ = atm%Zb(b) - atm%zb(b-1)
    r0 = Rc + atm%zb(b-1)
    If (z0 .NE. atm%zb(b-1)) Then
        zeta00 = zeta_upward(Rc+z0,Rc+z0-r0,zeta0)
    Else
        zeta00 = zeta0
    End If
    Write(*,'(A,ES24.16)')           '    zeta00 = ',zeta00
    Write(*,'(A,ES24.16,A,F6.2,A)')  '         S = ', &
                                  &  Smax, &
                                  &  ' (', &
                                  &  100._dp*( Smax /  &
                                  &          & (dZ*(2._dp * r0 + dZ) / (zeta00*r0 + Sqrt((zeta00*r0)**2 + dZ*(2._dp*r0 + dZ)))) & 
                                             & ), &
                                  &  '% of full layer path)'
#   if GFORT
        Call abort  !<--GFORT implementation
#   else
        ERROR STOP
#   endif   
Contains
    Function EPL_Integrand(s) Result(f)
        Use Kinds, Only: dp
        Use Utilities, Only: Smaller_Quadratic_Root
        Use Atmospheres, Only: inv_rho_SL
        Implicit None
        Real(dp) :: f
        Real(dp), Intent(In) :: s
        Real(dp) :: deltaZ

        deltaZ = Smaller_Quadratic_root(r0,s*(2._dp*r0*zeta0 + s))
        f = inv_rho_SL * atm%rho(z0 + deltaZ)
    End Function EPL_Integrand
End Subroutine
# endif

Subroutine EPL_straight_upward_layers(atm,r0,zeta0,z0,nb,bb,Lb,z1_in)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Utilities, Only: Bisection_Search
    Implicit None
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: r0
    Real(dp), Intent(In) :: zeta0
    Real(dp), Intent(In) :: z0
    Integer, Intent(Out) :: nb
    Integer, Intent(Out), Allocatable :: bb(:)
    Real(dp), Intent(Out), Allocatable :: Lb(:)
    Real(dp), Intent(In), Optional :: z1_in
    Integer :: b0,b1,b
    Real(dp) :: zeta1
    Logical :: partial_last_layer
    Real(dp) :: z1
    Integer :: i

    b0 = (atm%iZb(1)-1) + Bisection_Search(z0,atm%Zb(atm%iZb(1):atm%iZb(2)-1),atm%iZb(3)-1)
    If (Present(z1_in)) Then  !not going all the way to top of atm
        z1 = z1_in
        b1 = (atm%iZb(1)-1) + Bisection_Search(z1,atm%Zb(atm%iZb(1):atm%iZb(2)-1),atm%iZb(3)-1)
        partial_last_layer = .TRUE.
    Else
        z1 = atm%Z_top
        b1 = atm%iZb(2)
        partial_last_layer = .FALSE.
    End If
    nb = b1 - b0 + 1
    Allocate(bb(1:nb))
    bb = -1
    Allocate(Lb(1:nb))
    Lb = 0._dp
    If (b0 .EQ. b1) Then !single partial layer
        If (zeta0 .GE. 0.1_dp) Then
            Lb = EPL_Z_partial_layer(r0,z0,z1,zeta0,b0,atm%EPL_lay(b0)%nZ,atm)
        Else
            Lb = EPL_S_partial_layer(r0,z0,z1,zeta0,b0,atm%EPL_lay(b0)%nS,atm)
        End If
#       if CHECK_L
            Call Check_EPL(Lb(1),z0,z1,zeta0,b0,atm)
#       endif
        bb(1) = b0
        Return
    Else
        !first layer may or may not be partial
        bb(1) = b0
        If (z0 .NE. atm%Zb(b0-1)) Then  !partial first layer
            If (zeta0 .GE. 0.1_dp) Then
                Lb(1) = EPL_Z_partial_layer(r0,z0,atm%Zb(b0),zeta0,b0,atm%EPL_lay(b0)%nZ,atm)
            Else
                Lb(1) = EPL_S_partial_layer(r0,z0,atm%Zb(b0),zeta0,b0,atm%EPL_lay(b0)%nS,atm)
            End If
#           if CHECK_L
                Call Check_EPL(Lb(1),z0,atm%Zb(b0),zeta0,b0,atm)
#           endif
        Else  !full first layer
            If (zeta0 .GE. 0.1_dp) Then
                Lb(1) = EPL_Z_known_layer(atm%EPL_lay(b0),zeta0)
            Else
                Lb(1) = EPL_S_known_layer(zeta0,b0,atm)
            End If
#           if CHECK_L
                Call Check_EPL(Lb(1),atm%Zb(b0-1),atm%Zb(b0),zeta0,b0,atm)
#           endif
        End If
        !full layers from b0+1 to b1-1
        i = 2
        Do b = b0+1,b1-1
            bb(i) = b
            zeta1 = zeta_upward(r0,atm%Zb(b-1)-z0,zeta0)
            If (zeta1 .GE. 0.1_dp) Then
                Lb(i) = EPL_Z_known_layer(atm%EPL_lay(b),zeta1)
            Else
                Lb(i) = EPL_S_known_layer(zeta1,b,atm)
            End If
#           if CHECK_L
                Call Check_EPL(Lb(i),atm%Zb(b-1),atm%Zb(b),zeta1,b,atm)
#           endif
            i = i + 1
        End Do
        !last layer may or may not be partial
        bb(nb) = b1
        zeta1 = zeta_upward(r0,atm%Zb(b1-1)-z0,zeta0)
        If (partial_last_layer) Then
            If (zeta1 .GE. 0.1_dp) Then
                Lb(nb) = EPL_Z_partial_layer(atm%Rb(b1-1),atm%Zb(b1-1),z1,zeta1,b1,atm%EPL_lay(b1)%nZ,atm)
            Else
                Lb(nb) = EPL_S_partial_layer(atm%Rb(b1-1),atm%Zb(b1-1),z1,zeta1,b1,atm%EPL_lay(b1)%nS,atm)
            End If
#           if CHECK_L
                Call Check_EPL(Lb(nb),atm%Zb(b1-1),z1,zeta1,b1,atm)
#           endif
        Else
            If (zeta1 .GE. 0.1_dp) Then
                Lb(nb) = EPL_Z_known_layer(atm%EPL_lay(b),zeta1)
            Else
                Lb(nb) = EPL_S_known_layer(zeta1,b,atm)
            End If
#           if CHECK_L
                Call Check_EPL(Lb(nb),atm%Zb(b1-1),atm%Zb(b1),zeta1,b1,atm)
#           endif
        End If
    End If
End Subroutine EPL_straight_upward_layers

Subroutine EPL_orbit_upward_layers(atm,z0,zeta0,h,xi,p,e,rp,nb,bb,Lb,z1_in)
    Use Kinds, Only: dp
    Use Global, Only: mu => grav_param
    !Use Global, Only: Rc => R_center
    Use Atmospheres, Only: Atmosphere_Type
    Use Utilities, Only: Bisection_Search
    Implicit None
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: z0
    Real(dp), Intent(In) :: zeta0
    Real(dp), Intent(In) :: h,xi,p,e,rp
    Integer, Intent(Out) :: nb
    Integer, Intent(Out), Allocatable :: bb(:)
    Real(dp), Intent(Out), Allocatable :: Lb(:)
    Real(dp), Intent(In), Optional :: z1_in
    Integer :: b0,b1,b
    Real(dp) :: zeta1
    Logical :: partial_last_layer
    Real(dp) :: z1
    Integer :: i
    !Real(dp) :: k0
    Real(dp) :: vp

    b0 = (atm%iZb(1)-1) + Bisection_Search(z0,atm%Zb(atm%iZb(1):atm%iZb(2)-1),atm%iZb(3)-1)
    If (Present(z1_in)) Then  !not going all the way to top of atm
        z1 = z1_in
        b1 = (atm%iZb(1)-1) + Bisection_Search(z1,atm%Zb(atm%iZb(1):atm%iZb(2)-1),atm%iZb(3)-1)
        partial_last_layer = .TRUE.
    Else
        z1 = atm%Z_top
        b1 = atm%iZb(2)
        partial_last_layer = .FALSE.
    End If
    nb = b1 - b0 + 1
    Allocate(bb(1:nb))
    bb = -1
    Allocate(Lb(1:nb))
    Lb = 0._dp
    If (b0 .EQ. b1) Then !single partial layer
        If (zeta0 .GE. 0.1_dp) Then
            Lb = EPL_Rk_partial_layer(z0,z1,p,e,rp,b0,atm%EPL_lay(b0)%nRk,atm)
        Else
            Lb = EPL_Tk_partial_layer(z0,z1,p,e,b0,atm%EPL_lay(b0)%nTk,atm)
        End If
        bb(1) = b0
        Return
    Else
        !first layer may or may not be partial
        bb(1) = b0
        If (z0 .NE. atm%Zb(b0-1)) Then  !partial first layer
            If (zeta0 .GE. 0.1_dp) Then
                Lb(1) = EPL_Rk_partial_layer(z0,atm%Zb(b0),p,e,rp,b0,atm%EPL_lay(b0)%nRk,atm)
            Else
                Lb(1) = EPL_Tk_partial_layer(z0,atm%Zb(b0),p,e,b0,atm%EPL_lay(b0)%nTk,atm)
            End If
        Else  !full first layer
            If (zeta0 .GE. 0.1_dp) Then
                Lb(1) = EPL_Rk_known_layer(atm%EPL_lay(b0),p,e,rp)
            Else
                Lb(1) = EPL_Tk_known_layer(p,e,b0,atm)
            End If
        End If
        !full layers from b0+1 to b1-1
        i = 2
        !k0 = xi + mu / (Rc + z0)
        vp = Sqrt(mu * (2._dp / rp - (1._dp - e**2) / p))
        Do b = b0+1,b1-1
            bb(i) = b
            zeta1 = zeta_orbit_upward(rp,vp,p,e,atm%Rb(b-1))
            If (zeta1 .GE. 0.1_dp) Then
                Lb(i) = EPL_Rk_known_layer(atm%EPL_lay(b),p,e,rp)
            Else
                Lb(i) = EPL_Tk_known_layer(p,e,b,atm)
            End If
            i = i + 1
        End Do
        !last layer may or may not be partial
        bb(nb) = b1
        zeta1 = zeta_orbit_upward(rp,vp,p,e,atm%Rb(b1-1))
        If (partial_last_layer) Then
            If (zeta1 .GE. 0.1_dp) Then
                Lb(nb) = EPL_Rk_partial_layer(atm%Zb(b1-1),z1,p,e,rp,b1,atm%EPL_lay(b1)%nRk,atm)
            Else
                Lb(nb) = EPL_Tk_partial_layer(atm%Zb(b1-1),z1,p,e,b1,atm%EPL_lay(b1)%nTk,atm)
            End If
        Else
            If (zeta1 .GE. 0.1_dp) Then
                Lb(nb) = EPL_Rk_known_layer(atm%EPL_lay(b1),p,e,rp)
            Else
                Lb(nb) = EPL_Tk_known_layer(p,e,b1,atm)
            End If
        End If
    End If
End Subroutine EPL_orbit_upward_layers

Function EPL_S_partial_layer(r1,z1,z2,zeta1,b,n,atm) Result(L)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Atmospheres, Only: inv_rho_SL
    Implicit None
    Real(dp) :: L
    Real(dp), Intent(In) :: r1
    Real(dp), Intent(In) :: z1,z2
    Real(dp), Intent(In) :: zeta1
    Integer, Intent(In) :: b,n
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp) :: dZ
    Real(dp) :: S
    Real(dp) :: Ss(1:n)
    Integer :: i

    dZ = z2 - z1
    If (dZ .LE. 0._dp) Then
        L = 0._dp
        Return
    End If
    S = S_upward(r1,dZ,zeta1)
    Ss = atm%EPL_lay(b)%uS * S
    L = 0._dp
    Do i = 1,n
        L = L + atm%EPL_lay(b)%wS(i) * atm%rho(z1 + deltaZ(r1,zeta1,Ss(i)),b)
    End Do
    L = S * 0.5_dp * inv_rho_SL * L
End Function EPL_S_partial_layer

Function EPL_S_known_layer(zeta,b,atm) Result(L)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Implicit None
    Real(dp) :: L
    Real(dp), Intent(In) :: zeta
    Integer, Intent(In) :: b
    Type(Atmosphere_Type), Intent(In) :: atm

    L = EPL_S_partial_layer(atm%Rb(b-1),atm%Zb(b-1),atm%Zb(b),zeta,b,atm%EPL_lay(b)%nS,atm)
End Function EPL_S_known_layer

Function EPL_Z_partial_layer(r1,z1,z2,zeta1,b1,n,atm) Result(L)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Atmospheres, Only: inv_rho_SL
    Implicit None
    Real(dp) :: L
    Real(dp), Intent(In) :: r1
    Real(dp), Intent(In) :: z1,z2
    Real(dp), Intent(In) :: zeta1
    Integer, Intent(In) :: b1,n
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp) :: dZ,A,C
    Real(dp) :: z(1:n)
    Real(dp) :: B(1:n),D(1:n)
    Integer :: i

    dZ = z2 - z1
    If (dZ .LE. 0._dp) Then
        L = 0._dp
        Return
    End If
    A = 0.5_dp * inv_rho_SL * dZ
    z = atm%EPL_lay(b1)%uZ * dZ
    Do i = 1,n
        B(i) = atm%EPL_lay(b1)%wZ(i) * (r1 + z(i)) * atm%rho(z1 + z(i),b1)
    End Do
    C = r1
    D = z * (2._dp * r1 + z)
    L = EPL_Z_full_layer(zeta1,n,A,B,C,D)
End Function EPL_Z_partial_layer

Function EPL_Z_known_layer(lay,zeta) Result(L)
    Use Kinds, Only: dp
    Use Atmospheres, Only: EPL_Layer_Data
    Implicit None
    Real(dp) :: L
    Type(EPL_Layer_Data), Intent(In) :: lay
    Real(dp), Intent(In) :: zeta

    L = EPL_Z_full_layer(zeta,lay%nZ,lay%A,lay%B,lay%C,lay%D)
End Function EPL_Z_known_layer

Function EPL_Z_full_layer(zeta,n,A,B,C,D) Result(L)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: L
    Real(dp), Intent(In) :: zeta
    Integer, Intent(In) :: n
    Real(dp), Intent(In) :: A
    Real(dp), Intent(In) :: B(1:n)
    Real(dp), Intent(In) :: C
    Real(dp), Intent(In) :: D(1:n)

    L = A * Sum( B / Sqrt( (C*zeta)**2 + D) )
End Function EPL_Z_full_layer

Function EPL_Tk_partial_layer(z1,z2,p,e,b,n,atm) Result(L)
    Use Kinds, Only: dp
    Use Global, Only: Rc => R_center
    Use Global, Only: mu => grav_param
    Use Atmospheres, Only: Atmosphere_Type
    Use Atmospheres, Only: inv_rho_SL
    Implicit None
    Real(dp) :: L
    Real(dp), Intent(In) :: z1,z2
    Real(dp), Intent(In) :: p,e
    Integer, Intent(In) :: b,n
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp) :: r1,r2
    Real(dp) :: t1,t2,dt
    Real(dp) :: x1,x2
    Real(dp) :: cos_ts(1:n)
    Real(dp) :: dZs(1:n)
    Integer :: i

    !UNDONE small delta theta will result in poorly conditioned limits of integration
    !UNDONE need to find well conditioned formulae for limits or change variables
    r1 = Rc + z1
    r2 = Rc + z2
    x1 = (p - r1) / (e * r1)
    If (x1 .LT. 1._dp) Then
        t1 = ACos(x1)
    Else
        t1 = 0._dp
    End If
    x2 = (p - r2) / (e * r2)
    If (x2 .LT. 1._dp) Then
        t2 = ACos(x2)
    Else
        t2 = 0._dp
    End If
    dt = t2 - t1
    If (dt .LE. 0._dp) Then
        L = 0._dp
        Return
    End If
    cos_ts = Cos(t1 + atm%EPL_lay(b)%uTk * dt)
    dZs = (p / (1._dp + e*cos_ts)) - r1
    L = 0._dp
    Do i = 1,n
        L = L + atm%EPL_lay(b)%wTk(i)*atm%rho(z1 + dZs(i),b)*Sqrt(1._dp + e*e + 2._dp*e*cos_ts(i)) / ((1._dp + e*cos_ts(i))**2)
    End Do
    L = 0.5_dp * dt * p * inv_rho_SL * L
End Function EPL_Tk_partial_layer

Function EPL_Tk_known_layer(p,e,b,atm) Result(L)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Implicit None
    Real(dp) :: L
    Real(dp), Intent(In) :: p,e
    Integer, Intent(In) :: b
    Type(Atmosphere_Type), Intent(In) :: atm

    L = EPL_tk_partial_layer(atm%Zb(b-1),atm%Zb(b),p,e,b,atm%EPL_lay(b)%nTk,atm)
End Function EPL_Tk_known_layer

!NEEDSDONE Switch Kepler EPL calculations from p-e-rp formulations to p-xi form to reduce arithmetic and improve conditioning
Function EPL_Rk_partial_layer(z1,z2,p,e,rp,b1,n,atm) Result(L)
    Use Kinds, Only: dp
    Use Global, Only: Rc => R_center
    Use Atmospheres, Only: Atmosphere_Type
    Use Atmospheres, Only: inv_rho_SL
    Implicit None
    Real(dp) :: L
    Real(dp), Intent(In) :: z1,z2
    Real(dp), Intent(In) :: p,e,rp
    Integer, Intent(In) :: b1,n
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp) :: dZ,A
    Real(dp) :: z(1:n)
    Real(dp) :: B(1:n),C(1:n)
    Integer :: i

    dZ = z2 - z1
    If (dZ .LE. 0._dp) Then
        L = 0._dp
        Return
    End If
    A = 0.5_dp * inv_rho_SL * dZ
    z = atm%EPL_lay(b1)%uRk * dZ
    Do i = 1,n
        B(i) = atm%EPL_lay(b1)%wRk(i) * atm%rho(z1 + z(i),b1)
    End Do
    C = Rc + z1 + z
    L = EPL_Rk_full_layer(p,e,rp,n,A,B,C)
End Function EPL_Rk_partial_layer

Function EPL_Rk_known_layer(lay,p,e,rp) Result(L)
    Use Kinds, Only: dp
    Use Atmospheres, Only: EPL_Layer_Data
    Implicit None
    Real(dp) :: L
    Type(EPL_Layer_Data), Intent(In) :: lay
    Real(dp), Intent(In) :: p,e,rp

    L = EPL_Rk_full_layer(p,e,rp,lay%nRk,lay%A,lay%Bk,lay%Ck)
End Function EPL_Rk_known_layer

Function EPL_Rk_full_layer(p,e,rp,n,A,B,C) Result(L)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: L
    Real(dp), Intent(In) :: p,e,rp
    Integer, Intent(In) :: n
    Real(dp), Intent(In) :: A
    Real(dp), Intent(In) :: B(1:n)
    Real(dp), Intent(In) :: C(1:n)
    Real(dp) :: x(1:n),ra
    Real(dp) :: perpC(1:n)

    If (e .GT. 1._dp) Then
        x = (p - C * (1._dp - e)) * ((1._dp + e) * (C - rp))
    Else
        ra = p / (1._dp - e)
        x = ((1._dp - e)*(ra - C)) * ((1._dp + e) * (C - rp))
    End If
    perpC = Sqrt( (x + p*p) / x )
    L = A * Sum( B * perpC )
End Function EPL_Rk_full_layer

Function L_to_S_exact(atm,L,r0,z0,zeta0,nb,bb,Lb,db) Result(S)
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Global, Only: Rc => R_center
    Implicit None
    Real(dp) :: s
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp), Intent(In) :: L
    Real(dp), Intent(In) :: r0
    Real(dp), Intent(In) :: z0
    Real(dp), Intent(In) :: zeta0
    Integer, Intent(In) :: nb
    Integer, Intent(In) :: bb(1:nb)
    Real(dp), Intent(In) :: Lb(1:nb)
    Logical, Intent(In) :: db(1:nb)
    Real(dp) :: r1,z1,zeta1  !defines bottom (upward looking) end of ray segment
    Real(dp) :: r2,z2,zeta2  !defines top (downward looking) end of ray segment
    Real(dp) :: S0,Sb
    Real(dp) :: L0
    Integer :: i
    Real(dp) :: r_ca

    !find the segment on which the desired S falls
    Do i = 1,nb-1
        If (Sum(Lb(1:i)) .GT. L) Exit
    End Do
    If (db(i)) Then !downward segment
        !upward looking end of segment
        If (i .EQ. nb) Then
            !upward looking end of segment is at lower layer boundary
            r1 = atm%Rb(bb(i)-1)
            z1 = atm%Zb(bb(i)-1)
            zeta1 = -zeta_downward(r0,r1,zeta0)
        Else  !check for shallow (down-up) ray
            If (db(i+1)) Then !following segment is also downward, no down-up condition
                !upward looking end of segment is at lower layer boundary
                r1 = atm%Rb(bb(i)-1)
                z1 = atm%Zb(bb(i)-1)
                zeta1 = -zeta_downward(r0,r1,zeta0)
            Else
                !upward looking end of segment is at point of closest approach
                r1 = R_close_approach(r0,zeta0)
                z1 = r1 - Rc
                zeta1 = 0._dp
            End If
        End If
        !downward looking end of segment
        If (i .EQ. 1) Then
            !downward looking end of segment begins at r0,z0,zeta0
            r2 = r0
            z2 = z0
            zeta2 = zeta0
            !
            S0 = 0._dp
            L0 = 0._dp
        Else
            !downward looking end of segment begins at upper layer boundary
            r2 = atm%Rb(bb(i))
            z2 = atm%Zb(bb(i))
            zeta2 = zeta_downward(r0,r2,zeta0)
            !
            S0 = S_upward(r2,r0-r2,-zeta2)
            L0 = Sum(Lb(1:i-1))
        End If
        Sb = S_upward(r1,z2-z1,zeta1)
        s = S0 + L_to_S_from_top_of_segment(atm,bb(i),L-L0,Lb(i),Sb,r2,z2,zeta2)
    Else !upward segment
        !downward looking end of segment always at upper layer boundary
        r2 = atm%Rb(bb(i))
        z2 = atm%Zb(bb(i))
        zeta2 = -zeta_upward(r0,r2-r0,zeta0)
        !upward looking end of segment
        If (i .EQ. 1) Then
            !upward looking end of segment begins at r0,z0,zeta0
            r1 = r0
            z1 = z0
            zeta1 = zeta0
            !
            S0 = S_upward(r0,r2-r0,zeta0)
        Else If (Any(db)) Then !check for shallow (down-up) ray
            r_ca = R_close_approach(r0,zeta0)
            If (db(i-1)) Then  !preceeding segment is downward, down-up condition
                !upward looking end of segment is at point of closest approach
                r1 = r_ca
                z1 = r1 - Rc
                zeta1 = 0._dp
            Else  !preceeding segment is also upward, no down-up condition on segement, but still down-up on ray
                !upward looking end of segment begins at lower layer boundary
                r1 = atm%Rb(bb(i)-1)
                z1 = atm%Zb(bb(i)-1)
                zeta1 = zeta_upward(r0,z1-z0,zeta0)            
            End If
            S0 = -zeta0 * r0 + S_upward(r_ca,r2-r_ca,0._dp)
        Else
            !upward looking end of segment begins at lower layer boundary
            r1 = atm%Rb(bb(i)-1)
            z1 = atm%Zb(bb(i)-1)
            zeta1 = zeta_upward(r0,z1-z0,zeta0)
            !
            S0 = S_upward(r0,r2-r0,zeta0)
        End If
        !
        L0 = Sum(Lb(1:i))
        Sb = S_upward(r1,r2-r1,zeta1)
        s = S0 - L_to_S_from_top_of_segment(atm,bb(i),L0-L,Lb(i),Sb,r2,z2,zeta2)
    End If
End Function L_to_S_exact

Function L_to_S_from_top_of_segment(atm,b,L,Lmax,Smax,r1,z1,zeta1) Result(S)
    !Given a DOWNWARD ray in layer b with zenith cosine zeta1 at r1 and z1, 
    ! finds the distance S along the ray at which optical thickness L
    ! is traversed
    Use Kinds, Only: dp
    Use Atmospheres, Only: Atmosphere_Type
    Use Atmospheres, Only: rho_SL
    Use Utilities, Only: Converged
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Real(dp) :: S
    Type(Atmosphere_Type), Intent(In) :: atm
    Integer, Intent(In) :: b
    Real(dp), Intent(In) :: L !from TOP of segment
    Real(dp), Intent(In) :: Lmax
    Real(dp), Intent(In) :: Smax
    Real(dp), Intent(In) :: r1
    Real(dp), Intent(In) :: z1
    Real(dp), Intent(In) :: zeta1 !always negative
    Real(dp) :: L_high,L_low
    Real(dp) :: S_high,S_low
    Real(dp) :: S_old,L_this_S
    Real(dp) :: dZ,zeta
    Integer :: n_fix
    Integer :: i

    !initialize bracketing variables
    S_high = Smax
    S_low = 0._dp
    L_high = Lmax
    L_low = 0._dp
    n_fix = 0
    !initial guess for S
    S = 0._dp
    Do i = 1,100
        S_old = S
        dZ = deltaZ(r1,zeta1,S)  !change in height from z1 after traveling S, always negative
        zeta = -zeta_downward(r1,r1+dZ,zeta1)  !zenith cosine looking back UP ray from S, always positive
        !Compute EPL for this guess of S
        If (zeta .GE. 0.1_dp) Then
            L_this_S = EPL_Z_partial_layer(r1+dZ,z1+dZ,z1,zeta,b,atm%EPL_lay(b)%nZ,atm)  !EPL from zero to S
        Else
            L_this_S = EPL_S_partial_layer(r1+dZ,z1+dZ,z1,zeta,b,atm%EPL_lay(b)%nS,atm)  !EPL from zero to S
        End If
        !Newton's method for next s
        S = S_old - rho_SL * (L_this_S - L) / atm%rho(z1+dZ,b)
        !L is monotonically increasing, so we can refine the brackets in which the solution lies
        If (L_this_S.LT.L .AND. S_old.GT.S_low) Then  !solution is greater than this S_up
            S_low = S_old
            L_low = L_this_S
        Else If (L_this_S.GT.L .AND. S_old.LT.S_high) Then  !solution is less than this S_up
            S_high = S_old
            L_high = L_this_S
        End If
        !If this iteration landed outside known bounds, Newton's method has diverged
        If (S.GT.S_high .OR. S.LT.S_low) Then
            If (n_fix .GT. 10) Exit  !Newton w/ false position fix has failed 10 times, UNCONVERGED EXIT
            !Fix Newton with a False Position guess for this iteration
            S = S_old - (L_this_S - L) * (S_high - S_low) / (L_high - L_low)
            n_fix = n_fix + 1
        End If
        !Check for convergence
        If (Converged(S,S_old,rTol = 1.E-9_dp,aTol = 1.E-9_dp)) Return  !NORMAL EXIT for convergence on S
        If (Converged(L,L_this_S,rTol = 1.E-12_dp,aTol = 1.E-12_dp)) Return  !NORMAL EXIT for guessing correct L
    End Do
    !If we get this far, Newton with False Position helpler did not converge
    ! Use bisection instead
    Do i = 1,100
        S = 0.5_dp * (S_low + S_high)
        If (Converged(S,S_old,rTol = 1.E-9_dp,aTol = 1.E-9_dp)) Return  !NORMAL EXIT for convergence on S_up
        dZ = deltaZ(r1,zeta1,S)
        zeta = -zeta_downward(r1,r1+dZ,zeta1)
        !Compute EPL for this guess of s
        If (zeta .GE. 0.1_dp) Then
            L_this_S = EPL_Z_partial_layer(r1+dZ,z1+dZ,z1,zeta,b,atm%EPL_lay(b)%nZ,atm)
        Else
            L_this_S = EPL_S_partial_layer(r1+dZ,z1+dZ,z1,zeta,b,atm%EPL_lay(b)%nS,atm)
        End If
        If (L_this_S .GT. L) Then !guess for S_up was too high
            S_high = S
        Else If (L_this_S .LT. L) Then  !guess for S_up was too low
            S_low = S
            L_low = L_this_S
        Else !guess for S_up gave exact L_up
            Return  !NORMAL EXIT for guessing correct L
        End If
        S_old = S
    End Do
    !If we get this far, we failed to converge using Newton AND bisection
    Call Output_Message('ERROR:  Pathlengths: L_to_S_from_top_of_segment:  Failed to converge on a root for S.',kill=.TRUE.)
End Function L_to_S_from_top_of_segment

End Module Pathlengths
