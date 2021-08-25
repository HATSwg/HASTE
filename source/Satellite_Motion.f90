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
Module Satellite_Motion

    Use Kinds, Only: dp
    Implicit none

! Cartesian coordinate basis vectors are inertial
!   and are aligned with the rotational position of the Earth at time t = 0:
!   XHat = (1,0,0)  Points toward the equator at 0 deg Longitude (Greenwich meridian)
!   YHat = (0,1,0)  Points toward the equator at 90 deg E Longitude
!   ZHat = (0,0,1)  Points toward the north pole
! This system is used so that all the coordinate transformations
!   that are needed can be done once,
!   in establishing the satellites orbital parameters in this system,
!   and are not required while tracking neutrons in the simulation.

    Private
    Public :: Satellite_Position_Type
    Public :: Initialize_Satellite_Motion

    Type Satellite_Position_Type
        Real(dp) :: r0(1:3)  !position vector at t=0
        Real(dp) :: v0(1:3)  !velocity vector at t=0
        Real(dp) :: rp
        Real(dp) :: vp
        Logical :: is_stationary
        Logical :: is_conic
        Logical :: is_smooth
        Integer :: nt
        Real(dp), Allocatable :: ts(:)
        Real(dp), Allocatable :: rs_vs(:,:)
    Contains
        Procedure, Pass :: R_and_V => R_V_Satellite
        Procedure, Pass :: R => R_Satellite
        Procedure, Pass :: Cleanup => Cleanup_Satellite_position
    End Type

Contains

Subroutine Initialize_Satellite_Motion(resources_dir,motion_type,sat)
    Use Kinds, Only: dp
    Use Global, Only: mu => grav_param
    Use Global, Only: rot_Earth
    Use Global, Only: Z_hat
    Use Astro_Utilities, Only: Radius_of_Periapsis
    Use Astro_Utilities, Only: Velocity_of_Periapsis
    Use Astro_Utilities, Only: Period
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Use Utilities, Only: Cube_Root
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Character(*), Intent(In) :: resources_dir
    Character(10), Intent(In) :: motion_type
    Type(Satellite_Position_Type), Intent(InOut) :: sat

    Select Case(motion_type)
        Case('Stationary')
            sat%is_stationary = .TRUE.
            sat%is_conic = .FALSE.
            sat%is_smooth = .TRUE.
            sat%rp = Vector_Length(sat%r0)
            sat%vp = Vector_Length(sat%v0)
        Case('Linear')
            sat%is_stationary = .FALSE.
            sat%is_conic = .FALSE.
            sat%is_smooth = .TRUE.
            If (Dot_Product(sat%r0,sat%v0) .LT. 0._dp) Then !downward path
                sat%rp = -Vector_Length(sat%r0) * Dot_Product(sat%r0,sat%v0)
            Else !upward path
                sat%rp = Vector_Length(sat%r0)
            End If
            sat%vp = Vector_Length(sat%v0)
        Case('Conic')
            sat%is_stationary = .FALSE.
            sat%is_conic = .TRUE.
            sat%is_smooth = .TRUE.
            sat%rp = Radius_of_Periapsis(sat%r0,sat%v0)
            sat%vp = Velocity_of_Periapsis(sat%r0,sat%v0)
        Case ('Conic_tab')
            sat%is_stationary = .FALSE.
            sat%is_conic = .TRUE.
            sat%is_smooth = .FALSE.
            Call Read_sat_trace(resources_dir,sat%nt,sat%ts,sat%rs_vs)
            sat%r0 = sat%rs_vs(1:3,1)
            sat%v0 = sat%rs_vs(4:6,1)
        Case('GeoStat')
#           if LUNA
                ! A lunar-stationary orbit above the surface of the moon is not feasible (necessary altitude for this type of lunar 
                ! orbit places the satellite outside the moon's gravitational sphere of influence)
                Call Output_Message('ERROR:  Satellite_Motion: Initialize_Satellite_Motion:  Stationary lunar orbit specified.', & 
                                   & kill=.TRUE.)
#           else
                sat%is_stationary = .FALSE.
                sat%is_conic = .TRUE.
                sat%is_smooth = .TRUE.
                sat%r0(3) = 0._dp  !make sure r0 is in equatorial plane
                sat%r0 = Unit_Vector(sat%r0) * Cube_Root(mu / rot_Earth**2)  !make sure r0 has correct magnitude
                sat%v0 = Unit_Vector(Cross_Product(Z_hat,sat%r0)) * Cube_Root(mu / rot_Earth**2) * rot_Earth  !set v0
                sat%rp = Vector_Length(sat%r0)
                sat%vp = Velocity_of_Periapsis(sat%r0,sat%v0)
#           endif
        Case('StatGeo')  !special undocumented case:
            ! STATIONARY satellite in geostationary position with geostationary velocity, but NO MOTION
            sat%is_stationary = .TRUE.
            sat%is_conic = .FALSE.
            sat%is_smooth = .TRUE.
            sat%r0(3) = 0._dp  !make sure r0 is in equatorial plane
            sat%r0 = Unit_Vector(sat%r0) * Cube_Root(mu / rot_Earth**2)  !make sure r0 has correct magnitude
            sat%v0 = Unit_Vector(Cross_Product(Z_hat,sat%r0)) * Cube_Root(mu / rot_Earth**2) * rot_Earth  !set v0
            sat%rp = Vector_Length(sat%r0)
            sat%vp = Vector_Length(sat%v0)
        Case Default
            Call Output_Message('ERROR:  Satellite_Motion: Initialize_Satellite_Motion:  Unknown motion type.',kill=.TRUE.)
    End Select
End Subroutine Initialize_Satellite_Motion

Subroutine R_V_Satellite(s,t,r,v)
    Use Kinds, Only: dp
    Use Astro_Utilities, Only: Kepler_Gooding
    Implicit none
    Class(Satellite_Position_Type), Intent(In) :: s
    Real(dp), Intent(In) :: t
    Real(dp), Intent(Out) :: r(1:3),v(1:3)

    If (s%is_conic) Then
        If (s%is_smooth) Then
            Call Kepler_Gooding(s%r0,s%v0,t,r,v)
        Else
            Call Smooth_tabular_conic(s,t,r,v)
        End If
    Else If (s%is_stationary) Then
        r = s%r0
        v = s%v0
    Else !s%is_linear
        r = s%r0 + t*s%v0
        v = s%v0
    End If
End Subroutine R_V_Satellite

Function R_Satellite(s,t) Result(r)
    Use Kinds, Only: dp
    Use Astro_Utilities, Only: Kepler_Gooding
    Implicit none
    Real(dp) :: r(1:3)
    Class(Satellite_Position_Type), Intent(In) :: s
    Real(dp), Intent(In) :: t
    Real(dp) :: v(1:3)

    If (s%is_conic) Then
        If (s%is_smooth) Then
            Call Kepler_Gooding(s%r0,s%v0,t,r,v)
        Else
            Call Smooth_tabular_conic(s,t,r,v)
        End If
    Else If (s%is_stationary) Then
        r = s%r0
    Else !s%is_linear
        r = s%r0 + t*s%v0
    End If
End Function R_Satellite

Subroutine Smooth_tabular_conic(s,t,r,v)
    Use Kinds, Only: dp
    Use Astro_Utilities, Only: Kepler_Gooding
    Use Utilities, Only: Bisection_Search
    Implicit None
    Class(Satellite_Position_Type), Intent(In) :: s
    Real(dp), Intent(In) :: t
    Real(dp), Intent(Out) :: r(1:3),v(1:3)
    Integer :: i
    Real(dp) :: dt,w,time_from_last,time_till_next
    Real(dp) :: ra(1:3),va(1:3)
    Real(dp) :: rb(1:3),vb(1:3)

    !If the time of interest is beyond the extent of the list, use the first or last state vector as appropriate
    If (t .LE. 0._dp) Then
        Call Kepler_Gooding(s%rs_vs(1:3,1),s%rs_vs(4:6,1),t,r,v)
    Else If (t .GE. s%ts(s%nt)) Then
        Call Kepler_Gooding(s%rs_vs(1:3,s%nt),s%rs_vs(4:6,s%nt),t-s%ts(s%nt),r,v)
    Else !Otherwise, smoothly blend solutions from the two nearest (in time) state vectors
        i = Bisection_Search(t,s%ts,s%nt)
        dt = s%ts(i) - s%ts(i-1)  !time between state vectors
        time_from_last = t - s%ts(i-1)  !seconds since the prev state vector (positive)
        time_till_next = time_from_last - dt  !seconds until the next state vector (negative)
        Call Kepler_Gooding(s%rs_vs(1:3,i-1),s%rs_vs(4:6,i-1),time_from_last,ra,va) !solution from the prev state vector
        Call Kepler_Gooding(s%rs_vs(1:3,i  ),s%rs_vs(4:6,i  ),time_till_next,rb,vb) !solution from the next state vector
        w = Blend3(time_from_last / dt)  !weight factor for the next state vector's solution (prev state vector weight is 1-w)
        r = w*rb + (1._dp - w)*ra
        v = w*vb + (1._dp - w)*va
    End If
End Subroutine Smooth_tabular_conic

Pure Function Blend1(x) Result(w)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: w
    Real(dp), Intent(In) :: x

    w = x
End Function Blend1

Pure Function Blend3(x) Result(w)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: w
    Real(dp), Intent(In) :: x

    w = x * x * (3._dp - 2._dp * x)
End Function Blend3

Pure Function Blend5(x) Result(w)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: w
    Real(dp), Intent(In) :: x

    w = x * x * (5._dp + 2._dp*x * (x * (5._dp - 2._dp*x) - 5._dp))
End Function Blend5

Subroutine Cleanup_Satellite_Position(s)
    Use Kinds, Only: dp
    Implicit None
    Class(Satellite_Position_Type), Intent(InOut) :: s

    s%r0 = 0._dp
    s%v0 = Huge(s%v0)
    s%rp = 0._dp
    s%vp = Huge(s%vp)
    s%nt = 0
    If (Allocated(s%ts)) Deallocate(s%ts)
    If (Allocated(s%rs_vs)) Deallocate(s%rs_vs)
End Subroutine Cleanup_Satellite_Position

Subroutine Read_sat_trace(resources_dir,n,ts,rs_vs)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Character(*), Intent(In) :: resources_dir
    Integer, Intent(Out) :: n
    Real(dp), Intent(Out), Allocatable :: ts(:)
    Real(dp), Intent(Out), Allocatable :: rs_vs(:,:)
    Integer :: trace_unit,stat
    Integer :: i

    Open( NEWUNIT = trace_unit , FILE = resources_dir//'sat'//slash//'sat_trace.txt' , STATUS = 'OLD' , ACTION = 'READ' , & 
        & IOSTAT = stat )
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Satellite Motion: Read_sat_trace:  File open error, '//resources_dir// & 
        & 'sat'//slash//'sat_trace.txt'//', IOSTAT=',stat,kill=.TRUE.)
    !Count the number of lines in the file
    n = 0
    Do
        Read(trace_unit,*,IOSTAT=stat)
        If (stat .LT. 0) Exit
        n = n + 1
    End Do
    Rewind(trace_unit)
    Allocate(ts(1:n))
    ts = 0._dp
    Allocate(rs_vs(1:6,1:n))
    rs_vs = 0._dp
    Do i = 1,n
        Read(trace_unit,'(7ES27.16E3)') ts(i),rs_vs(1,i),rs_vs(2,i),rs_vs(3,i),rs_vs(4,i),rs_vs(5,i),rs_vs(6,i)
    End Do
    Close(trace_unit)
    !Set t(1) as epoch, t=0
    ts = ts - ts(1)
End Subroutine Read_sat_trace

End Module Satellite_Motion
