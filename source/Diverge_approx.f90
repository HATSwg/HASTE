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
Module Diverge_approx

    Implicit None
    Private
    Public :: Div_Fact_Straight
    Public :: Div_Fact_by_shooting

Contains

Function Div_Fact_Straight(r1,v1,r2,u,tof,vS2) Result(D)
    Use Kinds, Only: dp
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Unit_Vector
    Implicit None
    Real(dp) :: D
    Real(dp), Intent(In) :: r1(1:3),v1(1:3),r2(1:3)
    Real(dp), Intent(In) :: u(1:3)
    Real(dp), Intent(In) :: tof
    Real(dp), Intent(In) :: vS2(1:3)
    Real(dp) :: r1_to_r2(1:3),dist_r1_to_r2
    Real(dp) :: cos_theta

    r1_to_r2 = r2 - (r1 + u*tof)
    dist_r1_to_r2 = Vector_Length(r1_to_r2)
    D = 1._dp / ( dist_r1_to_r2**2 )
    cos_theta = Max( Abs( Dot_Product(Unit_Vector(vS2 - v1),r1_to_r2 / dist_r1_to_r2) ) , 0.1_dp )
    D = D / cos_theta
End Function Div_Fact_Straight

Function Div_Fact_by_shooting(r1,Omega_hat1cm,s1cm,u,tof,vS2,v2) Result(D)
    Use Kinds, Only: dp
    Use Utilities, Only: Cross_Product
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Quadrilateral_Area
    Use Utilities, Only: Normal_Vector_4
    Use Global, Only: Z_hat,X_hat
    Use Astro_Utilities, Only: Kepler_Gooding
    Implicit None
    Real(dp) :: D
    Real(dp), Intent(In) :: r1(1:3)
    Real(dp), Intent(In) :: Omega_hat1cm(1:3)
    Real(dp), Intent(In) :: s1cm
    Real(dp), Intent(In) :: u(1:3)
    Real(dp), Intent(In) :: tof
    Real(dp), Intent(In) :: vS2(1:3)
    Real(dp), Intent(In) :: v2(1:3)
    Real(dp) :: R_hat(1:3),N_hat(1:3),T_hat(1:3)
    Integer :: i
    Real(dp) :: Omega_hat1cm_shots(1:3,1:4)
    Real(dp) :: r2_shots(1:3,1:4)
    Real(dp) :: v2_shots(1:3,1:4)
    Real(dp) :: Acm,Aar
    Real(dp) :: cos_theta
    Real(dp) :: Surf_Normal(1:3)
    Real(dp) :: mag_Surf_Normal
    Real(dp), Parameter :: eps = 1.E-5_dp

    !construct basis in which to perturb launch direction (Radial, Normal, Transverse directions IN CM FRAME)
    R_hat = Unit_Vector(r1)
    If (Abs(Dot_Product(Omega_hat1cm,R_hat)) .LT. 1._dp) Then
        N_hat = Unit_Vector(Cross_Product(R_hat,Omega_hat1cm))
        T_hat = Unit_Vector(Cross_Product(N_hat,Omega_hat1cm))
    Else !emission direction is straight up, arbitrarily choose normal and transverse directions
        If (Abs(Dot_Product(Z_hat,R_hat)) .LT. 0.5_dp) Then
            T_hat = Unit_Vector(Cross_Product(R_hat,Z_hat))
        Else
            T_hat = Unit_Vector(Cross_Product(R_hat,X_hat))
        End If
        N_hat = Cross_Product(T_hat,Omega_hat1cm)
    End If
    !create four new launch directions by pertubing normal and radial/transverse components of Omega_hat1cm
    Omega_hat1cm_shots(:,1) = Unit_Vector( Omega_hat1cm + eps*N_hat )
    Omega_hat1cm_shots(:,3) = Unit_Vector( Omega_hat1cm - eps*N_hat )
    Omega_hat1cm_shots(:,2) = Unit_Vector( Omega_hat1cm + eps*T_hat )
    Omega_hat1cm_shots(:,4) = Unit_Vector( Omega_hat1cm - eps*T_hat )
    !propegate trajectories along the four test directions
    Do CONCURRENT (i = 1:4)
        Call Kepler_Gooding( r1,Omega_hat1cm_shots(:,i)*s1cm + u,tof,r2_shots(:,i),v2_shots(:,i) )
    End Do
    !Compute area as appoximation of solid angle at emission
    Acm = Quadrilateral_Area(Omega_hat1cm_shots)
    !Compute area as approximation of solid angle at rendezvous, accounting for incidence angle relative to detector motion
    Surf_Normal = Normal_Vector_4(r2_shots)
    mag_Surf_Normal = Vector_Length(Surf_Normal)
    Aar = 0.5_dp * mag_Surf_Normal !same as Quadrilateral_Area(r2_shots), but preserves intermediate quatity for later use
    !Divergence factor is the ratio of these areas adjusted for angle of incidence
    D = Acm / Aar
    cos_theta = Max( Abs( Dot_Product(Unit_Vector(vS2 - v2),Surf_Normal / mag_Surf_Normal) ) , 0.1_dp )
    D = D / cos_theta
End Function Div_Fact_by_shooting

End Module Diverge_approx
