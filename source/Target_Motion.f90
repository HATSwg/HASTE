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
Module Target_Motion

    Use Kinds, Only: dp
    Use Global, Only: k_Boltzmann
    Use Global, Only: keV_per_Joule
    Implicit None
    Private
    Public :: Random_Thermal_Velocity_AF
    Public :: Atm_Rotation_Velocity_AF

    Interface Random_Thermal_Velocity_AF
        Module Procedure Random_Thermal_Velocity_Monatomic  ! (T, An)
        Module Procedure Random_Thermal_Velocity_Diatomic   ! (T, An, AnPrime, Mn)
    End Interface Random_Thermal_Velocity_AF

    Real(dp), Parameter :: k_B = k_Boltzmann * keV_per_Joule  ![keV/K] Boltzmann constant with convenient units locally

        ! Note, CO2 can be approximated as diatomic (less accurately than for N2 and O2)
        !   with An and AnPrime having values for the O atoms (for rotation)
        !   but with Mn including the C as well (for translation).
        !   For this reason, Mn is an argument of Random_Thermal_Velocity_Diatomic.
        !   This lets us include all the usual constituents except H2O.
        !   We exclude H2O because it is significant only in the troposphere,
        !   where it depends on weather.

Contains

Function Atm_Rotation_Velocity_AF(big_r,xi) Result(V)
    Use Kinds, Only: dp
    Use Global, Only: rot_Earth
    Implicit None
    Real(dp) :: V(1:3)
    Real(dp), Intent(In) :: big_r  !distance from center of Earth
    Real(dp), Intent(In) :: xi  !cosine of lattitude

    V = (/ rot_Earth * big_r * Sqrt(1._dp - xi**2), &
         & 0._dp, &
         & 0._dp /)
End Function Atm_Rotation_Velocity_AF

Function Random_Thermal_Velocity_Monatomic(T,An,RNG) Result(v)
        ! Returns the thermal velocity of an atom in a diatomic molecule
        ! in [km/s] in the rest frame of the air, for air at the temperature at the altitude of
        ! point R (in ECI coordinates)
    Use Kinds, Only: dp
    Use Random_Numbers, Only: RNG_Type
    Use Random_Directions, Only: Isotropic_Omega_Hat
    Implicit None
    Real(dp):: v(1:3)   ! [km/s] Thermal velocity
    Real(dp), Intent(In) :: T   ! [K] Temperature
    Real(dp), Intent(In) :: An  ! [] mass of target atom / neutron mass
    Type(RNG_Type), Intent(InOut) :: RNG
    Real(dp) :: st   ! [km/s] speed of atom in rest frame of air due to thermal translation of molecule

    st = Random_Thermal_Translational_Speed(T,An,RNG)
    v= st * Isotropic_Omega_Hat(RNG)
End Function Random_Thermal_Velocity_Monatomic

Function Random_Thermal_Velocity_Diatomic(T, An, AnPrime, Mn, RNG) Result(v)
        ! Returns the thermal velocity of an atom in a diatomic molecule
        ! in [km/s] in the rest frame of the air, for air at the temperature at the altitude of
        ! point R (in ECI coordinates)
    Use Kinds, Only: dp
    Use Random_Numbers, Only: RNG_Type
    Use Random_Directions, Only: Isotropic_Omega_Hat
    Implicit None
    Real(dp):: v(1:3)   ! [km/s] Thermal velocity
    Real(dp), Intent(In) :: T   ! [K] Temperature
    Real(dp), Intent(In):: An       ! [] mass of target atom / neutron mass
    Real(dp), Intent(In):: AnPrime  ! [] mass of other atom in diatomic molecule
    Real(dp), Intent(In):: Mn       ! [] mass of molecule
    Type(RNG_Type), Intent(InOut) :: RNG
    Real(dp):: sr   ! [km/s] speed of atom in CM frame due to thermal rotation of molecule
    Real(dp):: st   ! [km/s] speed of atom in rest frame of air due to thermal translation of molecule
    Real(dp):: mu   ! [] cosine of angle between directions of rotational velocity of atom
                    !   and translational velocity of atom uniformly distributed in [-1,+1]
    Real(dp):: s    ! [km/s] speed of target atom resulting from thermal translation
                    !   combined with thermal rotation of diatomic molecule

    st = Random_Thermal_Translational_Speed(T, Mn, RNG)
    sr = Random_Thermal_Rotational_Speed(T, An, AnPrime, Mn, RNG)
    mu = 2._dp * RNG%Get_Random() - 1._dp
    s = Sqrt(st**2 + sr**2 + 2._dp * mu * st * sr)
    v = s * Isotropic_Omega_Hat(RNG)
End Function Random_Thermal_Velocity_Diatomic

Function Random_Thermal_Rotational_Speed(T, An, AnPrime, Mn, RNG) Result(speed)
    ! Samples random speed of one atom of a diatomic molecule (or CO2) 
    !   due to thermally-excited rotation 
    ! Result is in the CM frame of the molecule
    Use Kinds, Only: dp
    Use Random_Numbers, Only: RNG_Type
    Use Random_Directions, Only: Isotropic_Omega_Hat
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp):: speed            ! [km/s] random thermal motion speed
    Real(dp), Intent(In):: T    ! [degree_Kelvin] temperature
    Real(dp), Intent(In):: An   ! [] mass of atom of interest / neutron mass
    Real(dp), Intent(In):: AnPrime  ! [] An of the other atom in diatomic molecule
    Real(dp), Intent(In):: Mn   ! [] mass of molecule / neutron mass (An+AnPrime)
    Type(RNG_Type), Intent(InOut) :: RNG
    Real(dp):: Er               ! [keV] kinetic energy of the atom of interest
                                ! due to rotation of the diatomic molecule
    Real(dp):: zr               ! [] nondimensionalized rotational energy of molecule

    zr = -Log(1._dp - RNG%Get_Random())
    Er = (AnPrime / Mn) * zr * k_B * T
    speed = Neutron_Speed(Er / An)
End Function Random_Thermal_Rotational_Speed

Function Random_Thermal_Translational_Speed(T, Mn, RNG) Result(speed)
    ! Samples random speed of translation of CM of molecule of perfect gas
    !   in thermodynamic equilibrium at temperature T [Kelvin]
    ! Result is in the local rest frame of the air
    Use Kinds, Only: dp
    Use Random_Numbers, Only: RNG_Type
    Use Random_Directions, Only: Isotropic_Omega_Hat
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp):: speed            ! [km/s] random thermal motion speed
    Real(dp), Intent(In):: T    ! [Kelvin] temperature
    Real(dp), Intent(In):: Mn   ! [] mass of molecule / neutron mass
    Type(RNG_Type), Intent(InOut) :: RNG
    Real(dp):: Et               ! [keV] kinetic energy of translation of molecule's CM
                                !   in local rest frame of air
    Real(dp):: zt               ! [] Nondimensionalized speed = Sqrt(Et/(k*T))

    zt = Sample_Nondimensional_Translational_Speed(RNG)
    Et = zt**2 * k_B * T
    speed = Neutron_Speed(Et / Mn)
End Function Random_Thermal_Translational_Speed

Function Sample_Nondimensional_Translational_Speed(RNG) Result(zt)
    Use Kinds, Only: dp
    Use Random_Numbers, Only: RNG_Type
    Implicit None
    Real(dp):: zt               ! [] Nondimensionalized speed = Sqrt(Et/(k*T))
    Type(RNG_Type), Intent(InOut) :: RNG
    Real(dp):: r0, r(1:2)       ! [] random numbers from U[0,1)
    Real(dp), Parameter:: sp1 = 0.03901746166235435_dp   ! Probability of zt < 48/125
    Real(dp), Parameter:: sp2 = 0.29185870792747644_dp   ! Probability of zt < 5/6
    Real(dp), Parameter:: sp3 = 0.55827516574364871_dp   ! Probability of zt < 29/25
    Real(dp), Parameter:: sp4 = 0.90954539963811080_dp   ! Probability of zt < 9/5

    ! Use sampling by rejection with piecewise sampling function (s1 through s5)
    ! Note that probabilities p1 through p4, (p5 = 1-p1-p2-p3-p4 being implicit)
    !   are probabilites that zt is in the corresponding interval, 
    !   obtained by normalizing g(zt) to make f(zt) and integrating f(zt) over each interval.
    !   These are summed: sp1 = p1, sp2 = p1+p2, and so on.
    ! Therefore, the interval need be randomly chosen only once. 
    ! Then the rejection sampling for that interval is repeated until a sample is accepted.
    ! The efficiencies of the schemes, to the nearest 0.1%, are
    ! s1: 91.6%
    ! s2: 96.8%
    ! s3: 98.2%
    ! s4: 98.5%
    ! s5: 78.6%
    ! Weighting each by the probability that it is used (p1 through p5),
    !   the overall efficiency of this scheme is 95.9%.
    ! This scheme was developed by Professor Kirk Mathews, author of this code.
    r0 = RNG%Get_Random()
    If (r0 == 0.0_dp) Then
        zt = 0.0_dp
    Else if (r0 < sp1) then  ! 0 < zt < 48/125
        Do
            r = RNG%Get_Randoms(2)
            zt = S1Inverse(r(1))
            If (r(2) * s1(zt) <= g(zt)) Exit  ! Accept zt
        End do
    Else if (r0 < sp2) then  ! 48/125 <= zt < 5/6
        Do
            r = RNG%Get_Randoms(2)
            zt = S2Inverse(r(1))
            If (r(2) * s2(zt) <= g(zt)) Exit  ! Accept zt
        End do
    Else if (r0 < sp3) then  ! 5/6 <= zt < 29/25
        Do
            r = RNG%Get_Randoms(2)
            zt = S3Inverse(r(1))
            If (r(2) <= g(zt)) Exit  ! Accept zt 
            !   Note s3(zt) = 1.0_dp so no need to include it in the logical test
        End do
    Else if (r0 < sp4) then  ! 29/25 <= zt < 9/5
        Do
            r = RNG%Get_Randoms(2)
            zt = S4Inverse(r(1))
            If (r(2) * s4(zt) <= g(zt)) Exit  ! Accept zt
        End do
    Else                    ! 9/5 <= zt < Infinity
        Do
            r = RNG%Get_Randoms(2)
            zt = S5Inverse(r(1))
            If (r(2) * s5(zt) <= g(zt)) Exit  ! Accept zt
        End do
    End if
    Return

    Contains
        Function s1(zt) Result(s)
            Real(dp):: s
            Real(dp), Intent(In):: zt
            Real(dp), Parameter:: e = Exp(1.0_dp)
            s = e * zt**2
        End Function s1

        Function s2(zt) Result(s)
            Real(dp):: s
            Real(dp), Intent(In):: zt
            s = (2.0_dp - zt) * zt
        End Function s2

        Function s4(zt) Result(s)
            Real(dp):: s
            Real(dp), Intent(In):: zt
            s = 2.1809255541887791_dp - 1.02_dp * zt
        End Function s4

        Function s5(zt) Result(s)
            Real(dp):: s
            Real(dp), Intent(In):: zt
            Real(dp), Parameter:: a = 3.24_dp * Exp(2.24_dp)
            Real(dp), Parameter:: b = 112.0_dp / 45.0_dp
            s = a * Exp(-b * zt)
        End Function s5

        Function S1Inverse(r) Result(zt)
            Real(dp):: zt
            Real(dp), Intent(In):: r
            Real(dp), Parameter:: a = 48.0_dp / 125.0_dp
            Real(dp), Parameter:: b = 1.0_dp / 3.0_dp
            zt = a * r**b
        End Function S1Inverse

        Function S2Inverse(r) Result(zt)
            Real(dp):: zt
            Real(dp), Intent(In):: r
            Real(dp), Parameter:: a = 1.0_dp / 12.0_dp
            Real(dp), Parameter:: b = 0.203401_dp / 2.25_dp
            Real(dp), Parameter:: c = 0.66389_dp / 1.40625_dp
            zt = a + Sqrt(b + c * r)
        End Function S2Inverse

        Function S3Inverse(r) Result(zt)
            Real(dp):: zt
            Real(dp), Intent(In):: r
            Real(dp), Parameter:: a = 5.0_dp / 6.0_dp
            Real(dp), Parameter:: b = 49.0_dp / 150.0_dp
            zt = a + b * r
        End Function S3Inverse

        Function S4Inverse(r) Result(zt)
            Real(dp):: zt
            Real(dp), Intent(In):: r
            Real(dp), Parameter:: a = 5.7612645570087825_dp
            Real(dp), Parameter:: b = 1.3426416520733378_dp
            Real(dp), Parameter:: c = 3.4076721779841483_dp
            Real(dp), Parameter:: d = 0.00060876202977868390_dp
            Real(dp), Parameter:: e = 0.00053600480812538490_dp
            zt = (a + b * r) / (c + Sqrt(d - e * r))
        End Function S4Inverse

        Function S5Inverse(r) Result(zt)
            Real(dp):: zt
            Real(dp), Intent(In):: r
            Real(dp), Parameter:: a = 1.8_dp
            Real(dp), Parameter:: b = 45.0_dp / 112.0_dp
            zt = a - b * Log(1.0_dp - r)
            ! Note, most pseudo-random generators draw r in the interval 0 <= r < 1.
            ! For these, 1 - r is never zero, so the Log function does not fail.
        End Function S5Inverse

        Function g(zt)
            Use Kinds, Only: dp
            Implicit None
            Real(dp):: g                ! Standardized Boltzmann density function (peak is 1 at zt = 1)
            Real(dp), Intent(In):: zt   ! Nondimensionalized speed
            g = zt**2 * Exp(1.0_dp - zt**2)
        End Function g
End Function Sample_Nondimensional_Translational_Speed

End Module Target_Motion
