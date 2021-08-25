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
Module Diverge_exact

    Implicit None
    Integer, Parameter:: dp = Selected_Real_Kind(p=15)

    ! Mathematical Constants
    Integer, Parameter:: nMax = 10, nf = 6
    Real(dp), Parameter:: pi = 3.15159265359897932
    Real(dp), Parameter:: twoPi = 6.28318530717958650
    Real(dp):: BinomialCoefficient(0:nMax, 0:nMax)
    Real(dp):: Factorial(0:2*nMax), InverseFactorial(0:2*nMax)
    Real(dp):: Natural_Number(0:2*nMax)

    ! Units
    !   Distance [km]
    !   Time [s]
    !   Speed [km/s]
    !   Angle [radian]

    ! Physical Constants
    Real(dp), Parameter:: omegaEarth = 7.29212E-5_dp    ! radians / second
    Real(dp), Parameter:: muEarth = 3.986004418E+5_dp   ! km**3 / s**2
    Real(dp), Parameter:: neutronMass = 1.04541E-5_dp   ! keV / (km/s)**2

    ! Real(dp):: cp2(1:2), cp3(1:3), cp3a(1:3)   


    Interface Polynomial
        Module Procedure Polynomial_Scalar
        Module Procedure Polynomial_Vector
        Module Procedure Polynomial_Vector_Matrix
    End Interface Polynomial

    Interface Delta_Polynomial
        Module Procedure Delta_Polynomial_1
        Module Procedure Delta_Polynomial_2
    End Interface Delta_Polynomial

    Interface LegendreP
        Module Procedure LegendreP_Scalar   ! Inefficient: Wastes lower orders computed
        Module Procedure LegendreP_Vector   ! Efficient: Sequential but returns all orders computed
        Module Procedure LegendreP_Matrix   ! Most efficient: Vectorizes arithmentic for multiple values of x
    End Interface

    Interface deltaCos
        Module Procedure deltaCos2
        Module Procedure deltaCos3
    End Interface deltaCos

    Interface deltaSin
        Module Procedure deltaSin2
        Module Procedure deltaSin3
    End Interface deltaSin

    Contains

    ! This implementation of Divergence 
    !   assumes that the satellite is in a Circular Orbit

Subroutine Initialize_Diverge
    Integer:: i
    Do i = 1, 2 * nMax
        Natural_Number(i) = Real(i, dp)
    End do
    Call Initialize_Binomial_Coefficients(nMax, BinomialCoefficient)
    Call Initialize_Factorials(2*nMax, Factorial, InverseFactorial)
End Subroutine Initialize_Diverge

Function Effective_Detector_Radius(vRel, vRelXSat, vRelYSat, vRelZSat) Result(rEff)
    ! Currently, this is a dummy: rEff is set to a constant value.
    ! Hence, compiler warnings that input variables are not used can be ignored for now.
    ! Structure is here for future use.

    Real(dp), Intent(In):: vRel ! All speeds in [km/s]
        ! velocity components of neutron relative to detector (x,y,z right-handed)
    Real(dp), Intent(In), Optional:: vRelXSat 
        ! x in direction of motion of satellite in circular orbit
    Real(dp), Intent(In), Optional:: vRelYSat ! y in direction normal to orbital plane
    Real(dp), Intent(In), Optional:: vRelZSat ! z in radial direction
    Real(dp):: rEff             ! [km]
    ! For a prograde, equatorial, circular orbit: x is east; y is north; z is up.
    ! Positive vRelXSat implies that neutron is approaching the satellite from behind.
    ! Positive vRelZSat implies that neutron is approaching from below.
    ! This is a place-holder for a function
    !    that accounts for actual detector size and design.

    !   To account for cross section variation for spherical detector, use the next line:
    ! rEff = rDetector * fEff(2.0_dp * rDetector * sigmaTotal(KE(vRel))
    !   where KE converts neutron speed in km/s to neutron energy in units used by cross section routine
    !   and fEff is provided next
    !   Also, rDetector could be a parameter variable or an argument

    rEff = 1.0E-4_dp
End Function Effective_Detector_Radius

Function fEff(tau)
    Real(dp), Intent(In):: tau  ! optical thickness through diameter of sphere
                                !   tau = sigmaT(vRel) * Rd
    Real(dp):: fEff             ! ratio of effective radius to the actual radius
                                !   of detector sphere
    If (tau >= 0.01_dp) then
        fEff = Sqrt(1.0_dp - 2.0_dp * (1.0_dp - Exp(-tau) * (1.0_dp + tau)) / tau**2)
    Else
        fEff = Sqrt(tau / 6.0_dp) * (2.0_dp + tau * (-0.375_dp + tau * 0.06484375))
    End if
End Function fEff

    
Function Detection_Solid_Angle_CM &
    & (r0Vec, v0Vec, rIntVec, vIntVec, &    ! intercept solution
    & vSat, AxisHatSat, &                   ! satellite circular orbit
    & Emit, &                               ! .True. for emission of neutron; .False. for scattering
    & vEmitterVec, &                        ! Optional: must be present if Emit = .True.
    & vPrimeVec, An, QAbsorbed) &           ! Optional: must be present if Emit = .False.
    & Result(DeltaOmegaCM)

! NOTE: returns DeltaOmegaCM = 0.0_dp if scatter varaibles and solution variables
!       yield inconsistent values for vCM
!       Calling routine should test for this, 
!       then continue iterating for more precise intercept solution
!           (or code error should be found)

! Emission is implemented by setting some variables otherwise obtained from scattering variables
!   to have values consistent with emission.
!   vCM is held constant (fixed emission energy in emitter rest frame) 
!       as direction of emission varies.

! Intercept solution
Real(dp), Intent(In):: r0Vec(1:3)
Real(dp), Intent(In):: v0Vec(1:3)
Real(dp), Intent(In):: rIntVec(1:3)
Real(dp), Intent(In):: vIntVec(1:3)

! Satellite circular orbit parameters
Real(dp), Intent(In):: vSat
Real(dp), Intent(In):: AxisHatSat(1:3)  
    ! For a prograde equatorial orbit: 
    ! = Earth's axis = ZHat(1:3) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
    ! In general, equals rVec cross vVec for satellite at any point in its orbit

! Process Indicator
Logical, Intent(In):: Emit

! Emission variable MUST BE PRESENT IF EMIT = .TRUE.
Real(dp), Intent(In), Optional:: vEmitterVec(1:3)

! Scattering variables  MUST BE PRESENT IF EMIT = .FALSE.
Real(dp), Intent(In), Optional:: vPrimeVec(1:3)
Real(dp), Intent(In), Optional:: An
Real(dp), Intent(In), Optional:: QAbsorbed

! Result variable
Real(dp):: DeltaOmegaCM

! Flags
Logical:: Near_Vertical
Logical:: Inbound, Outbound
Logical:: Direct    ! ascertained by vertical component of arrival velocity
Logical:: Hyperbolic, Elliptical, Near_Parabolic
Logical:: Intercept_Error
Logical:: Collinear

! Tolerances
Real(dp), Parameter:: Near_Vertical_Tol = 0.001_dp    ! about 42 miles at GEO 
! *** Do the others work to this tolerance?
Real(dp), Parameter:: Near_Parabolic_Tol = 0.0001_dp
! *** where does accuracy of hyperbolic and elliptical formulas break down?
Real(dp), Parameter:: Intercept_RelTol = 0.00001_dp
! *** What does Whit expect?
Real(dp), Parameter:: s_RelTol = 1.0E-8_dp
Real(dp), Parameter:: Collinear_Tol = 0.01_dp
Real(dp), Parameter:: Use_ZHatCM_Tol = 0.1_dp
Real(dp), Parameter:: dxi_RelTol = 1.0E-8_dp
! See Real(dp), Parameter:: f_RelTol = 1.0E-8_dp in Subroutine Find_xi 

! Reference directions
Real(dp), Parameter:: XHat(1:3) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
Real(dp), Parameter:: YHat(1:3) = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
Real(dp), Parameter:: ZHat(1:3) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)

! Neutron quantities
Real(dp):: r0
Real(dp):: vAirVec(1:3)
Real(dp):: uVec(1:3)
Real(dp):: vPrimeCMVec(1:3)
Real(dp):: vPrimeCM
Real(dp):: OmegaHatCMPrime(1:3), OmegaHatCM(1:3)
Real(dp):: vCM
Real(dp):: rInt, kappaInt, vInt
Real(dp):: v0
Real(dp):: Omega0Hat(1:3)
Real(dp):: vCMVec(1:3)
Real(dp):: mu0CM, omega0CM
Real(dp):: xHatCM(1:3), yHatCM(1:3), zHatCM(1:3)
Real(dp):: x0Hat(1:3), y0Hat(1:3), z0Hat(1:3)
Real(dp):: xHatSat(1:3), yHatSat(1:3), zHatSat(1:3)
Real(dp):: zeta0
Real(dp):: v0x, v0y, v0z
Real(dp):: kappa0
Real(dp):: tau
Real(dp):: h
Real(dp):: p
Real(dp):: delta_e, e
Real(dp):: rP, rApogee
Real(dp):: vSatVec(1:3)
Real(dp):: vRelVec(1:3)
Real(dp):: phi0

! Derivatives
    ! Name convention: the partial derivative of y with respect to x is "dy_dx"
Real(dp):: ChiInt(1:2)

    ! Non-Collinear Derivatives
Real(dp):: dOmegaHatCM_dmu0CM(1:3), dOmegaHatCM_domega0CM(1:3)
Real(dp):: dkappa0_dmu0CM, dkappa0_domega0CM
Real(dp):: dv0_dmu0CM, dv0_domega0CM, dv0Vec_dmu0CM(1:3), dv0Vec_domega0CM(1:3)
Real(dp):: dv0x_dmu0CM, dv0x_domega0CM, dv0y_dmu0CM, dv0y_domega0CM
Real(dp):: dzeta0_dmu0CM, dzeta0_domega0CM
Real(dp):: dOmega0Hat_dmu0CM(1:3), dOmega0Hat_domega0CM(1:3)
Real(dp):: dphi0_dmu0CM, dphi0_domega0CM
Real(dp):: dp_dmu0CM, dp_domega0CM
Real(dp):: dh_dmu0CM, dh_domega0CM
Real(dp):: de_dmu0CM, de_domega0CM
Real(dp):: drPerigee_dmu0CM, drPerigee_domega0CM
Real(dp):: drApogee_dmu0CM, drApogee_domega0CM
Real(dp):: dvCMVec_dmu0CM(1:3), dvCMVec_domega0CM(1:3)

    ! Collinear Derivatives
Real(dp):: nuX, nuY
Real(dp):: dOmegaHatCM_dnuX(1:3), dOmegaHatCM_dnuY(1:3)
Real(dp):: dkappa0_dnuX, dkappa0_dnuY
Real(dp):: dv0_dnuX, dv0_dnuY, dv0x_dnuX, dv0x_dnuY, dv0y_dnuX, dv0y_dnuY
Real(dp):: dv0Vec_dnuX(1:3), dv0Vec_dnuY(1:3)
Real(dp):: dzeta0_dnuX, dzeta0_dnuY
Real(dp):: dOmega0Hat_dnuX(1:3), dOmega0Hat_dnuY(1:3)
Real(dp):: dphi0_dnuX, dphi0_dnuY
Real(dp):: dp_dnuX, dp_dnuY
Real(dp):: dh_dnuX, dh_dnuY
Real(dp):: de_dnuX, de_dnuY
Real(dp):: drPerigee_dnuX, drPerigee_dnuY
Real(dp):: drApogee_dnuX, drApogee_dnuY
Real(dp):: dvCMVec_dnuX(1:3), dvCMVec_dnuY(1:3)

! True Anomalies and derivatives
Real(dp):: theta0, thetaInt
Real(dp):: dtheta0_dmu0CM, dtheta0_domega0CM, dtheta0_dnuX, dtheta0_dnuY
Real(dp):: dthetaInt_dmu0CM, dthetaInt_domega0CM, dthetaInt_dnuX, dthetaInt_dnuY

! Time of Flight and derivatives
Real(dp):: DeltaTheta
Real(dp):: dDeltaTheta_dmu0CM, dDeltaTheta_domega0CM
Real(dp):: dDeltaTheta_dnuX, dDeltaTheta_dnuY

    ! Flight Times
Real(dp):: t0, dt0_dmu0CM, dt0_domega0CM, dt0_dnuX, dt0_dnuY
Real(dp):: tInt, dtInt_dmu0CM, dtInt_domega0CM, dtInt_dnuX, dtInt_dnuY
Real(dp):: DeltaT, dDeltaT_dmu0CM, dDeltaT_domega0CM, dDeltaT_dnuX, dDeltaT_dnuY

    ! Near_Parabolic
! Real(dp):: xi0, xiInt (not used)

    ! Hyperbolic & Elliptical
Real(dp):: a, da_dmu0CM, da_dOmega0CM, da_dnuX, da_dnuY
Real(dp):: E0, dE0_dmu0CM, dE0_dOmega0CM, dE0_dnuX, dE0_dnuY
Real(dp):: EInt, dEInt_dmu0CM, dEInt_dOmega0CM, dEInt_dnuX, dEInt_dnuY

    ! Elliptical
Real(dp):: Period, dPeriod_dmu0CM, dPeriod_domega0CM, dPeriod_dnuX, dPeriod_dnuY

! First Derivative Method
Real(dp):: draVec_dmu0CM(1:3), draVec_domega0CM(1:3)
Real(dp):: draVec_dnuX(1:3), draVec_dnuY(1:3)
Real(dp):: ScaleFactor, ScaleFactorVec(1:3)
Real(dp):: OmegaHatPerpendicular(1:3), OmegaHatRel(1:3)
Real(dp):: cos_IncidenceAngle               ! Determines whether finite difference method is used
Real(dp):: cos_IncidenceAngle_Tol = 0.2_dp  ! *** Needs tuning for balance of efficiency and accuracy
Real(dp):: Interaction_Area
Real(dp):: vRel, rEff
Real(dp):: cp(1:3)  ! holds result of a cross product (whence "cp")
                    ! that is used within the next few lines

! EXECUTION STARTS HERE

! Intercept Type
Direct = (Dot_Product(v0Vec, r0Vec) >= 0.0_dp)
Outbound = (Dot_Product(vIntVec, rIntVec) >= 0.0_dp)
Inbound = .not. Inbound

! Compute Neutron Parameters
r0 = Norm(r0Vec)
vAirVec = omegaEarth * (/ -r0Vec(2), r0Vec(1), 0.0_dp /)
rInt = Norm(rIntVec)
vInt = Norm(vIntVec)
kappaInt = vInt**2 / 2.0_dp
v0 = Norm(v0Vec)
Omega0Hat = v0Vec / v0
If (emit) then
    uVec = vEmitterVec
    vCMVec = v0Vec - uVec
    vCM = Norm(vcmVec)
    OmegaHatCMPrime = Hat(vCMVec) ! Used to compute mu0CM and set coordinate basis directions
Else
    uVec = (An * vAirVec + vPrimeVec) / (An + 1.0_dp)
    vPrimeCMVec = vPrimeVec - uVec
    vPrimeCM = Norm(vPrimeCMVec)
    OmegaHatCMPrime = vPrimeCMVec / vPrimeCM
    vCM = Sqrt(vPrimeCM**2 - (An / (An + 1.0_dp)) * 2.0_dp * QAbsorbed / neutronMass)
    ! Check consistency of scatter and intercept solution
    vCMVec = v0Vec - uVec
    Intercept_Error = (Abs(vCM - Norm(vCMVec)) > vCM * Intercept_RelTol)
    If (.not. Intercept_Error) then
        vCM = Norm(vcmVec)  ! for scatter, resets to be exactly consistent with intercept solution as given
    Else
        DeltaOmegaCM = 0.0_dp
        RETURN
    End if
End if
OmegaHatCM = vCMVec / vCM
mu0CM = Dot_Product(OmegaHatCM, OmegaHatCMPrime)
Collinear = ((1-Abs(mu0CM)) <= Collinear_Tol)
z0Hat = Hat(r0Vec)
zeta0 = Dot_Product(Omega0Hat,z0Hat)
Near_Vertical = ((1.0_dp - Abs(zeta0)) <= Near_Vertical_Tol)
If (Near_Vertical) then
    y0Hat = Hat(Cross_Product(vSatVec,z0Hat))
Else
    y0Hat = Hat(Cross_Product(z0Hat, Omega0Hat))
End if
x0Hat = Cross_Product(y0Hat, z0Hat)
v0x = Dot_Product(v0Vec, x0Hat)
v0y = Dot_Product(v0Vec, y0Hat)
v0z = Dot_Product(v0Vec, z0Hat)
kappa0 = v0**2 / 2.0_dp
tau = kappa0 - muEarth / r0
h = r0 * v0 * Sqrt(1.0_dp - zeta0**2)
p = h**2 / muEarth
delta_e = (2.0_dp * tau * h**2 / muEarth**2) &
    & / (1.0_dp + Sqrt(1.0_dp + 2.0_dp * tau * h**2 / muEarth**2))
Hyperbolic = delta_e > Near_Parabolic_Tol
Elliptical = delta_e < -Near_Parabolic_Tol
Near_Parabolic = (Abs(delta_E) <= Near_Parabolic_Tol)
e = 1 + delta_e
rP = p / (1.0_dp + e)
If (delta_e < 0) then
    rApogee = p / (-delta_e)
Else
    rApogee = Huge(rApogee)
End if

! Satellite Parameters
zHatSat = Hat(rIntVec)
xHatSat = Hat(Cross_Product(AxisHatSat, zHatSat))
yHatSat = Cross_Product(zHatSat, xHatSat)
vSatVec = vSat * xHatSat                        ! assumes circular orbit
vRelVec = vIntVec - vSatVec

! Derivatives of Neutron Parameters
    ! Non-collinear
If (.not. collinear) then
    zHatCM = OmegaHatCMPrime
    yHatCM = Hat(Cross_Product(zHatCM, OmegaHatCM))
    xHatCM = Cross_product(yHatCM, zHatCM)
    mu0cm = Dot_Product(OmegaHatCM, zHatCM)
    omega0CM = atan2(Dot_Product(OmegaHatCM, yHatCM), &
        & Dot_Product(OmegaHatCM, xHatCM))
    chiInt = (/ mu0cm, omega0cm /)
    dOmegaHatCM_dmu0CM = zHatCM - (mu0CM / Sqrt(1.0_dp - mu0CM**2)) &
        & * (cos(omega0CM) * xHatCM + sin(omega0CM) * yHatCM)
    dOmegaHatCM_domega0CM = Sqrt(1.0_dp - mu0CM**2) &
        & * (cos(omega0CM) * yHatCM - sin(omega0CM) * xHatCM)
    dkappa0_dmu0CM = vCM * Dot_Product(uVec, dOmegaHatCM_dmu0CM)
    dkappa0_domega0CM = 0.0_dp
    dv0_dmu0CM = dkappa0_dmu0cm / v0
    dv0_domega0CM = dkappa0_domega0cm / v0
    dzeta0_dmu0CM = (vCM / v0) &
        & * Dot_Product((z0Hat - (zeta0 / v0) * uVec), dOmegaHatCM_dmu0CM)
    dzeta0_domega0CM = (vCM / v0) &
        & * Dot_Product((z0Hat - (zeta0 / v0) * uVec), dOmegaHatCM_domega0CM)
    dOmega0Hat_dmu0CM = (vCM / v0) * dOmegaHatCM_dmu0CM &
        & - (dv0_dmu0CM / v0) * Omega0Hat
    dOmega0Hat_domega0CM = (vCM / v0) * dOmegaHatCM_domega0CM &
        & - (dv0_domega0CM / v0) * Omega0Hat
    dv0Vec_dmu0CM = dv0_dmu0CM * Omega0Hat + v0 * dOmega0Hat_dmu0CM
    dv0Vec_domega0CM = dv0_domega0CM * Omega0Hat + v0 * dOmega0Hat_domega0CM
    dv0x_dmu0CM = Dot_Product(dv0Vec_dmu0CM, x0Hat)
    dv0x_domega0CM = Dot_Product(dv0Vec_domega0CM, x0Hat)
    dv0y_dmu0CM = Dot_Product(dv0Vec_dmu0CM, y0Hat)
    dv0y_domega0CM = Dot_Product(dv0Vec_domega0CM, y0Hat)
    ! z components of derivatives of v0vec are not needed
    cp = Cross_Product(z0Hat, Omega0Hat)
    dphi0_dmu0CM = Dot_Product(cp, dOmega0Hat_dmu0CM) &
        & / (1.0_dp - zeta0**2)
    dphi0_domega0CM = Dot_Product(cp, dOmega0Hat_domega0CM) &
        & / (1.0_dp - zeta0**2)
    dp_dmu0CM = p * (dkappa0_dmu0CM / kappa0 &
        & - 2.0_dp * zeta0 * dzeta0_dmu0CM / (1.0_dp - zeta0**2))
    dp_domega0CM = p * (dkappa0_domega0CM / kappa0 &
        & - 2.0_dp * zeta0 * dzeta0_domega0CM / (1.0_dp - zeta0**2))
    dh_dmu0CM = muEarth * dp_dmu0CM / (2.0_dp * h)  ! fails if Near-Vertical
    dh_domega0CM = muEarth * dp_domega0CM / (2.0_dp * h)  ! fails if Near-Vertical
    de_dmu0CM = (tau * dp_dmu0CM + p * dkappa0_dmu0cm) / (e * muEarth)
    de_domega0CM = (tau * dp_domega0CM + p * dkappa0_domega0cm) / (e * muEarth)
    drPerigee_dmu0CM = rP * (dp_dmu0CM / p - de_dmu0CM / (1.0_dp + e))
    drPerigee_domega0CM = rP * (dp_domega0CM / p - de_domega0CM / (1.0_dp + e))
    If (Elliptical) then 
        drApogee_dmu0CM = rApogee * (dp_dmu0CM / p + de_dmu0CM / (-delta_e))
        drApogee_domega0CM = rApogee * (dp_domega0CM / p + de_domega0CM / (-delta_e))
    Else
        drApogee_dmu0CM = Huge(drApogee_dmu0CM)
        drApogee_domega0CM = Huge(drApogee_domega0CM)
    End if
     dvCMVec_dmu0CM = vCM * dOmegaHatCM_dmu0CM
     dvCMVec_domega0CM = vCM * dOmegaHatCM_domega0CM
Else
    ! Collinear
    zHatCM = OmegaHatCMPrime
    yHatCM = Cross_Product(ZHat, zHatCM)  ! intermediate quantity
    If (Norm(yHatCM) >= Use_ZHatCM_Tol) then
        yHatCM = Hat(yHatCM)
    Else ! Use XHat
        yHatCM = Hat(Cross_Product(XHat, zHatCM))
    End if
    xHatCM = Cross_Product(yHatcm, zHatCM)
    nuX = Dot_Product(OmegaHatCM, xHatCM)
    nuY = Dot_Product(OmegaHatCM, yHatCM)
    chiInt = (/ nuX, nuY /)
    mu0CM = Sqrt(1.0_dp - nuX**2 - nuY**2)
    dOmegaHatCM_dnuX = xHatCM - (nuX / mu0CM) * zHatCM
    dOmegaHatCM_dnuY = yHatCM - (nuY / mu0CM) * zHatCM
    dkappa0_dnuX = vCM * Dot_Product(uVec, dOmegaHatCM_dnuX)
    dkappa0_dnuY = vCM * Dot_Product(uVec, dOmegaHatCM_dnuY)
    dv0_dnuX = dKappa0_dnuX / v0
    dv0_dnuY = dKappa0_dnuY / v0
    dOmega0Hat_dnuX = (vCM / v0) * dOmegaHatCM_dnuX &
        & - (dv0_dnuX / v0) * Omega0Hat
    dOmega0Hat_dnuY = (vCM / v0) * dOmegaHatCM_dnuY &
        & - (dv0_dnuY / v0) * Omega0Hat
    dv0Vec_dnuX = dv0_dnuX * Omega0Hat + v0 * dOmega0Hat_dnuX
    dv0Vec_dnuY = dv0_dnuY * Omega0Hat + v0 * dOmega0Hat_dnuY
    dv0x_dnuX = Dot_Product(dv0Vec_dnuX, x0Hat)
    dv0x_dnuY = Dot_Product(dv0Vec_dnuY, x0Hat)
    dv0y_dnuX = Dot_Product(dv0Vec_dnuX, y0Hat)
    dv0y_dnuY = Dot_Product(dv0Vec_dnuY, y0Hat)
    ! z components of derivatives of v0vec are not needed
    dzeta0_dnuX = Dot_Product(z0Hat, dOmega0Hat_dnuX)
    dzeta0_dnuY = Dot_Product(z0Hat, dOmega0Hat_dnuY)
    cp = Cross_Product(z0Hat, Omega0Hat)  ! intermediate quantity
    dphi0_dnuX = Dot_Product(cp,dOmega0Hat_dnuX) &
        & / (1.0_dp - zeta0**2)
    dphi0_dnuY = Dot_Product(cp,dOmega0Hat_dnuY) &
        & / (1.0_dp - zeta0**2)
    dp_dnuX = (p / h) * (dkappa0_dnuX / kappa0 &
        & - (2.0_dp * zeta0 / (1.0_dp - zeta0**2)) * dzeta0_dnuX)
    dp_dnuY = (p / h) * (dkappa0_dnuY / kappa0 &
        & - (2.0_dp * zeta0 / (1.0_dp - zeta0**2)) * dzeta0_dnuY)
    dh_dnuX = muEarth * dp_dnuX / (2.0_dp * h)
    dh_dnuY = muEarth * dp_dnuY / (2.0_dp * h)
    de_dnuX = (tau * dp_dnuX + p * dkappa0_dnuX) / (e * muEarth)
    de_dnuY = (tau * dp_dnuY + p * dkappa0_dnuY) / (e * muEarth)
    drPerigee_dnuX = rP * (dp_dnuX / p - de_dnuX / (1.0_dp + e))
    drPerigee_dnuY = rP * (dp_dnuY / p - de_dnuY / (1.0_dp + e))
    If (delta_e < 0.0_dp) then
        drApogee_dnuX = rApogee * (dp_dnuX / p - de_dnuX / (-delta_e))
        drApogee_dnuY = rApogee * (dp_dnuY / p - de_dnuY / (-delta_e))
    Else
        drApogee_dnuX = Huge(drApogee_dnuX)
        drApogee_dnuY = Huge(drApogee_dnuY)
    End if
    dvCMVec_dnuX = vCM * dOmegaHatCM_dnuX
    dvCMVec_dnuY = vCM * dOmegaHatCM_dnuY
End if

! DeltaT and related variables, and their derivatives
If ( .not. Near_Vertical) then
    ! True Anomalies
    theta0 = p / (e * r0) + delta_e / e ! intermediate quantity
    theta0 = ASin(Sqrt((2.0_dp - theta0) * theta0))
    thetaInt = p / (e * rInt) + delta_e / e ! intermediate quantity
    thetaInt = ASin(Sqrt((2.0_dp - thetaInt) * thetaInt))
    dtheta0_dmu0CM = (cos(theta0) * de_dmu0CM - dp_dmu0CM / r0) &
        & / (e * sin(theta0))
    dtheta0_domega0CM = (cos(theta0) * de_domega0CM - dp_domega0CM / r0) &
        & / (e * sin(theta0))
    dthetaInt_dmu0CM = (cos(thetaInt) * de_dmu0CM - dp_dmu0CM / rInt) &
        & / (e * sin(thetaInt))
    dthetaInt_domega0CM = (cos(thetaInt) * de_domega0CM - dp_domega0CM / rInt) &
        & / (e * sin(thetaInt))
    If (Outbound) then
        If (Direct) then
            DeltaTheta = thetaInt - theta0
            dDeltaTheta_dmu0CM = dthetaInt_dmu0CM - dtheta0_dmu0CM
            dDeltaTheta_domega0CM = dthetaInt_domega0CM - dtheta0_domega0CM
        Else ! Indirect
            DeltaTheta = thetaInt + theta0
            dDeltaTheta_dmu0CM = dthetaInt_dmu0CM + dtheta0_dmu0CM
            dDeltaTheta_domega0CM = dthetaInt_domega0CM + dtheta0_domega0CM
        End if
    Else ! Inbound
        If (Direct) then
            DeltaTheta = twoPi - thetaInt - theta0
            dDeltaTheta_dmu0CM = -dthetaInt_dmu0CM - dtheta0_dmu0CM
            dDeltaTheta_domega0CM = -dthetaInt_domega0CM - dtheta0_domega0CM
        Else ! Indirect
            DeltaTheta = twoPi - thetaInt + theta0
            dDeltaTheta_dmu0CM = -dthetaInt_dmu0CM + dtheta0_dmu0CM
            dDeltaTheta_domega0CM = -dthetaInt_domega0CM + dtheta0_domega0CM
        End if
    End if
    ! Time of Flight
    If (Hyperbolic) Call Hyperbolic_DeltaT
    If (Elliptical) Call Elliptical_DeltaT
    If (Near_Parabolic) Call Near_Parabolic_DeltaT
Else
    If (Hyperbolic) then
        CALL Hyperbolic_Near_Vertical_DeltaT
    Else if (Elliptical) then
        CALL Elliptical_Near_Vertical_DeltaT
    Else if (Near_Parabolic) then
        CALL Parabolic_Near_Vertical_DeltaT
    End if
End if

! Detector Interaction Solid Angle by First Derivatives
If (.not. Near_Vertical) then
    If (.not. Collinear) then
        phi0 = atan2(Dot_Product(Omega0Hat,y0Hat), Dot_Product(Omega0Hat,x0Hat))
        draVec_dmu0CM = rInt * (((cos(DeltaTheta) * &
            & (cos(phi0) * x0Hat + sin(phi0) * y0Hat) &
            & - sin(DeltaTheta)) * z0Hat) * dDeltaTheta_dmu0CM &
            & + sin(DeltaTheta) * (-sin(phi0) * x0Hat &
            & + cos(phi0) * y0Hat) * dphi0_dmu0CM)
        draVec_domega0CM = rInt * (((cos(DeltaTheta) * &
            & (cos(phi0) * x0Hat + sin(phi0) * y0Hat) &
            & - sin(DeltaTheta)) * z0Hat) * dDeltaTheta_domega0CM &
            & + sin(DeltaTheta) * (-sin(phi0) * x0Hat &
            & + cos(phi0) * y0Hat) * dphi0_domega0CM)
    Else ! Collinear
        draVec_dnuX = rInt * (((cos(DeltaTheta) * &
            & (cos(phi0) * x0Hat + sin(phi0) * y0Hat) &
            & - sin(DeltaTheta)) * z0Hat) * dDeltaTheta_dnuX &
            & + sin(DeltaTheta) * (-sin(phi0) * x0Hat &
            & + cos(phi0) * y0Hat) * dphi0_dnuX)
        draVec_dnuY = rInt * (((cos(DeltaTheta) * &
            & (cos(phi0) * x0Hat + sin(phi0) * y0Hat) &
            & - sin(DeltaTheta)) * z0Hat) * dDeltaTheta_dnuY &
            & + sin(DeltaTheta) * (-sin(phi0) * x0Hat &
            & + cos(phi0) * y0Hat) * dphi0_dnuY)
    End if
Else ! Near Vertical
    If (.not. Collinear) then
        draVec_dmu0CM = (DeltaT * dv0x_dmu0CM + v0x * dDeltaT_dmu0CM) * x0Hat &
            & + (DeltaT * dv0y_dmu0CM + v0y * dDeltaT_dmu0CM) * y0Hat
        draVec_domega0CM = (DeltaT * dv0x_domega0CM + v0x * dDeltaT_domega0CM) * x0Hat &
            & + (DeltaT * dv0y_domega0CM + v0y * dDeltaT_domega0CM) * y0Hat
    Else ! Collinear
        draVec_dnuX = (DeltaT * dv0x_dnuX + v0x * dDeltaT_dnuX) * x0Hat &
            & + (DeltaT * dv0y_dnuX + v0y * dDeltaT_dnuX) * y0Hat
        draVec_dnuY = (DeltaT * dv0x_dnuY + v0x * dDeltaT_dnuY) * x0Hat &
            & + (DeltaT * dv0y_dnuY + v0y * dDeltaT_dnuY) * y0Hat
    End if
End if
If (.not. collinear) then
    draVec_dmu0CM = draVec_dmu0CM - vIntVec * dDeltaT_dmu0CM
    draVec_domega0CM = draVec_domega0CM - vIntVec * dDeltaT_domega0CM
    ScaleFactorVec = Cross_Product(draVec_dmu0CM, draVec_domega0CM)
    ScaleFactor = Norm(ScaleFactorVec)
Else ! collinear
    draVec_dnuX = draVec_dnuX - vIntVec * dDeltaT_dnuX
    draVec_dnuY = draVec_dnuY - vIntVec * dDeltaT_dnuY
    ScaleFactorVec = Cross_Product(draVec_dnuX, draVec_dnuY)
    ScaleFactor = Norm(ScaleFactorVec) * Sqrt(1.0_dp - nuX**2 - nuY**2) 
End if
OmegaHatPerpendicular = Hat(ScaleFactorVec)
OmegaHatRel = Hat(vRelVec)
vRel = Norm(vRelVec)
cos_IncidenceAngle = Abs(Dot_Product(OmegaHatPerpendicular, OmegaHatRel))

! A kludge to use until finite difference method is ready
If (cos_IncidenceAngle < cos_IncidenceAngle_Tol) Then
    cos_IncidenceAngle = cos_IncidenceAngle_Tol
End if

rEff = Effective_Detector_Radius(vRel)
Interaction_Area = pi * rEff**2 / cos_IncidenceAngle
DeltaOmegaCM = Interaction_Area / ScaleFactor

RETURN

Contains

    Subroutine Hyperbolic_DeltaT()
        If (Inbound) then
            Print *, "Error:: Hyperbolic_DeltaT : Inbound Intercept"
            Stop
        End if
        ! Quantities
        a = p / (delta_e * (1.0_dp + e))
        E0 = acosh((a + r0) / (a * e))
        EInt = acosh((a + rInt) / (a * e))
        t0 = a * Sqrt(a/muEarth) * (e * sinh(E0) - E0)
        tInt = a * Sqrt(a/muEarth) * (e * sinh(EInt) - EInt)
        ! Derivatives
        If (.not. collinear) then
            da_dmu0CM = (a / p) * (dp_dmu0CM - 2.0_dp * a * e * de_dmu0CM)
            da_domega0CM = (a / p) * (dp_domega0CM - 2.0_dp * a * e * de_domega0CM)
            dE0_dmu0CM = -((r0 / a) * da_dmu0CM + ((a + r0) / e) * de_dmu0CM)
            dE0_domega0CM = -((r0 / a) * da_domega0CM &
                & + ((a + r0) / e) * de_domega0CM)
            dEInt_dmu0CM = -((rInt / a) * da_dmu0CM + ((a + rInt) / e) * de_dmu0CM)
            dEInt_domega0CM = -((rInt / a) * da_domega0CM &
                & + ((a + rInt) / e) * de_domega0CM)
            dt0_dmu0CM = 1.5_dp * t0 * da_dmu0CM / a &
                & + (a * Sqrt(a / muEarth)) / (1.0_dp + e * cosh(E0)) &
                & * (sinh(E0) * de_dmu0CM  &
                & + (sinh(E0)**2 + delta_e * (1.0_dp + e) * cosh(E0)**2) &
                & * dE0_dmu0CM)
            dt0_domega0CM = 1.5_dp * t0 * da_domega0CM / a &
                & + (a * Sqrt(a / muEarth)) / (1.0_dp + e * cosh(E0)) &
                & * (sinh(E0) * de_domega0CM  &
                & + (sinh(E0)**2 + delta_e * (1.0_dp + e) * cosh(E0)**2) &
                & * dE0_domega0CM)
            dtInt_dmu0CM = 1.5_dp * tInt * da_dmu0CM / a &
                & + (a * Sqrt(a / muEarth)) / (1.0_dp + e * cosh(EInt)) &
                & * (sinh(EInt) * de_dmu0CM  &
                & + (sinh(EInt)**2 + delta_e * (1.0_dp + e) * cosh(EInt)**2) &
                & * dEInt_dmu0CM)
            dtInt_domega0CM = 1.5_dp * tInt * da_domega0CM / a &
                & + (a * Sqrt(a / muEarth)) / (1.0_dp + e * cosh(EInt)) &
                & * (sinh(EInt) * de_domega0CM  &
                & + (sinh(EInt)**2 + delta_e * (1.0_dp + e) * cosh(EInt)**2) &
                & * dEInt_domega0CM)
            ! Combine results depending on intercept class
            !   hyperbolic must be outbound
            If (Direct) then
                DeltaT = tInt - t0
                dDeltaT_dmu0CM = dtInt_dmu0CM - dt0_dmu0CM
                dDeltaT_domega0CM = dtInt_domega0CM - dt0_domega0CM
            Else
                DeltaT = tInt + t0
                dDeltaT_dmu0CM = dtInt_dmu0CM + dt0_dmu0CM
                dDeltaT_domega0CM = dtInt_domega0CM + dt0_domega0CM
            End if
        Else ! collinear
            da_dnuX = (a / p) * (dp_dnuX - 2.0_dp * a * e * de_dnuX)
            da_dnuY = (a / p) * (dp_dnuY - 2.0_dp * a * e * de_dnuY)
            dE0_dnuX = -((r0 / a) * da_dnuX + ((a + r0) / e) * de_dnuX)
            dE0_dnuY = -((r0 / a) * da_dnuY + ((a + r0) / e) * de_dnuY)
            dEInt_dnuX = -((rInt / a) * da_dnuX + ((a + rInt) / e) * de_dnuX)
            dEInt_dnuY = -((rInt / a) * da_dnuY + ((a + rInt) / e) * de_dnuY)
            dt0_dnuX = 1.5_dp * t0 * da_dnuX / a &
                & + (a * Sqrt(a / muEarth)) / (1.0_dp + e * cosh(E0)) &
                & * (sinh(E0) * de_dnuX  &
                & + (sinh(E0)**2 + delta_e * (1.0_dp + e) * cosh(E0)**2) * dE0_dnuX)
            dt0_dnuY = 1.5_dp * t0 * da_dnuY / a &
                & + (a * Sqrt(a / muEarth)) / (1.0_dp + e * cosh(E0)) &
                & * (sinh(E0) * de_dnuY  &
                & + (sinh(E0)**2 + delta_e * (1.0_dp + e) * cosh(E0)**2) * dE0_dnuY)
            dtInt_dnuX = 1.5_dp * tInt * da_dnuX / a &
                & + (a * Sqrt(a / muEarth)) / (1.0_dp + e * cosh(EInt)) &
                & * (sinh(EInt) * de_dnuX  &
                & + (sinh(EInt)**2 + delta_e * (1.0_dp + e) * cosh(EInt)**2)&
                & * dEInt_dnuX)
            dtInt_dnuY = 1.5_dp * tInt * da_dnuY / a &
                & + (a * Sqrt(a / muEarth)) / (1.0_dp + e * cosh(EInt)) &
                & * (sinh(EInt) * de_dnuY  &
                & + (sinh(EInt)**2 + delta_e * (1.0_dp + e) * cosh(EInt)**2) &
                & * dEInt_dnuY)
            ! Combine results depending on intercept class
            !   hyperbolic must be outbound
            If (Direct) then
                DeltaT = tInt - t0
                dDeltaT_dnuX = dtInt_dnuX - dt0_dnuX
                dDeltaT_dnuY = dtInt_dnuY - dt0_dnuY
            Else
                DeltaT = tInt + t0
                dDeltaT_dnuX = dtInt_dnuX + dt0_dnuX
                dDeltaT_dnuY = dtInt_dnuY + dt0_dnuY
            End if
        End if 
    End Subroutine Hyperbolic_DeltaT
  
    Subroutine Elliptical_DeltaT
        ! Quantities
        a = p / (-delta_e * (1.0_dp + e))
        E0 = acos((a - r0) / (a * e))
        EInt = acos((a - rInt) / (a * e))
        t0 = a * Sqrt(a/muEarth) * (e * sin(E0) - E0)
        tInt = a * Sqrt(a/muEarth) * (e * sin(EInt) - EInt)
        ! Derivatives
        If (.not. collinear) then
            da_dmu0CM = (dp_dmu0CM + 2.0_dp * e * a * de_dmu0CM) &
                & / (-delta_e * (1.0_dp +e))
            da_domega0CM = (dp_domega0CM + 2.0_dp * e * a * de_domega0CM) &
                & / (-delta_e * (1.0_dp +e))
            dE0_dmu0CM = ((r0/a**2) * da_dmu0CM - cos(theta0) * de_dmu0CM) / e
            dE0_domega0CM = ((r0/a**2) * da_domega0CM - cos(theta0) * de_domega0CM) / e
            dEInt_dmu0CM = ((rInt/a**2) * da_dmu0CM - cos(thetaInt) * de_dmu0CM) / e
            dEInt_domega0CM = ((rInt/a**2) * da_domega0CM - cos(thetaInt) * de_domega0CM) / e
            dt0_dmu0CM = 1.5 * t0 * da_dmu0CM / a + a * Sqrt(a/muEarth) &
                & * (r0 * dE0_dmu0CM / a - sin(E0) * de_dmu0CM)
            dt0_domega0CM = 1.5 * t0 * da_domega0CM / a + a * Sqrt(a/muEarth) &
                & * (r0 * dE0_domega0CM / a - sin(E0) * de_domega0CM)
            dtInt_dmu0CM = 1.5 * tInt * da_dmu0CM / a + a * Sqrt(a/muEarth) &
                & * (rInt * dEInt_dmu0CM / a - sin(EInt) * de_dmu0CM)
            dtInt_domega0CM = 1.5 * tInt * da_domega0CM / a + a * Sqrt(a/muEarth) &
                & * (rInt * dEInt_domega0CM / a - sin(EInt) * de_domega0CM)
            ! Combine to get time of flight and its derivatives
            If (Outbound) then
                If (Direct) then
                    DeltaT = tInt - t0
                    dDeltaT_dmu0CM = dtInt_dmu0CM - dt0_dmu0CM
                    dDeltaT_domega0CM = dtInt_domega0CM - dt0_domega0CM
                Else ! Indirect
                    DeltaT = tInt + t0
                    dDeltaT_dmu0CM = dtInt_dmu0CM + dt0_dmu0CM
                    dDeltaT_domega0CM = dtInt_domega0CM + dt0_domega0CM
                End if
            Else ! Inbound
                Period = twoPi * a * Sqrt(a / muEarth)
                dPeriod_dmu0CM = 1.5_dp * Period * da_dmu0CM / a
                dPeriod_domega0CM = 1.5_dp * Period * da_domega0CM / a
                If (Direct) then
                    DeltaT = Period - thetaInt - theta0
                    dDeltaT_dmu0CM = dPeriod_dmu0CM - dtInt_dmu0CM - dt0_dmu0CM
                    dDeltaT_domega0CM = dPeriod_domega0CM - dtInt_domega0CM - dt0_domega0CM
                Else ! Indirect
                    DeltaT = Period - thetaInt + theta0
                    dDeltaT_dmu0CM = dPeriod_dmu0CM - dtInt_dmu0CM + dt0_dmu0CM
                    dDeltaT_domega0CM = dPeriod_domega0CM - dtInt_domega0CM + dt0_domega0CM
                End if
            End if
        Else !collinear
            da_dnuX = (dp_dnuX + 2.0_dp * e * a * de_dnuX) &
                & / (-delta_e * (1.0_dp +e))
            da_domega0CM = (dp_domega0CM + 2.0_dp * e * a * de_domega0CM) &
                & / (-delta_e * (1.0_dp +e))
            dE0_dnuX = ((r0/a**2) * da_dnuX - cos(theta0) * de_dnuX) / e
            dE0_domega0CM = ((r0/a**2) * da_domega0CM - cos(theta0) * de_domega0CM) / e
            dEInt_dnuX = ((rInt/a**2) * da_dnuX - cos(thetaInt) * de_dnuX) / e
            dEInt_domega0CM = ((rInt/a**2) * da_domega0CM - cos(thetaInt) * de_domega0CM) / e
            dt0_dnuX = 1.5 * t0 * da_dnuX / a + a * Sqrt(a/muEarth) &
                & * (r0 * dE0_dnuX / a - sin(E0) * de_dnuX)
            dt0_dnuY = 1.5 * t0 * da_dnuY / a + a * Sqrt(a/muEarth) &
                & * (r0 * dE0_dnuY / a - sin(E0) * de_dnuY)
            dtInt_dnuX = 1.5 * tInt * da_dnuX / a + a * Sqrt(a/muEarth) &
                & * (rInt * dEInt_dnuX / a - sin(EInt) * de_dnuX)
            dtInt_dnuY = 1.5 * tInt * da_dnuY / a + a * Sqrt(a/muEarth) &
                & * (rInt * dEInt_dnuY / a - sin(EInt) * de_dnuY)
            ! Combine to get time of flight and its derivatives
            If (Outbound) then
                If (Direct) then
                    DeltaT = tInt - t0
                    dDeltaT_dnuX = dtInt_dnuX - dt0_dnuX
                    dDeltaT_dnuY = dtInt_dnuY - dt0_dnuY
                Else ! Indirect
                    DeltaT = tInt + t0
                    dDeltaT_dnuX = dtInt_dnuX + dt0_dnuX
                    dDeltaT_dnuY = dtInt_dnuY + dt0_dnuY
                End if
            Else ! Inbound
                Period = twoPi * a * Sqrt(a / muEarth)
                dPeriod_dnuX = 1.5_dp * Period * da_dnuX / a
                dPeriod_dnuY = 1.5_dp * Period * da_dnuY / a
                If (Direct) then
                    DeltaT = Period - thetaInt - theta0
                    dDeltaT_dnuX = dPeriod_dnuX - dtInt_dnuX - dt0_dnuX
                    dDeltaT_dnuY = dPeriod_dnuY - dtInt_dnuY - dt0_dnuY
                Else ! Indirect
                    DeltaT = Period - thetaInt + theta0
                    dDeltaT_dnuX = dPeriod_dnuX - dtInt_dnuX + dt0_dnuX
                    dDeltaT_dnuY = dPeriod_dnuY - dtInt_dnuY + dt0_dnuY
                End if
            End if
        End if
    End Subroutine Elliptical_DeltaT
 
    Subroutine Near_Parabolic_DeltaT
        Real(dp):: xi0
        Real(dp):: x0, f0, df0_dx0
        Real(dp):: g0, a0, b0, c0, d0
        Real(dp):: xiInt
        Real(dp):: xInt, fInt, dfInt_dxInt
        Real(dp):: gInt, aInt, bInt, cInt, dInt
        If (Inbound) then
            Print *, "Error:: Near_Parabolic_DeltaT : Inbound Intercept"
            ! Can't be inbound: Has (or nearly has) escape speed
            Stop
        End if
        xi0 = Find_xi(delta_e, tan(theta0/2.0_dp))
        Call g_NPDT(delta_e, xi0, g0, a0, b0)
        t0 = p * h * g0 / muEarth
        xiInt = Find_xi(delta_e, tan(thetaInt/2.0_dp))
        Call g_NPDT(delta_e, xiInt, gInt, aInt, bInt)
        tInt = p * h * gInt / muEarth
        If (.not. collinear) then
            Call Near_Parabolic_T (xi0, t0, delta_e, de_dmu0CM, de_domega0CM, &
                & theta0, dtheta0_dmu0CM, dtheta0_domega0CM, &
                & p, dp_dmu0CM, dp_domega0CM, h, dh_dmu0CM, dh_domega0CM,&
                & dt0_dmu0CM, dt0_domega0CM)
            Call Near_Parabolic_T (xiInt, tInt, delta_e, de_dmu0CM, de_domega0CM, &
                & thetaInt, dthetaInt_dmu0CM, dthetaInt_domega0CM, &
                & p, dp_dmu0CM, dp_domega0CM, h, dh_dmu0CM, dh_domega0CM, &
                & dtInt_dmu0CM, dtInt_domega0CM)
            If (Direct) then
                DeltaT = tInt - t0
                dDeltaT_dmu0CM = dtInt_dmu0CM - dt0_dmu0CM
                dDeltaT_domega0CM = dtInt_domega0CM - dt0_domega0CM
            Else ! Indirect
                DeltaT = tInt + t0
                dDeltaT_dmu0CM = dtInt_dmu0CM + dt0_dmu0CM
                dDeltaT_domega0CM = dtInt_domega0CM + dt0_domega0CM
            End if
        Else
            Call Near_Parabolic_T (xi0, t0, delta_e, de_dnuX, de_dnuY, &
                & theta0, dtheta0_dnuX, dtheta0_dnuY, &
                & p, dp_dnuX, dp_dnuY, h, dh_dnuX, dh_dnuY,&
                & dt0_dnuX, dt0_dnuY)
            Call Near_Parabolic_T (xiInt, tInt, delta_e, de_dnuX, de_dnuY, &
                & thetaInt, dthetaInt_dnuX, dthetaInt_dnuY, &
                & p, dp_dnuX, dp_dnuY, h, dh_dnuX, dh_dnuY, &
                & dtInt_dnuX, dtInt_dnuY)
            If (Direct) then
                DeltaT = tInt - t0
                dDeltaT_dnuX = dtInt_dnuX - dt0_dnuX
                dDeltaT_dnuY = dtInt_dnuY - dt0_dnuY
            Else ! Indirect
                DeltaT = tInt + t0
                dDeltaT_dnuX = dtInt_dnuX + dt0_dnuX
                dDeltaT_dnuY = dtInt_dnuY + dt0_dnuY
            End if
        End if
        ! combine to get DeltaT and its derivatives
    End Subroutine Near_Parabolic_DeltaT

    Subroutine Hyperbolic_Near_Vertical_DeltaT
        Real(dp):: u0, du0_dmu0CM, du0_domega0CM, du0_dnuX, du0_dnuY
        Real(dp):: uInt, duInt_dmu0CM, duInt_domega0CM, duInt_dnuX, duInt_dnuY
        Real(dp):: a, b0, bInt
        If (.not. Direct .or. .not. Outbound) then
            Print *, "Error:: Hyperbolic_Near_Vertical_DeltaT: wrong trajectory type"
            Stop
        End if
        u0 = tau * r0 / muEarth
        uInt = tau * rInt / muEarth
        t0 = (muEarth / Sqrt(2.0_dp * tau**3)) &
            & * (Sqrt(u0 * (1.0_dp + u0)) - Log(Sqrt(u0)+Sqrt(1.0_dp + u0)))
        tInt = (muEarth / Sqrt(2.0_dp * tau**3)) &
            & * (Sqrt(uInt * (1.0_dp + uInt)) - Log(Sqrt(uInt)+Sqrt(1.0_dp + uInt)))
        DeltaT = tInt - t0
        a = muEarth / (tau**2 * Sqrt(2.0_dp * tau))
        b0 = (u0 * Sqrt(u0 / (1.0_dp + u0)) - 1.5_dp * t0)
        bInt = (uInt * Sqrt(uInt / (1.0_dp + uInt)) - 1.5_dp * tInt)
        ! Note: derivatives of tau are also derivatives of kappa0 because r0 is constant
        If(.not. Collinear) then
            du0_dmu0CM = u0 * dkappa0_dmu0CM / tau
            du0_domega0CM = u0 * dkappa0_domega0CM / tau
            dt0_dmu0cm = a * b0 * dkappa0_dmu0CM
            dt0_domega0CM = a * b0 * dkappa0_domega0CM
            duInt_dmu0CM = uInt * dkappa0_dmu0CM / tau
            duInt_domega0CM = uInt * dkappa0_domega0CM / tau
            dtInt_dmu0cm = a * bInt * dkappa0_dmu0CM
            dtInt_domega0CM = a * bInt * dkappa0_domega0CM
            dDeltaT_dmu0cm = dtInt_dmu0cm - dt0_dmu0cm
            dDeltaT_domega0CM = dtInt_domega0CM - dt0_domega0CM
        Else
            du0_dnuX = u0 * dkappa0_dnuX / tau
            du0_dnuY = u0 * dkappa0_dnuY / tau
            dt0_dnuX = a * b0 * dkappa0_dnuX
            dt0_dnuY = a * b0 * dkappa0_dnuY
            duInt_dnuX = uInt * dkappa0_dnuX / tau
            duInt_dnuY = uInt * dkappa0_dnuY / tau
            dtInt_dnuX = a * bInt * dkappa0_dnuX
            dtInt_dnuY = a * bInt * dkappa0_dnuY
            dDeltaT_dnuX = dtInt_dnuX - dt0_dnuX
            dDeltaT_dnuY = dtInt_dnuY - dt0_dnuY
        End if   
    End Subroutine Hyperbolic_Near_Vertical_DeltaT

    Subroutine Elliptical_Near_Vertical_DeltaT 
        Real(dp):: u0, du0_dmu0CM, du0_domega0CM, du0_dnuX, du0_dnuY
        Real(dp):: uInt, duInt_dmu0CM, duInt_domega0CM, duInt_dnuX, duInt_dnuY
        Real(dp):: b0, bInt
        Real(dp):: Period
        If (.not. Direct) then
            Print *, "Error:: Elliptical_Near_Vertical_DeltaT: must be Direct"
            Stop
        End if
        u0 = -tau * r0 / muEarth
        uInt = -tau * rInt / muEarth
        t0 = (muEarth / Sqrt(2.0_dp * (-tau)**3)) &
            & * (ASin(Sqrt(u0)) - Sqrt(u0 * (1.0_dp + u0)))
        tInt = (muEarth / Sqrt(2.0_dp * (-tau)**3)) &
            & * (ASin(Sqrt(uInt)) - Sqrt(uInt * (1.0_dp + uInt)))
        a = muEarth / (tau**2 * Sqrt(-2.0_dp * tau))
        b0 = (1.5_dp * t0 - u0 * Sqrt(u0 / (1.0_dp - u0)))
        bInt = (1.5_dp * tInt - uInt * Sqrt(uInt / (1.0_dp - uInt)))
        If(.not. Collinear) then
            du0_dmu0CM = u0 * dkappa0_dmu0CM / tau
            du0_domega0CM = u0 * dkappa0_domega0CM / tau
            dt0_dmu0cm = a * b0 * dkappa0_dmu0CM
            dt0_domega0CM = a * b0 * dkappa0_domega0CM
            duInt_dmu0CM = uInt * dkappa0_dmu0CM / tau
            duInt_domega0CM = uInt * dkappa0_domega0CM / tau
            dtInt_dmu0cm = a * bInt * dkappa0_dmu0CM
            dtInt_domega0CM = a * bInt * dkappa0_domega0CM
            If (Outbound) then
                DeltaT = tInt - t0
                dDeltaT_dmu0CM = dtInt_dmu0CM - dt0_dmu0CM
                dDeltaT_domega0CM = dtInt_domega0CM - dt0_domega0CM
            Else ! Inbound
                Period = twoPi * a * Sqrt(a / muEarth)
                dPeriod_dmu0CM = 1.5_dp * Period * da_dmu0CM / a
                dPeriod_domega0CM = 1.5_dp * Period * da_domega0CM / a
                DeltaT = Period - thetaInt - theta0
                dDeltaT_dmu0CM = dPeriod_dmu0CM - dtInt_dmu0CM - dt0_dmu0CM
                dDeltaT_domega0CM = dPeriod_domega0CM - dtInt_domega0CM - dt0_domega0CM
            End if
        Else ! Collinear
            du0_dnuX = u0 * dkappa0_dnuX / tau
            du0_dnuY = u0 * dkappa0_dnuY / tau
            dt0_dnuX = a * b0 * dkappa0_dnuX
            dt0_dnuY = a * b0 * dkappa0_dnuY
            duInt_dnuX = uInt * dkappa0_dnuX / tau
            duInt_dnuY = uInt * dkappa0_dnuY / tau
            dtInt_dnuX = a * bInt * dkappa0_dnuX
            dtInt_dnuY = a * bInt * dkappa0_dnuY
            If (Outbound) then
                DeltaT = tInt - t0
                dDeltaT_dnuX = dtInt_dnuX - dt0_dnuX
                dDeltaT_dnuY = dtInt_dnuY - dt0_dnuY
            Else ! Inbound
                Period = twoPi * a * Sqrt(a / muEarth)
                dPeriod_dnuX = 1.5_dp * Period * da_dnuX / a
                dPeriod_dnuY = 1.5_dp * Period * da_dnuY / a
                DeltaT = Period - thetaInt - theta0
                dDeltaT_dnuX = dPeriod_dnuX - dtInt_dnuX - dt0_dnuX
                dDeltaT_dnuY = dPeriod_dnuY - dtInt_dnuY - dt0_dnuY
            End if
        End if
    End Subroutine Elliptical_Near_Vertical_DeltaT

    Subroutine Parabolic_Near_Vertical_DeltaT
        Real(dp):: u0
        Real(dp):: uInt
        Real(dp):: dPoly_du0, dPoly_duInt
        Real(dp):: tVP0, tVPInt
        Real(dp), Parameter:: a(0:6) = &
            & (/ 1.0_dp, 3.0_dp/10.0_dp, 9.0_dp/56.0_dp, 5.0_dp/48.0_dp, &
            & 105.0_dp/1408.0_dp, 189.0_dp/3328.0_dp, 231.0_dp/5120.0_dp /)
        If (.not. Direct .or. .not. Outbound) then
            Print *, &
                & "Error:: Near_Parabolic_Near_Vertical_DeltaT: wrong trajectory type"
            Stop
        End if
        u0 = -tau * r0 / muEarth
        tVP0 = Sqrt(2.0_dp * r0**3 / (9.0_dp * muEarth))
        t0 = tVP0 * Polynomial(u0, a, 6, dpdx = dPoly_du0)
        uInt = -tau * rInt / muEarth
        tVPInt = Sqrt(2.0_dp * rInt**3 / (9.0_dp * muEarth))
        tInt = tVPInt * Polynomial(uInt, a, 6, dpdx = dPoly_duInt)
        DeltaT = tInt - t0
        If (.not. collinear) then
            dt0_dmu0CM = tVP0 * dPoly_du0 * u0 * dkappa0_dmu0CM / tau
            dt0_domega0CM = tVP0 * dPoly_du0 *  u0 * dkappa0_domega0CM / tau
            dtInt_dmu0CM = tVPInt * dPoly_duInt *  uInt * dkappa0_dmu0CM / tau
            dtInt_domega0CM = tVPInt * dPoly_duInt *  uInt * dkappa0_domega0CM / tau
            dDeltaT_dmu0CM = dtInt_dmu0CM - dt0_dmu0CM
            dDeltaT_domega0CM = dtInt_domega0CM - dt0_domega0CM    
        Else ! Collinear
            dt0_dnuX = tVP0 * dPoly_du0 * u0 * dkappa0_dnuX / tau
            dt0_dnuY = tVP0 * dPoly_du0 *  u0 * dkappa0_dnuY / tau
            dtInt_dnuX = tVPInt * dPoly_duInt *  uInt * dkappa0_dnuX / tau
            dtInt_dnuY = tVPInt * dPoly_duInt *  uInt * dkappa0_dnuY / tau
            dDeltaT_dnuX = dtInt_dnuX - dt0_dnuX
            dDeltaT_dnuY = dtInt_dnuY - dt0_dnuY    
        End if
    End Subroutine Parabolic_Near_Vertical_DeltaT

    

    Subroutine g_NPDT(xi, dxi, g, a, b)
        Real(dp), Intent(In):: xi, dxi
        Real(dp), Intent(Out):: g
        Real(dp), Intent(Out), Optional:: a, b
        Real(dp):: y, dy, yPowers(0:nMax), p1, p2, dp1_dy, dp2_dy
        y = delta_e * xi**2
        yPowers = Powers_Vector(y, nMax)
        p1 = Polynomial(yPowers, &
            & InverseFactorial(1:2*nMax:2), nMax, dpdx = dp1_dy)
        p2 = Polynomial(yPowers(0:nMax-1), &
            & InverseFactorial(3:2*nMax:2), nMax-1, dpdx = dp2_dy)
        g = xi * (p1 + xi**2 * p2)
        If (Present(a) .and. Present(b)) then
            a = p1 + 3.0_dp * xi**2 * p2
            b = xi * (dp1_dy + xi**2 * dp2_dy)
        End if
    End Subroutine g_NPDT

    
    ! Subroutine Numerical_Quadrature_DeltaT_DeltaTheta
        ! Add later for testing and possible use ?
    ! End Subroutine Numerical_Quadrature_DeltaT_DeltaTheta  
    
End Function Detection_Solid_Angle_CM

! DivergenceFactor Helper Routines (not contained)

Subroutine Near_Parabolic_T (xi, t, delta_e, de_dchi1, de_dchi2, &
    & theta, dtheta_dchi1, dtheta_dchi2, p, dp_dchi1, dp_dchi2, &
    & h, dh_dchi1, dh_dchi2, dt_dchi1, dt_dchi2)
    Real(dp), Intent(In):: xi
    Real(dp), Intent(In):: t
    Real(dp), Intent(In):: delta_e, de_dchi1, de_dchi2
    Real(dp), Intent(In):: theta, dtheta_dchi1, dtheta_dchi2
    Real(dp), Intent(In):: p, dp_dchi1, dp_dchi2
    Real(dp), Intent(In):: h, dh_dchi1, dh_dchi2
    Real(dp), Intent(Out):: dt_dchi1, dt_dchi2
    Real(dp):: x, f, df_dx
    Real(dp):: a, b, c, d, e
    Real(dp):: dxi_dchi1, dxi_dchi2
    Real(dp):: dg_dchi1, dg_dchi2
    e = 1.0_dp + delta_e
    x = -delta_e * xi**2 / 4.0_dp
    Call f_NPDT(x, f, df_dx)
    c = xi**4 * df_dx / 4.0_dp &
        & - tan(theta/2.0_dp)  / ((1.0_dp + e) * Sqrt(1.0_dp + e))
    d = cos(theta/2.0_dp)**2 * Sqrt(1.0_dp + e)
    dxi_dchi1 = (c * de_dchi1 + dtheta_dchi1 / d) &
        & / (f + 3.0_dp * x * df_dx)
    dxi_dchi2 = (c * de_dchi2 + dtheta_dchi2 / d) &
        & / (f + 3.0_dp * x * df_dx)
    dg_dchi1 = a * dxi_dchi1 + b * &
        & (xi**2 * de_dchi1 + 2.0_dp * delta_e * xi * dxi_dchi1)
    dg_dchi2 = a * dxi_dchi2 + b * &
        & (xi**2 * de_dchi2 + 2.0_dp * delta_e * xi * dxi_dchi2)
    dt_dchi1 = t * (dp_dchi1 / p + dh_dchi1 / h &
        & - 1.5_dp * de_dchi1 / (1.0_dp + e)) &
        & + p * h * dg_dchi1 / (muEarth * (1.0_dp + e) * Sqrt(1.0_dp + e))
    dt_dchi2 = t * (dp_dchi2 / p + dh_dchi2 / h &
        & - 1.5_dp * de_dchi2 / (1.0_dp + e)) &
        & + p * h * dg_dchi2 / (muEarth * (1.0_dp + e) * Sqrt(1.0_dp + e))
End Subroutine Near_Parabolic_T

Function Find_xi(de, tanHalfTheta) Result (xi)
    Real(dp), Intent(In):: de   ! delta_e
    Real(dp), Intent(In):: tanHalfTheta
    Real(dp), Parameter:: f_RelTol = 1.0E-8_dp
    Real(dp):: xi, xiOld
    Real(dp):: a, f
    a = 2.0_dp * tanHalfTheta / Sqrt(2.0_dp + de)
    xi = Sqrt(2.0_dp) * tanHalfTheta
    Do
        xiOld = xi
        Call f_NPDT(-de * (xi / 2.0)**2, f)
        xi = a / f
        If (abs(xi-xiOld) <= f_RelTol * xi) Return
    End do
End Function Find_xi


Subroutine f_NPDT(x, f, df_dx)
    Real(dp), Intent(In):: x
    Real(dp), Intent(Out):: f
    Real(dp), Intent(Out), Optional:: df_dx
    Real(dp):: xPowers(0:nMax)
    Integer, Parameter:: nm = 9
    Real(dp):: f_coef(0:nm)
    Real(dp), Parameter:: f_coef_Numerators(0:nm) = &
        & (/ 1.0_dp, 1.0_dp, 2.0_dp, 17.0_dp, 62.0_dp, 1382.0_dp, 21844.0_dp, &
        & 929569.0_dp, 443861162.0_dp, 18888466084.0_dp /)
    Real(dp), Parameter:: f_coef_Denominators(0:nm) = &
        & (/ 1.0_dp, 3.0_dp, 15.0_dp, 315.0_dp, 2835.0_dp, 155925.0_dp, 6081075.0_dp, &
        & 638512875.0_dp, 1856156927625.0_dp, 194896477400625.0_dp /)
    f_coef = f_coef_Numerators / f_coef_Denominators
    If (Present (df_dx)) then
        f = Polynomial(x, f_coef, nMax, xPowers, df_dx)
    Else
        f = Polynomial(x, f_coef, nMax, xPowers)
    End if
End Subroutine f_NPDT


! UTILITIES

    ! Vector Utilities

Function Hat(x)
    Real(dp):: x(1:3)
    Real(dp):: Hat(1:3)
    Hat = x / Sqrt(Dot_Product(x,x))
End Function Hat

Function Norm(x)
    Real(dp):: x(:)
    Real(dp):: Norm
    Norm = Sqrt(Dot_Product(x,x))
End Function Norm


Function Cross_Product(x, y) Result (cp)
    Real(dp), Intent(In):: x(1:3), y(1:3)
    Real(dp):: cp(1:3)
    cp(1) = x(2) * y(3) - x(3) * y(2)
    cp(2) = x(3) * y(1) - x(1) * y(3)
    cp(3) = x(1) * y(2) - x(2) * y(1)
End Function Cross_Product

    ! Polynomial Utilities

Function Powers_Vector(x, n) Result (xPowers)
    Real(dp), Intent(In):: x
    Integer, Intent(In) :: n
    Real(dp):: xPowers(0:n)
    Integer:: k
    xPowers(0) = 1._dp
    Do k = 1, n
        xPowers(k) = x * xPowers(k-1)
    End do
End Function Powers_Vector

Function Polynomial_Scalar(x, coefficients, n, xPowers, dpdx) Result(p)
    Integer, Intent(In) :: n
    Real(dp), Intent(In):: x
    Real(dp), Intent(In):: coefficients(0:n)
    Real(dp), Intent(Out), Optional:: xPowers(0:n)
    REal(dp), Intent(Out), Optional:: dpdx
    Real(dp):: p, xp(0:n)
    xp = Powers_Vector(x, n)
    If (Present(xPowers)) xPowers = xp
    p = Dot_Product(xp, coefficients)
    If (Present(dpdx)) then
        dpdx = Dot_Product(xp(0:nMax-1), &
            & Natural_Number(1:nMax) * coefficients(1:nMax))
    End if
End Function Polynomial_Scalar

Function Polynomial_Vector(xPowers, coefficients, n, dpdx) Result(p)
    Integer, Intent(In) :: n
    Real(dp), Intent(In):: xPowers(0:n)
    Real(dp), Intent(In):: coefficients(0:n)
    Real(dp), Intent(Out), Optional:: dpdx
    Real(dp):: p
    p = Dot_Product(xPowers, coefficients)
    If (Present(dpdx)) dpdx = Dot_Product(xPowers(0:nMax-1), &
        & Natural_Number(1:nMax) * coefficients(1:nMax))
End Function Polynomial_Vector

Function Polynomial_Vector_Matrix(xPowers, coefficients, n, m) Result(p)
    Integer, Intent(In) :: n, m
    Real(dp), Intent(In):: xPowers(0:n)
    Real(dp), Intent(In):: coefficients(0:n,1:m)
    Real(dp):: p(1:m)
    p = MatMul(xPowers, coefficients)
End Function Polynomial_Vector_Matrix

Function Delta_Polynomial_Coefficients(xPowers, a, n) Result(b)
    Integer, Intent(In) :: n
    Real(dp), Intent(In):: xPowers(0:n)
    Real(dp), Intent(In):: a(0:n)
    Real(dp):: b(0:n)
    Integer:: j, k
    b = 0.0_dp
    Do k = 1, n
        Do j = k, n
            b(k) = b(k) + BinomialCoefficient(j,k) * a(j) * xPowers(j-k)
        End do
    End do
End Function Delta_Polynomial_Coefficients

Function Delta_Polynomial_1(xPowers, a, dx, n) Result (Delta_P)
    Integer, Intent(In) :: n
    Real(dp), Intent(In):: a(0:n)
    Real(dp), Intent(In):: xPowers(0:n)
    Real(dp), Intent(In):: dx
    Real(dp):: b(0:n)
    Real(dp):: Delta_P
    Real(dp):: dxPowers(0:n) 
    dxPowers = Powers_Vector(dx, n)
    b = Delta_Polynomial_Coefficients(xPowers, a, n)
    Delta_P = Polynomial(xPowers, b, n)
End Function Delta_Polynomial_1

Function Delta_Polynomial_2 (dx, b, n) Result(Delta_P)
    Integer, Intent(In) :: n
    Real(dp), Intent(In):: dx
    Real(dp), Intent(In):: b(0:n)
    Real(dp):: Delta_P
    Delta_P = Polynomial(Powers_Vector(dx, n), b, n)
End Function Delta_Polynomial_2

! Legendre Utilities

    ! Whit: Thanks for suggesting reordering the cases for efficiency. It is done in these routines.

Subroutine LegendreP_Scalar(x, P, L)
    ! This routine computes LegendreP(L,x) recursively
    !   discarding all the Legendre polynomials of lower order than L
    ! It is very inefficient to call it repeatedly to get P(0,x) through P(L,X)
    !   Use the vector routine which does not waste work.
    ! Better, use the matrix routine,
    !   which vectorizes the arithmetic for multiple arguments x(1:N)
    Real(dp), Intent(In):: x            ! Argument
    Real(dp), Intent(Out):: P           ! P = LegendreP(L, x)
    Integer, Intent(In):: L             ! Order of Legendre Polynomial
    Integer:: j
    Real(dp):: P_prior
    Real(dp):: P_current
    Select Case (L)
        Case(2:)
            P = 1.0_dp
            P_current = P
            P = x
            Do j = 1, L-1
                P_prior = P_current
                P_current = P
                P = (Real(2*j+1,dp)/Real(j+1,dp)) * x * P_current - (Real(j,dp)/Real(j+1,dp)) * P_prior
            End do
        Case (1)
            P = x
        Case(0)
            P = 1.0_dp
        End Select
End Subroutine LegendreP_Scalar


Subroutine LegendreP_Vector(x, P, L)
    ! More efficient to use the matrix version rather than calling this repeatedly,
    !   because it vectorizes and this is sequential.
    Integer, Intent(In):: L             ! Max order of Legendre Polynomials
    Real(dp), Intent(In):: x            ! Argument
    Real(dp), Intent(Out):: P(0:L)      ! P(j) = LegendreP(j, x)
    Integer:: j
    Select Case (L)
        Case (2:)
            P(0) = 1.0_dp
            P(1) = x
            Do j = 1, L-1
                P(j+1) = (Real(2*j+1,dp)/Real(j+1,dp)) * x * P(j) - (Real(j,dp)/Real(j+1,dp)) * P(j-1)
            End do
        Case(1)
            P(0) = 1.0_dp
            P(1) = x
        Case(0)
            P(0) = 1.0_dp
        End Select
    
End Subroutine LegendreP_Vector

Subroutine LegendreP_Matrix(x, P, L, N)
    Integer, Intent(In):: L     ! Max order of Legendre Polynomials
    Integer, Intent(In):: N     ! Number of arguments
    Real(dp), Intent(In):: x(1:N)       ! Arguments
    Real(dp), Intent(Out):: P(1:N,0:L)  ! P(k,j) = LegendreP(j, x(k))
        ! The inconvenient subscript order is required so that the vector operations are on P(:,j), not P(j,:),
        !     for sequential storage of vector arithmetic arguments.
    Integer:: j
    Select Case (L)
        Case (2:)
            P(:,0) = 1.0_dp
            P(:,1) = x
            Do j = 1, L-1
                P(:,j+1) = (Real(2*j+1,dp)/Real(j+1,dp)) * x(:) * P(:,j) - (Real(j,dp)/Real(j+1,dp)) * P(:,j-1)
            End do
        Case(1)
            P(:,0) = 1.0_dp
            P(:,1) = x
        Case(0)
            P(:,0) = 1.0_dp
        End Select        
End Subroutine LegendreP_Matrix



! Combinatoric Utilities

Subroutine Initialize_Binomial_Coefficients(n, BC)
    Integer, Intent(In):: n
    Real(dp), Intent(Out):: BC(0:n, 0:n)
    Integer :: i, j, k
    BC = 0.0_dp
    Do k = 0, n
        Do j = k, n
            Do i = 1, k
                BC(j,k) = BC(j,k) * (Real(j+1-i,dp) / Real(i,dp))
            End do ! i
            BC(j, j-k) = BC(j,k)
        End do ! j
    End do ! k
End Subroutine Initialize_Binomial_Coefficients

Subroutine Initialize_Factorials (n, f, f_inv)
    Integer, Intent(In):: n
    Real(dp), Intent(Out):: f(0:2*nMax), f_inv(0:2*nMax)
    Integer:: k
    f(0) = 1.0_dp
    Do k = 1, n
        f(k) = f(k-1) * Real(k,dp)
    End do
    f_inv = 1.0_dp / f
End Subroutine Initialize_Factorials

    ! Conversions among deltaSin, deltaCos, deltaX
    ! Well-conditioned for small abs(deltaX)

Function deltaSin2(X, deltaX)
    Real(dp):: X, deltaX
    Real(dp):: deltaSin2
    deltaSin2 = cos(X) * sin(deltaX) - sin(X) * sin(deltaX)**2 / (1 + cos(deltaX))
End Function deltaSin2

Function deltaCos2(X, deltaX)
    Real(dp):: X, deltaX
    Real(dp):: deltaCos2
    deltaCos2 = -sin(X) * sin(deltaX) + cos(X) * sin(deltaX)**2 / (1 + cos(deltaX))
End Function deltaCos2

Function deltaSin3(sinX, cosX, deltaCosX)
    Real(dp):: sinX, cosX, deltaCosX
    Real(dp):: deltaSin3
    deltaSin3 = -(2.0_dp * (cosX+deltaCosX) * deltaCosX) &
        & / (sinX + Sqrt(sinX**2 - (2.0_dp * cosX + deltaCosX) * deltaCosX))
End Function deltaSin3

Function deltaCos3(sinX, cosX, deltaSinX)
    Real(dp):: sinX, cosX, deltaSinX
    Real(dp):: deltaCos3
    deltaCos3 = -(2.0_dp * (sinX+deltaSinX) * deltaSinX) &
        & / (cosX + Sqrt(cosX**2 - (2.0_dp * sinX + deltaSinX) * deltaSinX))
End Function deltaCos3

Function deltaX(sinX, cosX, deltaSinX, relTol)
    Real(dp):: sinX, cosX, deltaSinX
    Real(dp), Optional:: relTol
    Real(dp):: deltaX
    Real(dp):: y, yOld
    Real(dp):: rTol
    If (Present(relTol)) then
        rTol = relTol
    Else
        rTol = 1.0E-15_dp
    End If
    y = 0.0_dp
    Do
        yOld = y
        y = 2.0_dp * deltaSinX &
            & / (cosX + Sqrt(cosX**2 - 4.0_dp * sinX * deltaSinX / (1.0_dp + Sqrt(1.0_dp - y**2))))
        If (Abs(y-yOld) <= rTol * abs(y)) Exit
    End do
    deltaX = asin(y)
End Function deltaX

End Module Diverge_exact
