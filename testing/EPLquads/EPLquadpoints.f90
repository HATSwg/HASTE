Program EPLquadpoints

Use Kinds, Only: dp
Use Global, Only: R_earth
Use US_Std_Atm_1976, Only: rho
Use US_Std_Atm_1976, Only: Zb
Use Quadratures, Only: Romberg_Quad
Use Quadratures, Only: GaussLegendreN
Use Utilities, Only: Prec

Implicit None

Real(dp) :: z0,z1,r0,zeta0,rho0,inv_rho0
Integer :: b,c
Integer, Parameter :: n_zeta = 10
Integer, Parameter :: n_max = 36
Real(dp) :: dZ,Smax
Real(dp) :: Ls
Real(dp) :: Lz
Real(dp) :: L0,L0p
Real(dp) :: reltol,abstol
Integer :: unit
Character(1) :: sub

rho0 = rho(0._dp)
inv_rho0 = 1._dp / rho0
reltol = 1.E-12_dp
abstol = 0._dp
Open(NEWUNIT=unit,FILE='EPLquadPrecs.tst',ACTION='WRITE',STATUS='REPLACE')

!Layers 1-7 (b = 0-6) are all single layers in the lower atmosphere
Do b = 0,6
    sub = ''
    z0 = Zb(b)
    z1 = Zb(b+1)
    r0 = R_earth + z0
    dZ = z1 - z0
    Call CheckEPLprecs()
End Do

!Layer 8 (b = 7) goes from 86 to 91 km
!  There is a discontinuity in the USSA76 model at 86km, so the EPL integrand 
!  functions have a special case for this to use the right endpoint
b = 7
sub = ''
z0 = Zb(b)
z1 = Zb(b+1)
r0 = R_earth + z0
dZ = z1 - z0
Call CheckEPLprecs()

!Layer 9 (b=8) has 4 sublayers
b = 8
sub = ''
z0 = Zb(b)
z1 = Zb(b+1)
r0 = R_earth + z0
dZ = z1 - z0
Call CheckEPLprecs()
Do c = 1,4
    Select Case (c)
        Case (1)
            sub = 'a'
            z0 = Zb(b)
            z1 = 95._dp
        Case (2)
            sub = 'b'
            z0 = 95._dp
            z1 = 97._dp
        Case (3)
            sub = 'c'
            z0 = 97._dp
            z1 = 100._dp
        Case (4)
            sub = 'd'
            z0 = 100._dp
            z1 = Zb(b+1)
    End Select
    r0 = R_earth + z0
    dZ = z1 - z0
    Call CheckEPLprecs()
End Do

!Layer 10 (b=9) has 2 sublayers
b = 9
sub = ''
z0 = Zb(b)
z1 = Zb(b+1)
r0 = R_earth + z0
dZ = z1 - z0
Call CheckEPLprecs()
Do c = 1,2
    Select Case (c)
        Case (1)
            sub = 'a'
            z0 = Zb(b)
            z1 = 115._dp
        Case (2)
            sub = 'b'
            z0 = 115._dp
            z1 = Zb(b+1)
    End Select
    r0 = R_earth + z0
    dZ = z1 - z0
    Call CheckEPLprecs()
End Do

!Layer 11 (b=10) has 2 sublayers
b = 10
sub = ''
z0 = Zb(b)
z1 = Zb(b+1)
r0 = R_earth + z0
dZ = z1 - z0
Call CheckEPLprecs()
Do c = 1,2
    Select Case (c)
        Case (1)
            sub = 'a'
            z0 = Zb(b)
            z1 = 500._dp
        Case (2)
            sub = 'b'
            z0 = 500._dp
            z1 = Zb(b+1)
    End Select
    r0 = R_earth + z0
    dZ = z1 - z0
    Call CheckEPLprecs()
End Do

Close(unit)

Contains

Subroutine CheckEPLprecs()
    Implicit None

    Write(*,'(A)') '------------------------------------------------------------------------'
    Write(*,'(A,I0,A,F0.3,A,F0.3,A)') 'LAYER ',b+1,Trim(sub)//',  z = ',z0,' to ',z1,' km'
    Write(*,'(A)') '------------------------------------------------------------------------'
    Write(unit,'(A)') '------------------------------------------------------------------------'
    Write(unit,'(A,I0,A,F0.3,A,F0.3,A)') 'LAYER ',b+1,Trim(sub)//',  z = ',z0,' to ',z1,' km'
    Write(unit,'(A)') '------------------------------------------------------------------------'
    !Shallow zetas
    Call CheckEPLprecsShallow()
    !Steep zetas
    Call CheckEPLprecsSteep()
End Subroutine CheckEPLprecs

Subroutine CheckEPLprecsShallow()
    Use Kinds, Only: dp
    Implicit None
    Integer :: i,j,p
    Character(2) :: p_char
    Character(11) :: fmt_char

    Write(*,'(A)') 'zeta0 (shallow)'
    Write(*,'(A)') '---------------'
    Write(unit,'(A)') 'zeta0 (shallow)'
    Write(unit,'(A)') '---------------'
    Do i = 0,n_zeta-1
        !Compute exact EPL
        zeta0 = Real(i,dp) * 0.01_dp
        Smax = dZ * (2._dp * r0 + dZ) / ( zeta0 * r0 + Sqrt( (zeta0 * r0)**2 + dZ * (2._dp * r0 + dZ) ) )
        L0 = Romberg_Quad(EPL_Integrand_dS,0._dp,Smax,atol=abstol,rtol=reltol,p=L0p)
        Write(*,'(F4.2,A,ES23.15,A,F0.5,A)') zeta0,'   exact: ',L0,' (~',L0p,' digits of precision)'
        Write(unit,'(F4.2,A,ES23.15,A,F0.5,A)') zeta0,'   exact: ',L0,' (~',L0p,' digits of precision)'
        !Compute approximate EPL, stopping when desired precision is achieved
        Do p = 12,6,-3
            If (L0p .LT. Real(p,dp)) Cycle
            Write(p_char,'(I2.2)') p
            Write(fmt_char,'(A,I0,A,I0,A,I0)') 'ES',(p-1)+8,'.',p-1,',A',16-p
            Do j = 3,n_max
                Ls = GaussLegendreN(j,EPL_Integrand_dS,0._dp,Smax)
                Write(*,'(A,'//fmt_char//',A,I2)',ADVANCE='NO') ACHAR(13)//'      dS-p'//p_char//': ',Ls,'',': ',j
                If (Floor(Prec(Ls,L0)) .GE. p) Then
                    If (Floor(Prec(GaussLegendreN(j+1,EPL_Integrand_dS,0._dp,Smax),L0)) .GE. Floor(Prec(Ls,L0))) Then
                        Write(*,'(A,I0)') ' gauss-pts for prec > ',p
                        Write(unit,'(A,'//fmt_char//',A,I2,A,I0)') '      dS-p'//p_char//': ',Ls,'',': ',j, & 
                                                                 & ' gauss-pts for prec > ',p
                        Exit
                    End If
                End If
                If (j .EQ. n_max) Then
                    Write(*,'(A,'//fmt_char//',A,I0,A,I0,A)') ACHAR(13)//'      dS-p'//p_char//': ',Ls,'',': Prec ',p, & 
                                                            & ' NOT MET w/ ',j,' gauss-pts'
                    Write(unit,'(A,'//fmt_char//',A,I0,A,I0,A)') '      dS-p'//p_char//': ',Ls,'',': Prec ',p, & 
                                                               & ' NOT MET w/ ',j,' gauss-pts'
                End If
            End Do
        End Do
    End Do
    Write(*,*)
    Write(Unit,*)
End Subroutine CheckEPLprecsShallow

Subroutine CheckEPLprecsSteep
    Use Kinds, Only: dp
    Implicit None
    Integer :: i,j,p
    Character(2) :: p_char
    Character(11) :: fmt_char

    Write(*,'(A)') 'zeta0 (steep)'
    Write(*,'(A)') '---------------'
    Write(unit,'(A)') 'zeta0 (steep)'
    Write(unit,'(A)') '---------------'
    Do i = 1,n_zeta
        zeta0 = Real(i,dp) * 0.1_dp
        Smax = dZ * (2._dp * r0 + dZ) / ( zeta0 * r0 + Sqrt( (zeta0 * r0)**2 + dZ * (2._dp * r0 + dZ) ) )
        L0 = Romberg_Quad(EPL_Integrand_dS,0._dp,Smax,atol=abstol,rtol=reltol,p=L0p)
        Write(*,'(F4.2,A,ES23.15,A,F0.5,A)') zeta0,'   exact: ',L0,' (~',L0p,' digits of precision)'
        Write(unit,'(F4.2,A,ES23.15,A,F0.5,A)') zeta0,'   exact: ',L0,' (~',L0p,' digits of precision)'
        !Compute approximate EPL, stopping when desired precision is achieved
        Do p = 12,6,-3
            If (L0p .LT. Real(p,dp)) Cycle
            Write(p_char,'(I2.2)') p
            Write(fmt_char,'(A,I0,A,I0,A,I0)') 'ES',(p-1)+8,'.',p-1,',A',16-p
            Do j = 3,n_max
                Ls = GaussLegendreN(j,EPL_Integrand_dS,0._dp,Smax)
                Write(*,'(A,'//fmt_char//',A,I2)',ADVANCE='NO') ACHAR(13)//'      dS-p'//p_char//': ',Ls,'',': ',j
                If (Floor(Prec(Ls,L0)) .GE. p) Then
                    If (Floor(Prec(GaussLegendreN(j+1,EPL_Integrand_dS,0._dp,Smax),L0)) .GE. Floor(Prec(Ls,L0))) Then
                        Write(*,'(A,I0)') ' gauss-pts for prec > ',p
                        Write(unit,'(A,'//fmt_char//',A,I2,A,I0)') '      dS-p'//p_char//': ',Ls,'',': ',j, & 
                                                                 & ' gauss-pts for prec > ',p
                        Exit
                    End If
                End If
                If (j .EQ. n_max) Then
                    Write(*,'(A,'//fmt_char//',A,I0,A,I0,A)') ACHAR(13)//'      dS-p'//p_char//': ',Ls,'',': Prec ',p, & 
                                                            & ' NOT MET w/ ',j,' gauss-pts'
                    Write(unit,'(A,'//fmt_char//',A,I0,A,I0,A)') '      dS-p'//p_char//': ',Ls,'',': Prec ',p, & 
                                                               & ' NOT MET w/ ',j,' gauss-pts'
                End If
            End Do
            Do j = 3,n_max
                Lz = GaussLegendreN(j,EPL_Integrand_dZ,0._dp,dZ)
                Write(*,'(A,'//fmt_char//',A,I2)',ADVANCE='NO') ACHAR(13)//'      dZ-p'//p_char//': ',Lz,'',': ',j
                If (Floor(Prec(Lz,L0)) .GE. p) Then
                    If (Floor(Prec(GaussLegendreN(j+1,EPL_Integrand_dZ,0._dp,dZ),L0)) .GE. Floor(Prec(Lz,L0))) Then
                        Write(*,'(A,I0)') ' gauss-pts for prec > ',p
                        Write(unit,'(A,'//fmt_char//',A,I2,A,I0)') '      dZ-p'//p_char//': ',Lz,'',': ',j, & 
                                                                 & ' gauss-pts for prec > ',p
                        Exit
                    End If
                End If
                If (j .EQ. n_max) Then
                    Write(*,'(A,'//fmt_char//',A,I0,A,I0,A)') ACHAR(13)//'      dZ-p'//p_char//': ',Lz,'',': Prec ',p, & 
                                                            & ' NOT MET w/ ',j,' gauss-pts'
                    Write(unit,'(A,'//fmt_char//',A,I0,A,I0,A)') '      dZ-p'//p_char//': ',Lz,'',': Prec ',p, & 
                                                               & ' NOT MET w/ ',j,' gauss-pts'
                End If
            End Do
        End Do
    End Do
    Write(*,*)
    Write(Unit,*)
End Subroutine CheckEPLprecsSteep

Function EPL_Integrand_dS(s) Result(f)
    Use Kinds, Only: dp
    Use Utilities, Only: Smaller_Quadratic_Root
    Use US_Std_Atm_1976, Only: rho
    Implicit None
    Real(dp) :: f
    Real(dp), Intent(In) :: s
    Real(dp) :: deltaZ
    
    deltaZ = Smaller_Quadratic_root(r0,s*(2._dp*r0*zeta0 + s))
    If (b.EQ.6 .OR. b.EQ.7) Then
        f = rho(z0 + deltaZ, layer = b+1)  !special case for b=6 and b=7, for 86km discontinuity
    Else
        f = rho(z0 + deltaZ)
    End If
    f = f * inv_rho0
End Function EPL_Integrand_dS

Function EPL_Integrand_dZ(deltaZ) Result(f)
    Use Kinds, Only: dp
    Use US_Std_Atm_1976, Only: rho
    Implicit None
    Real(dp) :: f
    Real(dp), Intent(In) :: deltaZ
    
    If (b.EQ.6 .OR. b.EQ.7) Then
        f = rho(z0 + deltaZ, layer = b+1)  !special case for b=6 and b=7, for 86km discontinuity
    Else
        f = rho(z0 + deltaZ)
    End If
    f = inv_rho0 * f * (r0 + deltaZ) / Sqrt((r0*zeta0)**2 + 2._dp*r0*deltaZ + deltaZ**2)
End Function EPL_Integrand_dZ

Function rho_Isotherm(Z) Result(d)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: d
    Real(dp), Intent(In) :: Z
    !Isothermal atmosphere parameter
    Real(dp), Parameter :: scale_Height_conv = 0.02927176966650182_dp
    
    d = rho0 * Exp(-Z / (scale_Height_conv * 273.15_dp))
End Function rho_Isotherm

End Program EPLquadpoints
