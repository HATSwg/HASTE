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
Module Atmospheres

    Use Kinds, Only: dp
    Implicit None
    Private
    Public :: Atmosphere_Type
    Public :: EPL_Layer_Data
    Public :: Setup_Atmosphere
    Public :: Write_Atmosphere
    Public :: rho_SL,inv_rho_SL

    Type :: EPL_Layer_Data
        !STRAIGHT PATHS
        !Q-points and precomputed values for zeta larger than 0.1 (integrate dZ)
        Integer :: nZ  !number of q-points
        Real(dp), Allocatable :: uZ(:)  !abscissa
        Real(dp), Allocatable :: wZ(:)  !weights
        Real(dp) :: A  !precomputed for full layer
        Real(dp), Allocatable :: B(:)  !precomputed for full layer
        Real(dp) :: C  !precomputed for full layer
        Real(dp), Allocatable :: D(:)  !precomputed for full layer
        !Q-points for zeta smaller than 0.1 (integrate dS)
        Integer :: nS  !number of q-points
        Real(dp), Allocatable :: uS(:)  !abscissa
        Real(dp), Allocatable :: wS(:)  !weights
        !ORBITAL PATHS
        !Q-points and precomputed values for zeta larger than 0.1 (integrate dR)
        Integer :: nRk  !number of q-points
        Real(dp), Allocatable :: uRk(:)  !abscissa
        Real(dp), Allocatable :: wRk(:)  !weights
        Real(dp), Allocatable :: Bk(:)  !precomputed for full layer
        Real(dp), Allocatable :: Ck(:)  !precomputed for full layer    
        !Q-points for zeta smaller than 0.1 (integrate dTheta)
        Integer :: nTk  !number of q-points
        Real(dp), Allocatable :: uTk(:)  !abscissa
        Real(dp), Allocatable :: wTk(:)  !weights
    End Type

    Type :: Atmosphere_Type
        Integer :: model_index
        Integer :: comp_index
        Integer :: n_el  !number of elements in the atmospheric representation
        Integer, Allocatable :: iso_ind(:)  !index of the first isotope of each element in iso_frac array, has dim 1:n_el+1 (the last entry is a ghost)
        Integer, Allocatable :: iso_map(:)  !maps each isotope to an element index, has dim 1:n_isotopes
        Real(dp), Allocatable :: iso_frac(:)  !isotopic fractions, has dim 1:number_of_isotopes
        Logical, Allocatable :: diatomic(:)  !flag indicating whether this isotope forms a diatomic atmospheric constituent, has dim 1:number_of_isotopes
        Real(dp) :: ref_density_ratio
        Real(dp) :: Isothermal_temp  ![K]
        Real(dp) :: z_top  ![km]
        Real(dp) :: z_bot  ![km]
        Real(dp) :: R_top  ![km]
        Real(dp) :: R_bot  ![km]
        Real(dp) :: wind_AF(1:3)  ![km/s]
        Logical :: discontinuous  !flag indicates an atmosphere model with discontinuities (layers)
        Real(dp), Allocatable :: Zb(:) !a local set of layer boundaries from atmospheres module, augmented with any natural breaks not explicit in the atmosphere model, used for finding layer boundaries
        Integer, Allocatable :: bZb(:) !list of base indexes mapping to the atmosphere model's base indexes
        Real(dp), Allocatable :: Rb(:)
        Integer :: iZb(1:3)  !indexes for bottom and top layers to be included (layers for z_bot and z_top), and number of layers
        Integer :: iZb_map(1:3)  !indexes for bottom and top layers and number of layers IN THE ACTUAL ATMOSPHERE MODEL
        Type(EPL_Layer_Data), Allocatable :: EPL_Lay(:)
#       if CHECK_L
            Integer :: EPL_prec
#       endif
    Contains
        Procedure, Pass :: T => Atm_Temperature  !returns atmosphere temperature at specified height
        Procedure, Pass :: rho => Atm_Density  !returns atmospheric density at specified heigh
    End Type

    !Integer designator for atmosphere model choice, can have value equal to one of the following parameters
    Integer, Parameter :: atm_mod_IsoTherm = 31
    Integer, Parameter :: atm_mod_USstd1976 = 33
    Integer, Parameter :: atm_mod_None = -1

    !Integer designator for atmosphere composition choice, can have value equal to one of the following parameters
    Integer, Parameter :: atm_comp_ALL = 121
    Integer, Parameter :: atm_comp_N14_O16 = 123
    Integer, Parameter :: atm_comp_N14_O16_Ar40 = 125
    Integer, Parameter :: atm_comp_N14 = 127
    Integer, Parameter :: atm_comp_N15 = 129
    Integer, Parameter :: atm_comp_O16 = 131
    Integer, Parameter :: atm_comp_O17 = 133
    Integer, Parameter :: atm_comp_O18 = 135
    Integer, Parameter :: atm_comp_Ar40 = 137
    Integer, Parameter :: atm_comp_None = -1

    !Uniform, Isothermal, and shared parameter
    !atmospheric density at seal level on a standard day according to US Standard Atmosphere 1976
    Real(dp), Parameter :: rho_SL = 1225.0002517923321962_dp  ![g/m^3]
    Real(dp), Parameter :: inv_rho_SL = 1._dp / rho_SL

    !Isothermal atmosphere parameter
    !converts temperature to scale height for isothermal atmosphere model
    Real(dp), Parameter :: scale_Height_conv = 0.02927176966650182_dp  ![km/K]

Contains

Function Setup_Atmosphere(setup_file_name,resources_dir,run_file_name,cs_file_name) Result(atm)
    Use Kinds, Only: dp
    Use Global, Only: Rc => R_center
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Output_Message
    Use FileIO_Utilities, Only: Worker_Index
    Use Utilities, Only: Bisection_Search
    Use US_Std_Atm_1976, Only: Zb_1976 => Zb
    Implicit None
    Type(Atmosphere_Type) :: atm
    Character(*), Intent(In) :: setup_file_name
    Character(*), Intent(In) :: resources_dir
    Character(*), Intent(In), Optional :: run_file_name
    Character(:), Allocatable, Intent(InOut), Optional :: cs_file_name
    Integer :: setup_unit,stat
    Character(15) :: atmosphere_model,composition
    Real(dp) :: ref_density,isothermal_temp,Z_top_atm,Z_bot_atm
    Real(dp) :: wind_N,wind_E
    Integer :: n_elements
    Integer, Allocatable :: n_isotopes(:)
    Real(dp), Allocatable :: el_fractions(:)
    Logical, Allocatable :: diatomic(:)
    Integer :: i,j,k
    Character(4), Allocatable :: isotope_names(:)
    Real(dp), Allocatable :: iso_fractions(:)
    Character(:), Allocatable :: f_name
    Real(dp) :: iso_fraction
    Integer :: n_absorption_modes,n_inelastic_lev
    Logical :: has_resonance
    Integer :: n_iso
    Logical :: has_atmosphere

    Real(dp), Parameter :: Zb_1976_extended(0:16) = (/  Zb_1976(0), & !adds the sublayers present in USSA76
                                                     &  Zb_1976(1), &
                                                     &  Zb_1976(2), &
                                                     &  Zb_1976(3), &
                                                     &  Zb_1976(4), &
                                                     &  Zb_1976(5), &
                                                     &  Zb_1976(6), &
                                                     &  Zb_1976(7), &
                                                     &  Zb_1976(8), &
                                                     &      95._dp, &
                                                     &      97._dp, &
                                                     &     100._dp, &
                                                     &  Zb_1976(9), &
                                                     &     115._dp, &
                                                     & Zb_1976(10), &
                                                     &     500._dp, &
                                                     & Zb_1976(11)  /)
    Integer, Parameter :: bZb_1976_extended(1:16) = (/ 0,1,2,3,4,5,6,7,8,8,8,8,9,9,10,10 /) !base indexes for each sublayer

    NameList /AtmosphereList/  atmosphere_model,ref_density,isothermal_temp, & 
                             & Z_top_atm,Z_bot_atm,wind_N,wind_E,composition
    NameList /csSetupList1/ n_elements
    NameList /csSetupList2/ el_fractions,n_isotopes
    NameList /csSetupList3/ isotope_names,diatomic,iso_fractions
    NameList /isoSetupList1/ iso_fraction,n_absorption_modes,n_inelastic_lev,has_resonance

    Open(NEWUNIT = setup_unit , FILE = setup_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Atmospheres: Setup_Atmosphere:  File open error, '//setup_file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Read(setup_unit,NML = AtmosphereList)
    Close(setup_unit)
    has_atmosphere = .TRUE.
    Select Case (atmosphere_model)
        Case ('USstd1976')
            atm%model_index = atm_mod_USstd1976
            atm%discontinuous = .TRUE.
            !find indexes in local layer list for bottom and top of atmosphere
            atm%iZb(1) = Bisection_Search(Z_bot_atm,Zb_1976_extended(1:15),15)
            atm%iZb(2) = Bisection_Search(Z_top_atm,Zb_1976_extended(1:15),15)
            If (Z_top_atm .EQ. Zb_1976_extended(atm%iZb(2)-1)) atm%iZb(2) = atm%iZb(2)-1
            atm%iZb(3) = atm%iZb(2) - atm%iZb(1) + 1
            !fill layer altitudes array
            Allocate(atm%Zb(atm%iZb(1)-1:atm%iZb(2)))
            atm%Zb = Zb_1976_extended(atm%iZb(1)-1:atm%iZb(2))
            atm%Zb(atm%iZb(1)-1) = Z_bot_atm
            atm%Zb(atm%iZb(2)) = Z_top_atm
            !fill base index arrays
            Allocate(atm%bZb(atm%iZb(1)-1:atm%iZb(2)))
            atm%bZb = bZb_1976_extended(atm%iZb(1):atm%iZb(2))
            atm%iZb_map(1) = bZb_1976_extended(atm%iZb(1))
            atm%iZb_map(2) = bZb_1976_extended(atm%iZb(2))
            atm%iZb_map(3) = atm%iZb_map(2) - atm%iZb_map(1) + 1
        Case ('IsoThermal')
            atm%model_index = atm_mod_IsoTherm
            atm%discontinuous = .FALSE.
            atm%iZb = 1
            Allocate(atm%Zb(0:1))
            atm%Zb(0) = Z_bot_atm
            atm%Zb(1) = Z_top_atm
            Allocate(atm%bZb(0:1))
            atm%bZb = 1
            atm%iZb_map = 1
        Case ('None')
            atm%model_index = atm_mod_None
            atm%discontinuous = .FALSE.
            atm%iZb = 1
            Allocate(atm%Zb(0:1))
            atm%Zb(0) = 0._dp
            atm%Zb(1) = 0._dp
            Allocate(atm%bZb(0:1))
            atm%bZb = 1
            atm%iZb_map = 1
            Allocate(atm%Rb(0:Size(atm%Zb)-1))
            atm%Rb = atm%Zb + Rc
            atm%comp_index = atm_comp_None
            atm%ref_density_ratio = 0._dp
            atm%isothermal_temp = 0._dp
            atm%z_top = 0._dp
            atm%z_bot = 0._dp
            atm%R_top = Rc + atm%z_top
            atm%R_bot = Rc + atm%z_bot
            atm%wind_AF = 0._dp
            has_atmosphere = .FALSE.
        Case Default
            Call Output_Message('ERROR:  Atmospheres: Setup_Atmosphere:  Undefined atmosphere model',kill=.TRUE.)
    End Select
    If (has_atmosphere) Then
        Allocate(atm%Rb(0:Size(atm%Zb)-1))
        atm%Rb = atm%Zb + Rc
        Allocate(Character(max_path_len) :: f_name)
        Select Case (composition)
            Case ('All')
                f_name = 'n_CS_setup_All.txt'
                atm%comp_index = atm_comp_ALL
            Case ('N14_O16')
                f_name = 'n_CS_setup_N14_O16.txt'
                atm%comp_index = atm_comp_N14_O16
            Case ('N14_O16_Ar40')
                f_name = 'n_CS_setup_N14_O16_Ar40.txt'
                atm%comp_index = atm_comp_N14_O16_Ar40
            Case ('N14')
                f_name = 'n_CS_setup_N14.txt'
                atm%comp_index = atm_comp_N14
            Case ('N15')
                f_name = 'n_CS_setup_N15.txt'
                atm%comp_index = atm_comp_N15
            Case ('O16')
                f_name = 'n_CS_setup_O16.txt'
                atm%comp_index = atm_comp_O16
            Case ('O17')
                f_name = 'n_CS_setup_O17.txt'
                atm%comp_index = atm_comp_O17
            Case ('O18')
                f_name = 'n_CS_setup_O18.txt'
                atm%comp_index = atm_comp_O18
            Case ('Ar40')
                f_name = 'n_CS_setup_Ar40.txt'
                atm%comp_index = atm_comp_Ar40
            Case Default
                Call Output_Message('ERROR:  Atmospheres: Setup_Atmosphere:  Undefined composition',kill=.TRUE.)
        End Select
        If (Present(cs_file_name)) cs_file_name = f_name
        !read number of elements, isotopes, and fractions from cross sections setup file
        Open( NEWUNIT = setup_unit , FILE = resources_dir//'cs'//slash//'n_cs'//slash//f_name , STATUS = 'OLD' , & 
            & ACTION = 'READ' , IOSTAT = stat)
        If (stat .NE. 0) Call Output_Message( 'ERROR:  Atmospheres: Setup_Atmosphere:  File open error, '//resources_dir// & 
                                            & 'cs'//slash//'n_cs'//slash//f_name//', IOSTAT=',stat,kill=.TRUE.)
        Read(setup_unit,NML = csSetupList1)
        atm%n_el = n_elements
        Allocate(el_fractions(1:n_elements))
        Allocate(n_isotopes(1:n_elements))
        Allocate(atm%iso_ind(1:n_elements+1))
        Read(setup_unit,NML = csSetupList2)
        n_iso = Sum(n_isotopes)
        Allocate(iso_fractions(1:n_iso))
        Allocate(isotope_names(1:n_iso))
        Allocate(diatomic(1:n_iso))
        Allocate(atm%iso_frac(1:n_iso))
        Allocate(atm%diatomic(1:n_iso))
        Allocate(atm%iso_map(1:n_iso))
        atm%iso_ind = 1
        Do i = 2,n_elements+1
            atm%iso_ind(i) = Sum(n_isotopes(1:i-1)) + 1
        End Do
        k = 1
        Do i = 1,n_elements
            Do j = 1,n_isotopes(i)
                atm%iso_map(k) = i
                k = k + 1
            End Do
        End Do
        Read(setup_unit,NML = csSetupList3)
        atm%iso_frac = iso_fractions
        atm%diatomic = diatomic
        Close(setup_unit)
        j = 1
        Do i = 1,n_elements
            atm%iso_frac(j:Sum(n_isotopes(1:i))) = atm%iso_frac(j:Sum(n_isotopes(1:i))) / & 
                                                 & Sum( atm%iso_frac(j:Sum(n_isotopes(1:i))) )
            j = Sum(n_isotopes(1:i)) + 1
        End Do
        atm%ref_density_ratio = ref_density / rho_SL
        atm%isothermal_temp = isothermal_temp
        atm%z_top = Z_top_atm
        atm%z_bot = Z_bot_atm
        atm%R_top = Rc + atm%z_top
        atm%R_bot = Rc + atm%z_bot
        atm%wind_AF = (/ wind_E, & 
                    & wind_N, &
                    & 0._dp /)
        Call Define_EPL_Layers(atm,resources_dir)
    End If
    If (Worker_Index() .EQ. 1) Then
        If (Present(run_file_name)) Then
            Open( NEWUNIT = setup_unit , FILE = run_file_name , STATUS = 'OLD' , ACTION = 'WRITE' , POSITION = 'APPEND' , & 
                & IOSTAT = stat)
            If (stat .NE. 0) Call Output_Message( 'ERROR:  Atmospheres: Setup_Atmosphere:  File open error, '//run_file_name// & 
                                                & ', IOSTAT=',stat,kill=.TRUE.)
            Write(setup_unit,NML = AtmosphereList)
            Write(setup_unit,*)
            Close(setup_unit)
        End If
    End If
End Function Setup_Atmosphere

Subroutine Write_Atmosphere(a,file_name)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Use FileIO_Utilities, Only: half_dash_line
    Implicit None
    Type(Atmosphere_Type), Intent(In) :: a
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat
    Integer :: i
    Real(dp) :: z

    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Atmospheres: Write_Atmosphere:  File open error, '//file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') 'ATMOSPHERE INFORMATION'
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A,ES24.16E3,A,ES24.16E3,A)') '  Extent: ',a%z_bot,' km to ',a%z_top,' km geometric altitude'
    Select Case (a%model_index)
        Case (atm_mod_None)
            Write(unit,'(A)')   '  Atmosphere Model:  NONE' 
            Write(unit,*)
            Write(unit,*)
            Close(unit)
            Return
        Case (atm_mod_IsoTherm)
            Write(unit,'(A,ES24.16E3,A,ES24.16E3,A,ES24.16E3,A)')   '  Atmosphere Model:  Isothermal, ', & 
                                                              & a%Isothermal_temp, & 
                                                                  & ' K,', & 
                                                                  & a%ref_density_ratio*rho_SL, & 
                                                                  & ' g/m^3 @ sea-level,', & 
                                                                  & scale_Height_conv*a%Isothermal_temp, & 
                                                                  & ' km scale height'
            Write(unit,'(2A27)') 'Alt [km]','Density [g/m^3]'
            Write(unit,'(2A27)') '-----------------------','-----------------------'
            i = 0
            Do
                z = a%z_bot + Real(i,dp)*0.1_dp  !100 meter increments
                If (z .GT. a%z_top) Exit
                Write(unit,'(2ES27.16E3)') z,a%rho(z)
                i = i + 1
            End Do
            Write(unit,*)
        Case (atm_mod_USstd1976)
            Write(unit,'(A)') '  Atmosphere Model:  US Std Atm 1976'
            Write(unit,'(3A27)') 'Alt [km]','Temperature [K]','Density [g/m^3]'
            Write(unit,'(3A27)') '-----------------------','-----------------------','-----------------------'
            i = 0
            Do
                z = a%z_bot + Real(i,dp)*0.1_dp  !100 meter increments
                If (z .GT. a%z_top) Exit
                Write(unit,'(3ES27.16E3)') z,a%T(z),a%rho(z)
                i = i + 1
            End Do
            Write(unit,*)
        Case Default
            Call Output_Message('ERROR:  Atmospheres: Write_Atmosphere: Undefined atmosphere model',kill=.TRUE.)
    End Select
    Select Case (a%comp_index)
        Case (atm_comp_ALL)
            Write(unit,'(A)') '  Composition:  N14,N15,O16,O17,O18,Ar40'
        Case (atm_comp_N14_O16)
            Write(unit,'(A)') '  Composition:  N14,O16'
        Case (atm_comp_N14_O16_Ar40)
            Write(unit,'(A)') '  Composition:  N14,O16,Ar40'
        Case (atm_comp_N14)
            Write(unit,'(A)') '  Composition:  N14'
        Case (atm_comp_N15)
            Write(unit,'(A)') '  Composition:  N15'
        Case (atm_comp_O16)
            Write(unit,'(A)') '  Composition:  O16'
        Case (atm_comp_O17)
            Write(unit,'(A)') '  Composition:  O17'
        Case (atm_comp_O18)
            Write(unit,'(A)') '  Composition:  O18'
        Case (atm_comp_Ar40)
            Write(unit,'(A)') '  Composition:  Ar40'
        Case Default
            Call Output_Message('ERROR:  Atmospheres: Write_Atmosphere: Undefined atmosphere composition',kill=.TRUE.)
    End Select
    Write(unit,'(A,ES24.16E3,A)') '  Global Winds:  ',a%wind_AF(1),' [km/s] northerly'
    Write(unit,'(A,ES24.16E3,A)') '                 ',a%wind_AF(2),' [km/s] westerly'
    Write(unit,*)
    Write(unit,*)
    Close(unit)
End Subroutine Write_Atmosphere

Function Atm_Temperature(atm,z,lay) Result(T)
    Use Kinds, Only: dp
    Use US_Std_Atm_1976, Only: T_1976 => T
    Implicit None
    Real(dp) :: T
    Class(Atmosphere_Type), Intent(In) :: atm  !atmosphere model info
    Real(dp), Intent(In) :: z  ![km] geometric height
    Integer, Intent(In), Optional :: lay

    Select Case (atm%model_index)
        Case (atm_mod_USstd1976)
            If (Present(lay)) Then
                T = T_1976(z,layer=atm%bZb(lay)+1)
            Else
                T = T_1976(z,layer_range=atm%iZb_map)
            End If
        Case (atm_mod_IsoTherm)
            T = atm%isothermal_temp
        Case (atm_mod_None)
            T = 0._dp
    End Select
End Function Atm_Temperature

Function Atm_Density(atm,z,lay) Result(rho)
    Use Kinds, Only: dp
    Use US_Std_Atm_1976, Only: rho_1976 => rho
    Implicit None
    Real(dp) :: rho
    Class(Atmosphere_Type), Intent(In) :: atm  !atmosphere model info
    Real(dp), Intent(In) :: z  ![km] geometric height
    Integer, Intent(In), Optional :: lay

    Select Case (atm%model_index)
        Case (atm_mod_USstd1976)
            If (Present(lay)) Then
                rho = rho_1976(z,layer=atm%bZb(lay)+1)
            Else
                rho = rho_1976(z,layer_range=atm%iZb_map)
            End If
        Case (atm_mod_IsoTherm)
            rho = atm%ref_density_ratio * rho_SL * Exp(-z / (scale_Height_conv * atm%isothermal_temp))
        Case (atm_mod_None)
            rho = 0._dp
    End Select
End Function Atm_Density

Subroutine Define_EPL_Layers(atm,resources_dir)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    !number of quadrature points for 6, 9 or 12 digits of precision on STRAIGHT paths
    !                                          USSA76 base index: 0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 9, 9, 10, 10 /)
    Integer, Parameter :: EPL_Quad_n1976_smallZeta_p6(1:16) =  (/ 4, 5, 5, 6, 4, 5, 5, 4, 4, 4, 4, 8, 5, 5, 14,  7 /)
    Integer, Parameter :: EPL_Quad_n1976_smallZeta_p9(1:16) =  (/ 5, 5, 6, 6, 4, 7, 7, 6, 6, 5, 5,14, 7, 6, 23, 10 /)
    Integer, Parameter :: EPL_Quad_n1976_smallZeta_p12(1:16) = (/ 7, 8, 9,10, 6, 9, 9, 8, 7, 6, 7,22,14, 7, 33, 13 /)
    Integer, Parameter :: EPL_Quad_n1976_largeZeta_p6(1:16) =  (/ 4, 4, 5, 5, 3, 5, 4, 4, 3, 3, 3, 6, 4, 4, 22, 15 /)
    Integer, Parameter :: EPL_Quad_n1976_largeZeta_p9(1:16) =  (/ 5, 7, 7, 7, 5, 6, 6, 5, 4, 4, 4,11, 6, 5, 28, 21 /)
    Integer, Parameter :: EPL_Quad_n1976_largeZeta_p12(1:16) = (/ 6, 6, 7, 8, 5, 8, 8, 6, 6, 5, 6,16,10, 6, 36, 28 /)
    Integer, Parameter :: EPL_Quad_nIsoT_smallZeta_p6 =   8
    Integer, Parameter :: EPL_Quad_nIsoT_smallZeta_p9 =  13!TODO
    Integer, Parameter :: EPL_Quad_nIsoT_smallZeta_p12 = 13
    Integer, Parameter :: EPL_Quad_nIsoT_largeZeta_p6 =   9
    Integer, Parameter :: EPL_Quad_nIsoT_largeZeta_p9 =  15!TODO
    Integer, Parameter :: EPL_Quad_nIsoT_largeZeta_p12 = 15
    !TODO Evaluate whether separate q-points numbers are needed for orbital paths... so far, results have always been the same.
    Type(Atmosphere_Type), Intent(InOut) :: atm
    Character(*), Intent(In) :: resources_dir
    Integer :: i,b
    Real(dp) :: dZ
    Real(dp), Allocatable :: z(:)

#   if CHECK_L
        atm%EPL_prec = 12
#   endif
    Allocate(atm%EPL_lay(atm%iZb(1):atm%iZb(2)))
    Do b = atm%iZb(1),atm%iZb(2)
        !get number of quad points
        Select Case (atm%model_index)
            Case (atm_mod_USstd1976)
                !number of quadrature points for 12 digits of precision on STRAIGHT paths
                atm%EPL_lay(b)%nZ = EPL_Quad_n1976_largeZeta_p12(b)
                atm%EPL_lay(b)%nS = EPL_Quad_n1976_smallZeta_p12(b)
                !number of quadrature points for 12 digits of precision on ORBITAL paths
                atm%EPL_lay(b)%nRk = EPL_Quad_n1976_largeZeta_p12(b)
                atm%EPL_lay(b)%nTk = EPL_Quad_n1976_smallZeta_p12(b)
            Case (atm_mod_IsoTherm)
                !number of quadrature points for 12 digits of precision on STRAIGHT paths
                atm%EPL_lay(b)%nZ = EPL_Quad_nIsoT_largeZeta_p12
                atm%EPL_lay(b)%nS = EPL_Quad_nIsoT_smallZeta_p12
                !number of quadrature points for 12 digits of precision on ORBITAL paths
                atm%EPL_lay(b)%nRk = EPL_Quad_nIsoT_largeZeta_p12
                atm%EPL_lay(b)%nTk = EPL_Quad_nIsoT_smallZeta_p12
            Case Default
                Call Output_Message('ERROR:  Amospheres: Define_EPL_Layers: Quad points for this atmosphere model are not &
                                   &implemented',kill=.TRUE.)
        End Select
        !Read in weights and abscissa
        Call Get_Q_points(resources_dir,atm%EPL_lay(b)%nZ,atm%EPL_lay(b)%uZ,atm%EPL_lay(b)%wZ)
        Call Get_Q_points(resources_dir,atm%EPL_lay(b)%nS,atm%EPL_lay(b)%uS,atm%EPL_lay(b)%wS)
        Call Get_Q_points(resources_dir,atm%EPL_lay(b)%nRk,atm%EPL_lay(b)%uRk,atm%EPL_lay(b)%wRk)
        Call Get_Q_points(resources_dir,atm%EPL_lay(b)%nTk,atm%EPL_lay(b)%uTk,atm%EPL_lay(b)%wTk)
        !rescale abscissa to (0,1)
        atm%EPL_lay(b)%uZ = 0.5_dp * atm%EPL_lay(b)%uZ + 0.5_dp
        atm%EPL_lay(b)%uS = 0.5_dp * atm%EPL_lay(b)%uS + 0.5_dp
        atm%EPL_lay(b)%uRk = 0.5_dp * atm%EPL_lay(b)%uRk + 0.5_dp
        atm%EPL_lay(b)%uTk = 0.5_dp * atm%EPL_lay(b)%uTk + 0.5_dp
        !precompute values for full layers
        dZ = atm%Zb(b) - atm%Zb(b-1)
        atm%EPL_lay(b)%A = 0.5_dp * inv_rho_SL * dZ
        Allocate(z(1:atm%EPL_lay(b)%nZ))
        z = atm%EPL_lay(b)%uZ * dZ
        Allocate(atm%EPL_lay(b)%B(1:atm%EPL_lay(b)%nZ))
        atm%EPL_lay(b)%B = 0._dp
        Do i = 1,atm%EPL_lay(b)%nZ
            atm%EPL_lay(b)%B(i) = atm%EPL_lay(b)%wZ(i) * (atm%Rb(b-1) + z(i)) * atm%rho(atm%Zb(b-1) + z(i))
        End Do
        atm%EPL_lay(b)%C = atm%Rb(b-1)
        Allocate(atm%EPL_lay(b)%D(1:atm%EPL_lay(b)%nZ))
        atm%EPL_lay(b)%D = z * (2._dp * atm%Rb(b-1) + z)
        Deallocate(z)
        Allocate(z(1:atm%EPL_lay(b)%nRk))
        z = atm%EPL_lay(b)%uRk * dZ
        Allocate(atm%EPL_lay(b)%Bk(1:atm%EPL_lay(b)%nRk))
        atm%EPL_lay(b)%Bk = 0._dp
        Do i = 1,atm%EPL_lay(b)%nRk
            atm%EPL_lay(b)%Bk(i) = atm%EPL_lay(b)%wRk(i) * atm%rho(atm%Zb(b-1) + z(i))
        End Do
        Allocate(atm%EPL_lay(b)%Ck(1:atm%EPL_lay(b)%nRk))
        atm%EPL_lay(b)%Ck = z + atm%Rb(b-1)
        Deallocate(z)
    End Do
End Subroutine Define_EPL_Layers

Subroutine Get_Q_points(dir,n,a,w)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Character(*), Intent(In) :: dir
    Integer, Intent(In) :: n
    Real(dp), Allocatable, Intent(Out) :: a(:)
    Real(dp), Allocatable, Intent(Out) :: w(:)
    Character(3) :: n_char
    Character(:), Allocatable :: file_name
    Integer :: unit,stat
    Integer :: i

    Allocate(Character(max_path_len) :: file_name)
    file_name = dir
    Write(n_char,'(I3.3)') n
    file_name = dir//'Qpoints'//slash//'Weights_Abscissa-GaussLegendre_n'//n_char//'.txt'
    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Atmospheres: Define_EPL_Layer:  File open error, '//file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Allocate(w(1:n))
    w = 0._dp
    Allocate(a(1:n))
    a = 0._dp
    Do i = 1,n
        Read(unit,*) w(i),a(i)
    End Do
    Close(unit)
End Subroutine Get_Q_points

End Module Atmospheres
