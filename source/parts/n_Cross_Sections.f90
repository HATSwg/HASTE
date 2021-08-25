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
Module n_Cross_Sections

    Use Kinds, Only: dp
    Implicit None
    Private
    Public :: CS_Type
    Public :: sig_Composite
    Public :: sig_Composite_list
    Public :: sig_Resonance
    Public :: Setup_Cross_Sections
    Public :: Write_Cross_Sections

    Type :: sig_Type
        Integer :: n_sig
        Real(dp), Allocatable :: sig(:) ![barns]  has dimension 1:n_sig
        Real(dp), Allocatable :: lnsig(:) ![barns]  has dimension 1:n_sig
        Integer, Allocatable :: E_map(:)  !has dimension 1:n_E_uni, indexes for each value in unified energy grid to indexes in sig
        Integer, Allocatable :: E_key(:)  !has dimension 1:n_sig, indexes for each value in sig to an energy in unified energy list
        Integer :: n_interp_r  !number of interpolation ranges
        Integer, Allocatable :: interp(:,:)  !has dimension 1:n_interp_r and 1:2, dim 1 is sig index up to which to use the 
                                             !interpolation method specified in dim 2
    End Type

    Type :: da_List_Type
        Integer :: n_a  !number of coeffs/points
        Logical :: is_Legendre
        Logical :: is_tab
        Real(dp), Allocatable :: a(:)  !has dim 0:n_a, list of legendre coefficients
        Real(dp), Allocatable :: ua(:,:)  !has dim 1:2,1:n_a, for tabulated cosine PDFs, dim 1 is 1=cosine value, 2=ln(pdf value)
    End Type

    Type :: da_Type
        Integer :: n_da  !number of energies at which da is given (0 indicates isotropic distribution)
        Logical :: is_iso
        Type(da_List_Type), Allocatable :: da(:) ![barns]  has dimension 1:n_da
        Integer, Allocatable :: E_map(:)  !has dimension 1:n_E_uni, indexes for each value in unified energy grid to indexes in da
        Integer, Allocatable :: E_key(:)  !has dimension 1:n_da, indexes for each value in da to an energy in unified energy list
    End Type

    Type :: lev_sig_Type
        Integer :: n_lev
        Real(dp), Allocatable :: Q(:)  ![positive keV]  has dimension 0:n_lev, Q-value for each level
        Integer, Allocatable :: thresh(:)  !has dimension 0:n_lev, index of level threshold energy in unified energy grid
        Type(sig_Type), Allocatable :: sig(:)  !has dimension 0:n_lev, cross sections for each level
        Type(da_Type), Allocatable :: da(:)  !has dimension 0:n_lev, angular cross sections for each level
    End Type

    Type :: abs_sig_Type
        Integer :: n_modes
        Integer, Allocatable :: MT_modes(:)
        Integer, Allocatable :: thresh(:)  !has dimension 1:modes, index of mode threshold energy in unified energy grid
        Type(sig_Type), Allocatable :: sig(:)  !has dim 1:n_modes, cross sections for each absorption mode
    End Type

    Type ::  res_sig_level_Type
        Integer :: n_r  !number of resonance energies in this level
        Real(dp) :: a  ![1E-12 cm] channel radius
        Real(dp) :: ap ![1E-12 cm] scattering radius
        Real(dp), Allocatable :: EraG(:,:)  !Has dimension 1:n_r and 1:b where b depends on the formalism used to define resonances
                                            !For all formalisms, in dim 2, 1=Er[eV] and 2=statistical factor gJ
                                            !For Riech-Moore formalism:  b=6, in dim 2, 3=neutron width, 4=radiation width
                                            !For BW formalisms: b=7, in dim 2, 3=total width, 4=neutron width, 5=radiation width
                                            !For all formalisms, in dim 2, the last 2 columns are shift and penetrability at Er
    End Type

    Type :: res_sig_Type
        Logical :: is_RM  !indicates use of Riech-Moore formalism for resonance representation
        Logical :: is_BW  !indicates use of one of the Breit-Wigner (BW) formalisms for resonance representation
            Logical :: is_SLBW  !indicates use of Single Level BW (SLBW) for resonance representation
            Logical :: is_MLBW  !indicates use of Mulit Level BW (MLBW) for resonance representation
        !Other formalisms (Adler-Adler, Limited-R) are not supported at this time
        Real(dp) :: E_range(1:2)  ![keV]  low and high energy bounds over which resonant contribution is computed
        Real(dp) :: k0  ![1/(1E-12 cm) (when multiplied by Sqrt(eV))] gives k (neutron wave #) when multiplied by Sqrt(E[eV])
        Real(dp) :: spi  !target spin
        Integer :: n_L  !number of levels in which resonances are grouped
        Integer :: s(1:2)  !for RM formalism, minumum and maximum values for 'spin' sums
        Type(res_sig_level_Type), Allocatable :: L(:)   !has dimension 1:n_L, resonance parameters for each level
    End Type

    Type :: CS_Type
        Integer :: n_E_uni  !number of energies in the unified energy list
        Real(dp), Allocatable :: E_uni(:) ![keV]  has dimension 1:n_E_uni, list of energies in the unified energy grid
        Real(dp), Allocatable :: lnE_uni(:) ![keV]  has dimension 1:n_E_uni, list of energies in the unified energy grid
        Integer :: n_iso  !number of isotopes in cross sections structure
        Real(dp), Allocatable :: iso_Fractions(:)  !has dimension 1:n_iso, fractional abundance of each isotope in the total atm
        Real(dp), Allocatable :: An(:)  ![neutron masses] has dimension 1:n_iso, mass (in neutron masses) of isotope nucleus
        Real(dp) :: Mn  ![kg] mean mass (in KILOGRAMS) of nuclei in the atmosphere
        Integer :: n_a_max,n_a_tab_max  !max number of coefficients in angular distribution cross sections
        Type(lev_sig_Type), Allocatable :: lev_cs(:)  !has dimension 1:n_iso, cross sections and angular distributions for each 
                                                      !inelastic level for each isotope
        Type(abs_sig_Type), Allocatable :: abs_cs(:)  !has dimension 1:n_iso, cross sections for each absorption type for each 
                                                      !isotope
        Logical, Allocatable :: has_res_cs(:)  !has dimension 1:n_iso, TRUE indicates resonance cross sections are included in an 
                                               !isotope cross section representation
        Type(res_sig_Type), Allocatable :: res_cs(:)  !has dimension 1:n_iso, resonance cross sections for each isotope
    Contains
        Procedure, Pass :: sig_T_all  !given energy, returns total microscopic cross section for the atmosphere
        Procedure, Pass :: sig_T_iso  !given energy, returns total microscopic cross section for a specific isotope
        Procedure, Pass :: sig_T_all_db
        Procedure, Pass :: sig_T_iso_db
        Procedure, Pass :: sig_S_all  !given energy, returns microscopic scatter cross section for the atmosphere
        Procedure, Pass :: sig_S_iso  !given energy, returns microscopic scatter cross section for a specific isotope
        Procedure, Pass :: sig_S_all_db
        Procedure, Pass :: sig_S_iso_db
        Procedure, Pass :: sig_A_all  !given energy, returns microscopic absorption cross section for the atmosphere
        Procedure, Pass :: sig_A_iso  !given energy, returns microscopic absorption cross section for a specific isotope
        Procedure, Pass :: sig_A_all_db
        Procedure, Pass :: sig_A_iso_db
        Procedure, Pass :: sig_T_A_all  !given energy, returns microscopic total and absorp cross sections for the atmosphere
        Procedure, Pass :: sig_T_A_iso  !given energy, returns microscopic total and absorp cross sections for a specific isotope
        Procedure, Pass :: sig_T_A_db_all
        Procedure, Pass :: sig_T_A_db_iso
        GENERIC :: sig_T_A => sig_T_A_all , sig_T_A_iso , sig_T_A_db_all , sig_T_A_db_iso
        GENERIC :: sig_T => sig_T_all , sig_T_iso , sig_T_all_db , sig_T_iso_db
        GENERIC :: sig_S => sig_S_all , sig_S_iso , sig_S_all_db , sig_S_iso_db
        GENERIC :: sig_A => sig_A_all , sig_A_iso , sig_A_all_db , sig_A_iso_db
    End Type

    !Constants and parameters for Doppler Broadening procedures
    Integer, Parameter :: db_choose_total_sig = 0
    Integer, Parameter :: db_choose_scatter_sig = 1
    Integer, Parameter :: db_choose_absorp_sig = -1
    Integer, Parameter :: Tmax = 21  !limit for number of Romberg extrapolation stages in broadening quadrature
    Real(dp), Parameter :: rTol = 5.E-5_dp  !relative tolerance for convergence of broadening integrals, cross section data has 
                                            !about 5 good digits...

    !MT_disappearance lists the ENDF MT numbers corresponding to neutron interactions WITHOUT a neutron in the exit channel.
    !MTs matching this list are read and stored as absorption modes for the isotope of interest.
    Integer, Parameter :: MT_disappearance(1:15) = (/ 102, &  !n,g
                                                    & 103, &  !n,p
                                                    & 104, &  !n,d
                                                    & 105, &  !n,t
                                                    & 106, &  !n,He3
                                                    & 107, &  !n,a
                                                    & 108, &  !n,2a
                                                    & 109, &  !n,3a
                                                    & 111, &  !n,2p
                                                    & 112, &  !n,p+a
                                                    & 113, &  !n,t+2a
                                                    & 114, &  !n,d+2a
                                                    & 115, &  !n,p+d
                                                    & 116, &  !n,p+t
                                                    & 117  /) !n,d+a
    Integer, Parameter :: MT_inelastic(1:40) = (/ 51, &  !n,n' of the 1-st level
                                                & 52, &  !n,n' of the 2-nd level
                                                & 53, &  !n,n' of the 3-rd level
                                                & 54, &  !n,n' of the 4-th level
                                                & 55, &  !n,n' of the 5-th level
                                                & 56, &  !n,n' of the 6-th level
                                                & 57, &  !n,n' of the 7-th level
                                                & 58, &  !n,n' of the 8-th level
                                                & 59, &  !n,n' of the 9-th level
                                                & 60, &  !n,n' of the 10-th level
                                                & 61, &  !n,n' of the 11-th level
                                                & 62, &  !n,n' of the 12-th level
                                                & 63, &  !n,n' of the 13-th level
                                                & 64, &  !n,n' of the 14-th level
                                                & 65, &  !n,n' of the 15-th level
                                                & 66, &  !n,n' of the 16-th level
                                                & 67, &  !n,n' of the 17-th level
                                                & 68, &  !n,n' of the 18-th level
                                                & 69, &  !n,n' of the 19-th level
                                                & 70, &  !n,n' of the 20-th level
                                                & 71, &  !n,n' of the 21-th level
                                                & 72, &  !n,n' of the 22-th level
                                                & 73, &  !n,n' of the 23-th level
                                                & 74, &  !n,n' of the 24-th level
                                                & 75, &  !n,n' of the 25-th level
                                                & 76, &  !n,n' of the 26-th level
                                                & 77, &  !n,n' of the 27-th level
                                                & 78, &  !n,n' of the 28-th level
                                                & 79, &  !n,n' of the 29-th level
                                                & 80, &  !n,n' of the 30-th level
                                                & 81, &  !n,n' of the 31-th level
                                                & 82, &  !n,n' of the 32-th level
                                                & 83, &  !n,n' of the 33-th level
                                                & 84, &  !n,n' of the 34-th level
                                                & 85, &  !n,n' of the 35-th level
                                                & 86, &  !n,n' of the 36-th level
                                                & 87, &  !n,n' of the 37-th level
                                                & 88, &  !n,n' of the 38-th level
                                                & 89, &  !n,n' of the 39-th level
                                                & 90  /) !n,n' of the 40-th level
    !MT_excluded lists the ENDF MT numbers corresponding to neutron interactions WITH a neutron in the exit channel other than 
    !n,n and n,n'.  MTs matching this list are currently IGNORED by the program.
    Integer, Parameter :: MT_excluded(1:22) = (/  5, &  !interactions not included in any other MT
                                               & 11, &  !n,2n+d
                                               & 16, &  !n,2n
                                               & 17, &  !n,3n
                                               & 22, &  !n,n+a
                                               & 23, &  !n,n+3a
                                               & 24, &  !n,2n+a
                                               & 25, &  !n,3n+a
                                               & 28, &  !n,n+p
                                               & 29, &  !n,n+2a
                                               & 30, &  !n,2n+2a
                                               & 32, &  !n,n+d
                                               & 33, &  !n,n+t
                                               & 34, &  !n,n+He3
                                               & 35, &  !n,n+d+2a
                                               & 36, &  !n,n+t+2a
                                               & 37, &  !n,4n
                                               & 41, &  !n,2n+p
                                               & 42, &  !n,3n+p
                                               & 44, &  !n,n+2p
                                               & 45, &  !n,n+p+a
                                               & 91  /) !n,n' not included in MTs 51-90

Contains

Function Setup_Cross_Sections( resources_directory, & 
                             & cs_setup_file,       & 
                             & cs_summary_file,     & 
                             & elastic_only,        & 
                             & aniso_dist,          & 
                             & E_min,E_max,         & 
                             & verbosity            ) Result(CS)
    Use Kinds, Only: dp
    Use Sorting, Only: Union_Sort
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Output_Message
    Use FileIO_Utilities, Only: half_dash_line
    Use FileIO_Utilities, Only: to_str
    Use Global, Only: neutron_mass
    Implicit None
    Type(CS_Type) :: CS
    Character(*), Intent(In) :: resources_directory
    Character(*), Intent(In) :: cs_setup_file
    Character(*), Intent(In) :: cs_summary_file
    Logical, Intent(In) :: elastic_only
    Logical, Intent(In) :: aniso_dist
    Real(dp), Intent(In) :: E_min,E_max
    Integer, Intent(In), Optional :: verbosity  ! >0 causes comprehensive output of cross section data to be generated during the 
                                                ! setup process.  Includes summary of ENDF data files read in and detailed output 
                                                ! of interpreted cross section data.  This output is placed in the folder specified 
                                                ! by resources_directory.
                                                !    1 = basic information on the content of the ENDF data file is output
                                                !    2 = (includes 1) ENDF data (as stored) is output to file
                                                !    3 = (includes 1 & 2) Cross sections (total, scatter, absorption) for each 
                                                !        isotope are written to file for plotting
                                                !    4 = (includes 1,2, & 3) Cross sections (total, scatter, absorption) doppler 
                                                !        broadened to 300K for each isotope are written to file for plotting

    Integer :: n_elements
    Character(:), Allocatable :: file_name_start,ENDF_file_name
    Integer, Allocatable :: n_isotopes(:)
    Character(4), Allocatable :: isotope_names(:)
    Real(dp), Allocatable :: el_fractions(:),iso_fractions(:)
    Integer :: n_energies
    Integer :: n_p,n_p_2,n_r,n_start
    Integer :: n_a,n_a_lines
    Integer :: h,i,j,k
    Real(dp) :: Q_scratch
    Real(dp) :: An_scratch
    Real(dp), Allocatable :: E_uni_scratch(:)
    Real(dp), Allocatable :: E_scratch(:)
    Real(dp), Allocatable :: sig_scratch(:)
    Integer, Allocatable :: Interp_scratch(:,:)
    Type(da_List_type), Allocatable :: Ang_Dist_scratch(:)
    Integer :: setup_unit,ENDF_unit,v_unit,stat
    Integer, Allocatable :: n_abs_modes(:),n_inel_lev(:)
    Integer, Allocatable :: abs_modes(:,:),inel_levs(:,:),inel_da_mf(:,:)
    Real(dp), Allocatable :: abs_thresh(:,:),inel_thresh(:,:)
    Logical, Allocatable :: diatomic(:)
    Integer :: LTT,LRP,LIP,LAW
    Integer :: MT,MF,line_num
    Real(dp) :: ZAP,AWP
    Character(80) :: trash_c
    Real(dp) :: trash_r
    Logical :: first_time
    Logical :: v_basic,v_full,v_burn0,v_burn300
    Character(15) :: v_string
    Integer, Parameter :: inter_pts = 5
    Real(dp) :: En,sT,sS,sA
    Real(dp) :: T_broad
    Logical :: skip_section
    Real(dp) :: thresh

    NameList /csSetupList1/ n_elements
    NameList /csSetupList2/ el_fractions,n_isotopes
    NameList /csSetupList3/ isotope_names,iso_fractions,diatomic

    v_basic = .FALSE. !default is NOT verbose
    v_full = .FALSE. !default is NOT verbose
    v_burn0 = .FALSE. !default is NOT verbose
    v_burn300 = .FALSE. !default is NOT verbose
    If (Present(verbosity)) Then
        If (verbosity .GE. 1) v_basic = .TRUE.
        If (verbosity .GE. 2) v_full = .TRUE.
        If (verbosity .GE. 3) v_burn0 = .TRUE.
        If (verbosity .GE. 4) v_burn300 = .TRUE.
    End If
    !read namelists from cross sections setup file
    Allocate(Character(max_path_len) :: file_name_start)
    file_name_start = resources_directory//'cs'//slash//'n_cs'//slash
    Open(NEWUNIT = setup_unit , FILE = file_name_start//cs_setup_file , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  File open error, ' & 
                                        & //file_name_start//cs_setup_file//', IOSTAT=',stat,kill=.TRUE. )
    Read(setup_unit,NML = csSetupList1)
    Allocate(el_fractions(1:n_elements))
    Allocate(n_isotopes(1:n_elements))
    Read(setup_unit,NML = csSetupList2)
    el_fractions = el_fractions / Sum(el_fractions) !make sure el_fractions is normalized to 1
    CS%n_iso = Sum(n_isotopes)
    Allocate(isotope_names(1:CS%n_iso))
    Allocate(iso_fractions(1:CS%n_iso))
    Allocate(diatomic(1:CS%n_iso))
    Read(setup_unit,NML = csSetupList3)
    j = 1
    Do i = 1,n_elements
        !make sure iso_fractions are normalized to 1 within each isotope
        iso_fractions(j:Sum(n_isotopes(1:i))) =  iso_fractions(j:Sum(n_isotopes(1:i))) / Sum(iso_fractions(j:Sum(n_isotopes(1:i))))
        !combine element and isotope fractions for total atmospheric fraction of each isotope
        iso_fractions(j:Sum(n_isotopes(1:i))) = iso_fractions(j:Sum(n_isotopes(1:i))) * el_fractions(i)
        j = Sum(n_isotopes(1:i)) + 1
    End Do
    Close(setup_unit)
    !Initialize the CS structure
    Allocate(CS%iso_fractions(1:CS%n_iso))
    CS%iso_fractions = iso_fractions
    Allocate(CS%An(1:CS%n_iso))
    CS%An = -1._dp
    Allocate(CS%lev_cs(1:CS%n_iso))
    Allocate(CS%abs_cs(1:CS%n_iso))
    Allocate(CS%has_res_cs(1:CS%n_iso))
    CS%has_res_cs = .FALSE.
    Allocate(CS%res_cs(1:CS%n_iso))
    If (v_basic) Then  !create file for verbose output
        Open(NEWUNIT = v_unit , FILE = cs_summary_file , STATUS = 'REPLACE' , ACTION = 'WRITE' , IOSTAT = stat)
        If (stat .NE. 0) Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  File open error, ' & 
                                            & //cs_summary_file//', IOSTAT=',stat,kill=.TRUE. )
    End If
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  CREATE A UNIFIED ENERGY LIST FOR THE CROSS SECTION DATA
    !count the number of energies (including duplicates) for MF 3 and 4 (interaction cross sections and angular distributions)
    !also count number of disappearance (absorption) modes (MF=3 and MTs appearing in MT_disappearance)
    !also count number of inelastic levels (MF=3 and MTs 51-90)
    !read in ALL (including duplicates) the energies and generate a unified list
    Allocate(Character(max_path_len) :: ENDF_file_name)
    first_time = .TRUE.
    n_start = 0
    Do h = 1,2
    !FIRST TIME: Count interactions and total number of energies for unified list
    !SECOND TIME: Read in ALL (including duplicates) the energies and generate unified list
        If (first_time) Then  !FIRST TIME
            If (v_basic) Then
                Write(v_unit,'(A)') half_dash_line
                Write(v_unit,'(A)') 'ENDF TAPE FILE INVENTORY'
                Write(v_unit,'(A)') half_dash_line
                Write(v_unit,*)
            End If
            n_energies = 0
            Allocate(n_abs_modes(1:CS%n_iso))
            n_abs_modes = 0
            Allocate(n_inel_lev(1:CS%n_iso))
            n_inel_lev = 0
        Else  !SECOND TIME
            n_start = 0
            !Allocate scratch arrays for ordering absorption modes and inelastic levels
            Allocate(abs_modes(1:MaxVal(n_abs_modes),1:CS%n_iso))
            abs_modes = -1
            Allocate(abs_thresh(1:MaxVal(n_abs_modes),1:CS%n_iso))
            abs_thresh = Huge(abs_thresh)
            Allocate(inel_levs(1:MaxVal(n_inel_lev),1:CS%n_iso))
            inel_levs = -1
            Allocate(inel_da_mf(1:MaxVal(n_inel_lev),1:CS%n_iso))
            inel_da_mf = -1
            Allocate(inel_thresh(1:MaxVal(n_inel_lev),1:CS%n_iso))
            inel_thresh = Huge(abs_thresh)
            !Allocate scratch arrays for energy
            Allocate(E_uni_scratch(1:n_energies))
            E_uni_scratch = -1._dp
        End If
        Do i = 1,CS%n_iso
            !create file name string and open the ENDF tape file for this isotope
            ENDF_file_name = file_name_start//Trim(isotope_names(i))//'.txt'
            Open(NEWUNIT = ENDF_unit , FILE = ENDF_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
            If (stat .NE. 0) Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  File open error, ' & 
                                                & //ENDF_file_name//', IOSTAT=',stat,kill=.TRUE. )
            !count energies in absorption, elastic, and inelastic files (MFs 3, 4, and 6)
            If (first_time .AND. v_basic) v_string = Trim(isotope_names(i))//' ENDF file:'  !FIRST TIME
            DO_SECTIONS: Do
                !check next section type
                Read(ENDF_unit,'(A71,I1,I3,I5)') trash_c,MF,MT,line_num
                If (MF .EQ. 0) Then  !this indicates a change in MF type, need to look for next section or end of file
                    NEXT_MF: Do
                        Read(ENDF_unit,'(A71,I1,I3,I5)',IOSTAT=stat) trash_c,MF,MT,line_num
                        If (stat .LT. 0) Exit DO_SECTIONS !end of file
                        If (MF .GT. 0) Exit NEXT_MF
                    End Do NEXT_MF
                End If
                If (first_time .AND. v_basic) Then  !FIRST TIME
                    Write(v_unit,'(A15,A4,I0,A4,I0,A1)',ADVANCE='NO') v_string,' MF=',MF,' MT=',MT,' '
                    v_string = ''
                End If
                If (MF.EQ.3 .OR. MF.EQ.4 .OR. MF.EQ.6) Then
                    skip_section = .FALSE.
                    If ( Any(MT_excluded .EQ. MT) .OR. & 
                       & (MF.EQ.6 .AND. Any(MT_disappearance.EQ.MT)) ) Then !this type/file is excluded, advance past it
                        If (first_time .AND. v_basic) Write(v_unit,'(A)',ADVANCE='YES') ' - excluded'  !FIRST TIME
                        skip_section = .TRUE.
                    Else If (MF .EQ. 3) Then !also need to check if range of energies has values between (E_min,E_max)
                        Read(ENDF_unit,*)
                        Read(ENDF_unit,*)
                        Read(ENDF_unit,'(E11.6E1)') thresh
                        Backspace(ENDF_unit)
                        Backspace(ENDF_unit)
                        Backspace(ENDF_unit)
                        If (thresh/1000._dp .GE. E_max) Then !energy out of range
                            If (first_time .AND. v_basic) Write(v_unit,'(A)',ADVANCE='YES') ' - E out of range'  !FIRST TIME
                            skip_section = .TRUE.
                        End If
                    End If
                    If (skip_section) Then
                        !advance to end of section
                        Call Find_MFMT_end(ENDF_unit)
                        !cycle to start next section
                        Cycle DO_SECTIONS
                    Else
                        If (first_time .AND. v_basic) Write(v_unit,'(A)',ADVANCE='NO') ' - counted'  !FIRST TIME
                    End If
                    !get number of energies in this MF 3, 4, or 6 section
                    Select Case (MF)
                        Case (3)
                            If (first_time) Then  !FIRST TIME
                                If (Any(MT_disappearance .EQ. MT)) Then
                                    n_abs_modes(i) = n_abs_modes(i) + 1
                                    If (v_basic) Write(v_unit,'(A)',ADVANCE='NO') ', disappearance'
                                Else If (Any(MT_inelastic .EQ. MT)) Then
                                    n_inel_lev(i) = n_inel_lev(i) + 1
                                    If (v_basic) Write(v_unit,'(A)',ADVANCE='NO') ', inelastic'
                                End If
                            End If
                            Read(ENDF_unit,'(A55,I11)') trash_c, n_p
                            If (first_time) Then  !FIRST TIME
                                If (v_basic) Write(v_unit,'(A,I0,A)',ADVANCE='YES') ', ',n_p,' E-points'
                            Else  !SECOND TIME
                                Read(ENDF_unit,*)
                                Do j = 1,n_p,3
                                    !Each line in the file has 3 pairs of (eV,barns)
                                    If (n_p-j .GT. 1) Then
                                        Read(ENDF_unit,'(5E11.6E1)') E_uni_scratch(n_start+j), &
                                                                   & trash_r, &
                                                                   & E_uni_scratch(n_start+j+1), &
                                                                   & trash_r, &
                                                                   & E_uni_scratch(n_start+j+2)
                                    Else If (n_p-j .EQ. 1) Then
                                        Read(ENDF_unit,'(3E11.6E1)') E_uni_scratch(n_start+j), &
                                                                   & trash_r, &
                                                                   & E_uni_scratch(n_start+j+1)
                                    Else
                                        Read(ENDF_unit,'(1E11.6E1)') E_uni_scratch(n_start+j)
                                    End If
                                    If (j .EQ. 1) Then !check if this is an absorption mode
                                        If (Any(MT_disappearance .EQ. MT)) Then  !this is an absorption mode
                                            !add it to the list of absorption modes, inserting in order of acending threshold energy
                                            Do k = 1,n_abs_modes(i)
                                                If (E_uni_scratch(n_start+j) .LT. abs_thresh(k,i)) Then 
                                                !this k is where this mode goes in the list
                                                    If (k .LT. n_abs_modes(i)) Then 
                                                    !make room by shifting the rest of the list down
                                                        abs_thresh(k+1:n_abs_modes(i),i) = abs_thresh(k:n_abs_modes(i)-1,i)
                                                        abs_modes(k+1:n_abs_modes(i),i) = abs_modes(k:n_abs_modes(i)-1,i)
                                                    End If
                                                    !insert
                                                    abs_thresh(k,i) = E_uni_scratch(n_start+j)
                                                    abs_modes(k,i) = MT
                                                    Exit
                                                End If
                                            End Do
                                        Else If (Any(MT_inelastic .EQ. MT)) Then  !this is an inelastic mode
                                            !add it to the list of inelastic modes, inserting in order of acending threshold energy
                                            Do k = 1,n_inel_lev(i)
                                                If (E_uni_scratch(n_start+j) .LT. inel_thresh(k,i)) Then 
                                                !this k is where this mode goes in the list
                                                    If (k .LT. n_inel_lev(i)) Then 
                                                    !make room by shifting the rest of the list down
                                                        inel_thresh(k+1:n_inel_lev(i),i) = inel_thresh(k:n_inel_lev(i)-1,i)
                                                        inel_levs(k+1:n_inel_lev(i),i) = inel_levs(k:n_inel_lev(i)-1,i)
                                                    End If
                                                    !insert
                                                    inel_thresh(k,i) = E_uni_scratch(n_start+j)
                                                    inel_levs(k,i) = MT
                                                    Exit
                                                End If
                                            End Do
                                        End If
                                    End If
                                End Do
                            End If
                        Case (4)
                            If (first_time .AND. v_basic) Then  !FIRST TIME
                                If (Any(MT_disappearance.EQ.MT)) Then
                                    Write(v_unit,'(A)',ADVANCE='NO') ', disappearance'
                                Else If (Any(MT_inelastic.EQ.MT)) Then
                                    Write(v_unit,'(A)',ADVANCE='NO') ', inelastic'
                                End If
                            End If
                            !need to read the first line again to get LTT
                            Backspace(ENDF_unit)
                            Read(ENDF_unit,'(A33,I11)') trash_c, LTT
                            Read(ENDF_unit,*)
                            Read(ENDF_unit,'(A55,I11)') trash_c,n_p
                            Read(ENDF_unit,*)
                            If (first_time) Then  !FIRST TIME
                                If (LTT .EQ. 0) Then !distribution is assumed isotropic over all energies
                                    n_p = 0
                                    n_p_2 = 0
                                    If (v_basic) Write(v_unit,'(A)',ADVANCE='YES') ', isotropic'
                                Else If (LTT .EQ. 3) Then !there is a second range of energies later in the section
                                    !advance in the file to the end of the Legendre section
                                    Do j = 1,n_p
                                        !The first line in each energy contains the energy in eV in the second position and the 
                                        !number of Legendre coefficients in the 5th position
                                        Read(ENDF_unit,'(A44,I11)') trash_c, n_a
                                        n_a_lines = (n_a / 6)  !integer divide
                                        If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
                                        !advance to the next energy
                                        Do k = 1,n_a_lines
                                            Read(ENDF_unit,*)
                                        End Do
                                    End Do
                                    Read(ENDF_unit,*)
                                    !first entry of next line is number of additional energy points
                                    Read(ENDF_unit,'(I11)') n_p_2
                                    If (v_basic) Write(v_unit,'(A,I0,A,I0,A)',ADVANCE='YES') ', ',n_p,' + ',n_p_2,' E-points'
                                Else
                                    n_p_2 = 0
                                    If (v_basic) Write(v_unit,'(A,I0,A)',ADVANCE='YES') ', ',n_p,' E-points'
                                End If
                                n_p = n_p + n_p_2
                            Else  !SECOND TIME
                                ! If (LTT .EQ. 0) Then !da is isotropic, no energies to read
                                If (LTT .EQ. 1) Then  !da is lists of legendre coeffs
                                    Do j = 1,n_p
                                        !The first line in each energy contains the energy in eV in the second position and the 
                                        !number of Legendre coefficients in the 5th position
                                        Read(ENDF_unit,'(A11,E11.6E1,A22,I11)') trash_c, E_uni_scratch(n_start+j), trash_c, n_a
                                        n_a_lines = (n_a / 6)  !integer divide
                                        If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
                                        !advance to the next energy
                                        Do k = 1,n_a_lines
                                            Read(ENDF_unit,*)
                                        End Do
                                    End Do
                                Else If (LTT .EQ. 2) Then  !da is tabulated
                                    Do j = 1,n_p
                                        !The first line in each energy contains the energy in eV in the second position
                                        Read(ENDF_unit,'(A11,E11.6E1)') trash_c, E_uni_scratch(n_start+j)
                                        !the next line contains the number of tabulation points in the first position
                                        Read(ENDF_unit,'(I11)') n_a
                                        n_a_lines = (n_a / 3)  !integer divide
                                        If (Mod(n_a,3) .GT. 0) n_a_lines = n_a_lines + 1
                                        !advance to the next energy
                                        Do k = 1,n_a_lines
                                            Read(ENDF_unit,*)
                                        End Do
                                    End Do
                                Else If (LTT .EQ. 3) Then  !da is tabulated for high energies but legendre for low energies
                                    !Read in low energy Legendre points
                                    Do j = 1,n_p
                                        !The first line in each energy contains the energy in eV in the second position and the 
                                        !number of Legendre coefficients in the 5th position
                                        Read(ENDF_unit,'(A11,E11.6E1,A22,I11)') trash_c, E_uni_scratch(n_start+j), trash_c, n_a
                                        n_a_lines = (n_a / 6)  !integer divide
                                        If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
                                        !advance to the next energy
                                        Do k = 1,n_a_lines
                                            Read(ENDF_unit,*)
                                        End Do
                                    End Do
                                    Read(ENDF_unit,*)
                                    !first entry of next line is number of additional energy points
                                    Read (ENDF_unit,'(I11)') n_p_2
                                    !Read in high energy tabulated cosine points
                                    Do j = n_p+1,n_p+n_p_2
                                        !The first line in each energy contains the energy in eV in the second position
                                        Read(ENDF_unit,'(A11,E11.6E1)') trash_c, E_uni_scratch(n_start+j)
                                        !the next line contains the number of tabulation points in the first position
                                        Read(ENDF_unit,'(I11)') n_a
                                        n_a_lines = (n_a / 3)  !integer divide
                                        If (Mod(n_a,3) .GT. 0) n_a_lines = n_a_lines + 1
                                        !advance to the next energy
                                        Do k = 1,n_a_lines
                                            Read(ENDF_unit,*)
                                        End Do
                                    End Do
                                    !update total n_p
                                    n_p = n_p + n_p_2
                                End If
                            End If
                        Case(6)
                            If (first_time .AND. v_basic) Then  !FIRST TIME
                                If (Any(MT_disappearance.EQ.MT)) Then
                                    Write(v_unit,'(A)',ADVANCE='NO') ', disappearance'
                                Else If (Any(MT_inelastic.EQ.MT)) Then
                                    Write(v_unit,'(A)',ADVANCE='NO') ', inelastic'
                                End If
                            End If
                            !find the start of the section describing the angular distribution of the neutron on the exit channel
                            Do
                                Read(ENDF_unit,'(2E11.6E1,2I11)') ZAP,AWP,LIP,LAW
                                If (ZAP.EQ.1._dp .AND. AWP.EQ.1._dp .AND. LAW.EQ.2) Exit !this is first line in the correct section
                            End Do
                            !the next 2 lines are always discarded
                            Read(ENDF_unit,*)
                            Read(ENDF_unit,*)
                            Read(ENDF_unit,'(A55,I11)') trash_c,n_p
                            Read(ENDF_unit,*)
                            If (first_time) Then  !FIRST TIME
                                n_p_2 = 0
                                If (v_basic) Write(v_unit,'(A,I0,A)',ADVANCE='YES') ', ',n_p,' E-points'
                                n_p = n_p + n_p_2
                            Else  !SECOND TIME
                                Do j = 1,n_p
                                    !The first line in each energy contains the energy in eV in the second position and the 
                                    !number of Legendre coefficients in the 5th position
                                    Read(ENDF_unit,'(A11,E11.6E1,A22,I11)') trash_c, E_uni_scratch(n_start+j), trash_c, n_a
                                    n_a_lines = (n_a / 6)  !integer divide
                                    If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
                                    !advance to the next energy
                                    Do k = 1,n_a_lines
                                        Read(ENDF_unit,*)
                                    End Do
                                End Do
                            End If
                    End Select
                    If (first_time) Then  !FIRST TIME
                        n_energies = n_energies + n_p
                    Else  !SECOND TIME
                        n_start = n_start + n_p
                    End If
                Else If (MF.EQ.1 .AND. MT.EQ.451) Then
                    If (first_time .AND. v_basic) Write(v_unit,'(A)',ADVANCE='YES') ' - uncounted, resonances'  !FIRST TIME
                Else
                    If (first_time .AND. v_basic) Write(v_unit,'(A)',ADVANCE='YES') ' - uncounted'  !FIRST TIME
                End If
                !advance to end of section
                Call Find_MFMT_end(ENDF_unit)
            End Do DO_SECTIONS
            Close(ENDF_unit)
        End Do
        If (first_time) Then  !FIRST TIME
            first_time = .FALSE.
        Else  !SECOND TIME
            If (v_basic) Then
                Write(v_unit,*)
                Do i = 1,CS%n_iso
                    Write(v_unit,'(A,I0,A,I0,A)') Trim(isotope_names(i))//' ENDF tape contains ',n_abs_modes(i), & 
                                                & ' absorption modes and ',n_inel_lev(i),' inelastic levels'
                    Write(v_unit,'(A10,A16,A11,A16)') ' Abs MT ','  E-thresh [keV]',' Inel MT ','  E-thresh [keV]'
                    Write(v_unit,'(A10,A16,A11,A16)') '--------','  --------------','---------','  --------------'
                    Do j = 1,Max(n_abs_modes(i),n_inel_lev(i))
                        If (j .LE. n_abs_modes(i)) Then
                            Write(v_unit,'(I10,ES16.7E2)',ADVANCE='NO') abs_modes(j,i),abs_thresh(j,i)/1000._dp
                        Else
                            Write(v_unit,'(A26)',ADVANCE='NO') ''
                        End If
                        If (j .LE. n_inel_lev(i)) Then
                            Write(v_unit,'(I11,ES16.7E2)') inel_levs(j,i),inel_thresh(j,i)/1000._dp
                        Else
                            Write(v_unit,*)
                        End If
                    End Do
                End Do
                Write(v_unit,*)
                Write(v_unit,'(I0,A)') n_energies,' total energies to map on unified list.'
            End If
            n_p = n_energies !stash n_energies before it gets overwritten
            E_uni_scratch = E_uni_scratch / 1000._dp  !convert to kev
            !E_uni_scratch is now a HUGE list of all the energies in all the files to be used, sort and eliminate duplicates
            Call Union_Sort(E_uni_scratch,n_energies,E_min,E_max)
            If (n_energies .GT. Huge(CS%n_E_uni)) Call Output_Message('ERROR:  Cross_Sections: Setup_Cross_Sections:  Length of &
                                                                      &unified energy grid exceeds available index',kill=.TRUE.)
            !Allocate and fill the unified energy list
            CS%n_E_uni = n_energies
            Allocate(CS%E_uni(1:n_energies))
            CS%E_uni = E_uni_scratch(1:n_energies)
            Deallocate(E_uni_scratch)
            Allocate(CS%lnE_uni(1:n_energies))
            CS%lnE_uni = Log(CS%E_uni)
            If (v_full) Then
                Write(v_unit,*)
                Write(v_unit,'(A)') half_dash_line
                Write(v_unit,'(A)') 'UNIFIED ENERGY LIST'
                Write(v_unit,'(A)') half_dash_line
                Write(v_unit,'(I0,A,I0,A,F0.2,A)') CS%n_E_uni,' unique points in unified grid out of ',n_p,', (', & 
                                                 & 100._dp*Real(CS%n_E_uni,dp)/Real(n_p,dp),'%)'
                Write(v_unit,'(A9,2A16)') ' Index ',' E [keV]      ',' ln(E)        '
                Write(v_unit,'(A9,2A16)') '-------','--------------','--------------'
                Do i = 1,n_energies
                    Write(v_unit,'(I9,2ES16.7E2)') i,CS%E_uni(i),CS%lnE_uni(i)
                End Do
            End If
        End If
    End Do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  READ IN INTERACTION CROSS SECTIONS AND ANGULAR DISTRIBUTIONS
    !Now read each file (yes, again) into its place in the cross section structure, constructing the integer maps and keys as we go
    If (v_full) Then
        Write(v_unit,*)
        Write(v_unit,'(A)') half_dash_line
        Write(v_unit,'(A)') 'CROSS SECTION DATA STORED AND INDEXED'
        Write(v_unit,'(A)') half_dash_line
        Write(v_unit,*)
    End If
    CS%n_a_max = 0
    CS%n_a_tab_max = 0
    Do i = 1,CS%n_iso
        !create file name string and open the ENDF tape for this isotope
        ENDF_file_name = file_name_start//Trim(isotope_names(i))//'.txt'
        Open(NEWUNIT = ENDF_unit , FILE = ENDF_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
        If (stat .NE. 0) Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  File open error, ' & 
                                            & //ENDF_file_name//', IOSTAT=',stat,kill=.TRUE. )
        !read in absorption, resonance, elastic, and inelastic interaction cross sections and angular distributions
        CS%abs_cs(i)%n_modes = n_abs_modes(i)
        Allocate(CS%abs_cs(i)%MT_modes(1:n_abs_modes(i)))
        CS%abs_cs(i)%MT_modes = -1
        Allocate(CS%abs_cs(i)%thresh(1:n_abs_modes(i)))
        CS%abs_cs(i)%thresh = -1
        Allocate(CS%abs_cs(i)%sig(1:n_abs_modes(i)))
        !!  ABSORPTION INTERACTION CROSS SECTIONS
        Do j = 1,n_abs_modes(i)
            !Find this interaction in the ENDF tape (MF=3, MT=abs_modes(j,i))
            MF = 3
            MT = abs_modes(j,i)
            If (.NOT.Find_MFMT(ENDF_unit,MF,MT)) Then
                Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  Section not found, MF='// & 
                                   & Trim(Adjustl(to_str(MF)))//', MT='//Trim(Adjustl(to_str(MT))) , kill=.TRUE. )
            End If
            !the next read statement on ENDF_unit will read the first line of MF=3, MT=abs_modes(j,i)
            Call Read_sig_sect(ENDF_unit,Q_scratch,An_scratch,E_scratch,sig_scratch,Interp_scratch,n_p,n_r)
            If (Trim_CS_for_E(n_p,E_scratch,sig_scratch,n_r,Interp_scratch,E_min,E_max)) Then
                CS%abs_cs(i)%MT_modes(j) = MT
                Call Map_and_Store_CS( CS%n_E_uni,            & 
                                     & CS%E_uni,              & 
                                     & n_p,                   & 
                                     & E_scratch,             & 
                                     & sig_scratch,           & 
                                     & n_r,                   & 
                                     & Interp_Scratch,        & 
                                     & CS%abs_cs(i)%sig(j),   & 
                                     & CS%abs_cs(i)%thresh(j) )
            End If
            Deallocate(E_scratch,sig_scratch,Interp_scratch)
            If (v_full) Then  !write the stored values for this absorption mode
                Write(v_unit,'(A,I0,A,I0,A)') Trim(isotope_names(i))//' MF=',3,', MT=',abs_modes(j,i),' (absorption)'
                Call Write_stored_sig(v_unit,CS%abs_cs(i)%sig(j),CS%n_E_uni,CS%E_uni)
            End If
        End Do
        !N2H For verbose output, write a summary of absorption modes for this isotope
        !!  SCATTERING INTERACTION CROSS SECTIONS
        CS%lev_cs(i)%n_lev = n_inel_lev(i)
        Allocate(CS%lev_cs(i)%Q(0:n_inel_lev(i)))
        Allocate(CS%lev_cs(i)%thresh(0:n_inel_lev(i)))
        Allocate(CS%lev_cs(i)%sig(0:n_inel_lev(i)))
        If (aniso_dist) Allocate(CS%lev_cs(i)%da(0:n_inel_lev(i)))
        !check if resonance cross sections are present
        MF = 1
        MT = 451
        If (.NOT.Find_MFMT(ENDF_unit,MF,MT)) Then
            Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  Section not found, MF='// & 
                               & Trim(Adjustl(to_str(MF)))//', MT='//Trim(Adjustl(to_str(MT))) , kill=.TRUE. )
        End If
        !the next read statement on ENDF_unit will read the first line of MF=1, MT=451
        Read(ENDF_unit,'(A22,I11)') trash_c,LRP
        If (LRP .EQ. 1) Then  !resonance parameters are included in the ENDF tape
            CS%has_res_cs(i) = .TRUE.
            !Find this interaction in the ENDF tape (MF=2, MT=151)
            MF = 2
            MT = 151
            If (.NOT.Find_MFMT(ENDF_unit,MF,MT)) Then
                Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  Section not found, MF='// & 
                                   & Trim(Adjustl(to_str(MF)))//', MT='//Trim(Adjustl(to_str(MT))) , kill=.TRUE. )
            End If
            !the next read statement on ENDF_unit will read the first line of MF=2, MT=151
            Call Read_res_sect(ENDF_unit,CS%res_cs(i))
            If (v_full) Then  !write the stored values for this resonance representation
                Write(v_unit,'(A,I0,A,I0,A)') Trim(isotope_names(i))//' MF=',2,', MT=',151,' (resonance)'
                Call Write_stored_res(v_unit,CS%res_cs(i))
            End If
        Else  !no resonance parameters, all interaction cross sections are tabulated in MF=3
            CS%has_res_cs(i) = .FALSE.
            If (v_full) Then  !write the stored values for this resonance representation
                Write(v_unit,'(A)') Trim(isotope_names(i))//' has no ENDF resonance data.'
                Write(v_unit,*)
            End If
        End If
        !!  ELASTIC SCATTERING INTERACTION CROSS SECTION
        !Find this interaction in the ENDF tape (MF=3, MT=2)
        MF = 3
        MT = 2
        If (.NOT.Find_MFMT(ENDF_unit,MF,MT)) Then
            Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  Section not found, MF='// & 
                               & Trim(Adjustl(to_str(MF)))//', MT='//Trim(Adjustl(to_str(MT))) , kill=.TRUE. )
        End If
        !the next read statement on ENDF_unit will read the first line of MF=3, MT=2
        Call Read_sig_sect(ENDF_unit,Q_scratch,CS%An(i),E_scratch,sig_scratch,Interp_scratch,n_p,n_r)
        If (Trim_CS_for_E(n_p,E_scratch,sig_scratch,n_r,Interp_scratch,E_min,E_max)) Then
            Call Map_and_Store_CS( CS%n_E_uni,            & 
                                 & CS%E_uni,              & 
                                 & n_p,                   & 
                                 & E_scratch,             & 
                                 & sig_scratch,           & 
                                 & n_r,                   & 
                                 & Interp_Scratch,        & 
                                 & CS%lev_cs(i)%sig(0),   & 
                                 & CS%lev_cs(i)%thresh(0) )
        End If
        Deallocate(E_scratch,sig_scratch,Interp_scratch)
        If (v_full) Then  !write the stored values for elastic scatter interaction cross section
            Write(v_unit,'(A,I0,A,I0,A)') Trim(isotope_names(i))//' MF=',3,', MT=',2,' (sig, elastic)'
            Call Write_stored_sig(v_unit,CS%lev_cs(i)%sig(0),CS%n_E_uni,CS%E_uni)
        End If
        !!  ELASTIC SCATTERING ANGULAR DISTRIBUTION
        If (aniso_dist) Then  !need elastic ang dist file
            !Find this interaction in the ENDF tape (MF=4, MT=2)
            MF = 4
            MT = 2
            If (.NOT.Find_MFMT(ENDF_unit,MF,MT)) Then
                Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  Section not found, MF='// & 
                                   & Trim(Adjustl(to_str(MF)))//', MT='//Trim(Adjustl(to_str(MT))) , kill=.TRUE. )
            End If
            !the next read statement on ENDF_unit will read the first line of MF=4, MT=2
            Call Read_da_sect_MF4(ENDF_unit,E_scratch,Ang_dist_scratch,n_p,LTT)
            If (LTT .EQ. 0) Then
                CS%lev_cs(i)%da(0)%n_da = 0
                CS%lev_cs(i)%da(0)%is_iso = .TRUE.
            Else
                CS%lev_cs(i)%da(0)%is_iso = .FALSE.
                If (Trim_AD_for_E(n_p,E_scratch,Ang_dist_scratch,E_min,E_max)) Then
                    Call Map_and_Store_AD( CS%n_E_uni,        & 
                                         & CS%E_uni,          & 
                                         & n_p,               & 
                                         & E_scratch,         & 
                                         & Ang_dist_scratch,  & 
                                         & CS%lev_cs(i)%da(0) )
                End If
                Deallocate(E_scratch,Ang_dist_scratch)
            End If
            If (v_full) Then  !write the stored values for elastic scatter angular distribution
                Write(v_unit,'(A,I0,A,I0,A)') Trim(isotope_names(i))//' MF=',4,', MT=',2,' (da, elastic)'
                Call Write_stored_AD(v_unit,CS%lev_cs(i)%da(0),CS%n_E_uni,CS%E_uni)
            End If
            If (LTT .EQ. 1) Then
                If (MaxVal(CS%lev_cs(i)%da(0)%da(:)%n_a) .GT. CS%n_a_max) CS%n_a_max = MaxVal(CS%lev_cs(i)%da(0)%da(:)%n_a)
            Else If (LTT .EQ. 2) Then
                If (MaxVal(CS%lev_cs(i)%da(0)%da(:)%n_a) .GT. CS%n_a_tab_max) CS%n_a_tab_max = MaxVal(CS%lev_cs(i)%da(0)%da(:)%n_a)
            Else If (LTT .EQ. 3) Then
                Do k = 1,CS%lev_cs(i)%da(0)%n_da
                    If (CS%lev_cs(i)%da(0)%da(k)%is_legendre) Then
                        If (CS%lev_cs(i)%da(0)%da(k)%n_a .GT. CS%n_a_max) CS%n_a_max = CS%lev_cs(i)%da(0)%da(k)%n_a
                    Else If (CS%lev_cs(i)%da(0)%da(k)%is_tab) Then
                        If (CS%lev_cs(i)%da(0)%da(k)%n_a .GT. CS%n_a_tab_max) CS%n_a_tab_max = CS%lev_cs(i)%da(0)%da(k)%n_a
                    End If
                End Do
            End If
        End If
        !!  INELASTIC SCATTER
        CS%lev_cs(i)%Q = -1._dp
        CS%lev_cs(i)%Q(0) = 0._dp
        CS%lev_cs(i)%thresh = -1
        CS%lev_cs(i)%thresh(0) = 0
        Do j = 1,n_inel_lev(i)
            !!  INELASTIC SCATTER INTERACTION CROSS SECTION
            !Find this interaction in the ENDF tape (MF=3, MT=50+j)
            MF = 3
            MT = 50 + j
            If (.NOT.Find_MFMT(ENDF_unit,MF,MT)) Then
                Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  Section not found, MF='// & 
                                   & Trim(Adjustl(to_str(MF)))//', MT='//Trim(Adjustl(to_str(MT))) , kill=.TRUE. )
            End If
            !the next read statement on ENDF_unit will read the first line of MF=3, MT=50+j
            Call Read_sig_sect(ENDF_unit,Q_scratch,An_scratch,E_scratch,sig_scratch,Interp_scratch,n_p,n_r)
            If (Trim_CS_for_E(n_p,E_scratch,sig_scratch,n_r,Interp_scratch,E_min,E_max)) Then
                CS%lev_cs(i)%Q(j) = Q_scratch
                Call Map_and_Store_CS( CS%n_E_uni,               & 
                                        & CS%E_uni,              & 
                                        & n_p,                   & 
                                        & E_scratch,             & 
                                        & sig_scratch,           & 
                                        & n_r,                   & 
                                        & Interp_Scratch,        & 
                                        & CS%lev_cs(i)%sig(j),   & 
                                        & CS%lev_cs(i)%thresh(j) )
            End If
            Deallocate(E_scratch,sig_scratch,Interp_scratch)
            If (v_full) Then  !write the stored values for this inelastic scatter interaction cross section
                Write(v_unit,'(A,I0,A,I0,A)') Trim(isotope_names(i))//' MF=',3,', MT=',50+j,' (sig, inelastic)'
                Call Write_stored_sig(v_unit,CS%lev_cs(i)%sig(j),CS%n_E_uni,CS%E_uni)
            End If
            !!  INELASTIC SCATTER ANGULAR DISTRIBUTION
            If (aniso_dist) Then  !need inelastic level ang dist files
                !Find this interaction in the ENDF tape (MF=4 or MF=6, MT=50+j)
                MF = 4
                MT = 50 + j
                If (.NOT.Find_MFMT(ENDF_unit,MF,MT)) Then
                    MF = 6
                    If (.NOT.Find_MFMT(ENDF_unit,MF,MT)) Then
                        Call Output_Message( 'ERROR:  Cross_Sections: Setup_Cross_Sections:  Section not found, MF='// & 
                                           & '4,6'//', MT='//Trim(Adjustl(to_str(MT))) , kill=.TRUE. )
                    End If
                End If
                LTT = -HUGE(LTT)
                LAW = -HUGE(LAW)
                Select Case (MF)
                    Case(4)
                        !the next read statement on ENDF_unit will read the first line of MF=4, MT=50+j
                        Call Read_da_sect_MF4(ENDF_unit,E_scratch,Ang_dist_scratch,n_p,LTT)
                    Case(6)
                        !the next read statement on ENDF_unit will read the first line of MF=6, MT=50+j
                        Call Read_da_sect_MF6(ENDF_unit,E_scratch,Ang_dist_scratch,n_p,LAW)
                End Select
                If (LTT.EQ.0 .OR. LAW.EQ.3) Then  !da is isotropic
                    CS%lev_cs(i)%da(j)%n_da = 0
                    CS%lev_cs(i)%da(j)%is_iso = .TRUE.
                Else  !da has been read, map and store the distribution data
                    If (Trim_AD_for_E(n_p,E_scratch,Ang_dist_scratch,E_min,E_max)) Then
                        Call Map_and_Store_AD( CS%n_E_uni, & 
                                                & CS%E_uni, & 
                                                & n_p, & 
                                                & E_scratch, & 
                                                & Ang_dist_scratch, & 
                                                & CS%lev_cs(i)%da(j), & 
                                                & CS%lev_cs(i)%thresh(j) )
                    End If
                    Deallocate(E_scratch,Ang_dist_scratch)
                End If
                If (v_full) Then  !write the stored values for this inelastic scatter angular distribution
                    Write(v_unit,'(A,I0,A,I0,A)') Trim(isotope_names(i))//' MF=',MF,', MT=',MT,' (da, inelastic)'
                    Call Write_stored_AD(v_unit,CS%lev_cs(i)%da(j),CS%n_E_uni,CS%E_uni)
                End If
                If (LTT.EQ.1 .OR. LAW.EQ.2) Then
                    Do k = 1,CS%lev_cs(i)%da(j)%n_da
                        If (CS%lev_cs(i)%da(j)%da(k)%n_a .GT. CS%n_a_max) CS%n_a_max = CS%lev_cs(i)%da(j)%da(k)%n_a
                    End Do
                Else If (LTT .EQ. 2) Then
                    Do k = 1,CS%lev_cs(i)%da(j)%n_da
                        If (CS%lev_cs(i)%da(j)%da(k)%n_a .GT. CS%n_a_tab_max) CS%n_a_tab_max = CS%lev_cs(i)%da(j)%da(k)%n_a
                    End Do
                Else If (LTT .EQ. 3) Then
                    Do k = 1,CS%lev_cs(i)%da(j)%n_da
                        If (CS%lev_cs(i)%da(j)%da(k)%is_legendre) Then
                            If (CS%lev_cs(i)%da(j)%da(k)%n_a .GT. CS%n_a_max) CS%n_a_max = CS%lev_cs(i)%da(j)%da(k)%n_a
                        Else If (CS%lev_cs(i)%da(j)%da(k)%is_tab) Then
                            If (CS%lev_cs(i)%da(j)%da(k)%n_a .GT. CS%n_a_tab_max) CS%n_a_tab_max = CS%lev_cs(i)%da(j)%da(k)%n_a
                        End If
                    End Do
                End If
            End If
        End Do
        !N2H For verbose output, write a summary of scattering modes for this isotope
    End Do
    CS%Mn = neutron_mass * Sum(CS%An) / CS%n_iso
    If (v_burn0) Then
        !write cross section traces for total, scatter, and absorption for each isotope
        Write(v_unit,*)
        Write(v_unit,'(A)') half_dash_line
        Write(v_unit,'(A)') 'CROSS SECTION TRACES FOR VERIFICATION (T = 0 K)'
        Write(v_unit,'(A)') half_dash_line
        Write(v_unit,*)
        Write(v_unit,'(A6,4A16)',ADVANCE='NO') ' Index',' E [keV]      ','Atmo sig-T [b]', & 
                                             &         'Atmo sig-S [b]','Atmo sig-A [b]'
        Do k = 1,CS%n_iso
            Write(v_unit,'(3A16)',ADVANCE='NO') isotope_names(k)//' sig-T [b]', & 
                                              & isotope_names(k)//' sig-S [b]', & 
                                              & isotope_names(k)//' sig-A [b]'
        End Do
        Write(v_unit,*)
        Write(v_unit,'(A6,4A16)',ADVANCE='NO') '------','--------------','--------------','--------------','--------------'
        Do k = 1,CS%n_iso
            Write(v_unit,'(3A16)',ADVANCE='NO') '--------------','--------------','--------------'
        End Do
        Write(v_unit,*)
        Do i = 1,CS%n_E_uni-1
            Do j = 0,inter_pts+1
                En = CS%E_uni(i) + Real(j,dp) * (CS%E_uni(i+1)-CS%E_uni(i)) / Real(inter_pts+1,dp)
                If (j .EQ. 0) Then
                    Write(v_unit,'(I6)',ADVANCE='NO') i
                Else If (j .EQ. inter_pts+1) Then
                    Write(v_unit,'(I6)',ADVANCE='NO') i + 1
                Else
                    Write(v_unit,'(A6)',ADVANCE='NO') ''
                End If
                sT = CS%sig_T(En)
                sS = CS%sig_S(En)
                sA = CS%sig_A(En)
                Write(v_unit,'(4ES16.7E2)',ADVANCE='NO') En , sT , sS , sA
                Do k = 1,CS%n_iso
                    sT = CS%sig_T(k,En)
                    sS = CS%sig_S(k,En)
                    sA = CS%sig_A(k,En)
                    Write(v_unit,'(3ES16.7E2)',ADVANCE='NO') sT , sS , sA
                End Do
                Write(v_unit,*)
                If (j.EQ.inter_pts .AND. i.LT.CS%n_E_uni-1) Exit
            End Do
        End Do
    End If
    If (v_burn300) Then
        !write cross section traces for total, scatter, and absorption for each isotope broadened to 300K
        T_broad = 300._dp
        Write(v_unit,*)
        Write(v_unit,'(A)') half_dash_line
        Write(v_unit,'(A)') 'CROSS SECTION TRACES FOR VERIFICATION (T = 300 K)'
        Write(v_unit,'(A)') half_dash_line
        Write(v_unit,*)
        Write(v_unit,'(A6,4A16)',ADVANCE='NO') 'Index',' E [keV]      ','Atmo sig-T [b]', & 
                                             &         'Atmo sig-S [b]','Atmo sig-A [b]'
        Do k = 1,CS%n_iso
            Write(v_unit,'(3A16)',ADVANCE='NO') isotope_names(k)//' sig-T [b]', & 
                                              & isotope_names(k)//' sig-S [b]', & 
                                              & isotope_names(k)//' sig-A [b]'
        End Do
        Write(v_unit,*)
        Write(v_unit,'(A6,4A16)',ADVANCE='NO') '------','--------------','--------------','--------------','--------------'
        Do k = 1,CS%n_iso
            Write(v_unit,'(3A16)',ADVANCE='NO') '--------------','--------------','--------------'
        End Do
        Write(v_unit,*)
        Do i = 1,CS%n_E_uni-1
            Do j = 0,inter_pts+1
                En = CS%E_uni(i) + Real(j,dp) * (CS%E_uni(i+1)-CS%E_uni(i)) / Real(inter_pts+1,dp)
                If (j .EQ. 0) Then
                    Write(v_unit,'(I6)',ADVANCE='NO') i
                Else If (j .EQ. inter_pts+1) Then
                    Write(v_unit,'(I6)',ADVANCE='NO') i + 1
                Else
                    Write(v_unit,'(A6)',ADVANCE='NO') ''
                End If
                sT = CS%sig_T(En,T_broad)
                sS = CS%sig_S(En,T_broad)
                sA = CS%sig_A(En,T_broad)
                Write(v_unit,'(4ES16.7E2)',ADVANCE='NO') En , sT , sS , sA
                Do k = 1,CS%n_iso
                    sT = CS%sig_T(k,En,T_broad)
                    sS = CS%sig_S(k,En,T_broad)
                    sA = CS%sig_A(k,En,T_broad)
                    Write(v_unit,'(3ES16.7E2)',ADVANCE='NO') sT , sS , sA
                End Do
                Write(v_unit,*)
                If (j.EQ.inter_pts .AND. i.LT.CS%n_E_uni-1) Exit
            End Do
        End Do
        Close(v_unit)
    End If
End Function Setup_Cross_Sections

Subroutine Write_stored_sig(v_unit,sig,n_E_uni,E_uni)
    Use Kinds, Only: dp
    Implicit None
    Integer, Intent(In) :: v_unit
    Type(sig_Type), Intent(In) :: sig
    Integer, Intent(In) :: n_E_uni
    Real(dp), Intent(In) :: E_uni(1:n_E_uni)
    Integer :: k
    Integer :: map_gap,m

    Write(v_unit,'(A8,A5)') '   Up to','  law'
    Write(v_unit,'(A8,A5)') '  ------','  ---'
    Do k = 1,sig%n_interp_r
        Write(v_unit,'(I8,I4)') sig%interp(k,1),sig%interp(k,2)
    End Do
    If (Any(sig%interp(:,2).EQ.4) .OR. Any(sig%interp(:,2).EQ.5)) Then !lnsigs are present
        Write(v_unit,'(3A16,3A9)') ' E-keyed [keV]',' sig [b]      ',' ln(sig)      ',' key   ',' index ',' map(s)'
        Write(v_unit,'(3A16,3A9)') '--------------','--------------','--------------','-------','-------','-------'
    Else  !ln(sig) not stored
        Write(v_unit,'(2A16,3A9)') ' E-keyed [keV]',' sig [b]      ',' key   ',' index ',' map(s)'
        Write(v_unit,'(2A16,3A9)') '--------------','--------------','-------','-------','-------'
    End If
    Do k = 1,sig%n_sig
        If (Any(sig%interp(:,2).EQ.4) .OR. Any(sig%interp(:,2).EQ.5)) Then !lnsigs are present
            Write(v_unit,'(3ES16.7E2,3I9)',ADVANCE='NO') E_uni(sig%E_key(k)), &
                                                       & sig%sig(k), &
                                                       & sig%lnsig(k), &
                                                       & sig%E_key(k), &
                                                       & k, &
                                                       & sig%E_map(sig%E_key(k))
        Else  !ln(sig) not stored
            Write(v_unit,'(2ES16.7E2,3I9)',ADVANCE='NO') E_uni(sig%E_key(k)), &
                                                       & sig%sig(k), &
                                                       & sig%E_key(k), &
                                                       & k, &
                                                       & sig%E_map(sig%E_key(k))
        End If
        If (k .EQ. 1) Then
            map_gap = sig%E_key(k) - 1
        Else
            map_gap = sig%E_key(k) - sig%E_key(k-1) - 1
        End If
        If (k.GT.1 .AND. map_gap.GT.0) Then !gap between cs points exists in the unified list
            !write the map values backwards through the gap
            Do m = 1,map_gap
                Write(v_unit,'(A1,I0)',ADVANCE='NO') ',' , sig%E_map(sig%E_key(k)-m)
            End Do
        End If
        Write(v_unit,*)
    End Do
    Write(v_unit,*)
End Subroutine Write_stored_sig

Subroutine Write_stored_AD(v_unit,d,n_E_uni,E_uni)
    Use Kinds, Only: dp
    Implicit None
    Integer, Intent(In) :: v_unit
    Type(da_Type), Intent(In) :: d
    Integer, Intent(In) :: n_E_uni
    Real(dp), Intent(In) :: E_uni(1:n_E_uni)
    Integer :: k
    Integer :: map_gap,m
    Integer :: a

    If (d%is_iso) Then  !isotropic distribution, no data stored
        Write(v_unit,'(A,I0,A)') 'Angular distribution is isotropic, ',d%n_da,' da points listed/stored.'
        Write(v_unit,*)
        RETURN
    End If
    Write(v_unit,'(A16,3A9)') ' E-keyed [keV]',' key   ',' index ',' map(s)'
    Write(v_unit,'(A16,3A9)') '--------------','-------','-------','-------'
    Do k = 1,d%n_da
       Write(v_unit,'(ES16.7E2,3I9)',ADVANCE='NO') E_uni(d%E_key(k)) , d%E_key(k) , k , d%E_map(d%E_key(k))
        If (k .EQ. 1) Then
            map_gap = d%E_key(k) - 1
        Else
            map_gap = d%E_key(k) - d%E_key(k-1) - 1
        End If
        If (k.GT.1 .AND. map_gap.GT.0) Then !gap between cs points exists in the unified list
            !write the map values backwards through the gap
            Do m = 1,map_gap
                Write(v_unit,'(A1,I0)',ADVANCE='NO') ',' , d%E_map(d%E_key(k)-m)
            End Do
        End If
        Write(v_unit,*)
        If (d%da(k)%is_Legendre) Then
            Write(v_unit,'(A9,I5,A16)') '',d%da(k)%n_a,' Legendre Coeffs'
            Write(v_unit,'(A9,A21)')    '',         '-------------------'
        Else !tabulated cosine pdf
            Write(v_unit,'(A9,I5,A11,2A16)') '',d%da(k)%n_a,' Cosine    ',' PDF          ',' ln(PDF)      '
            Write(v_unit,'(A9,3A16)') '',                '--------------','--------------','--------------'
        End If
        Do a = 1,d%da(k)%n_a
            If (d%da(k)%is_Legendre) Then
                Write(v_unit,'(A9,ES16.7E2)')  '',d%da(k)%a(a)
            Else !tabulated cosine pdf
                Write(v_unit,'(A9,3ES16.7E2)') '',d%da(k)%ua(1,a),Exp(d%da(k)%ua(2,a)),d%da(k)%ua(2,a)
            End If
        End Do
    End Do
    Write(v_unit,*)
End Subroutine Write_stored_AD

Subroutine Write_stored_res(v_unit,res)
    Implicit None
    Integer, Intent(In) :: v_unit
    Type(res_sig_Type), Intent(In) :: res
    Integer :: l,r

    Write(v_unit,'(2A16)') ' E-low [keV]  ',' E-high [keV] '
    Write(v_unit,'(2A16)') '--------------','--------------'
    Write(v_unit,'(2ES16.7E2)') res%E_range(1),res%E_range(2)
    If (res%is_RM) Then
        Write(v_unit,'(A)') '   Formalism:  Riech-Moore'
    Else If (res%is_SLBW) Then
        Write(v_unit,'(A)') '   Formalism:  Single Level Breit-Wigner'
    Else If (res%is_MLBW) Then
        Write(v_unit,'(A)') '   Formalism:  Multi Level Breit-Wigner'
    End If
    Do l = 1,res%n_L
        Write(v_unit,'(A,I0,2(A,ES16.7E2))') '  Level ',l,' -- Channel Rad [1E-12cm]:',res%L(l)%a, & 
                                                       & ', Scattering Rad [1E-12cm]:',res%L(l)%ap
        Write(v_unit,'(A4,2A16)',ADVANCE='NO') '',' E-res [eV]   ',' aJ           '
        If (res%is_RM) Then
            Write(v_unit,'(2A16)')       ' Neutron Width',' Rad Width    '
            Write(v_unit,'(A4,4A16)') '','--------------','--------------','--------------','--------------'
        Else If (res%is_BW) Then
            Write(v_unit,'(3A16)')       ' Total Width  ',' Neutron Width',' Rad Width    '
            Write(v_unit,'(A4,5A16)') '','--------------','--------------','--------------','--------------','--------------'
        End If
        Do r = 1,res%L(l)%n_r
            If (res%is_RM) Then
                Write(v_unit,'(A4,4ES16.7E2)') '',res%L(l)%EraG(r,1:4)
            Else If (res%is_MLBW) Then
                Write(v_unit,'(A4,5ES16.7E2)') '',res%L(l)%EraG(r,1:5)
            End If
        End Do
    End Do
    Write(v_unit,*)
End Subroutine Write_stored_res

Subroutine Find_MFMT_end(ENDF_unit)
    Implicit None
    Integer, Intent(In) :: ENDF_unit
    Character(77) :: trash_c
    Integer :: line_num

    Do
        Read(ENDF_unit,'(A75,I5)') trash_c,line_num
        If (line_num .EQ. 99999) Return
    End Do
End Subroutine Find_MFMT_end

Function Find_MFMT(ENDF_unit,MFi,MTi) Result(bingo)
    Implicit None
    Logical :: bingo
    Integer, Intent(In) :: ENDF_unit
    Integer, Intent(In) :: MFi,MTi
    Integer :: MF,MT
    Character(80) :: trash_c
    Integer :: stat

    bingo = .FALSE.
    !Find this interaction (MFi,MTi) in the ENDF tape
    Rewind(ENDF_unit)  !return to the start of the file
    Do
        Read(ENDF_unit,'(A71,I1,I3)',IOSTAT=stat) trash_c,MF,MT
        If (stat .LT. 0) Return
        If (MF.EQ.MFi .AND. MT.EQ.MTi) Then
            bingo = .TRUE.
            Backspace(ENDF_unit)  !go back one line
            Exit
        End If
    End Do
    !the next read statement on ENDF_unit will read the first line of MFi, MTi
End Function Find_MFMT

Subroutine Read_sig_sect(cs_unit,Q,An,E_list,sig_list,Int_list,n_p,n_r)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Integer, Intent(In) :: cs_unit  !unit number for the ENDF tape file
    !cursor must already be positioned in the file so that the next read will be the first line of the appropriate section
    Real(dp), Intent(Out) :: Q
    Real(dp), Intent(Out) :: An
    Real(dp), Allocatable, Intent(Out) :: E_list(:)
    Real(dp), Allocatable, Intent(Out) :: sig_list(:)
    Integer, Allocatable, Intent(Out) :: Int_list(:,:)
    Integer, Intent(Out) :: n_p,n_r
    Integer :: i
    Real(dp) :: trash

    !the second entry of the first line is the mass of the target in neutron masses
    Read(cs_unit,'(2E11.6E1)') trash, An
    !the second entry of the next line is the Q-value for the reaction, the fifth and sixth entries are the number of interpolation 
    !ranges and energy levels in the file
    Read(cs_unit,'(4E11.6E1,2I11)') trash, Q, trash, trash, n_r, n_p
    If (n_r .GT. 3) Call Output_Message('ERROR:  Cross_Sections: Read_sig_sect:  Number of interpolation ranges &
                                        &greater than 3:  n_r=',n_r,kill=.TRUE.)
    !The next line has interpolation methods and ranges, read and fill interpolation list
    Allocate(Int_list(1:n_r,1:2))
    Int_list = 0
    If (n_r .EQ. 1) Then
        Read(cs_unit,'(2I11)') Int_list(1,1), Int_list(1,2)
    Else If (n_r .EQ. 2) Then
        Read(cs_unit,'(4I11)') Int_list(1,1), Int_list(1,2), Int_list(2,1), Int_list(2,2)
    Else !n_r=3
        Read(cs_unit,'(6I11)') Int_list(1,1), Int_list(1,2), Int_list(2,1), Int_list(2,2), Int_list(3,1), Int_list(3,2)
    End If
    If (Any(Int_list(:,2).GT.5) .OR. Any(Int_list(:,2).LT.1)) Then
        Call Output_Message('ERROR:  Cross_Sections: Read_CS_file:  Unknown interpolation scheme.',kill=.TRUE.)
    End If
    !Allocate lists
    Allocate(E_list(1:n_p))
    Allocate(sig_list(1:n_p))
    !Read the rest of the file
    Do i = 1,n_p,3
        !Each line in the file has 3 pairs of (eV,barns)
        If (n_p-i .GT. 1) Then
            Read(cs_unit,'(6E11.6E1)') E_list(i), sig_list(i), E_list(i+1), sig_list(i+1), E_list(i+2), sig_list(i+2)
        Else If (n_p-i .EQ. 1) Then
            Read(cs_unit,'(4E11.6E1)') E_list(i), sig_list(i), E_list(i+1), sig_list(i+1)
        Else
            Read(cs_unit,'(2E11.6E1)') E_list(i), sig_list(i)
        End If
    End Do
    !Convert E from eV to keV
    Q = Abs(Q) / 1000._dp  !Also store Q as positive value, negativity of Q for inelastic scatter is handled by energy 
                           !book-keeeping formulations during transport
    E_list = E_list / 1000._dp
End Subroutine Read_sig_sect

Subroutine Read_da_sect_MF4(da_unit,E_list,da_list,n_p,LTT)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Integer, Intent(In) :: da_unit  !unit number for the ENDF tape file
    !cursor must already be positioned in the file so that the next read will be the first line of the appropriate section
    Real(dp), Allocatable, Intent(Out) :: E_list(:)
    Type(da_List_type), Allocatable, Intent(Out) :: da_list(:)
    Integer, Intent(Out):: n_p,LTT
    Integer :: i,j
    Real(dp) :: trash
    Character(80) :: trash_c
    Integer :: n_a,n_a_lines,n_skipped_lines
    Integer :: LCT
    Integer :: n_p_2
    Logical :: new_line

    !Check the LTT value on the first line to ensure Legendre Coeffs format
    Read(da_unit,'(A33,I11)') trash_c, LTT
    Read(da_unit,'(A33,I11)') trash_c, LCT
    !values must be in CM frame, check and throw error if they are not
    If (LCT .NE. 2) Call Output_Message('ERROR:  Cross_Sections: Read_da_sect:  Incorrectly formatted file, LCT=',LCT,kill=.TRUE.)
    If (LTT .EQ. 0) Return  !LTT=0 indicates isotropic distribution, special handling in calling routine
    !The 6th entry on the next line energy levels in the file (for LTT=1 or 2, if LTT=3 it is the number of energies for the 
    !legendre section)
    Read (da_unit,'(A55,I11)') trash_c,n_p
    !Skip the next line
    Read (da_unit,*)
    If (LTT .EQ. 3) Then !there is a second range of energies later in the section
        !Advance in the file to the end of the Legendre section
        n_skipped_lines = 0
        Do i = 1,n_p
            !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficients 
            !in the 5th position
            Read(da_unit,'(A44,I11)') trash_c, n_a
            n_skipped_lines = n_skipped_lines + 1
            n_a_lines = (n_a / 6)  !integer divide
            If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
            !advance to the next energy
            Do j = 1,n_a_lines
                Read(da_unit,*)
                n_skipped_lines = n_skipped_lines + 1
            End Do
        End Do
        !the 6th entry on the next line energy levels in the tabular section of the file
        Read (da_unit,'(A55,I11)') trash_c, n_p_2
        n_skipped_lines = n_skipped_lines + 1
        !need to rewind the file to the start of the legendre section
        Do i = 1,n_skipped_lines
            Backspace(da_unit)
        End Do
    Else
        n_p_2 = 0
    End If
    !Allocate the data array to the number of provided energy levels
    Allocate(E_list(1:n_p+n_p_2))
    Allocate(da_list(1:n_p+n_p_2))
    If (LTT .EQ. 1) Then  !da is lists of legendre coeffs
        da_list%is_legendre = .TRUE.
        da_list%is_tab = .FALSE.
        Do i = 1,n_p
            !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficients 
            !in the 5th position
            Read(da_unit,'(4E11.6E1,I11)') trash, E_list(i), trash, trash, da_list(i)%n_a
            Allocate(da_list(i)%a(0:da_list(i)%n_a))
            da_list(i)%a(0) = 0.5_dp
            Do j = 1,da_list(i)%n_a
                If (mod(j,6).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                    Read(da_unit,'(E11.6E1)', ADVANCE = 'YES') da_list(i)%a(j)
                Else !read the entry without advancing
                    Read(da_unit,'(E11.6E1)', ADVANCE = 'NO') da_list(i)%a(j)
                End If
                !multiply the coeff by (2k+1)/2
                da_list(i)%a(j) = 0.5_dp * da_list(i)%a(j) * Real(2*j + 1,dp)
            End Do
        End Do
    Else If (LTT .EQ. 2) Then  !da is tabulated
        da_list%is_legendre = .FALSE.
        da_list%is_tab = .TRUE.
        Do i = 1,n_p
            !The first line in each energy contains the energy in eV in the second position
            Read(da_unit,'(2E11.6E1)') trash, E_list(i)
            !the next line contains the number of tabulation points int he first position
            Read(da_unit,'(I11)') da_list(i)%n_a
            Allocate(da_list(i)%a(0:da_list(i)%n_a))
            da_list(i)%a = 0._dp
            Allocate(da_list(i)%ua(1:2,1:da_list(i)%n_a))
            da_list(i)%ua = 0._dp
            Do j = 1,da_list(i)%n_a
                If (mod(j,3).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                    Read(da_unit,'(2E11.6E1)', ADVANCE = 'YES') da_list(i)%ua(1,j), da_list(i)%ua(2,j)
                Else !read the entry without advancing
                    Read(da_unit,'(2E11.6E1)', ADVANCE = 'NO') da_list(i)%ua(1,j), da_list(i)%ua(2,j)
                End If
            End Do
        End Do
        da_list(i)%ua(2,:) = Log(da_list(i)%ua(2,:)) !Store logarithm of probability density (reduces cost of interpolation)
    Else If (LTT .EQ. 3) Then  !da is tabulated for high energies but legendre for low energies
        new_line = .FALSE.
        !Read in low energy Legendre points
        Do i = 1,n_p
            da_list(i)%is_legendre = .TRUE.
            da_list(i)%is_tab = .FALSE.
            !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficients 
            !in the 5th position
            Read(da_unit,'(4E11.6E1,I11)') trash, E_list(i), trash, trash, da_list(i)%n_a
            Allocate(da_list(i)%a(0:da_list(i)%n_a))
            da_list(i)%a(0) = 0.5_dp
            Do j = 1,da_list(i)%n_a
                If (mod(j,6).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                    Read(da_unit,'(E11.6E1)', ADVANCE = 'YES') da_list(i)%a(j)
                    new_line = .TRUE.
                Else !read the entry without advancing
                    Read(da_unit,'(E11.6E1)', ADVANCE = 'NO') da_list(i)%a(j)
                    new_line = .FALSE.
                End If
                !multiply the coeff by (2k+1)/2
                da_list(i)%a(j) = 0.5_dp * da_list(i)%a(j) * Real(2*j + 1,dp)
            End Do
        End Do
        !advance one or two lines based on where Legendre Coeff reading finished
        If (.NOT. new_line) Read(da_unit,*)
        Read(da_unit,*)
        !first entry of next line is number of additional energy points
        Read (da_unit,'(I11)') n_p_2
        !Read in high energy tabulated cosine points
        Do i = n_p+1,n_p+n_p_2
            da_list(i)%is_legendre = .FALSE.
            da_list(i)%is_tab = .TRUE.
            !The first line in each energy contains the energy in eV in the second position
            Read(da_unit,'(2E11.6E1)') trash, E_list(i)
            !the next line contains the number of tabulation points in the first position
            Read(da_unit,'(I11)') da_list(i)%n_a
            Allocate(da_list(i)%a(0:da_list(i)%n_a))
            da_list(i)%a = 0._dp
            Allocate(da_list(i)%ua(1:2,1:da_list(i)%n_a))
            da_list(i)%ua = 0._dp
            Do j = 1,da_list(i)%n_a
                If (mod(j,3).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                    Read(da_unit,'(2E11.6E1)', ADVANCE = 'YES') da_list(i)%ua(1,j), da_list(i)%ua(2,j)
                Else !read the entry without advancing
                    Read(da_unit,'(2E11.6E1)', ADVANCE = 'NO') da_list(i)%ua(1,j), da_list(i)%ua(2,j)
                End If
            End Do
            da_list(i)%ua(2,:) = Log(da_list(i)%ua(2,:)) !Store logarithm of probability density (reduces cost of interpolation)
        End Do
        !update total n_p
        n_p = n_p + n_p_2
    End If
    E_list = E_list / 1000._dp
End Subroutine Read_da_sect_MF4

Subroutine Read_da_sect_MF6(da_unit,E_list,da_list,n_p,LAW)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Integer, Intent(In) :: da_unit  !unit number for the ENDF tape file
    !cursor must already be positioned in the file so that the next read will be the first line of the appropriate section
    Real(dp), Allocatable, Intent(Out) :: E_list(:)
    Type(da_List_type), Allocatable, Intent(Out) :: da_list(:)
    Integer, Intent(Out):: n_p,LAW
    Integer :: i,j
    Real(dp) :: trash
    Integer :: trash_i
    Character(80) :: trash_c
    Integer :: LCT,LIP,LANG
    Real(dp) :: ZAP,AWP

    !Check the LTT value on the first line to ensure Legendre Coeffs format
    Read(da_unit,'(A33,I11)') trash_c, LCT
    !values must be in CM frame, check and throw error if they are not
    If (LCT .NE. 2) Call Output_Message('ERROR:  Cross_Sections: Read_da_sect:  Incorrectly formatted file, LCT=',LCT,kill=.TRUE.)
    !find the start of the section describing the angular distribution of the neutron on the exit channel
    Do
        Read(da_unit,'(2E11.6E1,2I11)') ZAP,AWP,LIP,LAW
        If (ZAP.EQ.1._dp .AND. AWP.EQ.1._dp .AND. LAW.EQ.2) Exit !this is first line in the correct section
    End Do
    If (LAW .EQ. 3) Return  !LAW=3 indicates isotropic distribution, special handling in calling routine
    !the next 2 lines are always discarded
    Read(da_unit,*)
    Read(da_unit,*)
    !The 6th entry on the next line energy levels in the section
    Read(da_unit,'(A55,I11)') trash_c,n_p
    !Skip the next line
    Read (da_unit,*)
    !Allocate the data array to the number of provided energy levels
    Allocate(E_list(1:n_p))
    Allocate(da_list(1:n_p))
    !Assume Legendre coeffs description of angular distribution, otherwise throw error and stop
    da_list%is_legendre = .TRUE.
    da_list%is_tab = .FALSE.
    Do i = 1,n_p
        !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficients 
        !in the 5th position
        Read(da_unit,'(2E11.6E1,4I11)') trash, E_list(i), LANG, trash_i, da_list(i)%n_a
        If (LANG .NE. 0) Then
            Call Output_Message('ERROR:  Cross_Sections: Read_da_sect:  Incorrectly formatted da, LANG=',LANG,kill=.TRUE.)
        End If
        Allocate(da_list(i)%a(0:da_list(i)%n_a))
        da_list(i)%a(0) = 0.5_dp
        Do j = 1,da_list(i)%n_a
            If (mod(j,6).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                Read(da_unit,'(E11.6E1)', ADVANCE = 'YES') da_list(i)%a(j)
            Else !read the entry without advancing
                Read(da_unit,'(E11.6E1)', ADVANCE = 'NO') da_list(i)%a(j)
            End If
            !multiply the coeff by (2k+1)/2
            da_list(i)%a(j) = 0.5_dp * da_list(i)%a(j) * Real(2*j + 1,dp)
        End Do
    End Do
    E_list = E_list / 1000._dp
End Subroutine Read_da_sect_MF6

Subroutine Read_res_sect(res_unit,res_List)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Use Global, Only: mn => neutron_mass ![kg]
    Use Global, Only: h_bar => h_bar_Planck ![J*s]
    Use Global, Only: keV_per_Joule ![keV/J]
    Use Sorting, Only: Union_Sort
    Implicit None
    Integer, Intent(In) :: res_unit
    Type(res_sig_Type), Intent(Out) :: res_list
    Integer :: l,r,d
    Real(dp) :: trash
    Integer :: LRF,NRO,NAPS
    Real(dp) :: SPI,AP0,AP
    Integer :: nL,nR
    Real(dp) :: AWRI
    Real(dp), Parameter :: one_third = 1._dp / 3._dp
    Logical :: found_fission

    !check the LRF value on the third line to ensure MLBW or Reich-Moore format
    !check the NRO value on the third line to ensure energy independent scattering radius
    Read(res_unit,*)
    Read(res_unit,*)
    Read(res_unit,'(3E11.6E1,3I11)') res_list%E_range(1), res_list%E_range(2), trash, LRF, NRO, NAPS
    res_list%E_range = res_list%E_range / 1000._dp  !convert to keV
    If (LRF.EQ.1 .OR. LRF.EQ.2) Then  !BW formalism
        res_list%is_RM = .FALSE.
        res_list%is_BW = .TRUE.
        If (LRF .EQ. 1) Then  !SLBW
            res_list%is_SLBW = .TRUE.
            res_list%is_MLBW = .FALSE.
        Else  !MLBW
            res_list%is_SLBW = .FALSE.
            res_list%is_MLBW = .TRUE.
        End If
    Else If (LRF .EQ. 3) Then  !Reich-Moore formalism
        res_list%is_BW = .FALSE.
        res_list%is_SLBW = .FALSE.
        res_list%is_MLBW = .FALSE.
        res_list%is_RM = .TRUE.
    Else  !this formalism is not supported
        Call Output_Message('ERROR:  Cross_Sections: Read_res_sect:  Incorrectly formatted file, LRF=',LRF,kill=.TRUE.)
    End If
    If (NRO .NE. 0) Then
        Call Output_Message('ERROR:  Cross_Sections: Read_res_sect:  Incorrectly formatted file, NRO=',NRO,kill=.TRUE.)
    End If
    !read AP and n_L from the next line
    Read(res_unit,'(4E11.6E1,I11)') SPI, AP0, trash, trash, nL
    res_list%spi = SPI
    If (res_list%is_RM) Then
        res_list%s = (/ NINT(Abs(SPI-0.5_dp)) , NINT(SPI+0.5_dp) /)
    Else
        res_list%s = 0
    End If
    !read AWRI from next line
    Read(res_unit,'(E11.6E1)') AWRI
    !compute k0 (when multiplied by Sqrt(E[eV]), gives neutron wave number)
    res_List%k0 = (Sqrt(2._dp*mn) / h_bar) * (AWRI / (AWRI + 1._dp)) * Sqrt(1._dp / (1000._dp*keV_per_Joule)) * 1.E-14_dp
    Backspace(res_unit)  !that line will need to be read again
    !allocate levels
    res_List%n_L = nL
    Allocate(res_List%L(1:nL))
    Do l = 1,nL
        res_List%L(l)%ap = AP0  !default value
        If (NAPS .EQ. 0) Then
            res_list%L(l)%a = (0.08_dp + 0.123_dp * AWRI**one_third)
        Else
            res_list%L(l)%a = AP0
        End If
        !read AWRI, APL and number of resonances in first layer from next line
        Read(res_unit,'(5E11.6E1,I11)') AWRI, AP, trash, trash, trash, nR
        res_list%L(l)%n_r = nR
        If (AP .NE. 0._dp) res_List%L(l)%ap = AP  !level specific scattering radius
        If (res_list%is_BW) Then
            d = 7
        Else
            d = 6
        End If
        Allocate(res_list%L(l)%EraG(1:nR,1:d))
        res_list%L(l)%EraG = 0._dp
        !read resonance parameters for this level
        Do r = 1,res_list%L(l)%n_r
            Read(res_unit,'(6E11.6E1)') res_list%L(l)%EraG(r,1), res_list%L(l)%EraG(r,2), res_list%L(l)%EraG(r,3), & 
                                      & res_list%L(l)%EraG(r,4), res_list%L(l)%EraG(r,5), res_list%L(l)%EraG(r,6)
        End Do
        !check for fission resonances
        found_fission = .FALSE.
        If (res_list%is_BW) Then
            If (Any(res_list%L(l)%EraG(:,6).NE.0._dp)) found_fission = .TRUE.
        Else If (res_list%is_RM) Then
            If (Any(res_list%L(l)%EraG(:,5:6).NE.0._dp)) found_fission = .TRUE.
        End If
        If (found_fission) Call Output_Message( 'ERROR:  Cross_Sections: Read_res_sect:  Dude, a fissionable '// &
                                              & 'atmosphere is just ridiculous. ',kill=.TRUE.)
        !precompute shift and penetrability factors at each stored resonance
        Do r = 1,res_list%L(l)%n_r
            Call PSphi( l-1 , res_List%L(l)%a * res_list%k0*Sqrt( Abs(res_list%L(l)%EraG(r,1)) ) , & 
                        & res_list%L(l)%EraG(r,d) , res_list%L(l)%EraG(r,d-1) )
        End Do
    End Do
End Subroutine Read_res_sect

Recursive Pure Subroutine PSphi(l,rho,P,S,phi)
    Use Kinds, Only: dp
    Implicit None
    Integer, Intent(In) :: l
    Real(dp), Intent(In) :: rho
    Real(dp), Intent(Out) :: P
    Real(dp), Intent(Out) :: S
    Real(dp), Intent(Out), Optional :: phi
    Real(dp) :: s0,p0,phi0
    Real(dp) :: c

    If (l .GT. 0) Then
        If (Present(phi)) Then
            Call PSphi(l-1,rho,p0,s0,phi0)
        Else
            Call PSphi(l-1,rho,p0,s0)
        End If
        s0 = Real(l,dp) - s0
        c = rho**2 / ( s0**2 + p0**2 )
        P = c * p0
        S = c * s0 - Real(l,dp)
        If (Present(phi)) phi = phi0 - ATAN( p0 / s0 )
    Else
        p = rho
        s = 0._dp
        If (Present(phi)) phi = rho
    End If
End Subroutine PSphi

Pure Function Phi_hard(l,rho) Result(phi)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: phi
    Integer, Intent(In) :: l
    Real(dp), Intent(In) :: rho
    Real(dp) :: p0,s0,phi0

    If (l .GT. 0) Then
        Call PSphi(l-1,rho,p0,s0,phi0)
        phi = phi0 - ATAN( p0 / (Real(l,dp) - s0) )
    Else
        phi = rho
    End If
End Function Phi_hard

Function Trim_CS_for_E(n_p,E_list,CS_list,n_r,Int_list,E_min,E_max) Result(bingo)
    Use Kinds, Only: dp
    Use Interpolation, Only: Linear_Interp
    Implicit None
    Logical :: bingo
    Integer, Intent(InOut) :: n_p
    Real(dp), Allocatable, Intent(InOut) :: E_list(:)
    Real(dp), Allocatable, Intent(InOut) :: CS_list(:)
    Integer, Intent(InOut) :: n_r
    Integer, Allocatable, Intent(InOut) :: Int_list(:,:)
    Real(dp), Intent(In) :: E_min,E_max
    Real(dp), Allocatable :: E_swap(:)
    Real(dp), Allocatable :: CS_swap(:)
    Integer, Allocatable :: Int_swap(:,:)
    Integer :: i,j

    bingo = .TRUE.
    If (E_list(1).GE.E_min .AND. E_list(n_p).LE.E_max) Return  !no trimming required
    If (E_list(1).GT.E_max .OR. E_list(n_p).LT.E_min) Then !there are no values in range
        bingo = .FALSE.
        n_p = 0
        Return
    End If
    !otherwise, need to find the range of E_list to keep
    i = 1
    Do
        If (E_list(i) .LT. E_min) Then
            i = i + 1
        Else
            Exit
        End If
    End Do
    j = n_p
    Do
        If (E_list(j) .GT. E_max) Then
            j = j - 1
        Else
            Exit
        End If
    End Do
    !compute new number of elements
    n_p = j - i + 1
    !resize E_list
    Allocate(E_swap(1:n_p))
    E_swap = E_list(i:j)
    Deallocate(E_list)
    Allocate(E_list(1:n_p))
    E_list = E_swap
    Deallocate(E_swap)
    !resize CS_list
    Allocate(CS_swap(1:n_p))
    CS_swap = CS_list(i:j)
    Deallocate(CS_list)
    Allocate(CS_list(1:n_p))
    CS_list = CS_swap
    Deallocate(CS_swap)
    !adjust Int_list for new indexes
    Int_list(:,1) = Int_list(:,1) - (i-1)
    If (Any(Int_list(:,1) .LT. 2)) Then !interpolation ranges at the beginning of the list need to be removed
        Do i = 2,n_r
            If (Int_list(i,1) .GE. 2) Then !remove the ranges before this one
                Allocate(Int_swap(1:n_r-(i-1),1:2))
                Int_swap = Int_list(i:n_r,:)
                Deallocate(Int_list)
                n_r = n_r - (i-1)
                Allocate(Int_list(1:n_r,1:2))
                Int_list = Int_swap
                Deallocate(Int_swap)
                Exit
            End If
        End Do
    End If
    If (Any(Int_list(:,1) .GT. n_p)) Then !interpolation ranges at the end of the list need to be removed
        Do i = 1,n_r-1
            If (Int_list(i,1) .GE. n_p) Then !remove the ranges after this one
                Allocate(Int_swap(1:i,1:2))
                Int_swap = Int_list(1:i,:)
                Deallocate(Int_list)
                n_r = i
                Allocate(Int_list(1:n_r,1:2))
                Int_list = Int_swap
                Int_list(n_r,1) = n_p
                Deallocate(Int_swap)
                Exit
            End If
        End Do
    End If
End Function Trim_CS_for_E

Elemental Subroutine copy_da(da_get,da_put)
    Implicit None
    Type(da_list_Type), Intent(In) :: da_get
    Type(da_list_Type), Intent(Out) :: da_put

    da_put%n_a = da_get%n_a
    da_put%is_Legendre = da_get%is_legendre
    da_put%is_Tab = da_get%is_Tab
    If (da_put%is_Legendre) Then
        Allocate(da_put%a(0:da_put%n_a))
        da_put%a = da_get%a
    Else If (da_put%is_Tab) Then
        Allocate(da_put%ua(1:2,1:da_put%n_a))
        da_put%ua = da_get%ua
    End If
End Subroutine copy_da

Function Trim_AD_for_E(n_p,E_list,AD_list,E_min,E_max) Result(bingo)
    Use Kinds, Only: dp
    Implicit None
    Logical :: bingo
    Integer, Intent(InOut) :: n_p
    Real(dp), Allocatable, Intent(InOut) :: E_list(:)
    Type(da_List_Type), Allocatable, Intent(InOut) :: AD_list(:)
    Real(dp), Intent(In) :: E_min,E_max
    Real(dp), Allocatable :: E_swap(:)
    Type(da_List_Type), Allocatable :: AD_swap(:)
    Integer :: i,j

    bingo = .TRUE.
    If (E_list(1).GE.E_min .AND. E_list(n_p).LE.E_max) Return  !no trimming required
    If (E_list(1).GT.E_max .OR. E_list(n_p).LT.E_min) Then !there are no values in range
        bingo = .FALSE.
        n_p = 0
        Return
    End If
    !otherwise, need to find the range of E_list to keep
    i = 1
    Do
        If (E_list(i) .LT. E_min) Then
            i = i + 1
        Else
            Exit
        End If
    End Do
    j = n_p
    Do
        If (E_list(j) .GT. E_max) Then
            j = j - 1
        Else
            Exit
        End If
    End Do
    !compute new number of elements
    n_p = j - i + 1
    !resize E_list
    Allocate(E_swap(1:n_p))
    E_swap = E_list(i:j)
    Deallocate(E_list)
    Allocate(E_list(1:n_p))
    E_list = E_swap
    Deallocate(E_swap)
    !resize CS_list
    Allocate(AD_swap(1:n_p))
    Call copy_da(AD_list(i:j),AD_swap(:))
    Deallocate(AD_list)
    Allocate(AD_list(1:n_p))
    Call copy_da(AD_swap(:),AD_list(:))
End Function Trim_AD_for_E

Subroutine Map_and_Store_CS(n_E_uni,E_uni,n_p,E_list,CS_list,n_r,Int_list,cs,i_thresh)
    Use Kinds, Only: dp
    Implicit None
    Integer, Intent(In) :: n_E_uni
    Real(dp), Intent(In) :: E_uni(1:n_E_uni)
    Integer, Intent(In) :: n_p
    Real(dp), Intent(In) :: E_list(1:n_p)
    Real(dp), Intent(In) :: CS_list(1:n_p)
    Integer, Intent(In) :: n_r
    Integer, Intent(In) :: Int_list(1:n_r,1:2)
    Type(sig_Type), Intent(Out) :: cs
    Integer, Optional, Intent(Out) :: i_thresh
    Integer :: i,j,t

    !Determine index of threshold energy
    t = 1
    If (Present(i_thresh)) Then  !determine index of the threshold energy, index of map will start at the threshold energy
        Do i = 1,n_E_uni
            If (E_list(1) .LE. E_uni(i)) Then
                i_thresh = i
                Exit
            End If
        End Do
        t = i_thresh
    End If
    !Create map
    Allocate(cs%E_map(t:n_E_uni))
    cs%E_map = -1
    j = 1
    Do i = t,n_E_uni
        If (E_uni(i) .GT. E_list(j)) j = j + 1
        If (j .GE. n_p) Then
            cs%E_map(i:n_E_uni) = j
            Exit
        End If
        cs%E_map(i) = j
        Do !increment the index past any equal or duplicate points
            If (E_uni(i) .EQ. E_list(j)) j = j + 1
            If (E_uni(i) .NE. E_list(j)) Exit
        End Do
    End Do
    !Create key
    Allocate(cs%E_key(1:n_p))
    cs%E_key = -1
    Do i = 1,n_p
        Do j = 1,n_E_uni
            If (E_uni(j) .EQ. E_list(i)) Then
                cs%E_key(i) = j
            End If
        End Do
    End Do
    !Store cross sections
    cs%n_sig = n_p
    Allocate(cs%sig(1:n_p))
    cs%sig = CS_list
    !store interpolation ranges
    cs%n_interp_r = n_r
    Allocate(cs%interp(1:n_r,1:2))
    cs%interp = Int_list
    !Store natural logarithms of cross sections if needed based on interpolation schemes
    If (Any(Int_list(:,2).EQ.4) .OR. Any(Int_list(:,2).EQ.5)) Then
        Allocate(cs%lnsig(1:n_p))
        cs%lnsig = -Huge(cs%lnsig)
        Where (CS_list .GT. 0._dp) cs%lnsig = Log(CS_list)
    End If
End Subroutine Map_and_Store_CS

Subroutine Map_and_Store_AD(n_E_uni,E_uni,n_p,E_list,AD_list,ad,i_thresh)
    Use Kinds, Only: dp
    Implicit None
    Integer, Intent(In) :: n_E_uni
    Real(dp), Intent(In) :: E_uni(1:n_E_uni)
    Integer, Intent(In) :: n_p
    Real(dp), Intent(In) :: E_list(1:n_p)
    Type(da_List_Type), Intent(In) :: AD_list(1:n_p)
    Type(da_Type), Intent(Out) :: ad
    Integer, Optional, Intent(In) :: i_thresh
    Integer :: i,j,t

    !Determine index of threshold energy
    t = 1
    If (Present(i_thresh)) t = i_thresh
    !Create map
    Allocate(ad%E_map(t:n_E_uni))
    ad%E_map = -1
    j = 1
    Do i = t,n_E_uni
        If (E_uni(i) .GT. E_list(j)) j = j + 1
        If (j .GE. n_p) Then
            ad%E_map(i:n_E_uni) = j
            Exit
        End If
        ad%E_map(i) = j
        Do !increment the index past any equal or duplicate points
            If (E_uni(i) .EQ. E_list(j)) j = j + 1
            If (j .GE. n_p) Exit
            If (E_uni(i) .NE. E_list(j)) Exit
        End Do
    End Do
    !Create key
    Allocate(ad%E_key(1:n_p))
    ad%E_key = -1
    Do i = 1,n_p
        Do j = 1,n_E_uni
            If (E_uni(j) .EQ. E_list(i)) Then
                ad%E_key(i) = j
            End If
        End Do
    End Do
    !Store cross sections
    ad%n_da = n_p
    Allocate(ad%da(1:n_p))
    Call copy_da(AD_list(:),ad%da(:))
End Subroutine Map_and_Store_AD

Subroutine sig_T_A_all(CS,E,sT,sA,iE_get,iE_put)
!returns total and absorption cross section (sT,sA) for the total atmosphere
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E
    Real(dp), Intent(Out) :: sT,sA
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Real(dp) :: sS
    Real(dp) :: resT,resS
    Integer :: i

    If (Present(iE_put)) Then
        E_index = Max(iE_put,2)
    Else
        E_index = Min(Bisection_Search(E,CS%E_uni,CS%n_E_uni),CS%n_E_uni)
    End If
    sS = 0._dp
    sA = 0._dp
    Do i = 1,CS%n_iso
        !resonant contribution
        If (CS%has_res_cs(i)) Then
            Call sig_Resonance(CS%res_cs(i),E,resT,resS)
            sS = sS + CS%iso_Fractions(i) * resS
            sA = sA + CS%iso_fractions(i) * (resT - resS)
        End If
        !background contribution
        sS = sS + CS%iso_Fractions(i) * sig_Composite( E, & 
                                                     & CS%n_E_uni,          & 
                                                     & CS%E_uni,            & 
                                                     & CS%lnE_uni,          & 
                                                     & E_index,             & 
                                                     & 0,                   & 
                                                     & CS%lev_cs(i)%n_lev,  & 
                                                     & CS%lev_cs(i)%thresh, & 
                                                     & CS%lev_cs(i)%sig     )
        sA = sA + CS%iso_Fractions(i) * sig_Composite( E,                    & 
                                                     & CS%n_E_uni,           & 
                                                     & CS%E_uni,             & 
                                                     & CS%lnE_uni,           & 
                                                     & E_index,              & 
                                                     & 1,                    & 
                                                     & CS%abs_cs(i)%n_modes, & 
                                                     & CS%abs_cs(i)%thresh,  & 
                                                     & CS%abs_cs(i)%sig      )
    End Do
    sT = sA + sS
    If (Present(iE_get)) iE_get = E_index
End Subroutine sig_T_A_all

Subroutine sig_T_A_iso(CS,iso,E,sT,sA,iE_get,iE_put)
!returns total and absorption cross section (sT,sA) for a single isotope (iso)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E
    Real(dp), Intent(Out) :: sT,sA
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Real(dp) :: sS
    Real(dp) :: resT,resS

    If (Present(iE_put)) Then
        E_index = Max(iE_put,2)
    Else
        E_index = Min(Bisection_Search(E,CS%E_uni,CS%n_E_uni),CS%n_E_uni)
    End If
    sS = 0._dp
    sA = 0._dp
    !resonant contribution
    If (CS%has_res_cs(iso)) Then
        Call sig_Resonance(CS%res_cs(iso),E,resT,resS)
        sS = sS + resS
        sA = sA + resT - resS
    End If
    !background contribution
    sS = sS + sig_Composite( E, & 
                           & CS%n_E_uni,            & 
                           & CS%E_uni,              & 
                           & CS%lnE_uni,            & 
                           & E_index,               & 
                           & 0,                     & 
                           & CS%lev_cs(iso)%n_lev,  & 
                           & CS%lev_cs(iso)%thresh, & 
                           & CS%lev_cs(iso)%sig     )
    sA = sA + sig_Composite( E,                      & 
                           & CS%n_E_uni,             & 
                           & CS%E_uni,               & 
                           & CS%lnE_uni,             & 
                           & E_index,                & 
                           & 1,                      & 
                           & CS%abs_cs(iso)%n_modes, & 
                           & CS%abs_cs(iso)%thresh,  & 
                           & CS%abs_cs(iso)%sig      )
    sT = sA + sS
    If (Present(iE_get)) iE_get = E_index
End Subroutine sig_T_A_iso

Function sig_T_all(CS,E,iE_get,iE_put) Result(sT)
!returns total cross section (sig_T_all) for the total atmosphere
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sT
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Real(dp) :: sA,sS
    Real(dp) :: resT,resS
    Integer :: i

    If (Present(iE_put)) Then
        E_index = Max(iE_put,2)
    Else
        E_index = Min(Bisection_Search(E,CS%E_uni,CS%n_E_uni),CS%n_E_uni)
    End If
    sS = 0._dp
    sA = 0._dp
    Do i = 1,CS%n_iso
        !resonant contribution
        If (CS%has_res_cs(i)) Then
            Call sig_Resonance(CS%res_cs(i),E,resT,resS)
            sS = sS + CS%iso_Fractions(i) * resS
            sA = sA + CS%iso_fractions(i) * (resT - resS)
        End If
        !background contribution from scatter
        sS = sS + CS%iso_Fractions(i) * sig_Composite( E,                   & 
                                                     & CS%n_E_uni,          & 
                                                     & CS%E_uni,            & 
                                                     & CS%lnE_uni,          & 
                                                     & E_index,             & 
                                                     & 0,                   & 
                                                     & CS%lev_cs(i)%n_lev,  & 
                                                     & CS%lev_cs(i)%thresh, & 
                                                     & CS%lev_cs(i)%sig     ) 
        !background contribution from absorption
        sA = sA + CS%iso_Fractions(i) * sig_Composite( E,                    & 
                                                     & CS%n_E_uni,           & 
                                                     & CS%E_uni,             & 
                                                     & CS%lnE_uni,           & 
                                                     & E_index,              & 
                                                     & 1,                    & 
                                                     & CS%abs_cs(i)%n_modes, & 
                                                     & CS%abs_cs(i)%thresh,  & 
                                                     & CS%abs_cs(i)%sig      )
    End Do
    sT = sA + sS
    If (Present(iE_get)) iE_get = E_index
End Function sig_T_all

Function sig_T_iso(CS,iso,E,iE_get,iE_put) Result(sT)
!returns total cross section (sig_T_iso) for a single isotope (iso)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sT
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Real(dp) :: resT,resS

    If (Present(iE_put)) Then
        E_index = Max(iE_put,2)
    Else
        E_index = Min(Bisection_Search(E,CS%E_uni,CS%n_E_uni),CS%n_E_uni)
    End If
    sT = 0._dp
    !resonant contribution
    If (CS%has_res_cs(iso)) Then
        Call sig_Resonance(CS%res_cs(iso),E,resT,resS)
        sT = sT + resT
    End If
    !background contribution from scatter
    sT = sT + sig_Composite( E,                     & 
                           & CS%n_E_uni,            & 
                           & CS%E_uni,              & 
                           & CS%lnE_uni,            & 
                           & E_index,               & 
                           & 0,                     & 
                           & CS%lev_cs(iso)%n_lev,  & 
                           & CS%lev_cs(iso)%thresh, & 
                           & CS%lev_cs(iso)%sig     )
    !background contribution from absorption
    sT = sT + sig_Composite( E,                      & 
                           & CS%n_E_uni,             & 
                           & CS%E_uni,               & 
                           & CS%lnE_uni,             & 
                           & E_index,                & 
                           & 1,                      & 
                           & CS%abs_cs(iso)%n_modes, & 
                           & CS%abs_cs(iso)%thresh,  & 
                           & CS%abs_cs(iso)%sig      )
    If (Present(iE_get)) iE_get = E_index
End Function sig_T_iso

Function sig_S_all(CS,E,iE_get,iE_put) Result(sS)
!returns scattering cross section (sig_S_all) for the total atmosphere
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sS
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Integer :: i
    Real(dp) :: resT,resS

    If (Present(iE_put)) Then
        E_index = Max(iE_put,2)
    Else
        E_index = Min(Bisection_Search(E,CS%E_uni,CS%n_E_uni),CS%n_E_uni)
    End If
    sS = 0._dp
    Do i = 1,CS%n_iso
        !resonant contribution
        If (CS%has_res_cs(i)) Then
            Call sig_Resonance(CS%res_cs(i),E,resT,resS)
            sS = sS + CS%iso_Fractions(i) * resS
        End If
        !background contribution
        sS = sS + CS%iso_Fractions(i) * sig_Composite( E,                   & 
                                                     & CS%n_E_uni,          & 
                                                     & CS%E_uni,            & 
                                                     & CS%lnE_uni,          & 
                                                     & E_index,             & 
                                                     & 0,                   & 
                                                     & CS%lev_cs(i)%n_lev,  & 
                                                     & CS%lev_cs(i)%thresh, & 
                                                     & CS%lev_cs(i)%sig     )
    End Do
    If (Present(iE_get)) iE_get = E_index
End Function sig_S_all

Function sig_S_iso(CS,iso,E,iE_get,iE_put) Result(sS)
!returns scattering cross section (sig_S_iso) for a single isotope (iso)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sS
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Real(dp) :: resT,resS

    If (Present(iE_put)) Then
        E_index = Max(iE_put,2)
    Else
        E_index = Min(Bisection_Search(E,CS%E_uni,CS%n_E_uni),CS%n_E_uni)
    End If
    sS = 0._dp
    !resonant contribution
    If (CS%has_res_cs(iso)) Then
        Call sig_Resonance(CS%res_cs(iso),E,resT,resS)
        sS = sS + resS
    End If
    !background contribution
    sS = sS + sig_Composite( E,                     & 
                           & CS%n_E_uni,            & 
                           & CS%E_uni,              & 
                           & CS%lnE_uni,            & 
                           & E_index,               & 
                           & 0,                     & 
                           & CS%lev_cs(iso)%n_lev,  & 
                           & CS%lev_cs(iso)%thresh, & 
                           & CS%lev_cs(iso)%sig     )
    If (Present(iE_get)) iE_get = E_index
End Function sig_S_iso

Function sig_A_all(CS,E,iE_get,iE_put) Result(sA)
!returns absorption cross section (sig_A_all) for the total atmosphere
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sA
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Integer :: i
    Real(dp) :: resT,resS

    If (Present(iE_put)) Then
        E_index = Max(iE_put,2)
    Else
        E_index = Min(Bisection_Search(E,CS%E_uni,CS%n_E_uni),CS%n_E_uni)
    End If
    sA = 0._dp
    Do i = 1,CS%n_iso
        !resonant contribution
        If (CS%has_res_cs(i)) Then
            Call sig_Resonance(CS%res_cs(i),E,resT,resS)
            sA = sA + CS%iso_Fractions(i) * (resT - resS)
        End If
        !background contribution
        sA = sA + CS%iso_Fractions(i) * sig_Composite( E,                    & 
                                                     & CS%n_E_uni,           & 
                                                     & CS%E_uni,             & 
                                                     & CS%lnE_uni,           & 
                                                     & E_index,              & 
                                                     & 1,                    & 
                                                     & CS%abs_cs(i)%n_modes, & 
                                                     & CS%abs_cs(i)%thresh,  & 
                                                     & CS%abs_cs(i)%sig      )
    End Do
    If (Present(iE_get)) iE_get = E_index
End Function sig_A_all

Function sig_A_iso(CS,iso,E,iE_get,iE_put) Result(sA)
!returns absorption cross section (sig_A_iso) for a single isotope (iso)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sA
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Real(dp) :: resT,resS

    If (Present(iE_put)) Then
        E_index = Max(iE_put,2)
    Else
        E_index = Min(Bisection_Search(E,CS%E_uni,CS%n_E_uni),CS%n_E_uni)
    End If
    sA = 0._dp
    !resonant contribution
    If (CS%has_res_cs(iso)) Then
        Call sig_Resonance(CS%res_cs(iso),E,resT,resS)
        sA = sA + resT - resS
    End If
    !background contribution
    sA = sA + sig_Composite( E,                       & 
                            & CS%n_E_uni,             & 
                            & CS%E_uni,               & 
                            & CS%lnE_uni,             & 
                            & E_index,                & 
                            & 1,                      & 
                            & CS%abs_cs(iso)%n_modes, & 
                            & CS%abs_cs(iso)%thresh,  & 
                            & CS%abs_cs(iso)%sig      )
    If (Present(iE_get)) iE_get = E_index
End Function sig_A_iso

Pure Subroutine sig_Resonance(res,E,sT,sS)
    Use Kinds, Only: dp
    Use Global, Only: Pi
    Use Global, Only: TwoPi
    Implicit None
    Type(res_sig_Type), Intent(In) :: res
    Real(dp), Intent(In) :: E ![keV]
    Real(dp), Intent(Out) :: sT
    Real(dp), Intent(Out) :: sS
    Real(dp) :: E_eV  ![eV] energy in eV for local calcs
    Real(dp) :: k  ![1/m] neutron wave number
    Real(dp) :: sA
    Real(dp) :: Sh,Pe,phi
    Real(dp) :: Sin_sq_phi,One_m_Cos2phi,Sin_2phi
    Complex(dp) :: Exp_m2iphi
    Complex(dp) :: ImK,rho,U
    Integer :: l,s,j,r
    Real(dp) :: gJ
    Real(dp) :: Gn
    Real(dp) :: Epr
    Real(dp) :: Gxr
    Real(dp) :: Gt
    Real(dp) :: c
    Real(dp) :: Smin,Smax,Sval
    Real(dp) :: Jmin,Jmax,Jval
    Integer :: nS,nJ
    Complex(dp), Parameter :: z0 = CMPLX(0._dp,Kind=dp)
    Complex(dp), Parameter :: z1 = CMPLX(1._dp,Kind=dp)

    sT = 0._dp
    sS = 0._dp
    sA = 0._dp
    !check if energy is in resonance range
    If (E.LT.res%E_range(1) .OR. E.GT.res%E_range(2)) Return
    E_eV = 1000._dp * E !local energy units in eV
    k = res%k0 * Sqrt(E_eV)
    Do l = 1,res%n_L !Performs the sum over levels for Brookhaven 2018:
                     !  eq D.1 and D.3 for SLBW
                     !  eq D.19 and D.3 for MLBW
                     !  eq D.30 and D.32 for RM
    !Compute Shift, Penetrability, and phase shift
    If (res%L(l)%a .EQ. res%L(l)%ap) Then
        !Sh,Pe,phi are computed from recursive forms from Brookhaven 2018 Table D.1
        Call PSphi(l-1,k*res%L(l)%a,Pe,Sh,phi)
    Else
        !Sh,Pe are computed from recursive forms from Brookhaven 2018 Table D.1
        Call PSphi(l-1,k*res%L(l)%a,Pe,Sh)
        !phi is computed from recursive forms from Brookhaven 2018 Table D.1 w/ a different value of rho
        phi = Phi_hard(l-1,k*res%L(l)%ap)
    End If
    !set up loops and precompute/initialize values based on formalism
    If (res%is_RM) Then
        !sums run through s AND J
        Smin = Abs(res%spi - 0.5_dp)
        Smax = res%spi + 0.5_dp
        nS = NINT(Smax-Smin) + 1
        !Initial and precomuted values depend formalism
        Exp_m2iphi = Exp( CMPLX(0._dp,-2._dp*phi,Kind=dp) )
        One_m_Cos2phi = 1._dp - Cos(2._dp*phi)
    Else !If (res%is_BW) Then
        !sums run through J
        Jmin = Abs(Abs(res%spi - Real(l,dp)) - 0.5_dp)
        Jmax = res%spi + Real(l,dp) + 0.5_dp
        nJ = NINT(Jmax-Jmin) + 1
        !Initial and precomuted values depend formalism
        If (res%is_MLBW) Then
            Exp_m2iphi = Exp( CMPLX(0._dp,-2._dp*phi,Kind=dp) )
        Else !If (res%is_SLBW) Then
            Sin_sq_phi = Sin(phi)**2
            Sin_2phi = Sin(2._dp*phi)
            sS = sS + Real(4*(2*(l-1)+1),dp) * Sin_sq_phi !first term of Brookhaven 2018 eq D.2
        End If
    End If
    If (res%is_RM) Then
        Do s = 1,nS
            Sval = Smin + Real(s-1,dp)
            Jmin = Abs( Real(l-1,dp) - Sval )
            Jmax = Real(l-1,dp) + Sval
            nJ = NINT(Jmax-Jmin) + 1
            Do j = 1,nJ
                Jval = Jmin + Real(j-1,dp)
                gJ = (Jval + 0.5_dp) / (2._dp * res%spi + 1._dp)
                ImK = z0
                Do r = 1,res%L(l)%n_r
                    If (res%L(l)%EraG(r,2) .NE. Jval) Cycle
                    !ASSOCIATE construct creates short-hand variable names for the duration of the block
                    ASSOCIATE ( Er  => res%L(l)%EraG(r,1), & !resonance energy
                              & Gnr => res%L(l)%EraG(r,3), & !neutron width at resonance energy
                              & Ggr => res%L(l)%EraG(r,4), & !radiation width
                              & Pr  => res%L(l)%EraG(r,6)  ) !penetration factor at resonance energy
                    Gn = Pe * Gnr / Pr !neutron width, Brookhaven 2018 eq D.7
                    ImK = ImK + & !sum term for Brookhaven 2018 eq D.28
                                & CMPLX(0._dp,0.5_dp*Gn,Kind=dp) / ( CMPLX(Er-E_eV,Kind=dp) - CMPLX(0._dp,0.5_dp*Ggr,Kind=dp) ) 
                    End ASSOCIATE
                End Do !r
                Imk = z1 - ImK !Brookhaven 2018 eq D.28
                !HACK The following implementation of Brookhaven eqs D.29 - D.31 is caoutioned against due to numerical stability
                rho = z1 - z1/ImK !Brookhaven 2018 eq D.29
                c = One_m_Cos2phi + 2._dp * Real(Exp_m2iphi * rho) !precomputed quantity
                sT = sT + gJ * c !Brookhaven 2018 eq D.30
                sS = sS + gJ * ( c + 2._dp * (Abs(rho)**2 - Real(rho)) ) !Brookhaven 2018 eq D.31
            End Do !j
        End Do !s
    Else If (res%is_BW) Then
        Do j = 1,nJ
            Jval = Jmin + Real(j-1,dp)
            gJ = (Jval + 0.5_dp) / (2._dp * res%spi + 1._dp)
            Do r = 1,res%L(l)%n_r
                If (res%L(l)%EraG(r,2) .NE. Jval) Cycle
                !ASSOCIATE construct creates short-hand variable names for the duration of the block
                ASSOCIATE ( Er  => res%L(l)%EraG(r,1), & !resonance energy
                          & Gtr => res%L(l)%EraG(r,3), & !total width at resonance energy
                          & Gnr => res%L(l)%EraG(r,4), & !neutron width at resonance energy
                          & Ggr => res%L(l)%EraG(r,5), & !radiation width
                          & Sr  => res%L(l)%EraG(r,6), & !shift factor at resonance energy
                          & Pr  => res%L(l)%EraG(r,7)  ) !penetration factor at resonance energy
                Gn = Pe * Gnr / Pr !neutron width, Brookhaven 2018 eq D.7
                Epr = Er + Gnr * (Sr - Sh) / (2._dp * Pr) !Primed resonance energy, Brookhaven 2018 eq D.9
                Gxr = Gtr - Gnr - Ggr !competitive width, Brookhaven 2018 eq D.8
                Gt = Gn + Ggr + Gxr !total width, Brookhaven 2018 eq p.340 (between D.7 and D.8)
                c = gJ * Gn / ((E_eV - Epr)**2 + 0.25_dp * Gt**2) !precomputed value
                sA = sA + c * Ggr !sum term for Brookhaven 2018 eq D.4
                If (res%is_SLBW) Then
                    sS = sS + c * ( Gn - 2._dp*(Gt*Sin_sq_phi + (E_eV-Epr)*Sin_2phi) ) !sum term for Brookhaven 2018 eq D.2
                Else !if(res%is_MLBW) Then
                    U = Exp_m2iphi * & !Brookhaven 2018 eq D.21
                        & ( z1 + CMPLX(0._dp,Gn,Kind=dp) / (CMPLX(Epr-E_eV,Kind=dp) - CMPLX(0._dp,0.5_dp*Gt,Kind=dp)) )
                    sS = sS + gJ * Abs(z1 - U)**2 !sum term for Brookhaven 2018 eq D.20
                End If
                End ASSOCIATE
            End Do !r
        End Do !j
    End If
    End Do !l
    !finish multiplying in scaling factors
    If (res%is_BW) Then
        If (res%is_SLBW) Then
            If (sS .LE. 0._dp) sS = Tiny(sS) !if SLBW returns negative cross section, set to tiny positive value instead
        End If
        sS = Pi * sS / k**2
        sA = Pi * sA / k**2
        sT = sS + sA
    Else !If (res%is_RM) Then
        sT = TwoPi * sT / k**2
        sS = TwoPi * sS / k**2
    End If
End Subroutine sig_Resonance

Pure Function sig_Composite(E,n_E,E_list,lnE_list,E_index,n1,n2,t_list,sig_list) Result(sig)
    Use Kinds, Only: dp
    Use Interpolation, Only: Linear_Interp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Real(dp) :: sig
    Real(dp), Intent(In) :: E
    Integer, Intent(In) :: n_E
    Real(dp), Intent(In) :: E_list(1:n_E)
    Real(dp), Intent(In) :: lnE_list(1:n_E)
    Integer, Intent(In) :: E_index
    Integer, Intent(In) :: n1,n2
    Integer, Intent(In) :: t_list(n1:n2)
    Type(sig_Type), Intent(In) :: sig_list(n1:n2)

    sig = Sum( sig_composite_list(E,n_E,E_list,lnE_list,E_index,n1,n2,t_list,sig_list) )
    ! Integer :: i,j
    ! Integer :: index1,index2
    ! Real(dp) :: E1,E2
    ! Real(dp) :: sig1,sig2,sig_n
    ! Integer :: method
    ! Integer, Parameter :: Hist_interpolation   = 1  !y is constant in x (constant, histogram)
    ! Integer, Parameter :: LinLin_interpolation = 2  !y is linear in x (linear-linear)
    ! Integer, Parameter :: LinLog_interpolation = 3  !y is linear in ln(x) (linear-log)
    ! Integer, Parameter :: LogLin_interpolation = 4  !ln(y) is linear in x (log-linear)
    ! Integer, Parameter :: LogLog_interpolation = 5  !ln(y) is linear in ln(x) (log-log)

    ! sig = 0._dp
    ! method = LinLin_interpolation
    ! i = n1 !initialize loop counter: loop runs at least once, from n1 to n2
    ! Do
    !     If (E_index .LE. t_list(i)) Exit
    !     index1 = sig_list(i)%E_map(E_index) - 1
    !     index2 = sig_list(i)%E_map(E_index)
    !     Do j = 1,sig_list(i)%n_interp_r
    !         If (index2 .LE. sig_list(i)%interp(j,1)) Then
    !             method = sig_list(i)%interp(j,2)
    !             Exit
    !         End If
    !     End Do
    !     Select Case (method)
    !         Case(Hist_interpolation)
    !             sig_n = sig_list(i)%sig( index1 )
    !         Case(LinLin_interpolation)
    !             E1 = E_list( sig_list(i)%E_key( index1 ) )
    !             E2 = E_list( sig_list(i)%E_key( index2 ) )
    !             sig1 = sig_list(i)%sig( index1 )
    !             sig2 = sig_list(i)%sig( index2 )
    !             sig_n = Linear_Interp(E,E1,E2,sig1,sig2)
    !         Case(LinLog_interpolation)
    !             E1 = lnE_list( sig_list(i)%E_key( index1 ) )
    !             E2 = lnE_list( sig_list(i)%E_key( index2 ) )
    !             sig1 = sig_list(i)%sig( index1 )
    !             sig2 = sig_list(i)%sig( index2 )
    !             sig_n = Linear_Interp(Log(E),E1,E2,sig1,sig2)
    !         Case(LogLin_interpolation)
    !             E1 = E_list( sig_list(i)%E_key( index1 ) )
    !             E2 = E_list( sig_list(i)%E_key( index2 ) )
    !             sig1 = sig_list(i)%lnsig( index1 )
    !             sig2 = sig_list(i)%lnsig( index2 )
    !             sig_n = Exp(Linear_Interp(E,E1,E2,sig1,sig2))
    !         Case(LogLog_interpolation)
    !             E1 = lnE_list( sig_list(i)%E_key( index1 ) )
    !             E2 = lnE_list( sig_list(i)%E_key( index2 ) )
    !             sig1 = sig_list(i)%lnsig( index1 )
    !             sig2 = sig_list(i)%lnsig( index2 )
    !             sig_n = Exp(Linear_Interp(Log(E),E1,E2,sig1,sig2))
    !     End Select
    !     If (sig_n .GT. 0._dp) sig = sig + sig_n
    !     i = i + 1
    !     If (i .GT. n2) Exit
    ! End Do
End Function sig_Composite

Pure Function sig_Composite_list(E,n_E,E_list,lnE_list,E_index,n1,n2,t_list,sig_list) Result(sig)
    Use Kinds, Only: dp
    Use Interpolation, Only: Linear_Interp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Real(dp), Intent(In) :: E
    Integer, Intent(In) :: n_E
    Real(dp), Intent(In) :: E_list(1:n_E)
    Real(dp), Intent(In) :: lnE_list(1:n_E)
    Integer, Intent(In) :: E_index
    Integer, Intent(In) :: n1,n2
    Integer, Intent(In) :: t_list(n1:n2)
    Type(sig_Type), Intent(In) :: sig_list(n1:n2)
    Real(dp) :: sig(n1:n2)
    Integer :: i,j
    Integer :: index1,index2
    Real(dp) :: E1,E2
    Real(dp) :: sig1,sig2,sig_n
    Integer :: method
    Integer, Parameter :: Hist_interpolation   = 1  !y is constant in x (constant, histogram)
    Integer, Parameter :: LinLin_interpolation = 2  !y is linear in x (linear-linear)
    Integer, Parameter :: LinLog_interpolation = 3  !y is linear in ln(x) (linear-log)
    Integer, Parameter :: LogLin_interpolation = 4  !ln(y) is linear in x (log-linear)
    Integer, Parameter :: LogLog_interpolation = 5  !ln(y) is linear in ln(x) (log-log)

    sig = 0._dp
    method = LinLin_interpolation
    Do i = n1,n2 !loop runs at least once, from n1 to n2
        If (E_index .LE. t_list(i)) Exit
        sig_n = -1._dp
        index1 = sig_list(i)%E_map(E_index) - 1
        index2 = sig_list(i)%E_map(E_index)
        Do j = 1,sig_list(i)%n_interp_r
            If (index2 .LE. sig_list(i)%interp(j,1)) Then
                method = sig_list(i)%interp(j,2)
                Exit
            End If
        End Do
        Select Case (method)
            Case(Hist_interpolation)
                sig_n = sig_list(i)%sig( index1 )
            Case(LinLin_interpolation)
                E1 = E_list( sig_list(i)%E_key( index1 ) )
                E2 = E_list( sig_list(i)%E_key( index2 ) )
                sig1 = sig_list(i)%sig( index1 )
                sig2 = sig_list(i)%sig( index2 )
                sig_n = Linear_Interp(E,E1,E2,sig1,sig2)
            Case(LinLog_interpolation)
                E1 = lnE_list( sig_list(i)%E_key( index1 ) )
                E2 = lnE_list( sig_list(i)%E_key( index2 ) )
                sig1 = sig_list(i)%sig( index1 )
                sig2 = sig_list(i)%sig( index2 )
                sig_n = Linear_Interp(Log(E),E1,E2,sig1,sig2)
            Case(LogLin_interpolation)
                E1 = E_list( sig_list(i)%E_key( index1 ) )
                E2 = E_list( sig_list(i)%E_key( index2 ) )
                sig1 = sig_list(i)%lnsig( index1 )
                sig2 = sig_list(i)%lnsig( index2 )
                sig_n = Exp(Linear_Interp(E,E1,E2,sig1,sig2))
            Case(LogLog_interpolation)
                E1 = lnE_list( sig_list(i)%E_key( index1 ) )
                E2 = lnE_list( sig_list(i)%E_key( index2 ) )
                sig1 = sig_list(i)%lnsig( index1 )
                sig2 = sig_list(i)%lnsig( index2 )
                sig_n = Exp(Linear_Interp(Log(E),E1,E2,sig1,sig2))
        End Select
        sig(i) = Max(sig_n,0._dp)
    End Do
End Function sig_Composite_list

Pure Subroutine Broad_sig_start(E,M,T,v,gamma,vRmin,vRmax)
    Use Kinds, Only: dp
    Use Global, Only: k_Boltzmann
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp), Intent(In) :: E,M,T
    Real(dp), Intent(Out) :: v,gamma,vRmin,vRmax
    Real(dp) :: vT
    Real(dp), Parameter :: twok_B = 2._dp * k_Boltzmann / 1000._dp**2  ![J/K]*[1 km^2 / 1000^2 m^2] Two multiplied by the Boltzmann 
                                                                       !constant with convenient units locally

    v = Neutron_Speed(E)
    gamma = Sqrt(M / (twok_b * T))
    vT = 4._dp / gamma  !same cutoff as SIGMA1 algorithm used by NJOY
    If (vT .LE. v) Then
        vRmin = v - vT
    Else
        !NOTE: Use small VRmin instead of zero to prevent Log(E=0) from appearing in log interpolation cases... 1.E-16 equates to 
        !an energy of about 1.E-38 keV
        vRmin = 1.E-16_dp  !0._dp
    End If
    vRmax = v + vT
End Subroutine Broad_sig_start

Subroutine sig_T_A_db_all(CS,E,T,sigT,sigA)
    Use Kinds, Only: dp
    Implicit None
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]
    Real(dp), Intent(Out) :: sigT,sigA

    Call sig_T_A_db(CS,E,T,sigT,sigA)
End Subroutine sig_T_A_db_all

Subroutine sig_T_A_db_iso(CS,iso,E,T,sigT,sigA)
    Use Kinds, Only: dp
    Implicit None
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]
    Real(dp), Intent(Out) :: sigT,sigA

    Call sig_T_A_db(CS,E,T,sigT,sigA,iso)
End Subroutine sig_T_A_db_iso

Subroutine sig_T_A_db(CS,E,T,sigT,sigA,iso)
    Use Kinds, Only: dp
    Use Global, Only: SqrtPi
    Use Utilities, Only: Bisection_Search
    Use Neutron_Utilities, Only: Neutron_Energy
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Type(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]
    Real(dp), Intent(Out) :: sigT,sigA
    Integer, Intent(In), Optional :: iso
    Real(dp) :: v  ![km/s] neutron velocity
    Real(dp) :: gamma  ![km/s]
    Real(dp) :: vR_min,vR_max
    Integer :: iE_min,iE_max
    Real(dp) :: sig_T_A(1:2)
    Integer :: i
    Integer, Parameter :: total = 1
    Integer, Parameter :: absor = 2

    Call Broad_sig_start(E,CS%Mn,T,v,gamma,vR_min,vR_max)
    iE_min = Bisection_Search(Neutron_Energy(vR_min),CS%E_uni,CS%n_E_uni)
    If (Neutron_Energy(vR_max) .LE. CS%E_uni(iE_min)) Then  !only a single interval spanned in energy grid
        If (Present(iso)) Then
            sig_T_A = db_Romberg_T_A(CS,iE_min,vR_min,vR_max,gamma,v,iso) * gamma / (SqrtPi * v**2)
        Else
            sig_T_A = db_Romberg_T_A(CS,iE_min,vR_min,vR_max,gamma,v) * gamma / (SqrtPi * v**2)
        End If
        sigT = sig_T_A(total)
        sigA = sig_T_A(absor)
        Return
    End If
    iE_max = Bisection_Search(Neutron_Energy(vR_max),CS%E_uni,CS%n_E_uni) - 1
    If (Present(iso)) Then
        !first interval (from E_min to first indexed E) is partial
        sig_T_A = db_Romberg_T_A(CS,iE_min,vR_min,Neutron_Speed(CS%E_uni(iE_min)),gamma,v,iso)
        !middle intervals
        Do i = iE_min,iE_max-1
            sig_T_A = sig_T_A + db_Romberg_T_A(CS,i+1,Neutron_Speed(CS%E_uni(i)),Neutron_Speed(CS%E_uni(i+1)),gamma,v,iso)
        End Do
        !last interval (from second to last indexed E to Emax) is partial, also apply normalization
        sig_T_A = sig_T_A + db_Romberg_T_A(CS,iE_max,Neutron_Speed(CS%E_uni(iE_max)),vR_max,gamma,v,iso)
    Else
        !first interval (from E_min to first indexed E) is partial
        sig_T_A = db_Romberg_T_A(CS,iE_min,vR_min,Neutron_Speed(CS%E_uni(iE_min)),gamma,v)
        !middle intervals
        Do i = iE_min,iE_max-1
            sig_T_A = sig_T_A + db_Romberg_T_A(CS,i+1,Neutron_Speed(CS%E_uni(i)),Neutron_Speed(CS%E_uni(i+1)),gamma,v)
        End Do
        !last interval (from second to last indexed E to Emax) is partial, also apply normalization
        sig_T_A = sig_T_A + db_Romberg_T_A(CS,iE_max,Neutron_Speed(CS%E_uni(iE_max)),vR_max,gamma,v)
    End If
    sig_T_A = sig_T_A * gamma / (SqrtPi * v**2)
    sigT = sig_T_A(total)
    sigA = sig_T_A(absor)
End Subroutine sig_T_A_db

Function db_Romberg_T_A(CS,iE,vR1,vR2,gamma,v,iso) Result(sig_T_A)
    Use Kinds, Only: dp
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Real(dp) :: sig_T_A(1:2)
    Type(CS_type), Intent(In) :: CS
    Integer, Intent(In) :: iE
    Real(dp), Intent(In) :: vR1,vR2
    Real(dp), Intent(In) :: gamma,v
    Integer, Intent(In), Optional :: iso
    Real(dp) :: sig_vR1(1:2),sig_vR2(1:2),sig_vR(1:2)
    Real(dp) :: s1(1:2),s2(1:2)
    Real(dp) :: vR
    Integer, Parameter :: total = 1
    Integer, Parameter :: absor = 2
    Real(dp) :: T(1:2,0:Tmax)  !Extrapolation table previous row
    Real(dp) :: Tk0(1:2),Tk(1:2)  !Extrapolation table current row values
    Integer :: i,j,k  !counters: i for table row, j for quadrature ordinates, k for table column
    Integer :: n      !number of intervals
    Real(dp) :: h0,h  !spacing between quadrature ordinates
    Real(dp) :: fk    !multiplier for extrapolation steps

    !Initial trapezoid estimate
    If (Present(iso)) Then
        Call CS%sig_T_A(iso,Neutron_Energy(vR1),sig_vR1(total),sig_vR1(absor),iE_put=iE)
        Call CS%sig_T_A(iso,Neutron_Energy(vR2),sig_vR2(total),sig_vR2(absor),iE_put=iE)
    Else
        Call CS%sig_T_A(Neutron_Energy(vR1),sig_vR1(total),sig_vR1(absor),iE_put=iE)
        Call CS%sig_T_A(Neutron_Energy(vR2),sig_vR2(total),sig_vR2(absor),iE_put=iE)
    End If
    s1 = 0.5_dp * (Broad_Integrand( vR1,sig_vR1,gamma,v) + Broad_Integrand( vR2,sig_vR2,gamma,v))
    s2 = 0.5_dp * (Broad_Integrand(-vR1,sig_vR1,gamma,v) + Broad_Integrand(-vR2,sig_vR2,gamma,v))
    h0 = vR2 - vR1
    T(:,0) = h0 * (s1 - s2)
    n = 1
    Do i = 1,Tmax !up to Tmax rows in the table
        !Trapezoid estimate for the 0-th column of the i-th row of table
        n = n * 2
        h = h0 / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate the integrand
            vR = vR1 + Real(j,dp)*h
            If (Present(Iso)) Then
                Call CS%sig_T_A(iso,Neutron_Energy(vR),sig_vR(total),sig_vR(absor),iE_put=iE)
            Else
                Call CS%sig_T_A(Neutron_Energy(vR),sig_vR(total),sig_vR(absor),iE_put=iE)
            End If
            s1 = s1 + Broad_Integrand( vR,sig_vR,gamma,v)
            s2 = s2 + Broad_Integrand(-vR,sig_vR,gamma,v)
        End Do
        Tk0 = h * (s1 - s2)
        !Fill i-th row, columns k = 1:i, with extrapolated estimates
        fk = 1._dp
        Do k = 1,i  !up to i columns this row
            fk = fk * 4._dp
            Tk = (fk * Tk0 - T(:,k-1)) / (fk - 1._dp)
            If (k .LT. i) Then
                T(:,k-1) = Tk0  !store Tk0 for next i
                Tk0 = Tk  !store Tk for next k
            End If !otherwise, skip storage steps if working final column
        End Do
        !check for convergence
        If ( All( Abs(T(:,i-1) - Tk) .LE. rTol * Abs(Tk) ) ) Then
            sig_T_A = Tk  !Tk is the highest precision converged value
            Return  !Normal exit
        Else  !prep for the next time though the loop
            !store Tk0 and Tk for next i
            T(:,i-1) = Tk0
            T(:,i) = Tk
        End If
    End Do
    !If we get this far, we did not converge
    Write(*,*)
    Write(*,'(A,I0,A)')        'ERROR:  n_Cross_sections: db_Romberg_T_A:  Failed to converge in ',Tmax,' extrapolations.'
    Write(*,'(A,2ES23.15)')    '        Final estimated value: ',Tk
    Write(*,'(A,2ES23.15)')    '        Prior estimated value: ',Tk0
    ERROR STOP
End Function db_Romberg_T_A

Function sig_T_all_db(CS,E,T) Result(sigT)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: sigT
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]

    sigT = sig_db(CS,E,T,db_choose_total_sig)
End Function sig_T_all_db

Function sig_T_iso_db(CS,iso,E,T) Result(sigT)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: sigT
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]

    sigT = sig_db(CS,E,T,db_choose_total_sig,iso)
End Function sig_T_iso_db

Function sig_S_all_db(CS,E,T) Result(sigS)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: sigS
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]

    sigS = sig_db(CS,E,T,db_choose_scatter_sig)
End Function sig_S_all_db

Function sig_S_iso_db(CS,iso,E,T) Result(sigS)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: sigS
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]

    sigS = sig_db(CS,E,T,db_choose_scatter_sig,iso)
End Function sig_S_iso_db

Function sig_A_all_db(CS,E,T) Result(sigA)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: sigA
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]

    sigA = sig_db(CS,E,T,db_choose_absorp_sig)
End Function sig_A_all_db

Function sig_A_iso_db(CS,iso,E,T) Result(sigA)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: sigA
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]

    sigA = sig_db(CS,E,T,db_choose_absorp_sig,iso)
End Function sig_A_iso_db

Function sig_db(CS,E,T,sig_type,iso) Result(sig)
    Use Kinds, Only: dp
    Use Global, Only: SqrtPi
    Use Utilities, Only: Bisection_Search
    Use Neutron_Utilities, Only: Neutron_Energy
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp) :: sig
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]
    Integer, Intent(In) :: sig_type
    Integer, Intent(In), Optional :: iso
    Real(dp) :: v  ![km/s] neutron velocity
    Real(dp) :: gamma  ![km/s]
    Real(dp) :: vR_min,vR_max
    Integer :: iE_min,iE_max
    Integer :: i

    Call Broad_sig_start(E,CS%Mn,T,v,gamma,vR_min,vR_max)
    iE_min = Bisection_Search(Neutron_Energy(vR_min),CS%E_uni,CS%n_E_uni)
    If (Neutron_Energy(vR_max) .LE. CS%E_uni(iE_min)) Then  !only a single interval spanned in energy grid
        If (Present(iso)) Then
            sig = db_Romberg(CS,iE_min,vR_min,vR_max,gamma,v,sig_type,iso) *  gamma / (SqrtPi * v**2)
        Else
            sig = db_Romberg(CS,iE_min,vR_min,vR_max,gamma,v,sig_type) *  gamma / (SqrtPi * v**2)
        End If
        Return
    End If
    iE_max = Bisection_Search(Neutron_Energy(vR_max),CS%E_uni,CS%n_E_uni) - 1
    If (Present(iso)) Then
        !first interval (from E_min to first indexed E) is partial
        sig = db_Romberg(CS,iE_min,vR_min,Neutron_Speed(CS%E_uni(iE_min)),gamma,v,sig_type,iso)
        !middle intervals
        Do i = iE_min,iE_max-1
            sig = sig + db_Romberg(CS,i+1,Neutron_Speed(CS%E_uni(i)),Neutron_Speed(CS%E_uni(i+1)),gamma,v,sig_type,iso)
        End Do
        !last interval (from second to last indexed E to Emax) is partial
        sig = sig + db_Romberg(CS,iE_max,Neutron_Speed(CS%E_uni(iE_max)),vR_max,gamma,v,sig_type,iso)
    Else
        !first interval (from E_min to first indexed E) is partial
        sig = db_Romberg(CS,iE_min,vR_min,Neutron_Speed(CS%E_uni(iE_min)),gamma,v,sig_type)
        !middle intervals
        Do i = iE_min,iE_max-1
            sig = sig + db_Romberg(CS,i+1,Neutron_Speed(CS%E_uni(i)),Neutron_Speed(CS%E_uni(i+1)),gamma,v,sig_type)
        End Do
        !last interval (from second to last indexed E to Emax) is partial
        sig = sig + db_Romberg(CS,iE_max,Neutron_Speed(CS%E_uni(iE_max)),vR_max,gamma,v,sig_type)
    End If
    sig = sig * gamma / (SqrtPi * v**2)
End Function sig_db

Function db_Romberg(CS,iE,vR1,vR2,gamma,v,sig_type,iso) Result(sT)
    Use Kinds, Only: dp
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Real(dp) :: sT
    Type(CS_type), Intent(In) :: CS
    Integer, Intent(In) :: iE
    Real(dp), Intent(In) :: vR1,vR2
    Real(dp), Intent(In) :: gamma,v
    Integer, Intent(In) :: sig_type
    Integer, Intent(In), Optional :: iso
    Real(dp) :: sig_vR1,sig_vR2,sig_vR
    Real(dp) :: s1,s2
    Real(dp) :: vR
    Real(dp) :: T(0:Tmax)  !Extrapolation table previous row
    Real(dp) :: Tk0,Tk  !Extrapolation table current row values
    Integer :: i,j,k  !counters: i for table row, j for quadrature ordinates, k for table column
    Integer :: n      !number of intervals
    Real(dp) :: h0,h  !spacing between quadrature ordinates
    Real(dp) :: fk    !multiplier for extrapolation steps

    !Initial trapezoid estimate
    If (sig_type .EQ. db_choose_total_sig) Then
        If (Present(iso)) Then
            sig_vR1 = CS%sig_T(iso,Neutron_Energy(vR1),iE_put=iE)
            sig_vR2 = CS%sig_T(iso,Neutron_Energy(vR2),iE_put=iE)
        Else
            sig_vR1 = CS%sig_T(Neutron_Energy(vR1),iE_put=iE)
            sig_vR2 = CS%sig_T(Neutron_Energy(vR2),iE_put=iE)
        End If
    Else If (sig_type .EQ. db_choose_scatter_sig) Then
        If (Present(iso)) Then
            sig_vR1 = CS%sig_S(iso,Neutron_Energy(vR1),iE_put=iE)
            sig_vR2 = CS%sig_S(iso,Neutron_Energy(vR2),iE_put=iE)
        Else
            sig_vR1 = CS%sig_S(Neutron_Energy(vR1),iE_put=iE)
            sig_vR2 = CS%sig_S(Neutron_Energy(vR2),iE_put=iE)
        End If
    Else If (sig_type .EQ. db_choose_absorp_sig) Then
        If (Present(iso)) Then
            sig_vR1 = CS%sig_A(iso,Neutron_Energy(vR1),iE_put=iE)
            sig_vR2 = CS%sig_A(iso,Neutron_Energy(vR2),iE_put=iE)
        Else
            sig_vR1 = CS%sig_A(Neutron_Energy(vR1),iE_put=iE)
            sig_vR2 = CS%sig_A(Neutron_Energy(vR2),iE_put=iE)
        End If
    End If
    s1 = 0.5_dp * (Broad_Integrand( vR1,sig_vR1,gamma,v) + Broad_Integrand( vR2,sig_vR2,gamma,v))
    s2 = 0.5_dp * (Broad_Integrand(-vR1,sig_vR1,gamma,v) + Broad_Integrand(-vR2,sig_vR2,gamma,v))
    h0 = vR2 - vR1
    T(0) = h0 * (s1 - s2)
    n = 1
    Do i = 1,Tmax !up to Tmax rows in the table
        !Trapezoid estimate for the 0-th column of the i-th row of table
        n = n * 2
        h = h0 / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate the integrand
            vR = vR1 + Real(j,dp)*h
            If (sig_type .EQ. db_choose_total_sig) Then
                If (Present(iso)) Then
                    sig_vR = CS%sig_T(iso,Neutron_Energy(vR),iE_put=iE)
                Else
                    sig_vR = CS%sig_T(Neutron_Energy(vR),iE_put=iE)
                End If
            Else If (sig_type .EQ. db_choose_scatter_sig) Then
                If (Present(iso)) Then
                    sig_vR = CS%sig_S(iso,Neutron_Energy(vR),iE_put=iE)
                Else
                    sig_vR = CS%sig_S(Neutron_Energy(vR),iE_put=iE)
                End If
            Else If (sig_type .EQ. db_choose_absorp_sig) Then
                If (Present(iso)) Then
                    sig_vR = CS%sig_A(iso,Neutron_Energy(vR),iE_put=iE)
                Else
                    sig_vR = CS%sig_A(Neutron_Energy(vR),iE_put=iE)
                End If
            End If
            s1 = s1 + Broad_Integrand( vR,sig_vR,gamma,v)
            s2 = s2 + Broad_Integrand(-vR,sig_vR,gamma,v)
        End Do
        Tk0 = h * (s1 - s2)
        !Fill i-th row, columns k = 1:i, with extrapolated estimates
        fk = 1._dp
        Do k = 1,i  !up to i columns this row
            fk = fk * 4._dp
            Tk = (fk * Tk0 - T(k-1)) / (fk - 1._dp)
            If (k .LT. i) Then
                T(k-1) = Tk0  !store Tk0 for next i
                Tk0 = Tk  !store Tk for next k
            End If !otherwise, skip storage steps if working final column
        End Do
        !check for convergence
        If ( Abs(T(i-1) - Tk) .LE. rTol * Abs(Tk) ) Then
            sT = Tk  !Tk is the highest precision converged value
            Return  !Normal exit
        Else  !prep for the next time though the loop
            !store Tk0 and Tk for next i
            T(i-1) = Tk0
            T(i) = Tk
        End If
    End Do
    !If we get this far, we did not converge
    Write(*,*)
    Write(*,'(A,I0,A)')    'ERROR:  n_Cross_sections: db_Romberg:  Failed to converge in ',Tmax,' extrapolations.'
    Write(*,'(A,ES23.15)') '        Final estimated value: ',Tk
    Write(*,'(A,ES23.15)') '        Prior estimated value: ',Tk0
    ERROR STOP
End Function db_Romberg

Elemental Function Broad_Integrand(vR,sig,gamma,v)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: Broad_Integrand
    Real(dp), Intent(In) :: vR,sig,gamma,v

    Broad_Integrand = vR**2 * sig * Exp( -((gamma * (v-vR))**2) )
End Function Broad_Integrand

Subroutine Write_Cross_Sections(CS,Broadened_CS,file_name)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: half_dash_line
    Implicit None
    Type(CS_Type), Intent(In) :: CS
    Logical, Intent(In) :: Broadened_CS
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat
    Integer :: i
    Real(dp) :: sT,sA,sTb,sAb

    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Then
        Print *,'ERROR:  Cross_Sections: Write_Cross_Sections:  File open error, '//file_name//', IOSTAT=',stat
        ERROR STOP
    End If
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') 'CROSS SECTIONS INFORMATION'
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') '  Microscopic cross sections for total atmosphere:'
    Write(unit,*)
    If (Broadened_CS) Then
        Write(unit,'(5A27)') 'Energy [keV]','Total CS [barns]','Total CS [barns]','Absorption CS [barns]','Absorption CS [barns]'
        Write(unit,'(5A27)') '','(0 deg K)','(273.15 deg K)','(0 deg K)','(273.15 deg K)'
        Write(unit,'(5A27)') '-----------------------','-----------------------','-----------------------', &
                           & '-----------------------','-----------------------'
        Do i = 1,CS%n_E_uni
            Call CS%sig_T_A(CS%E_uni(i),sT,sA)
            Call CS%sig_T_A(CS%E_uni(i),273.15_dp,sTb,sAb)
            Write(unit,'(5ES27.16E3)') CS%E_uni(i),sT,sTb,sA,sAb
        End Do
    Else
        Write(unit,'(3A27)') 'E [keV]','Total CS [barns]','Absorption CS [barns]'
        Write(unit,'(3A27)') '-----------------------','-----------------------','-----------------------'
        Do i = 1,CS%n_E_uni
            Call CS%sig_T_A(CS%E_uni(i),sT,sA)
            Write(unit,'(3ES27.16E3)') CS%E_uni(i),sT,sA
        End Do
    End If
    Write(unit,*)
    Write(unit,*)
    Close(unit)
End Subroutine Write_Cross_Sections

End Module n_Cross_Sections
