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
Module Neutron_Scatter
    
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Implicit None

    Private
    Public :: Neutron_Type
    Public :: Scatter_Data_Type
    Public :: Scatter_Model_Type
    Public :: Setup_Scatter_Model
    Public :: Finish_Scatter_Model_Setup
    Public :: Write_Scatter_Model
    Public :: Air_Velocity
    Public :: Apparent_Energy
    Public :: Scattered_Direction
    Public :: Scattered_Angles

    Type Neutron_Type
        Real(dp) :: r(1:3)  ![km] cartesian coordinates (x,y,z) of neutron
        Real(dp) :: big_r  ![km] magnitude of r
        Real(dp) :: Z  ![km]  geometric altitude
        Real(dp) :: zeta  !cosine of the zenith angle
        Real(dp) :: Omega_hat(1:3)  ![1 km/s] unit vector in direction of neutron travel (unit velocity vector)
        Real(dp) :: E0ef  ![keV]energy at emission in emission frame
        Real(dp) :: s0ef  ![km/s] speed in emission frame at emission
        Real(dp) :: E  ![keV] neutron energy
        Real(dp) :: weight  !weight (probability of neutron existence in this state), reduced for things like absorption supression
        Real(dp) :: t  ![s] simulation time since neutron emission
    End Type

    Type Scatter_Data_Type
        Real(dp) :: sig_A  ![km^-1]  absorption and total macroscopic (seal-level) cross sections
        Real(dp) :: sig_T
        Integer :: target_index  !index of the target isotope
        Real(dp) :: An  !mass of the target isotope in neutron masses
        Real(dp) :: vAir(1:3)  !side effect of computing vA, saved in case all materials/mechanims are computed
        Real(dp) :: vA(1:3)  ![km/s] velocity vector of the target nucleus
        Integer :: level  !sampled level of scatter
        Real(dp) :: Q  !Q-value for chosen scatter level
        Logical :: da_is_iso
        Logical :: da_is_Legendre
        Logical :: da_is_tab
        Integer :: n_a  !number of scatter coeffiecents
        Real(dp), Allocatable :: a(:)  !anisotropic scatter coeffs
        Integer :: n_a1,n_a2  !number of tabulated cosine pdf pairs
        Real(dp) :: a_tab_Econv
        Real(dp), Allocatable :: a_tab1(:,:),a_tab2(:,:)  !tabulated scatter cosine pdf values at energies +/- the incident energy
        Real(dp), Allocatable :: iso_cs(:)  !relative scatter cross section of each isotope at last computed energy
                                            !saved as a side effect in case all materials/mechanims are computed
        Real(dp), Allocatable :: lev_cs(:)  !relative scatter cross section of each level at last computed energy
                                            !saved as a side effect in case all materials/mechanims are computed
        Real(dp) :: Omega_hat0_cm(1:3)
        Real(dp) :: u(1:3),u_speed
        Real(dp) :: s1cm  !speed after the scatter in CM frame is computed as a side effect of the next-event procedure,
                          !saving it reduces operations in following scatter simulation
        Real(dp) :: A_hat(1:3),B_hat(1:3),C_hat(1:3)  !cartesian basis in CM frame computed as a side effect of next-event procedure
                                                      !saving it reduces operations in following scatter simulation
    End Type

    Type Scatter_Model_Type
        Integer :: n_scatters  !number of scatters to follow for history
        Logical :: direct_contribution
        Logical :: estimate_each_scatter
        Logical :: elastic_only
        Logical :: aniso_dist
        Real(dp) :: E_min,E_max
        Logical :: suppress_absorption
        Logical :: suppress_leakage
        Logical :: all_mat_mech
        Logical :: fast_analog
        Logical :: roulette
        Real(dp) :: roulette_weight
        Real(dp) :: roulette_rate
        Real(dp) :: roulette_mult
        Integer(id) :: n_kills(1:6)  !kills due to (1)Time, (2)Energy, (3)Weight, (4)Leakage, (5)Absorption. (6)Collision limit
        Integer(id) :: next_events(1:3)  !next-events at detector (1)Attempted, (2)Found, (3)Tallied
        Integer(id) :: n_no_tally(1:3)  !next-events at detector NOT tallied due to (1)Time & Energy, (2)Time, (3)Energy
        Integer(id) :: n_uncounted  !histories implicitly leaked due to exoatmospheric source geometries
        Type(Scatter_Data_Type) :: scat  !parameters defining the next scatter
        Logical :: Gravity
        Logical :: Neutron_Decay
        !N2H The logical 'Target_Motion' is currently unused except for informational purposes, consider removing
        Logical :: Target_Motion
            Logical :: Doppler_Broaden
            Logical :: Thermal_Motion
            Logical :: Diatomic_Atm
            Logical :: Rotating_Earth
            Logical :: Wind
        Logical :: cs_loaded
        Logical :: inbound_trajectories
    Contains
        Procedure, Pass :: Sample_Scatter
        Procedure, Pass :: Set_Scatter_prep
        Procedure, Pass :: Set_Scatter_iso
        Procedure, Pass :: Set_Scatter_lev
    End Type

Contains

Function Setup_Scatter_Model(setup_file_name,run_file_name,atm_model_i) Result(ScatMod)
    Use Kinds, Only: id
    Use n_Cross_Sections, Only: Setup_Cross_Sections
    Use Global, Only: n_kill_weight
    Use FileIO_Utilities, Only: Worker_Index
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Type(Scatter_Model_Type) :: ScatMod
    Character(*), Intent(In) :: setup_file_name
    Character(*), Intent(In) :: run_file_name
    Integer, Intent(In) :: atm_model_i
    Integer :: n_scatters
    Logical :: direct_contribution,estimate_each_scatter,elastic_only
    Logical :: suppress_absorption,suppress_leakage,all_mat_mech
    Logical :: Gravity,neutron_decay,doppler_broaden
    Logical :: thermal_motion,diatomic_Atm,Rotating_Earth,wind
    Logical :: inbound_trajectories
    Integer :: setup_unit,stat
    Character(10) :: scatter_model  !IsoCM, AnIsoCM
    Character(10) :: dist_to_next_event  !An-exact,An-fast
    Logical :: roulette
    Real(dp) :: roulette_weight
    Integer :: roulette_ratio
    Real(dp) :: E_min,E_max
    
    NameList /NeutronScatterList/ n_scatters,scatter_model,elastic_only,suppress_absorption, &
                                  & suppress_leakage,all_mat_mech,dist_to_next_event, &
                                  & roulette,roulette_weight,roulette_ratio, &
                                  & direct_contribution,estimate_each_scatter,E_min,E_max, &
                                  & Gravity,neutron_decay,doppler_broaden,thermal_motion, &
                                  & Diatomic_Atm,Rotating_Earth,Wind,inbound_trajectories
    
    Open(NEWUNIT = setup_unit , FILE = setup_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  File open error, '//setup_file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Read(setup_unit,NML = NeutronScatterList)
    Close(setup_unit)
    If (n_scatters.EQ.-1 .AND. .NOT.estimate_each_scatter) Then
        Call Output_Message('ERROR:  Neutron_Scatter: Setup_Scatter_Model:  & 
                           &estimate_each_scatter must be .TRUE. for n_scatters = -1 (unlimited)',kill=.TRUE.)
    End If
    If (n_scatters.EQ.0 .AND. .NOT.direct_contribution) Then
        Call Output_Message('ERROR:  Neutron_Scatter: Setup_Scatter_Model:  direct_contribution must be .TRUE. for &
                           &n_scatters = 0',kill=.TRUE.)
    End If
    If (direct_contribution) Then
        If (n_scatters.GT.0 .AND. .NOT.estimate_each_scatter) Then
            Call Output_Message('ERROR:  Neutron_Scatter: Setup_Scatter_Model:  direct_contribution cannot be TRUE w/ &
                               &either n_scatters.GT.0 or estimate_each_scatter=FALSE ',kill=.TRUE.)
        End If
    End If
    Select Case (scatter_model)
        Case('IsoCM')
            ScatMod%aniso_dist = .FALSE.
        Case('AnIsoCM')
            ScatMod%aniso_dist = .TRUE.
        Case Default
            Call Output_Message('ERROR:  Neutron_Scatter: Setup_Scatter_Model:  Undefined scatter model',kill=.TRUE.)
    End Select
    ScatMod%n_scatters = n_scatters
    ScatMod%direct_contribution = direct_contribution
    If (atm_model_i .EQ. -1) Then !scattering should be disabled by no atmosphere model selected
        If (n_scatters.NE.0 .OR. .NOT.direct_contribution) Then
            Call Output_Message('ERROR:  Neutron_Scatter: Setup_Scatter_Model:  With no atmosphere model n_scatters must be 0 &
                               &and direct_contribution must be TRUE ',kill=.TRUE.)
        End If
    End If
    ScatMod%estimate_each_scatter = estimate_each_scatter
    ScatMod%elastic_only = elastic_only
    ScatMod%E_min = E_min
    ScatMod%E_max = E_max
    ScatMod%suppress_absorption = suppress_absorption
    ScatMod%suppress_leakage = suppress_leakage
    ScatMod%all_mat_mech = all_mat_mech
    Select Case (dist_to_next_event)
        Case('An-exact')
            ScatMod%fast_analog = .FALSE.
        Case('An-fast')
            !UNDONE Fast analog option not yet implemented
            !ScatMod%fast_analog = .TRUE.
            Call Output_Message('ERROR:  Neutron_Scatter: Setup_Scatter_Model:  Fast Analog Monte-Carlo pathlength game not &
                               &yet implemented',kill=.TRUE.)
            ERROR STOP
        Case Default
            Call Output_Message('ERROR:  Neutron_Scatter: Setup_Scatter_Model:  Undefined Monte-Carlo pathlength game',kill=.TRUE.)
    End Select
    ScatMod%roulette = roulette
    If (roulette) Then
        ScatMod%roulette_weight = roulette_weight * n_kill_weight
        ScatMod%roulette_rate = 1._dp / Real(roulette_ratio,dp)
        ScatMod%roulette_mult = Real(roulette_ratio,dp)
    Else
        ScatMod%roulette_weight = 0._dp
        ScatMod%roulette_rate = 1._dp
        ScatMod%roulette_mult = 1._dp
    End If
    ScatMod%Gravity = Gravity
    If (ScatMod%Gravity) Then
        ScatMod%inbound_trajectories = inbound_trajectories
    Else
        ScatMod%inbound_trajectories = .FALSE.
    End If
    ScatMod%Neutron_Decay = Neutron_Decay
    ScatMod%Doppler_Broaden = Doppler_Broaden
    If (Any( (/ Thermal_Motion , Rotating_Earth , Wind /) )) Then
        ScatMod%Target_Motion = .TRUE.
    Else
        ScatMod%Target_Motion = .FALSE.
    End If
    ScatMod%Thermal_Motion = Thermal_Motion
    ScatMod%Rotating_Earth = Rotating_Earth
    ScatMod%Wind = Wind
    ScatMod%n_kills = 0_id
    ScatMod%next_events = 0_id
    ScatMod%n_no_tally = 0_id
    ScatMod%n_uncounted = 0_id
    If (Worker_Index() .EQ. 1) Then
        Open(NEWUNIT = setup_unit , FILE = run_file_name , STATUS = 'OLD' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
        If (stat .NE. 0) Call Output_Message( 'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  File open error, '//run_file_name// & 
                                            & ', IOSTAT=',stat,kill=.TRUE.)
        Write(setup_unit,NML = NeutronScatterList)
        Write(setup_unit,*)
        Close(setup_unit)
    End If
End Function Setup_Scatter_Model

Subroutine Finish_Scatter_Model_Setup(ScatMod,CS_loaded,n_a_max,n_a_tab_max,n_iso)
    Implicit None
    Type(Scatter_Model_Type), Intent(InOut) :: ScatMod
    Logical, Intent(In) :: CS_loaded
    Integer, Intent(In), Optional :: n_a_max
    Integer, Intent(In), Optional :: n_a_tab_max
    Integer, Intent(In), Optional :: n_iso
    
    If (CS_loaded) Then !cross sections were loaded
        ScatMod%cs_loaded = .TRUE.
        !initialize scatter CS arrays for sampled scatters, except the lev_cs array (used for sampled scatters)
        Allocate(ScatMod%scat%a(0:n_a_max))
        ScatMod%scat%a = 0._dp
        Allocate(ScatMod%scat%a_tab1(1:2,1:n_a_tab_max))
        ScatMod%scat%a_tab1 = 0._dp
        Allocate(ScatMod%scat%a_tab2(1:2,1:n_a_tab_max))
        ScatMod%scat%a_tab2 = 0._dp
        Allocate(ScatMod%scat%iso_cs(1:n_iso))
        ScatMod%scat%iso_cs = 0._dp
    Else !cross sections are not needed if atmosphere is disabled
        ScatMod%cs_loaded = .FALSE.
    End If
End Subroutine Finish_Scatter_Model_Setup

Subroutine Sample_Scatter(ScatMod,n,CS,atm,RNG)
    Use Kinds, Only: dp
    Use Global, Only: mfp_per_barn_per_km_at_seaLevel
    Use Global, Only: Z_hat,X_hat
    Use Atmospheres, Only: Atmosphere_Type
    Use Random_Numbers, Only: RNG_Type
    Use n_Cross_Sections, Only: CS_type
    Use n_Cross_Sections, Only: sig_Composite_list
    Use n_Cross_Sections, Only: sig_Resonance
    Use Utilities, Only: Bisection_Search
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Use Interpolation, Only: Linear_Interp
    Use Target_Motion, Only: Random_Thermal_Velocity_AF
    Use Neutron_Utilities, Only: Neutron_Speed
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Class(Scatter_Model_Type), Intent(InOut) :: ScatMod
    Type(Neutron_Type), Intent(In):: n
    Type(CS_Type), Intent(In) :: CS
    Type(Atmosphere_Type), Intent(In) :: atm
    Type(RNG_Type), Intent(InOut) :: RNG
    Integer :: E_index,i,index1,index2
    Real(dp) :: r
    Real(dp) :: E_cm
    Real(dp) :: E1,E2
    Real(dp), Allocatable :: a1(:),a2(:),level_cs(:)
    Integer :: iso,lev
    Integer :: target_el
    Real(dp) :: An_prime,Mn
    Real(dp) :: E_apparent
    Real(dp) :: T
    Real(dp) :: v0(1:3),v0cm(1:3)
    Real(dp) :: resT,resS
    
    !Get atmosphere total/absorption and isotope scatter cross sections at the point of scatter
    If (ScatMod%Rotating_Earth .OR. ScatMod%Wind) Then
        ScatMod%scat%vAir = Air_Velocity(ScatMod%Rotating_Earth,ScatMod%Wind,n,atm)
        E_apparent = Apparent_Energy(n%E,n%Omega_hat,ScatMod%scat%vAir)
    Else
        ScatMod%scat%vAir = 0._dp
        E_apparent = n%E
    End If
    If (ScatMod%Doppler_Broaden) Then
        T = atm%T(n%Z)
        Call CS%sig_T_A(E_apparent,T,ScatMod%scat%sig_T,ScatMod%scat%sig_A)
        Do i = 1,CS%n_iso
            ScatMod%scat%iso_cs(i) = CS%sig_S(i,E_apparent,T)
        End Do
    Else
        Call CS%sig_T_A(E_apparent,ScatMod%scat%sig_T,ScatMod%scat%sig_A,iE_get=E_index)
        Do i = 1,CS%n_iso
            ScatMod%scat%iso_cs(i) = CS%sig_S(i,E_apparent,iE_put=E_index)
        End Do
    End If
    ScatMod%scat%sig_A = ScatMod%scat%sig_A * mfp_per_barn_per_km_at_seaLevel  !convert to macroscopic cross section at sea-level
    ScatMod%scat%sig_T = ScatMod%scat%sig_T * mfp_per_barn_per_km_at_seaLevel  !convert to macroscopic cross section at sea-level
    !Choose a target isotope based on relative SCATTER cross sections
    ScatMod%scat%iso_cs = ScatMod%scat%iso_cs * CS%iso_Fractions
    r = RNG%Get_Random() * Sum(ScatMod%scat%iso_cs)
    Do i = 1,CS%n_iso
        If (r .LT. Sum(ScatMod%scat%iso_cs(1:i))) Then
            ScatMod%scat%target_index = i
            ScatMod%scat%An = CS%An(i)
            Exit
        End If
    End Do
    !Add thermal motion
    If (ScatMod%Thermal_Motion) Then
        If (ScatMod%Diatomic_Atm) Then
            If (atm%diatomic(ScatMod%scat%target_index)) Then !Sample isotope of other nucleus in diatomic molecule
                r = RNG%Get_Random()
                target_el = atm%iso_map(ScatMod%scat%target_index)
                Do i = atm%iso_ind(target_el),atm%iso_ind(target_el+1)-1
                    If (r .LT. Sum(atm%iso_frac(atm%iso_ind(target_el):i))) Then
                        An_prime = CS%An(i)
                        Mn = ScatMod%scat%An + An_prime
                        Exit
                    End If
                End Do
                ScatMod%scat%vA = ScatMod%scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),ScatMod%scat%An,An_prime,Mn,RNG)
            Else !monatomic atmospheric constituent
                ScatMod%scat%vA = ScatMod%scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),ScatMod%scat%An,RNG)
            End If
        Else !monotomic approximation of atmosphere
            ScatMod%scat%vA = ScatMod%scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),ScatMod%scat%An,RNG)
        End If
    Else
        ScatMod%scat%vA = ScatMod%scat%vAir
    End If
    !convert energy to the center of mass of the collision frame, storing side effects for later use
    v0 = Neutron_Speed(n%E) * n%Omega_hat
    ScatMod%scat%u = (v0 + ScatMod%scat%An * ScatMod%scat%vA) / (ScatMod%scat%An + 1._dp)
    v0cm = v0 - ScatMod%scat%u
    E_cm = Neutron_Energy(v0cm)
    !find E_cm in the unified energy grid
    E_index = Bisection_Search(E_cm,CS%E_uni,CS%n_E_uni)
    !compute a few more predetermined quantities that are known before the scatter
    ScatMod%scat%u_speed = Vector_Length(ScatMod%scat%u)
    ScatMod%scat%Omega_hat0_cm = Unit_Vector(v0cm)
    ScatMod%scat%A_hat = ScatMod%scat%Omega_hat0_cm
    If (Abs(Dot_Product(ScatMod%scat%A_hat,Z_Hat)) .LT. 0.5_dp) Then
        ScatMod%scat%B_hat = Unit_Vector(Cross_Product(Z_Hat,ScatMod%scat%A_hat))
    Else
        ScatMod%scat%B_hat = Unit_Vector(Cross_Product(X_Hat,ScatMod%scat%A_hat))
    End If
    ScatMod%scat%C_hat = Cross_Product(ScatMod%scat%A_hat,ScatMod%scat%B_hat)
    !choose scatter level
    iso = ScatMod%scat%target_index
    If (ScatMod%elastic_only) Then
        ScatMod%scat%level = 0
        ScatMod%scat%Q = 0._dp
    Else  !sample scatter level (including elastic as level 0)
        If (E_index .LT. CS%lev_cs(iso)%thresh(1)) Then !not enough energy for lowest inelastic level, scatter is elastic
            ScatMod%scat%level = 0
            ScatMod%scat%Q = 0._dp
        Else !sample for scatter level
            Allocate(level_cs(0:CS%lev_cs(iso)%n_lev))
            level_cs = sig_Composite_list( E_cm,                  & 
                                         & CS%n_E_uni,            & 
                                         & CS%E_uni,              & 
                                         & CS%lnE_uni,            & 
                                         & E_index,               & 
                                         & 0,                     & 
                                         & CS%lev_cs(iso)%n_lev,  & 
                                         & CS%lev_cs(iso)%thresh, & 
                                         & CS%lev_cs(iso)%sig     )
            If (CS%has_res_cs(iso)) Then  !resonance contribution needs to be added to level 0 (elastic)
                Call sig_Resonance(CS%res_cs(iso),E_cm,resT,resS)
                level_cs(0) = level_cs(0) + resS
            End If
            r = RNG%Get_Random() * Sum(level_cs)
            Do i = 0,CS%lev_cs(iso)%n_Lev
                If (r .LT. Sum(level_cs(0:i))) Then
                    ScatMod%scat%level = i
                    ScatMod%scat%Q = CS%lev_cs(ScatMod%scat%target_index)%Q(i)
                    Exit
                End If
            End Do
        End If
    End If
    lev = ScatMod%scat%level
    If (ScatMod%aniso_dist) Then
        If (CS%lev_cs(iso)%da(lev)%is_iso) Then
            ScatMod%scat%da_is_iso = .TRUE.
            ScatMod%scat%da_is_legendre = .FALSE.
            ScatMod%scat%da_is_tab = .FALSE.
        Else
            index1 = CS%lev_cs(iso)%da(lev)%E_map(E_index) - 1
            index2 = CS%lev_cs(iso)%da(lev)%E_map(E_index)
            E1 = CS%E_uni( CS%lev_cs(iso)%da(lev)%E_key( index1 ) )
            E2 = CS%E_uni( CS%lev_cs(iso)%da(lev)%E_key( index2 ) )
            ScatMod%scat%da_is_iso = .FALSE.
            If (CS%lev_cs(iso)%da(lev)%da(index1)%is_Legendre) Then  !angular dist is expressed in legendre coeffs
                ScatMod%scat%da_is_legendre = .TRUE.
                ScatMod%scat%da_is_tab = .FALSE.
                ScatMod%scat%n_a = MaxVal(CS%lev_cs(iso)%da(lev)%da(index1:index2)%n_a)
                Allocate(a1(0:ScatMod%scat%n_a))
                Allocate(a2(0:ScatMod%scat%n_a))
                a1 = 0._dp
                a2 = 0._dp
                a1(0:CS%lev_cs(iso)%da(lev)%da(index1)%n_a) = CS%lev_cs(iso)%da(lev)%da(index1)%a
                a2(0:CS%lev_cs(iso)%da(lev)%da(index2)%n_a) = CS%lev_cs(iso)%da(lev)%da(index2)%a
                ScatMod%scat%a(0:ScatMod%scat%n_a) = Linear_Interp(E_cm,E1,E2,a1,a2)
            Else !angular dist is expressed in tabulated cosine pdf values
                ScatMod%scat%da_is_legendre = .FALSE.
                ScatMod%scat%da_is_tab = .TRUE.
                ScatMod%scat%n_a1 = CS%lev_cs(iso)%da(lev)%da(index1)%n_a
                ScatMod%scat%n_a2 = CS%lev_cs(iso)%da(lev)%da(index2)%n_a
                ScatMod%scat%a_tab1(:,1:ScatMod%scat%n_a1) = CS%lev_cs(iso)%da(lev)%da(index1)%ua
                ScatMod%scat%a_tab2(:,1:ScatMod%scat%n_a2) = CS%lev_cs(iso)%da(lev)%da(index2)%ua
                ScatMod%scat%a_tab_Econv = (E_cm - E1) / (E2 - E1)
            End If
        End If
    End If
    !speed of neutron after scatter is independent of direction in CM frame
    ScatMod%scat%s1cm = Neutron_Speed( ( E_cm + ScatMod%scat%An*Neutron_Energy(ScatMod%scat%vA-ScatMod%scat%u) - ScatMod%scat%Q) * &
                                       & ScatMod%scat%An / (ScatMod%scat%An+1._dp) )
End Subroutine Sample_Scatter

Subroutine Set_Scatter_prep(ScatMod,CS,scat)
    Use Kinds, Only: dp
    Use n_Cross_Sections, Only: CS_Type
    Implicit None
    Class(Scatter_Model_Type), Intent(In) :: ScatMod
    Type(CS_Type), Intent(In) :: CS
    Type(Scatter_Data_Type), Intent(Out) :: scat
    
    !Get total and absorption apparent cross sections from sampled scatter
    scat%sig_A = ScatMod%scat%sig_A
    scat%sig_T = ScatMod%scat%sig_T
    !Get air velocity from sampled scatter
    scat%vAir = ScatMod%scat%vAir
    !Set up arrays for angular distributions
    Allocate(scat%a(0:CS%n_a_max))
    scat%a = 0._dp
    Allocate(scat%a_tab1(1:2,1:CS%n_a_tab_max))
    Allocate(scat%a_tab2(1:2,1:CS%n_a_tab_max))
    scat%a_tab1 = 0._dp
    scat%a_tab2 = 0._dp
    !Get precomputed isotope cross sections from sampled scatter
    Allocate(scat%iso_cs(1:CS%n_iso))
    scat%iso_cs = ScatMod%scat%iso_cs / Sum(ScatMod%scat%iso_cs)
End Subroutine Set_Scatter_prep

Subroutine Set_Scatter_iso(ScatMod,n,CS,atm,RNG,scat,iso,n_lev,E_cm,i_E_cm)
    Use Kinds, Only: dp
    Use Global, Only: X_hat,Z_hat
    Use Atmospheres, Only: Atmosphere_Type
    Use Random_Numbers, Only: RNG_Type
    Use n_Cross_Sections, Only: CS_type
    Use n_Cross_Sections, Only: sig_Composite_list
    Use n_Cross_Sections, Only: sig_Resonance
    Use Utilities, Only: Bisection_Search
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Use Target_Motion, Only: Random_Thermal_Velocity_AF
    Use Neutron_Utilities, Only: Neutron_Speed
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Class(Scatter_Model_Type), Intent(In) :: ScatMod
    Type(Neutron_Type), Intent(In):: n
    Type(CS_Type), Intent(In):: CS
    Type(Atmosphere_Type), Intent(In) :: atm
    Type(RNG_Type), Intent(InOut) :: RNG
    Type(Scatter_Data_Type), Intent(InOut) :: scat
    Integer, Intent(In) :: iso
    Integer, Intent(Out) :: n_lev
    Real(dp), Intent(Out) :: E_cm
    Integer, Intent(Out) :: i_E_cm
    Integer :: i
    Real(dp) :: r
    Integer :: target_el
    Real(dp) :: An_prime,Mn
    Real(dp) :: v0(1:3),v0cm(1:3)
    Real(dp) :: resT,resS
    
    !Set target isotope
    scat%target_index = iso
    scat%An = CS%An(iso)
    !Add thermal motion
    If (ScatMod%Thermal_Motion) Then
        If (ScatMod%scat%target_index .EQ. iso) Then !the sampled scatter already has values we can use here
            scat%vA = ScatMod%scat%vA
        Else
            If (ScatMod%Diatomic_Atm) Then
                If (atm%diatomic(iso)) Then !Sample isotope of other nucleus in diatomic molecule
                    r = RNG%Get_Random()
                    target_el = atm%iso_map(iso)
                    Do i = atm%iso_ind(target_el),atm%iso_ind(target_el+1)-1
                        If (r .LT. Sum(atm%iso_frac(atm%iso_ind(target_el):i))) Then
                            An_prime = CS%An(i)
                            Mn = scat%An + An_prime
                            Exit
                        End If
                    End Do
                    scat%vA = scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),scat%An,An_prime,Mn,RNG)
                Else !monatomic atmospheric constituent
                    scat%vA = scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),scat%An,RNG)
                End If
            Else !monotomic approximation of atmosphere
                scat%vA = scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),scat%An,RNG)
            End If
        End If
    Else
        scat%vA = scat%vAir
    End If
    !convert energy to the center of mass of the collision frame, storing side effects for later use
    v0 = Neutron_Speed(n%E) * n%Omega_hat
    scat%u = (v0 + scat%An * scat%vA) / (scat%An + 1._dp)
    v0cm = v0 - scat%u
    E_cm = Neutron_Energy(v0cm)
    !find E_cm in the unified energy grid
    i_E_cm = Bisection_Search(E_cm,CS%E_uni,CS%n_E_uni)
    !compute a few more predetermined quantities that are known before the scatter
    scat%u_speed = Vector_Length(scat%u)
    scat%Omega_hat0_cm = Unit_Vector(v0cm)
    scat%A_hat = scat%Omega_hat0_cm
    If (Abs(Dot_Product(scat%A_hat,Z_Hat)) .LT. 0.5_dp) Then
        scat%B_hat = Unit_Vector(Cross_Product(Z_Hat,scat%A_hat))
    Else
        scat%B_hat = Unit_Vector(Cross_Product(X_Hat,scat%A_hat))
    End If
    scat%C_hat = Cross_Product(scat%A_hat,scat%B_hat)
    !Get level cross sections
    Allocate(scat%lev_cs(0:CS%lev_cs(iso)%n_lev))
    scat%lev_cs = sig_Composite_list( E_cm,                  & 
                                    & CS%n_E_uni,            & 
                                    & CS%E_uni,              & 
                                    & CS%lnE_uni,            & 
                                    & i_E_cm,                & 
                                    & 0,                     & 
                                    & CS%lev_cs(iso)%n_lev,  & 
                                    & CS%lev_cs(iso)%thresh, & 
                                    & CS%lev_cs(iso)%sig     )
    If (CS%has_res_cs(iso)) Then  !resonance contribution needs to be added to level 0 (elastic)
        Call sig_Resonance(CS%res_cs(iso),E_cm,resT,resS)
        scat%lev_cs(0) = scat%lev_cs(0) + resS
    End If
    n_lev = Count(scat%lev_cs .GT. 0._dp) - 1
    If (n_lev .GT. 0) Then
        scat%lev_cs = scat%lev_cs / Sum(scat%lev_cs)
    Else
        scat%lev_cs = 1._dp
    End If
End Subroutine Set_Scatter_iso

Subroutine Set_Scatter_lev(ScatMod,scat,CS,lev,E_cm,i_E_cm)
    Use Kinds, Only: dp
    Use n_Cross_Sections, Only: CS_type
    Use Interpolation, Only: Linear_Interp
    Use Neutron_Utilities, Only: Neutron_Speed
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Class(Scatter_Model_Type), Intent(In) :: ScatMod
    Type(Scatter_Data_Type), Intent(InOut) :: scat
    Type(CS_type), Intent(In) :: CS
    Integer, Intent(In) :: lev
    Real(dp), Intent(In) :: E_cm
    Integer, Intent(In) :: i_E_cm
    Integer :: index1,index2
    Real(dp) :: E1,E2
    Real(dp), Allocatable :: a1(:),a2(:)
    
    !Set scatter level
    scat%level = lev
    scat%Q = CS%lev_cs(scat%target_index)%Q(lev)
    If (ScatMod%aniso_dist) Then  !get angular distribution coeffs
        If (CS%lev_cs(scat%target_index)%da(lev)%is_iso) Then
            scat%da_is_iso = .TRUE.
            scat%da_is_legendre = .FALSE.
            scat%da_is_tab = .FALSE.
        Else
            index1 = CS%lev_cs(scat%target_index)%da(lev)%E_map(i_E_cm) - 1
            index2 = CS%lev_cs(scat%target_index)%da(lev)%E_map(i_E_cm)
            E1 = CS%E_uni( CS%lev_cs(scat%target_index)%da(lev)%E_key( index1 ) )
            E2 = CS%E_uni( CS%lev_cs(scat%target_index)%da(lev)%E_key( index2 ) )
            scat%da_is_iso = .FALSE.
            If (CS%lev_cs(scat%target_index)%da(lev)%da(index1)%is_Legendre) Then  !angular dist is expressed in legendre coeffs
                scat%da_is_legendre = .TRUE.
                scat%da_is_tab = .FALSE.
                scat%n_a = MaxVal(CS%lev_cs(scat%target_index)%da(lev)%da(index1:index2)%n_a)
                Allocate(a1(0:scat%n_a))
                Allocate(a2(0:scat%n_a))
                a1 = 0._dp
                a2 = 0._dp
                a1(0:CS%lev_cs(scat%target_index)%da(lev)%da(index1)%n_a) = CS%lev_cs(scat%target_index)%da(lev)%da(index1)%a
                a2(0:CS%lev_cs(scat%target_index)%da(lev)%da(index2)%n_a) = CS%lev_cs(scat%target_index)%da(lev)%da(index2)%a
                scat%a(0:scat%n_a) = Linear_Interp(E_cm,E1,E2,a1,a2)
            Else !angular dist is expressed in tabulated cosine pdf values
                scat%da_is_legendre = .FALSE.
                scat%da_is_tab = .TRUE.
                scat%n_a1 = CS%lev_cs(scat%target_index)%da(lev)%da(index1)%n_a
                scat%n_a2 = CS%lev_cs(scat%target_index)%da(lev)%da(index2)%n_a
                scat%a_tab1(:,1:scat%n_a1) = CS%lev_cs(scat%target_index)%da(lev)%da(index1)%ua
                scat%a_tab2(:,1:scat%n_a2) = CS%lev_cs(scat%target_index)%da(lev)%da(index2)%ua
                scat%a_tab_Econv = (E_cm - E1) / (E2 - E1)
            End If
        End If
    End If
    !speed of neutron after scatter is independent of direction in CM frame
    scat%s1cm = Neutron_Speed( ( E_cm + scat%An * Neutron_Energy(scat%vA - scat%u) - scat%Q) * scat%An / (scat%An + 1._dp) )
End Subroutine Set_Scatter_lev

Function Air_Velocity(Rotating_Earth,Wind,n,atm) Result(vA)
    Use Kinds, Only: dp
    Use Global, Only: Z_hat
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Use Atmospheres, Only: Atmosphere_Type
    Use Target_Motion, Only: Atm_Rotation_Velocity_AF
    Implicit None
    Real(dp) :: vA(1:3)
    Logical, Intent(In) :: Rotating_Earth,Wind
    Type(Neutron_Type), Intent(In) :: n
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp) :: vA_AF(1:3)
    Real(dp) :: E_hat(1:3),N_hat(1:3),U_hat(1:3)
    
    vA_AF = 0._dp
    !Add rotation of the atmosphere
    If (Rotating_Earth) vA_AF = vA_AF + Atm_Rotation_Velocity_AF(n%big_r,n%r(3)/n%big_r)
    !Add wind
    If (Wind) vA_AF = vA_AF + atm%wind_AF
    !Convert from Air Frame to interial frame
    U_hat = Unit_Vector(n%r)
    E_hat = Unit_Vector(Cross_Product(Z_hat,U_hat))
    N_hat = Cross_Product(U_hat,E_hat)
    vA = E_hat * vA_AF(1) + &
       & N_hat * vA_AF(2) + &
       & U_hat * vA_AF(3)
End Function Air_Velocity

Function Apparent_Energy(E,Omega_hat,vA) Result(E_app)
    Use Kinds, Only: dp
    Use Neutron_Utilities, Only: Neutron_Energy
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp) :: E_app
    Real(dp), Intent(In) :: E,Omega_hat(1:3),vA(1:3)
    
    E_app = Neutron_Energy( Neutron_Speed(E) * Omega_hat - vA )
End Function Apparent_Energy

Function Scattered_Direction(mu,omega,A_hat,B_hat,C_hat) Result(Omega_hat)
    Use Kinds, Only: dp
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Implicit None
    Real(dp) :: Omega_hat(1:3)  ! Direction after scatter
    Real(dp), Intent(In) :: mu  ! Deflection angle cosine
    Real(dp), Intent(In) :: omega  ! Orientation angle around Omega_hat0
    Real(dp), Intent(In) :: A_hat(1:3),B_hat(1:3),C_hat(1:3)
    
    Omega_hat = mu * A_hat + Sqrt(1._dp - mu**2) * (Cos(omega) * B_hat + Sin(omega) * C_hat)
End Function Scattered_Direction

Subroutine Scattered_Angles(Omega_hat0,Omega_hat1,mu,omega,B_hat,C_hat)
    Use Kinds, Only: dp
    Implicit None
    Real(dp), Intent(In) :: Omega_hat0(1:3) ! Direction before scatter
    Real(dp), Intent(In) :: Omega_hat1(1:3) ! Direction after scatter
    Real(dp), Intent(Out) :: mu  ! Deflection angle cosine
    Real(dp), Intent(Out) :: omega  ! Orientation angle around Omega_hat0
    Real(dp), Intent(In) :: B_hat(1:3),C_hat(1:3) 
    
    mu = Dot_Product(Omega_hat0,Omega_hat1)
    omega = Atan2( Dot_Product(Omega_hat1,B_hat) , Dot_Product(Omega_hat1,C_hat) )
End Subroutine Scattered_Angles

Subroutine Write_Scatter_Model(s,file_name)
    Use n_Cross_Sections, Only: Write_Cross_Sections
    Use FileIO_Utilities, Only: Output_Message
    Use FileIO_Utilities, Only: half_dash_line
    Implicit None
    Type(Scatter_Model_Type), Intent(In) :: s
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat
    
    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Neutron_Scatter: Write_Scatter_Model:  File open error, '//file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') 'SCATTER MODEL INFORMATION'
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') '  Scatter Model Fidelity:'
    If (s%aniso_dist) Then
        Write(unit,'(A)') '    Scatter Model:       Anisotropic Center-of-Mass'
    Else
        Write(unit,'(A)') '    Scatter Model:       Isotropic Center-of-Mass'
    End If
    If (s%n_scatters .LT. 0) Then
        Write(unit,'(A)') '    Number of scatters:  UNLIMITED'
    Else
        Write(unit,'(A,I4)') '    Number of scatters: ',s%n_scatters
    End If
    If (s%estimate_each_scatter) Then
        Write(unit,'(A)') '                         Contributions estimated at each scatter'
    Else
        Write(unit,'(A)') '                         Contributions estimated at FINAL scatter only'
    End If
    If (s%direct_contribution) Then
        Write(unit,'(A)') '                         Direct Contribution (first event) SAMPLED'
    Else
        Write(unit,'(A)') '                         Direct Contribution (first event) NOT included'
    End If
    If (s%elastic_only) Then
        Write(unit,'(A)') '    Collision Types:     Elastic only'
    Else
        Write(unit,'(A)') '    Collision Types:     Elastic & Inelastic'
    End If
    If (s%fast_analog) Then
        Write(unit,'(A)') '    Collision Distance Estimator:  Fast Approximation to Analog'
    Else
        Write(unit,'(A)') '    Collision Distance Estimator:  Analog'
    End If
    Write(unit,'(A,L1)') '    Absorption suppression: ',s%suppress_absorption
    Write(unit,'(A,L1)') '    Leakage suppression:    ',s%suppress_leakage
    Write(unit,'(A,L1)') '    All Materials & Mech:   ',s%all_mat_mech
    Write(unit,'(A,L1)') '    Gravity:                ',s%Gravity
    Write(unit,'(A,L1)') '    Inbound Trajectories:   ',s%inbound_trajectories
    Write(unit,'(A,L1)') '    Neutron Decay:          ',s%Neutron_Decay
    Write(unit,'(A,L1)') '    OTF Doppler Broaden:    ',s%Doppler_Broaden
    Write(unit,'(A,L1)') '    Target Motion:          ',s%Target_Motion
    Write(unit,'(A,L1)') '      Thermal Motion:         ',s%Thermal_Motion
    Write(unit,'(A,L1)') '      Diatomic Atm:           ',s%Diatomic_atm
    Write(unit,'(A,L1)') '      Rotating Earth:         ',s%Rotating_Earth
    Write(unit,'(A,L1)') '      Global Wind:            ',s%wind
    Write(unit,*)
    Write(unit,'(A,I15)') '  Histories killed on:  Time:     ',s%n_kills(1)
    Write(unit,'(A,I15)') '                        Energy:   ',s%n_kills(2)
    Write(unit,'(A,I15)') '                        Weight:   ',s%n_kills(3)
    Write(unit,'(A,I15)') '                        Leaked:   ',s%n_kills(4)
    Write(unit,'(A,I15)') '                        Absorbed: ',s%n_kills(5)
    Write(unit,'(A,I15)') '                        Coll Lim: ',s%n_kills(6)
    Write(unit,'(A,I15)') '  Next-Events to detector:  Attempted: ',s%next_events(1)
    Write(unit,'(A,I15)') '                            Found:     ',s%next_events(2)
    Write(unit,'(A,I15)') '                            Tallied:   ',s%next_events(3)
    Write(unit,'(A,I15)') '  Next-Events NOT tallied:  Time-Energy: ',s%n_no_tally(1)
    Write(unit,'(A,I15)') '                            Time:        ',s%n_no_tally(2)
    Write(unit,'(A,I15)') '                            Energy:      ',s%n_no_tally(3)
    Write(unit,'(A,I15)') '  Histories implicitly leaked due to EXO source:  ',s%n_uncounted
    Write(unit,*)
    Write(unit,*)
    Close(unit)
End Subroutine Write_Scatter_Model

End Module Neutron_Scatter
