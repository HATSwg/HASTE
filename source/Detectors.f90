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
Module Detectors

    Use Kinds, Only: dp
    Use Tallies, Only: Contrib_triplet
    Use Satellite_Motion, Only: Satellite_Position_Type
    Implicit None
    Private
    Public :: Detector_Type
    Public :: Grid_Info_Type
    Public :: Define_Grid_Info
    Public :: Setup_Detector
    Public :: Write_Detector
    Public :: Close_Slice_files

    Type :: Grid_Info_Type
        Real(dp) :: min
        Real(dp) :: max
        Real(dp) :: res
        Real(dp), Allocatable :: bounds(:)
        Integer :: log_min,log_max
        Integer :: n_decades
        Integer :: n_bins
        Logical :: log_spacing
        !variables for shape data collection
        Logical, Allocatable :: collect_shape(:)  !flag indication if collection is still going for this slice
        Integer, Allocatable :: slice_bin(:)  !bin in which to collect data for the shape of this slice
        Integer, Allocatable :: slice_c(:)  !number of contributions recorded in this slice
        Integer, Allocatable :: slice_unit(:)  !file unit of connected file to record contributions in the slice
    Contains
        Procedure, Pass :: Bin_number  !returns the index of the bin in which a value falls
        Procedure, Pass :: Bin_center  !returns the geometric 'center' of a bin based on bin spacing type
    End Type

    Type :: Detector_Type
        Type(Satellite_Position_Type) :: sat  !defines motion of the detector
        Logical :: exoatmospheric  !T indicates detector is outside atmosphere, F indicates detector IN the atmosphere
        Type(Grid_Info_Type) :: TE_grid(1:2)
        Type(Grid_Info_Type) :: Dir_grid(1:2)
        Integer :: TE_contrib_index  !number of time-energy contributions stored for the current history
        Integer :: Dir_contrib_index  !number of direction contributions stored for the current history
        Integer :: Contrib_size  !number of contributions spaces for the current history
        Type(Contrib_triplet), Allocatable :: TE_contribs_this_history(:)  !list of contribs to time-energy for the current hist
        Type(Contrib_triplet), Allocatable :: Dir_contribs_this_history(:)  !list of contribs to dir of arrival for the current hist
        Real(dp), Allocatable :: TE_contribs_t(:)  !list of contributions to time bins for the current history
        Real(dp), Allocatable :: TE_contribs_E(:)  !list of contributions to energy bins for the current history
        Real(dp), Allocatable :: Dir_contribs_mu(:)  !list of contributions to mu bins for the current history
        Real(dp), Allocatable :: Dir_contribs_omega(:)  !list of contributions to omega bins for the current history
        Logical :: shape_data
        Integer :: max_shape_data
        Integer :: n_slices
    Contains
        Procedure, Pass :: Tally_Scatter  !records contribution in TE contribs arrays
        Procedure, Pass :: Combine_Duplicate_Tallies  !scrubs list of contributions: combines duplicate bins, condenses the list
    End Type

Contains

Function Setup_Detector(setup_file_name,resources_dir,run_file_name,slice_file_name,R_top_atm) Result(d)
    Use Kinds, Only: dp
    Use Global, Only: Pi
    Use Global, Only: TwoPi
    Use Global, Only: deg2r
    Use Global, Only: Rc => R_center
    Use Global, Only: Z_hat
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Cross_Product
    Use Utilities, Only: Celest_to_XYZ
    Use Satellite_Motion, Only: Initialize_Satellite_Motion
    Use FileIO_Utilities, Only: Output_Message
    Use FileIO_Utilities, Only: Worker_Index
    Implicit None
    Type(Detector_Type) :: d
    Character(*), Intent(In) :: setup_file_name
    Character(*), Intent(In) :: resources_dir
    Character(*), Intent(In) :: run_file_name
    Character(*), Intent(In) :: slice_file_name
    Real(dp), Intent(In) :: R_top_atm
    Real(dp) :: x_detector,y_detector,z_detector  ![km]  x,y,z coordinates of detector
    Character(3) :: velocity_frame
    Real(dp) :: declination_detector,hour_angle_detector
    Real(dp) :: v_A_detector,v_B_detector,v_C_detector
    Real(dp) :: E_max,E_min  !specifies min and max energies for detector grid
    Real(dp) :: t_max,t_min  !specifies min and max times for detector grid
    Real(dp) :: E_res,t_res  !resolution of time and energy grids
    Integer :: E_bins_per_decade  !number of energy bins per decade, used for 'log' grid spacing 
    Integer :: t_bins_per_decade  !number of time bins per decade, used for 'log' grid spacing
    Logical :: log_sp
    Character(10) :: E_grid_spacing,t_grid_spacing  !specifies 'Log' or 'Linear' for grid spacing
    Integer :: n_mu_bins,n_omega_bins
    Integer :: setup_unit,stat,i
    Character(9) :: position_geometry
    Character(10) :: detector_motion
    Real(dp) :: E_hat(1:3),N_hat(1:3),U_hat(1:3)
    Logical :: collect_shape_data
    Integer :: shape_data_n_slices,shape_data_limit
    Integer :: j
    Character(3) :: j_char
    Character(9) :: slice_name_end
    NameList /NeutronDetectorList/ position_geometry,x_detector,y_detector,z_detector, &
                                   & declination_detector,hour_angle_detector, &
                                   & velocity_frame,v_A_detector,v_B_detector,v_C_detector, &
                                   & detector_motion, &
                                   & E_max,E_min,E_grid_spacing,E_res,E_bins_per_decade, &
                                   & t_max,t_min,t_grid_spacing,t_res,t_bins_per_decade, &
                                   & n_mu_bins,n_omega_bins,collect_shape_data,shape_data_n_slices, &
                                   & shape_data_limit

    Open(NEWUNIT = setup_unit , FILE = setup_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Detectors: Setup_Detector:  File open error, '//setup_file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Read(setup_unit,NML = NeutronDetectorList)
    Close(setup_unit)
    Select Case(position_geometry)
        Case('Celestial')
            d%sat%r0 = Celest_to_XYZ(Rc+z_detector,hour_angle_detector*deg2r,declination_detector*deg2r)
        Case('Cartesian')
            d%sat%r0 = (/ x_detector , y_detector , z_detector /)
        Case Default
            Call Output_Message('ERROR:  Detectors: Setup_Detector:  Unknown position geometry.',kill=.TRUE.)
    End Select
    !check for exoatmospheric detector
    If (Vector_Length(d%sat%r0) .LT. R_top_atm) Then
        d%exoatmospheric = .FALSE.
        Call Output_Message('ERROR:  Detectors: Setup_Detector:  Endo-atmospheric detector not supported.',kill=.TRUE.)
    Else
        d%exoatmospheric = .TRUE.
    End If
    If (Any( (/v_A_detector,v_B_detector,v_C_detector/)  .NE. 0._dp)) Then
        Select Case(velocity_frame)
            Case('XYZ')
                d%sat%v0 = (/ v_A_detector , v_B_detector , v_C_detector /)
            Case('ENU')
                U_hat = Unit_Vector(d%sat%r0)
                E_hat = Unit_Vector(Cross_Product(Z_hat,U_hat))
                N_hat = Cross_Product(U_hat,E_hat)
                d%sat%v0 = E_hat * v_A_detector + &
                         & N_hat * v_B_detector + &
                         & U_hat * v_C_detector
            Case Default
                Call Output_Message('ERROR:  Detectors: Setup_Detector:  Unknown velocity frame.',kill=.TRUE.)
        End Select
    Else
        d%sat%v0 = 0._dp
    End If
    Call Initialize_Satellite_Motion(resources_dir,detector_motion,d%sat)
    !create time grid
    Select Case(t_grid_spacing)
        Case('Log')
            log_sp = .TRUE.
        Case('Linear')
            log_sp = .FALSE.
        Case Default
            Call Output_Message('ERROR:  Detectors: Setup_Detector:  Undefined t grid spacing',kill=.TRUE.)
    End Select
    d%TE_grid(1) = Define_Grid_Info(log_sp,t_min,t_max,t_res,t_bins_per_decade)
    !Create energy grid
    Select Case(E_grid_spacing)
        Case('Log')
            log_sp = .TRUE.
        Case('Linear')
            log_sp = .FALSE.
        Case Default
            Call Output_Message('ERROR:  Detectors: Setup_Detector:  Undefined E grid spacing',kill=.TRUE.)
    End Select
    d%TE_grid(2) = Define_Grid_Info(log_sp,E_min,E_max,E_res,E_bins_per_decade)
    !create mu grid
    d%Dir_grid(1) = Define_Grid_Info(.FALSE.,-1._dp,1._dp,2._dp / Real(n_mu_bins,dp),1)
    !create omega grid
    d%Dir_grid(2) = Define_Grid_Info(.FALSE.,-Pi,Pi,TwoPi / Real(n_omega_bins,dp),1)
    !Initialize contribution lists
    d%TE_contrib_index = 0
    d%Contrib_size = 2**8
    Allocate(d%TE_contribs_this_history(1:d%Contrib_size))
    d%TE_contribs_this_history(:)%i1 = -1
    d%TE_contribs_this_history(:)%i2 = -1
    d%TE_contribs_this_history(:)%f = 0._dp
    Allocate(d%TE_contribs_t(1:d%TE_grid(1)%n_bins))
    d%TE_contribs_t = 0._dp
    Allocate(d%TE_contribs_E(1:d%TE_grid(2)%n_bins))
    d%TE_contribs_E = 0._dp
    !Setup arrival direcion lists
    d%Dir_contrib_index = 0
    Allocate(d%Dir_contribs_this_history(1:d%Contrib_size))
    d%Dir_contribs_this_history(:)%i1 = -1
    d%Dir_contribs_this_history(:)%i2 = -1
    d%Dir_contribs_this_history(:)%f = 0._dp
    Allocate(d%Dir_contribs_mu(1:d%Dir_grid(1)%n_bins))
    d%Dir_contribs_mu = 0._dp
    Allocate(d%Dir_contribs_omega(1:d%Dir_grid(2)%n_bins))
    d%Dir_contribs_omega = 0._dp
    If (Worker_Index() .EQ. 1) Then
        !set up shape data collection
        If (collect_shape_data) Then
            d%shape_data = .TRUE.
            d%max_shape_data = shape_data_limit
            d%n_slices = shape_data_n_slices
            Do i = 1,2
                Allocate(d%TE_grid(i)%collect_shape(1:d%n_slices))
                d%TE_grid(i)%collect_shape = .TRUE.
                Allocate(d%TE_grid(i)%slice_bin(1:d%n_slices))
                Do CONCURRENT (j = 1:d%n_slices)
                    d%TE_grid(i)%slice_bin(j) = j * d%TE_grid(i)%n_bins / d%n_slices
                End Do
                Allocate(d%TE_grid(i)%slice_c(1:d%n_slices))
                d%TE_grid(i)%slice_c = 0
                Allocate(d%TE_grid(i)%slice_unit(1:d%n_slices))
                d%TE_grid(i)%slice_unit = 0
                Do j = 1,d%n_slices
                    Write(j_char,'(I3.3)') j
                    If (i .EQ. 1) Then
                        slice_name_end = '-t'//j_char//'.txt'
                    Else
                        slice_name_end = '-e'//j_char//'.txt'
                    End If
                    Open( NEWUNIT = d%TE_grid(i)%slice_unit(j) , FILE = slice_file_name//slice_name_end , STATUS = 'REPLACE' , &
                        & ACTION = 'WRITE' , IOSTAT = stat)
                    Write(d%TE_grid(i)%slice_unit(j),'(ES27.16E3)') d%TE_grid(i)%Bin_Center(d%TE_grid(i)%slice_bin(j))
                End Do
            End Do
        Else
            d%shape_data = .FALSE.
        End If
        !Write setup info to file
        Open(NEWUNIT = setup_unit , FILE = run_file_name , STATUS = 'OLD' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
        If (stat .NE. 0) Call Output_Message( 'ERROR:  Detectors: Setup_Detector:  File open error, '//run_file_name// & 
                                            & ', IOSTAT=',stat,kill=.TRUE.)
        Write(setup_unit,NML = NeutronDetectorList)
        Write(setup_unit,*)
        Close(setup_unit)
    Else
        d%shape_data = .FALSE.
    End If
End Function Setup_Detector

Function Define_Grid_Info(log_sp,x_min,x_max,x_res,ipd) Result(g)
    Use Kinds, Only : dp
    Implicit None
    Type(Grid_Info_Type) :: g
    Logical, Intent(In) :: log_sp
    Real(dp), Intent(In) :: x_min
    Real(dp), Intent(In) :: x_max
    Real(dp), Intent(In) :: x_res
    Integer, Intent(In) :: ipd
    Integer :: i

    g%log_spacing = log_sp
    If (log_sp) Then !logarithmic spacing
        !adjust x_min, x_max to enclose a set of complete regular decades
        g%log_min = Floor(log10(x_min))
        g%min = 10._dp**g%log_min
        g%log_max = Ceiling(log10(x_max))
        g%max = 10._dp**g%log_max
        g%n_decades = g%log_max - g%log_min
        g%n_bins = g%n_decades * ipd
        Allocate(g%bounds(0:g%n_bins))
        Do CONCURRENT (i = 0:g%n_bins)
            g%bounds(i) = 10._dp**(Real(g%log_min,dp) + Real(i,dp) / Real(ipd,dp))
        End Do
    Else !linear spacing
        g%min = x_min
        g%max = x_max
        g%n_bins = Ceiling((x_max - x_min) / x_res)
        Allocate(g%bounds(0:g%n_bins))
        Do CONCURRENT (i = 0:g%n_bins)
            g%bounds(i) = x_min + Real(i,dp) * x_res
        End Do
        g%res = x_res
    End If
End Function Define_Grid_Info

Subroutine Tally_Scatter(d,E,Omega_Hat,t,weight)
    Use Kinds, Only: dp
    Use Global, Only: X_hat
    Use Global, Only: Y_hat
    Use Global, Only: Z_hat
    Use Global, Only: Pi
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Implicit None
    Class(detector_type), Intent(InOut) :: d
    Real(dp), Intent(In) :: E  ![keV] neutron energy at arrival at detector
    Real(dp), Intent(In) :: Omega_Hat(1:3)  ![1 km/s] neutron direction of arrival at detector
    Real(dp), Intent(In) :: t  ![s] neutron time of arrival at detector
    Real(dp), Intent(In) :: weight  !weight (probability of neutron existence in this state) of neutron at arrival at detector
    Integer :: E_bin,t_bin
    Real(dp) :: r_sat(1:3),v_sat(1:3)
    Real(dp) :: D_hat(1:3)  !DOWN
    Real(dp) :: N_hat(1:3)  !negative orbit-normal
    Real(dp) :: F_hat(1:3)  !FORWARD
    Real(dp) :: Arr_dir_DNF(1:3)
    Integer :: mu_bin,omega_bin
    Integer :: i

    If (weight .EQ. 0._dp) Return !don't bother tallying zeroes
    t_bin = d%TE_grid(1)%Bin_number(t)
    E_bin = d%TE_grid(2)%Bin_number(E)
    d%TE_contrib_index = d%TE_contrib_index + 1
    !record shape data if applicable
    If (d%shape_data) Then
        Do i = 1,d%n_slices
            If (t_bin .EQ. d%TE_grid(1)%slice_bin(i)) Then
                If (d%TE_grid(1)%collect_shape(i)) Then
                    Write(d%TE_grid(1)%slice_unit(i),'(2ES27.16E3)') E,weight
                    d%TE_grid(1)%slice_c(i) = d%TE_grid(1)%slice_c(i) + 1
                    If (d%TE_grid(1)%slice_c(i) .GE. d%max_shape_data) Then
                        Close(d%TE_grid(1)%slice_unit(i))
                        d%TE_grid(1)%collect_shape(i) = .FALSE.
                    End If
                End If
            End If
            If (E_bin .EQ. d%TE_grid(2)%slice_bin(i)) Then
                If (d%TE_grid(2)%collect_shape(i)) Then
                    Write(d%TE_grid(2)%slice_unit(i),'(2ES27.16E3)') t,weight
                    d%TE_grid(2)%slice_c(i) = d%TE_grid(2)%slice_c(i) + 1
                    If (d%TE_grid(2)%slice_c(i) .GE. d%max_shape_data) Then
                        Close(d%TE_grid(2)%slice_unit(i))
                        d%TE_grid(2)%collect_shape(i) = .FALSE.
                    End If
                End If
            End If
        End Do
        If ( .NOT.Any(d%TE_grid(1)%collect_shape) .AND. .NOT.Any(d%TE_grid(2)%collect_shape) ) d%shape_data = .FALSE.
    End If
    !resize the list if necessary
    If (d%TE_contrib_index .GT. d%Contrib_size) Then
        Call Resize_Contrib_Lists(d%TE_contribs_this_history,d%Dir_contribs_this_history,d%Contrib_size)
    End If
    d%TE_contribs_this_history(d%TE_contrib_index)%i1 = t_bin
    d%TE_contribs_this_history(d%TE_contrib_index)%i2 = E_bin
    d%TE_contribs_this_history(d%TE_contrib_index)%f = weight
    d%TE_contribs_t(t_bin) = d%TE_contribs_t(t_bin) + weight
    d%TE_contribs_E(E_bin) = d%TE_contribs_E(E_bin) + weight
    !convert incoming vector to DNF coordinates
    If (d%sat%is_stationary) Then
        !N2H For stationary detector, we can save arithmetic by computing and saving the orientation basis just once
        D_hat = -Unit_Vector(d%sat%r0)
        If (d%sat%vp .GT. 0._dp) Then !stationary detector has velocity to use as reference direction
            N_hat = Unit_Vector(Cross_Product(D_hat,d%sat%v0))
            F_hat = Cross_Product(N_hat,D_hat)
            !Otherwise, use "easterly" (counterclockwise as viewed from positive Z in ECI)
        Else If (Abs(Dot_Product(D_hat,Z_hat)) .LT. 1._dp) Then
            F_hat = Unit_Vector(Cross_Product(D_Hat,Z_hat))
            N_hat = Cross_Product(F_hat,D_hat)
        Else !special case for polar location
            N_hat = -X_hat
            F_hat = Y_hat
        End If
    Else !non-stationary satellite motions have a basis that can be referenced from position and direction of travel
        Call d%sat%R_and_V(t,r_sat,v_sat)
        D_hat = -Unit_Vector(r_sat)
        N_hat = Unit_Vector(Cross_Product(D_hat,v_sat))
        F_hat = Cross_Product(N_hat,D_hat)
    End If
    F_hat = Cross_Product(N_hat,D_hat)
    Arr_dir_DNF = -(/ Dot_Product(Omega_hat,D_hat), &
                    & Dot_Product(Omega_hat,N_hat), &
                    & Dot_Product(Omega_hat,F_hat) /)
    mu_bin = Ceiling((Arr_dir_DNF(1)+1._dp) / d%Dir_grid(1)%res)
    omega_bin = Floor((ATan2(Arr_dir_DNF(2),Arr_dir_DNF(3)) + Pi) / d%Dir_grid(2)%res) + 1
    d%Dir_contrib_index = d%Dir_contrib_index + 1
    d%Dir_contribs_this_history(d%Dir_contrib_index)%i1 = mu_bin
    d%Dir_contribs_this_history(d%Dir_contrib_index)%i2 = omega_bin
    d%Dir_contribs_this_history(d%Dir_contrib_index)%f = weight
    d%Dir_contribs_mu(mu_bin) = d%Dir_contribs_mu(mu_bin) + weight
    d%Dir_contribs_omega(omega_bin) = d%Dir_contribs_omega(omega_bin) + weight
End Subroutine Tally_Scatter

Function Bin_Number(g,x) Result(b)
    Use Kinds, Only: dp
    Implicit None
    Integer :: b
    Class(Grid_Info_Type), Intent(In) :: g
    Real(dp), Intent(In) :: x

    If (g%log_spacing) Then
        b = Ceiling(Real(g%n_bins,dp) * (log10(x) - Real(g%log_min,dp)) / Real(g%n_decades,dp))
    Else  !linear spacing
        b = Ceiling((x - g%min) / g%res)
    End If
End Function Bin_Number

Subroutine Resize_Contrib_Lists(list1,list2,size)
    Use Kinds, Only: dp
    Implicit None
    Type(Contrib_triplet), Intent(InOut), Allocatable :: list1(:)
    Type(Contrib_triplet), Intent(InOut), Allocatable :: list2(:)
    Integer, Intent(InOut) :: size
    Type(Contrib_triplet), Allocatable :: swap(:)

    !Prepare swap array
    Allocate(swap(1:size))
    !Resize first (time-energy) list
    swap = list1
    Deallocate(list1)
    Allocate(list1(1:2*size))
    list1(:)%i1 = -1
    list1(:)%i2 = -1
    list1(:)%f = 0._dp
    list1(1:size) = swap
    !resize second (arrival direction) list
    swap = list2
    Deallocate(list2)
    Allocate(list2(1:2*size))
    list2(:)%i1 = -1
    list2(:)%i2 = -1
    list2(:)%f = 0._dp
    list2(1:size) = swap
    !update new size
    size = 2 * size
    Deallocate(swap)
End Subroutine Resize_Contrib_Lists

Subroutine Combine_Duplicate_Tallies(d)
    Implicit None
    Class(Detector_Type), Intent(InOut) :: d

    Call Sort_Combine_Triplets(d%Contrib_size,d%TE_contribs_this_history,d%TE_contrib_index)
    Call Sort_Combine_Triplets(d%Contrib_size,d%Dir_contribs_this_history,d%Dir_contrib_index)
End Subroutine Combine_Duplicate_Tallies

Subroutine Sort_Combine_Triplets(size,list,index)
    Implicit None
    Integer, Intent(In) :: size
    Type(Contrib_Triplet), Intent(InOut) :: list(1:size)
    Integer, Intent(InOut) :: index
    Integer :: i,j,k,n_1,n_2

    !Sort the list of contributions for this history on 1st index
    Call Contribution_Cocktail_Sort_i1(index,list)
    !Within each 1st index bin: sort in 2nd index and then consolidate duplicate contributions
    i = 1
    Do While (i .LT. index)
        !check for duplicate bins
        If (list(i)%i1 .EQ. list(i+1)%i1) Then !at least two in this 1st index bin, sort and check for duplicate 2nd index bins
            !count total contributions in this time bin
            n_1 = 2
            Do j = i+1,index-1
                If (list(i)%i1 .EQ. list(j+1)%i1) Then  !still in same 1st index bin
                    n_1 = n_1 + 1
                Else  !moved on to next 1st index bin
                    Exit
                End If
            End Do
            !sort the list of contributions in this 1st index bin on 2nd index
            Call Contribution_Cocktail_Sort_i2(n_1,list(i:i+n_1-1))
            !Search for and sum duplicates
            j = i
            Do While (j .LT. i+n_1-1)
                If (list(j)%i2 .EQ. list(j+1)%i2) Then !at least two in this 2nd index bin
                    !count total contributions in this bin
                    n_2 = 1
                    Do k = j+1,i+n_1-1
                        If (list(j)%i2 .EQ. list(k)%i2) Then  !still in same 2nd index bin
                            !add the contribution to the first instance of this bin in the list
                            list(j)%f = list(j)%f + list(k)%f
                            !update the count of duplicates
                            n_2 = n_2 + 1
                        Else  !moved on to next 2nd index bin
                            Exit
                        End If
                    End Do
                    !discard the duplicates by shifting the rest of the list up
                    list(j+1:index-n_2+1) = list(j+n_2:index)
                    index = index - n_2 + 1
                    n_1 = n_1 - n_2 + 1
                    j = j + 1
                Else !no duplicate, advance j by 1 to next 2nd index bin
                    j = j + 1
                End If
            End Do
            i = i + n_1
        Else  !no duplicate, advance i by 1 to next 1st index bin
            i = i + 1
        End If
    End Do
End Subroutine Sort_Combine_Triplets

Subroutine Contribution_Cocktail_Sort_i1(n,c)
    Implicit None
    Integer, Intent(In) :: n  !length of contribution triplet list
    Type(Contrib_triplet), Intent(InOut) :: c(1:n)  !list of contribution triplets to sort
    Type(Contrib_triplet)  :: swap
    Logical :: no_swaps
    Integer :: i,j
    
    Do  j = 1,n  !Cocktail sort (bidirectional bubble sort)
        !Bubble sort down through the list
        no_swaps = .TRUE.
        Do i = j,n-j
            If (c(i)%i1 .GT. c(i+1)%i1) Then
                swap = c(i)
                c(i) = c(i+1)
                c(i+1) = swap
                no_swaps = .FALSE.
            End If
        End Do
        If (no_swaps) Exit  !if nothing moved, then the list is sorted
        !Bubble sort back up through the list
        no_swaps = .TRUE.
        Do i = n,j+1,-1
            If (c(i)%i1 .LT. c(i-1)%i1) Then
                swap = c(i)
                c(i) = c(i-1)
                c(i-1) = swap
                no_swaps = .FALSE.
            End If
        End Do
        If (no_swaps) Exit  !if nothing moved, then the list is sorted
    End Do
End Subroutine Contribution_Cocktail_Sort_i1

Subroutine Contribution_Cocktail_Sort_i2(n,c)
    Implicit None
    Integer, Intent(In) :: n  !length of contribution triplet list
    Type(Contrib_triplet), Intent(InOut) :: c(1:n)  !list of contribution triplets to sort
    Type(Contrib_triplet)  :: swap
    Logical :: no_swaps
    Integer :: i,j
    
    Do  j = 1,n  !Cocktail sort (bidirectional bubble sort)
        !Bubble sort down through the list
        no_swaps = .TRUE.
        Do i = j,n-j
            If (c(i)%i2 .GT. c(i+1)%i2) Then
                swap = c(i)
                c(i) = c(i+1)
                c(i+1) = swap
                no_swaps = .FALSE.
            End If
        End Do
        If (no_swaps) Exit  !if nothing moved, then the list is sorted
        !Bubble sort back up through the list
        no_swaps = .TRUE.
        Do i = n,j+1,-1
            If (c(i)%i2 .LT. c(i-1)%i2) Then
                swap = c(i)
                c(i) = c(i-1)
                c(i-1) = swap
                no_swaps = .FALSE.
            End If
        End Do
        If (no_swaps) Exit  !if nothing moved, then the list is sorted
    End Do
End Subroutine Contribution_Cocktail_Sort_i2

Subroutine Write_Detector(d,file_name)
    Use Kinds, Only: dp
    Use Global, Only: Pi
    Use Global, Only: Rc => R_center
    Use Utilities, Only: Vector_Length
    Use FileIO_Utilities, Only: Output_Message
    Use FileIO_Utilities, Only: half_dash_line
    Implicit None
    Type(Detector_Type), Intent(In) :: d
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat
    Integer :: i
    Logical :: write_full_position_trace
    Real(dp) :: t,r(1:3),v(1:3)

    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Detectors: Write_Detector:  File open error, '//file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') 'DETECTOR INFORMATION'
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') '  Position (t = 0) :'
    Write(unit,'(A,ES24.16E3,A)') '    x = ',d%sat%r0(1),' km'
    Write(unit,'(A,ES24.16E3,A)') '    y = ',d%sat%r0(2),' km'
    Write(unit,'(A,ES24.16E3,A)') '    z = ',d%sat%r0(3),' km'
    Write(unit,'(A,ES24.16E3,A)') '    Right Ascension = ',Acos( & 
                                                               & d%sat%r0(1)/(Vector_Length(d%sat%r0)* & 
                                                               & Cos(Asin(d%sat%r0(3)/Vector_Length(d%sat%r0))))) & 
                                                          & / (Pi/180._dp),' deg'
    Write(unit,'(A,ES24.16E3,A)') '    Declination     = ',Asin(d%sat%r0(3)/Vector_Length(d%sat%r0)) / (Pi/180._dp),' deg'
    Write(unit,'(A,ES24.16E3,A)') '    Altitude        = ',Vector_Length(d%sat%r0)-Rc,' km'
    If (d%sat%is_stationary) Then
        Write(unit,'(A)') '  Motion Type:  STATIONARY'
        write_full_position_trace = .FALSE.
    Else If (d%sat%is_conic) Then
        Write(unit,'(A)') '  Motion Type:  ORBITAL'
        write_full_position_trace = .TRUE.
    Else
        Write(unit,'(A)') '  Motion Type:  LINEAR'
        write_full_position_trace = .TRUE.
    End If
    Write(unit,'(A)') '  Velocity (t = 0) :'
    Write(unit,'(A,ES24.16E3,A)') '    x = ',d%sat%v0(1),' km/s'
    Write(unit,'(A,ES24.16E3,A)') '    y = ',d%sat%v0(2),' km/s'
    Write(unit,'(A,ES24.16E3,A)') '    z = ',d%sat%v0(3),' km/s'
    Write(unit,*)
    Write(unit,'(A)') '  Position Trace:'
    Write(unit,'(A7,7A27)') 'Bin #','t [sec]','x-pos [km]','y-pos [km]','z-pos [km]','x-vel [km/s]','y-vel [km/s]','z-vel [km/s]'
    Write(unit,'(A7,7A27)') '-----','-------------------------','-------------------------','-------------------------', & 
                          & '-------------------------','-------------------------','-------------------------', & 
                          & '-------------------------'
    Write(unit,'(I7,7ES27.16E3)') 0,0._dp,d%sat%r0(1),d%sat%r0(2),d%sat%r0(3),d%sat%v0(1),d%sat%v0(2),d%sat%v0(3)
    If (write_full_position_trace) Then
        Do i = 1,d%TE_grid(1)%n_bins
            t = d%TE_grid(1)%bounds(i)
            Call d%sat%R_and_V(t,r,v)
            Write(unit,'(I7,7ES27.16E3)') i,t,r(1),r(2),r(3),v(1),v(2),v(3)
        End Do
    End If
    Write(unit,*)
    Write(unit,'(A,ES24.16E3,A,ES24.16E3,A)') '  Time range: ',d%TE_grid(1)%min,' sec to ',d%TE_grid(1)%max,' sec'
    If (d%TE_grid(1)%log_spacing) Then
        Write(unit,'(A,I5,A,I5,A)') '  Logarithmically spaced bins: ', & 
                                  & d%TE_grid(1)%n_bins, & 
                                  & ' bins, with ', & 
                                  & d%TE_grid(1)%n_bins / d%TE_grid(1)%n_decades, & 
                                  & ' bins per decade'
    Else
        Write(unit,'(A,I7,A,ES24.16E3,A)') '  Linearly spaced bins: ', & 
                                         & d%TE_grid(1)%n_bins, & 
                                         & ' bins, with width ', & 
                                         & d%TE_grid(1)%res, & 
                                         & ' sec'
    End If
    Write(unit,*)
    Write(unit,'(A7,2A27)') 'Bin #','t-low [sec]','t-high [sec]'
    Write(unit,'(A7,2A27)') '-----','-------------------------','-------------------------'
    Do i = 1,d%TE_grid(1)%n_bins
        Write(unit,'(I7,2ES27.16E3)') i,d%TE_grid(1)%bounds(i-1),d%TE_grid(1)%bounds(i)
    End Do
    Write(unit,*)
    Write(unit,'(A,ES24.16E3,A,ES24.16E3,A)') '  Energy range: ',d%TE_grid(2)%min,' keV to ',d%TE_grid(2)%max,' keV'
    If (d%TE_grid(2)%log_spacing) Then
        Write(unit,'(A,I5,A,I5,A)') '  Logarithmically spaced bins: ', & 
                                  & d%TE_grid(2)%n_bins, & 
                                  & ' bins, with ', & 
                                  & d%TE_grid(2)%n_bins / d%TE_grid(2)%n_decades, & 
                                  & ' bins per decade'
    Else
        Write(unit,'(A,I7,A,ES24.16E3,A)') '  Linearly spaced bins: ', & 
                                         & d%TE_grid(2)%n_bins, & 
                                         & ' bins, with width ', & 
                                         & d%TE_grid(2)%res, & 
                                         & ' keV'
    End If
    Write(unit,'(A7,2A27)') 'Bin #','E-low [keV]','E-high [keV]'
    Write(unit,'(A7,2A27)') '-----','-------------------------','-------------------------'
    Do i = 1,d%TE_grid(2)%n_bins
        Write(unit,'(I7,2ES27.16E3)') i,d%TE_grid(2)%bounds(i-1),d%TE_grid(2)%bounds(i)
    End Do
    Write(unit,*)
    Write(unit,'(A)') '  Arrival Direction:  mu, cosine of angle from forward'
    Write(unit,'(A,I7,A,ES24.16E3)') '  Linearly spaced bins: ',d%Dir_grid(1)%n_bins,' bins, with width ',d%Dir_grid(1)%res
    Write(unit,'(A,ES24.16E3,A)') '                                       max width ', & 
                                & MaxVal( Abs( Acos(d%Dir_grid(1)%bounds(1:d%Dir_grid(1)%n_bins))- & 
                                             & Acos(d%Dir_grid(1)%bounds(0:d%Dir_grid(1)%n_bins-1)) ) & 
                                        & ), & 
                                & ' rad'
    Write(unit,'(A7,2A27)') 'Bin #','mu-low','mu-high'
    Write(unit,'(A7,2A27)') '-----','-------------------------','-------------------------'
    Do i = 1,d%Dir_grid(1)%n_bins
        Write(unit,'(I7,2ES27.16E3)') i,d%Dir_grid(1)%bounds(i-1),d%Dir_grid(1)%bounds(i)
    End Do
    Write(unit,*)
    Write(unit,'(A)') '  Arrival Direction:  omega, angle about downward'
    Write(unit,'(A,I7,A,ES24.16E3,A)') '  Linearly spaced bins: ',d%Dir_grid(2)%n_bins,' bins, with width ',d%Dir_grid(2)%res,' rad'
    Write(unit,'(A7,2A27)') 'Bin #','omega-low [rad]','omega-high [rad]'
    Write(unit,'(A7,2A27)') '-----','-------------------------','-------------------------'
    Do i = 1,d%Dir_grid(2)%n_bins
        Write(unit,'(I7,2ES27.16E3)') i,d%Dir_grid(2)%bounds(i-1),d%Dir_grid(2)%bounds(i)
    End Do
    Write(unit,*)
    Write(unit,'(A,I11)') '  Contrib Stack Size: ',d%Contrib_size
    If (d%shape_data) Then
        Write(unit,*)
        Write(unit,'(A,I5,A,I11,A)') '  Shape data collected in ', & 
                                   & d%n_slices, & 
                                   & ' time and energy slices, up to ', & 
                                   & d%max_shape_data, & 
                                   & ' contributions recorded'
        Write(unit,'(A13,A13,A19,A13)') 'Time Bin','Contribs','Energy Bin','Contribs'
        Write(unit,'(A13,A13,A19,A13)') '--------','--------','----------','--------'
        Do i = 1,d%n_slices
            Write(unit,'(I13,I13,I19,I13)') d%TE_grid(1)%slice_bin(i), & 
                                          & d%TE_grid(1)%slice_c(i), & 
                                          & d%TE_grid(2)%slice_bin(i), & 
                                          & d%TE_grid(2)%slice_c(i)
        End Do
    End If
    Write(unit,*)
    Write(unit,*)
    Close(unit)
End Subroutine Write_Detector

Function Bin_Center(g,b) Result(x)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: x
    Class(Grid_info_Type), Intent(In) :: g
    Integer, Intent(In) :: b

    If (g%log_spacing) Then
        x = 10._dp**( 0.5_dp * (Log10(g%bounds(b-1)) + Log10(g%bounds(b))) )
    Else
        x = 0.5_dp * (g%bounds(b-1) + g%bounds(b))
    End If
End Function Bin_Center

Subroutine Close_Slice_Files(n,flags1,units1,flags2,units2)
    Implicit None
    Integer, Intent(In) :: n
    Logical, Intent(In) :: flags1(1:n),flags2(1:n)
    Integer, Intent(In) :: units1(1:n),units2(1:n)
    Integer :: i

    Do i = 1,n
        If (flags1(i)) Close(units1(i))
        If (flags2(i)) Close(units2(i))
    End Do
End Subroutine

End Module Detectors
