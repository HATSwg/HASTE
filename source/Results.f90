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
Module Results

    Implicit None
    Private

    Public :: Write_Run_Summary
    Public :: Write_Tally_Grids
#   if CAF
        Public :: Image_Result_to_Disk
        Public :: Image_Results_from_Disk
#   endif

Contains

# if CAF
Subroutine Image_Result_to_Disk(run_time,wait_time,n_hit,n_run,t,d,n_kills,next_events,no_tally,uncounted)
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use FileIO_Utilities, Only: slash
    Use Tallies, Only: Contrib_array
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: Working_Directory
    Use FileIO_Utilities, Only: Var_to_File
    Implicit None
    Real(dp), Intent(In) :: run_time
    Real(dp), Intent(In) :: wait_time
    Integer(id), Intent(In) :: n_hit
    Integer(id), Intent(In) :: n_run
    Type(Contrib_Array), Intent(In) :: t
    Type(Contrib_Array), Intent(In) :: d
    Integer(id), Intent(In) :: n_kills(1:6)
    Integer(id), Intent(In) :: next_events(1:3)
    Integer(id), Intent(In) :: no_tally(1:3)
    Integer(id), Intent(In) :: uncounted
    Character(3) :: i_char
    Character(max_path_len) :: dir
    Character(:), Allocatable :: file_name,file_dir

    Call Working_Directory(GETdir = dir,s = slash)
    Allocate(Character(max_path_len) :: file_dir)
    file_dir = Trim(dir)//'temp'//slash
    Write(i_char,'(I3.3)') this_image()
    Allocate(Character(max_path_len) :: file_name)
    !write elapsed time to file
    file_name = file_dir//'img'//i_char//'_t.tmp'
    Call Var_to_File(run_time,file_name)
    !write wait time to file
    file_name = file_dir//'img'//i_char//'_tw.tmp'
    Call Var_to_File(wait_time,file_name)
    !write number of histories hitting TE grid run to file
    file_name = file_dir//'img'//i_char//'_nh.tmp'
    Call Var_to_File(n_hit,file_name)
    !write number of histories run to file
    file_name = file_dir//'img'//i_char//'_nr.tmp'
    Call Var_to_File(n_run,file_name)
    !write time-energy tallies arrays and values to file
    file_name = file_dir//'img'//i_char//'_t_s.tmp'
    Call Var_to_File(t%size,file_name)
    file_name = file_dir//'img'//i_char//'_t_i.tmp'
    Call Var_to_File(t%index,file_name)
    file_name = file_dir//'img'//i_char//'_t_ts.tmp'
    Call Var_to_File(t%contribs(:)%i1,file_name)
    file_name = file_dir//'img'//i_char//'_t_Es.tmp'
    Call Var_to_File(t%contribs(:)%i2,file_name)
    file_name = file_dir//'img'//i_char//'_t_fs.tmp'
    Call Var_to_File(t%contribs(:)%f,file_name)
    file_name = file_dir//'img'//i_char//'_t_f2s.tmp'
    Call Var_to_File(t%contribs(:)%f_sq,file_name)
    file_name = file_dir//'img'//i_char//'_t_totf2.tmp'
    Call Var_to_File(t%tot_f_sq,file_name)
    file_name = file_dir//'img'//i_char//'_t_tf2.tmp'
    Call Var_to_File(t%f_sq_1,file_name)
    file_name = file_dir//'img'//i_char//'_t_Ef2.tmp'
    Call Var_to_File(t%f_sq_2,file_name)
    !write arrival direction tallies arrays and values to file
    file_name = file_dir//'img'//i_char//'_d_s.tmp'
    Call Var_to_File(d%size,file_name)
    file_name = file_dir//'img'//i_char//'_d_i.tmp'
    Call Var_to_File(d%index,file_name)
    file_name = file_dir//'img'//i_char//'_d_us.tmp'
    Call Var_to_File(d%contribs(:)%i1,file_name)
    file_name = file_dir//'img'//i_char//'_d_os.tmp'
    Call Var_to_File(d%contribs(:)%i2,file_name)
    file_name = file_dir//'img'//i_char//'_d_fs.tmp'
    Call Var_to_File(d%contribs(:)%f,file_name)
    file_name = file_dir//'img'//i_char//'_d_f2s.tmp'
    Call Var_to_File(d%contribs(:)%f_sq,file_name)
    file_name = file_dir//'img'//i_char//'_d_totf2.tmp'
    Call Var_to_File(d%tot_f_sq,file_name)
    file_name = file_dir//'img'//i_char//'_d_uf2.tmp'
    Call Var_to_File(d%f_sq_1,file_name)
    file_name = file_dir//'img'//i_char//'_d_of2.tmp'
    Call Var_to_File(d%f_sq_2,file_name)
    !write n_kills to file
    file_name = file_dir//'img'//i_char//'_n_k.tmp'
    Call Var_to_File(n_kills,file_name)
    !write next_events to file
    file_name = file_dir//'img'//i_char//'_n_e.tmp'
    Call Var_to_File(next_events,file_name)
    !write no_tally to file
    file_name = file_dir//'img'//i_char//'_n_t.tmp'
    Call Var_to_File(no_tally,file_name)
    !write uncounted to file
    file_name = file_dir//'img'//i_char//'_n_u.tmp'
    Call Var_to_File(uncounted,file_name)
End Subroutine Image_Result_to_Disk
# endif

# if CAF
Subroutine Image_Results_from_Disk(nt,nE,nm,no,t_runs,t_waits,n_hist_run,n_hist_hit,t,d,n_kills,next_events,no_tally,uncounted)
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use FileIO_Utilities, Only: slash
    Use Tallies, Only: Contrib_array
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: Working_Directory
    Use FileIO_Utilities, Only: Var_from_File
    Implicit None
    Integer, Intent(In) :: nt,nE
    Integer, Intent(In) :: nm,no
    Real(dp), Allocatable, Intent(Out) :: t_runs(:)
    Real(dp), Allocatable, Intent(Out) :: t_waits(:)
    Integer(id), Allocatable, Intent(Out) :: n_hist_run(:)
    Integer(id), Allocatable, Intent(Out) :: n_hist_hit(:)
    Type(Contrib_Array), Intent(InOut) :: t
    Type(Contrib_Array), Intent(InOut) :: d
    Integer(id), Intent(Out) :: n_kills(1:6)
    Integer(id), Intent(Out) :: next_events(1:3)
    Integer(id), Intent(Out) :: no_tally(1:3)
    Integer(id), Intent(Out) :: uncounted
    Integer :: i,j
    Integer :: size,index
    Real(dp) :: totf2
    Real(dp) :: tf2_1(1:nt),tf2_2(1:ne)
    Real(dp) :: df2_1(1:nm),df2_2(1:no)
    Real(dp), Allocatable :: fs(:),f2s(:)
    Real(dp), Allocatable :: TE_tmp(:,:,:)
    Real(dp), Allocatable :: Dir_tmp(:,:,:)
    Integer, Allocatable :: i1s(:),i2s(:)
    Integer(id) :: n_k(1:6),n_e(1:3),n_t(1:3),n_u
    Character(3) :: i_char
    Character(max_path_len) :: dir
    Character(:), Allocatable :: file_name,file_dir

    Call Working_Directory(GETdir = dir,s = slash)
    Allocate(Character(max_path_len) :: file_dir)
    file_dir = Trim(dir)//'temp'//slash
    !initialize time variables
    Allocate(t_runs(1:num_images()))
    t_runs = 0._dp
    Allocate(t_waits(1:num_images()))
    t_waits = 0._dp
    !initialize number of histories run
    Allocate(n_hist_hit(1:num_images()))
    n_hist_hit = -1_id
    Allocate(n_hist_run(1:num_images()))
    n_hist_run = -1_id
    !initialize time_energy tallies array and temp TE array
    t%size = 0
    t%index = 0
    Deallocate(t%contribs)
    t%tot_f_sq = 0._dp
    t%f_sq_1 = 0._dp
    t%f_sq_2 = 0._dp
    Allocate(TE_tmp(1:nt,1:ne,1:2))
    TE_tmp = 0._dp
    !initialize arrival direction tallies array and temp TE array
    d%size = 0
    d%index = 0
    Deallocate(d%contribs)
    d%tot_f_sq = 0._dp
    d%f_sq_1 = 0._dp
    d%f_sq_2 = 0._dp
    Allocate(Dir_tmp(1:nm,1:no,1:2))
    Dir_tmp = 0._dp
    !initialize counters
    n_kills = 0_id
    next_events = 0_id
    no_tally = 0_id
    uncounted = 0_id
    !read variables from files and accumulate appropriately
    Allocate(Character(max_path_len) :: file_name)
    Do i = 1,num_images()
        Write(i_char,'(I3.3)') i
        !read times
        file_name = file_dir//'img'//i_char//'_t.tmp'
        Call Var_from_File(t_runs(i),file_name,delete_file = .TRUE.)
        file_name = file_dir//'img'//i_char//'_tw.tmp'
        Call Var_from_File(t_waits(i),file_name,delete_file = .TRUE.)
        !read n_hist_run
        file_name = file_dir//'img'//i_char//'_nr.tmp'
        Call Var_from_File(n_hist_run(i),file_name,delete_file = .TRUE.)
        !read n_hist_hit
        file_name = file_dir//'img'//i_char//'_nh.tmp'
        Call Var_from_File(n_hist_hit(i),file_name,delete_file = .TRUE.)
        !read time-energy tallies arrays
        file_name = file_dir//'img'//i_char//'_t_s.tmp'
        Call Var_from_File(size,file_name,delete_file = .TRUE.)
        If (size .GT. t%size) t%size = size
        file_name = file_dir//'img'//i_char//'_t_i.tmp'
        Call Var_from_File(index,file_name,delete_file = .TRUE.)
        Allocate(i1s(1:index))
        i1s = -1
        Allocate(i2s(1:index))
        i2s = -1
        Allocate(fs(1:index))
        fs = 0._dp
        Allocate(f2s(1:index))
        f2s = 0._dp
        file_name = file_dir//'img'//i_char//'_t_ts.tmp'
        Call Var_from_File(i1s,file_name,delete_file = .TRUE.)
        file_name = file_dir//'img'//i_char//'_t_Es.tmp'
        Call Var_from_File(i2s,file_name,delete_file = .TRUE.)
        file_name = file_dir//'img'//i_char//'_t_fs.tmp'
        Call Var_from_File(fs,file_name,delete_file = .TRUE.)
        file_name = file_dir//'img'//i_char//'_t_f2s.tmp'
        Call Var_from_File(f2s,file_name,delete_file = .TRUE.)
        Do j = 1,index
            TE_tmp(i1s(j),i2s(j),1) = TE_tmp(i1s(j),i2s(j),1) + fs(j)
            TE_tmp(i1s(j),i2s(j),2) = TE_tmp(i1s(j),i2s(j),2) + f2s(j)
        End Do
        Deallocate(i1s)
        Deallocate(i2s)
        Deallocate(fs)
        Deallocate(f2s)
        file_name = file_dir//'img'//i_char//'_t_totf2.tmp'
        Call Var_from_File(totf2,file_name,delete_file = .TRUE.)
        t%tot_f_sq = t%tot_f_sq + totf2
        file_name = file_dir//'img'//i_char//'_t_tf2.tmp'
        Call Var_from_File(tf2_1,file_name,delete_file = .TRUE.)
        t%f_sq_1 = t%f_sq_1 + tf2_1
        file_name = file_dir//'img'//i_char//'_t_Ef2.tmp'
        Call Var_from_File(tf2_2,file_name,delete_file = .TRUE.)
        t%f_sq_2 = t%f_sq_2 + tf2_2
        !read arrival direction tallies arrays
        file_name = file_dir//'img'//i_char//'_d_s.tmp'
        Call Var_from_File(size,file_name,delete_file = .TRUE.)
        If (size .GT. d%size) d%size = size
        file_name = file_dir//'img'//i_char//'_d_i.tmp'
        Call Var_from_File(index,file_name,delete_file = .TRUE.)
        Allocate(i1s(1:index))
        i1s = -1
        Allocate(i2s(1:index))
        i2s = -1
        Allocate(fs(1:index))
        fs = 0._dp
        Allocate(f2s(1:index))
        f2s = 0._dp
        file_name = file_dir//'img'//i_char//'_d_us.tmp'
        Call Var_from_File(i1s,file_name,delete_file = .TRUE.)
        file_name = file_dir//'img'//i_char//'_d_os.tmp'
        Call Var_from_File(i2s,file_name,delete_file = .TRUE.)
        file_name = file_dir//'img'//i_char//'_d_fs.tmp'
        Call Var_from_File(fs,file_name,delete_file = .TRUE.)
        file_name = file_dir//'img'//i_char//'_d_f2s.tmp'
        Call Var_from_File(f2s,file_name,delete_file = .TRUE.)
        Do j = 1,index
            Dir_tmp(i1s(j),i2s(j),1) = Dir_tmp(i1s(j),i2s(j),1) + fs(j)
            Dir_tmp(i1s(j),i2s(j),2) = Dir_tmp(i1s(j),i2s(j),2) + f2s(j)
        End Do
        Deallocate(i1s)
        Deallocate(i2s)
        Deallocate(fs)
        Deallocate(f2s)
        file_name = file_dir//'img'//i_char//'_d_totf2.tmp'
        Call Var_from_File(totf2,file_name,delete_file = .TRUE.)
        d%tot_f_sq = d%tot_f_sq + totf2
        file_name = file_dir//'img'//i_char//'_d_uf2.tmp'
        Call Var_from_File(df2_1,file_name,delete_file = .TRUE.)
        d%f_sq_1 = d%f_sq_1 + df2_1
        file_name = file_dir//'img'//i_char//'_d_of2.tmp'
        Call Var_from_File(df2_2,file_name,delete_file = .TRUE.)
        d%f_sq_2 = d%f_sq_2 + df2_2
        !read counter values
        file_name = file_dir//'img'//i_char//'_n_k.tmp'
        Call Var_from_File(n_k,file_name,delete_file = .TRUE.)
        n_kills = n_kills + n_k
        file_name = file_dir//'img'//i_char//'_n_e.tmp'
        Call Var_from_File(n_e,file_name,delete_file = .TRUE.)
        next_events = next_events + n_e
        file_name = file_dir//'img'//i_char//'_n_t.tmp'
        Call Var_from_File(n_t,file_name,delete_file = .TRUE.)
        no_tally = no_tally + n_t
        file_name = file_dir//'img'//i_char//'_n_u.tmp'
        Call Var_from_File(n_u,file_name,delete_file = .TRUE.)
        uncounted = uncounted + n_u
        n_hist_run(i) = n_hist_run(i) + n_u  !includes implicity leaked histories from exatmospheric sources
    End Do
    !fill the time-energy tallies contribution array
    index = COUNT(TE_tmp(:,:,1) .GT. 0._dp)
    If (index .GT. t%size) t%size = index
    Allocate(t%contribs(1:t%size))
    Do i = 1,nt
        Do j = 1,ne
            If (TE_tmp(i,j,1) .GT. 0._dp) Then
                t%index = t%index + 1
                t%contribs(t%index)%i1 = i
                t%contribs(t%index)%i2 = j
                t%contribs(t%index)%f = TE_tmp(i,j,1)
                t%contribs(t%index)%f_sq = TE_tmp(i,j,2)
            End If
        End Do
    End Do
    Deallocate(TE_tmp)
    !fill the arrival direction tallies contribution array
    index = COUNT(Dir_tmp(:,:,1) .GT. 0._dp)
    If (index .GT. d%size) d%size = index
    Allocate(d%contribs(1:d%size))
    Do i = 1,nm
        Do j = 1,no
            If (Dir_tmp(i,j,1) .GT. 0._dp) Then
                d%index = d%index + 1
                d%contribs(d%index)%i1 = i
                d%contribs(d%index)%i2 = j
                d%contribs(d%index)%f = Dir_tmp(i,j,1)
                d%contribs(d%index)%f_sq = Dir_tmp(i,j,2)
            End If
        End Do
    End Do
    Deallocate(Dir_tmp)
End Subroutine Image_Results_from_Disk
# endif

Subroutine Write_Run_Summary(n_img,t_runs,t_waits,n_h_hit,n_h_run,RNG,paths_files,a,sm,s,d,TE_tallies,Dir_tallies,file_name)
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use Setups, Only: paths_files_type
    Use Setups, Only: Write_Setup_Information
    Use Random_Numbers, Only: RNG_type
    Use Atmospheres, Only: Atmosphere_Type
    Use Atmospheres, Only: Write_Atmosphere
    Use Neutron_Scatter, Only: Scatter_Model_Type
    Use Neutron_Scatter, Only: Write_Scatter_Model
    Use Sources, Only: Source_Type
    Use Sources, Only: Write_Source
    Use Detectors, Only: Detector_Type
    Use Detectors, Only: Write_Detector
    Use Tallies, Only: Contrib_array
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Integer, Intent(In) :: n_img
    Real(dp), Intent(In) :: t_runs(1:n_img)
    Real(dp), Intent(In) :: t_waits(1:n_img)
    Integer(id), Intent(In) :: n_h_hit(1:n_img)
    Integer(id), Intent(In) :: n_h_run(1:n_img)
    Type(RNG_Type), Intent(In) :: RNG
    Type(Paths_Files_Type), Intent(In) :: paths_files
    Type(Atmosphere_Type), Intent(In) :: a
    Type(Scatter_Model_Type), Intent(In) :: sm
    Type(Source_Type), Intent(In) :: s
    Type(Detector_Type), Intent(In) :: d
    Type(Contrib_array), Intent(In) :: TE_tallies
    Type(Contrib_array), Intent(In) :: Dir_tallies
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat

    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Run_Summary:  File open error, '//file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Close(unit)
    Call Write_Setup_Information(n_img,t_runs,t_waits,n_h_hit,n_h_run,RNG,paths_files,file_name)
    Call Write_Atmosphere(a,file_name)
    Call Write_Scatter_Model(sm,file_name)
    Call Write_Source(s,file_name)
    Call Write_Detector(d,file_name)
    Call Write_Results_summary(TE_tallies,Dir_tallies,d,Sum(n_h_run),file_name)
End Subroutine Write_Run_Summary

Subroutine Write_Tally_Grids( TE_list,Dir_list,d,n_h, & 
                            & F_file_name,TE_file_name,t_file_name,E_file_name,d_file_name,m_file_name,o_file_name)
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use Tallies, Only: Contrib_array
    Use Detectors, Only: Detector_Type
    Use Statistics, Only: Std_Err
    Use FileIO_Utilities, Only: Output_Message
   Implicit None
    Type(Contrib_array), Intent(In) :: TE_list
    Type(Contrib_array), Intent(In) :: Dir_list
    Type(Detector_Type), Intent(In) :: d
    Integer(id), Intent(In) :: n_h  !number of histories
    Character(*), Intent(In), Optional :: F_file_name
    Character(*), Intent(In), Optional :: TE_file_name
    Character(*), Intent(In), Optional :: t_file_name
    Character(*), Intent(In), Optional :: E_file_name
    Character(*), Intent(In), Optional :: d_file_name
    Character(*), Intent(In), Optional :: m_file_name
    Character(*), Intent(In), Optional :: o_file_name
    Integer :: unit,stat
    Real(dp) :: N_hist

    N_hist = Real(n_h,dp)
    If (Sum(TE_list%contribs(1:TE_list%index)%f) .GT. 0._dp) Then
        !write the total fluence file
        If (Present(F_file_name)) Then
            Open(NEWUNIT = unit , FILE = F_file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , IOSTAT = stat)
            If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Tally_Grids:  File open error, '//F_file_name// & 
                                                & ', IOSTAT=',stat,kill=.TRUE.)
            Call Write_tot_Tallies(unit,N_hist,TE_list%index,TE_list%contribs(1:TE_list%index),TE_list%tot_f_sq)
            Close(unit)
        End If
        !write the complete TE file
        If (Present(TE_file_name)) Then
            Open(NEWUNIT = unit , FILE = TE_file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , IOSTAT = stat)
            If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Tally_Grids:  File open error, '//TE_file_name// & 
                                                & ', IOSTAT=',stat,kill=.TRUE.)
            !write each time-energy bin's contribution and its standard error
            Call Write_2D_Tallies(unit,N_hist,TE_list%index,TE_list%contribs(1:TE_list%index),d%TE_grid(1),d%TE_grid(2))
            Close(unit)
        End If
        If (Present(t_file_name) .OR. Present(E_file_name)) Then
            If (Present(t_file_name)) Then
                !write the time integrated results to file
                Open(NEWUNIT = unit , FILE = t_file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , IOSTAT = stat)
                If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Tally_Grids:  File open error, '//t_file_name// & 
                                                    & ', IOSTAT=',stat,kill=.TRUE.)
                Call Write_1D_Tallies( unit, & 
                                     & N_hist, & 
                                     & TE_list%index, & 
                                     & TE_list%contribs(1:TE_list%index), & 
                                     & d%TE_grid(1), & 
                                     & d%TE_grid(2), & 
                                     & 1, & 
                                     & d%TE_grid(1)%n_bins, & 
                                     & TE_list%f_sq_1 )
                Close(unit)
            End If
            If (Present(E_file_name)) Then
                !write the energy integrated results to file
                Open(NEWUNIT = unit , FILE = E_file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , IOSTAT = stat)
                If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Tally_Grids:  File open error, '//E_file_name// & 
                                                    & ', IOSTAT=',stat,kill=.TRUE.)
                Call Write_1D_Tallies( unit, & 
                                     & N_hist, & 
                                     & TE_list%index, & 
                                     & TE_list%contribs(1:TE_list%index), & 
                                     & d%TE_grid(1), & 
                                     & d%TE_grid(2), & 
                                     & 2, & 
                                     & d%TE_grid(2)%n_bins, & 
                                     & TE_list%f_sq_2 )
                Close(unit)
            End If
        End If
    End If
    If (Sum(Dir_list%contribs(1:Dir_list%index)%f) .GT. 0._dp) Then
        !write the complete direction file
        If (Present(d_file_name)) Then
            Open(NEWUNIT = unit , FILE = d_file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , IOSTAT = stat)
            If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Tally_Grids:  File open error, '//d_file_name// & 
                                                & ', IOSTAT=',stat,kill=.TRUE.)
            !write each arrival direction bin's contribution and its standard error
            Call Write_2D_Tallies(unit,N_hist,Dir_list%index,Dir_list%contribs(1:Dir_list%index),d%Dir_grid(1),d%Dir_grid(2))
            Close(unit)
        End If
        If (Present(m_file_name) .OR. Present(o_file_name)) Then
            If (Present(m_file_name)) Then
                !write the time integrated results to file
                Open(NEWUNIT = unit , FILE = m_file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , IOSTAT = stat)
                If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Tally_Grids:  File open error, '//m_file_name// & 
                                                    & ', IOSTAT=',stat,kill=.TRUE.)
                Call Write_1D_Tallies( unit, & 
                                     & N_hist, & 
                                     & Dir_list%index, & 
                                     & Dir_list%contribs(1:Dir_list%index), & 
                                     & d%Dir_grid(1), & 
                                     & d%Dir_grid(2), & 
                                     & 1, & 
                                     & d%Dir_grid(1)%n_bins, & 
                                     & Dir_list%f_sq_1 )
                Close(unit)
            End If
            If (Present(o_file_name)) Then
                !write the energy integrated results to file
                Open(NEWUNIT = unit , FILE = o_file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , IOSTAT = stat)
                If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Tally_Grids:  File open error, '//o_file_name// & 
                                                    & ', IOSTAT=',stat,kill=.TRUE.)
                Call Write_1D_Tallies( unit, & 
                                     & N_hist, & 
                                     & Dir_list%index, & 
                                     & Dir_list%contribs(1:Dir_list%index), & 
                                     & d%Dir_grid(1), & 
                                     & d%Dir_grid(2), & 
                                     & 2, & 
                                     & d%Dir_grid(2)%n_bins, & 
                                     & Dir_list%f_sq_2 )
                Close(unit)
            End If
        End If
    End If
End Subroutine Write_Tally_Grids

Subroutine Write_tot_Tallies(unit,N_hist,index,contribs,tot_sq)
    Use Kinds, Only: dp
    Use Tallies, Only: Contrib_Quadruplet
    Use Statistics, Only: Std_Err
    Use Statistics, Only: std_devs_for_95CI
    Implicit None
    Integer, Intent(In) :: unit
    Real(dp), Intent(In) :: N_hist
    Integer, Intent(In) :: index
    Type(Contrib_Quadruplet), Intent(In) :: contribs(1:index)
    Real(dp), Intent(In) :: tot_sq
    Real(dp) :: f,err

    f = Sum(contribs(:)%f) / N_hist
    err = Std_Err(N_hist,Sum(contribs(:)%f),tot_sq)
    Write(unit,'(5ES27.16E3)') f, err, err/f, f-std_devs_for_95CI*err, f+std_devs_for_95CI*err
End Subroutine Write_tot_Tallies

Subroutine Write_2D_Tallies(unit,N_hist,index,contribs,grid1,grid2)
    Use Kinds, Only: dp
    Use Detectors, Only: Grid_Info_Type
    Use Tallies, Only: Contrib_Quadruplet
    Use Statistics, Only: Std_Err
    Use Statistics, Only: std_devs_for_95CI
    Implicit None
    Integer, Intent(In) :: unit
    Real(dp), Intent(In) :: N_hist
    Integer, Intent(In) :: index
    Type(Contrib_Quadruplet), Intent(In) :: contribs(1:index)
    Type(GRid_Info_Type), Intent(In) :: grid1,grid2
    Integer :: i,i1,i2
    Real(dp) :: min1,max1,mid1
    Real(dp) :: min2,max2,mid2
    Real(dp) :: f,err

    !write each bin's details, contribution, standard error, and confidence interval
    Do i = 1,index
        i1 = contribs(i)%i1
        min1 = grid1%bounds(i1-1)
        max1 = grid1%bounds(i1)
        mid1 = grid1%Bin_Center(i1)
        i2 = contribs(i)%i2
        min2 = grid2%bounds(i2-1)
        max2 = grid2%bounds(i2)
        mid2 = grid2%Bin_Center(i2)
        f = contribs(i)%f / N_hist
        If (f .LT. Tiny(f)) Cycle  !f has underflowed to zero, the tallied contribution is too small to report
        err = Std_Err(N_hist,contribs(i)%f,contribs(i)%f_sq)
        Write(unit,'(I11,3ES27.16E3,I14,8ES27.16E3)') i1, & 
                                                    & min1, & 
                                                    & max1, & 
                                                    & mid1, & 
                                                    & i2, & 
                                                    & min2, & 
                                                    & max2, & 
                                                    & mid2, & 
                                                    & f, & 
                                                    & err, & 
                                                    & err/f, & 
                                                    & f-std_devs_for_95CI*err, & 
                                                    & f+std_devs_for_95CI*err
    End Do
End Subroutine Write_2D_Tallies

Subroutine Write_1D_Tallies(unit,N_hist,index,contribs,grid1,grid2,dim,n,sq_list)
    Use Kinds, Only: dp
    Use Detectors, Only: Grid_Info_Type
    Use Tallies, Only: Contrib_Quadruplet
    Use Statistics, Only: Std_Err
    Use Statistics, Only: std_devs_for_95CI
    Implicit None
    Integer, Intent(In) :: unit
    Real(dp), Intent(In) :: N_hist
    Integer, Intent(In) :: index
    Type(Contrib_Quadruplet), Intent(In) :: contribs(1:index)
    Type(Grid_Info_Type), Intent(In) :: grid1,grid2
    Integer, Intent(In) :: dim
    Integer, Intent(In) :: n
    Real(dp), Intent(In) :: sq_list(1:n)
    Integer :: i,i1,i2
    Real(dp) :: min1,max1,mid1
    Real(dp) :: min2,max2,mid2
    Real(dp) :: f,err
    Real(dp), Allocatable :: tmp_grid(:,:)

    Allocate(tmp_grid(1:grid1%n_bins,1:grid2%n_bins))
    tmp_grid = 0._dp
    Do i = 1,index
        i1 = contribs(i)%i1
        i2 = contribs(i)%i2
        tmp_grid(i1,i2) = contribs(i)%f
    End Do
    If (dim .EQ. 1) Then
        Do i = 1,grid1%n_bins
            If ( Any(tmp_grid(i,:) .NE. 0._dp) ) Then
                min1 = grid1%bounds(i-1)
                max1 = grid1%bounds(i)
                mid1 = grid1%Bin_Center(i)
                f = Sum(tmp_grid(i,:)) / N_hist
                If (f .LT. Tiny(f)) Cycle  !f has underflowed to zero, the tallied contribution is too small to report
                err = Std_Err(N_hist,Sum(tmp_grid(i,:)),sq_list(i))
                Write(unit,'(I11,8ES27.16E3)') i, min1, max1, mid1, f, err, err/f, f-std_devs_for_95CI*err, f+std_devs_for_95CI*err
            End If
        End Do
    Else If (dim .EQ. 2) Then
        Do i = 1,grid2%n_bins
            If ( Any(tmp_grid(:,i) .NE. 0._dp) ) Then
                min2 = grid2%bounds(i-1)
                max2 = grid2%bounds(i)
                mid2 = grid2%Bin_Center(i)
                f = Sum(tmp_grid(:,i)) / N_hist
                If (f .LT. Tiny(f)) Cycle  !f has underflowed to zero, the tallied contribution is too small to report
                err = Std_Err(N_hist,Sum(tmp_grid(:,i)),sq_list(i))
                Write(unit,'(I11,8ES27.16E3)') i, min2, max2, mid2, f, err, err/f, f-std_devs_for_95CI*err, f+std_devs_for_95CI*err
            End If
        End Do
    End If
    Deallocate(tmp_grid)
End Subroutine Write_1D_Tallies

Subroutine Write_Results_summary(TE_list,Dir_List,d,n_h,file_name)
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use Tallies, Only: Contrib_array
    Use Detectors, Only: Detector_Type
    Use Statistics, Only: Std_Err
    Use FileIO_Utilities, Only: n_Workers
    Use FileIO_Utilities, Only: Output_Message
    Use FileIO_Utilities, Only: half_dash_line
    Implicit None
    Type(Contrib_array), Intent(In) :: TE_list
    Type(Contrib_array), Intent(In) :: Dir_list
    Type(Detector_Type), Intent(In) :: d
    Integer(id), Intent(In) :: n_h  !number of histories
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat
    Real(dp) :: N_hist
    Integer :: n_img

    N_hist = Real(n_h,dp)
    n_img = n_Workers()
    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Results: Write_Results:  File open error, '//file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') 'TALLIES'
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A,I11)')          '  Number of TE bins:             ',d%TE_grid(1)%n_bins * d%TE_grid(2)%n_bins
    Write(unit,'(A,I11,A,F6.2,A)') '  Number of TE bins w/ tallies:  ', & 
                                 & TE_list%index, & 
                                 & ' (', & 
                                 & 100._dp*Real(TE_list%index,dp)/Real(d%TE_grid(1)%n_bins*d%TE_grid(2)%n_bins,dp), & 
                                 & '% of detector grid)'
    If (n_img .EQ. 1) Write(unit,'(A,I11,A,F6.2,A)') '  TE bins list size (% used):    ', & 
                                                   & TE_list%size, & 
                                                   & ' (', & 
                                                   & 100._dp*Real(TE_list%index,dp)/Real(TE_list%size,dp), & 
                                                   & '%)'
    Write(unit,*)
    Write(unit,'(A,I11)')          '  Number of Dir bins:            ',d%Dir_grid(1)%n_bins * d%Dir_grid(2)%n_bins
    Write(unit,'(A,I11,A,F6.2,A)') '  Number of Dir bins w/ tallies: ', & 
                                 & Dir_list%index, & 
                                 & ' (', & 
                                 & 100._dp*Real(Dir_list%index,dp)/Real(d%Dir_grid(1)%n_bins*d%Dir_grid(2)%n_bins,dp), & 
                                 & '% of detector grid)'
    If (n_img .EQ. 1) Write(unit,'(A,I11,A,F6.2,A)') '  Dir bins list size (% used):   ', & 
                                                   & Dir_list%size, & 
                                                   & ' (', & 
                                                   & 100._dp*Real(Dir_list%index,dp)/Real(Dir_list%size,dp), & 
                                                   & '%)'
    Write(unit,*)
    Write(unit,*)
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') 'RESULTS'
    Write(unit,'(A)') half_dash_line
    If (Sum(TE_list%contribs(1:TE_list%index)%f).EQ.0._dp .AND. Sum(Dir_list%contribs(1:Dir_list%index)%f).EQ.0._dp) Then
        Write(unit,'(A)') 'No contributions recorded at detector.'
    Else
        Write(unit,'(A)') 'Total Fluence:'
        Write(unit,'(5A27)') 'F [n/km^2 per source n]','Std Err','Rel Std Err','95% CI Low','95% CI High'
        Write(unit,'(5A27)') '-----------------------','-----------------------','-----------------------', & 
                           & '-----------------------','-----------------------'
        Call Write_tot_Tallies(unit,N_hist,TE_list%index,TE_list%contribs(1:TE_list%index),TE_list%tot_f_sq)
        Write(unit,*)
        Write(unit,'(A)') 'Fluence by Time-Energy bin:'
        Write(unit,'(A11,3A27,A14,8A27)') 'Time bin','t-low [s]','t-high [s]','t-mid [s]', & 
                                        & 'Energy bin','E-low [keV]','E-high [keV]','E-mid [keV]', & 
                                        & 'F [n/km^2 per source n]','Std Err','Rel Std Err','95% CI Low','95% CI High'
        Write(unit,'(A11,3A27,A14,8A27)') '----------','-----------------------','-----------------------', & 
                                        & '-----------------------','----------','-----------------------', & 
                                        & '-----------------------','-----------------------','-----------------------', & 
                                        & '-----------------------','-----------------------','-----------------------', & 
                                        & '-----------------------'
        Call Write_2D_Tallies(unit,N_hist,TE_list%index,TE_list%contribs(1:TE_list%index),d%TE_grid(1),d%TE_grid(2))
        Write(unit,*)
        Write(unit,'(A)') 'Fluence by Time bin:'
        Write(unit,'(A11,8A27)') 'Time bin','t-low [s]','t-high [s]','t-mid [s]', & 
                               & 'F [n/km^2 per source n]','Std Err','Rel Std Err','95% CI Low','95% CI High'
        Write(unit,'(A11,8A27)') '----------','-----------------------','-----------------------','-----------------------', & 
                               & '-----------------------','-----------------------','-----------------------', & 
                               & '-----------------------','-----------------------'
        Call Write_1D_Tallies( unit, & 
                             & N_hist, & 
                             & TE_list%index, & 
                             & TE_list%contribs(1:TE_list%index), & 
                             & d%TE_grid(1), & 
                             & d%TE_grid(2), & 
                             & 1, & 
                             & d%TE_grid(1)%n_bins, & 
                             & TE_list%f_sq_1 )
        Write(unit,*)
        Write(unit,'(A)') 'Fluence by Energy bin:'
        Write(unit,'(A11,8A27)') 'Energy bin','E-low [keV]','E-high [keV]','E-mid [keV]', & 
                               & 'F [n/km^2 per source n]','Std Err','Rel Std Err','95% CI Low','95% CI High'
        Write(unit,'(A11,8A27)') '----------','-----------------------','-----------------------','-----------------------', & 
                               & '-----------------------','-----------------------','-----------------------', & 
                               & '-----------------------','-----------------------'
        Call Write_1D_Tallies( unit, & 
                             & N_hist, & 
                             & TE_list%index, & 
                             & TE_list%contribs(1:TE_list%index), & 
                             & d%TE_grid(1), & 
                             & d%TE_grid(2), & 
                             & 2, & 
                             & d%TE_grid(2)%n_bins, & 
                             & TE_list%f_sq_2 )
        Write(unit,*)
        Write(unit,'(A)') 'Fluence by Arrival Direction Bin:'
        Write(unit,'(A11,3A27,A14,8A27)') 'Mu bin','mu-low','mu-high','mu-mid', & 
                                        & 'Omega bin','omega-low [rad]','omega-high [rad]','omega-mid [rad]', & 
                                        & 'F [n/km^2 per source n]','Std Err','Rel Std Err','95% CI Low','95% CI High'
        Write(unit,'(A11,3A27,A14,8A27)') '----------','-----------------------','-----------------------', & 
                                        & '-----------------------','----------','-----------------------', & 
                                        & '-----------------------','-----------------------','-----------------------', & 
                                        & '-----------------------','-----------------------','-----------------------', & 
                                        & '-----------------------'
        Call Write_2D_Tallies(unit,N_hist,Dir_list%index,Dir_list%contribs(1:Dir_list%index),d%Dir_grid(1),d%Dir_grid(2))
        Write(unit,*)
        Write(unit,'(A)') 'Fluence by Arrival Direction Mu bin:'
        Write(unit,'(A11,8A27)') 'Mu bin','mu-low','mu-high','mu-mid', & 
                               & 'F [n/km^2 per source n]','Std Err','Rel Std Err','95% CI Low','95% CI High'
        Write(unit,'(A11,8A27)') '----------','-----------------------','-----------------------','-----------------------', & 
                               & '-----------------------','-----------------------','-----------------------', & 
                               & '-----------------------','-----------------------'
        Call Write_1D_Tallies( unit, & 
                             & N_hist, & 
                             & Dir_list%index, & 
                             & Dir_list%contribs(1:Dir_list%index), & 
                             & d%Dir_grid(1), & 
                             & d%Dir_grid(2), & 
                             & 1, & 
                             & d%Dir_grid(1)%n_bins, & 
                             & Dir_list%f_sq_1 )
        Write(unit,*)
        Write(unit,'(A)') 'Fluence by Arrival Direction Omega bin:'
        Write(unit,'(A11,8A27)') 'Omega bin','omega-low [rad]','omega-high [rad]','omega-mid [rad]', & 
                               & 'F [n/km^2 per source n]','Std Err','Rel Std Err','95% CI Low','95% CI High'
        Write(unit,'(A11,8A27)') '----------','-----------------------','-----------------------','-----------------------', & 
                               & '-----------------------','-----------------------','-----------------------', & 
                               & '-----------------------','-----------------------'
        Call Write_1D_Tallies( unit, & 
                             & N_hist, & 
                             & Dir_list%index, & 
                             & Dir_list%contribs(1:Dir_list%index), & 
                             & d%Dir_grid(1), & 
                             & d%Dir_grid(2), & 
                             & 2, & 
                             & d%Dir_grid(2)%n_bins, & 
                             & Dir_list%f_sq_2 )
        Write(unit,*)
        Write(unit,*)
    End If
    Close(unit)
End Subroutine Write_Results_summary

End Module Results
