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
Module Setups

    Implicit None
    Private
    Public :: Setup_HASTE
    Public :: Paths_Files_Type
    Public :: Write_Setup_Information
    Public :: Create_Output_File_names
#   if CAF
        Public :: Setup_Info_to_disk
        Public :: Setup_Info_from_disk
        Public :: Cleanup_temp_files
#   endif

    Type :: Paths_Files_Type
        Character(:), Allocatable :: app_title !name of program
        Character(:), Allocatable :: app_ver !version number of program
        Character(:), Allocatable :: program_exe !path and name of the running executable
        Character(:), Allocatable :: setup_file !specifies setup file, default is 'HASTE_Setup.txt' in the current working directory
        Character(:), Allocatable :: log_file_name !specifies log file, default is 'HASTE_log.txt' in the default results directory
        Character(:), Allocatable :: run_file_name !specifies run file, default is 'HASTE-setup.txt' in the default directory
        Character(:), Allocatable :: resources_directory !specifies resources dir, default is 'Resources' in current working dir
        Character(:), Allocatable :: cs_setup_file  !specifies setup file for cross sections
                                                    !this file sets atmospheric constituents and specifies what data is available
        Character(:), Allocatable :: results_directory !specifies results dir, default is 'Results' in current working dir
        Character(:), Allocatable :: file_suffix !specifies suffix to be appended to file names
        Character(:), Allocatable :: output_folder !specifies sub-folder in results folder,
                                                   !this variable is an interim holder, and is appended to results path during setup
        Character(:), Allocatable :: TE_file_name
        Character(:), Allocatable :: t_file_name
        Character(:), Allocatable :: E_file_name
        Character(:), Allocatable :: f_file_name
        Character(:), Allocatable :: d_file_name
        Character(:), Allocatable :: m_file_name
        Character(:), Allocatable :: o_file_name
        Character(:), Allocatable :: ss_file_name
        Character(:), Allocatable :: s_file_name
        Character(:), Allocatable :: cs_summary_file_name
    Contains
        Procedure, Pass :: Initialize => Initialize_Paths_Files
    End Type

Contains

Subroutine Setup_HASTE(prompt_for_exit,screen_progress,paths_files,n_neutron_histories,absolute_n_histories,setup_file)
    !Reads in problem data, method data, physics data
    !Initializes processes and variables
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use FileIO_Utilities, Only: Worker_Index
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Working_Directory
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Logical, Intent(Out) :: prompt_for_exit
    Logical, Intent(Out) :: screen_progress
    Type(paths_files_type), Intent(InOut) :: paths_files
    Integer(id), Intent(Out) :: n_neutron_histories
    Logical, Intent(Out) :: absolute_n_histories
    Character(*), Intent(In), Optional :: setup_file  !if specified, overrides default and command line setup file names
    Character(:), Allocatable :: file_suffix  !specifies suffix to be appended to file names
    Character(:), Allocatable :: output_folder  !specifies sub-folder in results folder,
                                                !this variable is an interim holder, and is appended to results path during setup
    Logical :: force_overwrite
    Integer :: setup_unit,stat,num_cmd_args
    Integer :: pathlen, status
    Character(max_path_len) :: path
    Logical :: interim_outputs
    Character(4) :: interim_unit
    Integer :: interim_period
    
    NameList /ProgramSetupList/   prompt_for_exit,screen_progress, &
                                & output_folder,file_suffix,force_overwrite, &
                                & interim_outputs,interim_unit,interim_period
    
    Call paths_files%Initialize()
    Allocate(Character(max_path_len) :: output_folder)
    Allocate(Character(max_path_len) :: file_suffix)
    !Get path and name of executable
    Call GET_COMMAND_ARGUMENT(0, path, pathlen, status)
    If (status .NE. 0) Call Output_Message('ERROR:  Setups: Setup_HASTE:  Read command line argument 0 failed',kill=.TRUE.)
    paths_files%program_exe = Trim(path)
    !default values for files and directories
    Call Working_Directory(GETdir=path,s=slash)
    paths_files%setup_file = Trim(path)//'HASTE_Setup.txt'
    paths_files%resources_directory = Trim(path)//'Resources'//slash
    paths_files%results_directory = Trim(path)//'Results'//slash
    If (Present(setup_file)) paths_files%setup_file = setup_file
    !open setup file and read namelist
    Open(NEWUNIT = setup_unit , FILE = paths_files%setup_file , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Setups: Setup_HASTE:  File open error, '//paths_files%setup_file// & 
                                        & ', IOSTAT=',stat,kill=.TRUE.)
    Read(setup_unit,NML = ProgramSetupList)
    Close(setup_unit)
    !Trim output folder and file suffix character strings
    output_folder = Trim(output_folder)
    file_suffix = Trim(file_suffix)
    If (Len(Trim(output_folder)) .GT. 0) Then
        paths_files%output_folder = Trim(output_folder)//slash
    Else
        paths_files%output_folder = Trim(output_folder)
    End If
    paths_files%file_suffix = Trim(file_suffix)
    !Check if files exist, if they do overwrite or prompt for revised file/folders from user
    Call Check_files_exist(force_overwrite,num_cmd_args,paths_files)
    !Check for correct directories, create them if they don't exist
    Call Check_folders_exist(paths_files)
    !append output folder to results directory
    paths_files%results_directory = Trim(paths_files%results_directory//paths_files%output_folder)
    !Construct file names
    paths_files%cs_setup_file = ''
    paths_files%log_file_name = ''
    paths_files%run_file_name = ''
    paths_files%TE_file_name = ''
    paths_files%t_file_name = ''
    paths_files%E_file_name = ''
    paths_files%f_file_name = ''
    paths_files%d_file_name = ''
    paths_files%m_file_name = ''
    paths_files%o_file_name = ''
    paths_files%ss_file_name = ''
    paths_files%s_file_name = ''
    Call Create_Output_File_names( paths_files%results_directory, & 
                                 & paths_files%file_suffix, & 
                                 & paths_files%log_file_name, & 
                                 & paths_files%TE_file_name, & 
                                 & paths_files%t_file_name, & 
                                 & paths_files%E_file_name, & 
                                 & paths_files%f_file_name, & 
                                 & paths_files%d_file_name, & 
                                 & paths_files%m_file_name, & 
                                 & paths_files%o_file_name, & 
                                 & paths_files%ss_file_name, & 
                                 & paths_files%s_file_name, & 
                                 & paths_files%cs_summary_file_name, & 
                                 & paths_files%run_file_name )
    !Create backup setup file in results folder and write namelist
    If (Worker_Index() .EQ. 1) Then
        !the backup setup file will not contain any continuation or study set configuration information
        Open( NEWUNIT = setup_unit , FILE = paths_files%run_file_name , STATUS = 'REPLACE' , ACTION = 'WRITE' , & 
            & POSITION = 'APPEND' , IOSTAT = stat )
        If (stat .NE. 0) Call Output_Message( 'ERROR:  Setups: Setup_HASTE:  File open error, '//paths_files%run_file_name// & 
                                            & ', IOSTAT=',stat,kill=.TRUE.)
        Write(setup_unit,NML = ProgramSetupList)
        Write(setup_unit,*)
        Close(setup_unit)
    End If
    !Read in the next setup namelist
    Call Setup_Estimator(paths_files%setup_file,paths_files%run_file_name,n_neutron_histories,absolute_n_histories)
End Subroutine Setup_HASTE

Subroutine Create_Output_File_names( dir,suff, & 
                                   & log_name, & 
                                   & TE_name, & 
                                   & T_name, & 
                                   & E_name, & 
                                   & F_name, & 
                                   & D_name, & 
                                   & M_name, & 
                                   & O_name, & 
                                   & SS_name, & 
                                   & S_name, & 
                                   & CS_name, & 
                                   & run_name )
    Implicit None
    Character(*), Intent(In) :: dir,suff
    Character(:), Allocatable, Intent(InOut) :: log_name
    Character(:), Allocatable, Intent(InOut) :: TE_name
    Character(:), Allocatable, Intent(InOut) :: T_name
    Character(:), Allocatable, Intent(InOut) :: E_name
    Character(:), Allocatable, Intent(InOut) :: F_name
    Character(:), Allocatable, Intent(InOut) :: D_name
    Character(:), Allocatable, Intent(InOut) :: M_name
    Character(:), Allocatable, Intent(InOut) :: O_name
    Character(:), Allocatable, Intent(InOut) :: SS_name
    Character(:), Allocatable, Intent(InOut) :: S_name
    Character(:), Allocatable, Intent(InOut) :: CS_name
    Character(:), Allocatable, Intent(InOut), Optional :: run_name
    
    !Construct file names
    !Construct log file name
    log_name = dir//'HASTE-Log'//suff//'.txt'
    If (Present(run_name)) run_name = dir//'HASTE-setup'//suff//'.txt'
    !Construct neutron output file names
    TE_name = dir//'TE_fluence'//suff//'.txt'
    T_name = dir//'T_fluence'//suff//'.txt'
    E_name = dir//'E_fluence'//suff//'.txt'
    F_name = dir//'tot_fluence'//suff//'.txt'
    D_name = dir//'ArrDirs'//suff//'.txt'
    M_name = dir//'ArrDirs_mu'//suff//'.txt'
    O_name = dir//'ArrDirs_omega'//suff//'.txt'
    SS_name = dir//'SliceShape'//suff
    S_name = dir//'Source'//suff//'.txt'
    CS_name = dir//'n_CS_summary'//suff//'.txt'
End Subroutine Create_Output_File_names

Subroutine Check_files_exist(overwrite,n_cmd_args,paths_files)
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Implicit None
    Logical, Intent(In) :: overwrite
    Integer, Intent(In) :: n_cmd_args
    Type(paths_files_type), Intent(InOut) :: paths_files
    Integer, Parameter :: n_names = 10
    Logical :: file_exists(1:n_names)
    Character(max_path_len) :: file_name(1:n_names)
    Character(64) :: user_inp
    Integer :: i

    Do_CheckFiles: Do
        !construct file names to check if they already exist
        file_name(1) = paths_files%results_directory//paths_files%output_folder//'TE_Fluence'//paths_files%file_suffix//'.txt'
        file_name(2) = paths_files%results_directory//paths_files%output_folder//'E_Fluence'//paths_files%file_suffix//'.txt'
        file_name(3) = paths_files%results_directory//paths_files%output_folder//'t_Fluence'//paths_files%file_suffix//'.txt'
        file_name(4) = paths_files%results_directory//paths_files%output_folder//'tot_Fluence'//paths_files%file_suffix//'.txt'
        file_name(5) = paths_files%results_directory//paths_files%output_folder//'ArrDirs'//paths_files%file_suffix//'.txt'
        file_name(6) = paths_files%results_directory//paths_files%output_folder//'ArrDirs_mu'//paths_files%file_suffix//'.txt'
        file_name(7) = paths_files%results_directory//paths_files%output_folder//'ArrDirs_omega'//paths_files%file_suffix//'.txt'
        file_name(8) = paths_files%results_directory//paths_files%output_folder//'HASTE-Log'//paths_files%file_suffix//'.txt'
        file_name(9) = paths_files%results_directory//paths_files%output_folder//'Source'//paths_files%file_suffix//'.txt'
        file_name(10) = paths_files%results_directory//paths_files%output_folder//'n_CS_summary'//paths_files%file_suffix//'.txt'
        Do i = 1,n_names
            INQUIRE(FILE = file_name(i) , EXIST = file_exists(i))
        End Do
        If (Any(file_exists)) Then
            Write(*,'(A)') 'WARNING:  One or more files already exist.'
            Write(*,'(A)') '    Output Folder: '//paths_files%results_directory//paths_files%output_folder
            Write(*,'(A)') '    File suffix:   '//paths_files%file_suffix
            If (overwrite) Then !file overwrite is forced by input settings
                Write(*,'(A)') 'WARNING:  Automatic overwrite selected by setup inputs.'
                Exit Do_CheckFiles
            Else If (n_cmd_args .GT. 0) Then  !command line specified path cannot be changed, prompt for overwrite or abort
                Write(*,'(A)') 'WARNING:  Command line specified results files already exist.'
                Do_cmd_line_overwrite: Do
                    Write(*,'(A)') 'Press:  O to overwrite'
                    Write(*,'(A)') '        A to abort the simulation'
                    Read (*,*) user_inp
                    Select Case (user_inp)
                        Case ('O','o')
                            Exit Do_CheckFiles
                        Case ('A','a')
                            Write(*,'(A)') '  Confirm:  Abort simulation? (Y/N)'
                            Read (*,*) user_inp
                            Select Case (user_inp)
                                Case ('Y','y')
                                    Write(*,*) 
                                    Write(*,'(A)') 'Simulation aborted by user.'
                                    ERROR STOP
                            End Select
                        Case Default
                            Write(*,'(A)') 'Invalid response.'
                    End Select
                End Do Do_cmd_line_overwrite
            End If
            Do_GetInput: Do
                Write(*,'(A)') 'Press:  O to overwrite'
                Write(*,'(A)') '        D to enter a new directory for the new files'
                Write(*,'(A)') '        S to enter a different file suffix for the new files'
                Write(*,'(A)') '        A to abort the simulation'
                Read (*,*) user_inp
                Select Case (user_inp)
                    Case ('O','o')
                        Exit Do_CheckFiles
                    Case ('D','d')
                        Write(*,'(A)') '  Enter new directory (omit ''Results'//slash//''', omit trailing '''//slash//'''):'
                        Read (*,*) user_inp
                        paths_files%output_folder = Trim(user_inp)
                        Exit Do_GetInput
                    Case ('S','s')
                        Write(*,'(A)') '  Enter new file suffix (recommend 1st character delimiter (''_'','' '',''-'',etc):'
                        Read (*,*) user_inp
                        paths_files%file_suffix = Trim(user_inp)
                        Exit Do_GetInput
                    Case ('A','a')
                        Write(*,'(A)') '  Confirm:  Abort simulation? (Y/N)'
                        Read (*,*) user_inp
                        Select Case (user_inp)
                            Case ('Y','y')
                                Write(*,*) 
                                Write(*,'(A)') 'Simulation aborted by user.'
                                ERROR STOP
                        End Select
                    Case Default
                        Write(*,'(A)') 'Invalid response.'
                End Select
            End Do Do_GetInput
        Else
            Exit Do_CheckFiles
        End If
    End Do Do_CheckFiles
End Subroutine Check_files_exist

Subroutine Check_folders_exist(paths_files)
    Use FileIO_Utilities, Only: Check_Directory
    Use FileIO_Utilities, Only: Create_Directory
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Type(paths_files_type), Intent(In) :: paths_files
    
    !Check if resources directories exist
    If (.NOT. Check_Directory(paths_files%resources_directory)) Then
        Call Output_Message( 'ERROR:  Setups: Setup_HASTE:  Resources directory not found: '// & 
                           & paths_files%resources_directory,kill=.TRUE. )
    End If
    !Make sure results directories exist
    Call Create_Directory(paths_files%results_directory)
    If (paths_files%output_folder .NE. '')  Call Create_Directory(paths_files%results_directory//paths_files%output_folder)
End Subroutine Check_folders_exist

Subroutine Setup_Estimator(setup_file_name,run_file_name,n_neutron_histories,absolute_n_histories)
    Use Kinds, Only: id
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Character(*), Intent(In) :: setup_file_name,run_file_name
    Integer(id), Intent(Out) :: n_neutron_histories
    Logical, Intent(Out) :: absolute_n_histories
    Integer :: setup_unit,stat
    
    NameList /EstimatorSetupList/ n_neutron_histories,absolute_n_histories
    
    !open setup file and read namelist
    Open(NEWUNIT = setup_unit , FILE = setup_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Setups: Setup_Estimator:  File open error, '//setup_file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE. )
    Read(setup_unit,NML = EstimatorSetupList)
    Close(setup_unit)
    Open(NEWUNIT = setup_unit , FILE = run_file_name , STATUS = 'OLD' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Setups: Setup_Estimator:  File open error, '//run_file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE. )
    Write(setup_unit,NML = EstimatorSetupList)
    Write(setup_unit,*)
    Close(setup_unit)
End Subroutine Setup_Estimator

# if CAF
Subroutine Setup_Info_to_disk(n_histories,abs_n_histories,prompt_for_exit,screen_progress,paths_files)
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Working_Directory
    Use FileIO_Utilities, Only: Create_Directory
    Use FileIO_Utilities, Only: Var_to_File
    Implicit None
    Integer(id), Intent(In) :: n_histories
    Logical, Intent(In) :: abs_n_histories
    Logical, Intent(In) :: prompt_for_exit
    Logical, Intent(In) :: screen_progress
    Type(Paths_Files_Type), Intent(In) :: paths_files
    Logical :: folder_exists,success
    Character(max_path_len) :: dir
    Character(:), Allocatable :: file_name,file_dir

    Call Working_Directory(GETdir = dir,s = slash)
    Allocate(Character(max_path_len) :: file_dir)
    file_dir = Trim(dir)//'temp'//slash
    !Make sure temp results directory exists
    Call Create_Directory(file_dir)
    Allocate(Character(max_path_len) :: file_name)
    !write n_histories to file
    file_name = file_dir//'n.tmp'
    Call Var_to_File(n_histories,file_name)
    !write absolute_n_histories to file
    file_name = file_dir//'abs_n.tmp'
    Call Var_to_File(abs_n_histories,file_name)
    !write prompt_for_exit to file
    file_name = file_dir//'eb.tmp'
    Call Var_to_File(prompt_for_exit,file_name)
    !write screen_progress to file
    file_name = file_dir//'sp.tmp'
    Call Var_to_File(screen_progress,file_name)
    !write each element of paths_files to file
    file_name = file_dir//'pf_app_t.tmp'
    Call Var_to_File(paths_files%app_title,file_name)
    file_name = file_dir//'pf_app_v.tmp'
    Call Var_to_File(paths_files%app_ver,file_name)
    file_name = file_dir//'pf_prog_exe.tmp'
    Call Var_to_File(paths_files%program_exe,file_name)
    file_name = file_dir//'pf_set_f.tmp'
    Call Var_to_File(paths_files%setup_file,file_name)
    file_name = file_dir//'pf_log_f.tmp'
    Call Var_to_File(paths_files%log_file_name,file_name)
    file_name = file_dir//'pf_run_f.tmp'
    Call Var_to_File(paths_files%run_file_name,file_name)
    file_name = file_dir//'pf_resource_d.tmp'
    Call Var_to_File(paths_files%resources_directory,file_name)
    file_name = file_dir//'pf_cs_f.tmp'
    Call Var_to_File(paths_files%cs_setup_file,file_name)
    file_name = file_dir//'pf_result_d.tmp'
    Call Var_to_File(paths_files%results_directory,file_name)
    file_name = file_dir//'pf_f_suff.tmp'
    Call Var_to_File(paths_files%file_suffix,file_name)
    file_name = file_dir//'pf_out_fold.tmp'
    Call Var_to_File(paths_files%output_folder,file_name)
    file_name = file_dir//'pf_TE_f.tmp'
    Call Var_to_File(paths_files%TE_file_name,file_name)
    file_name = file_dir//'pf_t_f.tmp'
    Call Var_to_File(paths_files%t_file_name,file_name)
    file_name = file_dir//'pf_E_f.tmp'
    Call Var_to_File(paths_files%E_file_name,file_name)
    file_name = file_dir//'pf_f_f.tmp'
    Call Var_to_File(paths_files%f_file_name,file_name)
    file_name = file_dir//'pf_d_f.tmp'
    Call Var_to_File(paths_files%d_file_name,file_name)
    file_name = file_dir//'pf_m_f.tmp'
    Call Var_to_File(paths_files%m_file_name,file_name)
    file_name = file_dir//'pf_o_f.tmp'
    Call Var_to_File(paths_files%o_file_name,file_name)
    file_name = file_dir//'pf_ss_f.tmp'
    Call Var_to_File(paths_files%ss_file_name,file_name)
    file_name = file_dir//'pf_s_f.tmp'
    Call Var_to_File(paths_files%s_file_name,file_name)
    file_name = file_dir//'pf_cs_f.tmp'
    Call Var_to_File(paths_files%cs_summary_file_name,file_name)
End Subroutine Setup_Info_to_disk
# endif

# if CAF
Subroutine Setup_Info_from_disk(n_histories,abs_n_histories,prompt_for_exit,screen_progress,paths_files)
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Working_Directory
    Use FileIO_Utilities, Only: Var_from_File
    Implicit None
    Integer(id), Intent(Out) :: n_histories
    Logical, Intent(Out) :: abs_n_histories
    Logical, Intent(Out) :: prompt_for_exit
    Logical, Intent(Out) :: screen_progress
    Type(Paths_Files_Type), Intent(InOut) :: paths_files
    Character(max_path_len) :: dir
    Character(:), Allocatable :: file_name,file_dir
    Character(max_path_len) :: C_tmp
    
    Call Working_Directory(GETdir=dir,s=slash)
    Allocate(Character(max_path_len) :: file_dir)
    file_dir = Trim(dir)//'temp'//slash
    Allocate(Character(max_path_len) :: file_name)
    !read n_histories from file
    file_name = file_dir//'n.tmp'
    Call Var_from_File(n_histories,file_name)
    !read abs_n_histories from file
    file_name = file_dir//'abs_n.tmp'
    Call Var_from_File(abs_n_histories,file_name)
    !read prompt_for_exit from file
    file_name = file_dir//'eb.tmp'
    Call Var_from_File(prompt_for_exit,file_name)
    !read screen_progress from file
    file_name = file_dir//'sp.tmp'
    Call Var_from_File(screen_progress,file_name)
    !initialize paths_files character strings
    Call paths_files%Initialize()
    !read each element of paths_files from file
    file_name = file_dir//'pf_app_t.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%app_title = Trim(C_tmp)
    file_name = file_dir//'pf_app_v.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%app_ver = Trim(C_tmp)
    file_name = file_dir//'pf_prog_exe.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%program_exe = Trim(C_tmp)
    file_name = file_dir//'pf_set_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%setup_file = Trim(C_tmp)
    file_name = file_dir//'pf_log_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%log_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_run_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%run_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_resource_d.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%resources_directory = Trim(C_tmp)
    file_name = file_dir//'pf_cs_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%cs_setup_file = Trim(C_tmp)
    file_name = file_dir//'pf_result_d.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%results_directory = Trim(C_tmp)
    file_name = file_dir//'pf_f_suff.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%file_suffix = Trim(C_tmp)
    file_name = file_dir//'pf_out_fold.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%output_folder = Trim(C_tmp)
    file_name = file_dir//'pf_TE_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%TE_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_t_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%t_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_E_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%E_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_f_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%f_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_d_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%d_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_m_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%m_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_o_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%o_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_ss_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%ss_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_s_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%s_file_name = Trim(C_tmp)
    file_name = file_dir//'pf_cs_f.tmp'
    Call Var_from_File(C_tmp,file_name)
    paths_files%cs_summary_file_name = Trim(C_tmp)
End Subroutine Setup_Info_from_disk
# endif

# if CAF
Subroutine Cleanup_temp_files()
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Working_Directory
    Use FileIO_Utilities, Only: Delete_Directory
    Implicit None
    Character(max_path_len) :: dir
    
    Call Working_Directory(GETdir=dir,s=slash)
    !delete the temp directory
    Call Delete_Directory(Trim(dir)//'temp')
End Subroutine Cleanup_temp_files
# endif

Subroutine Initialize_Paths_Files(paths_files)
    Use FileIO_Utilities, Only: max_path_len
    Implicit None
    Class(Paths_Files_Type), Intent(InOut) :: paths_files
    Character(9), Parameter :: empty_string = '<<EMPTY>>'

    !allocate character variables with an arbitrary length, each assignment statement then reallocates them to the correct length
    Allocate(Character(max_path_len) :: paths_files%app_title)
    paths_files%app_title = empty_string
    Allocate(Character(max_path_len) :: paths_files%app_ver)
    paths_files%app_ver = empty_string
    Allocate(Character(max_path_len) :: paths_files%program_exe)
    paths_files%program_exe = empty_string
    Allocate(Character(max_path_len) :: paths_files%setup_file)
    paths_files%setup_file = empty_string
    Allocate(Character(max_path_len) :: paths_files%resources_directory)
    paths_files%resources_directory = empty_string
    Allocate(Character(max_path_len) :: paths_files%results_directory)
    paths_files%results_directory = empty_string
    Allocate(Character(max_path_len) :: paths_files%cs_setup_file)
    paths_files%cs_setup_file = empty_string
    Allocate(Character(max_path_len) :: paths_files%output_folder)
    paths_files%output_folder = empty_string
    Allocate(Character(max_path_len) :: paths_files%file_suffix)
    paths_files%file_suffix = empty_string
    Allocate(Character(max_path_len) :: paths_files%log_file_name)
    paths_files%log_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%run_file_name)
    paths_files%run_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%TE_file_name)
    paths_files%TE_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%t_file_name)
    paths_files%t_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%E_file_name)
    paths_files%E_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%f_file_name)
    paths_files%f_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%d_file_name)
    paths_files%d_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%m_file_name)
    paths_files%m_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%o_file_name)
    paths_files%o_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%ss_file_name)
    paths_files%ss_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%s_file_name)
    paths_files%s_file_name = empty_string
    Allocate(Character(max_path_len) :: paths_files%cs_summary_file_name)
    paths_files%s_file_name = empty_string
End Subroutine Initialize_Paths_Files

Subroutine Write_Setup_Information(n_img,t_runs,t_waits,n_h_hit,n_h_run,RNG,paths_files,file_name)
    Use Kinds, Only: dp
    Use Kinds, Only: id
    Use Random_Numbers, Only: RNG_Type
    Use FileIO_Utilities, Only: Date_Time_string
    Use FileIO_Utilities, Only: Get_Host_Name
    Use FileIO_Utilities, Only: Output_Message
    Use FileIO_Utilities, Only: full_dash_line
    Use FileIO_Utilities, Only: half_dash_line
    Implicit None
    Integer, Intent(In) :: n_img
    Real(dp), Intent(In) :: t_runs(1:n_img)
    Real(dp), Intent(In) :: t_waits(1:n_img)
    Integer(id), Intent(In) :: n_h_hit(1:n_img)
    Integer(id), Intent(In) :: n_h_run(1:n_img)
    Type(RNG_Type), Intent(In) :: RNG
    Type(Paths_Files_Type), Intent(In) :: paths_files
    Character(*), Intent(In) :: file_name
    Character(80) :: hostname
    Integer :: unit,stat
    Integer :: i
    
    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message( 'ERROR:  Setups: Write_Setup_Information:  File open error, '//file_name// & 
                                        & ', IOSTAT=',stat,kill=.TRUE. )
    Write(unit,'(A)') full_dash_line
    Write(unit,'(A)') paths_files%app_title
    Write(unit,'(A)') paths_files%app_ver
    Write(unit,'(A)') full_dash_line
    Write(unit,'(A)') '   Copyright (C) 2017  Whitman T. Dailey'
    Write(unit,*)
    Write(unit,'(A)') '   This program is free software: you can redistribute it and/or modify'
    Write(unit,'(A)') '   it under the terms of the GNU General Public License version 3 as'
    Write(unit,'(A)') '   published by the Free Software Foundation.'
    Write(unit,*)
    Write(unit,'(A)') '   This program is distributed in the hope that it will be useful,'
    Write(unit,'(A)') '   but WITHOUT ANY WARRANTY; without even the implied warranty of'
    Write(unit,'(A)') '   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    Write(unit,'(A)') '   GNU General Public License for more details.'
    Write(unit,*)
    Write(unit,'(A)') '   You should have received a copy of the GNU General Public License'
    Write(unit,'(A)') '   along with this program.  If not, see <http://www.gnu.org/licenses/>.'
    Write(unit,'(A)') full_dash_line
    Write(unit,*)
    Write(unit,'(2A)') '  Run Complete: ',Date_Time_string()
    Write(unit,'(A,F11.3,A)') '  Total Compute Time: ',Sum(t_runs),' sec'
    Write(unit,'(A,F11.3,A)') '  Min Compute Time:   ',MinVal(t_runs),' sec'
    Write(unit,'(A,F11.3,A)') '  Max Compute Time:   ',MaxVal(t_runs),' sec'
    !spin time is computed as time spent waiting (sources of waiting are different run end times and random number generation)
    Write(unit,'(A,F7.2,A)') '  Spin Fraction:  ',100._dp*(Sum(t_waits)+Sum(MaxVal(t_runs)-t_runs)) / Sum(t_runs),'%'
    hostname = Get_Host_Name()
    If (n_img .GT. 1) Then
        Write(unit,'(A,I0,A)') '  Host: '//Trim(hostname)//', ',n_img,' coarray images'
    Else
        Write(unit,'(A)') '  Host: '//Trim(hostname)
    End If
    Write(unit,*)
    Write(unit,*)
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') 'SETUP INFORMATION'
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') '  Paths & Files:'
    Write(unit,'(2A)') '    Executable:     ',paths_files%program_exe
    Write(unit,'(2A)') '    Setup File:     ',paths_files%setup_file
    Write(unit,'(2A)') '    Resources Dir:  ',paths_files%resources_directory
    Write(unit,'(2A)') '    Cross Sect set: ',paths_files%resources_directory//paths_files%cs_setup_file
    Write(unit,'(2A)') '    Results Dir:    ',paths_files%results_directory
    If (Trim(paths_files%file_suffix) .EQ. '') Then
        Write(unit,'(2A)') '    File Suffix:    ','<<none>>'
    Else
        Write(unit,'(2A)') '    File Suffix:    ',paths_files%file_suffix
    End If
    Write(unit,'(2A)') '    Log File:       ',paths_files%log_file_name
    Write(unit,'(2A)') '    Setup File cpy: ',paths_files%run_file_name
    Write(unit,'(2A)') '    TE file:        ',paths_files%TE_file_name
    Write(unit,'(2A)') '    t file:         ',paths_files%t_file_name
    Write(unit,'(2A)') '    E file:         ',paths_files%E_file_name
    Write(unit,'(2A)') '    F file:         ',paths_files%F_file_name
    Write(unit,'(2A)') '    Arr Dirs file:        ',paths_files%d_file_name
    Write(unit,'(2A)') '    Arr Dirs mu file:     ',paths_files%m_file_name
    Write(unit,'(2A)') '    Arr Dirs omega file:  ',paths_files%o_file_name
    Write(unit,'(2A)') '    Slice Shape files:    ',paths_files%ss_file_name//'<<...>>.txt'
    Write(unit,'(2A)') '    Source Data file:     ',paths_files%s_file_name
    Write(unit,'(2A)') '    CS Summary file:      ',paths_files%cs_summary_file_name
    Write(unit,*)
    Write(unit,'(A,I11)') '  RNG Seed: ',RNG%seed
    Write(unit,'(A,I11)') '  RNG Array Length:    ',RNG%q_size
    Write(unit,'(A,I11)') '  RNG Array Used:      ',RNG%q_index - 1
    Write(unit,'(A,I11)') '  RNG Array Refreshes: ',RNG%q_refreshed
    Write(unit,'(A,I18)') '  Total R used:        ',RNG%q_refreshed * Int(RNG%q_size,id) + Int(RNG%q_index - 1,id)
    Write(unit,'(A,I0,A,I0,A,F6.2,A)') '  Number of Histories:  ', & 
                                     & Sum(n_h_hit), & 
                                     & ' contributing, ', & 
                                     & Sum(n_h_run), & 
                                     & ' total run, (', & 
                                     & 100._dp*Real(Sum(n_h_hit),dp)/Real(Sum(n_h_run),dp), & 
                                     & '% efficiency)'
    Write(unit,'(A)') '  Histories per image/thread:'
    Write(unit,'(A11,2A17)') 'Image','Contributing','Total Run'
    Write(unit,'(A11,2A17)') '-----','---------------','---------------'
    Do i = 1,n_img
        Write(unit,'(I11,2I17,A,F6.2,A)') i,n_h_hit(i),n_h_run(i),' (',100._dp*Real(n_h_hit(i),dp)/Real(n_h_run(i),dp),'%)'
    End Do
    Write(unit,*)
    Write(unit,*)
    Close(unit)
End Subroutine Write_Setup_Information

End Module Setups
