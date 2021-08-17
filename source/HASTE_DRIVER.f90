!-------------------------------------------------------------------------------
!   High-Altitude to Space Transport Estimator for Neutrons (HASTEn)
!   
!   The High-Altitude to Space Transport Estimator (HASTE) is a high fidelity 
!   Monte Carlo code written in modern Fortran for estimating the radiation 
!   field seen by a space-based detector from a point source in or above the 
!   atmosphere.
!   
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
!   Version History
!   0.0    Initial version from HASTE-n by K.Mathews and HASTE-n-TE by W.Dailey
!   0.1    Initial release candidate, debugged & ready for studies.  Ready for
!          multithreading or coarray implementations.
!   0.2    Initial coarray implementation.
!   0.3    Revised coarray implementation for memory stability.  Revised next-
!          event orbital trajectory selection and supporting astro routines.
!   0.4    Converted from QuickWin to Console application.
!   0.5    Speed improvements.
!   0.6    Added direct (first-flight by sampling) computations. Revised cross 
!          section implementation for speed.
!   0.7    Revised EPL quadrature routines, Kepler Problem solver, & gravity 
!          divergence approach.
!   0.8    Opened development to community WG, migrated project to Github.
!   0.9    Major portability revision: gFortran compiler, IA32/IA64/ARM 
!          architectures now supported
!   0.10   Add albedo sources and limited-function adjoint solver routines.
!   0.11   Revise neutron cross sections data and procedures.
!-------------------------------------------------------------------------------
Program HASTE

Use Kinds, Only: dp  !double precision real
Use Kinds, Only: id  !double integer
Use Random_Numbers, Only: RNG_Type
Use Random_Numbers, Only: Setup_RNG
Use Setups, Only: Setup_HASTE
Use Setups, Only: Paths_Files_Type
Use Setups, Only: Create_Output_File_names
Use Atmospheres, Only: Atmosphere_Type
Use Atmospheres, Only: Setup_Atmosphere
Use Sources, Only: Source_Type
Use Sources, Only: Setup_Source
Use Detectors, Only: Detector_Type
Use Detectors, Only: Setup_Detector
Use Detectors, Only: Close_Slice_Files
Use Neutron_Scatter, Only: Scatter_Model_Type
Use Neutron_Scatter, Only: Setup_Scatter_Model
Use MC_Neutron, Only: Do_Neutron
Use Tallies, Only: Contrib_array
Use Tallies, Only: Setup_Tallies
Use Results, Only: Write_Run_Summary
Use Results, Only: Write_Tally_Grids
Use Statistics, Only: Std_Err
Use FileIO_Utilities, Only: max_line_len
Use FileIO_Utilities, Only: Worker_Index
Use FileIO_Utilities, Only: n_Workers
Use FileIO_Utilities, Only: Make_Boom
Use FileIO_Utilities, Only: creturn
Use FileIO_Utilities, Only: full_dash_line
Use FileIO_Utilities, Only: Delta_Time
# if CAF
    Use Setups, Only: Setup_info_to_disk
    Use Setups, Only: Setup_info_from_disk
    Use Setups, Only: Cleanup_temp_files
    Use Results, Only: Image_Result_to_Disk
    Use Results, Only: Image_Results_from_Disk
# endif

Implicit None

Character(max_line_len), Parameter :: title = 'High-Altitude to Space Transport Estimator for Neutrons (HASTEn)'
Character(max_line_len), Parameter :: ver =   '    v0.11.00, 18 Aug 2021'
Integer(id) :: n_histories
Logical :: absolute_n_histories  !specifies whether number of histories is an absolute limit or a target number of contributions   
Logical :: prompt_exit  !specifies whether the simulation will wait for user unput before exiting
Logical :: screen_progress  !determines if progress will be updated on screen during run
Type(Source_Type) :: source  !contains data defining neutron source
Type(Detector_Type) :: detector  !contains data defining detector
Type(Scatter_Model_Type) :: ScatterModel  !contains data defining the neutron scattering model
Type(Contrib_array) :: TE_Tallies  !list of contributions tallied to time-energy bins
Type(Contrib_array) :: Dir_Tallies  !list of contributions tallied to arrival direction bins
Type(Paths_Files_Type) :: paths_files  !contains path and file name character strings
Type(Atmosphere_Type) :: atmosphere  !contains data defining atmospheric representation
Type(RNG_Type) :: RNG  !contains data defining the random number generator
Integer(id) :: c_start,c_last,c_now  !clock counters for timing data
Real(dp) :: dt,ETTC  !times (in seconds) for computing estimated time to completion and recording run time
Integer :: HH,MM,SS  !integer times for computing estimated time to completion
Integer(id) :: n_p,p  !loop counter for histories
Logical :: contributed  !flag indicating that at least one contribution was accumulated for the current history
Real(dp) :: t_tot
Real(dp), Allocatable :: t_runs(:)
Real(dp), Allocatable :: t_waits(:)
Integer :: n_img,i_img
Integer(id) :: n_done
Integer(id), Allocatable :: n_hist_run(:),n_hist_hit(:)

n_img = n_Workers()
i_img = Worker_Index()
If (i_img .EQ. 1) Then
    Write(*,'(A)') full_dash_line
    Write(*,'(A)') title
    Write(*,'(A)') ver
    Write(*,'(A)') full_dash_line
    Write(*,'(A)') 'Setting up... '
    Call Setup_HASTE(prompt_exit,screen_progress,paths_files,n_histories,absolute_n_histories)
    paths_files%app_title = title
    paths_files%app_ver = ver
#   if CAF
        Write(*,'(I6,A)') n_img,' images sharing histories'
        !Write processed setup info to disk for other images
        Call Setup_Info_to_disk(n_histories,absolute_n_histories,prompt_exit,screen_progress,paths_files)
#   endif
End If
# if CAF
    !Sync so that setup info is available on disk for images other than 1
    SYNC ALL
    !Read setup info processed by image 1 from disk
    If (i_img .NE. 1) Then
        Call Setup_Info_from_disk(n_histories,absolute_n_histories,prompt_exit,screen_progress,paths_files)
        screen_progress = .FALSE.
    End If
# endif
RNG = Setup_RNG(paths_files%setup_file,paths_files%run_file_name)
atmosphere = Setup_Atmosphere( paths_files%setup_file,          & 
                             & paths_files%resources_directory, & 
                             & paths_files%run_file_name,       & 
                             & paths_files%cs_setup_file        )
ScatterModel = Setup_Scatter_Model( paths_files%setup_file,          & 
                                  & paths_files%resources_directory, & 
                                  & paths_files%cs_setup_file,       & 
                                  & paths_files%run_file_name,       & 
                                  & atmosphere%model_index           )
source = Setup_Source( paths_files%setup_file,          & 
                     & paths_files%resources_directory, & 
                     & paths_files%run_file_name,       & 
                     & paths_files%s_file_name,         & 
                     & atmosphere%R_top                 )
detector = Setup_Detector( paths_files%setup_file,          & 
                         & paths_files%resources_directory, & 
                         & paths_files%run_file_name,       & 
                         & paths_files%ss_file_name,        & 
                         & atmosphere%R_top                 )
TE_Tallies = Setup_Tallies(detector%TE_grid(1)%n_bins,detector%TE_grid(2)%n_bins)
Dir_Tallies = Setup_Tallies(detector%Dir_grid(1)%n_bins,detector%Dir_grid(2)%n_bins)
# if CAF
    !divide n_histories equally among images
    n_p = n_histories / n_img  !integer divide
    If (i_img .EQ. n_img) n_p = n_p + Mod(n_histories,INT(n_img,id))  !remainder added to last image
# else
    !all histories will be executed on a single image
    n_p = n_histories
# endif
!run a set of histories
If (i_img .EQ. 1) Then
    Write(*,'(A)') 'Running histories:  '
    Write(*,'(4A)') '  ',paths_files%results_directory,'...',paths_files%file_suffix
    If (screen_Progress) Then !prep screen for progress updates
        Write(*,'(A)') '    % Comp    ETTC     NE/min      H/min     % Eff'
        Write(*,'(A)') '    ------  --------  ---------  ---------  -------'
        !Initialize progress to screen
        Write(*,'(F9.2,A11,2ES11.3E2,F8.2,A1,A)',ADVANCE='NO') & 
              & 0._dp,'%  **:**:**',0._dp,0._dp,0._dp,'%',creturn
    End If
End If
p = 0_id
n_done = 0_id
!start a clock to track processing time
Call SYSTEM_CLOCK(c_start)
Do !run histories, periodically updating progress to the display if required
    If (screen_progress) Then  !update progress on screen
        If (Delta_Time(clock_then = c_last,clock_now = c_now) .GT. 1._dp) Then  !only update ETTC if elapsed time >1 second
            c_last = c_now
            If (p .GT. 0_id) Then
                dt = Delta_Time(clock_then = c_start)
                ETTC = dt / Real(p,dp) * Real(n_p-p,dp)  !time per particle so far multiplied by particles remaining
                HH = Floor(ETTC/3600._dp)
                MM = Floor((ETTC - Real(3600*HH,dp))/60._dp)
                SS = Floor(ETTC - Real(3600*HH,dp) - Real(60*MM,dp))
                Write(*,'(F9.2,A2,I3.2,A1,I2.2,A1,I2.2,2ES11.3E2,F8.2,A1,A)',ADVANCE='NO') & 
                      & 100._dp * Real(p,dp) / Real(n_p,dp),'% ', & !Percent Complete
                      & HH,':',MM,':',SS, & !ETTC
                      & Real(n_img*ScatterModel%next_events(1),dp) / (dt / 60._dp), & !Next-Events per minute
                      & Real(n_img*p,dp) / (dt / 60._dp), & !Histories per minute
                      & 100._dp * Real(p,dp) / Real(n_done,dp),'%', & !Percent efficiency
                      & creturn
            End If 
        End If
    End If
    !run a history
    Call Do_Neutron(source,detector,atmosphere,ScatterModel,RNG,contributed)
    n_done = n_done + 1_id
    If (contributed .OR. absolute_n_histories) Then
        p = p + 1_id
        !add the contributions from this history to the main time-energy and arrival direction contribution lists
        Call TE_Tallies%Tally_History( detector%TE_contrib_index,                                      & 
                                     & detector%TE_contribs_this_history(1:detector%TE_contrib_index), & 
                                     & detector%TE_grid(1)%n_bins,                                     & 
                                     & detector%TE_contribs_t,                                         & 
                                     & detector%TE_grid(2)%n_bins,                                     & 
                                     & detector%TE_contribs_E                                          )
        Call Dir_Tallies%Tally_History( detector%Dir_contrib_index,                                       & 
                                      & detector%Dir_contribs_this_history(1:detector%Dir_contrib_index), & 
                                      & detector%Dir_grid(1)%n_bins,                                      & 
                                      & detector%Dir_contribs_mu,                                         & 
                                      & detector%Dir_grid(2)%n_bins,                                      & 
                                      & detector%Dir_contribs_omega                                       )
    End If
    If (p .GE. n_p) Exit
End Do
!record final completion time
t_tot = Delta_Time(clock_then = c_start)
If (i_img.EQ.1 .AND. detector%shape_data) Call Close_Slice_Files( detector%n_slices,                 & 
                                                                & detector%TE_grid(1)%collect_shape, & 
                                                                & detector%TE_grid(1)%slice_unit,    & 
                                                                & detector%TE_grid(2)%collect_shape, & 
                                                                & detector%TE_grid(2)%slice_unit     )
# if CAF
    !Write image results to disk
    Call Image_Result_to_Disk( t_tot,                    & 
                             & RNG%t_wait,               & 
                             & n_p,                      & 
                             & n_done,                   & 
                             & TE_Tallies,               & 
                             & Dir_Tallies,              & 
                             & ScatterModel%n_kills,     & 
                             & ScatterModel%next_events, & 
                             & ScatterModel%n_no_tally,  & 
                             & ScatterModel%n_uncounted  )
    SYNC ALL
    If (i_img .EQ. 1) Then
        !Gather image results from temporary files to image 1
        Call Image_Results_from_Disk( detector%TE_grid(1)%n_bins,  & 
                                    & detector%TE_grid(2)%n_bins,  & 
                                    & detector%Dir_grid(1)%n_bins, & 
                                    & detector%Dir_grid(2)%n_bins, & 
                                    & t_runs,                      & 
                                    & t_waits,                     & 
                                    & n_hist_run,                  & 
                                    & n_hist_hit,                  & 
                                    & TE_Tallies,                  & 
                                    & Dir_Tallies,                 & 
                                    & ScatterModel%n_kills,        & 
                                    & ScatterModel%next_events,    & 
                                    & ScatterModel%n_no_tally,     & 
                                    & ScatterModel%n_uncounted     )
    End If
# else
    Allocate(t_runs(1:1))
    t_runs = t_tot
    Allocate(t_waits(1:1))
    t_waits = RNG%t_wait
    Allocate(n_hist_run(1:1))
    n_hist_run = n_done + ScatterModel%n_uncounted  !includes implicity leaked histories from exatmospheric sources
    Allocate(n_hist_hit(1:1))
    n_hist_hit = n_p
# endif
If (i_img .EQ. 1) Then
    If (screen_progress) Then !finalize progress to screen
        Write(*,'(F9.2,A2,I3.2,A1,I2.2,A1,I2.2,2ES11.3E2,F8.2,A1)') &
              & 100._dp,'% ', & !Percent Complete
              & 0,':',0,':',0, & !ETTC
              & Real(ScatterModel%next_events(1),dp) / (Sum(t_runs) / 60._dp), & !Next-Events per minute
              & Real(Sum(n_hist_run),dp) / (Sum(t_runs) / 60._dp), & !Histories per minute
              & 100._dp * Real(Sum(n_hist_hit),dp) / Real(Sum(n_hist_run),dp),'%' !Percent efficiency
    End If
    !Write results
    Write(*,'(A)', ADVANCE = 'NO') 'Writing output... '
    Call Write_Run_Summary( n_img,                    & 
                          & t_runs,                   & 
                          & t_waits,                  & 
                          & n_hist_hit,               & 
                          & n_hist_run,               & 
                          & RNG,                      & 
                          & paths_files,              & 
                          & atmosphere,               & 
                          & ScatterModel,             & 
                          & source,                   & 
                          & detector,                 & 
                          & TE_Tallies,               & 
                          & Dir_Tallies,              & 
                          & paths_files%log_file_name )
    Call Write_Tally_Grids( TE_Tallies,               & 
                          & Dir_Tallies,              & 
                          & detector,                 & 
                          & Sum(n_hist_run),          & 
                          & paths_files%F_file_name,  & 
                          & paths_files%TE_file_name, & 
                          & paths_files%t_file_name,  & 
                          & paths_files%E_file_name,  & 
                          & paths_files%d_file_name,  & 
                          & paths_files%m_file_name,  & 
                          & paths_files%o_file_name   )
#   if CAF
        !cleanup temp files
        Call Cleanup_temp_files()
#   endif
    Write(*,'(A)') 'Done.'
    Write(*,*)
End If
If (i_img .EQ. 1) Then
    Call Make_Boom()
    Write(*,'(A)') full_dash_line
    If (prompt_exit) Then
        Write(*,'(A)',ADVANCE='NO') 'Finished.  Press RETURN to continue...'
        Read(*,*)
    End If
End If

End Program HASTE
