Program testCS

Use Kinds, Only: dp
Use FileIO_utilities, Only: slash
Use n_Cross_Sections, Only:  CS_type
Use n_Cross_Sections, Only:  Setup_Cross_Sections

Implicit None

Type(CS_type) :: CS

CS = Setup_Cross_Sections( & 
                           & resources_directory = 'Resources'//slash, & 
                           & cs_setup_file = 'n_CS_setup_All.txt', & 
                           & cs_summary_file = 'n_CS_summary.txt', & 
                           & elastic_only = .FALSE., & 
                           & aniso_dist = .TRUE., & 
                           & E_min = 0._dp , E_max = 20000._dp, &
                           & verbosity = 4)

End Program testCS
