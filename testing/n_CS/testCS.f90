Program testCS

Use Kinds, Only: dp
Use FileIO_utilities, Only: slash
Use n_Cross_Sections, Only:  CS_type
Use n_Cross_Sections, Only:  Setup_Cross_Sections

Implicit None

Type(CS_type) :: CS

CS = Setup_Cross_Sections( & 
                           & resources_directory = 'Resources/', & 
                           & cs_setup_file = 'n_CS_setup_All.txt', & 
                           & elastic_only = .FALSE., & 
                           & aniso_dist = .TRUE., & 
                           & E_min = 0._dp , E_max = 20000._dp, &
                           & verbose = .TRUE.)

End Program testCS
