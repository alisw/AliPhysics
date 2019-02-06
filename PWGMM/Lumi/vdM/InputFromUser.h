
//-------------------------------------------------------
// here the user defines the variables
// to be used as input for the analysis
//-------------------------------------------------------



void Set_input_file_names(Int_t Fill)
{
  if (Fill == 4937) {
    // set fill and number of scans in the fill
    g_vdm_Fill = Fill;
    g_n_Scans_in_Fill = 2; 
    // set name of input files
    sprintf(g_Input_vdm_File,"../Fill-%d/vdm_time_4937_6m11_12p17_1_v3.root",g_vdm_Fill);
    sprintf(g_Input_vdm_DDL2_File,"../Fill-%d/vdm_DDL2_4937-6m11_12p17_1_v3.root",g_vdm_Fill);
    sprintf(g_Input_vdm_BPTX_File,"../Fill-%d/vdm_time_4937_6m11_12p17_1_v3-BPTX.root",g_vdm_Fill);
    // charge of beams
    gBeamA = 1; // proton
    gBeamB = 1; // lead    
  } else if (Fill == 5533) {
    // set fill and number of scans in the fill
    g_vdm_Fill = Fill;
    g_n_Scans_in_Fill = 2; 
    // set name of input files
    sprintf(g_Input_vdm_File,"../Fill-%d/vdm_time_5533_4m12_10p18_1_v3.root",g_vdm_Fill);
    sprintf(g_Input_vdm_DDL2_File,"../Fill-%d/vdm_DDL2_5533-4m12_10p18.root",g_vdm_Fill);
    sprintf(g_Input_vdm_BPTX_File,"../Fill-%d/vdm_time_5533_4m12_10p18_1_v3-BPTX.root",g_vdm_Fill);
    // charge of beams
    gBeamA = 1; // proton
    gBeamB = 82; // lead    
  } else {
    cout << " Fill " << Fill << " not know. Bye " << endl;
    exit(-100);
  }
}

