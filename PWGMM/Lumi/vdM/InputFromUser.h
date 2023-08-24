#ifndef __INPUT_FROM_USER_HH__
#define __INPUT_FROM_USER_HH__

//-------------------------------------------------------
// here the user defines the variables
// to be used as input for the analysis
//-------------------------------------------------------

void Set_input_file_names(Int_t Fill)
{
	if (Fill == 4937) // pp @ 13 TeV May 17, 2016
	{
		// set fill and number of scans in the fill
		g_vdm_Fill = Fill;
		g_n_Scans_in_Fill = 2;

		const char* CUT = "5.5m11.5_11.5p17.5";
		//cout <<" Following V0/T0 timing cut is being used: " <<CUT <<endl;
		// set name of input files
		sprintf(g_Input_vdm_File,      "../Fill-%d/vdm_time_4937_%s_1_v3.root",      g_vdm_Fill, CUT);
		sprintf(g_Input_vdm_DDL2_File, "../Fill-%d/vdm_DDL2_4937-%s_1_v3.root",      g_vdm_Fill, CUT);
		sprintf(g_Input_vdm_BPTX_File, "../Fill-%d/vdm_time_4937_%s_1_v3-BPTX.root", g_vdm_Fill, CUT);

		// charge of beams
		gBeamA = 1; // proton
		gBeamB = 1; // proton
	}
	else if (Fill == 5533) // p-Pb @ 8.16 TeV, Nov 23, 2016
	{
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
	}
	else if (Fill == 6012) // pp@13 TeV, July 27, 2017
	{
		// set fill and number of scans in the fill
		g_vdm_Fill = Fill;
		g_n_Scans_in_Fill = 2; 

		const char* CUT = "5.5m11.5_11.5p17.5";
		//cout <<" Following V0/T0 timing cut is being used: " <<CUT <<endl;
		sprintf(g_Input_vdm_File,      "../Fill-%d/vdm_time_6012_%s_1_v3.root",      g_vdm_Fill, CUT);
		sprintf(g_Input_vdm_DDL2_File, "../Fill-%d/vdm_DDL2_6012-%s_1_v3.root",      g_vdm_Fill, CUT);
		sprintf(g_Input_vdm_BPTX_File, "../Fill-%d/vdm_time_6012_%s_1_v3-BPTX.root", g_vdm_Fill, CUT);

		// charge of beams
		gBeamA = 1; // proton
		gBeamB = 1; // proton    
	}
	else if (Fill == 6864) // pp@13 June 2018
	{
		// set fill and number of scans in the fill
		g_vdm_Fill = Fill;
		g_n_Scans_in_Fill = 2; 

		const char* CUT1 = "5.5m11.5"; const char* CUT2 = "11.5p17.5";
		//cout <<"Following V0/T0 timing cut is being used: " <<CUT1 <<", " <<CUT2 <<endl;
		sprintf(g_Input_vdm_File,      "../Fill-%d/vdm_time_6864_%s_%s_1_v3.root",      g_vdm_Fill, CUT1, CUT2);
		sprintf(g_Input_vdm_DDL2_File, "../Fill-%d/vdm_DDL2_6864-%s.root",              g_vdm_Fill, CUT1);
		sprintf(g_Input_vdm_BPTX_File, "../Fill-%d/vdm_time_6864_%s_%s_1_v3-BPTX.root", g_vdm_Fill, CUT1, CUT2);

		// charge of beams
		gBeamA = 1; // proton
		gBeamB = 1; // proton    
	}
	else if (Fill == 7483) // Pb-Pb @ 5.02 TeV, Nov 29, 2018
	{
		// set fill and number of scans in the fill
		g_vdm_Fill = Fill;
		g_n_Scans_in_Fill = 2; 
		// set name of input files
		sprintf(g_Input_vdm_File,"../Fill-%d/vdm_time_7483_vstC-25feb19_1_v3.root",g_vdm_Fill);
		sprintf(g_Input_vdm_DDL2_File,"../Fill-%d/vdm_DDL2_7483-vstC-25feb19_1_v3.root",g_vdm_Fill);
		sprintf(g_Input_vdm_BPTX_File,"../Fill-%d/vdm_time_7483_vstC-25feb19_1_v3-BPTX.root",g_vdm_Fill);
		// charge of beams
		gBeamA = 82; // lead
		gBeamB = 82; // lead    
	}
	else
	{
		cout << " Fill " << Fill << " not know. Bye " << endl;
		exit(-100);
	}

	return;
}
#endif

