#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Macro to print the middle timestamp in each step of a vdm scan
//-------------------------------------------------------

void QA_print_timestamps(Int_t Fill, Int_t scan)
{
	// initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// first get the files and trees
	TFile *sepx_file = new TFile(Form("../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan));
	TFile *sepy_file = new TFile(Form("../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan));
	TTree *sepx_tree = (TTree *) sepx_file->Get("SepInfo");
	TTree *sepy_tree = (TTree *) sepy_file->Get("SepInfo");

	// Next step: prepare variables
	Long64_t time_separation_start_x;
	Long64_t time_separation_end_x;
	Long64_t time_separation_start_y;
	Long64_t time_separation_end_y;
	// --> set branches
	sepx_tree->ResetBranchAddresses();
	sepx_tree->SetBranchAddress("time_separation_start",&time_separation_start_x);
	sepx_tree->SetBranchAddress("time_separation_end",&time_separation_end_x);
	sepy_tree->ResetBranchAddresses();
	sepy_tree->SetBranchAddress("time_separation_start",&time_separation_start_y);
	sepy_tree->SetBranchAddress("time_separation_end",&time_separation_end_y);

	// Final step: printout
	// --> open file
	ofstream ofs;
	ofs.open(Form("../Fill-%d/Timestamps_Scan_%d.txt",g_vdm_Fill,scan));
	// --> loop over tree entries
	const Int_t nx = sepx_tree->GetEntries();
	for (Int_t i = 0; i<nx; i++)
	{
		sepx_tree->GetEntry(i);
		Long_t time = time_separation_start_x + 0.5*(time_separation_end_x-time_separation_start_x);
		ofs << time << endl;
		if(i==0) cout << " start time x-scan " << time_separation_start_x << endl;
		if(i==(nx-1)) cout << " end time x-scan " << time_separation_end_x << endl;    
	}
	const Int_t ny = sepy_tree->GetEntries();
	for (Int_t i = 0; i<ny; i++)
	{
		sepy_tree->GetEntry(i);
		Long_t time = time_separation_start_y + 0.5*(time_separation_end_y-time_separation_start_y);
		ofs << time << endl;
		if(i==0) cout << " start time y-scan " << time_separation_start_y << endl;
		if(i==(nx-1)) cout << " end time y-scan " << time_separation_end_y << endl;    
	}
	// --> close file
	ofs.close();

	return;
}
