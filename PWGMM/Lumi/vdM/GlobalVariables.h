//-------------------------------------------------------
// all global variables used in the code are defined here
// names of global variables start with g_
//-------------------------------------------------------

// the fill number to be analysed
Int_t g_vdm_Fill = -1;
// number of scans in the analysed fill
Int_t g_n_Scans_in_Fill = -1;

//-------------------------------------------------------
// indices at start and end of scans and steps
//-------------------------------------------------------
Int_t *g_Idx_Start_Scan_x = NULL;
Int_t *g_Idx_Start_Scan_y = NULL;
Int_t *g_Idx_End_Scan_x = NULL;
Int_t *g_Idx_End_Scan_y = NULL;

//-------------------------------------------------------
// To store names of input files
const Int_t kg_string_size = 1024;
char g_Input_vdm_File[kg_string_size];
char g_Input_vdm_DDL2_File[kg_string_size];
char g_Input_vdm_BPTX_File[kg_string_size]; 

//-------------------------------------------------------
// pointers to trees and files

TFile *g_vdm_File = NULL;
TFile *g_vdm_DDL2_File = NULL;
TFile *g_vdm_BPTX_File = NULL;

TTree *g_vdm_Tree = NULL;
TTree *g_vdm_DDL2_Tree = NULL;
TTree *g_vdm_BPTX_Tree = NULL;

//-------------------------------------------------------
// LHC frequency

const Double_t g_kLHCFreq = 11245.558; // Hz

//-------------------------------------------------------
// charge of beams (to be set to the real values in
// InputFromUser.h )

Double_t gBeamA = 0;
Double_t gBeamB = 0;

//-------------------------------------------------------
// The reference bunch intensities are not taken at one
// point in time but are the average of intensities in a
// range of half width given by
const Int_t gk_half_range = 0; // seconds

