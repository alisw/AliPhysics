#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the rate
//-------------------------------------------------------

void QA_rate_vs_sep(Int_t Fill, const char *rate_name, const char *rate_type,
		    const char *sep_type, Int_t scan_type, Int_t scan, Int_t bc)
// scan_type: 1 => x-scan; 2 => y-scan
{
  
  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // first get the files and trees
  // --> rate file
  char *file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_type,rate_name,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_type,rate_name,scan);
  TFile *rate_file = new TFile(file_name);
  // --> separation file
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sSep_x_Scan_%d.root",g_vdm_Fill,sep_type,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sSep_y_Scan_%d.root",g_vdm_Fill,sep_type,scan);
  TFile *sep_file = new TFile(file_name);
    // --> get the trees
  TTree *rate_tree = (TTree *) rate_file->Get("Rate");
  TTree *sep_tree = (TTree *) sep_file->Get("Separations");

  // Next step: prepare variables to be fitted
  // --> get number of separations
  Int_t n_sep = FindNumberSeparations(scan_type, scan);
  // --> reserve space for rates and separations
  Double_t *rate = new Double_t[n_sep];
  Double_t *rate_error = new Double_t[n_sep];  
  Double_t *sep = new Double_t[n_sep];
  // --> set branches
  rate_tree->ResetBranchAddresses();
  rate_tree->SetBranchAddress("rate",rate);
  rate_tree->SetBranchAddress("rate_error",rate_error);
  sep_tree->ResetBranchAddresses();
  sep_tree->SetBranchAddress("separation",sep);

  // get the right bunch crossing
  rate_tree->GetEntry(bc);
  sep_tree->GetEntry(bc);

  // fill graph
  TGraphErrors *gr = new TGraphErrors(n_sep,sep,rate,NULL,rate_error);
  gr->SetMarkerStyle(20);

  // define the limits for the plot
  // --> separation
  Double_t sep_min = 0;
  Double_t sep_max = 0;
  for(Int_t i=0;i<n_sep;i++) {
    if(sep[i]<sep_min) sep_min=sep[i];
    if(sep[i]>sep_max) sep_max=sep[i];
  }
  sep_min*=1.2;
  sep_max*=1.2;  
  // --> find maximum rate
  Double_t rate_max = 0;
  for(Int_t i=0;i<n_sep;i++) {
    if ((rate[i]+rate_error[i])>rate_max) rate_max = rate[i]+rate_error[i];
  }
  rate_max *= 1.2; //consider 20% larger limit

  // plot TGraph
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *rvss_C = new TCanvas("rvss_C","rate versus separation",600,400);
  TH1F* frame = gPad->DrawFrame(sep_min,0.001,sep_max,rate_max);
  frame->SetTitle(";separation (mm); rate (Hz)");
  gr->Draw("p,e1,same");
  
}
