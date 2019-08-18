#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the cross section
//-------------------------------------------------------

void QA_xs(Int_t Fill, const char *rate_name, const char *rate_type,
	   const char *sep_type, const char *intensity_type, Int_t fit_type, Int_t scan)
{
  
  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // get xs file and tree
  char *xs_file_name = new char[kg_string_size];
  sprintf(xs_file_name,"../Fill-%d/xs_%sRate_%s_%sSep_%s_Scan_%d_Fit_%s.root",
	  g_vdm_Fill,rate_type,rate_name,sep_type,intensity_type,scan,g_fit_model_name[fit_type]);
  TFile *xs_file = new TFile(xs_file_name);
  TTree *xs_tree = (TTree *) xs_file->Get("XS");

  // prepare tree branches
  Double_t xs=0;
  Double_t xs_error=0;  
  xs_tree->ResetBranchAddresses();
  xs_tree->SetBranchAddress("xs",&xs);
  xs_tree->SetBranchAddress("xs_error",&xs_error);    

  // get number of bunch crossings
  Int_t nIBC = GetNumberInteractingBunchCrossings();

  // reserve space
  Double_t *xs_all = new Double_t [nIBC];
  Double_t *xse_all = new Double_t [nIBC];
  TH1D *xs_H = new TH1D("xs_H","cross section",nIBC,-0.5,nIBC-0.5);
  
  // fill info
  for(Int_t i=0;i<nIBC;i++) {
    xs_tree->GetEntry(i);
    if (xs<0) continue;
    xs_H->SetBinContent(i+1,xs);
    xs_H->SetBinError(i+1,xs_error);
  }

  // plot histo

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *bi_C = new TCanvas("xs_histo","xs_histo",1200,800);
  xs_H->SetTitle("#sigma (mb);bunch crossing; #sigma (mb)");
  xs_H->SetMarkerStyle(20);
  xs_H->Draw("p,e");

  // clean
  delete [] xs_all;
  delete [] xse_all;  

  
}
