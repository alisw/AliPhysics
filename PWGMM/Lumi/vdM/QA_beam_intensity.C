#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the beam intensity
//-------------------------------------------------------

//-------------------------------------------------------
// produce the canvas with the intensities at the
// normalisation point

void Canvas_intensity(Int_t Fill, Int_t opt, Int_t scan)
{
  // get number of bunch crossings
  Int_t nIBC = GetNumberInteractingBunchCrossings();

  // open the correct file
  char filename[120];
  if (opt==0) sprintf(filename,"../Fill-%d/FBCT_Scan_%d.root",Fill,scan);
  if (opt==1) sprintf(filename,"../Fill-%d/BPTX_Scan_%d.root",Fill,scan);
  TFile *IntensityFile = new TFile(filename);
  
  // get the info from the tree
  Double_t *bunch_intensity_1 = new Double_t[nIBC];
  Double_t *bunch_intensity_2 = new Double_t[nIBC];  
  TTree *intensity_tree = (TTree *) IntensityFile->Get("Bunch_Intensity");
  intensity_tree->SetBranchAddress("bunch_intensity_1",bunch_intensity_1);
  intensity_tree->SetBranchAddress("bunch_intensity_2",bunch_intensity_2);  intensity_tree->GetEntry(0);

  // fill histograms
  TH1D *bi1_H = new TH1D("bi1_H","bunch intensities beam 1",nIBC,-0.5,nIBC-0.5);
  TH1D *bi2_H = new TH1D("bi2_H","bunch intensities beam 2",nIBC,-0.5,nIBC-0.5);
  for(Int_t i=0;i<nIBC;i++){
    bi1_H->SetBinContent(i+1,bunch_intensity_1[i]);
    bi1_H->SetBinError(i+1,0);
    bi2_H->SetBinContent(i+1,bunch_intensity_2[i]);
    bi2_H->SetBinError(i+1,0);    
  }

  
  // plot histos
  char txt[120];
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  if (opt == 0) {
    sprintf(txt,"FBCT");
    bi1_H->SetTitle("beam 1;bunch crossing; beam 1 FBCT");
    bi2_H->SetTitle("beam 2;bunch crossing; beam 2 FBCT");    
  } else  {
    sprintf(txt,"BPTX");
    bi1_H->SetTitle("beam 1;bunch crossing; beam 1 BPTX");
    bi2_H->SetTitle("beam 2;bunch crossing; beam 2 BPTX");    
  }
  TCanvas *bi_C = new TCanvas(txt,txt,1200,800);

  bi_C->Divide(1,2);
  bi_C->cd(1);
  bi1_H->SetMarkerStyle(20);
  bi1_H->Draw("p");
  bi_C->cd(2);
  bi2_H->SetMarkerStyle(20);
  bi2_H->Draw("p");


  // clean
  delete [] bunch_intensity_1;
  delete [] bunch_intensity_2;  
}

//-------------------------------------------------------
Bool_t CheckInputs(Int_t opt, Int_t scan)
{
    // check inputs
  // --> opt
  if (opt != 0 && opt != 1) {
    cout << " Option has to be 0 (fBCT) or 1 (BPTX), try again!" << endl;
    return kFALSE;
  }
  // --> scan
  if ((scan <0) || (scan >=g_n_Scans_in_Fill)) {
    cout << " Invalid scan number. It has to be between zero and "
	 << g_n_Scans_in_Fill << " try again!" << endl;
    return kFALSE;
  }
  return kTRUE;
}

//-------------------------------------------------------

void QA_beam_intensity(Int_t Fill, Int_t opt, Int_t scan)
// opt = 0 => fBCT; opt = 1 => BPTX

{
  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // check
  if (!CheckInputs(opt, scan)) return;
  
  // make plot
  Canvas_intensity(Fill,opt,scan);
}
