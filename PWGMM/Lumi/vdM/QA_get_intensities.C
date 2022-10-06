#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// creates a canvas with the time dependence of the intensities of a given bc
//-------------------------------------------------------

void get_single_intensities(Int_t Fill, Int_t opt, Int_t bc)
// opt = 0 => fBCT; opt = 1 => BPTX

{
  // get the branches with intensity information
  // (also, set the name of the output file)
  char *file_name = new char[kg_string_size];
  TFile *NormFile = NULL;
  g_vdm_Tree->ResetBranchAddresses();
  Double_t *bunch1 = new Double_t[3564];  
  Double_t *bunch2 = new Double_t[3564];
  Double_t timerel;
  g_vdm_Tree->SetBranchAddress("timerel",&timerel);
  if (opt == 0) {
    g_vdm_Tree->SetBranchAddress("bunch1",bunch1);
    g_vdm_Tree->SetBranchAddress("bunch2",bunch2);
    sprintf(file_name,"../Fill-%d/FBCT_norm.root",g_vdm_Fill);    
  } if (opt == 1) {
    g_vdm_Tree->SetBranchAddress("bptx1",bunch1);
    g_vdm_Tree->SetBranchAddress("bptx2",bunch2);
    sprintf(file_name,"../Fill-%d/BPTX_norm.root",g_vdm_Fill);
  }
  
  // --> set up tree with info on normalisation
  NormFile = new TFile(file_name);      
  TTree *norm_tree = (TTree *) NormFile->Get("Beam_Normalisation");
  Double_t cf_dcct_1 = 0;
  Double_t cf_dcct_2 = 0;  
  norm_tree->ResetBranchAddresses();  
  norm_tree->SetBranchAddress("cf_dcct_1",&cf_dcct_1);
  norm_tree->SetBranchAddress("cf_dcct_2",&cf_dcct_2);
 
  // -- bunch indices
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  Int_t *bunches = new Int_t [nIBC];
  GetBunchIndices(bunches);
  // Int_t *BucketA = new Int_t [nIBC];
  // Int_t *BucketC = new Int_t [nIBC];
  Int_t BucketA[nIBC];
  Int_t BucketC[nIBC];  
  GetBucketInfo(BucketA,BucketC);

  Int_t idx1 = ((Int_t) (BucketA[bc]/10.0));
  Int_t idx2 = ((Int_t) (BucketC[bc]/10.0));

  // storage space
  Int_t ngr =  g_vdm_Tree->GetEntries();
  Double_t *beam1 = new Double_t[ngr];  
  Double_t *beam2 = new Double_t[ngr];
  Double_t *time = new Double_t[ngr];

  // loop over tree to print intensity
  Double_t imax = 0;
  Double_t imin = 1e20;
  for(Int_t i=0;i<ngr;i++) {
    g_vdm_Tree->GetEntry(i);
    norm_tree->GetEntry(i);

    // info on bunch intensities for the given bunch
    time[i] = timerel;
    beam1[i] = cf_dcct_1*bunch1[idx1]/gBeamA;
    beam2[i] = cf_dcct_2*bunch2[idx2]/gBeamB;
    if (imax<beam1[i]) imax=beam1[i];
    if (imax<beam2[i]) imax=beam2[i];
    if (imin>beam1[i]) imin=beam1[i];
    if (imin>beam2[i]) imin=beam2[i];
    /*
    if (i<10) cout << time[i] << " " << beam1[i] << " " << beam2[i]
		   << " " << (cf_dcct_1*beam1[i]) << " " << (cf_dcct_2*beam2[i]) << endl;
    else return;
    */
  }

  // graphs
  TGraph *gr1 = new TGraph(ngr,time,beam1);
  TGraph *gr2 = new TGraph(ngr,time,beam2);  
  gr1->SetMarkerStyle(20);  gr1->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(20);  gr2->SetMarkerColor(kBlue);  

  // plot
  TCanvas *cx = new TCanvas("cx","cx",1200,600);
  cx->Divide(1,1);
  cx->cd(1);
  TH1F* frame1 = gPad->DrawFrame(-1,0.1*imin,ngr+1,10*imax);
  frame1->GetYaxis()->SetTitleOffset(1.3);
  frame1->GetYaxis()->SetLabelSize(0.025);
  frame1->SetTitle("beam1 (red), beam2 (blue); timerel; beam intensity");
  gr1->Draw("p,same");
  gr2->Draw("p,same");
  cx->Print(Form("c1c_Fill%i_opt%i_bc%i.%s", Fill, opt, bc, FFormat));
  
   // clean
  delete [] beam1; 
  delete [] beam2; 
  delete [] bunch1; 
  delete [] bunch2; 
  delete [] bunches;
  // delete [] BucketA;
  // delete [] BucketC;

}

//-------------------------------------------------------
// Create normalisation histograms 
//-------------------------------------------------------

void QA_get_intensities(Int_t Fill, Int_t opt, Int_t bc)
// opt = 0 => fBCT; opt = 1 => BPTX

{
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create histograms for all scans
  get_single_intensities(Fill, opt, bc);
    
}
