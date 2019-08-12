#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Create the files with the information on the beam
//-------------------------------------------------------

void Create_single_normalisation_histogram(Int_t scan, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX

{
  // get the branches with intensity information
  // (also, set the name of the output file)
  char *file_name = new char[kg_string_size];
  TTree *tree = g_vdm_Tree;
  g_vdm_Tree->ResetBranchAddresses();
  Double_t *bunch1 = new Double_t[3564];  
  Double_t *bunch2 = new Double_t[3564];
  Double_t timerel;
  g_vdm_Tree->SetBranchAddress("timerel",&timerel);
  if (opt == 0) {
    g_vdm_Tree->SetBranchAddress("bunch1",bunch1);
    g_vdm_Tree->SetBranchAddress("bunch2",bunch2);
  } if (opt == 1) {
    g_vdm_Tree->SetBranchAddress("bptx1",bunch1);
    g_vdm_Tree->SetBranchAddress("bptx2",bunch2);
  }
  
  // get the histograms with DCCT information
  // two devices: A and B
  TH1D* DCCT1AH = (TH1D*) g_vdm_File->Get("Beam1A");
  TH1D* DCCT2AH = (TH1D*) g_vdm_File->Get("Beam2A");
  TH1D* DCCT1BH = (TH1D*) g_vdm_File->Get("Beam1B");
  TH1D* DCCT2BH = (TH1D*) g_vdm_File->Get("Beam2B");  

  // get indices for start and end of scan
  Int_t *indices = new Int_t [2];
  indices[0]=indices[1]=-1;
  FindIdicesOfScanStartEnd(scan, indices);

  // create histograms
  g_vdm_Tree->GetEntry(indices[0]);
  Double_t t0 = timerel;
  g_vdm_Tree->GetEntry(indices[1]);
  Double_t t1 = timerel;
  char *histo_name = new char[kg_string_size];
  sprintf(histo_name,"Norm1_%d_%d_%d",g_vdm_Fill,scan,opt);
  TH1D *h1 = new TH1D(histo_name,histo_name,((Int_t) (t1-t0)),t0,t1);
  sprintf(histo_name,"Norm2_%d_%d_%d",g_vdm_Fill,scan,opt);
  TH1D *h2 = new TH1D(histo_name,histo_name,((Int_t) (t1-t0)),t0,t1);
 
  // loop over tree to normalise each timerel
  for(Int_t i=indices[0];i<indices[1];i++) {
    // get relative time
    g_vdm_Tree->GetEntry(i);
    
    // varaibles to store correction factors
    Double_t cf_dcct_1a = 0;
    Double_t cf_dcct_2a = 0;  
    Double_t cf_dcct_1b = 0;
    Double_t cf_dcct_2b = 0;  

    // get the index of each histogram corresponding to the relative time
    Int_t idx_DCCT1AH = GetHistogramIndex(DCCT1AH,timerel);
    Int_t idx_DCCT1BH = GetHistogramIndex(DCCT1BH,timerel);  
    Int_t idx_DCCT2AH = GetHistogramIndex(DCCT2AH,timerel);
    Int_t idx_DCCT2BH = GetHistogramIndex(DCCT2BH,timerel);
    
    // get total intensity
    Double_t total_beam1 = 0;
    Double_t total_beam2 = 0;    
    for(Int_t j=0;j<3564;j++) { // get total intensity
      total_beam1 += bunch1[j];
      total_beam2 += bunch2[j];      
    }

    // get correction factor 
    cf_dcct_1a = ((DCCT1AH->GetBinContent(idx_DCCT1AH))/total_beam1);
    cf_dcct_1b = ((DCCT1BH->GetBinContent(idx_DCCT1BH))/total_beam1);
    cf_dcct_2a = ((DCCT2AH->GetBinContent(idx_DCCT2AH))/total_beam2);
    cf_dcct_2b = ((DCCT2BH->GetBinContent(idx_DCCT2BH))/total_beam2);    
    
    // correction factor is average over A and B factors
    Double_t cf_dcct_1 = 0.5*(cf_dcct_1a+cf_dcct_1b);
    Double_t cf_dcct_2 = 0.5*(cf_dcct_2a+cf_dcct_2b);

    // fill histos
    Int_t ibin=h1->FindBin(timerel);
    h1->SetBinContent(ibin,cf_dcct_1);
    h2->SetBinContent(ibin,cf_dcct_2);    

    /*
    // print info
    cout << " timerel " << timerel << endl;
    cout << "   cf_dcct_1 = " << cf_dcct_1 << " cf_dcct_2 = " << cf_dcct_2 << endl;
    */
  }

  // plot canvas
  sprintf(histo_name,"Canvas_%d_%d_%d",g_vdm_Fill,scan,opt);
  TCanvas *c = new TCanvas(histo_name,histo_name,900,600);
  c->Divide(1,2);
  c->cd(1);
  h1->SetTitle(";timerel;correction for beam 1");
  h1->Draw("p");
  c->cd(2);
  h2->SetTitle(";timerel;correction for beam 2");  
  h2->Draw("p");

  // clean
  delete [] histo_name;
  delete [] indices;
  delete [] bunch1; 
  delete [] bunch2; 


}

//-------------------------------------------------------
// Create normalisation histograms 
//-------------------------------------------------------

void Create_normalisation_histograms(Int_t Fill, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX

{
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create histograms for all scans
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
     Create_single_normalisation_histogram(i,opt);
  }
}
