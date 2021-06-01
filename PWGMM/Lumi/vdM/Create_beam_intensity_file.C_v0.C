#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// get indices between x and y scans
//-------------------------------------------------------
void GetIndices(Int_t *idx)
{
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    idx[i] = FindIdxBetweenScans(i);			 
  }
}

//-------------------------------------------------------
// Create the files with the information on the beam
//-------------------------------------------------------

void Create_single_beam_intensity_file(Int_t scan, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX

// first compute the normalization factors
// then compute the intensities at the reference point
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
    sprintf(file_name,"../Fill-%d/FBCT_Scan_%d.root",g_vdm_Fill,scan);
  } if (opt == 1) {
    g_vdm_Tree->SetBranchAddress("bptx1",bunch1);
    g_vdm_Tree->SetBranchAddress("bptx2",bunch2);
    sprintf(file_name,"../Fill-%d/BPTX_Scan_%d.root",g_vdm_Fill,scan);
  }


  // get the histograms with DCCT information
  // two devices: A and B
  TH1D* DCCT1AH = (TH1D*) g_vdm_File->Get("Beam1A");
  TH1D* DCCT2AH = (TH1D*) g_vdm_File->Get("Beam2A");
  TH1D* DCCT1BH = (TH1D*) g_vdm_File->Get("Beam1B");
  TH1D* DCCT2BH = (TH1D*) g_vdm_File->Get("Beam2B");  

  // varaibles to store correction factors
  Double_t cf_dcct_1a = 0;
  Double_t cf_dcct_2a = 0;  
  Double_t cf_dcct_1b = 0;
  Double_t cf_dcct_2b = 0;  

  // obtain the normalization at different times
  // around the index chosen to normalize and then average
  Int_t idx = FindIdxBetweenScans(scan);
  Double_t counter = 0.0;
  Int_t start = idx-gk_half_range;
  Int_t end = idx+gk_half_range;
  Int_t i = start;
  while(i<=end) {
    g_vdm_Tree->GetEntry(i);
    counter += 1.0;
    Double_t total_beam1 = 0;
    Double_t total_beam2 = 0;    
    for(Int_t j=0;j<3564;j++) { // get total intensity
      total_beam1 += bunch1[j];
      total_beam2 += bunch2[j];      
    }
    // get correction factor 
    cf_dcct_1a += ((DCCT1AH->GetBinContent(2*i+1))/total_beam1);
    cf_dcct_1b += ((DCCT1BH->GetBinContent(2*i+1))/total_beam1);
    cf_dcct_2a += ((DCCT2AH->GetBinContent(2*i+1))/total_beam2);
    cf_dcct_2b += ((DCCT2BH->GetBinContent(2*i+1))/total_beam2);    
    i++;
    cout << " time rel " << timerel << " bin center " << DCCT1AH->GetBinCenter(2*i+1) << endl;
  }
  // average
  cf_dcct_1a /= counter;
  cf_dcct_2a /= counter;
  cf_dcct_1b /= counter;
  cf_dcct_2b /= counter;  
  cout << " cf_dcct_1a = " << cf_dcct_1a << " cf_dcct_1b = " << cf_dcct_1b << endl
       << " cf_dcct_2a = " << cf_dcct_2a << " cf_dcct_2b = " << cf_dcct_2b << endl;
  // correction factor is average over A and B factors
  Double_t cf_dcct_1 = 0.5*(cf_dcct_1a+cf_dcct_1b);
  Double_t cf_dcct_2 = 0.5*(cf_dcct_2a+cf_dcct_2b);  

  // get bunch crossing information
  // -- number of bc
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  // -- bucket info
  Int_t *BucketA = new Int_t [nIBC];
  Int_t *BucketC = new Int_t [nIBC];  
  GetBucketInfo(BucketA,BucketC);
  
  // set up output tree for rates
  TFile *IntensityFile = new TFile(file_name,"recreate");
  TTree *intensity_tree = new TTree("Bunch_Intensity","Bunch_Intensity");
  Double_t *bunch_intensity_1 = new Double_t[nIBC];
  Double_t *bunch_intensity_2 = new Double_t[nIBC];
  for(Int_t i=0;i<nIBC;i++) bunch_intensity_1[i]=bunch_intensity_2[i]=0.0;
  char *txt_tmp = new char[kg_string_size];
  sprintf(txt_tmp,"bunch_intensity_1[%d]/D",nIBC);
  intensity_tree->Branch("bunch_intensity_1",bunch_intensity_1,txt_tmp);
  sprintf(txt_tmp,"bunch_intensity_2[%d]/D",nIBC);
  intensity_tree->Branch("bunch_intensity_2",bunch_intensity_2,txt_tmp);
  intensity_tree->Branch("cf_dcct_1a",&cf_dcct_1a,"cf_dcct_1a/D");
  intensity_tree->Branch("cf_dcct_1b",&cf_dcct_1b,"cf_dcct_1b/D");
  intensity_tree->Branch("cf_dcct_1",&cf_dcct_1,"cf_dcct_1/D");      
  intensity_tree->Branch("cf_dcct_2a",&cf_dcct_2a,"cf_dcct_2a/D");
  intensity_tree->Branch("cf_dcct_2b",&cf_dcct_2b,"cf_dcct_2b/D");
  intensity_tree->Branch("cf_dcct_2",&cf_dcct_2,"cf_dcct_2/D");      

  // get the intensities
  //average over n_half_range before and after the chosen index
  counter = 0.0;
  Int_t m = start;
  while(m<=end) {
    g_vdm_Tree->GetEntry(m);
    counter += 1.0;
    for(Int_t j=0;j<nIBC;j++) {
      Int_t idx1 = ((Int_t) (BucketA[j]/10.0));
      Int_t idx2 = ((Int_t) (BucketC[j]/10.0));
      bunch_intensity_1[j] += (bunch1[idx1]/gBeamA);
      bunch_intensity_2[j] += (bunch2[idx2]/gBeamB);
    }
    m++;
  } // end while
  for(Int_t j=0;j<nIBC;j++) { // average
    bunch_intensity_1[j] /= counter;
    bunch_intensity_2[j] /= counter;
  }
  intensity_tree->Fill();

  // save  tree
  IntensityFile->cd();
  intensity_tree->SetDirectory(IntensityFile);
  intensity_tree->Write();
  IntensityFile->Close();

  // free memory
  delete [] bunch1;
  delete [] bunch2;     
  delete [] bunch_intensity_1;
  delete [] bunch_intensity_2;     
  delete [] file_name;
  delete [] txt_tmp;  
  delete [] BucketA;
  delete [] BucketC;  

}

//-------------------------------------------------------
// Create the files with the information on the beam
// intensities for all bunch crossings 
//-------------------------------------------------------

void Create_beam_intensity_file(Int_t Fill, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX

{
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create files for all scans
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    Create_single_beam_intensity_file(i,opt);
  }
}
