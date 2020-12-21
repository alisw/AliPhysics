

#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Compute the separations for beam-beam deflection 
//-------------------------------------------------------

Double_t BBDcorrection(Double_t sep, Double_t d1, Double_t d2)
{
  Double_t Delta = (d1+d2)/1000.0; // factor of 1000 to convert to mm
  Double_t SepNew = sep+Delta;
  if ((TMath::Abs(SepNew)-TMath::Abs(sep))<0) SepNew = sep-Delta;
  return SepNew;
}

//-------------------------------------------------------
// Correct the separations for beam-beam deflection using the
// fit of orbit-drift data
//-------------------------------------------------------

void Get_BBD_separations(Int_t scan, Int_t opt)
// opt: 0 => nominal separations
// opt: 1 =>  separations corrected for orbit drift
{
  // -- get number of separations
  Int_t n_sep_x = FindNumberSeparations(1, scan);
  Int_t n_sep_y = FindNumberSeparations(2, scan);
  if (n_sep_x != n_sep_y) {
    cout << " number of x and y separations is not compatible! " << endl;
    exit(-1);
  }
  Int_t n_sep = n_sep_x;
  Double_t *sep_x = new Double_t[n_sep];
  Double_t *sep_y = new Double_t[n_sep];  
 
  // create names for BBD separation files created by Ivan
  // ** WARNING: Directory structure and nameing convention
  //    valid for fill 7483 ... the corresponding structure
  //    has to be added for other fills ... **
  char *bbd_file_name = new char[kg_string_size];
  if (opt == 0) 
    sprintf(bbd_file_name,"../Fill-%d/BBD_FromIvan/NOM-%d/ROOT/bbroot_%d_V0.root",
	    g_vdm_Fill,(scan+1),g_vdm_Fill);
  else
    sprintf(bbd_file_name,"../Fill-%d/BBD_FromIvan/ODC-%d/ROOT/bbroot_%d_V0.root",
	    g_vdm_Fill,(scan+1),g_vdm_Fill);
   
  // open  file and get the trees
  TFile *bbd_file = new TFile(bbd_file_name);
  TTree *sep_tree = (TTree *) bbd_file->Get("beamsep");
  sep_tree->ResetBranchAddresses();
  Int_t bc;
  Int_t nseps;
  Double_t *sepx = new Double_t [n_sep*2];
  Double_t *sepy = new Double_t [n_sep*2];  
  Double_t *dx1 = new Double_t [n_sep*2];
  Double_t *dy1 = new Double_t [n_sep*2];  
  Double_t *dx2 = new Double_t [n_sep*2];
  Double_t *dy2 = new Double_t [n_sep*2];  
  sep_tree->SetBranchAddress("bc",&bc);
  sep_tree->SetBranchAddress("nseps",&nseps);  
  sep_tree->SetBranchAddress("sepx",sepx);  
  sep_tree->SetBranchAddress("sepy",sepy);  
  sep_tree->SetBranchAddress("dx1",dx1);  
  sep_tree->SetBranchAddress("dy1",dy1);  
  sep_tree->SetBranchAddress("dx2",dx2);  
  sep_tree->SetBranchAddress("dy2",dy2);  

  // now prepare names of output
  char *sep_file_name_x = new char[kg_string_size];
  char *sep_file_name_y = new char[kg_string_size];    
  if (opt == 0) {
    sprintf(sep_file_name_x,"../Fill-%d/NomBBDSep_x_Scan_%d.root",
	    g_vdm_Fill,scan);
    sprintf(sep_file_name_y,"../Fill-%d/NomBBDSep_y_Scan_%d.root",
	    g_vdm_Fill,scan);
  } else {
    sprintf(sep_file_name_x,"../Fill-%d/ODCBBDSep_x_Scan_%d.root",
	    g_vdm_Fill,scan);
    sprintf(sep_file_name_y,"../Fill-%d/ODCBBDSep_y_Scan_%d.root",
	    g_vdm_Fill,scan);
  } 

  // next piece is repeated twice, because we want to use
  // the same name for the tree ...
  // -- create tree with separations
  TFile *BBDFile_x = new TFile(sep_file_name_x,"recreate");
  TTree *bbd_sep_tree_x = new TTree("Separations","Separations");
  char *txt_tmp = new char[kg_string_size];
  sprintf(txt_tmp,"separation[%d]/D",n_sep);
  bbd_sep_tree_x->Branch("separation",sep_x,txt_tmp);

  // loop over bunches
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  Int_t nEntries = sep_tree->GetEntries();
  if (nIBC != nEntries) {
    cout << " number of bc is not compatible! " << endl;
    exit(-1);
  }

  for(Int_t k=0;k<nIBC;k++) {
     sep_tree->GetEntry(k);
     Int_t nseps_half = nseps>>1;
     if (nseps_half != n_sep) {
       cout << " number of separations is not compatible! " << nseps << " - " << n_sep << " - " << nseps_half << endl;
       exit(-1);
     }
  
     for(Int_t isep = 0;isep<nseps_half;isep++) {
       
       sep_x[isep] = BBDcorrection(sepx[isep],dx1[isep],dx2[isep]);
       
       //     cout << " -x- " << isep << ": " << sep_x[isep]-sepx[isep] << endl;
     }
     BBDFile_x->cd();
     bbd_sep_tree_x->Fill();
   }

  // save
  BBDFile_x->cd();
  bbd_sep_tree_x->Write();
  BBDFile_x->Close();
  
  // here is the repetition for the y-coordinate
  // -- create tree with separations
  TFile *BBDFile_y = new TFile(sep_file_name_y,"recreate");
  TTree *bbd_sep_tree_y = new TTree("Separations","Separations");  
  bbd_sep_tree_y->Branch("separation",sep_y,txt_tmp);

  // loop over bunches
  for(Int_t k=0;k<nIBC;k++) {
     sep_tree->GetEntry(k);
     Int_t nseps_half = nseps>>1; 
     for(Int_t isep = 0;isep<nseps_half;isep++) {
       
       sep_y[isep] = BBDcorrection(sepy[isep+nseps_half],dy1[isep+nseps_half],dy2[isep+nseps_half]);
       
       // cout << " -y- " << sep_y[isep]-sepy[isep] << endl;
     }
     BBDFile_y->cd();     
     bbd_sep_tree_y->Fill();     
  }

  // save
  BBDFile_y->cd();     
  bbd_sep_tree_y->Write();     
  BBDFile_y->Close();  
  
  // clean up
  delete [] sep_file_name_x;
  delete [] sep_file_name_y;
  delete [] bbd_file_name;
  delete [] sep_x;
  delete [] sep_y;
  delete [] sepx;
  delete [] sepy;
  delete [] dx1;
  delete [] dx2;
  delete [] dy1;
  delete [] dy2;

}


//-------------------------------------------------------
// Create root files with the information of
// the separations corrected for beam-beam deflection (bbd)
// Note that as input you would need the files produced by Ivan
// You need files for the bbd correction to nominal *and* also
// to orbit-drift-corrected separations 
//-------------------------------------------------------

void Create_BBD_separation_file(Int_t Fill)
{

  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create separation files  corrected for
  // beam-beam deflection
  
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    // nominal separations
    Get_BBD_separations(i,0);
    // separations corrected for orbit drift
    Get_BBD_separations(i,1);
  }

}
