/******************************************************************************************************************************************
Contact person: Marco Leoncino (leoncino@to.infn.it)
Macro to run the ZDC QA trending by accessing the std QA output, to be mainly used with the automatic scripts to fill the QA repository.
Launch with aliroot -l -b -q "MakeTrendingZDCQA.C(\"${fullpath}/QAresults.root\", ${run}, ...) 
The macro produces a file containing the tree of trending variables and the main plots.
A feature that displays the plots in canvases must be enable when needed.
******************************************************************************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TClassTable.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLine.h>
#include <TGrid.h>
#include <TBits.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFileMerger.h>
#include <TGridResult.h>
#include <TSystem.h>
#include <TGaxis.h>

#endif

Int_t MakeTrendingZDCQA(TString qafilename,               //full path of the QA output; set IsOnGrid to prepend "alien://"
                        Int_t runNumber,                  //run number
                        Bool_t isMC = kFALSE,             //MC flag, to disable meaningless checks
                        Bool_t canvasE = kFALSE,          //enable display plots on canvas and save png
                        Bool_t IsOnGrid = kFALSE,         //set to kTRUE to access files on the grid
                        TString ocdbStorage = "raw://")   //set the default ocdb storage
{   

  // macro to generate tree with ZDC QA trending variables
  // access qa PWGPP output files  
  if (!qafilename) {
    printf("Error - Invalid input file");
    return 1;
  }

  /*set graphic style*/
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetStatColor(kWhite); 
  gStyle->SetStatBorderSize(1);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(10);

  char defaultQAoutput[30]="QAresults.root";
  char * treePostFileName="trending.root";  
  
  if (IsOnGrid) TGrid::Connect("alien://");
  TFile * fin = TFile::Open(qafilename,"r");
  if (!fin) {
    Printf("ERROR: QA output not found. Exiting...\n");
    return -1;
  } else {
    Printf("INFO: QA output file %s open. \n",fin->GetName());
  }
  
  //access histograms lists
  char zdcQAdirName[20]="ZDC_Performance";
  char genListName[20]="QAZDCHists";
  
  TDirectoryFile * zdcQAdir=(TDirectoryFile*)fin->Get(zdcQAdirName);
  if (!zdcQAdir) {
    Printf("ERROR: ZDC QA directory not present in input file.\n");
    return -1;
  }
  
  TList * generalList=(TList*)zdcQAdir->Get(genListName);

  if (!generalList) Printf("WARNING: general QA histograms absent or not accessible\n");

  TH1F    *fhTDCZNC          = (TH1F*)generalList->FindObject("fhTDCZNC");       //! TDC ZNC sum
  TH1F    *fhTDCZNA          = (TH1F*)generalList->FindObject("fhTDCZNA");       //! TDC DIFF sum
  TH1F    *fhTDCZNSum        = (TH1F*)generalList->FindObject("fhTDCZNSum");     //! TDC ZNC sum
  TH1F    *fhTDCZNDiff       = (TH1F*)generalList->FindObject("fhTDCZNDiff");    //! TDC DIFF sum
  TH1F    *fhZNCSumQ         = (TH1F*)generalList->FindObject("fhZNCSumQ");      //! ZNC sum 4Q
  TH1F    *fhZNASumQ         = (TH1F*)generalList->FindObject("fhZNASumQ");      //! ZNA sum 4Q
  TH1F    *fhZPCSumQ         = (TH1F*)generalList->FindObject("fhZPCSumQ");      //! ZPC sum 4Q
  TH1F    *fhZPASumQ         = (TH1F*)generalList->FindObject("fhZPASumQ");      //! ZPA sum 4Q
  TH1F    *fhZEM1Spectrum    = (TH1F*)generalList->FindObject("fhZEM1Spectrum"); //! ZEM1 spectra
  TH1F    *fhZEM2Spectrum    = (TH1F*)generalList->FindObject("fhZEM2Spectrum"); //! ZEM2 spectra
  TH1F    *fhZNCpmc          = (TH1F*)generalList->FindObject("fhZNCpmc");       //! ZNC PMCs
  TH1F    *fhZNApmc          = (TH1F*)generalList->FindObject("fhZNApmc");       //! ZNA PMCs
  TH1F    *fhZPCpmc          = (TH1F*)generalList->FindObject("fhZPCpmc");       //! ZPC PMCs
  TH1F    *fhZPApmc          = (TH1F*)generalList->FindObject("fhZPApmc");       //! ZPA PMCs
  TH2F    *fhZNCCentroid     = (TH2F*)generalList->FindObject("fhZNCCentroid");  //! ZNC centroid
  TH2F    *fhZNACentroid     = (TH2F*)generalList->FindObject("fhZNACentroid");  //! ZNA centroid
  TH1F    *fhPMCZNCemd       = (TH1F*)generalList->FindObject("fhPMCZNCemd");    //! ZNC PMC low gain chain
  TH1F    *fhPMCZNAemd       = (TH1F*)generalList->FindObject("fhPMCZNAemd");    //! ZNA PMC low gain chain
  TH2F    *fDebunch          = (TH2F*)generalList->FindObject("fDebunch");       //! TDC sum vs. diff
  TH1F    *fhTDCZNAcorr      = (TH1F*)generalList->FindObject("fhTDCZNAcorr");   //! ZNA corrected TDC
  TH1F    *fhTDCZNCcorr      = (TH1F*)generalList->FindObject("fhTDCZNCcorr");   //! ZNC corrected TDC
  
  Double_t ZNC_mean = fhZNCpmc->GetMean()/TMath::Sqrt(fhZNCpmc->GetEntries());
  Double_t ZNA_mean = fhZNApmc->GetMean()/TMath::Sqrt(fhZNApmc->GetEntries());
  Double_t ZPC_mean = fhZPCpmc->GetMean()/TMath::Sqrt(fhZPCpmc->GetEntries());
  Double_t ZPA_mean = fhZPApmc->GetMean()/TMath::Sqrt(fhZPApmc->GetEntries());  
  Double_t ZEM1_mean = fhZEM1Spectrum->GetMean()/TMath::Sqrt(fhZEM1Spectrum->GetEntries());
  Double_t ZEM2_mean = fhZEM2Spectrum->GetMean()/TMath::Sqrt(fhZEM2Spectrum->GetEntries());  
  Double_t ZNC_XCent = fhZNCCentroid->GetMean(1);
  Double_t ZNC_YCent = fhZNCCentroid->GetMean(2);    
  Double_t ZNA_XCent = fhZNACentroid->GetMean(1);
  Double_t ZNA_YCent = fhZNACentroid->GetMean(2);
  Double_t ZNC_XCent_err = fhZNCCentroid->GetRMS(1);
  Double_t ZNC_YCent_err = fhZNCCentroid->GetRMS(2);    
  Double_t ZNA_XCent_err = fhZNACentroid->GetRMS(1);
  Double_t ZNA_YCent_err = fhZNACentroid->GetRMS(2);     
  Double_t ZN_TDC_Sum = fhTDCZNSum->GetMean();
  Double_t ZN_TDC_Diff = fhTDCZNDiff->GetMean();  
  Double_t ZN_TDC_Sum_err = fhTDCZNSum->GetRMS();
  Double_t ZN_TDC_Diff_err = fhTDCZNDiff->GetRMS();    
  
  TTree * ttree=new TTree("trending","tree of trending variables");
  ttree->Branch("run",&runNumber,"run/I");
  ttree->Branch("ZNC_mean_value",&ZNC_mean,"ZNC_mean_value/D");
  ttree->Branch("ZNA_mean_value",&ZNA_mean,"ZNA_mean_value/D"); 
  ttree->Branch("ZPC_mean_value",&ZPC_mean,"ZPC_mean_value/D");
  ttree->Branch("ZPA_mean_value",&ZPA_mean,"ZPA_mean_value/D");
  ttree->Branch("ZEM1_mean_value",&ZEM1_mean,"ZEM1_mean_value/D");
  ttree->Branch("ZEM2_mean_value",&ZEM2_mean,"ZEM2_mean_value/D");  
  ttree->Branch("ZNC_X_Centroid",&ZNC_XCent,"ZNC_X_Centroid/D");
  ttree->Branch("ZNC_Y_Centroid",&ZNC_YCent,"ZNC_Y_Centroid/D");   
  ttree->Branch("ZNA_X_Centroid",&ZNA_XCent,"ZNA_X_Centroid/D");
  ttree->Branch("ZNA_Y_Centroid",&ZNA_YCent,"ZNA_Y_Centroid/D");  
  ttree->Branch("ZNC_X_Centroid_Err",&ZNC_XCent_err,"ZNC_X_Centroid_Err/D");
  ttree->Branch("ZNC_Y_Centroid_Err",&ZNC_YCent_err,"ZNC_Y_Centroid_Err/D");   
  ttree->Branch("ZNA_X_Centroid_Err",&ZNA_XCent_err,"ZNA_X_Centroid_Err/D");
  ttree->Branch("ZNA_Y_Centroid_Err",&ZNA_YCent_err,"ZNA_Y_Centroid_Err/D");    
  ttree->Branch("ZN_TDC_Sum",&ZN_TDC_Sum,"ZN_TDC_Sum/D");
  ttree->Branch("ZN_TDC_Diff",&ZN_TDC_Diff,"ZN_TDC_Diff/D");    
  ttree->Branch("ZN_TDC_Sum_Err",&ZN_TDC_Sum_err,"ZN_TDC_Sum_Err/D");
  ttree->Branch("ZN_TDC_Diff_Err",&ZN_TDC_Diff_err,"ZN_TDC_Diff_Err/D");  
    
  Printf(":::: Getting post-analysis info for run %i",runNumber);
  TFile * trendFile = new TFile(treePostFileName,"recreate");

  printf("==============  Saving histograms for run %i ===============\n",runNumber);
  
  fhTDCZNC->Write();      
  fhTDCZNA->Write();      
  fhTDCZNSum->Write();    
  fhTDCZNDiff->Write();   
  fhZNCSumQ->Write();     
  fhZNASumQ->Write();     
  fhZPCSumQ->Write();     
  fhZPASumQ->Write();     
  fhZEM1Spectrum->Write();
  fhZEM2Spectrum->Write();
  fhZNCpmc->Write();      
  fhZNApmc->Write();      
  fhZPCpmc->Write();      
  fhZPApmc->Write();      
  fhZNCCentroid->Write(); 
  fhZNACentroid->Write(); 
  fhPMCZNCemd->Write();   
  fhPMCZNAemd->Write();   
  fDebunch->Write();      
  fhTDCZNAcorr->Write();  
  fhTDCZNCcorr->Write();
    
  ttree->Fill();
  printf("==============  Saving trending quantities in tree for run %i ===============\n",runNumber);
  
  trendFile->cd();
  ttree->Write();
  trendFile->Close();

}