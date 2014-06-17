#include "TTree.h"
#include "TF1.h"
#include "TH2D.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TDatime.h"

#include <iostream>
#include "ProgressBar.h"

#define pidTypeAvailable

Int_t splitPbPbTreeMultiplicity(TString pathTree, Bool_t isMC, Int_t stepSize = 2000, Int_t maxMult = 14000,
                                TString fileNameTree = "bhess_PIDetaTree.root", TString treeName = "fTree") 
{ 
	// Extract the data tree		       
	TFile* fTree = TFile::Open(Form("%s/%s", pathTree.Data(), fileNameTree.Data()));
  if (!fTree)  {
    std::cout << "Failed to open file \"" << Form("%s/%s", pathTree.Data(), fileNameTree.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  TTree* tree = dynamic_cast<TTree*>(fTree->Get(treeName.Data()));
  if (!tree) {
    std::cout << "Failed to load data tree!" << std::endl;
    return -1;
  }
  
  Long64_t nTreeEntries = tree->GetEntriesFast();
  Long64_t nTreeEntriesRejected = 0;
  Long64_t nTreeEntriesSkipped = 0;
  
  Double_t pTPC = 0.; 
  Double_t dEdx = 0.; 
  Double_t dEdxExpected = 0.;
  Double_t tanTheta = 0.; 
  Int_t fMultiplicity = 0.; 
  UShort_t tpcSignalN = 0.;
#ifdef pidTypeAvailable
  UChar_t pidType = 0;
#endif
  
  tree->SetBranchStatus("*", 0); // Disable all branches
  tree->SetBranchStatus("pTPC", 1);
  tree->SetBranchStatus("dEdx", 1);
  tree->SetBranchStatus("dEdxExpected", 1);
  tree->SetBranchStatus("tanTheta", 1);
  tree->SetBranchStatus("fMultiplicity", 1);
  tree->SetBranchStatus("tpcSignalN", 1);
#ifdef pidTypeAvailable
  tree->SetBranchStatus("pidType", 1);
#endif
  
  tree->SetBranchAddress("pTPC", &pTPC);
  tree->SetBranchAddress("dEdx", &dEdx);
  tree->SetBranchAddress("dEdxExpected", &dEdxExpected);
  tree->SetBranchAddress("tanTheta", &tanTheta);
  tree->SetBranchAddress("tpcSignalN", &tpcSignalN);
  tree->SetBranchAddress("fMultiplicity", &fMultiplicity);
#ifdef pidTypeAvailable
  tree->SetBranchAddress("pidType", &pidType);
#endif
 
  /* OLD - For nContributorsToPrimaryVertex
  const Int_t stepSize = 600;
  
  const Int_t numMultBins = TMath::Ceil((4000. - 0.) / stepSize);
  */
  
  const Int_t numMultBins = TMath::Ceil((maxMult - 0.) / stepSize);
  
  TFile* fSave[numMultBins];
  TTree* treeMult[numMultBins];
  
  for (Int_t i = 0; i < numMultBins; i++)   {
    // Output files
    TString savefileName = Form("%s_mult%d_%d.root", fileNameTree.ReplaceAll(".root", "").Data(), i * stepSize, (i + 1) * stepSize - 1);
    
    fSave[i] = TFile::Open(Form("%s/%s", pathTree.Data(), savefileName.Data()), "recreate");
    if (!fSave[i]) {
      std::cout << "Failed to open save file \"" << Form("%s/%s", pathTree.Data(), savefileName.Data()) << "\"!" << std::endl;
      return -1;
    }
    
    fSave[i]->cd();
    treeMult[i] = new TTree("fTree", Form("Tree with multiplicity %d - %d", i * stepSize, (i + 1) * stepSize - 1));
    treeMult[i]->Write(0, TObject::kOverwrite);
    treeMult[i]->Branch("pTPC", &pTPC);
    treeMult[i]->Branch("dEdx", &dEdx);
    treeMult[i]->Branch("dEdxExpected", &dEdxExpected);
    treeMult[i]->Branch("tanTheta", &tanTheta);
    treeMult[i]->Branch("tpcSignalN", &tpcSignalN);
    treeMult[i]->Branch("fMultiplicity", &fMultiplicity);
#ifdef pidTypeAvailable
    treeMult[i]->Branch("pidType", &pidType);
#endif
  }
  
  TF1 cutFunc("cutFunc","50./TMath::Power(x,1.3)", 0.05, 6);
  
  Bool_t isProtonTree = kFALSE;
  if (treeName.CompareTo("fTree") == 0) {
    isProtonTree = kTRUE;
    printf("Applying cut for protons (if no MC)...\n");
  }
  else {
    printf("Species seems to be different from protons - no cut...\n");
  }
  
  progressbar(0.);
  
  for (Long64_t i = 0; i < nTreeEntries; i++) {
    tree->GetEntry(i);

    // Filter tracks according to shape recognition
    if (dEdx < 10 || dEdxExpected < 10 || (!isMC && isProtonTree && pTPC < 0.6 && dEdx < cutFunc.Eval(pTPC))) {
      nTreeEntriesRejected++;
      continue;
    }
    
    Int_t multBin = (Int_t)(fMultiplicity / stepSize);
    if (multBin >= numMultBins || multBin < 0) {
      nTreeEntriesSkipped++;
      //std::cout << std::endl << "Warning: Skipping entry with fMultiplicity " << fMultiplicity << std::endl;
      //printedSomething = kTRUE;
      continue;
    }
    
    treeMult[multBin]->Fill();
    
    if (i % 100000 == 0)
      progressbar(100. * (((Double_t)i) / nTreeEntries));
  }
  
  Long64_t nSplitTreeEntries = 0;
  for (Int_t i = 0; i < numMultBins; i++)   {
    fSave[i]->cd();
    treeMult[i]->Write(0, TObject::kOverwrite);
    nSplitTreeEntries += treeMult[i]->GetEntriesFast();
    fSave[i]->Close();
  }
  
  progressbar(100.);
  
  std::cout << std::endl << "Entries original tree: " << nTreeEntries << std::endl << "Entries split trees: " << nSplitTreeEntries << std::endl;
  std::cout << "Rejected entries (cuts): " << nTreeEntriesRejected << std::endl;
  std::cout << "Skipped entries (out of range): " << nTreeEntriesSkipped << std::endl;
  std::cout << "=> Difference (after subtracting rejected and skipped entries) is "
            << nTreeEntries - nTreeEntriesRejected - nTreeEntriesSkipped - nSplitTreeEntries << std::endl;
  
  return 0;
}
