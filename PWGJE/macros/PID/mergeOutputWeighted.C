#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TKey.h"
#include "TList.h"
#include "TProfile.h"
#include "TString.h"
#include "TSystem.h"
#include "AliCFContainer.h"

void mergeOutputWeighted(const Int_t numDirs, const TString* inputDirs, const TString fileWithWeighting,
                         const TString dirInFileWithWeighting, const TString listInFileWithWeighting, const TString fileToMerge,
                         const TString dirInFileToMerge, 
                         const TString listInFileToMerge/*leave empty if objects directly in file w/o list*/,
                         const TString outDir, const Int_t normScheme)
{
  // Based on PWG/Tools/AliAnalysisHelperJetTask.C:
  // This is used to merge the analysis-output from different 
  // data samples/pt_hard bins
  // in case the eventweigth was set to xsection/ntrials already, this
  // is not needed. Both methods only work in case we do not mix different 
  // pt_hard bins, and do not have overlapping bins

  
  Float_t xsection[numDirs];
  Float_t nTrials[numDirs];
  Float_t sf[numDirs];
  Float_t sfSum = 0;
  Float_t sfSmallest = -1;
  TList* lIn[numDirs];
  TDirectory* dIn[numDirs];
  TFile* fIn[numDirs];
  Bool_t isKeyList = kFALSE;
  Bool_t hasSubDir = kFALSE;
  
  TFile* fTemp = 0x0;
  TDirectory* dTemp = 0x0;
  TList* lTemp = 0x0;
  
  for (Int_t i = 0; i < numDirs; i++) {
    xsection[i] = 0;
    nTrials[i] = 0;
    sf[i] = 0;
    dIn[i] = 0x0;
    lIn[i] = 0x0;
    fIn[i] = 0x0;
  }
  
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  
  
  // Obtain nTrials and xsection histos for all pt_hard bins first
  // from file containing this information and also obtain the lists
  // and objects from the files to merge
  

  Int_t ibTotal = 0;
  
  for (Int_t iDir = 0; iDir < numDirs; iDir++) {
    TString cDir = Form("%s/%s", inputDirs[iDir].Data(), fileWithWeighting.Data());
    fTemp = TFile::Open(cDir.Data(), "READ");
    
    if (!fTemp) {
      Printf("%s:%d File %s not found, exiting...", __FILE__, __LINE__, cDir.Data());
      
      return;
    }
    
    TString dirInFileWithWeightingUsed = dirInFileWithWeighting;
    TString listInFileWithWeightingUsed = listInFileWithWeighting;
    
    if (strlen(dirInFileWithWeighting.Data()) == 0) {
      dTemp = gDirectory;
    }
    else {
      if (dirInFileWithWeighting == "default") {
        dTemp = (TDirectory*)fTemp->Get(fTemp->GetListOfKeys()->At(0)->GetName());
        if (dTemp)
          dirInFileWithWeightingUsed = Form("%s", dTemp->GetName());
      }
      else
        dTemp = (TDirectory*)fTemp->Get(dirInFileWithWeighting.Data());
    }
    
    if (!dTemp) {
      Printf("%s:%d No directory %s found, exiting...", __FILE__, __LINE__, dirInFileWithWeightingUsed.Data());
      fTemp->ls();
      
      return;
    }

    if (listInFileWithWeighting == "default") {
      lTemp = (TList*)dTemp->Get(dTemp->GetListOfKeys()->At(0)->GetName());
      if (lTemp)
        listInFileWithWeightingUsed = Form("%s", lTemp->GetName());
    }
    else
      lTemp = (TList*)dTemp->Get(listInFileWithWeighting.Data());
    
    if (!lTemp) {
      Printf("%s:%d No list %s found, exiting...", __FILE__, __LINE__, listInFileWithWeightingUsed.Data());
      fTemp->ls();
      
      return;
    }
    
    TH1* hTrials = (TH1F*)lTemp->FindObject("fh1Trials");
    if (!hTrials) {
      Printf("%s:%d fh1PtHard_Trials not found in list, exiting...", __FILE__, __LINE__);
      
      return;
    }
    
    TProfile* hXsec = (TProfile*)lTemp->FindObject("fh1Xsec");
    if (!hXsec) {
      Printf("%s:%d fh1Xsec  not found in list, exiting...", __FILE__, __LINE__);
      
      return;
    }
    
    xsection[ibTotal] = hXsec->GetBinContent(1);
    nTrials[ibTotal] = hTrials->Integral();
    sf[ibTotal] = xsection[ibTotal]/ nTrials[ibTotal];
    sfSum += sf[ibTotal];
    
    if (sf[ibTotal] < sfSmallest || sfSmallest < 0)
      sfSmallest = sf[ibTotal];
    
    fTemp->Close();
    
    // Obtain lists for merging
    TString fileNameToMerge = Form("%s/%s", inputDirs[iDir].Data(), fileToMerge.Data());
    fIn[ibTotal] = TFile::Open(fileNameToMerge.Data(), "READ");
    
    if (!fIn[ibTotal]) {
      Printf("%s:%d File %s not found, exiting...", __FILE__, __LINE__, fileNameToMerge.Data());
      
      return;
    }
    
    TString dirInFileToMergeUsed = dirInFileToMerge;
    TString listInFileToMergeUsed = listInFileToMerge;
    
    if (strlen(dirInFileToMerge.Data()) == 0) {
      dIn[ibTotal] = gDirectory;
      if (ibTotal == 0)
        hasSubDir = kFALSE;
    }
    else {
      if (dirInFileToMerge == "default") {
        dIn[ibTotal] = (TDirectory*)fIn[ibTotal]->Get(fIn[ibTotal]->GetListOfKeys()->At(0)->GetName());
        if (dIn[ibTotal])
          dirInFileToMergeUsed = Form("%s", dIn[ibTotal]->GetName());
      }
      else {
        dIn[ibTotal] = (TDirectory*)fIn[ibTotal]->Get(dirInFileToMerge.Data());
      }
      
      if (ibTotal == 0)
        hasSubDir = kTRUE;
    }
    
    if (!dIn[ibTotal]) {
      Printf("%s:%d No directory %s found, exiting...", __FILE__, __LINE__, dirInFileToMergeUsed.Data());
      fIn[ibTotal]->ls();
      
      return;
    }
    
    if (strlen(listInFileToMerge.Data()) == 0) {
      lIn[ibTotal] = dIn[ibTotal]->GetListOfKeys();
      if (ibTotal == 0)
        isKeyList = kTRUE;
    }
    else {
      if (listInFileToMerge == "default") {
        lIn[ibTotal] = (TList*)dIn[ibTotal]->Get(dIn[ibTotal]->GetListOfKeys()->At(0)->GetName());
        if (lIn[ibTotal])
          listInFileToMergeUsed = Form("%s", lIn[ibTotal]->GetName());
      }
      else
        lIn[ibTotal] = (TList*)dIn[ibTotal]->Get(listInFileToMerge.Data());
      
      if (ibTotal == 0)
        isKeyList = kFALSE;
    }
    
    Printf("Merging file %s %s\nWeight factor %e\n", fileNameToMerge.Data(), dirInFileToMergeUsed.Data(), sf[ibTotal]);
    if (!lIn[ibTotal]) {
      Printf("%s:%d No list %s found, exiting...", __FILE__, __LINE__, listInFileToMergeUsed.Data());
      lIn[ibTotal]->ls();
      
      return;
    }
    
    ibTotal++;
  }

  if (ibTotal == 0) {
    Printf("%s:%d No files found for merging, exiting", __FILE__, __LINE__);
    
    return;
  }
  if (ibTotal < numDirs) {
    Printf("%s:%d No files found for merging in at least one directory, exiting", __FILE__, __LINE__);
    
    return;
  }

  if (sfSum <= 0) {
    Printf("%s:%d Sum of scale factors <= 0, exiting", __FILE__, __LINE__);
    
    return;
  }
  
  if (normScheme > 0) {
    printf("\nNormalised scale factors:\n");
    for (Int_t ib = 0; ib < ibTotal; ib++) {
      if (normScheme == 1)
        sf[ib] = sf[ib] / sfSum;
      else
        sf[ib] = sf[ib] / sfSmallest;
      printf("Index %d:\t%f\n", ib, sf[ib]);
    }
    printf("\n");
  }
  
  // Open output file
  TFile* fOut = 0;
  gSystem->mkdir(outDir.Data(), kTRUE);
  fOut = TFile::Open(Form("%s/%s", outDir.Data(), fileToMerge.Data()), "RECREATE");
  
  if (!fOut) {
    Printf("%s:%d Failed to open output file %s, exiting", __FILE__, __LINE__, Form("%s/%s", outDir.Data(), fileToMerge.Data()));
    
    return;
  }
  fOut->cd();
  TDirectory* dOut = hasSubDir ? fOut->mkdir(dIn[0]->GetName()) : gDirectory;
  if (!dOut) {
    Printf("%s:%d Failed to create directory %s in output file, exiting", __FILE__, __LINE__, dIn[0]->GetName());
    
    return;
  }
  dOut->cd();
  TList* lOut = new TList();
  lOut->SetName(lIn[0]->GetName());
  // Merge lists / objects
  
  for(Int_t ie = 0; ie < lIn[0]->GetEntries(); ++ie) {
    TH1* h1Add = 0x0;
    THnSparse* hnAdd = 0x0;
    AliCFContainer* contAdd = 0x0;
    
    for(Int_t ib = 0; ib < ibTotal; ++ib) {
      TObject* object = 0x0;
      
      if (isKeyList) {
        TKey* lKey = (TKey*)lIn[ib]->At(ie);
        object = lKey->ReadObj();
      }
      else {
        object = lIn[ib]->At(ie);
      }
      
      fOut->cd();
      if (object->InheritsFrom("TH1")) {
        TH1* h1 = (TH1*)object;
        if (ib == 0) {
          h1Add = (TH1*)h1->Clone(h1->GetName());
          if (h1Add->GetSumw2N() == 0) {
            // If no errors have been set so far, calculate errors as sqrt(binContent) and set sumw2 before
            h1Add->Sumw2();

            // Also catch under and overflow bins!
            // -> Compare TH1::Add(...)
            Int_t nBinsX = h1Add->GetNbinsX();
            Int_t nBinsY = h1Add->GetNbinsY();
            Int_t nBinsZ = h1Add->GetNbinsZ();
            
            if (h1Add->GetDimension() < 2)
              nBinsY = -1;
            if (h1Add->GetDimension() < 3)
              nBinsZ = -1;
            
            for (Int_t xBin = 0; xBin <= nBinsX + 1; xBin++) {
              for (Int_t yBin = 0; yBin <= nBinsY + 1; yBin++) {
                for (Int_t zBin = 0; zBin <= nBinsZ + 1; zBin++) {
                  Double_t binContent = h1Add->GetBinContent(xBin, yBin, zBin);
                  h1Add->SetBinError(xBin, yBin, zBin, TMath::Sqrt(binContent));
                }
              }
            }
          }
          
          h1Add->Scale(sf[ib]);
        }
        else {
          h1Add->Add(h1, sf[ib]);
        }
      }
      else if (object->InheritsFrom("THnSparse")) {
        THnSparse* hn = (THnSparse*)object;
        if (ib == 0) {
          hnAdd = (THnSparse*)hn->Clone(hn->GetName());
          
          if (!hnAdd->GetCalculateErrors()) {
            // If no errors have been set so far, calculate errors as sqrt(binContent) and set sumw2 before
            hnAdd->Sumw2();
            
            for (Long64_t bin = 0; bin < hnAdd->GetNbins(); bin++) {
              Double_t binContent = hnAdd->GetBinContent(bin);
              hnAdd->SetBinError(bin, TMath::Sqrt(binContent));
            }  
          }
  
          hnAdd->Scale(sf[ib]);
        }
        else {
          hnAdd->Add(hn, sf[ib]);
        }
      }
      else if (object->InheritsFrom("AliCFContainer")) {
        AliCFContainer* cont = (AliCFContainer*)object;
        if (ib == 0) {
          contAdd = (AliCFContainer*)cont->Clone(cont->GetName());
          contAdd->Scale(sf[ib]);
        }
        else {
          contAdd->Add(cont, sf[ib]);
        }
      }
      
      printf("\r");
      printf("Done file %d/%d for item %d/%d                    ", ib + 1, ibTotal, ie + 1, lIn[0]->GetEntries());
      fflush(stdout);
    }// ib
    
    if (h1Add)
      lOut->Add(h1Add);
    else if (hnAdd)
      lOut->Add(hnAdd);
    else if (contAdd)
      lOut->Add(contAdd);
    
    //printf("\n\nDone item %d/%d\n\n", ie + 1, lIn[0]->GetEntries()); 
  }//ie
  dOut->cd();
  lOut->Write(isKeyList ? 0 : lOut->GetName(), isKeyList ? 0 : TObject::kSingleKey);
  fOut->Close();
  
  for (Int_t i = 0; i < numDirs; i++) {
    if (fIn[i])
      fIn[i]->Close();
  }
  
  delete lOut;
  lOut = 0x0;
  
  printf("\n\n");
  
  printf("WARNING: If AnalysisResults.root is taken to obtain the weight for merging and this is done for the jet analysis with pile-up rejection, the weights in AnalysisResults.root do NOT take into account the pile-up rejection, whereas this is true for the values in the PID tasks. Does not matter for current MC productions, which do not simulate pile-up (und false positives are negligable).\n\n");
}

void mergeOutputWeighted(const TString inputDirBase, const TString subDir, const TString dirInFileWithWeighting,
                         const TString listInFileWithWeighting, const TString fileToMerge, const TString dirInFileToMerge,
                         const TString listInFileToMerge/*leave empty if objects directly in file w/o list*/,
                         const Int_t normScheme)
{
  const Int_t numDirs = 10;
  TString inputDirs[numDirs];
  
  TString subDirExtended = "";
  if (subDir != "")
    subDirExtended = Form("/%s", subDir.Data());
  
  for (Int_t i = 0; i < numDirs; i++) 
    inputDirs[i] = Form("%s/LHC11a1%c%s", inputDirBase.Data(), 'a' + i, subDirExtended.Data());
  
  mergeOutputWeighted(numDirs, inputDirs, "AnalysisResults.root", dirInFileWithWeighting, listInFileWithWeighting,
                      fileToMerge, dirInFileToMerge, listInFileToMerge, Form("%s/allJetPtMergedWeighted", inputDirBase.Data()),
                      normScheme);
}


void mergeOutputWeighted(const Int_t mergeUpToHardBin, const TString inputDirBase, const TString outDirBase, const TString subDir,
                         const TString dirInFileWithWeighting,
                         const TString listInFileWithWeighting, const TString fileToMerge, const TString dirInFileToMerge,
                         const TString listInFileToMerge/*leave empty if objects directly in file w/o list*/,
                         const Int_t normScheme,
                         const TString period = "LHC11a1")
{
  const TString hardBin[10] = { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j" };
  
  
  const Int_t numDirs = mergeUpToHardBin;
  TString inputDirs[numDirs];
  
  TString outDir = Form("%s/%s%s", outDirBase.Data(), period.Data(), hardBin[mergeUpToHardBin - 1].Data());
  
  TString subDirExtended = "";
  if (subDir != "")
    subDirExtended = Form("/%s", subDir.Data());
  
  for (Int_t i = 0; i < numDirs; i++) 
    inputDirs[i] = Form("%s/%s%c%s", inputDirBase.Data(), period.Data(), 'a' + i, subDirExtended.Data());
  
  mergeOutputWeighted(numDirs, inputDirs, "AnalysisResults.root", dirInFileWithWeighting, listInFileWithWeighting,
                      fileToMerge, dirInFileToMerge, listInFileToMerge, outDir.Data(),
                      normScheme);
}



// a 'mergeOutputWeighted.C+("finalCuts/MC_pp/7TeV/LHC11a1/corrected/legoTrain/nclCut", "finalCuts/MC_pp/7TeV/LHC11a1/corrected/legoTrain/nclCut/mergedWeightedUntilHardBin", "", "default", "default", "AnalysisResults.root", "default", "default", 0, kTRUE, "LHC11a1")'
// a 'mergeOutputWeighted.C+("finalCuts/MC_pp/7TeV/LHC11a1/corrected/legoTrain/nclCut", "finalCuts/MC_pp/7TeV/LHC11a1/corrected/legoTrain/nclCut/mergedWeightedUntilHardBin", "", "default", "default", "PWGJE_taskPID_Jets_efficiency.root", "", "", 0, kTRUE, "LHC11a1")'
// a 'mergeOutputWeighted.C+("finalCuts/MC_pp/7TeV/LHC11a1/corrected/legoTrain/nclCut", "finalCuts/MC_pp/7TeV/LHC11a1/corrected/legoTrain/nclCut/mergedWeightedUntilHardBin", "", "default", "default", "PWGJE_taskPID_Jets.root", "", "default", 0, kTRUE, "LHC11a1")'

// New MC production
// a 'mergeOutputWeighted.C+("finalCuts/MC_pp/7TeV/LHC14b6/legoTrain/nclCut", "finalCuts/MC_pp/7TeV/LHC14b6/legoTrain/nclCut/allJetPtMergedWeighted", "", "PWGJE_IDFragmentationFunction_clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00_AODMC2b_AODMCb_cl0_trackRefs", "default", "AnalysisResults.root", "PWGJE_IDFragmentationFunction_clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00_AODMC2b_AODMCb_cl0_trackRefs", "default", 0, kFALSE, "LHC14b6")'
// a 'mergeOutputWeighted.C+("finalCuts/MC_pp/7TeV/LHC14b6/legoTrain/nclCut", "finalCuts/MC_pp/7TeV/LHC14b6/legoTrain/nclCut/allJetPtMergedWeighted", "", "PWGJE_IDFragmentationFunction_clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00_AODMC2b_AODMCb_cl0_trackRefs", "default", "PWGJE_taskPID_Jets_efficiency.root", "", "", 0, kFALSE, "LHC14b6")'
// a 'mergeOutputWeighted.C+("finalCuts/MC_pp/7TeV/LHC14b6/legoTrain/nclCut", "finalCuts/MC_pp/7TeV/LHC14b6/legoTrain/nclCut/allJetPtMergedWeighted", "", "PWGJE_IDFragmentationFunction_clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00_AODMC2b_AODMCb_cl0_trackRefs", "default", "PWGJE_taskPID_Jets.root", "", "default", 0, kFALSE, "LHC14b6")'
void mergeOutputWeighted(const TString inputDirBase, const TString outDirBase, const TString subDir,
                         const TString dirInFileWithWeighting,
                         const TString listInFileWithWeighting, const TString fileToMerge, const TString dirInFileToMerge,
                         const TString listInFileToMerge/*leave empty if objects directly in file w/o list*/,
                         const Int_t normScheme,
                         const Bool_t mergeUntilHardBin,
                         const TString period)
{
  Int_t upperIndex = 10;
  if (period == "LHC14b6")
    upperIndex = 4; // Stops with hard bin "d"
  
  for (Int_t i = 1; i <= upperIndex; i++) {
    if (!mergeUntilHardBin && i < upperIndex)
      continue;
    mergeOutputWeighted(i, inputDirBase, outDirBase, subDir, dirInFileWithWeighting, listInFileWithWeighting, fileToMerge, dirInFileToMerge,
                        listInFileToMerge, normScheme, period);
  }
}

  