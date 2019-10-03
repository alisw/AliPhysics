// $Id$
//
// Jet embedding from PYTHIA task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetEmbeddingFromPYTHIATask.h"

#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TRandom.h>
#include <TParameter.h>
#include <TH1I.h>
#include <TGrid.h>
#include <THashTable.h>
#include <TSystem.h>

#include "AliVEvent.h"
#include "AliLog.h"

ClassImp(AliJetEmbeddingFromPYTHIATask)

//________________________________________________________________________
AliJetEmbeddingFromPYTHIATask::AliJetEmbeddingFromPYTHIATask() : 
  AliJetEmbeddingFromAODTask("AliJetEmbeddingFromPYTHIATask"),
  fPYTHIAPath(),
  fPtHardBinScaling(),
  fLHC11hAnchorRun(kTRUE),
  fAnchorRun(-1),
  fFileTable(0),
  fUseAsVetoTable(kTRUE),
  fMinEntriesPerPtHardBin(1),
  fCurrentPtHardBin(-1),
  fPtHardBinParam(0),
  fPtHardBinCount(0),
  fDebugEmbedding(kFALSE),
  fHistPtHardBins(0)
{
  // Default constructor.
  SetSuffix("PYTHIAEmbedding");
  fTotalFiles = 140;
  fRandomAccess = kTRUE;
  SetAODMC(kTRUE);
}

//________________________________________________________________________
AliJetEmbeddingFromPYTHIATask::AliJetEmbeddingFromPYTHIATask(const char *name, Bool_t drawqa) : 
  AliJetEmbeddingFromAODTask(name, drawqa),
  fPYTHIAPath("alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/%04d/AliAOD.root"),
  fPtHardBinScaling(),
  fLHC11hAnchorRun(kTRUE),
  fAnchorRun(-1),
  fFileTable(0),
  fUseAsVetoTable(kTRUE),
  fMinEntriesPerPtHardBin(1),
  fCurrentPtHardBin(-1),
  fPtHardBinParam(0),
  fPtHardBinCount(0),
  fDebugEmbedding(kFALSE),
  fHistPtHardBins(0)
{
  // Standard constructor.
  SetSuffix("PYTHIAEmbedding");
  fTotalFiles = 140;
  fRandomAccess = kTRUE;
  SetAODMC(kTRUE);
}

//________________________________________________________________________
AliJetEmbeddingFromPYTHIATask::~AliJetEmbeddingFromPYTHIATask()
{
  // Destructor
}

//________________________________________________________________________
void AliJetEmbeddingFromPYTHIATask::UserCreateOutputObjects()
{
  if (!fQAhistos)
    return;

  AliJetModelBaseTask::UserCreateOutputObjects();

  fHistPtHardBins = new TH1F("fHistPtHardBins", "fHistPtHardBins", 11, 0, 11);
  fHistPtHardBins->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistPtHardBins->GetYaxis()->SetTitle("total events");
  fOutput->Add(fHistPtHardBins);

  const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  for (Int_t i = 1; i < 12; i++) 
    fHistPtHardBins->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));

  fHistEmbeddingQA = new TH1F("fHistEmbeddingQA", "fHistEmbeddingQA", 2, 0, 2);
  fHistEmbeddingQA->GetXaxis()->SetTitle("Event state");
  fHistEmbeddingQA->GetYaxis()->SetTitle("counts");
  fHistEmbeddingQA->GetXaxis()->SetBinLabel(1, "OK");
  fHistEmbeddingQA->GetXaxis()->SetBinLabel(2, "Not embedded");
  fOutput->Add(fHistEmbeddingQA);

  fHistRejectedEvents = new TH1F("fHistRejectedEvents", "fHistRejectedEvents", 500, -0.5, 499.5);
  fHistRejectedEvents->GetXaxis()->SetTitle("# of rejected events");
  fHistRejectedEvents->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistRejectedEvents);

  PostData(1, fOutput);
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromPYTHIATask::ExecOnce() 
{
  if (fPtHardBinScaling.GetSize() > 0) {
    Double_t sum = 0;
    for (Int_t i = 0; i < fPtHardBinScaling.GetSize(); i++) 
      sum += fPtHardBinScaling[i];
    
    if (sum == 0) {
      AliWarning("No hard pt bin scaling!");
      sum = fPtHardBinScaling.GetSize();
      for (Int_t i = 0; i < fPtHardBinScaling.GetSize(); i++) 
        fPtHardBinScaling[i] /= sum;
    }
    else {
      for (Int_t i = 0; i < fPtHardBinScaling.GetSize(); i++) 
        fPtHardBinScaling[i] /= sum;
    }
  }

  fPtHardBinParam = static_cast<TParameter<int>*>(InputEvent()->FindListObject("PYTHIAPtHardBin"));
  if (!fPtHardBinParam) {
    fPtHardBinParam = new TParameter<int>("PYTHIAPtHardBin", 0);
    AliDebug(3,"Adding pt hard bin param object to the event list...");
    InputEvent()->AddObject(fPtHardBinParam);
  }

  return AliJetEmbeddingFromAODTask::ExecOnce();
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromPYTHIATask::GetNextEntry()
{
  if (fCurrentPtHardBin < 0 || (fMinEntriesPerPtHardBin >= 0 && fPtHardBinCount >= fMinEntriesPerPtHardBin)) {
    fPtHardBinCount = 0;

    Int_t newPtHard = GetRandomPtHardBin();
    
    if (newPtHard != fCurrentPtHardBin) {
      fPtHardBinParam->SetVal(newPtHard);
      fCurrentPtHardBin = newPtHard;
      if (!OpenNextFile()) return kFALSE;
    }
  }

  fPtHardBinCount++;
  if (fHistPtHardBins) fHistPtHardBins->SetBinContent(fCurrentPtHardBin+1, fHistPtHardBins->GetBinContent(fCurrentPtHardBin+1)+1);

  return AliJetEmbeddingFromAODTask::GetNextEntry();
}

//________________________________________________________________________
Int_t AliJetEmbeddingFromPYTHIATask::GetRandomPtHardBin() 
{
  static Int_t order[20]={-1};

  if (order[0] == -1)
    TMath::Sort(fPtHardBinScaling.GetSize(), fPtHardBinScaling.GetArray(), order);

  Double_t rnd = gRandom->Rndm();
  Double_t sum = 0;
  Int_t ptHard = -1;
  for (Int_t i = 0; i < fPtHardBinScaling.GetSize(); i++) {
    sum += fPtHardBinScaling[order[i]];
    if (sum >= rnd) {
      ptHard = order[i];
      break;
    }
  }

  return ptHard;
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromPYTHIATask::UserNotify()
{
  if (!fLHC11hAnchorRun)
    return kTRUE;
  
  Int_t runNumber = InputEvent()->GetRunNumber();

  Int_t semiGoodRunList[32] = {169975, 169981, 170038, 170040, 170083, 170084, 170085, 
			       170088, 170089, 170091, 170152, 170155, 170159, 170163, 
			       170193, 170195, 170203, 170204, 170228, 170230, 170268, 
			       170269, 170270, 170306, 170308, 170309, 169238, 169160, 
			       169156, 169148, 169145, 169144};

  fAnchorRun = 169838; // Assume it is a good run

  for (Int_t i = 0; i < 32; i++) {
    if (runNumber == semiGoodRunList[i]) { // If it is semi good, change the anchor run
      fAnchorRun = 170040;
      break;
    }
  }

  return kTRUE;
}

//________________________________________________________________________
TFile* AliJetEmbeddingFromPYTHIATask::GetNextFile() 
{
  if (fDebugEmbedding == kTRUE) {
    // Embed files in order
    fCurrentAODFileID++;
    // Ensure only files up to the total number of files are used
    if (fCurrentAODFileID > fTotalFiles) {
      fCurrentAODFileID = 1;
    }
  }
  else {
    fCurrentAODFileID = TMath::Nint(gRandom->Rndm()*(fTotalFiles-1))+1;
  }

  if (fMinEntriesPerPtHardBin < 0) {
    fCurrentPtHardBin = GetRandomPtHardBin();
    fPtHardBinParam->SetVal(fCurrentPtHardBin);
  }

  TString fileName;

  if (fAnchorRun>0)
    fileName = Form(fPYTHIAPath.Data(), fAnchorRun, fCurrentPtHardBin, fCurrentAODFileID);
  else
    fileName = Form(fPYTHIAPath.Data(), fCurrentPtHardBin, fCurrentAODFileID);

  if (fFileTable && fFileTable->GetEntries() > 0) {
    TObject* obj = fFileTable->FindObject(fileName);
    if (obj != 0 && fUseAsVetoTable) {
      AliWarning(Form("File %s found in the vetoed file table. Skipping...", fileName.Data()));
      return 0;
    } 
    if (obj == 0 && !fUseAsVetoTable) {
      AliWarning(Form("File %s not found in the allowed file table. Skipping...", fileName.Data()));
      return 0;
    }
  }

  if (fileName.BeginsWith("alien://") && !gGrid) {
    AliInfo("Trying to connect to AliEn ...");
    TGrid::Connect("alien://");
  }

  TString baseFileName(fileName);
  if (baseFileName.Contains(".zip#")) {
    Ssiz_t pos = baseFileName.Last('#');
    baseFileName.Remove(pos);
  }
  
  if (gSystem->AccessPathName(baseFileName)) {
    AliError(Form("File %s does not exist!", baseFileName.Data()));
    return 0;
  }

  AliDebug(3,Form("Trying to open file %s...", fileName.Data()));
  TFile *file = TFile::Open(fileName);

  if (!file || file->IsZombie()) {
    AliError(Form("Unable to open file: %s!", fileName.Data()));
    return 0;
  }

  return file;
}
