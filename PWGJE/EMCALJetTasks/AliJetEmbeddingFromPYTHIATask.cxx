// $Id: AliJetEmbeddingFromPYTHIATask.cxx $
//
// Jet embedding from PYTHIA task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetEmbeddingFromPYTHIATask.h"

#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TRandom.h>

#include "AliVEvent.h"
#include "AliLog.h"

ClassImp(AliJetEmbeddingFromPYTHIATask)

//________________________________________________________________________
AliJetEmbeddingFromPYTHIATask::AliJetEmbeddingFromPYTHIATask() : 
  AliJetEmbeddingFromAODTask("AliJetEmbeddingFromPYTHIATask", kFALSE),
  fPYTHIAPath(),
  fPtHardBinScaling(),
  fLHC11hAnchorRun(kTRUE),
  fAnchorRun(-1),
  fCurrentPtHardBin(-1)
{
  // Default constructor.
  SetSuffix("PYTHIAEmbedding");
  fTotalFiles = 2000;
  fRandomAccess = kTRUE;
}

//________________________________________________________________________
AliJetEmbeddingFromPYTHIATask::AliJetEmbeddingFromPYTHIATask(const char *name) : 
  AliJetEmbeddingFromAODTask(name, kFALSE),
  fPYTHIAPath("/alice/sim/2012/LHC12a15e"),
  fPtHardBinScaling(),
  fLHC11hAnchorRun(kTRUE),
  fAnchorRun(-1),
  fCurrentPtHardBin(-1)
{
  // Standard constructor.
  SetSuffix("PYTHIAEmbedding");
  fTotalFiles = 2000;
  fRandomAccess = kTRUE;
}

//________________________________________________________________________
AliJetEmbeddingFromPYTHIATask::~AliJetEmbeddingFromPYTHIATask()
{
  // Destructor
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
    }
    
    for (Int_t i = 0; i < fPtHardBinScaling.GetSize(); i++) 
      fPtHardBinScaling[i] /= sum;
  }

  return AliJetEmbeddingFromAODTask::ExecOnce();
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromPYTHIATask::GetNextEntry() 
{
  Int_t newPtHard = GetRandomPtHardBin();

  if (newPtHard != fCurrentPtHardBin) {
    fCurrentPtHardBin = newPtHard;
    if (!OpenNextFile()) return kFALSE;
  }

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

  Int_t semiGoodRunList[28] = {169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 
			       170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 
			       170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 
			       170270, 170306, 170308, 170309};

  fAnchorRun = 169838; // Assume it is a good run

  for (Int_t i = 0; i < 28; i++) {
    if (runNumber == semiGoodRunList[i]) { // If it is semi good, change the anchor run
      fAnchorRun = 170040;
      break;
    }
  }

  return kTRUE;
}

//________________________________________________________________________
TString AliJetEmbeddingFromPYTHIATask::GetNextFileName() 
{
  fCurrentAODFileID = TMath::Nint(gRandom->Rndm()*(fTotalFiles-1))+1;

  TString fileName(fPYTHIAPath);
  if (!fileName.EndsWith("/"))
    fileName += "/";
  
  fileName += fAnchorRun;
  fileName += "/";
  fileName += fCurrentPtHardBin;
  fileName += "/";
  if (fCurrentAODFileID < 10)
    fileName += "00";
  else if (fCurrentAODFileID < 100)
    fileName += "0";
  fileName += fCurrentAODFileID;
  fileName += "/AliAOD.root";

  return fileName;
}
