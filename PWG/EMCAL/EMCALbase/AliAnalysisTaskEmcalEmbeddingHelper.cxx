/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>
#include <TGrid.h>
#include <TSystem.h>

#include <AliLog.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalEmbeddingHelper);
/// \endcond

AliAnalysisTaskEmcalEmbeddingHelper* AliAnalysisTaskEmcalEmbeddingHelper::fgInstance = 0;

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalEmbeddingHelper::AliAnalysisTaskEmcalEmbeddingHelper() :
  AliAnalysisTaskSE(),
  fTreeName(),
  fPYTHIAPath(),
  fAnchorRun(169838),
  fPtHardBin(-1),
  fAttempts(0),
  fRandomAccess(kFALSE),
  fExternalFile(0),
  fExternalTree(0),
  fCurrentEntry(-1),
  fFirstEntry(-1),
  fLastEntry(-1),
  fExternalEvent(0)
{
  if (fgInstance != 0) {
    AliError("An instance of AliAnalysisTaskEmcalEmbeddingHelper already exists: it will be deleted!!!");
    delete fgInstance;
  }

  fgInstance = this;
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalEmbeddingHelper::AliAnalysisTaskEmcalEmbeddingHelper(const char *name) :
  AliAnalysisTaskSE(name),
  fTreeName("aodTree"),
  fPYTHIAPath("alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/%04d/AliAOD.root"),
  fAnchorRun(169838),
  fPtHardBin(-1),
  fAttempts(10),
  fRandomAccess(kFALSE),
  fExternalFile(0),
  fExternalTree(0),
  fCurrentEntry(-1),
  fFirstEntry(-1),
  fLastEntry(-1),
  fExternalEvent(0)
{
  if (fgInstance != 0) {
    AliError("An instance of AliAnalysisTaskEmcalEmbeddingHelper already exists: it will be deleted!!!");
    delete fgInstance;
  }

  fgInstance = this;
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalEmbeddingHelper::~AliAnalysisTaskEmcalEmbeddingHelper()
{
  if (fgInstance == this) fgInstance = 0;
  if (fExternalEvent) delete fExternalEvent;
  if (fExternalFile) {
    fExternalFile->Close();
    delete fExternalFile;
  }
}

/**
 * Open the next file for embedding
 * This function attempts to open the next available file for embedding
 * by calling GetNextFile(). If a pointer to a file is successfully returned
 * it gets the tree from it, then resets current, first and last entry values.
 * @return kTRUE is the next file is successfully open
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::OpenNextFile()
{
  if (fExternalFile) {
    fExternalFile->Close();
    delete fExternalFile;
    fExternalFile = 0;
  }

  Int_t i = 0;

  while ((!fExternalFile || fExternalFile->IsZombie()) && i < fAttempts) {
    fExternalFile = GetNextFile();
    i++;
  }

  if (!fExternalFile || fExternalFile->IsZombie()) {
    AliError("Could not open AOD file to embed!");
    return kFALSE;
  }

  fExternalTree = static_cast<TTree*>(fExternalFile->Get(fTreeName));
  if (!fExternalTree) {
    AliError(Form("Could not get tree %s from file %s", fTreeName.Data(), fExternalFile->GetName()));
    return kFALSE;
  }

  Bool_t res = InitEvent();
  if (!res) return kFALSE;

  if (fRandomAccess) {
    fFirstEntry = TMath::Nint(gRandom->Rndm()*fExternalTree->GetEntries())-1;
  }
  else {
    fFirstEntry = -1;
  }

  fLastEntry = fExternalTree->GetEntries();
  fCurrentEntry = fFirstEntry;

  AliDebug(3,Form("Will start embedding from entry %d", fCurrentEntry+1));

  return kTRUE;
}

/**
 * Get the name of the next file available for embedding
 * @return A string containing the file name
 */
TString AliAnalysisTaskEmcalEmbeddingHelper::GetNextFileName()
{
  


  TString fname;
    if (fAnchorRun>0)
    fname = Form(fPYTHIAPath.Data(), fAnchorRun, fPtHardBin, fCurrentEntry);
  else
    fname = Form(fPYTHIAPath.Data(), fPtHardBin, fCurrentEntry);

 
  return fname;
}

/**
 * This function calls GetNextFileName() to obtain the name of the next
 * file available for embedding. Then it attempts to open it, opening the alien
 * connection if needed.
 * @return A pointer to a file
 */
TFile* AliAnalysisTaskEmcalEmbeddingHelper::GetNextFile()
{
  TString fileName = GetNextFileName();

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

/**
 * Get the next event to make it available for embedding.
 * If needed it calls OpenNextFile() to open a file.
 * @return kTRUE if successful
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::GetNextEntry()
{
  Int_t attempts = -1;

  do {
    if (fCurrentEntry+1 >= fLastEntry) { // in case it did not start from the first entry, it will go back
      fLastEntry = fFirstEntry;
      fFirstEntry = -1;
      fCurrentEntry = -1;
    }
     

    if (!fExternalFile || !fExternalTree || fCurrentEntry+1 >= fLastEntry) {
      if (!OpenNextFile()) {
        AliError("Could not open the next file!");
        return kFALSE;
      }
    }

    if (!fExternalTree) {
      AliError("Could not get the tree!");
      return kFALSE;
    }

    fCurrentEntry++;
    fExternalTree->GetEntry(fCurrentEntry);

    attempts++;
    if (attempts == 1000)
      AliWarning("After 1000 attempts no event has been accepted by the event selection (trigger, centrality...)!");

  } while (!IsEventSelected());

  if (!fExternalTree) return kFALSE;

  return kTRUE;
}

/**
 * Performs an event selection on the current external event
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::IsEventSelected()
{
  return kTRUE;
}

/**
 * Initialize the external event reading it from the current external file.
 * @return kTRUE if successful
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::InitEvent()
{
  if (!fExternalTree) return kFALSE;

  if (!fExternalEvent) {
    if (fTreeName == "aodTree") {
      fExternalEvent = new AliAODEvent();
    }
    else if (fTreeName == "esdTree") {
      fExternalEvent = new AliESDEvent();
    }
    else {
      AliError(Form("Tree name %s not recognized!", fTreeName.Data()));
      return kFALSE;
    }
  }

  fExternalEvent->ReadFromTree(fExternalTree, fTreeName);

  return kTRUE;
}


/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::UserCreateOutputObjects()
{

}




/**
 * Run analysis code here
 */
void AliAnalysisTaskEmcalEmbeddingHelper::UserExec(Option_t*)
{
  Bool_t res = GetNextEntry();

  if (!res) {
    AliError("Unable to get the AOD event to embed. Nothing will be embedded.");
    return;
  }
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::Terminate(Option_t*)
{
}
