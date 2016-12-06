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

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TChain.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TUUID.h>

#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliInputEventHandler.h>

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
  fAnchorRun(169838),
  fPtHardBin(-1),
  fRandomEventNumberAccess(kFALSE),
  fRandomFileAccess(kFALSE),
  fFilePattern(""),
  fFileListFilename(""),
  fFilenameIndex(-1),
  fFilenames(),
  fTriggerMask(AliVEvent::kAny),
  fZVertexCut(10),
  fMaxVertexDist(999),
  fExternalFile(0),
  fCurrentEntry(0),
  fLowerEntry(0),
  fUpperEntry(0),
  fOffset(0),
  fMaxNumberOfFiles(0),
  fFileNumber(0),
  fInitializedEmbedding(false),
  fInitializedNewFile(false),
  fWrappedAroundTree(false),
  fChain(0),
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
  fAnchorRun(169838),
  fPtHardBin(-1),
  fRandomEventNumberAccess(kFALSE),
  fRandomFileAccess(kFALSE),
  fFilePattern(""),
  fFileListFilename(""),
  fFilenameIndex(-1),
  fFilenames(),
  fTriggerMask(AliVEvent::kAny),
  fZVertexCut(10),
  fMaxVertexDist(999),
  fExternalFile(0),
  fCurrentEntry(0),
  fLowerEntry(0),
  fUpperEntry(0),
  fOffset(0),
  fMaxNumberOfFiles(0),
  fFileNumber(0),
  fInitializedEmbedding(false),
  fInitializedNewFile(false),
  fWrappedAroundTree(false),
  fChain(0),
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
 * Get the names of the files to embed.
 * The filenames will be used to initialize a TChain
 */
void AliAnalysisTaskEmcalEmbeddingHelper::GetFilenames()
{
  // Retrieve filenames if we don't have them yet.
  if (fFilenames.size() == 0)
  {
    // Handle if fPtHardBin or fAnchorRun are set
    // This will require formatting the file pattern in the proper way to support these substitutions
    if (fPtHardBin != -1) {
      if (fAnchorRun > 0) {
        fFilePattern = TString::Format(fFilePattern, fAnchorRun, fPtHardBin);
      }
      else {
        fFilePattern = TString::Format(fFilePattern, fPtHardBin);
      }
    }

    // Retrieve AliEn filenames
    if (fFilePattern.BeginsWith("alien://") && !gGrid) {
      AliInfo("Trying to connect to AliEn ...");
      TGrid::Connect("alien://");

      if (!gGrid) {
        AliFatal(Form("Cannot access AliEn to retrieve file list with pattern %s!", fFilePattern.Data()));
      }

      AliDebug(2,Form("Trying to retrieve file list from AliEn with pattern file %s...", fFilePattern.Data()));

      // Create a temporary filename based on a UUID to make sure that it doesn't overwrite any files
      if (fFileListFilename == "") {
        TUUID tempUUID;
        fFileListFilename = "fileList";
        fFileListFilename += tempUUID.AsString();
        fFileListFilename += ".txt";
      }

      // Retrieve filenames from alien using alien_find
      TString command = "alien_find -l";
      command += fFilePattern;
      command += " > ";
      command += fFileListFilename;

      // Execute the alien_find command to get the filenames
      AliDebug(2,Form("Trying to retrieve file list from AliEn with alien_find command \"%s\"", command.Data()));
      gSystem->Exec(command.Data());
    }

    std::ifstream inputFile(fFileListFilename);

    // Copy available files into the filenames vector
    // From:: https://stackoverflow.com/a/8365247
    std::copy(std::istream_iterator<std::string>(inputFile),
        std::istream_iterator<std::string>(),
        std::back_inserter(fFilenames));
  }

  if (fFilenames.size() == 0) {
    AliFatal(Form("Filenames from pattern \"%s\" and file list \"%s\" yielded an empty list!", fFilePattern.Data(), fFileListFilename.Data()));
  }

  // Add "#" to files in there are any zip files
  // It is require to open the proper root file within the zip
  for (auto filename : fFilenames)
  {
    if (filename.find(".zip") != std::string::npos && filename.find("#") == std::string::npos) {
      if (fTreeName == "aodTree") {
        filename += "#AliAOD.root";
      }
      else if (fTreeName == "esdTree") {
        filename += "#AliESDs.root";
      }
      else {
        AliError(TString::Format("Filename %s contains \".zip\" and not \"#\", but tree name %s is not recognized. Please check the file list to ensure that the proper path is set.", filename.c_str(), fTreeName.Data()));
      }
    }
  }

  // Determine the next filename index
  // This determines which file is added first to the TChain, thus determining the order of processing
  // Random file access. Only do this if the user has no set the filename index and request random file access
  if (fFilenameIndex == -1 && fRandomFileAccess) {
    // - 1 ensures that we it doesn't overflow
    fFilenameIndex = TMath::Nint(gRandom->Rndm()*fFilenames.size()) - 1;
    AliInfo(TString::Format("Starting with random file number %i!", fFilenameIndex));
  }
  // If not random file access, then start from the beginning
  if (fFilenameIndex >= fFilenames.size() || fFilenameIndex < 0) {
    // Skip notifying on -1 since it will likely be set there due to constructor.
    if (fFilenameIndex != -1) {
      AliWarning(Form("File index %i out of range from 0 to %lu! Resetting to 0!", fFilenameIndex, fFilenames.size()));
    }
    fFilenameIndex = 0;
  }

  AliInfo(TString::Format("Starting with file number %i out of %lu", fFilenameIndex, fFilenames.size()));
}

/**
 * Get the next event (entry) in the TChain to make it available for embedding.
 * If needed it calls InitTree() to setup the next tree within the TChain.
 * @return kTRUE if successful
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::GetNextEntry()
{
  Int_t attempts = -1;

  do {
    // Reset to start of tree
    if (fCurrentEntry == fUpperEntry) {
      fCurrentEntry = fLowerEntry;
      fWrappedAroundTree = true;
    }

    if ((fCurrentEntry < fLowerEntry + fOffset) || !fWrappedAroundTree) {
      // Continue with GetEntry as normal
    }
    else {
      // NOTE: On transition from one file to the next, this calls the next entry that would be expected.
      //       However, if it is for the last file, it tries to GetEntry() of one entry past the end of the last file.
      //       Normally, this would be a problem, however GetEntry() just doesn't fill the fields of an invalid index
      //       instead of throwing an error. So "invalid values" are filled for a file that doesn't exist, but then 
      //       they are immediately replaced by the lines below that reset the access values and re-init the tree.
      //       The benefit of this approach is it simplies file counting (we don't need to carefully increment here
      //       and in InitTree()) and preserves the desired behavior when we are not at the last file.
      InitTree();
    }

    // Load current event
    // Can be a simple less than, because fFileNumber counts from 0.
    if (fFileNumber < fMaxNumberOfFiles) {
      fChain->GetEntry(fCurrentEntry);
    }
    else {
      AliError("====================================================================================================");
      AliError("== No more files available to embed from the TChain! Restarting from the beginning of the TChain! ==");
      AliError("== Be careful to check that this is the desired action!                                           ==");
      AliError("====================================================================================================");

      // Reset the relevant access values
      // fCurrentEntry and fLowerEntry are automatically reset in InitTree()
      fFileNumber = 0;
      fUpperEntry = 0;

      // Re-init back to the start
      InitTree();

      // Access the relevant entry
      // We are certain that fFileNumber is less than fMaxNumberOfFiles, so we are resetting to start
      fChain->GetEntry(fCurrentEntry);
    }
    AliDebug(4, TString::Format("Loading entry %i between %i-%i, starting with offset %i from the lower bound of %i", fCurrentEntry, fLowerEntry, fUpperEntry, fOffset, fLowerEntry));

    // Increment current entry
    fCurrentEntry++;
    
    // Provide a check for number of attempts
    attempts++;
    if (attempts == 1000)
      AliWarning("After 1000 attempts no event has been accepted by the event selection (trigger, centrality...)!");

  } while (!IsEventSelected());

  if (!fChain) return kFALSE;

  return kTRUE;
}

/**
 * Performs an event selection on the current external event
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::IsEventSelected()
{
  // AOD Event selection.
  // TODO: Can we make centrality make sense here?
  /*if (!fEsdTreeMode && fAODHeader) {
    AliAODHeader *aodHeader = static_cast<AliAODHeader*>(fAODHeader);

    // Centrality selection
    if (fMinCentrality >= 0) {
      AliCentrality *cent = aodHeader->GetCentralityP();
      Float_t centVal = cent->GetCentralityPercentile("V0M");
      if (centVal < fMinCentrality || centVal >= fMaxCentrality) {
	AliDebug(2,Form("Event rejected due to centrality selection. Event centrality: %f, centrality range selection: %f to %f",
			centVal, fMinCentrality, fMaxCentrality));
	return kFALSE;
      }
    }
  }*/

  // Trigger selection
  if (fTriggerMask != AliVEvent::kAny) {
    UInt_t res = 0;
    const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(InputEvent());
    if (eev) {
      res = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    } else {
      const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
      if (aev) {
        res = ((AliVAODHeader*)aev->GetHeader())->GetOfflineTrigger();
      }
    }

    if ((res & fTriggerMask) == 0) {
      AliDebug(3, Form("Event rejected due to physics selection. Event trigger mask: %d, trigger mask selection: %d.",
                      res, fTriggerMask));
      return kFALSE;
    }
  }

  // Vertex selection
  Double_t externalVertex[3]={0};
  Double_t inputVertex[3]={0};
  const AliVVertex *externalVert = fExternalEvent->GetPrimaryVertex();
  const AliVVertex *inputVert = InputEvent()->GetPrimaryVertex();
  if (externalVert && inputVert) {
    externalVert->GetXYZ(externalVertex);
    inputVert->GetXYZ(inputVertex);

    if (TMath::Abs(externalVertex[2]) > fZVertexCut) {
      AliDebug(3, Form("Event rejected due to Z vertex selection. Event Z vertex: %f, Z vertex cut: %f",
       externalVertex[2], fZVertexCut));
      return kFALSE;
    }
    Double_t dist = TMath::Sqrt((externalVertex[0]-inputVertex[0])*(externalVertex[0]-inputVertex[0])+(externalVertex[1]-inputVertex[1])*(externalVertex[1]-inputVertex[1])+(externalVertex[2]-inputVertex[2])*(externalVertex[2]-inputVertex[2]));
    if (dist > fMaxVertexDist) {
      AliDebug(3, Form("Event rejected because the distance between the current and embedded verteces is > %f. "
       "Current event vertex (%f, %f, %f), embedded event vertex (%f, %f, %f). Distance = %f",
       fMaxVertexDist, inputVertex[0], inputVertex[1], inputVertex[2], externalVertex[0], externalVertex[1], externalVertex[2], dist));
      return kFALSE;
    }
  }

  // TODO: Can we do selection based on the contents of the external event input objects?
  //       The previous embedding task could do so by directly accessing the elemnts.
  //       Certainly can't do jets (say minPt of leading jet) :because this has to be embedded before them.
  //       See AliJetEmbeddingFromAODTask::IsAODEventSelected()

  return kTRUE;
}

/**
 * Initialize the external event reading it from the TChain.
 * @return kTRUE if successful
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::InitEvent()
{
  if (!fChain) return kFALSE;

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

  fExternalEvent->ReadFromTree(fChain, fTreeName);

  return kTRUE;
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::UserCreateOutputObjects()
{
  SetupEmbedding();
}

/**
 * Setup TChain
 * @return kTRUE if successful in setting up the input files.
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::SetupInputFiles()
{
  // Setup TChain
  fChain = new TChain(fTreeName);

  // Determine whether AliEn is needed
  bool requiresAlien = false;
  for (auto filename : fFilenames)
  {
    if (filename.find("alien://") != std::string::npos) {
      requiresAlien = true;
    }
  }

  if (requiresAlien && !gGrid) {
    AliInfo("Trying to connect to AliEn ...");
    TGrid::Connect("alien://");
  }

  // Add files for TChain
  // See: https://stackoverflow.com/a/8533198
  bool wrapped = false;
  TString baseFileName("");
  for (auto filename = fFilenames.begin() + fFilenameIndex; (filename != fFilenames.begin() + fFilenameIndex || !wrapped); filename++)
  {
    // Wraps the loop back around to the beginning
    if (filename == fFilenames.end())
    {
      // Explicit check is needed. Otherwise, an offset of 0 would load the 0th entry twice.
      if (fFilenameIndex == 0) {
        break;
      }
      filename = fFilenames.begin();
      wrapped = true;
    }

    // AccessPathName() cannot handle the "#", so we need to strip it to check that the file exists.
    baseFileName = filename->c_str();
    if (baseFileName.Contains(".zip#")) {
      Ssiz_t pos = baseFileName.Last('#');
      baseFileName.Remove(pos);
    }

    // Ensure that the file is accessible 
    if (gSystem->AccessPathName(baseFileName)) {
      AliError(Form("File %s does not exist! Skipping!", baseFileName.Data()));
      // Do not process the file if it is unaccessible, but continue processing
      continue;
    }

    // Add to the Chain
    fChain->Add(filename->c_str());
  }

  // Keep track of the total number of files in the TChain to ensure that we don't start repeating within the chain
  fMaxNumberOfFiles = fChain->GetListOfFiles()->GetEntries();

  if (fFilenames.size() > fMaxNumberOfFiles) {
    AliWarning(TString::Format("Number of input files (%lu) is larger than the number of available files (%i). Some filenames were likely invalid!", fFilenames.size(), fMaxNumberOfFiles));
  }

  // Setup input event
  Bool_t res = InitEvent();
  if (!res) return kFALSE;

  return kTRUE;
}

/**
 * Sets up the TChain needed to embed external files.
 *
 * Sets fInitializedEmbedding to kTRUE when the initialization is completed.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::SetupEmbedding()
{
  // Get file list
  GetFilenames();

  // Setup TChain
  Bool_t res = SetupInputFiles();
  if (!res) { return; }

  // Note if getting random event access
  if (fRandomEventNumberAccess) {
    AliInfo("Random event number access enabled!");
  }
  
  fInitializedEmbedding = kTRUE;
}

/**
 * Determines the limits of the current TTree within the TChain. By carefully keeping track of
 * the first and lest entry of the current tree in the opened file, we allow a random entry
 * point into the entries in the tree without jumping between files, which would not perform
 * well on the grid.
 *
 * This function sets up random entry point into the file if requested.
 *
 * Sets fInitializedNewFile to kTRUE when the new entries in the tree have been initialized.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::InitTree()
{
  // Load first entry of the (next) file so that we can query information about it
  // (it is unaccessible otherwise).
  // Since fUpperEntry is the total number of entries, loading it will retrieve the
  // next tree (in the next file) since entries are indexed starting from 0.
  fChain->GetEntry(fUpperEntry);

  // Determine tree size and current entry
  // Set the limits of the new tree
  fLowerEntry = fUpperEntry;
  // Fine to be += as long as we started at 0
  fUpperEntry += fChain->GetTree()->GetEntries();

  // Jump ahead at random if desired
  // Determines the offset into the tree
  if (fRandomEventNumberAccess) {
    fOffset = TMath::Nint(gRandom->Rndm()*(fUpperEntry-fLowerEntry))-1;
  }
  else {
    fOffset = 0;
  }

  // Sets which entry to start if the try
  fCurrentEntry = fLowerEntry + fOffset;

  // Keep track of the number of files that we have gone through
  // To start from 0, we only increment if fLowerEntry > 0
  if (fLowerEntry > 0) {
    fFileNumber++;
  }

  AliDebug(2, TString::Format("Will start embedding file %i beginning from entry %i (entry %i within the file). NOTE: This file number is not equal to the absolute file number in the file list!", fFileNumber, fCurrentEntry, fCurrentEntry - fLowerEntry));
  // NOTE: Cannot use this print message, as it is possible that fMaxNumberOfFiles != fFilenames.size() because
  //       invalid filenames may be included in the fFilenames count!
  //AliDebug(2, TString::Format("Will start embedding file %i as the %ith file beginning from entry %i.", (fFilenameIndex + fFileNumber) % fMaxNumberOfFiles, fFileNumber, fCurrentEntry));

  // (re)set whether we have wrapped the tree
  fWrappedAroundTree = false;

  // Note that the tree in the new file has been initialized
  fInitializedNewFile = kTRUE;
}

/**
 * Run analysis code here
 */
void AliAnalysisTaskEmcalEmbeddingHelper::UserExec(Option_t*)
{
  if (!fInitializedEmbedding) {
    AliError("Chain not initialized before running! Setting up now.");
    SetupEmbedding();
  }

  if (!fInitializedNewFile) {
    InitTree();
  }

  Bool_t res = GetNextEntry();

  if (!res) {
    AliError("Unable to get the event to embed. Nothing will be embedded.");
    return;
  }
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::Terminate(Option_t*)
{
}

AliAnalysisTaskEmcalEmbeddingHelper * AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper()
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalEmbeddingHelper", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalEmbeddingHelper", "This task requires an input event handler");
    return 0;
  }

  TString name = "AliAnalysisTaskEmcalEmbeddingHelper";

  AliAnalysisTaskEmcalEmbeddingHelper * mgrTask = static_cast<AliAnalysisTaskEmcalEmbeddingHelper *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  // Create the task that manages
  AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = new AliAnalysisTaskEmcalEmbeddingHelper(name.Data());

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(embeddingHelper);

  // Create containers for input/output
  AliAnalysisDataContainer* cInput = mgr->GetCommonInputContainer();

  /*TString outputContainerName(name);
  outputContainerName += "_histos";

  AliAnalysisDataContainer * cOutput = mgr->CreateContainer(outputContainerName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));*/

  mgr->ConnectInput(embeddingHelper, 0, cInput);
  //mgr->ConnectOutput(embeddingHelper, 1, cOutput);

  return embeddingHelper;
}
