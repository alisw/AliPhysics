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
#include <sstream>
#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <iostream>
#include <bitset>

#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TChain.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TSystem.h>
#include <TUUID.h>
#include <TKey.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TList.h>

#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliInputEventHandler.h>
#include <AliVHeader.h>
#include <AliAODMCHeader.h>
#include <AliAODMCParticle.h>
#include <AliGenPythiaEventHeader.h>

#include "AliYAMLConfiguration.h"
#include "AliEmcalList.h"
#include "AliEmcalContainerUtils.h"

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalEmbeddingHelper);
/// \endcond

/**
 * Helper function to connect to AliEn (so the code doesn't need to be duplicated).
 */
void ConnectToAliEn()
{
  if (!gGrid) {
    AliInfoGeneralStream("AliAnalysisTaskEmcalEmbeddingHelper") << "Trying to connect to AliEn ...\n";
    TGrid::Connect("alien://");
  }
  if (!gGrid) {
    AliFatalGeneral("AliAnalysisTaskEmcalEmbeddingHelper", "Cannot access AliEn!");
  }
}

/**
 * Helper function to check if a given filename is accessible.
 *
 * @param[in] filename Filename to be checked
 * @return true if the file exists.
 */
bool IsFileAccessible(std::string filename)
{
  // Connect to AliEn if necessary
  // Usually, `gGrid` will exist, so we won't need to waste time on find()
  if (!gGrid && filename.find("alien://") != std::string::npos) {
    ::ConnectToAliEn();
  }

  // AccessPathName() cannot handle the "#", so we need to strip it to check that the file exists.
  if (filename.find(".zip#") != std::string::npos) {
    std::size_t pos = filename.find_last_of("#");
    filename.erase(pos);
  }

  // AccessPathName() has an odd return value - false means that the file exists.
  // NOTE: This is extremely inefficienct for TAlienSystem. It calls
  //       -> gapi_access() -> gapi_stat(). gapi_stat() calls "ls -l" on the basename directory,
  //       which can cause a load on AliEn.
  bool res = gSystem->AccessPathName(filename.c_str());
  // Normalize the result to true if file exists (needed because of the odd return value)
  res = (res == false);
  if (res == false) {
    AliDebugGeneralStream("AliAnalysisTaskEmcalEmbeddingHelper", 4) << "File \"" << filename << "\" doees not exist!\n";
  }

  return res;
}

AliAnalysisTaskEmcalEmbeddingHelper* AliAnalysisTaskEmcalEmbeddingHelper::fgInstance = nullptr;

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalEmbeddingHelper::AliAnalysisTaskEmcalEmbeddingHelper() :
  AliAnalysisTaskSE(),
  fTriggerMask(0),
  fMCRejectOutliers(false),
  fPtHardJetPtRejectionFactor(4),
  fZVertexCut(10),
  fMaxVertexDist(999),
  fInitializedConfiguration(false),
  fInitializedNewFile(false),
  fInitializedEmbedding(false),
  fWrappedAroundTree(false),
  fTreeName(""),
  fNPtHardBins(1),
  fPtHardBin(-1),
  fRandomEventNumberAccess(kFALSE),
  fRandomFileAccess(kTRUE),
  fCreateHisto(true),
  fYAMLConfig(),
  fUseInternalEventSelection(false),
  fUseManualInternalEventCuts(false),
  fInternalEventCuts(),
  fEmbeddedEventUsed(true),
  fValidatedPhysicsSelection(false),
  fInternalEventTriggerMask(0),
  fCentMin(-999),
  fCentMax(-999),
  fRandomRejectionFactor(1.),
  fRandom(0),
  fAutoConfigurePtHardBins(false),
  fAutoConfigureBasePath(""),
  fAutoConfigureTrainTypePath(""),
  fAutoConfigureIdentifier(""),
  fFilePattern(""),
  fInputFilename(""),
  fFileListFilename(""),
  fFilenameIndex(-1),
  fFilenames(),
  fConfigurationPath(""),
  fEmbeddedRunlist(),
  fPythiaCrossSectionFilenames(),
  fExternalFile(nullptr),
  fChain(nullptr),
  fCurrentEntry(0),
  fLowerEntry(0),
  fUpperEntry(0),
  fOffset(0),
  fMaxNumberOfFiles(0),
  fFileNumber(0),
  fHistManager(),
  fOutput(nullptr),
  fExternalEvent(nullptr),
  fExternalMCEvent(nullptr),
  fExternalHeader(nullptr),
  fPythiaHeader(nullptr),
  fPythiaTrials(0),
  fPythiaTrialsFromFile(0),
  fPythiaCrossSection(0.),
  fPythiaCrossSectionFromFile(0.),
  fPythiaPtHard(0.),
  fPrintTimingInfoToLog(false),
  fTimer()
{
  if (fgInstance != nullptr) {
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
  fTriggerMask(0),
  fMCRejectOutliers(false),
  fPtHardJetPtRejectionFactor(4),
  fZVertexCut(10),
  fMaxVertexDist(999),
  fInitializedConfiguration(false),
  fInitializedNewFile(false),
  fInitializedEmbedding(false),
  fWrappedAroundTree(false),
  fTreeName("aodTree"),
  fNPtHardBins(1),
  fPtHardBin(-1),
  fRandomEventNumberAccess(kFALSE),
  fRandomFileAccess(kTRUE),
  fCreateHisto(true),
  fYAMLConfig(),
  fUseInternalEventSelection(false),
  fUseManualInternalEventCuts(false),
  fInternalEventCuts(),
  fEmbeddedEventUsed(true),
  fValidatedPhysicsSelection(false),
  fInternalEventTriggerMask(0),
  fCentMin(-999),
  fCentMax(-999),
  fRandomRejectionFactor(1.),
  fRandom(0),
  fAutoConfigurePtHardBins(false),
  fAutoConfigureBasePath("alien:///alice/cern.ch/user/a/alitrain/"),
  fAutoConfigureTrainTypePath("PWGJE/Jets_EMC_PbPb/"),
  fAutoConfigureIdentifier("autoConfigIdentifier"),
  fFilePattern(""),
  fInputFilename(""),
  fFileListFilename(""),
  fFilenameIndex(-1),
  fFilenames(),
  fConfigurationPath(""),
  fEmbeddedRunlist(),
  fPythiaCrossSectionFilenames(),
  fExternalFile(nullptr),
  fChain(nullptr),
  fCurrentEntry(0),
  fLowerEntry(0),
  fUpperEntry(0),
  fOffset(0),
  fMaxNumberOfFiles(0),
  fFileNumber(0),
  fHistManager(name),
  fOutput(nullptr),
  fExternalEvent(nullptr),
  fExternalMCEvent(nullptr),
  fExternalHeader(nullptr),
  fPythiaHeader(nullptr),
  fPythiaTrials(0),
  fPythiaTrialsFromFile(0),
  fPythiaCrossSection(0.),
  fPythiaCrossSectionFromFile(0.),
  fPythiaPtHard(0.),
  fPrintTimingInfoToLog(false),
  fTimer()
{
  if (fgInstance != 0) {
    AliError("An instance of AliAnalysisTaskEmcalEmbeddingHelper already exists: it will be deleted!!!");
    delete fgInstance;
  }

  fgInstance = this;

  if (fCreateHisto) {
    DefineOutput(1, AliEmcalList::Class());
  }
}

/**
 * Destructor
 *
 * Deletes the singleton instance and ensures that any open file is closed.
 */
AliAnalysisTaskEmcalEmbeddingHelper::~AliAnalysisTaskEmcalEmbeddingHelper()
{
  if (fgInstance == this) fgInstance = nullptr;
  if (fExternalEvent) delete fExternalEvent;
  if (fExternalMCEvent) delete fExternalMCEvent;
  if (fExternalFile) {
    fExternalFile->Close();
    delete fExternalFile;
  }
}

bool AliAnalysisTaskEmcalEmbeddingHelper::Initialize(bool removeDummyTask)
{
  // Initialize %YAML configuration, if one is given
  bool initializedYAML = InitializeYamlConfig();

  RetrieveTaskPropertiesFromYAMLConfig();
  
  // Get file list
  bool result = GetFilenames();

  if (result && initializedYAML) {
    fInitializedConfiguration = true;
  }

  if (removeDummyTask == true) {
    RemoveDummyTask();
  }

  // Initialize the YAML config object for streaming
  fYAMLConfig.Initialize();

  // Print the results of the initialization
  // Print outside of the ALICE Log system to ensure that it is always available!
  std::cout << *this;

  return result;
}

/**
 * Retrieve embedding helper properties from a %YAML configuration file.
 *
 * Unlike the yaml configuration used in the correction task, there is no "default" yaml configuration
 * file - nothing is required to exist in the yaml file
 */
void AliAnalysisTaskEmcalEmbeddingHelper::RetrieveTaskPropertiesFromYAMLConfig()
{
  // Following the variable blocks defined in the header
  // Embedded event properties
  std::vector<std::string> physicsSelection;
  bool res = fYAMLConfig.GetProperty("embeddedEventPhysicsSelection", physicsSelection, false);
  if (res) {
    fTriggerMask = AliEmcalContainerUtils::DeterminePhysicsSelectionFromYAML(physicsSelection);
  }
  res = fYAMLConfig.GetProperty("enableMCOutlierRejection", fMCRejectOutliers, false);
  res = fYAMLConfig.GetProperty("ptHardJetPtRejectionFactor", fPtHardJetPtRejectionFactor, false);
  res = fYAMLConfig.GetProperty("embeddedEventZVertexCut", fZVertexCut, false);
  res = fYAMLConfig.GetProperty("maxVertexDifferenceDistance", fMaxVertexDist, false);

  // Embedding helper properties
  res = fYAMLConfig.GetProperty("treeName", fTreeName, false);
  res = fYAMLConfig.GetProperty("nPtHardBins", fNPtHardBins, false);
  res = fYAMLConfig.GetProperty("ptHardBin", fPtHardBin, false);
  res = fYAMLConfig.GetProperty("randomEventNumberAccess", fRandomEventNumberAccess, false);
  res = fYAMLConfig.GetProperty("randomFileAccess", fRandomFileAccess, false);
  res = fYAMLConfig.GetProperty("createHisto", fCreateHisto, false);
  res = fYAMLConfig.GetProperty("printTimingInfoInLog", fPrintTimingInfoToLog, false);
  // More general embedding helper properties
  res = fYAMLConfig.GetProperty("filePattern", fFilePattern, false);
  res = fYAMLConfig.GetProperty("inputFilename", fInputFilename, false);
  res = fYAMLConfig.GetProperty("fileListFilename", fFileListFilename, false);
  res = fYAMLConfig.GetProperty("filenameIndex", fFilenameIndex, false);
  // Configuration path makes no sense, as we are already using the %YAML configuration
  res = fYAMLConfig.GetProperty("runlist", fEmbeddedRunlist, false);
  // Generally should not be set
  res = fYAMLConfig.GetProperty("filenames", fFilenames, false);
  res = fYAMLConfig.GetProperty("fPythiaCrossSectionFilenames", fPythiaCrossSectionFilenames, false);

  // Internal event selection properties
  // NOTE: Need to define the base name here so that the property path is not ambiguous (due to otherwise only being `const char *`)
  std::string baseName = "internalEventSelection";
  res = fYAMLConfig.GetProperty({baseName, "enabled"}, fUseInternalEventSelection, false);
  res = fYAMLConfig.GetProperty({baseName, "useManualCuts"}, fUseManualInternalEventCuts, false);
  // Centrality
  std::vector <double> centralityRange;
  res = fYAMLConfig.GetProperty({baseName, "centralityRange"}, centralityRange, false);
  if (res) {
    if (centralityRange.size() != 2) {
      AliErrorStream() << "Passed centrality range with " << centralityRange.size() << " entries, but 2 values are required. Ignoring values.\n";
    }
    else {
      AliDebugStream(1) << "Setting internal event centrality range to [" << centralityRange.at(0) << ", " << centralityRange.at(1) << "]\n";
      fCentMin = centralityRange.at(0);
      fCentMax = centralityRange.at(1);
    }
  }
  // Physics selection
  res = fYAMLConfig.GetProperty({baseName, "physicsSelection"}, physicsSelection, false);
  if (res) {
    fInternalEventTriggerMask = AliEmcalContainerUtils::DeterminePhysicsSelectionFromYAML(physicsSelection);
  }

  // Auto configure pt hard properties
  res = fYAMLConfig.GetProperty("autoConfigurePtHardBins", fAutoConfigurePtHardBins, false);
  res = fYAMLConfig.GetProperty("autoConfigureBasePath", fAutoConfigureBasePath, false);
  res = fYAMLConfig.GetProperty("autoConfigureTrainTypePath", fAutoConfigureTrainTypePath, false);
  res = fYAMLConfig.GetProperty("autoConfigureIdentifier", fAutoConfigureIdentifier, false);
  // Random rejection 
  res = fYAMLConfig.GetProperty("randomRejectionFactor", fRandomRejectionFactor, false);
}

/**
 * Get the names of the files to embed and determine which file to start from. The filenames will be
 * used to initialize a TChain. Filenames can either be specified in a local file with one filename per
 * line or found on AliEn by specifying any pattern that would work in alien_find.
 *
 * Notes on this function:
 * - In the case that the passed files contain ".zip", the root file name will be appended. Be certain to
 *   set the proper file type (AOD or ESD)!
 * - In the case of AliEn, the result will be copied into a local file.
 * - The file to start from can be randomized if requested.
 *
 * The file pattern can be specified by a number of different options:
 * - The simplest is pass a fully formed path to a file either locally or on AliEn.
 * - You could specify a pT hard bin, which will be added into the passed file pattern. For example,
 *   "/period/runNumber/%d/aod_archive.root", where the %d will the be the pT hard bin.
 * - You could specify a pT hard bin and anchor run, which will both be added into the passed file pattern.
 *   "/period/%d/%d/aod_archive.root", where the first %d will the be the anchor run number, and the second
 *   %d will be the pT hard bin.
 *
 * NOTE: The file pattern only makes sense in terms of AliEn. Passing a local pattern will not work!
 * NOTE: Exercise care if you set both the file pattern and the filename! Doing so will probably cause your
 * file to be overwritten!
 */
bool AliAnalysisTaskEmcalEmbeddingHelper::GetFilenames()
{
  // Determine the pattern filename if not yet set
  if (fInputFilename == "") {
    if (fTreeName == "aodTree") {
      fInputFilename = "AliAOD.root";
    }
    else if (fTreeName == "esdTree") {
      fInputFilename = "AliESDs.root";
    }
    else {
      AliFatal(TString::Format("Requested default (pattern) input filename, but could not determine file type from tree name \"%s\". Please check the tree name and set the pattern input filename", fTreeName.Data()));
    }
  }

  // Retrieve filenames if we don't have them yet.
  if (fFilenames.size() == 0)
  {
    // Handle pt hard bin auto configuration
    if (fAutoConfigurePtHardBins)
    {
      if (fPtHardBin > 0) {
        AliFatal("Requested both pt hard bin auto configuration and selected a non-zero pt hard bin. These are incompatible options. Please check your configuration.");
      }
      bool success = AutoConfigurePtHardBins();
      if (success == false) {
        AliFatal("Pt hard bin auto configuration requested, but it failed. Please check the logs.\n");
      }
    }

    // Handle if fPtHardBin is set
    // This will require formatting the file pattern in the proper way to support these substitutions
    if (fPtHardBin != -1 && fFilePattern != "") {
      fFilePattern = TString::Format(fFilePattern, fPtHardBin);
    }

    // Setup AliEn access if needed
    if (fFilePattern.Contains("alien://") || fFileListFilename.Contains("alien://")) {
      ::ConnectToAliEn();
    }

    // Retrieve AliEn filenames directly from AliEn
    bool usedFilePattern = false;
    if (fFilePattern.Contains("alien://")) {
      usedFilePattern = true;
      AliDebug(2, TString::Format("Trying to retrieve file list from AliEn with pattern file %s...", fFilePattern.Data()));

      // Create a temporary filename based on a UUID to make sure that it doesn't overwrite any files
      if (fFileListFilename == "") {
        fFileListFilename += GenerateUniqueFileListFilename();
      }

      // The query command cannot handle "alien://" in the file pattern, so we need to remove it for the command
      TString filePattern = fFilePattern;
      filePattern.ReplaceAll("alien://", "");

      // Execute the grid query to get the filenames
      AliDebug(2, TString::Format("Trying to retrieve file list from AliEn with pattern \"%s\" and input filename \"%s\"", filePattern.Data(), fInputFilename.Data()));
      auto result = gGrid->Query(filePattern.Data(), fInputFilename.Data());

      if (result) {

        // Loop over the result to store it in the fileList file
        std::ofstream outFile(fFileListFilename);
        for (int i = 0; i < result->GetEntries(); i++)
        {
          TString path = result->GetKey(i, "turl");
          // "turl" corresponds to the full AliEn url
          
          // If a runlist is specified for good embedded runs, only include the file if it is in this runlist
          if (IsRunInRunlist(path.Data())) {
            outFile << path << "\n";
          }
        }
        outFile.close();
      }
      else {
        AliErrorStream() << "Failed to run grid query\n";
        return false;
      }
    }

    // Handle a filelist on AliEn
    if (fFileListFilename.Contains("alien://")) {
      // Check if we already used the file pattern
      if (usedFilePattern) {
        AliErrorStream() << "You set both the file pattern and the file list filename! The file list filename will override the pattern! Pattern: \"" << fFilePattern << "\", filename: \"" << fFileListFilename << "\"\nPlease check that this is the desired behavior!\n";
      }

      // Determine the local filename and copy file to local directory
      std::string alienFilename = fFileListFilename.Data();
      fFileListFilename = gSystem->BaseName(alienFilename.c_str());
      fFileListFilename += GenerateUniqueFileListFilename();

      TFile::Cp(alienFilename.c_str(), fFileListFilename.Data());
    }

    std::ifstream inputFile(fFileListFilename);

    // Copy available files into the filenames vector
    // From:: https://stackoverflow.com/a/8365247
    std::copy(std::istream_iterator<std::string>(inputFile),
        std::istream_iterator<std::string>(),
        std::back_inserter(fFilenames));

    inputFile.close();
  }

  if (fFilenames.size() == 0) {
    AliFatal(TString::Format("Filenames from pattern \"%s\" and file list \"%s\" yielded an empty list!", fFilePattern.Data(), fFileListFilename.Data()));
  }

  // Add "#" to files in there are any zip files
  // It is require to open the proper root file within the zip
  for (auto filename : fFilenames)
  {
    if (filename.find(".zip") != std::string::npos && filename.find("#") == std::string::npos) {
      filename += "#";
      filename += fInputFilename;
      if (fTreeName == "aodTree") {
        filename += "#AliAOD.root";
      }
      else if (fTreeName == "esdTree") {
        filename += "#AliESDs.root";
      }
      else {
        AliError(TString::Format("Filename %s contains \".zip\" and not \"#\", but tree name %s is not recognized. Please check the file list to ensure that the proper path is set.", filename.c_str(), fTreeName.Data()));
        return false;
      }
    }
  }

  // Determine whether AliEn is needed
  // It is possible that this has not been determined up to this point
  for (auto filename : fFilenames)
  {
    if (filename.find("alien://") != std::string::npos) {
      ::ConnectToAliEn();
      // No point in continuing to search once we know that it is needed
      break;
    }
  }

  // Check if each filenames exists. If they do not exist, then remove them for fFilenames
  unsigned int initialSize = fFilenames.size();
  // NOTE: We invert the result of IsFileAccessible because we should return true for files that should be _removed_ (ie are inaccessible)
  fFilenames.erase(std::remove_if(fFilenames.begin(), fFilenames.end(), [](const std::string & filename) {return (::IsFileAccessible(filename) == false);} ), fFilenames.end());

  AliInfoStream() << "Found " << fFilenames.size() << " files to embed (" << (initialSize - fFilenames.size()) << " filename(s) inaccessible or invalid)\n";

  // Determine pythia filename
  DeterminePythiaXSecFilename();

  return true;
}

/**
 * Determine the Pythia cross section filename by checking for the existence of files with various possible filenames.
 * Note that it uses the first input filename as a proxy for all other input files following the same pattern.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::DeterminePythiaXSecFilename()
{
  // Get the initial filename. Use the first entry as a proxy for other input files
  std::string externalEventFilename = "";
  if (fFilenames.size() > 0) {
    externalEventFilename = fFilenames.at(0);
  }
  else {
    return;
  }

  std::vector <std::string> pythiaBaseFilenames = {"pyxsec.root", "pyxsec_hists.root"};
  AliInfoStream() << "Attempting to determine pythia cross section filename. It can be normal to see some TFile::Init() errors!\n";
  std::string pythiaXSecFilename = "";
  for (auto & name : pythiaBaseFilenames) {
    pythiaXSecFilename = ConstructFullPythiaXSecFilename(externalEventFilename, name, true);
    if (pythiaXSecFilename != "") {
      AliDebugStream(4) << "Found pythia cross section filename \"" << name << "\"\n";
      fPythiaXSecFilename = name;
      break;
    }
  }

  if (fPythiaXSecFilename == "") {
    // Failed entirely - just give up on this
    // We will use an empty filename as a proxy for whether the file has been found (empty is equivalent to not found)
    AliErrorStream() << "Failed to find pythia x sec file! Continuing with only the pythia header!\n";
  }
  else {
    AliInfoStream() << "Found pythia cross section file \"" << fPythiaXSecFilename << "\".\n";
  }
}


/**
 * Check if a given filename is from a run in the good embedded runlist. If no runlist was defined,
 * it will always return true.
 *
 * @param path path of a single filename
 * @return true if the path contains a run in the good embedded runlist.
 */
bool AliAnalysisTaskEmcalEmbeddingHelper::IsRunInRunlist(const std::string & path) const
{
  if (fEmbeddedRunlist.size() == 0) {
    return true;
  }
  
  for (auto run : fEmbeddedRunlist) {
    if (path.find(run) != std::string::npos) {
      return true;
    }
  }
  return false;
}

/**
 * Initialize the %YAML configuration with a potentially specified %YAML configuration file.
 *
 * @return true if no yaml file or it exists and was successfully accessed.
 */
bool AliAnalysisTaskEmcalEmbeddingHelper::InitializeYamlConfig()
{
  if (fConfigurationPath == "") {
    AliInfo("No Embedding YAML configuration was provided");
  }
  else {
    AliInfoStream() << "Embedding YAML configuration was provided: \"" << fConfigurationPath << "\".\n";

    int addedConfig = fYAMLConfig.AddConfiguration(fConfigurationPath, "yamlConfig");
    if (addedConfig < 0) {
      AliError("YAML Configuration not found!");
      return false;
    }
  }

  return true;
}

/**
 * Handle auto-configuration of pt hard bins on LEGO trains. It gets the train number from the LEGO train
 * environment, thereby assigning a pt hard bin to a particular train. This assignment is written out to a
 * %YAML file so all of the trains can determine which pt hard bins are available. The %YAML file that is written
 * contains a map of ptHardBin to train number.
 *
 * Note that when the number of pt hard bins is exhausted and all trains are assigned, the file is removed.
 *
 * Relevant directory information for an example train:
 *   - PWD=/home/alitrain/train-workdir/PWGJE/Jets_EMC_PbPb/2558_20170930-2042/config
 *   - Grid workdir relative to user $HOME: _________ /alice/cern.ch/user/a/alitrain/PWGJE/Jets_EMC_PbPb/2558_20170930-2042
 *   - Grid output directory relative to workdir: ___ $1/PWGJE/Jets_EMC_PbPb/2558_20170930-2042
 *
 * @return true if the pt hard bin was successfully extracted.
 */
bool AliAnalysisTaskEmcalEmbeddingHelper::AutoConfigurePtHardBins()
{
  bool returnValue = false;

  AliInfoStream() << "Attempting to auto configure pt hard bins.\n";
  // %YAML configuration containing pt hard bin to train number mapping
  PWG::Tools::AliYAMLConfiguration config;

  // Handle AliEn explicitly here since the default base path contains "alien://"
  if (fAutoConfigureBasePath.find("alien://") != std::string::npos && !gGrid) {
    ::ConnectToAliEn();
  }

  // Get train ID
  // Need to get the char * directly because it may be null.
  const char * trainNumberStr = gSystem->Getenv("TRAIN_RUN_ID");
  std::stringstream trainNumberSS;
  if (trainNumberStr) {
    trainNumberSS << trainNumberStr;
  }
  if (trainNumberSS.str() == "") {
    AliFatal("Cannot retrieve train ID.");
  }
  // Extract train number from the string
  int trainNumber;
  trainNumberSS >> trainNumber;

  // Determine the file path
  auto filename = RemoveTrailingSlashes(fAutoConfigureBasePath);
  filename += "/";
  filename += RemoveTrailingSlashes(fAutoConfigureTrainTypePath);
  filename += "/";
  filename += fAutoConfigureIdentifier;
  // Add ".yaml" if it is not already there
  std::string yamlExtension = ".yaml";
  if (filename.find(yamlExtension) == std::string::npos) {
    filename += yamlExtension;
  }

  // Check if file exists
  if (gSystem->AccessPathName(filename.c_str())) {
    // File _does not_ exist
    AliInfoStream() << "Train pt hard bin configuration file not available, so creating a new empty configuration named \"" << fAutoConfigureIdentifier << "\".\n";
    // Use an empty configuration
    config.AddEmptyConfiguration(fAutoConfigureIdentifier);
  }
  else {
    AliInfoStream() << "Opening configuration located at \"" << filename << "\".\n";
    // Use the existing configuration
    config.AddConfiguration(filename, fAutoConfigureIdentifier);
  }

  // Look for each pt hard bin, and then retrieve the corresponding train number
  // Once an open pt hard bin is found, add the current train number
  int tempTrainNumber = -1;
  bool getPropertyReturnValue = false;
  std::stringstream propertyName;
  for (int ptHardBin = 1; ptHardBin <= fNPtHardBins; ptHardBin++)
  {
    propertyName.str("");
    propertyName << ptHardBin;
    getPropertyReturnValue = config.GetProperty(propertyName.str(), tempTrainNumber, false);
    if (getPropertyReturnValue != true) {
      AliInfoStream() << "Train " << trainNumber << " will use pt hard bin " << ptHardBin << ".\n";
      // We have determine our pt hard bin!
      fPtHardBin = ptHardBin;

      // Write the train number back out to the %YAML configuration and save it
      config.WriteProperty(propertyName.str(), trainNumber, fAutoConfigureIdentifier);
      config.WriteConfiguration(filename, fAutoConfigureIdentifier);

      // NOTE: Cannot clean up the YAML file on the last pt hard bin because the train can be launched
      // multiple times due to tests, etc. Therefore, we have to accept that we are leaving around used
      // YAML config files.

      // We are done - continue on.
      returnValue = true;
      break;
    }
    else {
      AliDebugStream(2) << "Found pt hard bin " << ptHardBin << " corresponding to train number " << trainNumber << ".\n";
      // If train was already allocated (say, by a test train), then use that pt hard bin
      if (tempTrainNumber == trainNumber) {
        AliInfoStream() << "Train run number " << trainNumber << " was already found assigned to pt hard bin " << ptHardBin << ". That pt hard bin will be used.\n";
        fPtHardBin = ptHardBin;

        // We are done - continue on.
        returnValue = true;
        break;
      }
      // Otherwise, nothing to be done.
    }
  }

  return returnValue;
}

/**
 * Simple helper function to generate a unique file list filename. This filename will be used to store the
 * filelist locally. It will be of the form fFileListFilename + ".UUID.txt". If fFileListFilename is empty,
 * then the beginning of the filename will be "fileList".
 *
 * @return std::string containing the rest of the file to be appended to fFileListFilename
 */
std::string AliAnalysisTaskEmcalEmbeddingHelper::GenerateUniqueFileListFilename() const
{
  std::string tempStr = "";
  if (fFileListFilename == "") {
    tempStr = "fileList";
  }
  TUUID tempUUID;
  tempStr += ".";
  tempStr += tempUUID.AsString();
  tempStr += ".txt";

  return tempStr;
}

/**
 * Remove slashes at the end of strings. See: https://stackoverflow.com/a/14878124
 *
 * @param[in] filename String containing a filename with some number of extra trailing slashes.
 *
 * @return string without trailing slashes.
 */
std::string AliAnalysisTaskEmcalEmbeddingHelper::RemoveTrailingSlashes(std::string filename) const
{
  while (filename.rbegin() != filename.rend() && *(filename.rbegin()) == '/') {
    filename.pop_back();
  }

  return filename;
}

/**
 * Determine the first file to embed and store the index. The index will either be
 * random or the first file in the list, depending on the task configuration.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::DetermineFirstFileToEmbed()
{
  // This determines which file is added first to the TChain, thus determining the order of processing
  // Random file access. Only do this if the user has no set the filename index and request random file access
  if (fFilenameIndex == -1 && fRandomFileAccess) {
    // Floor ensures that we it doesn't overflow
    TRandom3 rand(0);
    fFilenameIndex = TMath::FloorNint(rand.Rndm()*fFilenames.size());
    // +1 to account for the fact that the filenames vector is 0 indexed.
    AliInfo(TString::Format("Starting with random file number %i!", fFilenameIndex+1));
  }
  // If not random file access, then start from the beginning
  if (fFilenameIndex < 0 || static_cast<UInt_t>(fFilenameIndex) >= fFilenames.size()) {
    // Skip notifying on -1 since it will likely be set there due to constructor.
    if (fFilenameIndex != -1) {
      AliWarning(TString::Format("File index %i out of range from 0 to %lu! Resetting to 0!", fFilenameIndex, fFilenames.size()));
    }
    fFilenameIndex = 0;
  }

  // +1 to account for the fact that the filenames vector is 0 indexed.
  AliInfo(TString::Format("Starting with file number %i out of %lu", fFilenameIndex+1, fFilenames.size()));
}

/**
 * Get the next event (entry) in the TChain to make it available for embedding. The event will be selected
 * according to the conditions determined in IsEventSelected(). If needed it calls InitTree() to setup the
 * next tree within the TChain. In the case of running of out files to embed, an error is thrown and embedding
 * begins again from the start of the file list.
 *
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
      //       However, if it is for the last file, it tries to GetEntry() of one entry past the end of the last
      //       file. Normally, this would be a problem, however GetEntry() just doesn't fill the fields of an
      //       invalid index instead of throwing an error. So "invalid values" are filled for a file that doesn't
      //       exist, but then they are immediately replaced by the lines below that reset the access values and
      //       re-init the tree. The benefit of this approach is it simplies file counting (we don't need to
      //       carefully increment here and in InitTree()) and preserves the desired behavior when we are not at
      //       the last file.
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

    // Set relevant event properties
    SetEmbeddedEventProperties();

    // Increment current entry
    fCurrentEntry++;
    
    // Provide a check for number of attempts
    attempts++;
    if (attempts == 1000)
      AliWarning("After 1000 attempts no event has been accepted by the event selection (trigger, centrality...)!");

    // Record event properties
    if (fCreateHisto) {
      RecordEmbeddedEventProperties();
    }

  } while (!IsEventSelected());

  if (fCreateHisto) {
    fHistManager.FillTH1("fHistEventCount", "Accepted");
    fHistManager.FillTH1("fHistEmbeddedEventsAttempted", attempts);
  }

  if (!fChain) return kFALSE;

  return kTRUE;
}

/**
 * Set some properties of the event that are not immediately available from the external event to make them
 * available to user tasks.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::SetEmbeddedEventProperties()
{
  AliDebug(4, "Set event properties");
  fExternalHeader = fExternalEvent->GetHeader();
  
  // Handle pythia header if AOD
  AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(fExternalEvent->FindListObject(AliAODMCHeader::StdBranchName()));
  if (aodMCH) {
    for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
      fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
      if (fPythiaHeader) break;
    }
    
    // Get MC Event                                                                        
    TClonesArray* mcParticles = static_cast<TClonesArray*> (fExternalEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( mcParticles ) {
      if ( !fExternalMCEvent ) fExternalMCEvent = new AliMCEvent();
      fExternalMCEvent->SetExternalHeader(aodMCH);
      fExternalMCEvent->SetParticleArray(mcParticles);
    }
  }

  if (fPythiaHeader)
  {
    fPythiaCrossSection = fPythiaHeader->GetXsection();
    fPythiaTrials = fPythiaHeader->Trials();
    fPythiaPtHard = fPythiaHeader->GetPtHard();
    // It is identically zero if the cross section is not available
    if (fPythiaCrossSection == 0.) {
      AliDebugStream(4) << "Taking the pythia cross section avg from the xsec file.\n";
      fPythiaCrossSection = fPythiaCrossSectionFromFile;
    }
    // It is identically zero if the number of trials is not available
    if (fPythiaTrials == 0.) {
      AliDebugStream(4) << "Taking the pythia trials avg from the xsec file.\n";
      fPythiaTrials = fPythiaTrialsFromFile;
    }
    // Pt hard is inherently event-by-event and cannot by taken as a avg quantity.

    AliDebugStream(4) << "Pythia header is defined!\n";
    AliDebugStream(4) << "fPythiaCrossSection: " << fPythiaCrossSection << "\n";
  }
}

/**
 * Record event properties
 */
void AliAnalysisTaskEmcalEmbeddingHelper::RecordEmbeddedEventProperties()
{
  // Fill trials, xsec, pt hard
  fHistManager.FillTH1("fHistTrials", fPtHardBin, fPythiaTrials);
  fHistManager.FillProfile("fHistXsection", fPtHardBin, fPythiaCrossSection);
  fHistManager.FillTH1("fHistPtHard", fPythiaPtHard);
}

/**
 * Handles (ie wraps) event selection and proper event counting.
 *
 * @return kTRUE if the event successfully passes all criteria.
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::IsEventSelected()
{
  if (CheckIsEmbeddedEventSelected()) {
    return kTRUE;
  }

  if (fCreateHisto) {
    // Keep count of number of rejected events
    fHistManager.FillTH1("fHistEventCount", "Rejected");
  }

  return kFALSE;
}

/**
 * Performs the embedded event selection on the current external event.
 *
 * @return kTRUE if the event successfully passes all criteria.
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::CheckIsEmbeddedEventSelected()
{
  // Check if pt hard bin is 0, indicating a problem with the event or the grid.
  // In such a case, the event should be rejected.
  // This condition should only be applied if we have a valid pythia header.
  // (pt hard should still be set even if the production wasn't done in pt hard bins).
  if (fPythiaPtHard == 0. && fPythiaHeader) {
    AliDebugStream(3) << "Event rejected due to pt hard = 0, indicating a problem with the external event.\n";
    if (fCreateHisto) {
      fHistManager.FillTH1("fHistEmbeddedEventRejection", "PtHardIs0", 1);
    }
    return kFALSE;
  }

  // Physics selection
  if (fTriggerMask != 0) {
    UInt_t res = 0;
    const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(fExternalEvent);
    if (eev) {
      AliFatal("Event selection is not implemented for embedding ESDs.");
      // Unfortunately, the normal method of retrieving the trigger mask (commented out below) doesn't work for the embedded event since we don't
      // create an input handler and I am not an expert on getting a trigger mask. Further, embedding ESDs is likely to be inefficient, so it is
      // probably best to avoid it if possible.
      //
      // Suggestions are welcome here!
      //res = (dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    } else {
      const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(fExternalEvent);
      if (aev) {
        res = (dynamic_cast<AliVAODHeader*>(aev->GetHeader()))->GetOfflineTrigger();
      }
    }

    if ((res & fTriggerMask) == 0) {
      AliDebug(3, Form("Event rejected due to physics selection. Event trigger mask: %d, trigger mask selection: %d.",
                      res, fTriggerMask));
      if (fCreateHisto) {
        fHistManager.FillTH1("fHistEmbeddedEventRejection", "PhysSel", 1);
      }
      return kFALSE;
    }
  }

  // Vertex selection
  Double_t externalVertex[3]={0};
  Double_t inputVertex[3]={0};
  const AliVVertex *externalVert = fExternalEvent->GetPrimaryVertex();
  const AliVVertex *inputVert = AliAnalysisTaskSE::InputEvent()->GetPrimaryVertex();
  if (externalVert && inputVert) {
    externalVert->GetXYZ(externalVertex);
    inputVert->GetXYZ(inputVertex);

    if (TMath::Abs(externalVertex[2]) > fZVertexCut) {
      AliDebug(3, Form("Event rejected due to Z vertex selection. Event Z vertex: %f, Z vertex cut: %f",
       externalVertex[2], fZVertexCut));
      if (fCreateHisto) {
        fHistManager.FillTH1("fHistEmbeddedEventRejection", "Vz", 1);
      }
      return kFALSE;
    }
    Double_t dist = TMath::Sqrt((externalVertex[0]-inputVertex[0])*(externalVertex[0]-inputVertex[0])+(externalVertex[1]-inputVertex[1])*(externalVertex[1]-inputVertex[1])+(externalVertex[2]-inputVertex[2])*(externalVertex[2]-inputVertex[2]));
    if (dist > fMaxVertexDist) {
      AliDebug(3, Form("Event rejected because the distance between the current and embedded vertices is > %f. "
       "Current event vertex (%f, %f, %f), embedded event vertex (%f, %f, %f). Distance = %f",
       fMaxVertexDist, inputVertex[0], inputVertex[1], inputVertex[2], externalVertex[0], externalVertex[1], externalVertex[2], dist));
      if (fCreateHisto) {
        fHistManager.FillTH1("fHistEmbeddedEventRejection", "VertexDist", 1);
      }
      return kFALSE;
    }
  }

  // Check for pt hard bin outliers
  if (fPythiaHeader && fMCRejectOutliers)
  {
    // Pythia jet / pT-hard > factor
    // This corresponds to "condition 1" in AliAnalysisTaskEmcal
    // NOTE: The other "conditions" defined there are not really suitable to define here, since they
    //       depend on the input objects of the event
    if (fPtHardJetPtRejectionFactor > 0.) {
      TLorentzVector jet;

      Int_t nTriggerJets =  fPythiaHeader->NTriggerJets();

      AliDebugStream(4) << "Pythia Njets: " << nTriggerJets << ", pT Hard: " << fPythiaPtHard << "\n";

      Float_t tmpjet[]={0,0,0,0};
      for (Int_t iJet = 0; iJet< nTriggerJets; iJet++) {
        fPythiaHeader->TriggerJet(iJet, tmpjet);

        jet.SetPxPyPzE(tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3]);

        AliDebugStream(5) << "Pythia jet " << iJet << ", pycell jet pT: " << jet.Pt() << "\n";

        //Compare jet pT and pt Hard
        if (jet.Pt() > fPtHardJetPtRejectionFactor * fPythiaPtHard) {
          AliDebugStream(3) << "Event rejected because of MC outlier removal. Pythia header jet with: pT Hard " << fPythiaPtHard << ", pycell jet pT " << jet.Pt() << ", rejection factor " << fPtHardJetPtRejectionFactor << "\n";
          fHistManager.FillTH1("fHistEmbeddedEventRejection", "MCOutlier", 1);
          return kFALSE;
        }
      }
    }
  }

  return kTRUE;
}

/**
 * Initialize the external event by creating an event and then reading the event info from the TChain.
 *
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
 * Performing run-independent initialization to setup embedding.
 *
 * If the setup is not successful here, it will be repeated in UserExec().
 */
void AliAnalysisTaskEmcalEmbeddingHelper::UserCreateOutputObjects()
{
  SetupEmbedding();

  // Reinitialize the YAML config after it was streamed so that it can be used properly.
  fYAMLConfig.Reinitialize();

  // Setup AliEventCuts
  if (fUseInternalEventSelection)
  {
    AliDebugStream(1) << "Configuring AliEventCuts for internal event selection.\n";
    // Handle manual cuts
    if (fUseManualInternalEventCuts) {
      fInternalEventCuts.SetManualMode();
      // Implement these cuts by retrieving the event cuts object and setting them manually.
    }

    // Trigger selection
    bool useEventCutsAutomaticTriggerSelection = false;
    bool res = fYAMLConfig.GetProperty(std::vector<std::string>({"internalEventSelection", "useEventCutsAutomaticTriggerSelection"}), useEventCutsAutomaticTriggerSelection, false);
    if (res && useEventCutsAutomaticTriggerSelection) {
      // Use the automatic selection. We will validate the trigger mask later because it
      // isn't yet set if we're using automatic mode.
      AliDebugStream(1) << "Using the automatic trigger selection from AliEventCuts.\n";
    }
    else {
      // Use the cuts selected by via YAML.
      std::bitset<32> triggerMask(fInternalEventTriggerMask);
      AliDebugStream(1) << "Using the trigger selection specified via YAML: " << triggerMask << ".\n";
      fInternalEventCuts.OverrideAutomaticTriggerSelection(fInternalEventTriggerMask);
    }
  }
  
  // Set up timer for logging purposes
  if (fPrintTimingInfoToLog) {
    fTimer = TStopwatch();
  }

  if (!fCreateHisto) {
    return;
  }

  // Create output list
  OpenFile(1);
  fOutput = new AliEmcalList();
  fOutput->SetOwner();

  // Get the histograms from AliEventCuts
  if (fUseInternalEventSelection) {
    // This list will be owned by fOutput, so it won't be leaked.
    TList * eventCutsOutput = new TList();
    eventCutsOutput->SetOwner(kTRUE);
    eventCutsOutput->SetName("EventCuts");

    // Add the event cuts to the output
    fInternalEventCuts.AddQAplotsToList(eventCutsOutput);
    fOutput->Add(eventCutsOutput);
  }

  // Create histograms
  TString histName;
  TString histTitle;

  // Cross section
  histName = "fHistXsection";
  histTitle = "Pythia Cross Section;p_{T} hard bin; XSection";
  fHistManager.CreateTProfile(histName, histTitle, fNPtHardBins, 0, fNPtHardBins);

  // Trials
  histName = "fHistTrials";
  histTitle = "Number of Pythia Trials;p_{T} hard bin;Trials";
  fHistManager.CreateTH1(histName, histTitle, fNPtHardBins, 0, fNPtHardBins);

  // Pt hard spectra
  histName = "fHistPtHard";
  histTitle = "p_{T} Hard Spectra;p_{T} hard;Counts";
  fHistManager.CreateTH1(histName, histTitle, 500, 0, 1000);

  // Count of accepted and rejected events
  histName = "fHistEventCount";
  histTitle = "Event count;Result;Count";
  auto histEventCount = fHistManager.CreateTH1(histName, histTitle, 2, 0, 2);
  histEventCount->GetXaxis()->SetBinLabel(1,"Accepted");
  histEventCount->GetXaxis()->SetBinLabel(2,"Rejected");

  // Event rejection reason
  histName = "fHistEmbeddedEventRejection";
  histTitle = "Reasons to reject embedded event";
  std::vector<std::string> binLabels = {"PhysSel", "MCOutlier", "Vz", "VertexDist", "PtHardIs0"};
  auto histEmbeddedEventRejection = fHistManager.CreateTH1(histName, histTitle, binLabels.size(), 0, binLabels.size());
  // Set label names
  for (unsigned int i = 1; i <= binLabels.size(); i++) {
    histEmbeddedEventRejection->GetXaxis()->SetBinLabel(i, binLabels.at(i-1).c_str());
  }
  histEmbeddedEventRejection->GetYaxis()->SetTitle("Counts");

  // Rejected events in embedded event selection
  histName = "fHistEmbeddedEventsAttempted";
  histTitle = "Number of embedded events rejected by event selection before success;Number of rejected events;Counts";
  fHistManager.CreateTH1(histName, histTitle, 200, 0, 200);

  // Number of files embedded
  histName = "fHistNumberOfFilesEmbedded";
  histTitle = "Number of files which contributed events to be embedded";
  fHistManager.CreateTH1(histName, histTitle, 1, 0, 2);

  // File number which was embedded
  histName = "fHistAbsoluteFileNumber";
  histTitle = "Number of times each absolute file number was embedded";
  fHistManager.CreateTH1(histName, histTitle, fMaxNumberOfFiles, 0, fMaxNumberOfFiles);

  if (fUseInternalEventSelection) {
    // Internal event cut statistics
    histName = "fHistInternalEventCutsStats";
    histTitle = "Number of events to pass each cut";
    binLabels = {"passedEventCuts", "centrality", "passedRandomRejection", "passedAllCuts"};
    auto histInternalEventCutsStats = fHistManager.CreateTH1(histName, histTitle, binLabels.size(), 0, binLabels.size());
    // Set label names
    for (unsigned int i = 1; i <= binLabels.size(); i++) {
      histInternalEventCutsStats->GetXaxis()->SetBinLabel(i, binLabels.at(i-1).c_str());
    }
    histInternalEventCutsStats->GetYaxis()->SetTitle("Number of selected events");
  }
  
  // Time to execute InitTree()
  if (fPrintTimingInfoToLog) {
    histName = "fInitTreeCPUtime";
    histTitle = "CPU time to execute InitTree() (s)";
    fHistManager.CreateTH1(histName, histTitle, 200, 0, 2000);
    
    histName = "fInitTreeRealtime";
    histTitle = "Real time to execute InitTree() (s)";
    fHistManager.CreateTH1(histName, histTitle, 200, 0, 2000);
  }

  // Add all histograms to output list
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }

  PostData(1, fOutput);
}

/**
 * Use the input files to setup the TChain. This involves accessing the files on AliEn and ensuring that
 * each file actually exists. If the TChain is successfully set up, then embedding is ready to proceed.
 *
 * @return kTRUE if successful in setting up the TChain.
 */
Bool_t AliAnalysisTaskEmcalEmbeddingHelper::SetupInputFiles()
{
  // Determine which file to start with
  DetermineFirstFileToEmbed();

  // Setup TChain
  fChain = new TChain(fTreeName);

  // Determine whether AliEn is needed
  for (auto filename : fFilenames)
  {
    if (filename.find("alien://") != std::string::npos) {
      ::ConnectToAliEn();
      // No point in continuing to search once we know that it is needed
      break;
    }
  }

  // Add files for TChain
  // See: https://stackoverflow.com/a/8533198
  bool wrapped = false;
  std::string fullPythiaXSecFilename = "";
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

    // Add to the Chain
    AliDebugStream(4) << "Adding file to the embedded input chain \"" << *filename << "\".\n";
    fChain->Add(filename->c_str());

    // Determine the full pythia cross section filename based on the previously determined filename
    // If we have determined that it doesn't exist in the initialization then we don't repeated attempt to open
    // the file (which will fail)
    if (fPythiaXSecFilename != "") {
      // Could check here again whether it exists here, but almost certainly unnecessary.
      // Further, we won't check to ensure that rapid, repeated file access on AliEn doesn't cause any problmes!
      fullPythiaXSecFilename = ConstructFullPythiaXSecFilename(*filename, fPythiaXSecFilename, false);

      AliDebugStream(4) << "Adding pythia cross section file \"" << fullPythiaXSecFilename << "\".\n";

      // They will automatically be ordered the same as the files to embed!
      fPythiaCrossSectionFilenames.push_back(fullPythiaXSecFilename);
    }
  }

  // Keep track of the total number of files in the TChain to ensure that we don't start repeating within the chain
  fMaxNumberOfFiles = fChain->GetListOfFiles()->GetEntries();

  if (fFilenames.size() > fMaxNumberOfFiles) {
    AliErrorStream() << "Number of input files (" << fFilenames.size() << ") is larger than the number of available files (" << fMaxNumberOfFiles << "). Something went wrong when adding some of those files to the TChain!\n";
  }

  // Setup input event
  Bool_t res = InitEvent();
  if (!res) return kFALSE;

  return kTRUE;
}

/**
 * Check if the file pythia base filename can be found in the folder or archive corresponding where
 * the external event input file is found.
 *
 * @param externalEventFilename Path to external event input file.
 * @param pythiaFilename Name of the pythia cross section file to try.
 * @param testIfExists If true, will check if the filename that it has determined actually exists.
 *
 * @return True if the file was found
 */
std::string AliAnalysisTaskEmcalEmbeddingHelper::ConstructFullPythiaXSecFilename(std::string externalEventFilename, const std::string & pythiaFilename, bool testIfExists) const
{
  std::string pythiaXSecFilename = "";

  // Remove "#*.root" if necessary
  if (externalEventFilename.find(".zip#") != std::string::npos) {
    std::size_t pos = externalEventFilename.find_last_of("#");
    externalEventFilename.erase(pos);
  }

  // Handle different file types
  if (externalEventFilename.find(".zip") != std::string::npos)
  {
    // Handle zip files
    pythiaXSecFilename = externalEventFilename;
    pythiaXSecFilename += "#";
    pythiaXSecFilename += pythiaFilename;

    // Check if the file is accessible
    if (testIfExists) {
      // Unfortunately, we cannot test for the existence of a file in an archive.
      // Instead, we have to tolerate TFile throwing an error (maximum of two).
      std::unique_ptr<TFile> fTemp(TFile::Open(pythiaXSecFilename.c_str(), "READ"));

      if (!fTemp) {
        AliDebugStream(4) << "File " << pythiaXSecFilename.c_str() << " does not exist!\n";
        pythiaXSecFilename = "";
      }
      else {
        AliDebugStream(4) << "Found pythia cross section file \"" << pythiaXSecFilename.c_str() << "\".\n";
      }
    }
  }
  else
  {
    // Handle normal root files
    pythiaXSecFilename = gSystem->DirName(externalEventFilename.c_str());
    pythiaXSecFilename += "/";
    pythiaXSecFilename += pythiaFilename;

    // Check if the file is accessible
    if (testIfExists) {
      if(::IsFileAccessible(pythiaXSecFilename)) {
        AliDebugStream(4) << "Found pythia cross section file \"" << pythiaXSecFilename.c_str() << "\".\n";
      }
      else {
        AliDebugStream(4) << "File " << pythiaXSecFilename.c_str() << " does not exist!\n";
        pythiaXSecFilename = "";
      }
    }
  }

  return pythiaXSecFilename;
}

/**
 * Setup embedding by retrieving the file list and then setting up the TChain based on that file list.
 *
 * Sets fInitializedEmbedding to kTRUE when the initialization is completed.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::SetupEmbedding()
{
  if (fInitializedConfiguration == false) {
    AliFatal("The configuration is not initialized. Check that Initialize() was called!");
  }

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
 * Initializes a new TTree within the TChain by determining the limits of the current TTree
 * within the TChain. By carefully keeping track of the first and lest entry of the current
 * tree in the opened file, we allow a random entry point into the event in the tree without
 * jumping between files, which would not perform well on the grid.
 *
 * This function sets up a random entry point into the file if requested.
 *
 * Sets fInitializedNewFile to kTRUE when the new entries in the tree have been initialized.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::InitTree()
{
  // Start the timer (for logging purposes)
  if (fPrintTimingInfoToLog) {
    fTimer.Start(kTRUE);
    std::cout << "InitTree() has started for file " << (fFilenameIndex + fFileNumber + 1) % fMaxNumberOfFiles << fChain->GetCurrentFile()->GetName() << "..." << std::endl;
  }
  
  // Load first entry of the (next) file so that we can query information about it
  // (it is inaccessible otherwise).
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
    TRandom3 rand(0);
    fOffset = TMath::Nint(rand.Rndm()*(fUpperEntry-fLowerEntry))-1;
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

  // Add to the count the number of files which were embedded
  fHistManager.FillTH1("fHistNumberOfFilesEmbedded", 1);
  fHistManager.FillTH1("fHistAbsoluteFileNumber", (fFileNumber + fFilenameIndex) % fMaxNumberOfFiles);

  // Check for pythia cross section and extract if possible
  // fFileNumber corresponds to the next file
  // If there are pythia filenames, the number of match the file number of the tree.
  // If we previously gave up on extracting then there should be no entires
  if (fPythiaCrossSectionFilenames.size() > 0) {
    // Need to check that fFileNumber is smaller than the size of the vector because we don't check if
    if (fFileNumber < fPythiaCrossSectionFilenames.size()) {
      bool success = PythiaInfoFromCrossSectionFile(fPythiaCrossSectionFilenames.at(fFileNumber));

      if (!success) {
        AliDebugStream(3) << "Failed to retrieve cross section from xsec file. Will still attempt to get the information from the header.\n";
      }
    }
    else {
      AliErrorStream() << "Attempted to read past the end of the pythia cross section filenames vector. File number: " << fFileNumber << ", vector size: " << fPythiaCrossSectionFilenames.size() << ".\nThis should only occur if we have run out of files to embed!\n";
    }
  }

  AliDebug(2, TString::Format("Will start embedding file %i beginning from entry %i (entry %i within the file). NOTE: This file number is not equal to the absolute file number in the file list!", fFileNumber, fCurrentEntry, fCurrentEntry - fLowerEntry));
  // NOTE: Cannot use this print message, as it is possible that fMaxNumberOfFiles != fFilenames.size() because
  //       invalid filenames may be included in the fFilenames count!
  //AliDebug(2, TString::Format("Will start embedding file %i as the %ith file beginning from entry %i.", (fFilenameIndex + fFileNumber) % fMaxNumberOfFiles, fFileNumber, fCurrentEntry));

  // (re)set whether we have wrapped the tree
  fWrappedAroundTree = false;

  // Note that the tree in the new file has been initialized
  fInitializedNewFile = kTRUE;
  
  // Stop timer (for logging purposes)
  if (fPrintTimingInfoToLog) {
    fTimer.Stop();
    std::cout << "InitTree() complete. CPU time: " << fTimer.CpuTime() << " (s). Real time: " << fTimer.RealTime() << " (s)." << std::endl;
    fHistManager.FillTH1("fInitTreeCPUtime", fTimer.CpuTime());
    fHistManager.FillTH1("fInitTreeRealtime", fTimer.RealTime());
  }

}

/**
 * Extract pythia information from a cross section file. Modified from AliAnalysisTaskEmcal::PythiaInfoFromFile().
 *
 * @param pythiaFileName Path to the pythia cross section file.
 *
 * @return True if the information has been successfully extracted.
 */
bool AliAnalysisTaskEmcalEmbeddingHelper::PythiaInfoFromCrossSectionFile(std::string pythiaFileName)
{
  std::unique_ptr<TFile> fxsec(TFile::Open(pythiaFileName.c_str()));

  if (fxsec)
  {
    int trials = 0;
    double crossSection = 0;
    double nEvents = 0;
    // Check if it's a tree
    TTree *xtree = dynamic_cast<TTree*>(fxsec->Get("Xsection"));
    if (xtree) {
      UInt_t ntrials  = 0;
      Double_t xsection  = 0;
      xtree->SetBranchAddress("xsection",&xsection);
      xtree->SetBranchAddress("ntrials",&ntrials);
      xtree->GetEntry(0);
      trials = ntrials;
      crossSection = xsection;
      // TODO: Test this on a file which has pyxsec.root!
      nEvents = 1.;
      AliFatal("Have no tested pyxsec.root files. Need to determine the proper way to get nevents!!");
    }
    else {
      // Check if it's instead the histograms
      // find the Tlist we want to be independent of the name so use the Tkey
      TKey* key = static_cast<TKey*>(fxsec->GetListOfKeys()->At(0));
      if (!key) return false;
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) return false;
      TProfile * crossSectionHist = static_cast<TProfile*>(list->FindObject("h1Xsec"));
      // check for failure
      if(!(crossSectionHist->GetEntries())) {
        // No cross section information available - fall back to raw
        AliErrorStream() << "No cross section information available in file \"" << fxsec->GetName() << "\". Will still attempt to extract cross section information from pythia header.\n";
      } else {
        // Cross section histogram filled - take it from there
        crossSection = crossSectionHist->GetBinContent(1);
        if(!crossSection) AliErrorStream() << GetName() << ": Cross section 0 for file " << pythiaFileName << std::endl;
      }
      TH1F * trialsHist = static_cast<TH1F*>(list->FindObject("h1Trials"));
      trials = trialsHist->GetBinContent(1);
      nEvents = trialsHist->GetEntries();
    }

    // If successful in retrieving the values, normalize the xsec and trials by the number of events
    // in the file. This way, we can use it as an approximate event-by-event value
    // We do not want to just use the overall value because some of the events may be rejected by various
    // event selections, so we only want that ones that were actually use. The easiest way to do so is by
    // filling it for each event.
    fPythiaTrialsFromFile = trials/nEvents;
    // Do __NOT__ divide by nEvents here! The value is already from a TProfile and therefore is already the mean!
    fPythiaCrossSectionFromFile = crossSection;

    return true;
  }
  else {
    AliDebugStream(3) << "Unable to open file \"" << pythiaFileName << "\". Will attempt to use values from the header.";
  }

  // Could not open file
  return false;
}

/**
 * Run the main analysis code here. If for some reason the embedding was not successfully set up
 * in UserCreateOutputObjects(), it is set up against before continuing. It also ensures that the
 * current tree is available before attempting to get the next entry.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::UserExec(Option_t*)
{
  if (!fInitializedEmbedding) {
    AliError("Chain not initialized before running! Setting up now.");
    SetupEmbedding();
  }

  // Apply internal event selection
  if (fUseInternalEventSelection) {
    fEmbeddedEventUsed = false;
    if (fInternalEventCuts.AcceptEvent(AliAnalysisTaskSE::InputEvent()) == true)
    {
      fEmbeddedEventUsed = true;
      fHistManager.FillTH1("fHistInternalEventCutsStats", "passedEventCuts", 1);

      // Validate the event selection. We will only do this once, but we must wait
      // until `AcceptEvent(...)` is called once in case automatic setup is used.
      if (fValidatedPhysicsSelection == false) {
        // Validate that the trigger selection in the other tasks is a subset of the internal event selection
        AliDebugStream(1) << "Validating physics selection.\n";
        ValidatePhysicsSelectionForInternalEventSelection();
        fValidatedPhysicsSelection = true;
        AliDebugStream(1) << "Successfully validated the physics selection!\n";
      }

      // The event was accepted by AliEventCuts. Now check for additional cuts.
      // Centrality
      // NOTE: If the centrality range is the same as AliEventCuts, then simply all will pass
      //       If a wider centrality range than in AliEventCuts is needed then it must be _entirely_
      //       configured through manual mode.
      if (fCentMin != -999 && fCentMax != -999) {
        if (fInternalEventCuts.GetCentrality() < fCentMin || fInternalEventCuts.GetCentrality() > fCentMax) {
          fEmbeddedEventUsed = false;
        }
        else {
          fHistManager.FillTH1("fHistInternalEventCutsStats", "centrality", 1);
        }
      }
      // Now reject events based on a rejection factor if set, where the fraction of events 
      // kept is equal to 1 / fRandomRejectionFactor
      if(fRandomRejectionFactor > 1.) {
        Double_t rand = fRandom.Rndm();
        if(fRandomRejectionFactor > 1./rand) {
          fEmbeddedEventUsed = false;
        }
        else {
          fHistManager.FillTH1("fHistInternalEventCutsStats", "passedRandomRejection", 1);
        }
      }
      if (fEmbeddedEventUsed) {
        // Record all cuts passed
        fHistManager.FillTH1("fHistInternalEventCutsStats", "passedAllCuts", 1);
      }
    }

    // If the internal event was rejected, then record and skip this internal event by asking the
    // analysis manager to break the execution.
    if (fEmbeddedEventUsed == false) {
      if (fCreateHisto) {
        PostData(1, fOutput);
      }
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      if (!mgr) {
        AliFatal("No analysis manager to connect to.");
      }
      // Ask the analysis manager to break execution, which will prevent any downstream tasks
      // from executing.
      AliDebugStream(3) << "Internal event rejected due to event selection. Breaking execution early.\n";
      mgr->BreakExecutionChain(kTRUE);
      return;
    }
  }

  if (!fInitializedNewFile) {
    InitTree();
  }

  Bool_t res = GetNextEntry();

  if (!res) {
    AliError("Unable to get the event to embed. Nothing will be embedded.");
    return;
  }

  if (fCreateHisto && fOutput) {
    PostData(1, fOutput);
  }
}

/**
 * Validate that the physics selection of the embedding helper is at least as broad
 * as all other tasks. Also checks that the internal event selection isn't broader
 * than the main physics selection. If either is not the case, it will throw an error.
 *
 * This ensures that other tasks are never accessing an old embedded event due to other
 * tasks seeing events that the embedding helper is not selected for.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::ValidatePhysicsSelectionForInternalEventSelection()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliFatal("No analysis manager to connect to.");
  }

  // Check that the internal event physics selection that has been applied is a subset of
  // the collision candidates that are selected. Otherwise, it could be quite misleading.
  // NOTE: We either set `fTriggerMask` to our internal event selection, or it was set automatically.
  //       In either case, `fInternalEventCuts.fTriggerMask` should have the correct internal
  //       event physics selection.
  // For information on the comparison method, see: https://stackoverflow.com/a/8639510
  UInt_t collisionCandidates = this->GetCollisionCandidates();
  bool res = (fInternalEventCuts.fTriggerMask | collisionCandidates) == collisionCandidates;
  if (res == false) {
    std::stringstream message;
    message << "Collision candidates selected for the embedding helper are more restrictive than"
        << " the internal event physics selection! You will not have access to all of the events"
        << " selected in the internal event physics selection. Please expand your trigger mask set"
        << " via SelectCollisionCandidates().\n"
        << std::bitset<32>(fInternalEventCuts.fTriggerMask) << " <- Embedding helper internal event physics selection\n"
        << std::bitset<32>(collisionCandidates) << " <- Collision candidates\n";
    AliFatal(message.str().c_str());
  }

  auto tasks = mgr->GetTasks();
  for (auto t : *tasks)
  {
    auto task = dynamic_cast<AliAnalysisTaskSE *>(t);
    if (!task || task->GetName() == GetName()) {
      // Skip the task if it's not an analysis task or if it's the embedding
      // helper, since it's allowed to have a broader physics selection than the
      // internal event physics selection.
      continue;
    }

    // Compare the selected collision candidates to the embedding helper physics selection.
    // Every subsequent task must be a subset or equal to the embedding helper physics selection.
    // Otherwise, subsequent tasks will be be selected for internal events where we haven't update
    // the embedded event. This will lead to double counting, double corrections, etc.
    // For information on the comparison method, see: https://stackoverflow.com/a/8639510
    UInt_t taskCollisionCandidates = task->GetCollisionCandidates();
    res = (taskCollisionCandidates | collisionCandidates) == collisionCandidates;
    if (res == false) {
      std::stringstream message;
      message << "The physics selection of all tasks must be a subset of the physics selection used"
          << " in the embedding helper.\n"
          << std::bitset<32>(collisionCandidates) << " <- Embedding helper internal event physics selection\n"
          << std::bitset<32>(taskCollisionCandidates) << " <- Task \"" << task->GetName() << "\" collision candidates\n";
      AliFatal(message.str().c_str());
    }
  }
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::Terminate(Option_t*)
{
}

/**
 * Remove the dummy task which had to be added by ConfigureEmcalEmbeddingHelperOnLEGOTrain()
 * from the Analysis Manager. This is the same function as in AliEmcalCorrectionTask.
 */
void AliAnalysisTaskEmcalEmbeddingHelper::RemoveDummyTask() const
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    AliErrorStream() << "No analysis manager to connect to.\n";
    return;
  }

  // Remove the dummy task
  std::string dummyTaskName = GetName();
  dummyTaskName += "_dummyTask";
  TObjArray * tasks = mgr->GetTasks();
  if (tasks) {
    AliAnalysisTaskSE * dummyTask = dynamic_cast<AliAnalysisTaskSE *>(tasks->FindObject(dummyTaskName.c_str()));
    if (!dummyTask) {
      AliErrorStream() << "Could not remove dummy task \"" << dummyTaskName << "\" from analysis manager! Was it added?\n";
    }
    // Actually remove the task
    tasks->Remove(dummyTask);
    AliDebugStream(1) << "Removed dummy task named \"" << dummyTaskName << "\".\n";
  }
  else {
    AliErrorStream() << "Could not retrieve tasks from the analysis manager.\n";
  }
}

AliEventCuts * AliAnalysisTaskEmcalEmbeddingHelper::GetInternalEventCuts()
{
  if (fUseManualInternalEventCuts) {
    return &fInternalEventCuts;
  }
  else {
    AliErrorStream() << "Enable manual mode in AliEventCuts (though the embedding helper) to access this object.\n";
  }
  return nullptr;
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

  TString outputContainerName(name);
  outputContainerName += "_histos";

  AliAnalysisDataContainer * cOutput = mgr->CreateContainer(outputContainerName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(embeddingHelper, 0, cInput);
  mgr->ConnectOutput(embeddingHelper, 1, cOutput);

  return embeddingHelper;
}

AliAnalysisTaskEmcalEmbeddingHelper * AliAnalysisTaskEmcalEmbeddingHelper::ConfigureEmcalEmbeddingHelperOnLEGOTrain()
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("ConfigureEmcalEmbeddingHelperOnLEGOTrain", "No analysis manager to connect to.");
    return nullptr;
  }

  // Retrieve the embedding helper
  auto embeddingHelperConst = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  // Cast away const-ness on the pointer since the underlying object is not const and we need to be able to modify it.
  auto embeddingHelper = const_cast<AliAnalysisTaskEmcalEmbeddingHelper *>(embeddingHelperConst);

  // Fatal if we can't find the task
  if (!embeddingHelper) {
    AliFatalClass("Could not find embedding helper, Did you remember to create it?");
  }

  AliInfoClassStream() << "Found embedding helper to configure.\n";

  // AliAnalysisTaskCfg will require a task to be returned, so we add a dummy task to the analysis manager,
  // which will be removed when the user calls Initialize(true) on the embedding helper.
  mgr->AddTask(new AliAnalysisTaskSE("AliAnalysisTaskEmcalEmbeddingHelper_dummyTask"));

  return embeddingHelper;
}

/**
 * Prints information about the embedding helper.
 *
 * @return std::string containing information about the task.
 */
std::string AliAnalysisTaskEmcalEmbeddingHelper::toString(bool includeFileList) const
{
  std::stringstream tempSS;

  // General embedding helper information
  tempSS << std::boolalpha;
  tempSS << GetName() << ": Embedding helper configuration:\n";
  tempSS << "Create histos: " << fCreateHisto << "\n";
  tempSS << "Pt Hard Bin: " << fPtHardBin << "\n";
  tempSS << "N Pt Hard Bins: " << fNPtHardBins << "\n";
  tempSS << "File pattern: \"" << fFilePattern << "\"\n";
  tempSS << "Input filename: \"" << fInputFilename << "\"\n";
  tempSS << "Pythia cross section filename: \"" << fPythiaXSecFilename << "\"\n";
  tempSS << "File list filename: \"" << fFileListFilename << "\"\n";
  tempSS << "Tree name: " << fTreeName << "\n";
  tempSS << "Print timing info to log: " << fPrintTimingInfoToLog << "\n";
  tempSS << "Random event number access: " << fRandomEventNumberAccess << "\n";
  tempSS << "Random file access: " << fRandomFileAccess << "\n";
  tempSS << "Starting file index: " << fFilenameIndex << "\n";
  tempSS << "Number of files to embed: " << fFilenames.size() << "\n";
  tempSS << "YAML configuration path: \"" << fConfigurationPath << "\"\n";
  tempSS << "Enable internal event selection: " << fUseInternalEventSelection << "\n";
  tempSS << "Use manual event cuts mode for internal event selection: " << fUseManualInternalEventCuts << "\n";
  if (fCentMin != -999 && fCentMax != -999) {
    tempSS << "Internal event selection centrality range: [" << fCentMin << ", " << fCentMax << "]\n";
  }
  else {
    tempSS << "Internal event selection centrality range disabled.\n";
  }
  tempSS << "Internal event physics selection via class (should propagate to AliEventCuts): " << std::bitset<32>(fInternalEventTriggerMask) << "\n";
  tempSS << "Internal event physics selection via AliEventCuts: " << std::bitset<32>(fInternalEventCuts.fTriggerMask) << "\n";

  std::bitset<32> triggerMask(fTriggerMask);
  tempSS << "\nEmbedded event settings:\n";
  tempSS << "Trigger mask (binary): " << triggerMask << "\n";
  tempSS << "Reject outliers: " << fMCRejectOutliers << "\n";
  tempSS << "Pt hard jet pt rejection factor: " << fPtHardJetPtRejectionFactor << "\n";
  tempSS << "Z vertex cut: " << fZVertexCut << "\n";
  tempSS << "Max difference between internal and embedded vertex: " << fMaxVertexDist << "\n";
  tempSS << "Random event rejection factor: " << fRandomRejectionFactor << "\n";

  if (includeFileList) {
    tempSS << "\nFiles to embed:\n";
    for (auto filename : fFilenames) {
      tempSS << "\t" << filename << "\n";
    }
  }

  return tempSS.str();
}

/**
 * Print embedding helper information on an output stream using the string representation provided by
 * AliAnalysisTaskEmcalEmbeddingHelper::toString(). Used by operator<<
 *
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream & AliAnalysisTaskEmcalEmbeddingHelper::Print(std::ostream & in) const {
  in << toString();
  return in;
}

/**
 * Implementation of the output stream operator for AliAnalysisTaskEmcalEmbeddingHelper. Printing
 * basic embedding helper information provided by function toString()
 *
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream & operator<<(std::ostream & in, const AliAnalysisTaskEmcalEmbeddingHelper & myTask)
{
  std::ostream & result = myTask.Print(in);
  return result;
}

/**
 * Print basic embedding helper information using the string representation provided by
 * AliAnalysisTaskEmcalEmbeddingHelper::toString()
 *
 * @param opt If "FILELIST" is passed, then the list of files to embed is also printed
 */
void AliAnalysisTaskEmcalEmbeddingHelper::Print(Option_t* opt) const
{
  std::string temp(opt);
  bool includeFileList = false;
  if (temp == "FILELIST") {
    includeFileList = true;
  }
  Printf("%s", toString(includeFileList).c_str());
}
