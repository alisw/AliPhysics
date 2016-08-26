// AliEmcalCorrectionTask
//
// Runs the cell and cluster level corrections for the EMCal
//

#include "AliEmcalCorrectionTask.h"
#include "AliEmcalCorrectionComponent.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include <yaml-cpp/yaml.h>

#include <TChain.h>
#include <TSystem.h>
#include <TGrid.h>

#include "AliVEventHandler.h"
#include "AliEMCALGeometry.h"
#include "AliVCaloCells.h"
#include "AliAODCaloCells.h"
#include "AliESDCaloCells.h"
#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMultSelection.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionTask);
/// \endcond

/**
 *
 *
 */
AliEmcalCorrectionTask::AliEmcalCorrectionTask() :
  AliAnalysisTaskSE("AliEmcalCorrectionTask"),
  fUserConfiguration(),
  fUserConfigurationFilename(""),
  fUserConfigurationString(""),
  fDefaultConfiguration(),
  fDefaultConfigurationFilename(""),
  fDefaultConfigurationString(""),
  fCorrectionComponents(),
  fIsEsd(false),
  fForceBeamType(kNA),
  fRunPeriod(""),
  fConfigurationInitialized(false),

  fCreateNewObjectBranches(false),
  fCreatedClusterBranchName(""),
  fCreatedTrackBranchName(""),
  fEventInitialized(false),
  fCent(0),
  fCentBin(-1),
  fMinCent(-999),
  fMaxCent(-999),
  fNcentBins(4),
  fCentEst("V0M"),
  fUseNewCentralityEstimation(kFALSE),
  fNVertCont(0),
  fBeamType(kNA),
  fCaloCellsName(),
  fNeedEmcalGeom(kTRUE),
  fGeom(0),
  fParticleCollArray(),
  fClusterCollArray(),
  fCaloCells(0),
  fCaloCellsFromInputEvent(0),
  fOutput(0)
{
  // Default constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);
}

/**
 *
 *
 */
AliEmcalCorrectionTask::AliEmcalCorrectionTask(const char * name) :
  AliAnalysisTaskSE(name),
  fUserConfiguration(),
  fUserConfigurationFilename(""),
  fUserConfigurationString(""),
  fDefaultConfiguration(),
  fDefaultConfigurationFilename(""),
  fDefaultConfigurationString(""),
  fCorrectionComponents(),
  fIsEsd(false),
  fForceBeamType(kNA),
  fRunPeriod(""),
  fConfigurationInitialized(false),

  fCreateNewObjectBranches(false),
  fCreatedClusterBranchName(""),
  fCreatedTrackBranchName(""),
  fEventInitialized(false),
  fCent(0),
  fCentBin(-1),
  fMinCent(-999),
  fMaxCent(-999),
  fNcentBins(4),
  fCentEst("V0M"),
  fUseNewCentralityEstimation(kFALSE),
  fNVertCont(0),
  fBeamType(kNA),
  fCaloCellsName(),
  fNeedEmcalGeom(kTRUE),
  fGeom(0),
  fParticleCollArray(),
  fClusterCollArray(),
  fCaloCells(0),
  fCaloCellsFromInputEvent(0),
  fOutput(0)
{
  // Standard constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

/** Checks if a file exists. This is done inline to make it efficient.
 * See: https://stackoverflow.com/a/19841704
 *
 * \param filename String containing the filename of the file to check.
 *
 * \return True if the file exists.
 */
inline bool AliEmcalCorrectionTask::doesFileExist(const std::string & filename)
{
  std::ifstream inFile(filename);
  return inFile.good();
}

/**
 * Handles expanding ALICE_PHYSICS and copying the file from the grid if necessary.
 *
 */
void AliEmcalCorrectionTask::SetupConfigurationFilePath(std::string & filename, bool userFile)
{
  if (filename != "")
  {
    // Handle if in AliPhysics
    // Check for and replace $ALICE_PHYSICS with the actual path if needed
    std::size_t alicePhysicsPathLocation = filename.find("$ALICE_PHYSICS");
    if (alicePhysicsPathLocation != std::string::npos)
    {
      TString alicePhysicsPath = gSystem->Getenv("ALICE_PHYSICS");
      // "$ALICE_PHYSICS "is 14 characters
      filename.replace(alicePhysicsPathLocation, alicePhysicsPathLocation + 14, alicePhysicsPath.Data());
    }

    // Handle grid
    if(filename.find("alien://") != std::string::npos)
    {
      AliDebug(2, TString::Format("Opening file \"%s\" on the grid!", filename.c_str()));
      // Init alien connection if needed
      if (!gGrid) {
        TGrid::Connect("alien://");
      }

      // Determine the loca filename and copy file to local directory
      std::string localFilename = gSystem->BaseName(filename.c_str());
      // Ensures that the default and user files do not conflict if both are taken from the grid and have the same filename
      if (userFile == true) {
        localFilename = "user" + localFilename;
      }
      TFile::Cp(filename.c_str(), localFilename.c_str());

      // yaml-cpp should only open the local file
      filename = localFilename;
    }
  }
}


/**
 *
 *
 */
void AliEmcalCorrectionTask::InitializeConfiguration()
{
  // Determine file path
  if (fDefaultConfigurationFilename == "")
  {
    // Use the default if nothing is set
    fDefaultConfigurationFilename = "$ALICE_PHYSICS/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml";
  }

  // Setup the YAML files
  // Default
  SetupConfigurationFilePath(fDefaultConfigurationFilename);

  if (doesFileExist(fDefaultConfigurationFilename) == true)
  {
    AliInfo(TString::Format("Using default EMCal corrections configuration located at %s", fDefaultConfigurationFilename.c_str()));

    fDefaultConfiguration = YAML::LoadFile(fDefaultConfigurationFilename);
    // Check for valid file
    if (fDefaultConfiguration.IsNull() == true)
    {
      AliFatal(TString::Format("Could not open the default configuration file \"%s\"!", fDefaultConfigurationFilename.c_str()));
    }
  }
  else
  {
    AliFatal(TString::Format("Default file located at \"%s\" does not exist!", fDefaultConfigurationFilename.c_str()));
  }

  // User
  SetupConfigurationFilePath(fUserConfigurationFilename, true);

  if (doesFileExist(fUserConfigurationFilename) == true)
  {
    AliInfo(TString::Format("Using user EMCal corrections configuration located at %s", fUserConfigurationFilename.c_str()));

    fUserConfiguration = YAML::LoadFile(fUserConfigurationFilename);
  }
  else
  {
    AliInfo(TString::Format("User file at \"%s\" does not exist! All settings will be from the default file!", fUserConfigurationFilename.c_str()));
  }

  // Ensure that there is a run period
  if (fRunPeriod == "")
  {
    AliFatal("Must pass a run period to the correction task!");
  }
  // Check the user provided run period
  TString userRunPeriod = "kNoUserFile";
  // Test if the period exists in the user file
  if (fUserConfiguration.IsNull() != true)
  {
    if (fUserConfiguration["period"])
    {
      userRunPeriod = fUserConfiguration["period"].as<std::string>();
    }
    else
    {
      AliFatal("User must specify a period. Leave the period as an empty string to apply to all periods.");
    }
  }
  
  AliDebug(3, TString::Format("userRunPeriod: %s", userRunPeriod.Data()));
  // Normalize the user run period to lower case to ensure that we don't miss any matches
  userRunPeriod.ToLower();
  // "" means the user wants their settings to apply to all periods
  if (userRunPeriod != "" && userRunPeriod != "knouserfile" && userRunPeriod != fRunPeriod)
  {
    AliFatal(TString::Format("User run period \"%s\" does not match the run period of \"%s\" passed to the correction task!", userRunPeriod.Data(), fRunPeriod.Data()));
  }

  // Ensure that the user is aware
  if (userRunPeriod == "")
  {
    AliWarning("User run period is an empty string. Settings apply to all run periods!");
  }

  // Save configuration into strings so that they can be streamed
  // Need the stringstream because yaml implements streamers
  std::stringstream tempConfiguration;
  tempConfiguration << fUserConfiguration;
  fUserConfigurationString = tempConfiguration.str();
  tempConfiguration.str("");
  tempConfiguration << fDefaultConfiguration;
  fDefaultConfigurationString = tempConfiguration.str();

  // Note that it is initialized properly so that the analysis can proceed
  fConfigurationInitialized = true;
}

/**
 * Writes the desired yaml configuration to a file.
 * 
 * \param filename The name of the file to write.
 * \param userCofig True to write the user configuration.
 * \return Whether writing the configuration to the file was successful.
 */
bool AliEmcalCorrectionTask::WriteConfigurationFile(std::string filename, bool userConfig)
{
  bool returnValue = false;
  if (filename != "")
  {
    if (fConfigurationInitialized == true)
    {
      std::ofstream outFile(filename);
      std::string stringToWrite = userConfig ? fUserConfigurationString : fDefaultConfigurationString;
      if (stringToWrite == "") {
        AliWarning(TString::Format("%s configuration is empty!", userConfig ? "User" : "Default"));
      }
      outFile << stringToWrite;
      outFile.close();

      returnValue = true;
    }
    else
    {
      AliWarning(TString::Format("Configuration not properly initialized! Cnanot print %s configuration!", userConfig ? "user" : "default"));
    }

  }
  else
  {
    AliWarning("Please pass a valid filename instead of empty qutoes!");
  }
  return returnValue;
}


/**
 *
 *
 */
void AliEmcalCorrectionTask::RetrieveExecutionOrder(std::vector <std::string> & executionOrder)
{
  AliEmcalCorrectionComponent::GetProperty("executionOrder", executionOrder, fUserConfiguration, fDefaultConfiguration, true);
  // Need to append "AliEmcalCorrection" to allow the tasks to be found!
  AliDebug(2, "Creating EMCal Correction Components: ");
  for (auto & component : executionOrder)
  {
    component = "AliEmcalCorrection" + component;
    AliDebug(2, TString::Format("%s", component.c_str()) );
  }
}

/**
 *
 *
 */
void AliEmcalCorrectionTask::InitializeComponents()
{
  // YAML Objects cannot be streamed, so we need to reinitialize them here.
  // They need reinitialize if they are null
  if (fUserConfiguration.IsNull() == true && fUserConfigurationString != "")
  {
    AliInfo(TString::Format("%s: Reinitializing user configuration from string. Expected if running on grid!", GetName()));
    fUserConfiguration = YAML::Load(fUserConfigurationString);
  }
  if (fDefaultConfiguration.IsNull() == true)
  {
    AliInfo(TString::Format("%s: Reinitializing default configuration from string. Expected if running on grid!", GetName()));
    fDefaultConfiguration = YAML::Load(fDefaultConfigurationString);
  }

  // Create a function to handle creation and configuration of the all of the created module
  std::vector<std::string> executionOrder;
  RetrieveExecutionOrder(executionOrder);
 
  // Allow for output files
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  
  // Iterate over the list
  AliEmcalCorrectionComponent * component = 0;
  bool componentEnabled = false;
  for (auto componentName : executionOrder)
  {
    componentEnabled = false;
    AliEmcalCorrectionComponent::GetProperty("enabled", componentEnabled, fUserConfiguration, fDefaultConfiguration, true, componentName);
    if (componentEnabled == false)
    {
      AliInfo(TString::Format("%s: Component %s is disabled and will not be run!", GetName(), componentName.c_str()));
      continue;
    }

    //Printf("Attempting to add task: %s", tempString->GetString().Data());
    component = AliEmcalCorrectionComponentFactory::createInstance(componentName);
    if (!component)
    {
      AliFatal(TString::Format("%s: Failed to create requested component %s!", GetName(), componentName.c_str()));
    }

    // For setting names of tasks to differentiate between tasks of the same class
    component->SetName(componentName.c_str());
    component->SetTitle(componentName.c_str());

    // Initialize the YAML configurations in each component
    component->SetUserConfiguration(fUserConfiguration);
    component->SetDefaultConfiguration(fDefaultConfiguration);
    
    // configure needed fields for components to properly initialize
    component->SetNcentralityBins(fNcentBins);
    
    // Initialize each component
    component->Initialize();

    // Attempt to retrieve any containers they may have created in initialization
    AliClusterContainer * tempClusterContainer = component->GetClusterContainer();
    if (tempClusterContainer)
    {
      AdoptClusterContainer(tempClusterContainer);
    }

    AliParticleContainer * tempParticleContainer = component->GetParticleContainer();
    if (tempParticleContainer)
    {
      AdoptParticleContainer(tempParticleContainer);
    }

    if (component)
    {
      AliInfo(TString::Format("Successfully added correction task: %s", componentName.c_str()));
      fCorrectionComponents.push_back(component);
    }

    if (component->GetOutputList() != 0)
    {
      // Adds a list to the list -- this doesn't work for some unknown reason
      //fOutput->Add(component->GetOutputList());
      
      // iterate through lists for each component, and fill in output
      TList* t = new TList();
      t->SetName(componentName.c_str());
      fOutput->AddLast(t);
      t = (TList*)fOutput->Last();
      TIter next(component->GetOutputList());
      while (TObject *obj = next()){
        t->Add(obj);
      }
      
      AliDebug(1, TString::Format("Added output list from task %s to output.", componentName.c_str()));
    }
  }

  PostData(1, fOutput);

  // Just printing in order is probably better
  // Available components
  /*AliEmcalCorrectionComponentFactory::map_type * availableComponents = AliEmcalCorrectionComponentFactory::getMap();

  // Need to sort by execution order
  Printf("Correction Components:");
  bool availableFlag = false;
  for (AliEmcalCorrectionComponentFactory::map_type::const_iterator it = availableComponents->begin(); it != availableComponents->end(); ++it)
  {
    if (fCorrectionComponents.FindObject(it->first.c_str()))
    {
      availableFlag = true;
    }

    Printf("%s:\t(%s)", it->first.c_str(), availableFlag ? "Enabled" : "Disabled");

    availableFlag = false;
  }*/
}

/**
 *
 *
 */
void AliEmcalCorrectionTask::UserCreateOutputObjects()
{
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));

  // Determine file type
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) {
    AliVEventHandler *evhand = mgr->GetInputEventHandler();
    if (evhand) {
      if (evhand->InheritsFrom("AliESDInputHandler")) {
        fIsEsd = true;
      }
      else {
        fIsEsd = false;
      }
    }
    else {
      AliError("Event handler not found!");
    }
  }
  else {
    AliError("Analysis manager not found!");
  }

  // Initialize the components
  //InitializeConfiguration();
  if (fConfigurationInitialized != true)
  {
    AliFatal(TString::Format("%s: YAML configuration must be initialized before running (ie. the AddTask, run macro or wagon)!", GetName()));
  }

  // Retrieve cells name from configruation
  std::string cellsName = "";
  AliEmcalCorrectionComponent::GetProperty("cellBranchName", cellsName, fUserConfiguration, fDefaultConfiguration, true, "");
  // In the case of fCreateNewObjectBranches, we need to get the default name to retrieve the normal cells
  // We will then retrieve the new cells name later
  if (cellsName == "usedefault" || fCreateNewObjectBranches) {
    cellsName = AliEmcalCorrectionComponent::DetermineUseDefaultName(AliEmcalCorrectionComponent::kCaloCells, fIsEsd);
  }
  fCaloCellsName = cellsName;

  if (fForceBeamType == kpp)
    fNcentBins = 1;

  InitializeComponents();
}

/**
 *
 *
 */
void AliEmcalCorrectionTask::ExecOnceComponents()
{
  // Run the initialization for all derived classes.
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  for (auto component : fCorrectionComponents)
  {
    component->ExecOnce();
  }
}

/**
 *
 *
 */
void AliEmcalCorrectionTask::CreateNewObjectBranches()
{
  // Create new cell branch
  // cellBranchName doesn't belong to any particular component
  AliEmcalCorrectionComponent::GetProperty("cellBranchName", fCreatedCellBranchName, fUserConfiguration, fDefaultConfiguration, true, "");
  // While it is wrong to use "usedefault" here, this will provide a more meaningful message error message
  if (fCreatedCellBranchName == "usedefault") {
    fCreatedCellBranchName = AliEmcalCorrectionComponent::DetermineUseDefaultName(AliEmcalCorrectionComponent::kCaloCells, fIsEsd);
  }

  // Check to ensure that we are not trying to create a new branch on top of the old branch
  TObject * existingObject = InputEvent()->FindListObject(fCreatedCellBranchName.c_str());
  if (existingObject) {
    AliFatal(TString::Format("Attempted to create a new cell branch, \"%s\", with the same name as an existing branch! Check your configuration! Perhaps \"usedefault\" was used incorrectly?", fCreatedCellBranchName.c_str()));
  }

  // Create new branch and add it to the event
  AliVCaloCells * newCells = 0;
  if (fIsEsd) {
    newCells = new AliESDCaloCells(fCreatedCellBranchName.c_str(), fCreatedCellBranchName.c_str(), AliVCaloCells::kEMCALCell);
  }
  else {
    newCells = new AliAODCaloCells(fCreatedCellBranchName.c_str(), fCreatedCellBranchName.c_str(), AliVCaloCells::kEMCALCell);
  }
  InputEvent()->AddObject(newCells);
  // Set fCaloCells here to avoid looking up every event
  // newCells is the same address as that reference by the TList in the InputEvent
  fCaloCells = newCells;

  // Create new cluster branch
  // Clusterizer is used since it is the first place that clusters are used.
  AliEmcalCorrectionComponent::GetProperty("clusterBranchName", fCreatedClusterBranchName, fUserConfiguration, fDefaultConfiguration, true, "Clusterizer");
  // While it is wrong to use "usedefault" here, this will provide a more meaningful message error message
  if (fCreatedClusterBranchName == "usedefault") {
    fCreatedClusterBranchName = AliEmcalCorrectionComponent::DetermineUseDefaultName(AliEmcalCorrectionComponent::kCluster, fIsEsd);
  }

  // Check to ensure that we are not trying to create a new branch on top of the old branch
  existingObject = InputEvent()->FindListObject(fCreatedClusterBranchName.c_str());
  if (existingObject) {
    AliFatal(TString::Format("Attempted to create a new cluster branch, \"%s\", with the same name as an existing branch! Check your configuration! Perhaps \"usedefault\" was used incorrectly?", fCreatedClusterBranchName.c_str()));
  }
  TClonesArray * existingArray = dynamic_cast<TClonesArray *>(existingObject);

  // Create new branch and add it to the event
  TClonesArray * newClusters = 0;
  if (fIsEsd) {
    newClusters = new TClonesArray("AliESDCaloCluster");
  }
  else {
    newClusters = new TClonesArray("AliAODCaloCluster");
  }
  newClusters->SetName(fCreatedClusterBranchName.c_str());
  InputEvent()->AddObject(newClusters);

  // Create new tracks branch
  // ClusterTrackMatcher is used since it is the first place that tracks are used.
  AliEmcalCorrectionComponent::GetProperty("trackBranchName", fCreatedTrackBranchName, fUserConfiguration, fDefaultConfiguration, true, "ClusterTrackMatcher");
  // While it is wrong to use "usedefault" here, this will provide a more meaningful message error message
  if (fCreatedTrackBranchName == "usedefault") {
    fCreatedTrackBranchName = AliEmcalCorrectionComponent::DetermineUseDefaultName(AliEmcalCorrectionComponent::kTrack, fIsEsd);
  }

  // Check to ensure that we are not trying to create a new branch on top of the old branch
  existingObject = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(fCreatedTrackBranchName.c_str()));
  if (existingObject) {
    AliFatal(TString::Format("Attempted to create a new track branch, \"%s\", with the same name as existing branch! Check your configuration! Perhaps \"usedefault\" was used incorrectly?", fCreatedTrackBranchName.c_str()));
  }
  existingArray = dynamic_cast<TClonesArray *>(existingObject);

  // Create new branch and add it to the event
  TClonesArray * newTracks = 0;
  if (fIsEsd) {
    newTracks = new TClonesArray("AliESDtrack");
  }
  else {
    newTracks = new TClonesArray("AliAODTrack");
  }
  newTracks->SetName(fCreatedTrackBranchName.c_str());
  InputEvent()->AddObject(newTracks);
}

/**
 *
 *
 */
void AliEmcalCorrectionTask::CopyBranchesToNewObjects()
{
  // Cells
  AliDebug(3, Form("Number of old cells: %d", fCaloCellsFromInputEvent->GetNumberOfCells()));
  // The function CopyCaloCells() does not work, so we use the assignment operator instead!
  if (fIsEsd)
  {
    AliESDCaloCells * currentCells = dynamic_cast<AliESDCaloCells *>(fCaloCellsFromInputEvent);
    AliESDCaloCells * newCells = dynamic_cast<AliESDCaloCells *>(fCaloCells);
    fCaloCells = dynamic_cast<AliESDCaloCells *>(new (newCells) AliESDCaloCells(*currentCells));
    // The name was changed to the currentCells name, but we want it to be newCells, so we restore it
    fCaloCells->SetName(fCreatedCellBranchName.c_str());
  }
  else
  {
    AliAODCaloCells * currentCells = dynamic_cast<AliAODCaloCells *>(fCaloCellsFromInputEvent);
    AliAODCaloCells * newCells = dynamic_cast<AliAODCaloCells *>(fCaloCells);
    fCaloCells = dynamic_cast<AliAODCaloCells *>(new (newCells) AliAODCaloCells(*currentCells));
    // The name was changed to the currentCells name, but we want it to be newCells, so we restore it
    fCaloCells->SetName(fCreatedCellBranchName.c_str());
  }
  AliDebug(3, Form("Number of old cells: %d \tNumber of new cells: %d", fCaloCellsFromInputEvent->GetNumberOfCells(), fCaloCells->GetNumberOfCells() ));

  // Clusters
  TClonesArray * newClusters = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(fCreatedClusterBranchName.c_str()));
  std::string currentClustersName = AliEmcalCorrectionComponent::DetermineUseDefaultName(AliEmcalCorrectionComponent::kCluster, fIsEsd);
  TClonesArray * currentClusters = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(currentClustersName.c_str()));
  AliDebug(3, Form("before copy:\t currentClusters->GetEntries(): %d \t newClusters->GetEntries(): %d", currentClusters->GetEntries(), newClusters->GetEntries()));
  CopyClusters(currentClusters, newClusters);
  AliDebug(3, Form("after  copy:\t currentClusters->GetEntries(): %d \t newClusters->GetEntries(): %d", currentClusters->GetEntries(), newClusters->GetEntries()));

  // Tracks
  TClonesArray * newTracks = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(fCreatedTrackBranchName.c_str()));
  std::string currentTracksName = AliEmcalCorrectionComponent::DetermineUseDefaultName(AliEmcalCorrectionComponent::kTrack, fIsEsd);
  TClonesArray * currentTracks = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(currentTracksName.c_str()));
  AliDebug(3, Form("before copy:\t currentTracks->GetEntries(): %d \t newTracks->GetEntries(): %d", currentTracks->GetEntries(), newTracks->GetEntries()));
  for (Int_t i = 0; i < currentTracks->GetEntriesFast(); i++)
  {
    if (fIsEsd)
    {
      AliESDtrack *currentTrack = dynamic_cast<AliESDtrack *>(currentTracks->At(i));
      // Calls copy constructor to create a new track at position i in newTracks
      AliESDtrack *newTrack = dynamic_cast<AliESDtrack *>(new ((*newTracks)[i]) AliESDtrack(*currentTrack));

      // Assign new process ID so that it is from the current process as opposed to the old one copied from the current track
      TProcessID::AssignID(newTrack);
    }
    else
    {
      AliAODTrack *currentTrack = dynamic_cast<AliAODTrack *>(currentTracks->At(i));
      // Calls copy constructor to create a new track at position i in newTracks
      AliAODTrack *newTrack = dynamic_cast<AliAODTrack *>(new ((*newTracks)[i]) AliAODTrack(*currentTrack));

      // TEMP
      /*if (i == 5 || i == 7)
      {
        std::cout << "Track properties:" << std::endl;
        currentTrack->Print();
        newTrack->Print();
        std::cout << "ProcessID for current track: " << TProcessID::GetProcessWithUID(currentTrack)->GetName() << "/" << TProcessID::GetProcessWithUID(currentTrack)->GetTitle() << ". Memory Address: " << TProcessID::GetProcessWithUID(currentTrack) << std::endl;
        std::cout << "ProcessID for new track: " << TProcessID::GetProcessWithUID(newTrack)->GetName() << "/" << TProcessID::GetProcessWithUID(newTrack)->GetTitle() << ". Memory Address: " << TProcessID::GetProcessWithUID(newTrack) << std::endl;
      }*/

      // Assign new process ID so that it is from the current process as opposed to the old one copied from the current track
      TProcessID::AssignID(newTrack);

      /*if (i == 5 || i == 7)
      {
        std::cout << "ProcessID for new track after assign ID: " << TProcessID::GetProcessWithUID(newTrack)->GetName() << "/" << TProcessID::GetProcessWithUID(newTrack)->GetTitle() << ". Memory Address: " << TProcessID::GetProcessWithUID(newTrack) << std::endl;
        std::cout << "ProcessID for newTracks: " << TProcessID::GetProcessWithUID(newTracks)->GetName() << "/" << TProcessID::GetProcessWithUID(newTracks)->GetTitle() << ". Memory Address: " << TProcessID::GetProcessWithUID(newTracks) << std::endl;
      }*/
    }
  }
  AliDebug(3, Form("after  copy:\t currentTracks->GetEntries(): %d \t newTracks->GetEntries(): %d", currentTracks->GetEntries(), newTracks->GetEntries()));
}

/**
 *
 *
 */
void AliEmcalCorrectionTask::CopyClusters(TClonesArray *orig, TClonesArray *dest)
{
  const Int_t Ncls = orig->GetEntries();

  for(Int_t i=0; i < Ncls; ++i) {
    AliVCluster *oc = static_cast<AliVCluster*>(orig->At(i));

    if (!oc)
      continue;

    // We likely want to include PHOS cells as well, so the copy we make should avoid this check
    //if (!oc->IsEMCAL())
    //    continue;
    
    AliVCluster *dc = static_cast<AliVCluster*>(dest->New(i));
    //dc->SetType(AliVCluster::kEMCALClusterv1);
    dc->SetType(oc->GetType());
    dc->SetE(oc->E());
    Float_t pos[3] = {0};
    oc->GetPosition(pos);
    dc->SetPosition(pos);
    dc->SetNCells(oc->GetNCells());
    dc->SetCellsAbsId(oc->GetCellsAbsId());
    dc->SetCellsAmplitudeFraction(oc->GetCellsAmplitudeFraction());
    dc->SetID(oc->GetID());
    dc->SetDispersion(oc->GetDispersion());
    dc->SetEmcCpvDistance(-1);
    dc->SetChi2(-1);
    dc->SetTOF(oc->GetTOF());     //time-of-flight
    dc->SetNExMax(oc->GetNExMax()); //number of local maxima
    dc->SetM02(oc->GetM02());
    dc->SetM20(oc->GetM20());
    dc->SetDistanceToBadChannel(oc->GetDistanceToBadChannel()); 
    dc->SetMCEnergyFraction(oc->GetMCEnergyFraction());

    //MC
    UInt_t nlabels = oc->GetNLabels();
    Int_t *labels = oc->GetLabels();

    if (nlabels == 0 || !labels)
      continue;

    AliESDCaloCluster *esdClus = dynamic_cast<AliESDCaloCluster*>(dc);
    if (esdClus) {
      TArrayI parents(nlabels, labels);
      esdClus->AddLabels(parents); 
    }
    else {
      AliAODCaloCluster *aodClus = dynamic_cast<AliAODCaloCluster*>(dc);
      if (aodClus) 
        aodClus->SetLabel(labels, nlabels);
    }
  }
}

/**
 * Retrieve objects from event.
 * @return
 */
Bool_t AliEmcalCorrectionTask::RetrieveEventObjects()
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fNVertCont = 0;

  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  if (vert) {
    vert->GetXYZ(fVertex);
    fNVertCont = vert->GetNContributors();
  }

  fBeamType = GetBeamType();

  if (fBeamType == kAA || fBeamType == kpA ) {
    if (fUseNewCentralityEstimation) {
      AliMultSelection *MultSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
      if (MultSelection) {
        fCent = MultSelection->GetMultiplicityPercentile(fCentEst.Data());
      }
      else {
        AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
      }
    }
    else { // old centrality estimation < 2015
      AliCentrality *aliCent = InputEvent()->GetCentrality();
      if (aliCent) {
        fCent = aliCent->GetCentralityPercentile(fCentEst.Data());
      }
      else {
        AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
      }
    }

    if (fNcentBins==4) {
      if      (fCent >=  0 && fCent <   10) fCentBin = 0;
      else if (fCent >= 10 && fCent <   30) fCentBin = 1;
      else if (fCent >= 30 && fCent <   50) fCentBin = 2;
      else if (fCent >= 50 && fCent <= 100) fCentBin = 3;
      else {
        AliWarning(Form("%s: Negative centrality: %f. Assuming 99", GetName(), fCent));
        fCentBin = fNcentBins-1;
      }
    }
    else if (fNcentBins==5) {  // for PbPb 2015
      if      (fCent >=  0 && fCent <   10) fCentBin = 0;
      else if (fCent >= 10 && fCent <   30) fCentBin = 1;
      else if (fCent >= 30 && fCent <   50) fCentBin = 2;
      else if (fCent >= 50 && fCent <= 90) fCentBin = 3;
      else if (fCent > 90) {
        fCent = 99;
        fCentBin = 4;
      }
      else {
        AliWarning(Form("%s: Negative centrality: %f. Assuming 99", GetName(), fCent));
        fCentBin = fNcentBins-1;
      }
    }
    else {
      Double_t centWidth = (fMaxCent-fMinCent)/(Double_t)fNcentBins;
      if(centWidth>0.) {
        fCentBin = TMath::FloorNint(fCent/centWidth);
      }
      else {
        fCentBin = 0;
      }
      if (fCentBin>=fNcentBins) {
        AliWarning(Form("%s: fCentBin too large: cent = %f fCentBin = %d. Assuming 99", GetName(),fCent,fCentBin));
        fCentBin = fNcentBins-1;
      }
    }
  }
  else {
    fCent = 99;
    fCentBin = 0;
  }

  AliEmcalContainer* cont = 0;

  TIter nextPartColl(&fParticleCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextPartColl()))) cont->NextEvent();

  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextClusColl()))) cont->NextEvent();

  return kTRUE;
}


/**
 * Perform steps needed to initialize the analysis.
 * This function relies on the presence of an input
 * event (ESD or AOD event). Consequently it is called
 * internally by UserExec for the first event.
 *
 * This function connects all containers attached to
 * this task to the corresponding arrays in the
 * input event. Furthermore it initializes the geometry.
 */
void AliEmcalCorrectionTask::ExecOnce()
{
  if (!InputEvent()) {
    AliError(Form("%s: Could not retrieve event! Returning!", GetName()));
    return;
  }

  if (fNeedEmcalGeom) {
    fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
    if (!fGeom) {
      AliFatal(Form("%s: Can not get EMCal geometry instance. If you do not need the EMCal geometry, disable it by setting task->SetNeedEmcalGeometry(kF    ALSE).", GetName()));
      return;
    }
  }

  // Create the new object branch here.
  // The cell pointer already exists. All we need to do is load it later.
  if (fCreateNewObjectBranches)
    CreateNewObjectBranches();

  // Load all requested track branches - each container knows name already
  for (Int_t i =0; i<fParticleCollArray.GetEntriesFast(); i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    cont->SetArray(InputEvent());
  }

  // Load all requested cluster branches - each container knows name already
  for (Int_t i =0; i<fClusterCollArray.GetEntriesFast(); i++) {
    AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
    cont->SetArray(InputEvent());
  }

  if (!fCaloCellsName.IsNull() && !fCaloCellsFromInputEvent) {
    fCaloCellsFromInputEvent =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCaloCellsName));
    if (!fCaloCellsFromInputEvent) {
      AliError(Form("%s: Could not retrieve cells %s!", GetName(), fCaloCellsName.Data())); 
      return;
    }
  }

  fEventInitialized = kTRUE;

  ExecOnceComponents();
}

/**
 * Create new container for MC particles and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new container for MC particles
 */
AliMCParticleContainer* AliEmcalCorrectionTask::AddMCParticleContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;

  AliMCParticleContainer* cont = new AliMCParticleContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

/**
 * Create new track container and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new track container
 */
AliTrackContainer* AliEmcalCorrectionTask::AddTrackContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;

  AliTrackContainer* cont = new AliTrackContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

/**
 * Create new particle container and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new particle container
 */
AliParticleContainer* AliEmcalCorrectionTask::AddParticleContainer(const char *n) 
{
  if (TString(n).IsNull()) return 0;

  AliParticleContainer* cont = new AliParticleContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

/**
 * Create new cluster container and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new cluster container
 */
AliClusterContainer* AliEmcalCorrectionTask::AddClusterContainer(const char *n) 
{
  if (TString(n).IsNull()) return 0;

  AliClusterContainer* cont = new AliClusterContainer(n);

  fClusterCollArray.Add(cont);

  return cont;
}

/**
 * Get \f$ i^{th} \f$ particle container attached to this task
 * @param[in] i Index of the particle container
 * @return Particle container found for the given index (NULL if no particle container exists for that index)
 */
AliParticleContainer* AliEmcalCorrectionTask::GetParticleContainer(Int_t i) const 
{
  if (i<0 || i>fParticleCollArray.GetEntriesFast()) return 0;
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
  return cont;
}

/**
 * Get \f$ i^{th} \f$ cluster container attached to this task
 * @param[in] i Index of the cluster container
 * @return Cluster container found for the given index (NULL if no cluster container exists for that index)
 */
AliClusterContainer* AliEmcalCorrectionTask::GetClusterContainer(Int_t i) const 
{
  if (i<0 || i>fClusterCollArray.GetEntriesFast()) return 0;
  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
  return cont;
}

/**
 * Find particle container attached to this task according to its name
 * @param[in] name Name of the particle container
 * @return Particle container found under the given name
 */
AliParticleContainer* AliEmcalCorrectionTask::GetParticleContainer(const char *name) const 
{
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.FindObject(name));
  return cont;
}

/**
 * Find cluster container attached to this task according to its name
 * @param[in] name Name of the cluster container
 * @return Cluster container found under the given name
 */
AliClusterContainer* AliEmcalCorrectionTask::GetClusterContainer(const char *name) const 
{
  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.FindObject(name));
  return cont;
}



/**
 *
 *
 */
void AliEmcalCorrectionTask::UserExec(Option_t *option)
{
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));

  // Initialize the event if not intialized
  if (!fEventInitialized)
    ExecOnce();

  // Only continue if we are initialized successfully
  if (!fEventInitialized)
    return;

  // Copy cells and cluster to new branches if desired
  if (fCreateNewObjectBranches)
  {
    CopyBranchesToNewObjects();
  }
  else
  {
    fCaloCells = fCaloCellsFromInputEvent;
  }

  // Get the objects for each event
  if (!RetrieveEventObjects())
    return;

  // TODO: Consider adding IsEventSelected()??

  // Call run for each correction
  if (!Run())
    return;
}

/**
 *
 *
 */
Bool_t AliEmcalCorrectionTask::Run()
{
  // Run the initialization for all derived classes.
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  for (auto component : fCorrectionComponents)
  {
    component->SetEvent(InputEvent());
    component->SetMCEvent(MCEvent());
    component->SetCaloCells(fCaloCells);
    component->SetEMCALGeometry(fGeom);
    component->SetCentralityBin(fCentBin);
    component->SetCentrality(fCent);
    component->SetNcentralityBins(fNcentBins);
    AliClusterContainer * tempClusterContainer = GetClusterContainer(component->GetName());
    if (tempClusterContainer)
    {
      component->SetClusterContainer(tempClusterContainer);
    }
    AliParticleContainer * tempParticleContainer = GetParticleContainer(component->GetName());
    if (tempParticleContainer)
    {
      // TEMP
      //Printf("Particle container name: %s, branch: %s", tempParticleContainer->GetName(), tempParticleContainer->GetArrayName().Data());
      component->SetParticleContainer(tempParticleContainer);
    }

    component->Run();
  }
  
  PostData(1, fOutput);
  
  // Need something more sophisticated here
  return kTRUE;
}

/**
 *
 *
 */
Bool_t AliEmcalCorrectionTask::UserNotify()
{
  // Run the initialization for all derived classes.
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  for (auto component : fCorrectionComponents)
  {
    component->UserNotify();
  }

  // Need something more sophisticated here
  return kTRUE;
}

/**
 *
 *
 */
AliEmcalCorrectionTask::~AliEmcalCorrectionTask()
{
  // Destructor
}

/**
 * Get beam type : pp-AA-pA
 * ESDs have it directly, AODs get it from hardcoded run number ranges
 * @return Beam type of the run.
 */
AliEmcalCorrectionTask::BeamType AliEmcalCorrectionTask::GetBeamType()
{
  if (fForceBeamType != kNA)
    return fForceBeamType;

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (esd) {
    const AliESDRun *run = esd->GetESDRun();
    TString beamType = run->GetBeamType();
    if (beamType == "p-p")
      return kpp;
    else if (beamType == "A-A")
      return kAA;
    else if (beamType == "p-A")
      return kpA;
    else
      return kNA;
  } else {
    Int_t runNumber = InputEvent()->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
        (runNumber >= 166529 && runNumber <= 170593)) {  // LHC11h
      return kAA;
    } else if ((runNumber>=188365 && runNumber <= 188366) ||   // LHC12g
        (runNumber >= 195344 && runNumber <= 196608)) { // LHC13b-f
      return kpA;
    } else {
      return kpp;
    }
  }  
}
