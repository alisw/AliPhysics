// AliEmcalCorrectionTask
//
// Runs the cell and cluster level corrections for the EMCal
//

#include "AliEmcalCorrectionTask.h"
#include "AliEmcalCorrectionComponent.h"

#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <iostream>

#include <TChain.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TFile.h>

#include "AliVEventHandler.h"
#include "AliEMCALGeometry.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliAnalysisManager.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionTask);
/// \endcond

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellContainer);
/// \endcond

/**
 * Default constructor.
 *
 */
AliEmcalCorrectionTask::AliEmcalCorrectionTask() :
  AliAnalysisTaskSE("AliEmcalCorrectionTask"),
  fSuffix(""),
  fUserConfiguration(),
  fUserConfigurationFilename(""),
  fUserConfigurationString(""),
  fDefaultConfiguration(),
  fDefaultConfigurationFilename(""),
  fDefaultConfigurationString(""),
  fCorrectionComponents(),
  fIsEsd(false),
  fForceBeamType(kNA),
  fRunPeriod("noSetRunPeriod"),
  fConfigurationInitialized(false),
  fOrderedComponentsToExecute(),
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
  fNeedEmcalGeom(kTRUE),
  fGeom(0),
  fParticleCollArray(),
  fClusterCollArray(),
  fCellCollArray(),
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
 * Standard constructor.
 *
 * If "_" is included in the name, then all characters after the underscore will
 * be taken as the task suffix (sometimes described as a "specialization") and used
 * to select settings in the YAML configuration file.
 *
 * @param[in] name Name of the correction task.
 */
AliEmcalCorrectionTask::AliEmcalCorrectionTask(const char * name) :
  AliAnalysisTaskSE(name),
  fSuffix(""),
  fUserConfiguration(),
  fUserConfigurationFilename(""),
  fUserConfigurationString(""),
  fDefaultConfiguration(),
  fDefaultConfigurationFilename(""),
  fDefaultConfigurationString(""),
  fCorrectionComponents(),
  fIsEsd(false),
  fForceBeamType(kNA),
  fRunPeriod("noSetRunPeriod"),
  fConfigurationInitialized(false),
  fOrderedComponentsToExecute(),
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
  fNeedEmcalGeom(kTRUE),
  fGeom(0),
  fParticleCollArray(),
  fClusterCollArray(),
  fCellCollArray(),
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

/**
 * Destructor
 *
 */
AliEmcalCorrectionTask::~AliEmcalCorrectionTask()
{
  // Destructor
}

/**
 * Initializes the Correction Task by initializing the YAML configuration and selected correction components,
 * including setting up the input objects (cells, clusters, and tracks).
 *
 * This function is the main function for initialization and should be called from a run macro!
 * Once called, most of the configuration of the correction task and the correction components is locked in,
 * so be certain to change any additional configuration before that!
 */
void AliEmcalCorrectionTask::Initialize()
{
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

  // Determine the suffix of the correction task
  std::string tempName = GetName();
  std::size_t foundSuffix = tempName.find("_");
  if (foundSuffix != std::string::npos) {
    // +1 to skip "_"
    fSuffix = tempName.substr(foundSuffix + 1).c_str();
  }

  // Initialize YAML configuration
  InitializeConfiguration();
  // Check that the configuration is initialized
  if (fConfigurationInitialized != true)
  {
    AliFatal("YAML configuration must be initialized before running (ie. the AddTask, run macro or wagon)!");
  }

  // Determine component execution order
  DetermineComponentsToExecute(fOrderedComponentsToExecute);

  // Check for user defined settings that are not in the default file
  CheckForUnmatchedUserSettings();

  // Setup input objects
  // Setup Cells
  // Cannot do this entirely yet because we need input objects
  CreateInputObjects(kCaloCells);
  PrintRequestedContainersInformation(kCaloCells);
  // Create cluster input objects
  CreateInputObjects(kCluster);
  PrintRequestedContainersInformation(kCluster);
  // Create track input objects
  CreateInputObjects(kTrack);
  PrintRequestedContainersInformation(kTrack);

  // Initialize components
  InitializeComponents();
}

/**
 * Initializes and sets up the user and default configuration files.
 * This includes opening the files, checking to ensure that the run period is valid, and
 * storing the contents of the user and default YAML files into strings so that they can
 * be streamed to the grid. (yaml-cpp objects do not work properly with ROOT streamers).
 *
 * NOTE: fConfigurationInitialized is set to true if the function is successful.
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

  if (DoesFileExist(fDefaultConfigurationFilename) == true)
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

  if (DoesFileExist(fUserConfigurationFilename) == true)
  {
    AliInfo(TString::Format("Using user EMCal corrections configuration located at %s", fUserConfigurationFilename.c_str()));

    fUserConfiguration = YAML::LoadFile(fUserConfigurationFilename);
  }
  else
  {
    AliInfo(TString::Format("User file at \"%s\" does not exist! All settings will be from the default file!", fUserConfigurationFilename.c_str()));
  }

  // Ensure that there is a run period
  if (fRunPeriod == "noSetRunPeriod")
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
  if (userRunPeriod != "knouserfile" && userRunPeriod != fRunPeriod)
  {
    AliFatal(TString::Format("User run period \"%s\" does not match the run period of \"%s\" passed to the correction task!", userRunPeriod.Data(), fRunPeriod.Data()));
  }

  // "" means the user wants their settings to apply to all periods
  // Ensure that the user is aware!
  if (userRunPeriod == "")
  {
    AliWarning("User run period is an empty string. Settings apply to all run periods!");
  }

  // Save configuration into strings so that they can be streamed
  // Need the stringstream because YAML implements streamers
  std::stringstream tempConfiguration;
  tempConfiguration << fUserConfiguration;
  fUserConfigurationString = tempConfiguration.str();
  tempConfiguration.str("");
  tempConfiguration << fDefaultConfiguration;
  fDefaultConfigurationString = tempConfiguration.str();

  //AliInfo(TString::Format("User configuration: %s", fUserConfigurationString.c_str()));
  //AliInfo(TString::Format("Default configuration: %s", fDefaultConfigurationString.c_str()));

  // Note that it is initialized properly so that the analysis can proceed
  fConfigurationInitialized = true;
}

/**
 * Determines which components to execute based on which are selected via specialization (ie suffix),
 * as well as which are enabled and in what order they should be executed. The returned components in
 * order of which they should be executed.
 *
 * It is recommended to store the result, as it requires a large number of string comparisons and lookups
 * in the YAML configuration, whose associated methods also utilize a number of string comparisons.
 *
 * @param[out] correctionComponents Names of the selected correction components in the order in which they should be executed.
 */
void AliEmcalCorrectionTask::DetermineComponentsToExecute(std::vector <std::string> & correctionComponents)
{
  std::vector <std::string> executionOrder;
  // executionOrder determines the order of tasks to execute, but it doesn't name the particular tasks
  AliEmcalCorrectionComponent::GetProperty("executionOrder", executionOrder, fUserConfiguration, fDefaultConfiguration, true);

  // Possible components to create from both the user and default configurations
  // Use set so that the possible components are not repeated
  std::set <std::string> possibleComponents;
  for (auto node : fUserConfiguration) {
    possibleComponents.insert(node.first.as<std::string>());
  }
  for (auto node : fDefaultConfiguration) {
    possibleComponents.insert(node.first.as<std::string>());
  }

  // Determine the correction names associated with the correction task
  std::string expectedComponentName = "";
  bool foundSuffixComponent = false;
  bool foundComponent = false;
  bool componentEnabled = true;

  // Execution order determines the order that corrections should be added to our execution list
  for (auto & execName : executionOrder)
  {
    // Construct the expected component name with the suffix
    expectedComponentName = TString::Format("%s_%s", execName.c_str(), fSuffix.c_str()).Data();
    foundComponent = false;
    componentEnabled = false;

    foundComponent = CheckPossibleNamesForComponentName(expectedComponentName, possibleComponents);
    if (foundComponent)
    {
      // Check if the component is enabled
      AliEmcalCorrectionComponent::GetProperty("enabled", componentEnabled, fUserConfiguration, fDefaultConfiguration, true, expectedComponentName);
      // If enabled, then store the name so that it can be executed
      if (componentEnabled == true) {
        foundSuffixComponent = true;
        correctionComponents.push_back(expectedComponentName);
      }
      else {
        AliInfo(TString::Format("Component %s is disabled and will not be run!", expectedComponentName.c_str()));
      }

      continue;
    }
    else
    {
      // Look for the normal component
      expectedComponentName = execName;
      foundComponent = CheckPossibleNamesForComponentName(expectedComponentName, possibleComponents);
      // Check if it is enabled
      AliEmcalCorrectionComponent::GetProperty("enabled", componentEnabled, fUserConfiguration, fDefaultConfiguration, true, expectedComponentName);

      if (componentEnabled == true) {
        if (foundSuffixComponent == true) {
          AliFatal(TString::Format("Found earlier component %s with suffix \"%s\", but could not found component %s with that same suffix!", correctionComponents.back().c_str(), fSuffix.c_str(), expectedComponentName.c_str()));
        }
        else {
          // Take the normal component and store it to be executed
          correctionComponents.push_back(expectedComponentName);
        }
      }
      else {
        AliInfo(TString::Format("Component %s is disabled and will not be run!", expectedComponentName.c_str()));
      }
    }
  }

  // Need to append "AliEmcalCorrection" to allow the tasks to be found!
  AliDebug(2, "Found EMCal Correction Components: ");
  for (auto & component : correctionComponents)
  {
    component = "AliEmcalCorrection" + component;
    AliDebug(2, TString::Format("%s", component.c_str()) );
  }
}

/**
 * Check each property defined in the user configuration file for a match to the properties in the default
 * configuration file. The default configuration file should have all possible properties defined, so it can
 * be treated as a reference. Thus, if there is a property that is defined in the user configuration that
 * cannot be found in the default configuration, then it must be a typo and we should alert the user.
 *
 * Utilizes the stored ordered correction component names.
 */
void AliEmcalCorrectionTask::CheckForUnmatchedUserSettings()
{
  // Names of properties for a particular component in the user and default configurations
  std::set <std::string> userPropertyNames;
  std::set <std::string> defaultPropertyNames;
  // Notes whether a match was found between user and default properties
  bool foundMatch = false;
  std::string tempComponentName = "";

  // Loop over all components
  for (const auto componentName : fOrderedComponentsToExecute)
  {
    // Reset for each loop
    userPropertyNames.clear();
    defaultPropertyNames.clear();
    // We need to remove "AliEmcalCorrection" so that the correction will be found in the configuration
    // "AliEmcalCorrection" is 18 characters
    tempComponentName = componentName.substr(componentName.find("AliEmcalCorrection")+18);

    AliDebugStream(2) << "Checking component " << componentName << " for unmatched user settings" << std::endl;

    // Get the user property names
    GetPropertyNamesFromNode(tempComponentName, fUserConfiguration, userPropertyNames, false);

    // Get the same from default
    // Not required here because the default configuration may not have the specialized component
    GetPropertyNamesFromNode(tempComponentName, fDefaultConfiguration, defaultPropertyNames, false);

    // We need to check the base correction as well to fill out the options
    if (tempComponentName.find("_") != std::string::npos) {
      // Get the base user component
      GetPropertyNamesFromNode(tempComponentName.substr(0, tempComponentName.find("_")), fUserConfiguration, userPropertyNames, false);

      // Required here because the default must have the base component!
      GetPropertyNamesFromNode(tempComponentName.substr(0, tempComponentName.find("_")), fDefaultConfiguration, defaultPropertyNames, true);
    }

    // Check each property defined in the user file for a match to the properties in the default file
    for (auto userPropertyName : userPropertyNames)
    {
      AliDebugStream(2) << "Checking property " << userPropertyName << std::endl;
      foundMatch = false;
      for (auto defaultPropertyName : defaultPropertyNames)
      {
        if (userPropertyName == defaultPropertyName) {
          AliDebugStream(2) << "Found match of " << userPropertyName << " with " << defaultPropertyName << std::endl;
          foundMatch = true;
        }
      }
      if (foundMatch == false) {
        AliFatal(TString::Format("Property \"%s:%s\" defined in the user configuration file cannot be found in the default configuration file! Check the spelling in your user file!", tempComponentName.c_str(), userPropertyName.c_str()));
      }
    }
  }
}

/**
 * Creates, configures, and initializes components based on the configuration described in the YAML files.
 * Configuration includes making the user and default configuration available to each component, as
 * well as setting up the input clusters and tracks (cells have to be handled during ExecOnce()).
 * Each component's individual initialization is also called.
 */
void AliEmcalCorrectionTask::InitializeComponents()
{
  // Iterate over the ordered components list and create the components
  AliEmcalCorrectionComponent * component = 0;
  for (auto componentName : fOrderedComponentsToExecute)
  {
    std::string noPrefixComponentName = componentName.substr(0, componentName.find("_" + fSuffix));
    component = AliEmcalCorrectionComponentFactory::createInstance(noPrefixComponentName);
    if (!component)
    {
      AliFatal(TString::Format("Failed to create requested component %s!", componentName.c_str()));
    }

    // For setting names of tasks to differentiate between tasks of the same class
    component->SetName(componentName.c_str());
    component->SetTitle(componentName.c_str());

    // Initialize the YAML configurations in each component
    component->SetUserConfiguration(fUserConfiguration);
    component->SetDefaultConfiguration(fDefaultConfiguration);

    // configure needed fields for components to properly initialize
    component->SetIsESD(fIsEsd);

    // Add the require containers to the component
    // Cells must be set during UserExec() because we need to add them as a pointer
    AddContainersToComponent(component, kCluster);
    AddContainersToComponent(component, kTrack);

    // Initialize each component
    component->Initialize();

    if (component)
    {
      AliInfo(TString::Format("Successfully added correction task: %s", componentName.c_str()));
      fCorrectionComponents.push_back(component);
    }
  }
}

/**
 * Steers the creation of all input objects of a selected type. In the case of clusters and tracks, 
 * AliEmcalContainer can be used, and a container is created for each requested object of the specified
 * input object type. In the case of cells, an EMCal Correction Cell Container is created to handle the
 * relevant information.
 *
 * Normally, when properties are retrieved from the YAML configuration, the entire configuration is considered.
 * However, here we only consider a subset in both the user and default configurations. In particular, the
 * "inputObjects" section is selected and then all properties are drawn from this subset (the shared parameters
 * are also retained so they can be used). Thus, it takes a bit more care to use this task, but it should be
 * entirely transparent to the users - no functionality or expected flexibility is lost.
 *
 * @param inputObjectType Type of the input object(s) to create
 */
void AliEmcalCorrectionTask::CreateInputObjects(InputObject_t inputObjectType)
{
  // Get container node
  std::string inputObjectName = GetInputFieldNameFromInputObjectType(inputObjectType);

  // Get the user and default input nodes for the object type
  YAML::Node userInputObjectNode;
  YAML::Node defaultInputObjectNode;
  GetNodeForInputObjects(userInputObjectNode, fUserConfiguration, inputObjectName, false);
  GetNodeForInputObjects(defaultInputObjectNode, fDefaultConfiguration, inputObjectName, true);

  AliDebugStream(3) << "userInputObjectNode: " << userInputObjectNode << std::endl;
  AliDebugStream(3) << "defaultInputObjectNode: " << defaultInputObjectNode << std::endl;

  // Determine which containers we need based on which are requested by the enabled correction tasks
  std::set <std::string> requestedContainers;
  std::vector <std::string> componentRequest;
  for ( const auto & componentName : fOrderedComponentsToExecute )
  {
    componentRequest.clear();
    // Not required because not all components will have all kinds of containers
    // "AliEmcalCorrection" is 18 characters
    AliEmcalCorrectionComponent::GetProperty(inputObjectName + "Names", componentRequest, fUserConfiguration, fDefaultConfiguration, false, componentName.substr(componentName.find("AliEmcalCorrection")+18));
    for ( auto & req : componentRequest )
    {
      AliDebugStream(3) << "Component " << componentName << " requested container name " << req << std::endl;
      requestedContainers.insert(req);
    }
  }

  AliDebugStream(2) << inputObjectName << " Containers requested by components: " << std::endl;
  for (auto & str : requestedContainers) {
    AliDebugStream(2) << "\t" << str << std::endl;;
  }

  // Create all requested containers
  AliDebug(2, TString::Format("Setting up requested containers!"));
  SetupContainersFromInputNodes(inputObjectType, userInputObjectNode, defaultInputObjectNode, requestedContainers);
}

/**
 * Adds the previously created input objects that are managed by the Correction Task into a correction
 * component based on which input objects are requested in the YAML configuration by the component.
 *
 * @param[in] component The correction component to which the input objects will be added
 * @param[in] inputObjectType The type of input object to add to the component
 */
void AliEmcalCorrectionTask::AddContainersToComponent(AliEmcalCorrectionComponent * component, InputObject_t inputObjectType)
{
  std::string inputObjectName = GetInputFieldNameFromInputObjectType(inputObjectType);
  // Need to be of the form "clusterContainersNames"
  inputObjectName = inputObjectName + "Names";

  std::vector <std::string> inputObjects;
  // Not required, because not all components need Clusters or Tracks
  AliEmcalCorrectionComponent::GetProperty(inputObjectName.c_str(), inputObjects, fUserConfiguration, fDefaultConfiguration, false, component->GetName());

  //AliDebugStream(4) << "inputObjects.size(): " << inputObjects.size() << std::endl;

  // If it is not found, then there will be nothing to iterate over, so we don't need to explicitly check the return value
  for (auto const & str : inputObjects)
  {
    // TODO: Generalize to arrays for clusters and tracks...
    // NOTE: The AliEmcalContainer derived objects operate differently than the cells. The containers should be added during initialization while the cells should be added during ExecOnce()!
    if (inputObjectType == kCluster)
    {
      AliEmcalContainer * cont = GetClusterContainer(str.c_str());
      AliDebugStream(2) << "Adding cluster container " << str << " of array " << cont->GetArrayName() << " to component " << component->GetName() << std::endl;
      component->SetClusterContainer(GetClusterContainer(str.c_str()));
    }
    else if (inputObjectType == kTrack)
    {
      AliEmcalContainer * cont = GetParticleContainer(str.c_str());
      AliDebugStream(2) << "Adding particle container " << str << " of array " << cont->GetArrayName() << " to component " << component->GetName() << std::endl;
      component->SetParticleContainer(GetParticleContainer(str.c_str()));
    }
    else if (inputObjectType == kCaloCells)
    {
      // NOTE: This operates different than the others. This should be executed during run time rather than during initialization!
      if (inputObjects.size() > 1) {
        AliFatal(TString::Format("Component %s requested more than one cell branch, but this is not supported! Check the configuration!", component->GetName()));
      }

      // If we've made it here, this must be at least one entry
      AliDebugStream(2) << "Adding calo cells " << GetCellContainer(str)->GetName() << " of branch name " << GetCellContainer(str)->GetBranchName() << "to component " << component->GetName() << std::endl;
      component->SetCaloCells(GetCellContainer(str)->GetCells());
      AliDebugStream(3) << "component GetNumberOfCells: " << component->GetCaloCells()->GetNumberOfCells() << std::endl;
    }
  }
}

/**
 * Creates the input objects containers requested by the components.
 *
 * @param[in] inputObjectType Type of the input objects to create
 * @param[in] userInputObjectNode YAML Node corresponding to the user input objects configuration
 * @param[in] defaultInputObjectNode YAML Node corresponding to the default input objects configuration
 * @param[in] requestedContainers Containers to be created
 */
void AliEmcalCorrectionTask::SetupContainersFromInputNodes(InputObject_t inputObjectType, YAML::Node & userInputObjectNode, YAML::Node & defaultInputObjectNode, std::set <std::string> & requestedContainers)
{
  // Our node contains all of the objects that we will want to create.
  for(auto & containerName : requestedContainers)
  {
    // The section is the container name
    //std::string containerName = it->first.as<std::string>();
    // Skip if the particle or cluster container already exists
    if (GetParticleContainer(containerName.c_str()) || GetClusterContainer(containerName.c_str())) {
      continue;
    }

    AliDebug(2, TString::Format("Processing container %s of inputType %d", containerName.c_str(), inputObjectType));
    if (inputObjectType == kCluster || inputObjectType == kTrack) {
      SetupContainer(inputObjectType, containerName, userInputObjectNode, defaultInputObjectNode);
    }
    else if (inputObjectType == kCaloCells) {
      SetupCellsInfo(containerName, userInputObjectNode, defaultInputObjectNode);
    }
  }
}

/**
 * Setup cell container with information from the YAML configuration nodes corresponding to the selected
 * input object.
 *
 * The created cell containers is stored by in the correction task.
 *
 * @param[in] containerName Name of the container to create (as defined in the YAML configuration)
 * @param[in] userNode YAML Node corresponding to the user input object's configuration
 * @param[in] defaultNode YAML Node corresponding to the default input object's configuration
 */
void AliEmcalCorrectionTask::SetupCellsInfo(std::string containerName, YAML::Node & userNode, YAML::Node & defaultNode)
{
  // Define cell info
  AliEmcalCorrectionCellContainer * cellObj = new AliEmcalCorrectionCellContainer();

  AliDebugStream(2) << "User: " << std::endl << userNode << std::endl << "default: " << std::endl << defaultNode << std::endl;

  // Set properties
  // Cells (object) name
  cellObj->SetName(containerName);
  // Branch name
  std::string tempString = "";
  AliEmcalCorrectionComponent::GetProperty("branchName", tempString, userNode, defaultNode, true, containerName);
  if (tempString == "usedefault") {
    tempString = DetermineUseDefaultName(kCaloCells, fIsEsd);
  }
  cellObj->SetBranchName(tempString);

  // IsEmbedding
  bool tempBool = false;
  AliEmcalCorrectionComponent::GetProperty("embedded", tempString, userNode, defaultNode, false, containerName);
  cellObj->SetIsEmbedding(tempBool);

  // Add to the array to keep track of it
  fCellCollArray.push_back(cellObj);
}

/**
 * Configures AliEmcalContainer derived tasks. Sets both general properties (such as min energy, etc),
 * as well as container type specific properties (such as Non-linearity energy cut for a cluster container).
 * Available options are enumerated in the general documentation related to the Correction Framework, available
 * [here](\ref READMEemcCorrections).
 *
 * Note for experts: This implementation explicitly calls particular functions, which is much more
 * straightforward to understand and work with, but is less flexible. This means that any new options need to
 * be implemented here by hand. An idea for a better implementation is mentioned in comments in the function
 * for those who are interested. Once it is implemented, AliEmcalContainer derived classes could be entirely
 * configured through YAML.
 *
 * @param[in] inputObjectType Type of the input object to configure
 * @param[in] containerName Name of the container to create (as defined in the YAML configuration)
 * @param[in] userNode YAML Node corresponding to the user input object's configuration
 * @param[in] defaultNode YAML Node corresponding to the default input object's configuration
 */
void AliEmcalCorrectionTask::SetupContainer(InputObject_t inputObjectType, std::string containerName, YAML::Node & userNode, YAML::Node & defaultNode)
{
  // Create container
  AliDebugStream(2) << "Adding container" << std::endl;
  AliEmcalContainer * cont = AddContainer(inputObjectType, containerName, userNode, defaultNode);
  AliDebugStream(2) << "Added container" << std::endl;

  // Set the container properties
  //
  // TODO: Consider if this can be converted to a map to function pointers. There are a number of details
  //       which can make it a bit complicated. Those details include inheritance, pointing to member
  //       functions, etc. It should all be possible, but may not be worth all of the extra work and code.
  //       Example ccode:
  //          SetValueInContainer("minPt", &cont::SetMinPt, tempDouble, userNode, defaultNode);
  //          SetValueInContainer("minE", &cont::SetMinE, tempDouble, userNode, defaultNode);

  // Temporary variables to store requested properties
  std::string tempString = "";
  Double_t tempDouble = 0;
  bool tempBool = false;

  // AliEmcalContainer properties
  // Min Pt
  bool result = AliEmcalCorrectionComponent::GetProperty("minPt", tempDouble, userNode, defaultNode, false, containerName);
  if (result) {
    AliDebugStream(2) << cont->GetName() << ": Setting minPt of " << tempDouble << std::endl;
    cont->SetMinPt(tempDouble);
  }
  // Min E
  result = AliEmcalCorrectionComponent::GetProperty("minE", tempDouble, userNode, defaultNode, false, containerName);
  if (result) {
    AliDebugStream(2) << cont->GetName() << ": Setting minE of " << tempDouble << std::endl;
    cont->SetMinE(tempDouble);
  }
  // Eta min, max
  result = AliEmcalCorrectionComponent::GetProperty("minEta", tempDouble, userNode, defaultNode, false, containerName);
  if (result) {
    // Only continue checking if the min is there, since we must set both together
    Double_t tempDouble2 = 0;
    result = AliEmcalCorrectionComponent::GetProperty("maxEta", tempDouble, userNode, defaultNode, false, containerName);
    if (result) {
      AliDebugStream(2) << cont->GetName() << ": Setting eta limits of " << tempDouble << " to " << tempDouble2 << std::endl;
      cont->SetEtaLimits(tempDouble, tempDouble2);
    }
  }
  // Phi min, max
  result = AliEmcalCorrectionComponent::GetProperty("minPhi", tempDouble, userNode, defaultNode, false, containerName);
  if (result) {
    // Only continue checking if the min is there, since we must set both together
    Double_t tempDouble2 = 0;
    result = AliEmcalCorrectionComponent::GetProperty("maxPhi", tempDouble, userNode, defaultNode, false, containerName);
    if (result) {
      AliDebugStream(2) << cont->GetName() << ": Setting phi limits of " << tempDouble << " to " << tempDouble2 << std::endl;
      cont->SetPhiLimits(tempDouble, tempDouble2);
    }
  }
  // Embedded
  // TODO: Enable embedded when that branch is committed!
  /*result = AliEmcalCorrectionComponent::GetProperty("IsEmbedded", tempBool, userNode, defaultNode, false, containerName);
  if (result) {
    AliDebugStream(2) << cont->GetName() << ": Setting embedding to " << (tempBool ? "enabled" : "disabled") << std::endl;
    cont->SetIsEmbedding(tempBool);
  }*/

  // Cluster specific properties
  AliClusterContainer * clusterContainer = dynamic_cast<AliClusterContainer *>(cont);
  if (clusterContainer) {
    // Default energy
    // Probably not needed for the corrections
    /*result = AliEmcalCorrectionComponent::GetProperty("defaultClusterEnergy", tempString, userNode, defaultNode, false, containerName);
    if (result) {
      // Need to get the enumeration
      AliVCluster::VCluUserDefEnergy_t clusterEnergyType = clusterEnergyTypeMap.at(tempString);
      AliDebugStream(2) << clusterContainer->GetName() << ": Setting cluster energy type to " << clusterEnergyType << std::endl;
      clusterContainer->SetDefaultClusterEnergy(clusterEnergyType);
    }*/

    // NonLinCorrEnergyCut
    result = AliEmcalCorrectionComponent::GetProperty("clusNonLinCorrEnergyCut", tempDouble, userNode, defaultNode, false, containerName);
    if (result) {
      AliDebugStream(2) << clusterContainer->GetName() << ": Setting clusNonLinCorrEnergyCut of " << tempDouble << std::endl;
      clusterContainer->SetClusNonLinCorrEnergyCut(tempDouble);
    }

    // HadCorrEnergyCut
    result = AliEmcalCorrectionComponent::GetProperty("clusHadCorrEnergyCut", tempDouble, userNode, defaultNode, false, containerName);
    if (result) {
      AliDebugStream(2) << clusterContainer->GetName() << ": Setting clusHadCorrEnergyCut of " << tempDouble << std::endl;
      clusterContainer->SetClusHadCorrEnergyCut(tempDouble);
    }

    // SetIncludePHOS
    result = AliEmcalCorrectionComponent::GetProperty("includePHOS", tempBool, userNode, defaultNode, false, containerName);
    if (result) {
      AliDebugStream(2) << clusterContainer->GetName() << ": Setting Include PHOS to " << (tempBool ? "enabled" : "disabled") << std::endl;
      clusterContainer->SetIncludePHOS(tempBool);
    }
  }

  // Track specific
  AliTrackContainer * trackContainer = dynamic_cast<AliTrackContainer *>(cont);
  if (trackContainer) {
    // Track selection
    // AOD Filter bits as a sequence
    std::vector <UInt_t> filterBitsVector;
    result = AliEmcalCorrectionComponent::GetProperty("aodFilterBits", filterBitsVector, userNode, defaultNode, false, containerName);
    if (result){
      UInt_t filterBits = 0;
      for (int filterBit : filterBitsVector) {
        filterBits += filterBit;
      }
      AliDebugStream(2) << trackContainer->GetName() << ": Setting filterBits of " << filterBits << std::endl;
      trackContainer->SetAODFilterBits(filterBits);
    }

    // SetTrackFilterType enum
    result = AliEmcalCorrectionComponent::GetProperty("trackFilterType", tempString, userNode, defaultNode, false, containerName);
    if (result) {
      // Need to get the enumeration
      AliEmcalTrackSelection::ETrackFilterType_t trackFilterType = trackFilterTypeMap.at(tempString);
      AliDebugStream(2) << trackContainer->GetName() << ": Setting trackFilterType of " << trackFilterType << std::endl;
      trackContainer->SetTrackFilterType(trackFilterType);
    }

    // Track cuts period
    result = AliEmcalCorrectionComponent::GetProperty("trackCutsPeriod", tempString, userNode, defaultNode, false, containerName);
    if (result) {
      // Need to get the enumeration
      AliDebugStream(2) << trackContainer->GetName() << ": Setting track cuts period to " << tempString << std::endl;
      trackContainer->SetTrackCutsPeriod(tempString.c_str());
    }
  }
}

/**
 * Creates a new AliEmcalContainer derived container based on the requested type and the branch name set in
 * the user and default YAML configuration and requested by a particular correction component. Supports the
 * "usedefault" pattern to simplify setting the proper branch name. Any track input objects are created as 
 * track containers unless the branch is named "mcparticles". If this is problematic for your analysis, the
 * track selection behavior of the track container can be disabled by setting the track selection to
 * "kNoTrackFilter".
 *
 * Note that the created container is adopted and managed by the Correction Task. 
 *
 * @param[in] contType Type of the input object to add
 * @param[in] containerName Name of the container to create (as defined in the YAML configuration)
 * @param[in] userNode YAML Node corresponding to the user input object's configuration
 * @param[in] defaultNode YAML Node corresponding to the default input object's configuration
 *
 * @return The created container
 */
AliEmcalContainer * AliEmcalCorrectionTask::AddContainer(InputObject_t contType, std::string & containerName, YAML::Node & userNode, YAML::Node & defaultNode)
{
  // Determine the type of branch to request
  std::string containerBranch = "";
  if (contType != kCluster && contType != kTrack){
    AliFatal("Must specify type of container when requesting branch.");
  }

  // Retrieve branch name
  // YAML::Node() is just an empty node
  AliDebugStream(2) << "User Node: " << userNode << std::endl;
  AliDebugStream(2) << "Default Node: " << defaultNode << std::endl;
  AliEmcalCorrectionComponent::GetProperty("branchName", containerBranch, userNode, defaultNode, true, containerName);
  // Should be unnecessary, since the user can only do this if done explicitly.
  /*if (containerBranch == "")
  {
    AliFatal(TString::Format("Request %i container, but the container branch is empty!", contType));
  }*/

  // Determine proper name if using "usedefault" pattern
  if (containerBranch == "usedefault") {
    containerBranch = DetermineUseDefaultName(contType, fIsEsd);
  }

  // Create containers and set them to the name of the component
  AliEmcalContainer * cont = 0;
  if (contType == kCluster)
  {
    cont = new AliClusterContainer(containerBranch.c_str());
    AdoptClusterContainer(dynamic_cast<AliClusterContainer *>(cont));
  }
  else if (contType == kTrack)
  {
    if (containerBranch == "mcparticles") {
      cont = new AliMCParticleContainer(containerBranch.c_str());
    }
    else {
      cont = new AliTrackContainer(containerBranch.c_str());
    }
    AdoptParticleContainer(dynamic_cast<AliParticleContainer *>(cont));
  }
  cont->SetName(containerName.c_str());

  return cont;
}

/**
 * Reinitializes the YAML configurations if necessary and sets up for output from the correction components.
 * The reinitialization is necessary if the object is streamed because yaml-cpp objects cannot be streamed.
 * Instead, the YAML configuration is stored in strings and the nodes are recreated here from the string.
 *
 * Note that the number of centrality bins is also set here in the case of a forced beam-type since this is
 * how it was done in AliAnalysisTaskEmcal.
 */
void AliEmcalCorrectionTask::UserCreateOutputObjects()
{
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));

  // Check that the configuration is initialized
  if (fConfigurationInitialized != true)
  {
    AliFatal("YAML configuration must be initialized before running (ie. the AddTask, run macro or wagon)!");
  }

  // Show the configurations info this is available
  AliDebugStream(4) << "User configuration string: " << fUserConfigurationString << std::endl;
  AliDebugStream(4) << "User configuration: " << fUserConfiguration << std::endl;
  AliDebugStream(4) << "Default configuration string: " << fDefaultConfigurationString << std::endl;
  AliDebugStream(4) << "Default configuration: " << fDefaultConfiguration << std::endl;

  // YAML Objects cannot be streamed, so we need to reinitialize them here.
  // They need reinitialize if they are null
  if (fUserConfiguration.IsNull() == true && fUserConfigurationString != "")
  {
    AliInfo("Reinitializing user configuration from string. Expected if running on grid!");
    fUserConfiguration = YAML::Load(fUserConfigurationString);
  }
  if (fDefaultConfiguration.IsNull() == true)
  {
    AliInfo("Reinitializing default configuration from string. Expected if running on grid!");
    fDefaultConfiguration = YAML::Load(fDefaultConfigurationString);
  }

  // Debug to check that the configuration has been (re)initiailzied has been completed correctly
  AliDebugStream(4) << "(Re)initialized user configuration: " << fUserConfigurationString << std::endl;
  AliDebugStream(4) << "(Re)initialized user configuration: " << fUserConfiguration << std::endl;
  AliDebugStream(4) << "(Re)initialized default configuration: " << fDefaultConfigurationString << std::endl;
  AliDebugStream(4) << "(Re)initialized default configuration: " << fDefaultConfiguration << std::endl;

  if (fForceBeamType == kpp)
    fNcentBins = 1;

  // Allow for output files
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

  UserCreateOutputObjectsComponents();

  PostData(1, fOutput);
}

/**
 * Calls UserCreateOutputObjects() for each component and ensures that the output from the correction
 * components is eventually output by the correction task.
 *
 * It also sets the number of centrality bins.
 */
void AliEmcalCorrectionTask::UserCreateOutputObjectsComponents()
{
  // Run the initialization for all derived classes.
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  for (auto component : fCorrectionComponents)
  {
    // Set cent bins (usually used for hist creation)
    // It cannot be set until now because it can be changed after initialization
    // For instance, by SetForceBeamType()
    component->SetNcentralityBins(fNcentBins);

    component->UserCreateOutputObjects();

    if (component->GetOutputList() != 0)
    {
      // Adds a list to the list -- this doesn't work for some unknown reason
      //fOutput->Add(component->GetOutputList());

      // iterate through lists for each component, and fill in output
      TList* t = new TList();
      t->SetName(component->GetName());
      fOutput->AddLast(t);
      t = (TList*)fOutput->Last();
      TIter next(component->GetOutputList());
      while (TObject *obj = next()){
        t->Add(obj);
      }

      AliDebug(1, TString::Format("Added output list from task %s to output.", component->GetName()));
    }
  }
}

/**
 * Steers each event. It enforces that the event is initialized before executing the main analysis of the event.
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

  // Get the objects for each event
  if (!RetrieveEventObjects())
    return;

  // TODO: Consider adding IsEventSelected()??

  // Call run for each correction
  if (!Run())
    return;
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
    AliError("Could not retrieve event! Returning!");
    return;
  }

  if (fNeedEmcalGeom) {
    fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
    if (!fGeom) {
      AliFatal("Can not get EMCal geometry instance. If you do not need the EMCal geometry, disable it by setting task->SetNeedEmcalGeometry(kFALSE).");
      return;
    }
  }

  // Load all requested track branches - each container knows name already
  for (Int_t i =0; i<fParticleCollArray.GetEntriesFast(); i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    CheckForContainerArray(cont, kTrack);
    cont->SetArray(InputEvent());
  }

  // Load all requested cluster branches - each container knows name already
  for (Int_t i =0; i<fClusterCollArray.GetEntriesFast(); i++) {
    AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
    CheckForContainerArray(cont, kCluster);
    cont->SetArray(InputEvent());
  }

  // Determine the proper pointer for each cell object and save them to the cell contianer
  // At this point, they should all be created
  for (auto cellObj : fCellCollArray)
  {
    SetCellsObjectInCellContainerBasedOnProperties(cellObj);
  }

  fEventInitialized = kTRUE;

  ExecOnceComponents();
}

/**
 * Calls ExecOnce() for each component. Additionally, the cells container is added to each component, as
 * this is the first time the pointer to the CaloCells object is available.
 */
void AliEmcalCorrectionTask::ExecOnceComponents()
{
  // Run the initialization for all derived classes.
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  for (auto component : fCorrectionComponents)
  {
    // Setup geomertry
    component->SetEMCALGeometry(fGeom);

    // Add the requested cells to the component
    //AliDebugStream(3) << "Adding CaloCells" << std::endl;
    AddContainersToComponent(component, kCaloCells);
    //AliDebugStream(3) << "Added CaloCells" << std::endl;

    // Component ExecOnce()
    component->ExecOnce();
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
        AliWarning("Could not retrieve centrality information! Assuming 99");
      }
    }
    else { // old centrality estimation < 2015
      AliCentrality *aliCent = InputEvent()->GetCentrality();
      if (aliCent) {
        fCent = aliCent->GetCentralityPercentile(fCentEst.Data());
      }
      else {
        AliWarning("Could not retrieve centrality information! Assuming 99");
      }
    }

    if (fNcentBins==4) {
      if      (fCent >=  0 && fCent <   10) fCentBin = 0;
      else if (fCent >= 10 && fCent <   30) fCentBin = 1;
      else if (fCent >= 30 && fCent <   50) fCentBin = 2;
      else if (fCent >= 50 && fCent <= 100) fCentBin = 3;
      else {
        AliWarning(Form("Negative centrality: %f. Assuming 99", fCent));
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
        AliWarning(Form("Negative centrality: %f. Assuming 99", fCent));
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
        AliWarning(Form("fCentBin too large: cent = %f fCentBin = %d. Assuming 99", fCent, fCentBin));
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
 * Executed each event. It sets run-by-run properties in the correction components and calls Run() for each
 * component.
 */
Bool_t AliEmcalCorrectionTask::Run()
{
  // Run the initialization for all derived classes.
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  for (auto component : fCorrectionComponents)
  {
    component->SetEvent(InputEvent());
    component->SetMCEvent(MCEvent());
    component->SetCentralityBin(fCentBin);
    component->SetCentrality(fCent);

    component->Run();
  }

  PostData(1, fOutput);

  return kTRUE;
}

/**
 * Executed when the file is changed. Also calls UserNotify() for each component.
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

  return kTRUE;
}

/**
 * Write the desired YAML configuration to a file.
 *
 * @param filename The name of the file to write
 * @param userConfig True to write the user configuration
 * @return True when writing the configuration to the file was successful
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
      AliWarning(TString::Format("Configuration not properly initialized! Cannot print %s configuration!", userConfig ? "user" : "default"));
    }

  }
  else
  {
    AliWarning("Please pass a valid filename instead of empty quotes!");
  }
  return returnValue;
}

/**
 * Given a container type, it returns the proper default branch name based on the "usedefault" pattern.
 * This is useful to properly handle creating input objects such as AliEmcalContainer derived objects.
 * If returnObjectType is true, it returns the "default" (unlikely to change) object type instead of the
 * branch name. This is useful to properly determine the type of an object for a TClonesArray.
 *
 * This function can also be very useful in places such as an AddTask(). Using it can significantly reduce
 * code duplication!
 *
 * @param[in] objType Type of the input object
 * @param[in] esdMode True if running with an ESD
 * @param[in] returnObjectType Returns the "default" type of the object rather than the branch name
 *
 * @return The name corresponding to the request branch name or object type.
 */
std::string AliEmcalCorrectionTask::DetermineUseDefaultName(InputObject_t objType, bool esdMode, bool returnObjectType)
{
  std::string returnValue = "";
  if (objType == kCluster) {
    if (esdMode == true) {
      if (returnObjectType == true) {
        returnValue = "AliESDCaloCluster";
      }
      else {
        returnValue = "CaloClusters";
      }
    }
    else {
      if (returnObjectType == true) {
        returnValue = "AliAODCaloCluster";
      }
      else {
        returnValue = "caloClusters";
      }
    }
  }
  else if (objType == kTrack) {
    if (esdMode == true) {
      if (returnObjectType == true) {
        returnValue = "AliESDtrack";
      }
      else {
        returnValue = "Tracks";
      }
    }
    else {
      if (returnObjectType == true) {
        returnValue = "AliAODTrack";
      }
      else {
        returnValue = "tracks";
      }
    }
  }
  else if (objType == kCaloCells) {
    if (esdMode == true) {
      if (returnObjectType == true) {
        returnValue = "AliESDCaloCells";
      }
      else {
        returnValue = "EMCALCells";
      }
    }
    else {
      if (returnObjectType == true) {
        returnValue = "AliAODCaloCells";
      }
      else {
        returnValue = "emcalCells";
      }
    }
  }
  else {
    // Default to empty if we are given an unrecognized type with "usedefault"
    returnValue = "";
  }

  return returnValue;
}

/**
 * Get the proper event based on whether embedding is enabled or not. Useful when determining from which event
 * an input object should be retrieved. It could either be the current input event or an embedded event.
 *
 * @param[in] inputEvent The input event of the analysis. Will be returned if nothing else is requested
 * @param[in] isEmbedding True if the event from embedding should be used.
 *
 * @return The input event to be used
 */
AliVEvent * AliEmcalCorrectionTask::GetEvent(AliVEvent * inputEvent, bool isEmbedding)
{
  AliVEvent * event = 0;
  if (isEmbedding) {
    // TODO: Enable embedded when that branch is committed!
    /*const AliAnalysisTaskEmcalEmbeddingHelper* embedding = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
    if (!embedding) return 0;

    event = embedding->GetExternalEvent();*/
  }
  else {
    event = inputEvent;
  }

  return event;
}

/**
 * Checks if a file exists. This is done inline to make it efficient.
 * See: https://stackoverflow.com/a/19841704
 *
 * @param filename String containing the filename of the file to check.
 *
 * @return True if the file exists.
 */
inline bool AliEmcalCorrectionTask::DoesFileExist(const std::string & filename)
{
  std::ifstream inFile(filename);
  return inFile.good();
}

/**
 * Handles setting up the configuration file to be opened, including in AliPhysics and on the grid.
 * Cannot just use TFile::Open() because the YAML file is just text as opposed to a root file.
 * In the case of a file on the grid, it is copied locally.
 *
 * @param[in] filename Name of the file to be open
 * @param[in] userFile True if the file to be open is a user YAML file
 */
void AliEmcalCorrectionTask::SetupConfigurationFilePath(std::string & filename, bool userFile)
{
  if (filename != "")
  {
    // Handle if in AliPhysics and includes $ALICE_PHYSICS
    filename = gSystem->ExpandPathName(filename.c_str());

    // Handle grid
    if(filename.find("alien://") != std::string::npos)
    {
      AliDebug(2, TString::Format("Opening file \"%s\" on the grid!", filename.c_str()));
      // Initialize alien connection if needed
      if (!gGrid) {
        TGrid::Connect("alien://");
      }

      // Determine the local filename and copy file to local directory
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
 * Uses the information in the cell container to properly set the pointer to the CaloCells object in the
 * cell container.
 *
 * @param[in,out] cellContainer Cell container to set the pointer in
 */
void AliEmcalCorrectionTask::SetCellsObjectInCellContainerBasedOnProperties(AliEmcalCorrectionCellContainer * cellContainer)
{
  AliDebugStream(2) << "Retrieving cells object " << cellContainer->GetName() << std::endl;
  // Check for embedding and return object
  AliVEvent * event = GetEvent(InputEvent(), cellContainer->GetIsEmbedding());

  cellContainer->SetCells(dynamic_cast<AliVCaloCells *>(event->FindListObject(cellContainer->GetBranchName().c_str())));
}

/**
 * Checks whether a container branch exists in the event. If it doesn't exist, then the branch is created
 * automatically. Which this approach requires some care since no fatal error will be thrown when a
 * nonexistent branch is requested, it also allows the creation of output branches just by selecting
 * an unused branch name. For instance, this allows the clusterizer to easily create a new output branch.
 *
 * @param cont The container which is requesting the branch
 * @param objectType The type of the input object
 */
void AliEmcalCorrectionTask::CheckForContainerArray(AliEmcalContainer * cont, InputObject_t objectType)
{
  // TODO: Enable embedded when that branch is committed!
  //AliVEvent * event = GetEvent(InputEvent(), cont->GetIsEmbedding());
  AliVEvent * event = GetEvent(InputEvent());

  TClonesArray *  array = dynamic_cast<TClonesArray *>(event->FindListObject(cont->GetArrayName()));
  if (!array) {
    AliWarning(TString::Format("Container %s requested branch %s, but it does not exist! Creating it for you! Please check that this is the proper action!", cont->GetName(), cont->GetArrayName().Data()));
    array = new TClonesArray(DetermineUseDefaultName(objectType, fIsEsd, true).c_str());
    array->SetName(cont->GetArrayName());
    event->AddObject(array);
  }
}

/**
 * Given the input object type, it return the name of the field in the YAML configuration where information
 * about it should be located. 
 *
 * @param inputObjectType The type of the input object
 * @return The name of the field of the requested input object in the YAML configuration file
 */
std::string AliEmcalCorrectionTask::GetInputFieldNameFromInputObjectType(InputObject_t inputObjectType)
{
  // Get container node
  std::string inputObjectName = "";
  if (inputObjectType == kCluster) {
    inputObjectName = "clusterContainers";
  }
  else if (inputObjectType == kTrack) {
    inputObjectName = "trackContainers";
  }
  else if (inputObjectType == kCaloCells) {
    inputObjectName = "cells";
  }
  else {
    AliFatal(TString::Format("Unrecognized input object type %d", inputObjectType));
  }

  return inputObjectName;
}

/**
 * Checks for a component name in a list of possible component names. This is necessary to search for the
 * components that are associated with a given correction task and it's associated suffix. The comparison is
 * done between strings, so some care is needed not to execute this function too often, especially for a 
 * large number of possible components
 *
 * @param name Name to search for in the possible names
 * @param possibleComponents Possible names of components that have been retrieved from the YAML file
 *
 * @return True name in possible components name
 */
bool AliEmcalCorrectionTask::CheckPossibleNamesForComponentName(std::string & name, std::set <std::string> & possibleComponents)
{
  bool foundComponent = false;
  for (auto & possibleComponent : possibleComponents)
  {
    if (possibleComponent == name) {
      foundComponent = true;
      break;
    }
  }

  return foundComponent;
}

/**
 * Get beam type : pp-AA-pA
 * ESDs have it directly, AODs get it from hardcoded run number ranges
 *
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

/**
 *
 */
void AliEmcalCorrectionTask::PrintRequestedContainersInformation(InputObject_t inputObjectType)
{
  if (inputObjectType == kCaloCells) {
    AliInfoStream() << "Cells info: " << std::endl;
    for (auto cellInfo : fCellCollArray) {
      AliInfoStream() << "\tName: " << cellInfo->GetName() << "\tBranch: " << cellInfo->GetBranchName() << "\tIsEmbedding: " << std::boolalpha << cellInfo->GetIsEmbedding() << std::endl;
    }
  }
  else if (inputObjectType == kCluster || inputObjectType == kTrack) {
    AliInfoStream() << (inputObjectType == kCluster ? "Cluster" : "Track") << " container info: " << std::endl;
    AliEmcalContainer * cont = 0;
    for (auto containerInfo : (inputObjectType == kCluster ? fClusterCollArray : fParticleCollArray) ) {
      cont = static_cast<AliEmcalContainer *>(containerInfo);
      AliInfoStream() << "\tName: " << cont->GetName() << "\tBranch: " << cont->GetArrayName() << "\tTitle: " << cont->GetTitle() << std::endl;
      // TODO: Enable in embedding branch
      //AliInfoStream() << "\tName: " << cont->GetName() << "\tBranch: " << cont->GetArrayName() << "\tTitle: " << cont->GetTitle() << "\tIsEmbedding:" << std::boolalpha << cont->GetIsEmbedding() << std::endl;
    }
  }
  else {
    AliErrorStream() << "Unrecognized input object type " << inputObjectType << std::endl;
  }
}

/**
 * Get the YAML node associated with the named input object. The shared parameters of the passed YAML file
 * is also retrieved and attached to the returned YAML node so that it can be used when retrieving properties.
 *
 * @param[out] inputNode The node that will contain the requested properties
 * @param[in] nodeToRetrieveFrom The node from which the input object properties should be retrieved. Usually the user or default configuration
 * @param[in] inputObjectName Name of the input object node to be retrieved
 * @param[in] requiredProperty True if the input object node is required to exist. It may not if it was not defined in the user configuration.
 */
void AliEmcalCorrectionTask::GetNodeForInputObjects(YAML::Node & inputNode, YAML::Node & nodeToRetrieveFrom, std::string & inputObjectName, bool requiredProperty)
{
  // Get the user input node
  AliEmcalCorrectionComponent::GetProperty(inputObjectName.c_str(), inputNode, YAML::Node(), nodeToRetrieveFrom, requiredProperty, "inputObjects");

  // Get the user shared node and add it back to the user node so that shared parameters are available
  if (nodeToRetrieveFrom["sharedParameters"]) {
    inputNode["sharedParameters"] = nodeToRetrieveFrom["sharedParameters"];
  }
}

/**
 * Utility function for CheckForUnmatchedUserSettings() which returns the names of all of the
 * properties defined in a YAML node. This can then be used to check for consistency in how properties
 * are defined in the user and default configurations.
 *
 * @param[in] componentName Name of the node from which properties are extracted
 * @param[in] node YAML Node of either the user or default configuration
 * @param[out] propertyNames Names of all of the properties that were in the desired node
 * @param[in] nodeRequired True if the node is required to exist
 */
void AliEmcalCorrectionTask::GetPropertyNamesFromNode(const std::string & componentName, const YAML::Node & node, std::set <std::string> & propertyNames, const bool nodeRequired)
{
  YAML::Node tempNode;
  AliEmcalCorrectionComponent::GetProperty(componentName, tempNode, YAML::Node(), node, nodeRequired, "");
  for (auto propertyName : tempNode)
  {
    propertyNames.insert(propertyName.first.as<std::string>());
  }
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
 * Finds the desired cell container by name.
 *
 * @param cellsContainerName Name of the desired cells container
 *
 * @return Pointer to the found cell container, or 0 if not found
 */
AliEmcalCorrectionCellContainer * AliEmcalCorrectionTask::GetCellContainer(const std::string & cellsContainerName) const
{
  for (auto cellContainer : fCellCollArray)
  {
    if (cellContainer->GetName() == cellsContainerName) {
      return cellContainer;
    }
  }

  return 0;
}

AliEmcalCorrectionTask * AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask(TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalCorrectionTask", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalCorrectionTask", "This task requires an input event handler");
    return 0;
  }

  TString name = "AliEmcalCorrectionTask";
  if (suffix != "") {
    name += TString::Format("_%s", suffix.Data());
  }

  AliEmcalCorrectionTask* mgrTask = static_cast<AliEmcalCorrectionTask *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  // Create the task that manages the corrections
  AliEmcalCorrectionTask* correctionTask = new AliEmcalCorrectionTask(name.Data());

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(correctionTask);

  // Create containers for input/output
  AliAnalysisDataContainer* cInput = mgr->GetCommonInputContainer();

  TString outputContainerName(name);
  outputContainerName += "_histos";

  AliAnalysisDataContainer * cOutput = mgr->CreateContainer(outputContainerName.Data(),
                               TList::Class(),
                               AliAnalysisManager::kOutputContainer,
                               Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(correctionTask, 0, cInput);
  mgr->ConnectOutput(correctionTask, 1, cOutput);

  //TObjArray* cnt = mgr->GetContainers();

  return correctionTask;
}
