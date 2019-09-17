/**
 * Implemenation for AliAnalysisTaskEmcalJetHPerformance
 */

#include "AliAnalysisTaskEmcalJetHPerformance.h"

#include <cmath>
#include <map>
#include <vector>
#include <iostream>
#include <bitset>

#include <TObject.h>
#include <TCollection.h>
#include <TAxis.h>
#include <THnSparse.h>

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliLog.h>

#include "yaml-cpp/yaml.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliTrackContainer.h"
#include "AliTLorentzVector.h"
#include "AliClusterContainer.h"
#include "AliEmcalContainerUtils.h"
#include "AliJetContainer.h"
#include "AliEmcalJet.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsQnVector.h"

#include "AliAnalysisTaskEmcalJetHUtils.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance);

namespace PWGJE {
namespace EMCALJetTasks {

/**
 * Default constructor
 */
AliAnalysisTaskEmcalJetHPerformance::AliAnalysisTaskEmcalJetHPerformance():
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetHPerformance", kFALSE),
  fYAMLConfig(),
  fConfigurationInitialized(false),
  fHistManager(),
  fEmbeddingQA(),
  fCreateQAHists(false),
  fCreateCellQAHists(false),
  fPerformJetMatching(false),
  fCreateResponseMatrix(false),
  fEmbeddedCellsName("emcalCells"),
  fPreviousEventTrigger(0),
  fPreviousEmbeddedEventSelected(false),
  fEfficiencyPeriodIdentifier(AliAnalysisTaskEmcalJetHUtils::kDisableEff),
  fFlowQnVectorManager(nullptr),
  fResponseMatrixFillMap(),
  fResponseFromThreeJetCollections(true),
  fMinFractionShared(0.),
  fLeadingHadronBiasType(AliAnalysisTaskEmcalJetHUtils::kCharged)
{
}

/**
 * Standard constructor
 */
AliAnalysisTaskEmcalJetHPerformance::AliAnalysisTaskEmcalJetHPerformance(const char * name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fYAMLConfig(),
  fConfigurationInitialized(false),
  fHistManager(name),
  fEmbeddingQA(),
  fCreateQAHists(false),
  fCreateCellQAHists(false),
  fPerformJetMatching(false),
  fCreateResponseMatrix(false),
  fEmbeddedCellsName("emcalCells"),
  fPreviousEventTrigger(0),
  fPreviousEmbeddedEventSelected(false),
  fEfficiencyPeriodIdentifier(AliAnalysisTaskEmcalJetHUtils::kDisableEff),
  fFlowQnVectorManager(nullptr),
  fResponseMatrixFillMap(),
  fResponseFromThreeJetCollections(true),
  fMinFractionShared(0.),
  fLeadingHadronBiasType(AliAnalysisTaskEmcalJetHUtils::kCharged)
{
  // Ensure that additional general histograms are created
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Copy constructor
 */
AliAnalysisTaskEmcalJetHPerformance::AliAnalysisTaskEmcalJetHPerformance(const AliAnalysisTaskEmcalJetHPerformance & other):
  fYAMLConfig(other.fYAMLConfig),
  fConfigurationInitialized(other.fConfigurationInitialized),
  fHistManager(other.fHistManager.GetName()),
  //fEmbeddingQA(other.fEmbeddingQA), // Cannot use because the THistManager (which is in the class) copy constructor is private.
  fCreateQAHists(other.fCreateQAHists),
  fCreateCellQAHists(other.fCreateCellQAHists),
  fPerformJetMatching(other.fPerformJetMatching),
  fCreateResponseMatrix(other.fCreateResponseMatrix),
  fEmbeddedCellsName(other.fEmbeddedCellsName),
  fPreviousEventTrigger(other.fPreviousEventTrigger),
  fPreviousEmbeddedEventSelected(other.fPreviousEmbeddedEventSelected),
  fEfficiencyPeriodIdentifier(other.fEfficiencyPeriodIdentifier),
  fFlowQnVectorManager(other.fFlowQnVectorManager),
  fResponseMatrixFillMap(other.fResponseMatrixFillMap),
  fResponseFromThreeJetCollections(other.fResponseFromThreeJetCollections),
  fMinFractionShared(other.fMinFractionShared),
  fLeadingHadronBiasType(other.fLeadingHadronBiasType)
{
  TIter next(other.fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fHistManager.SetObject(obj, obj->GetName());
  }
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
AliAnalysisTaskEmcalJetHPerformance & AliAnalysisTaskEmcalJetHPerformance::operator=(AliAnalysisTaskEmcalJetHPerformance other)
{
  swap(*this, other);
  return *this;
}

/**
 * Utility function to apply the determine the single track efficiency.
 *
 * @param trackEta Eta of the track
 * @param trackPt pT of the track
 *
 * @return Track efficiency of the track (the entry in a histogram should be weighted as 1/(return value))
 */
double AliAnalysisTaskEmcalJetHPerformance::DetermineTrackingEfficiency(double trackPt, double trackEta)
{
  return AliAnalysisTaskEmcalJetHUtils::DetermineTrackingEfficiency(
   trackPt, trackEta, fCentBin, fEfficiencyPeriodIdentifier,
   "PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance");
}

/**
 * Retrieve task properties from the YAML configuration.
 */
void AliAnalysisTaskEmcalJetHPerformance::RetrieveAndSetTaskPropertiesFromYAMLConfig()
{
  // Base class options
  // Task physics (trigger) selection.
  std::vector<std::string> physicsSelection;
  bool res = fYAMLConfig.GetProperty(std::vector<std::string>({"eventCuts", "physicsSelection"}), physicsSelection, false);
  if (res) {
    fOfflineTriggerMask = AliEmcalContainerUtils::DeterminePhysicsSelectionFromYAML(physicsSelection);
  }

  // Same ordering as in the constructor (for consistency)
  std::string baseName = "enable";
  // These are disabled by default.
  fYAMLConfig.GetProperty({baseName, "QAHists"}, fCreateQAHists, false);
  fYAMLConfig.GetProperty({baseName, "CellQAHists"}, fCreateCellQAHists, false);
  fYAMLConfig.GetProperty({baseName, "jetMatching"}, fPerformJetMatching, false);
  fYAMLConfig.GetProperty({baseName, "responseMatrix"}, fCreateResponseMatrix, false);

  // Event cuts
  baseName = "eventCuts";
  // If event cuts are enabled (which they exceptionally are by default), then we want to configure them here.
  // If the event cuts are explicitly disabled, then we invert that value to enable the AliAnylsisTaskEmcal
  // builtin event selection.
  bool tempBool;
  fYAMLConfig.GetProperty({baseName, "enabled"}, tempBool, false);
  fUseBuiltinEventSelection = !tempBool;
  if (fUseBuiltinEventSelection == false) {
    // Need to include the namespace so that AliDebug will work properly...
    std::string taskName = "PWGJE::EMCALJetTasks::";
    taskName += GetName();
    AliAnalysisTaskEmcalJetHUtils::ConfigureEventCuts(fAliEventCuts, fYAMLConfig, fOfflineTriggerMask, baseName, taskName);
  }

  // General task options
  baseName = "general";
  fYAMLConfig.GetProperty({baseName, "nCentBins"}, fNcentBins, false);

  // QA options
  baseName = "QA";
  fYAMLConfig.GetProperty({baseName, "cellsName"}, fCaloCellsName, false);
  // Defaults to "emcalCells" if not set.
  fYAMLConfig.GetProperty({baseName, "embeddedCellsName"}, fEmbeddedCellsName, false);
  // Efficiency
  std::string tempStr = "";
  baseName = "efficiency";
  res = fYAMLConfig.GetProperty({baseName, "periodIdentifier"}, tempStr, false);
  if (res) {
    fEfficiencyPeriodIdentifier = AliAnalysisTaskEmcalJetHUtils::fgkEfficiencyPeriodIdentifier.at(tempStr);
  }

  // Jet matching
  baseName = "jetMatching";
  fYAMLConfig.GetProperty({baseName, "maxMatchingDistance"}, fMaxJetMatchingDistance, false);

  // Response matrix properties
  baseName = "responseMatrix";
  fYAMLConfig.GetProperty({baseName, "useThreeJetCollections"}, fResponseFromThreeJetCollections, false);
  fYAMLConfig.GetProperty({baseName, "minFractionSharedPt"}, fMinFractionShared, false);
  std::string hadronBiasStr = "";
  baseName = "jets";
  res = fYAMLConfig.GetProperty({baseName, "leadingHadronBiasType"}, hadronBiasStr, false);
  // Only attempt to set the property if it is retrieved successfully
  if (res) {
    fLeadingHadronBiasType = AliAnalysisTaskEmcalJetHUtils::fgkLeadingHadronBiasMap.at(hadronBiasStr);
  }
}

/**
 * Create jet containers based on the jet values defined in the YAML config.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupJetContainersFromYAMLConfig()
{
  std::string baseName = "jets";
  std::vector <std::string> jetNames = {"hybridLevelJets", "detLevelJets", "partLevelJets", "analysisJets"};
  for (const auto & jetName : jetNames) {
    // Retrieve the node just to see if it is exists. If so, then we can proceed
    YAML::Node node;
    bool res = fYAMLConfig.GetProperty(std::vector<std::string>({baseName, jetName}), node, false);
    if (res) {
      // Retrieve jet container properties
      std::string collectionName = "";
      double R = -1;
      fYAMLConfig.GetProperty({baseName, jetName, "collection"}, collectionName, true);
      fYAMLConfig.GetProperty({baseName, jetName, "R"}, R, true);

      // Determine the jet acceptance.
      std::vector<std::string> jetAcceptances;
      fYAMLConfig.GetProperty({baseName, jetName, "acceptance"}, jetAcceptances, true);
      UInt_t jetAcceptance = AliAnalysisTaskEmcalJetHUtils::DetermineJetAcceptanceFromYAML(jetAcceptances);

      // Create jet container and set the name
      AliDebugStream(1) << "Creating jet from jet collection name " << collectionName << " with acceptance "
               << jetAcceptance << " and R=" << R << "\n";
      AliJetContainer* jetCont = AddJetContainer(collectionName.c_str(), jetAcceptance, R);
      jetCont->SetName(jetName.c_str());

      // Jet area cut percentage
      double jetAreaCut = -1;
      bool res = fYAMLConfig.GetProperty({ baseName, jetName, "areaCutPercentage" }, jetAreaCut, false);
      if (res) {
        AliDebugStream(1) << "Setting jet area cut percentage of " << jetAreaCut << " for jet cont " << jetName
                 << "\n";
        jetCont->SetPercAreaCut(jetAreaCut);
      }

      // Rho name
      std::string tempStr = "";
      res = fYAMLConfig.GetProperty({ baseName, jetName, "rhoName" }, tempStr, false);
      if (res) {
        AliDebugStream(1) << "Setting rho name of " << tempStr << " for jet cont " << jetName << "\n";
        jetCont->SetRhoName(tempStr.c_str());
      }

      // Leading hadron type
      int leadingHadronType = -1;
      res = fYAMLConfig.GetProperty({ baseName, jetName, "leadingHadronType" }, leadingHadronType, false);
      if (res) {
        AliDebugStream(1) << "Setting leading hadron type of " << leadingHadronType << " for jet cont "
                 << jetName << "\n";
        jetCont->SetLeadingHadronType(leadingHadronType);
      }

      // Leading hadron cut
      double leadingHadronBias = 0.;
      res = fYAMLConfig.GetProperty({ baseName, jetName, "leadingHadronBias" }, leadingHadronBias, false);
      if (res) {
        AliDebugStream(1) << "Setting leading hadron bias of " << leadingHadronBias << " for jet cont "
                 << jetName << "\n";
        jetCont->SetMinTrackPt(leadingHadronBias);
      }

      // Max track pt
      double maxTrackPt = 0.;
      res = fYAMLConfig.GetProperty({ baseName, jetName, "maxTrackPt" }, maxTrackPt, false);
      if (res) {
        AliDebugStream(1) << "Setting max track pt of " << maxTrackPt << " for jet cont "
                 << jetName << "\n";
        jetCont->SetMinTrackPt(maxTrackPt);
      }
    }
    else {
      AliInfoStream() << "Unable to find definition of jet container corresponding to \"" << jetName << "\"\n";
    }
  }
}

/**
 * Setup particle containers based on the provided YAML configuration.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupParticleContainersFromYAMLConfig()
{
  std::string baseName = "particles";
  // Retrieve the node just to see if it is exists. If so, then we can proceed
  YAML::Node node;
  fYAMLConfig.GetProperty(baseName, node);
  // Iterate over all of the particle and track containers
  for (const auto & n : node) {
    std::string containerName = n.first.as<std::string>();

    // Create the container.
    std::string branchName = "";
    fYAMLConfig.GetProperty({baseName, containerName, "branchName"}, branchName, true);
    if (branchName == "usedefault") {
      branchName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack);
    }
    AliParticleContainer * partCont = AliAnalysisTaskEmcalJetHUtils::CreateParticleOrTrackContainer(branchName);
    this->AdoptParticleContainer(partCont);

    // Configure the container
    // Need to include the namespace so that AliDebug will work properly...
    std::string taskName = "PWGJE::EMCALJetTasks::";
    taskName += GetName();
    std::vector<std::string> baseNameWithContainer = { baseName, containerName };
    AliAnalysisTaskEmcalJetHUtils::ConfigureEMCalContainersFromYAMLConfig(baseNameWithContainer, containerName,
                                       partCont, fYAMLConfig, taskName);

    // Track specific properties
    AliTrackContainer * trackCont = dynamic_cast<AliTrackContainer *>(partCont);
    if (trackCont) {
      AliAnalysisTaskEmcalJetHUtils::ConfigureTrackContainersFromYAMLConfig(baseNameWithContainer, trackCont,
                                         fYAMLConfig, taskName);
    }

    AliDebugStream(2) << "Particle/track container: " << partCont->GetName()
             << ", array class: " << partCont->GetClassName() << ", collection name: \""
             << partCont->GetArrayName() << "\", min pt: " << partCont->GetMinPt() << "\n";
  }
}

/**
 * Setup cluster containers based on the provided YAML configuration.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupClusterContainersFromYAMLConfig()
{
  std::string baseName = "clusters";
  // Retrieve the node just to see if it is exists. If so, then we can proceed
  YAML::Node node;
  fYAMLConfig.GetProperty(baseName, node);
  // Iterate over all of the cluster containers
  for (const auto & n : node) {
    std::string containerName = n.first.as<std::string>();

    // Create the container.
    std::string branchName = "";
    fYAMLConfig.GetProperty({baseName, containerName, "branchName"}, branchName, true);
    if (branchName == "usedefault") {
      branchName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster);
    }
    AliClusterContainer * clusterCont = this->AddClusterContainer(branchName.c_str());

    // Configure the container
    // Need to include the namespace so that AliDebug will work properly...
    std::string taskName = "PWGJE::EMCALJetTasks::";
    taskName += GetName();
    std::vector<std::string> baseNameWithContainer = { baseName, containerName };
    AliAnalysisTaskEmcalJetHUtils::ConfigureEMCalContainersFromYAMLConfig(baseNameWithContainer, containerName,
                                       clusterCont, fYAMLConfig, taskName);

    // Cluster specific properties
    AliAnalysisTaskEmcalJetHUtils::ConfigureClusterContainersFromYAMLConfig(baseNameWithContainer, clusterCont, fYAMLConfig, taskName);

    AliDebugStream(2) << "Cluster container: " << clusterCont->GetName()
             << ", array class: " << clusterCont->GetClassName() << ", collection name: \""
             << clusterCont->GetArrayName() << "\", min pt: " << clusterCont->GetMinPt() << "\n";
  }
}

/**
 * Initialize task.
 */
bool AliAnalysisTaskEmcalJetHPerformance::Initialize()
{
  fConfigurationInitialized = false;

  // Ensure that we have at least one configuration in the YAML config.
  if (fYAMLConfig.DoesConfigurationExist(0) == false) {
    // No configurations exist. Return immediately.
    return fConfigurationInitialized;
  }

  // Always initialize for streaming purposes
  fYAMLConfig.Initialize();

  // Setup task based on the properties defined in the YAML config
  AliDebugStream(2) << "Configuring task from the YAML configuration.\n";
  RetrieveAndSetTaskPropertiesFromYAMLConfig();
  SetupJetContainersFromYAMLConfig();
  SetupParticleContainersFromYAMLConfig();
  SetupClusterContainersFromYAMLConfig();
  AliDebugStream(2) << "Finished configuring via the YAML configuration.\n";

  // Print the results of the initialization
  // Print outside of the ALICE Log system to ensure that it is always available!
  std::cout << *this;

  fConfigurationInitialized = true;
  return fConfigurationInitialized;
}

/**
 * Create output objects
 */
void AliAnalysisTaskEmcalJetHPerformance::UserCreateOutputObjects()
{
  // First call the base class
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Check that the task was initialized
  if (!fConfigurationInitialized) {
    AliFatal("Task was not initialized. Please ensure that Initialize() was called!");
  }
  // Reinitialize the YAML configuration
  fYAMLConfig.Reinitialize();

  // Create histograms
  if (fCreateQAHists) {
    SetupQAHists();
  }
  if (fPerformJetMatching) {
    SetupJetMatchingQA();
  }
  if (fCreateResponseMatrix) {
    SetupResponseMatrixHists();
  }

  // Store hist manager output in the output list
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }

  // Initialize
  const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (embeddingHelper) {
    bool res = fEmbeddingQA.Initialize();
    if (res) {
      fEmbeddingQA.AddQAPlotsToList(fOutput);
    }
  }

  // Post the data when finished
  PostData(1, fOutput);
}

/**
 * Setup cell QA histograms at a specific prefix within the histogram manager. This allows us to avoid repeating
 * definitions for the same histograms in the embedded vs internal event. This function actually defines the histograms.
 *
 * @param[in] prefix Prefix under which the histograms will be created within the hist manager.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupCellQAHistsWithPrefix(const std::string & prefix)
{
  std::string name = prefix + "/fHistCellEnergy";
  std::string title = name + ";\\mathit{E}_{\\mathrm{cell}} (\\mathrm{GeV});counts";
  fHistManager.CreateTH1(name.c_str(), title.c_str(), 300, 0, 150);

  name = prefix + "/fHistCellTime";
  title = name + ";t (s);counts";
  fHistManager.CreateTH1(name.c_str(), title.c_str(), 1000, -10e-6, 10e-6);

  name = prefix + "/fHistNCells";
  title = name + ";\\mathrm{N}_{\\mathrm{cells}};counts";
  fHistManager.CreateTH1(name.c_str(), title.c_str(), 100, 0, 5000);

  name = prefix + "/fHistCellID";
  title = name + ";\\mathrm{N}_{\\mathrm{cell}};counts";
  fHistManager.CreateTH1(name.c_str(), title.c_str(), 20000, 0, 20000);

  // Histograms for embedding QA which use the cell timing to determine whether the
  // embedded event has been double corrected.
  if (prefix.find("embedding") != std::string::npos) {
    name = prefix + "/fHistEmbeddedEventUsed";
    title = name + ";Embedded event used";
    fHistManager.CreateTH1(name.c_str(), title.c_str(), 2, 0, 2);

    name = prefix + "/fHistInternalEventSelection";
    title = name + ";Embedded event used;Trigger bit";
    fHistManager.CreateTH2(name.c_str(), title.c_str(), 2, 0, 2, 32, -0.5, 31.5);
  }
}

/**
 * Directs the creation of Cell QA histograms for the internal and external events (if it exists).
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupCellQAHists()
{
  // Only create and fill the cells QA if we've setup the cells name.
  if (fCaloCellsName != "") {
    std::string prefix = "QA/";
    prefix += fCaloCellsName.Data();
    SetupCellQAHistsWithPrefix(prefix);
    auto embeddingInstance = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
    if (embeddingInstance) {
      prefix = "QA/embedding/";
      prefix += fEmbeddedCellsName;
      SetupCellQAHistsWithPrefix(prefix);
    }
  }
}

/**
 * Setup and allocate histograms related to QA histograms.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupQAHists()
{
  // Cell level QA
  if (fCreateCellQAHists) {
    SetupCellQAHists();
  }

  // Tracks
  AliTrackContainer * trackCont = nullptr;
  TIter nextTrackColl(&fParticleCollArray);
  while ((trackCont = static_cast<AliTrackContainer*>(nextTrackColl()))) {
    std::string name = "QA/%s/fHistTrackPtEtaPhi";
    std::string title = name + ";#it{p}_{T} (GeV);#eta;#phi";
    fHistManager.CreateTH3(TString::Format(name.c_str(), trackCont->GetName()),
                TString::Format(title.c_str(), trackCont-> GetName()), 50, 0, 25, 40, -1, 1, 72, 0,
                TMath::TwoPi());
    name = "QA/%s/fHistTrackPtEtaPhiEfficiencyCorrected";
    title = name + ";#it{p}_{T} (GeV);#eta;#phi";
    fHistManager.CreateTH3(TString::Format(name.c_str(), trackCont->GetName()),
                TString::Format(title.c_str(), trackCont-> GetName()), 50, 0, 25, 40, -1, 1, 72, 0,
                TMath::TwoPi());
  }

  // Clusters
  AliClusterContainer * clusterCont = nullptr;
  TIter nextClusColl(&fClusterCollArray);
  while ((clusterCont = static_cast<AliClusterContainer*>(nextClusColl()))) {
    // Cluster time vs energy
    std::string name = "QA/%s/fHistClusterEnergyVsTime";
    std::string title = name + ";#it{E}_{cluster} (GeV);t_{cluster} (s);counts";
    fHistManager.CreateTH2(TString::Format(name.c_str(), clusterCont->GetName()),
                TString::Format(title.c_str(), clusterCont->GetName()), 1000, 0, 100, 300, -300e-9,
                300e-9);
    // Hadronically corrected energy (which is usually what we're using)
    name = "QA/%s/fHistClusterHadCorrEnergy";
    title = name + ";#it{E}_{cluster}^{had.corr.} (GeV);counts";
    fHistManager.CreateTH1(TString::Format(name.c_str(), clusterCont->GetName()), TString::Format(title.c_str(), clusterCont->GetName()), 200, 0, 100);
    // Cluster eta-phi
    name = "QA/%s/fHistClusterEtaPhi";
    title = name + ";#eta_{cluster};#phi_{cluster};counts";
    fHistManager.CreateTH2(TString::Format(name.c_str(), clusterCont->GetName()), TString::Format(title.c_str(), clusterCont->GetName()), 28, -0.7, 0.7, 72, 0, TMath::TwoPi());
  }

  // Jets
  AliJetContainer * jetCont = nullptr;
  TIter nextJetColl(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(nextJetColl()))) {
    // Jet pT
    std::string name = "QA/%s/fHistJetPt";
    std::string title = name + ";p_{T} (GeV)";
    fHistManager.CreateTH1(TString::Format(name.c_str(), jetCont->GetName()), TString::Format(title.c_str(), jetCont->GetName()), 600, -50, 250);
  }

  // Event plane resolution
  // Only enable if the Qn vector corrections task is available.
  auto qnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*>(
   AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  if (qnVectorTask) {
    // Setup the Qn corrections manager.
    fFlowQnVectorManager = qnVectorTask->GetAliQnCorrectionsManager();
    // Then setup the resolution histograms.
    // Resolution are calculated via TProfiles as a function of centrality.
    std::map<std::string, std::string> resolutionNamesToTitles = {
      std::make_pair("VZERO", R"(VZERO: $<\cos(R * (\Psi_{\mathrm{TPC_{Pos}}} - \Psi_{\mathrm{TPC_{Neg}}}))>$)"),
      std::make_pair("TPC_Pos", R"(TPC $\eta$ > 0: $<\cos(R * (\Psi_{\mathrm{V0M}} - \Psi_{\mathrm{TPC_{Neg}}}))>$)"),
      std::make_pair("TPC_Neg", R"(TPC $\eta$ < 0: $<\cos(R * (\Psi_{\mathrm{V0M}} - \Psi_{\mathrm{TPC_{Pos}}}))>$)")
    };
    for (const auto& nameAndTitle: resolutionNamesToTitles) {
      for (unsigned int R = 2; R < 9; R++) {
        std::string name = "QA/eventPlaneRes/%s/R%i";
        std::string title = "%s event plane res;Centrality (%%);R_{%i}";
        fHistManager.CreateTProfile(TString::Format(name.c_str(), nameAndTitle.first.c_str(), R),
                      TString::Format(title.c_str(), nameAndTitle.second.c_str(), R),
                      10, 0, 100);
      }
    }
  }
}

/**
 * Setup and allocate histograms related to jet matching.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupJetMatchingQA()
{
  // Setup
  const int nBinsPt        = 40;
  const int nBinsDPhi      = 72;
  const int nBinsDEta      = 100;
  const int nBinsDR        = 50;
  const int nBinsFraction  = 101;

  const double minPt       = -50.;
  const double maxPt       = 150.;
  const double minDPhi     = -0.5;
  const double maxDPhi     =  0.5;
  const double minDEta     = -0.5;
  const double maxDEta     =  0.5;
  const double minDR       =  0.;
  const double maxDR       =  0.5;
  const double minFraction =  -0.005;
  const double maxFraction =  1.005;

  // Create the histograms
  // "s" ensures that Sumw2() is called
  std::string name = "";
  std::string title = "";
  std::vector<std::string> matchingTypes = {"hybridToDet", "detToPart"};
  for (auto matchingType : matchingTypes)
  {
    for (unsigned int centBin = 0; centBin < fNcentBins; centBin++)
    {
      // Jet 1 pt vs delta eta vs delta phi
      name = "jetMatching/%s/fh3PtJet1VsDeltaEtaDeltaPhi_%d";
      title = "fh3PtJet1VsDeltaEtaDeltaPhi_%d;#it{p}_{T,jet1};#it{#Delta#eta};#it{#Delta#varphi}";
      fHistManager.CreateTH3(TString::Format(name.c_str(), matchingType.c_str(), centBin),
                  TString::Format(title.c_str(), centBin), nBinsPt, minPt, maxPt, nBinsDEta,
                  minDEta, maxDEta, nBinsDPhi, minDPhi, maxDPhi, "s");

      // Jet pt 1 vs delta R
      name = "jetMatching/%s/fh2PtJet1VsDeltaR_%d";
      title = "fh2PtJet1VsDeltaR_%d;#it{p}_{T,jet1};#it{#Delta R}";
      fHistManager.CreateTH2(TString::Format(name.c_str(), matchingType.c_str(), centBin),
                  TString::Format(title.c_str(), centBin), nBinsPt, minPt, maxPt, nBinsDR, minDR,
                  maxDR, "s");

      // Jet pt 2 vs shared momentum fraction
      name = "jetMatching/%s/fh2PtJet2VsFraction_%d";
      title = "fh2PtJet2VsFraction_%d;#it{p}_{T,jet2};#it{f}_{shared}";
      fHistManager.CreateTH2(TString::Format(name.c_str(), matchingType.c_str(), centBin),
                  TString::Format(title.c_str(), centBin), nBinsPt, minPt, maxPt, nBinsFraction,
                  minFraction, maxFraction, "s");

      // Jet pt 1 vs lead pt before checking if the jet has a matched jet.
      name = "jetMatching/%s/fh2PtJet1VsLeadPtAllSel_%d";
      title = "fh2PtJet1VsLeadPtAllSel_%d;#it{p}_{T,jet1};#it{p}_{T,lead trk}";
      fHistManager.CreateTH2(TString::Format(name.c_str(), matchingType.c_str(), centBin),
                  TString::Format(title.c_str(), centBin), nBinsPt, minPt, maxPt, 20, 0., 20.,
                  "s");

      // Jet pt 1 vs lead pt after checking if the jet has a matched jet.
      name = "jetMatching/%s/fh2PtJet1VsLeadPtTagged_%d";
      title = "fh2PtJet1VsLeadPtTagged_%d;#it{p}_{T,jet1};#it{p}_{T,lead trk}";
      fHistManager.CreateTH2(TString::Format(name.c_str(), matchingType.c_str(), centBin),
                  TString::Format(title.c_str(), centBin), nBinsPt, minPt, maxPt, 20, 0., 20.,
                  "s");

      // Jet pt 1 vs jet pt 2.
      name = "jetMatching/%s/fh2PtJet1VsPtJet2_%d";
      title = "fh2PtJet1VsPtJet2_%d;#it{p}_{T,jet1};#it{p}_{T,jet2}";
      fHistManager.CreateTH2(TString::Format(name.c_str(), matchingType.c_str(), centBin),
                  TString::Format(title.c_str(), centBin), nBinsPt, minPt, maxPt, nBinsPt, minPt,
                  maxPt, "s");

      // Jet pt2 residual (aka. the jet energy scale).
      name = "jetMatching/%s/fh2PtJet2VsRelPt_%d";
      title = "fh2PtJet2VsRelPt_%d;#it{p}_{T,jet2};(#it{p}_{T,jet1}-#it{p}_{T,jet2})/#it{p}_{T,jet2}";
      fHistManager.CreateTH2(TString::Format(name.c_str(), matchingType.c_str(), centBin),
                  TString::Format(title.c_str(), centBin), nBinsPt, minPt, maxPt, 241, -2.41, 2.41,
                  "s");
    }

    // Number of jets accepted.
    name = "jetMatching/%s/fNAccJets";
    title = "fNAccJets;N/ev";
    fHistManager.CreateTH1(TString::Format(name.c_str(), matchingType.c_str()), title.c_str(), 11, -0.5, 9.5, "s");
  }

}

/**
 * Setup and allocate histograms related to creating a response matrix.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupResponseMatrixHists()
{
  // Main response matrix THnSparse
  std::string name = "response/fHistResponseMatrix";
  std::string title = name;

  // Retrieve binning from the YAML configuration
  std::vector<TAxis *> binning;
  // This structure is designed to preserve the order of the axis in the YAML by using a YAML sequence (decoded into
  // a vector), while defining a pair of the axis name and axis limts. Using this structure avoids the need to create
  // a new object and conversion to retrieve the data
  std::vector<std::pair<std::string, std::vector<double>>> sparseAxes;
  std::string baseName = "responseMatrix";
  fYAMLConfig.GetProperty({baseName, "axes"}, sparseAxes, true);
  for (auto axis : sparseAxes) {
    auto axisLimits = axis.second;
    AliDebugStream(3) << "Creating axis " << axis.first << " with nBins " << axisLimits.at(0) << ", min: " << axisLimits.at(1) << ", max: " << axisLimits.at(2) << "\n";
    binning.emplace_back(new TAxis(axisLimits.at(0), axisLimits.at(1), axisLimits.at(2)));
  }

  // "s" ensures that Sumw2() is called
  // The explicit const_cast is required
  THnSparse * hist = fHistManager.CreateTHnSparse(name.c_str(), title.c_str(), binning.size(), const_cast<const TAxis **>(binning.data()), "s");
  // Set the axis titles
  int axisNumber = 0;
  for (auto axis = sparseAxes.begin(); axis != sparseAxes.end(); axis++) {
    AliDebugStream(5) << "ResponseMatrix: Add axis " << axis->first << " to sparse\n";
    hist->GetAxis(axisNumber)->SetTitle(axis->first.c_str());
    axisNumber++;
  }

  // Define mapping from value name to value (including jet information)
  for (unsigned int iJet = 1; iJet < 3; iJet++)
  {
    // pt
    // ex: p_{T,1}
    fResponseMatrixFillMap.insert(std::make_pair(std::string("p_{T,") + std::to_string(iJet) + "}", std::make_pair(iJet, &ResponseMatrixFillWrapper::fPt)));
    // Area
    // ex: A_{jet,1}
    fResponseMatrixFillMap.insert(std::make_pair(std::string("A_{jet,") + std::to_string(iJet) + "}", std::make_pair(iJet, &ResponseMatrixFillWrapper::fArea)));
    // EP angle
    // ex: #theta_{jet,1}^{EP}
    fResponseMatrixFillMap.insert(std::make_pair(std::string("#theta_{jet,") + std::to_string(iJet) + "}^{EP}", std::make_pair(iJet, &ResponseMatrixFillWrapper::fRelativeEPAngle)));
    // Leading hadron
    // ex: p_{T,particle,1}^{leading} (GeV/c)
    fResponseMatrixFillMap.insert(std::make_pair(std::string("p_{T,particle,") + std::to_string(iJet) + "}^{leading} (GeV/c)", std::make_pair(iJet, &ResponseMatrixFillWrapper::fLeadingHadronPt)));
  }
  // Distance from one jet to another
  fResponseMatrixFillMap.insert(std::make_pair("distance", std::make_pair(1, &ResponseMatrixFillWrapper::fDistance)) );
  // Centrality
  fResponseMatrixFillMap.insert(std::make_pair("centrality", std::make_pair(1, &ResponseMatrixFillWrapper::fCentrality)) );

  // Shared momentum fraction
  name = "fHistFractionSharedPt";
  title = "Fraction of p_{T} shared between matched jets";
  // Need to include the bin from 1-1.01 to ensure that jets which shared all of their momentum
  // due not end up in the overflow bin!
  fHistManager.CreateTH1(name.c_str(), title.c_str(), 101, 0, 1.01);
}

/**
 * Main event loop
 */
Bool_t AliAnalysisTaskEmcalJetHPerformance::Run()
{
  // Only fill the embedding qa plots if:
  //  - We are using the embedding helper
  //  - The class has been initialized
  if (fEmbeddingQA.IsInitialized()) {
    fEmbeddingQA.RecordEmbeddedEventProperties();
  }

  // QA
  if (fCreateQAHists) {
    QAHists();
  }

  // Jet matching
  if (fPerformJetMatching) {
    // Setup
    AliDebugStream(1) << "Performing jet matching\n";
    // We don't perform any additional jet matching initialization because we will restrict
    // the jets accepted via the jet acceptance cuts (ie. EMCal fiducial cuts)
    // Retrieve the releveant jet collections
    AliJetContainer * jetsHybrid = GetJetContainer("hybridLevelJets");
    AliJetContainer * jetsDetLevel = GetJetContainer("detLevelJets");
    AliJetContainer * jetsPartLevel = GetJetContainer("partLevelJets");
    // Validation
    if (!jetsHybrid) {
      AliErrorStream() << "Could not retrieve hybrid jet collection.\n";
      return kFALSE;
    }
    if (!jetsDetLevel) {
      AliErrorStream() << "Could not retrieve det level jet collection.\n";
      return kFALSE;
    }
    if (!jetsPartLevel) {
      AliErrorStream() << "Could not retrieve part level jet collection.\n";
      return kFALSE;
    }

    // Now, begin the actual matching.
    // Hybrid <-> det first
    AliDebugStream(2) << "Matching hybrid to detector level jets.\n";
    // First, we reset the tagging
    ResetMatching(*jetsHybrid);
    ResetMatching(*jetsDetLevel);
    // Next, we perform the matching
    PerformGeometricalJetMatching(*jetsHybrid, *jetsDetLevel, fMaxJetMatchingDistance);
    // Now, begin the next matching stage
    // det <-> particle
    AliDebugStream(2) << "Matching detector level to particle level jets.\n";
    // First, we reset the tagging. We need to reset the det matching again to ensure
    // that it doesn't accidentally keep some latent matches to the hybrid jets.
    ResetMatching(*jetsDetLevel);
    ResetMatching(*jetsPartLevel);
    // Next, we perform the matching
    PerformGeometricalJetMatching(*jetsDetLevel, *jetsPartLevel, fMaxJetMatchingDistance);

    // Fill QA hists.
    FillJetMatchingQA(*jetsHybrid, *jetsDetLevel, "hybridToDet");
    FillJetMatchingQA(*jetsDetLevel, *jetsPartLevel, "detToPart");
  }

  // Response matrix
  if (fCreateResponseMatrix) {
    ResponseMatrix();
  }

  return kTRUE;
}

/**
 * Process and fill QA hists
 */
void AliAnalysisTaskEmcalJetHPerformance::QAHists()
{
  // Continue on to filling the histograms
  FillQAHists();
}

/**
 * Helper function to fill cell QA into the defined histograms. This actually fills the histograms at a given prefix
 * from the provided calo cells.
 *
 * @param[in] prefix Prefix under which the histograms will be created within the hist manager.
 * @param[in] cells Cells from which the information should be extracted.
 */
void AliAnalysisTaskEmcalJetHPerformance::FillCellQAHists(const std::string & prefix, AliVCaloCells * cells)
{
  AliDebugStream(4) << "Storing cells with prefix \"" << prefix << "\". N cells: " << cells->GetNumberOfCells() << "\n";
  short absId = -1;
  double eCell = 0;
  double tCell = 0;
  double eFrac = 0;
  int mcLabel = -1;

  std::string energyName = prefix + "/fHistCellEnergy";
  std::string timeName = prefix + "/fHistCellTime";
  std::string idName = prefix + "/fHistCellID";
  bool embeddedCellWithLateCellTime = false;
  bool fillingEmbeddedCells = (prefix.find("embedding") != std::string::npos);
  for (unsigned int iCell = 0; iCell < cells->GetNumberOfCells(); iCell++) {
    cells->GetCell(iCell, absId, eCell, tCell, mcLabel, eFrac);

    AliDebugStream(5) << "Cell " << iCell << ": absId: " << absId << ", E: " << eCell << ", t: " << tCell
             << ", mcLabel: " << mcLabel << ", eFrac: " << eFrac << "\n";
    fHistManager.FillTH1(energyName.c_str(), eCell);
    fHistManager.FillTH1(timeName.c_str(), tCell);
    fHistManager.FillTH1(idName.c_str(), absId);

    // We will record the event selection if the time is less than -400 ns
    // This corresponds to a doubly corrected embedded event, which shouldn't be possible, and therefore
    // indicates that something has gone awry
    // NOTE: We must also require that the time is greater than -1 because apparently some uncalibrated cells
    //       will report their time as -1. We don't want to include those cells.
    if (tCell < -400e-9 && tCell > -1 && fillingEmbeddedCells) {
      embeddedCellWithLateCellTime = true;
    }
  }

  // If we have one embedded cell with late cell time, then we want to fill out the QA to
  // help identify the event.
  std::string embeddedEventUsed = prefix + "/fHistEmbeddedEventUsed";
  std::string embeddedInteranlEventSelection = prefix + "/fHistInternalEventSelection";
  if (embeddedCellWithLateCellTime)
  {
    auto embeddingInstance = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
    if (embeddingInstance) {
      fHistManager.FillTH1(embeddedEventUsed.c_str(),
                 static_cast<double>(fPreviousEmbeddedEventSelected));

      // Determine the physics selection. This isn't quite a perfect way to store it, as it mingles the
      // selections between different events. But it is simple, which will let us investigate quickly.
      // Plus, it's a reasonable bet that the event selection when be the same when it goes wrong.
      std::bitset<sizeof(UInt_t) * 8> testBits = fPreviousEventTrigger;
      for (unsigned int i = 0; i < 32; i++) {
        if (testBits.test(i)) {
          fHistManager.FillTH2(embeddedInteranlEventSelection.c_str(),
                     static_cast<double>(fPreviousEmbeddedEventSelected), i);
        }
      }
    }
  }

  std::string nCellsName = prefix + "/fHistNCells";
  fHistManager.FillTH1(nCellsName.c_str(), cells->GetNumberOfCells());
}

/**
 * Directs the filling of cell QA into the defined histograms.
 */
void AliAnalysisTaskEmcalJetHPerformance::FillCellQAHists()
{
  // Fill standard cell QA
  std::string prefix = "QA/";
  prefix += fCaloCellsName.Data();
  FillCellQAHists(prefix, fCaloCells);

  // Fill embedded cell QA it if's available.
  auto embeddingInstance = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (embeddingInstance) {
    auto embeddedCells = dynamic_cast<AliVCaloCells*>(
     embeddingInstance->GetExternalEvent()->FindListObject(fEmbeddedCellsName.c_str()));
    if (embeddedCells) {
      prefix = "QA/embedding/";
      prefix += fEmbeddedCellsName;
      FillCellQAHists(prefix, embeddedCells);
    }
  }
}

/**
 * Fill QA histograms.
 */
void AliAnalysisTaskEmcalJetHPerformance::FillQAHists()
{
  if (fCreateCellQAHists) {
    FillCellQAHists();
  }

  // Tracks
  AliTrackContainer * trackCont = 0;
  TIter nextTrackColl(&fParticleCollArray);
  while ((trackCont = static_cast<AliTrackContainer*>(nextTrackColl()))) {
    for (auto track : trackCont->accepted())
    {
      fHistManager.FillTH3(TString::Format("QA/%s/fHistTrackPtEtaPhi", trackCont->GetName()), track->Pt(),
                 track->Eta(), track->Phi());
      fHistManager.FillTH3(TString::Format("QA/%s/fHistTrackPtEtaPhiEfficiencyCorrected", trackCont->GetName()),
                 track->Pt(), track->Eta(), track->Phi(),
                 DetermineTrackingEfficiency(track->Pt(), track->Eta()));
    }
  }

  // Clusters
  AliClusterContainer* clusCont = 0;
  TIter nextClusColl(&fClusterCollArray);
  AliVCluster * cluster = 0;
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusColl()))) {
    for (auto clusIter : clusCont->accepted_momentum())
    {
      AliTLorentzVector c = clusIter.first;
      cluster = clusIter.second;
      // Intentionally plotting against raw energy
      fHistManager.FillTH2(TString::Format("QA/%s/fHistClusterEnergyVsTime", clusCont->GetName()), cluster->E(), cluster->GetTOF());
      // Hadronically corrected energy (which is usually what we're using)
      fHistManager.FillTH1(TString::Format("QA/%s/fHistClusterHadCorrEnergy", clusCont->GetName()),
                 cluster->GetHadCorrEnergy());
      fHistManager.FillTH2(TString::Format("QA/%s/fHistClusterEtaPhi", clusCont->GetName()), c.Eta(),
                 c.Phi_0_2pi());
    }
  }

  // Jets
  AliJetContainer * jetCont = 0;
  double jetPt = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(nextJetColl()))) {
    double rhoVal = jetCont->GetRhoVal();
    for (auto jet : jetCont->accepted())
    {
      jetPt = AliAnalysisTaskEmcalJetHUtils::GetJetPt(jet, rhoVal);
      fHistManager.FillTH1(TString::Format("QA/%s/fHistJetPt", jetCont->GetName()), jetPt);
    }
  }

  // Event plane resolution
  // Only attempt to fill the histogram if the corrected Qn vectors are available.
  if (fFlowQnVectorManager) {
    // Retrieve the Qn vectors
    const AliQnCorrectionsQnVector * vzeroQnVector = fFlowQnVectorManager->GetDetectorQnVector("VZEROQoverM");
    const AliQnCorrectionsQnVector * tpcPosQnVector = fFlowQnVectorManager->GetDetectorQnVector("TPCPosEtaQoverM");
    const AliQnCorrectionsQnVector * tpcNegQnVector = fFlowQnVectorManager->GetDetectorQnVector("TPCNegEtaQoverM");
    if (vzeroQnVector == nullptr || tpcPosQnVector == nullptr || tpcNegQnVector == nullptr) {
      AliWarningStream() << "Q vector unavailable. Skipping.\n";
    }
    else {
      // We need the centrality bin, but the bin provided by the base class isn't that fine grained
      // (and may cause problems elsewhere if it is). So we just calculate the value here.
      std::string name = "QA/eventPlaneRes/%s/R%i";
      // The resolution is calculated by changing the leading term, but it is all with respect
      // to the same harmonic. We measure with respect to the second order harmonic, so we retrieve
      // that event plane.
      double vzero = vzeroQnVector->EventPlane(2);
      double tpcPos = tpcPosQnVector->EventPlane(2);
      double tpcNeg = tpcNegQnVector->EventPlane(2);
      for (unsigned int R = 2; R < 9; R++) {
        fHistManager.FillProfile(TString::Format(name.c_str(), "VZERO", R), fCent,
                     std::cos(R * (tpcPos - tpcNeg)));
        fHistManager.FillProfile(TString::Format(name.c_str(), "TPC_Pos", R), fCent,
                     std::cos(R * (vzero - tpcNeg)));
        fHistManager.FillProfile(TString::Format(name.c_str(), "TPC_Neg", R), fCent,
                     std::cos(R * (vzero - tpcPos)));
      }
    }
  }

  // Update the previous event trigger to the current event trigger so that it is available next event.
  // This is stored for keeping track of when embedded events are double corrected.
  // This must be updated after filling the relevant hists above!
  fPreviousEventTrigger =
   ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  auto embeddingInstance = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (embeddingInstance) {
    fPreviousEmbeddedEventSelected = embeddingInstance->EmbeddedEventUsed();
  }
}

/**
 * Reset matching for the jet container.
 *
 * @params[in] c Jet container for the matching to be reset.
 */
void AliAnalysisTaskEmcalJetHPerformance::ResetMatching(const AliJetContainer &c) const
{
  for(auto j : c.all()){
    j->ResetMatching();
  }
}

/**
 * Perform matching between jet collections.
 *
 * This code is basically combined from `AliAnalysisTaskEmcalJetTagger::MatchJetsGeo(...)`` and
 * ``AliEmcalJetTaggerTaskFast::MatchJetsGeo(...)``. I mainly took the appraoch from the fast tagger,
 * but I'm not using the KDTrees because they still need some additional work as of to properly wrap
 * around at 2pi in phi. So for now we just use brute force matching at O(n^2).
 *
 * Why copy this code? It's a less than ideal thing to do, but for some unidentified, but well known reason
 * (present from at least 2018), there is an issue that causes most jobs to fail when embedding into LHC15o
 * for full jets when using one of the tagger tasks. This particularly seems to occur when embedding LHC16j5.
 * This can be worked around by not running the standard or fast taggers. Thus, the code is copied.
 *
 * NOTE: This makes a number of assumptions about the jet containers that are specific to the jet-hadron analysis
 *       for full jets. Among other factors, it expects that the acceptances are restricted to the EMCal fiducial
 *       volume for detector level jets. If these assumptions are not true, matching may not work properly.
 *
 * @params[in] contBase Base jet container.
 * @params[in] contTag Second jet container.
 * @params[in] maxDist Max distance for matching.
 * @returns True if the matching was successful.
 */
bool AliAnalysisTaskEmcalJetHPerformance::PerformGeometricalJetMatching(AliJetContainer& contBase,
                                    AliJetContainer& contTag, double maxDist) const
{
  // Setup
  const Int_t kNacceptedBase = contBase.GetNAcceptedJets(), kNacceptedTag = contTag.GetNAcceptedJets();
  if (!(kNacceptedBase && kNacceptedTag)) {
    return false;
  }

  // Build up vectors of jet pointers to use when assigning the closest jets.
  // The storages are needed later for applying the tagging, in order to avoid multiple occurrence of jet selection
  std::vector<AliEmcalJet*> jetsBase(kNacceptedBase), jetsTag(kNacceptedTag);

  int countBase(0), countTag(0);
  for (auto jb : contBase.accepted()) {
    jetsBase[countBase] = jb;
    countBase++;
  }
  for (auto jt : contTag.accepted()) {
    jetsTag[countTag] = jt;
    countTag++;
  }

  TArrayI faMatchIndexTag(kNacceptedBase), faMatchIndexBase(kNacceptedTag);
  faMatchIndexBase.Reset(-1);
  faMatchIndexTag.Reset(-1);

  // find the closest distance to the base jet
  countBase = 0;
  for (auto jet1 : contBase.accepted()) {
    double distance = maxDist;

    // Loop over all accepted jets and brute force search for the closest jet.
    // NOTE: current_index() returns the jet index in the underlying array, not
    //       the index within the accepted jets that are returned.
    int contTagAcceptedIndex = 0;
    for (auto jet2 : contTag.accepted()) {
      double dR = jet1->DeltaR(jet2);
      if (dR < distance && dR < maxDist) {
        faMatchIndexTag[countBase] = contTagAcceptedIndex;
        distance = dR;
      }
      contTagAcceptedIndex++;
    }

    // Let us know whether a match was found successfully.
    if (faMatchIndexTag[countBase] >= 0 && distance < maxDist) {
      AliDebugStream(1) << "Found closest tag jet for " << countBase << " with match index "
               << faMatchIndexTag[countBase] << " and distance " << distance << "\n";
    } else {
      AliDebugStream(1) << "Not found closest tag jet for " << countBase << ".\n";
    }

    countBase++;
  }

  // other way around
  countTag = 0;
  for (auto jet1 : contTag.accepted()) {
    double distance = maxDist;

    // Loop over all accepted jets and brute force search for the closest jet.
    // NOTE: current_index() returns the jet index in the underlying array, not
    //       the index within the accepted jets that are returned.
    int contBaseAcceptedIndex = 0;
    for (auto jet2 : contBase.accepted()) {
      double dR = jet1->DeltaR(jet2);
      if (dR < distance && dR < maxDist) {
        faMatchIndexBase[countTag] = contBaseAcceptedIndex;
        distance = dR;
      }
      contBaseAcceptedIndex++;
    }

    // Let us know whether a match was found successfully.
    if (faMatchIndexBase[countTag] >= 0 && distance < maxDist) {
      AliDebugStream(1) << "Found closest base jet for " << countTag << " with match index "
               << faMatchIndexBase[countTag] << " and distance " << distance << "\n";
    } else {
      AliDebugStream(1) << "Not found closest base jet for " << countTag << ".\n";
    }

    countTag++;
  }

  // check for "true" correlations
  // these are pairs where the base jet is the closest to the tag jet and vice versa
  // As the lists are linear a loop over the outer base jet is sufficient.
  AliDebugStream(1) << "Starting true jet loop: nbase(" << kNacceptedBase << "), ntag(" << kNacceptedTag << ")\n";
  for (int ibase = 0; ibase < kNacceptedBase; ibase++) {
    AliDebugStream(2) << "base jet " << ibase << ": match index in tag jet container " << faMatchIndexTag[ibase]
             << "\n";
    if (faMatchIndexTag[ibase] > -1) {
      AliDebugStream(2) << "tag jet " << faMatchIndexTag[ibase] << ": matched base jet "
               << faMatchIndexBase[faMatchIndexTag[ibase]] << "\n";
    }
    // We have a true correlation where each jet points to the other.
    if (faMatchIndexTag[ibase] > -1 && faMatchIndexBase[faMatchIndexTag[ibase]] == ibase) {
      AliDebugStream(2) << "found a true match \n";
      AliEmcalJet *jetBase = jetsBase[ibase], *jetTag = jetsTag[faMatchIndexTag[ibase]];
      // We have a valid pair of matched jets, so set the closest jet properties.
      if (jetBase && jetTag) {
        Double_t dR = jetBase->DeltaR(jetTag);
        jetBase->SetClosestJet(jetTag, dR);
        jetTag->SetClosestJet(jetBase, dR);
      }
    }
  }
  return true;
}

/**
 * Fill jet matching QA histograms to describe the matching quality.
 *
 * @param[in] contBase Base jet container.
 * @param[in] contTag Second jet container.
 */
void AliAnalysisTaskEmcalJetHPerformance::FillJetMatchingQA(AliJetContainer& contBase, AliJetContainer& contTag,
                              const std::string& prefix)
{
  // Fill histograms.
  std::string name = "";
  std::string centBin = std::to_string(fCentBin);
  AliDebugStream(2) << "Filling matching hists with prefix: " << prefix << "\n";
  for (auto jet1 : contBase.accepted()) {
    // Setup
    Double_t ptJet1 = jet1->Pt() - contBase.GetRhoVal() * jet1->Area();
    // Record jet 1 only properties
    AliDebugStream(4) << "jet1: " << jet1->toString() << "\n";
    name = "jetMatching/" + prefix + "/fh2PtJet1VsLeadPtAllSel_" + centBin;
    fHistManager.FillTH2(name.c_str(), ptJet1, jet1->MaxTrackPt());

    // Retrieve jet 2
    AliEmcalJet * jet2 = jet1->ClosestJet();
    if (!jet2) { continue; }
    AliDebugStream(4) << "jet2: " << jet2->toString() << "\n";
    Double_t ptJet2 = jet2->Pt() - contTag.GetRhoVal() * jet2->Area();
    // This will retrieve the fraction of jet2's momentum in jet1.
    Double_t fraction = contBase.GetFractionSharedPt(jet1);

    name = "jetMatching/" + prefix + "/fh2PtJet2VsFraction_" + centBin;
    fHistManager.FillTH2(name.c_str(), ptJet2, fraction);
    AliDebugStream(5) << "Fraction: " << fraction << ", minimum: " << fMinFractionShared << "\n";
    if (fraction < fMinFractionShared) {
      continue;
    }
    name = "jetMatching/" + prefix + "/fh2PtJet1VsLeadPtTagged_" + centBin;
    fHistManager.FillTH2(name.c_str(), ptJet1,jet1->MaxTrackPt());
    name = "jetMatching/" + prefix + "/fh2PtJet1VsPtJet2_" + centBin;
    fHistManager.FillTH2(name.c_str(), ptJet1, ptJet2);
    if (ptJet2 > 0.) {
      name = "jetMatching/" + prefix + "/fh2PtJet2VsRelPt_" + centBin;
      fHistManager.FillTH2(name.c_str(), ptJet2, (ptJet1 - ptJet2) / ptJet2);
    }

    // Recall that the arguments are backwards of the expectation.
    Double_t dPhi = DeltaPhi(jet2->Phi(), jet1->Phi());
    if (dPhi > TMath::Pi()) {
      dPhi -= TMath::TwoPi();
    }
    if (dPhi < (-1. * TMath::Pi())) {
      dPhi += TMath::TwoPi();
    }

    name = "jetMatching/" + prefix + "/fh3PtJet1VsDeltaEtaDeltaPhi_" + centBin;
    fHistManager.FillTH3(name.c_str(), ptJet1, jet1->Eta() - jet2->Eta(), dPhi);
    name = "jetMatching/" + prefix + "/fh2PtJet1VsDeltaR_" + centBin;
    fHistManager.FillTH2(name.c_str(), ptJet1, jet1->DeltaR(jet2));
  }

  // Number of accepted jets
  name = "jetMatching/" + prefix + "/fNAccJets";
  fHistManager.FillTH1(name.c_str(), contBase.GetNAcceptedJets());
}

/**
 * Process the jets according to defined cuts and fill the response matrix.
 */
void AliAnalysisTaskEmcalJetHPerformance::ResponseMatrix()
{
  AliJetContainer * jetsHybrid = GetJetContainer("hybridLevelJets");
  AliJetContainer * jetsDetLevel = GetJetContainer("detLevelJets");
  AliJetContainer * jetsPartLevel = GetJetContainer("partLevelJets");
  // NOTE: Defaults to 0 if rho was not specified.
  double rhoHybrid = jetsHybrid->GetRhoVal();
  if (!jetsHybrid) {
    AliErrorStream() << "Could not retrieve hybrid jet collection.\n";
    return;
  }
  if (!jetsDetLevel) {
    AliErrorStream() << "Could not retrieve det level jet collection.\n";
    return;
  }
  if (fResponseFromThreeJetCollections && !jetsPartLevel) {
    AliErrorStream() << "Could not retrieve part level jet collection.\n";
    return;
  }

  // Handle matching of jets.
  for (auto jet1 : jetsHybrid->accepted())
  {
    AliDebugStream(4) << "jet1: " << jet1->toString() << "\n";
    AliDebugStream(4) << "jet1 address: " << jet1 << "\n";

    // Get jet the det level jet from the hybrid jet
    AliEmcalJet * jet2 = jet1->ClosestJet();
    if(!jet2) continue;

    AliDebugStream(4) << "jet2: " << jet2->toString() << "\n";
    AliDebugStream(4) << "jet2 address: " << jet2 << "\n";

    // Check shared fraction
    double sharedFraction = jetsHybrid->GetFractionSharedPt(jet1);
    fHistManager.FillTH1("fHistFractionSharedPt", sharedFraction);
    if (sharedFraction < fMinFractionShared) {
      AliDebugStream(4) << "Rejecting jet due to momentum fraction of " << sharedFraction << ", smaller than the minimum.\n";
      continue;
    }
    else {
      AliDebugStream(4) << "Jet passed momentum fraction cut with value of " << sharedFraction << "\n";
    }

    // NOTE: We apply no explicit event selection to jet 2 because the matching in the tagger
    // only matches jets which are accepted.

    // Determine the jet to use to fill the response
    // The jet that is passed may be the embedded detector level (for two stage matching), or it may
    // be the embedded particle level (for direct matching).
    AliEmcalJet * jetToPass = 0;
    if (fResponseFromThreeJetCollections) {
      // Retrieve the MC level jet
      AliEmcalJet * jet3 = jet2->ClosestJet();
      UInt_t rejectionReason = 0;
      if (!jetsPartLevel->AcceptJet(jet3, rejectionReason)) {
        // NOTE: This shouldn't ever happen because the tagger applies acceptance
        //       cuts when matching. However, we keep the check here for good measure.
        // NOTE: Could store rejection reasons if needed with below:
        //fHistRejectionReason2->Fill(jets2->GetRejectionReasonBitPosition(rejectionReason), jet2->Pt());
        continue;
      }

      AliDebugStream(4) << "jet3: " << jet3->toString() << "\n";
      AliDebugStream(4) << "jet3 address: " << jet3 << "\n";

      // Use for the response
      AliDebugStream(4) << "Using part level jet for response (ie. two stage matching)\n";
      jetToPass = jet3;
    }
    else {
      // Use for the response
      AliDebugStream(4) << "Using one stage matching for response\n";
      jetToPass = jet2;
    }

    // Fill response
    FillResponseMatrix(jet1, jetToPass, rhoHybrid);
  }

}

/**
 * If given multiple jet collections, handle creating a reasponse matrix
 */
void AliAnalysisTaskEmcalJetHPerformance::FillResponseMatrix(AliEmcalJet * jet1, AliEmcalJet * jet2, const double jet1Rho)
{
  if (!jet1 || !jet2) {
    AliErrorStream() << "Null jet passed to fill response matrix";
  }

  AliDebugStream(3) << "About to create ResponseMatrixFillWrappers\n";
  AliDebugStream(4) << "jet1: " << jet1->toString() << "\n";
  AliDebugStream(4) << "jet2: " << jet2->toString() << "\n";
  // Create map from jetNumber to jet and initialize the objects
  std::map<unsigned int, ResponseMatrixFillWrapper> jetNumberToJet = {
    std::make_pair(1, CreateResponseMatrixFillWrapper(jet1, jet1Rho)),
    // We would never want to rho subtract the second jet, regardless of whether it's external
    // detector level or particle level.
    std::make_pair(2, CreateResponseMatrixFillWrapper(jet2, 0))
  };

  // Fill histograms
  std::string histName = "response/fHistResponseMatrix";
  std::vector<double> values;
  THnSparse * response = static_cast<THnSparse*>(fHistManager.FindObject(histName.c_str()));
  AliDebugStream(3) << "About to fill response matrix values\n";
  AliDebugStream(4) << "jet1: " << jet1->toString() << "\n";
  AliDebugStream(4) << "jet2: " << jet2->toString() << "\n";
  for (unsigned int i = 0; i < response->GetNdimensions(); i++) {
    std::string title = response->GetAxis(i)->GetTitle();

    // Retrieve pair of jet and pointer to extract the fill value
    auto jetPair = fResponseMatrixFillMap.find(title);
    if (jetPair != fResponseMatrixFillMap.end()) {
      auto wrapper = jetNumberToJet.at(jetPair->second.first);
      auto member = jetPair->second.second;
      AliDebugStream(4) << "Filling value " << wrapper.*member << " into axis " << title << "\n";
      values.emplace_back(wrapper.*member);
    }
    else {
      AliWarningStream() << "Unable to fill dimension " << title << "!\n";
    }
  }

  fHistManager.FillTHnSparse(histName.c_str(), values.data());
}

/**
 *
 */
AliAnalysisTaskEmcalJetHPerformance::ResponseMatrixFillWrapper AliAnalysisTaskEmcalJetHPerformance::CreateResponseMatrixFillWrapper(AliEmcalJet * jet, const double rho) const
{
  ResponseMatrixFillWrapper wrapper;
  if (!jet) {
    AliErrorStream() << "Must pass valid jet to create object.\n";
    return wrapper;
  }
  wrapper.fPt = AliAnalysisTaskEmcalJetHUtils::GetJetPt(jet, rho);
  wrapper.fArea = jet->Area();
  wrapper.fPhi = jet->Phi();
  wrapper.fDistance = jet->ClosestJetDistance();
  wrapper.fCentrality = fCent;
  wrapper.fRelativeEPAngle = AliAnalysisTaskEmcalJetHUtils::RelativeEPAngle(jet->Phi(), fEPV0);
  wrapper.fLeadingHadronPt = AliAnalysisTaskEmcalJetHUtils::GetLeadingHadronPt(jet, fLeadingHadronBiasType);

  return wrapper;
}

/**
 * Add task to an existing analysis manager.
 *
 * @param[in] suffix Suffix string to attach to the task name
 */
AliAnalysisTaskEmcalJetHPerformance * AliAnalysisTaskEmcalJetHPerformance::AddTaskEmcalJetHPerformance(const char * suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    AliErrorClass("No analysis manager to connect to.");
    return nullptr;
  }

  // Setup task name
  std::string taskName = "AliAnalysisTaskEmcalJetHPerformance";
  std::string suffixName(suffix);
  if (suffixName != "") {
    taskName += "_";
    taskName += suffixName;
  }

  // Create task and configure as desired.
  AliAnalysisTaskEmcalJetHPerformance * task = new AliAnalysisTaskEmcalJetHPerformance(taskName.c_str());
  // Set a few general default.
  task->SetNCentBins(5);
  // Configuration is via YAML.
  mgr->AddTask(task);

  // Create containers for input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer* outputContainer =
   mgr->CreateContainer(task->GetName(), TList::Class(), AliAnalysisManager::kOutputContainer,
              Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectOutput(task, 1, outputContainer);

  return task;
}

/**
 * Prints information about the jet-hadron performance task.
 *
 * @return std::string containing information about the task.
 */
std::string AliAnalysisTaskEmcalJetHPerformance::toString() const
{
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  tempSS << "Particle collections:\n";
  TIter nextParticleCont(&fParticleCollArray);
  AliParticleContainer * particleCont;
  while ((particleCont = static_cast<AliParticleContainer *>(nextParticleCont()))) {
    tempSS << "\t" << particleCont->GetName() << ": " << particleCont->GetArrayName() << "\n";
  }
  tempSS << "Cluster collections:\n";
  TIter nextClusterCont(&fClusterCollArray);
  AliClusterContainer * clusterCont;
  while ((clusterCont = static_cast<AliClusterContainer *>(nextClusterCont()))) {
    tempSS << "\t" << clusterCont->GetName() << ": " << clusterCont->GetArrayName() << "\n";
  }
  tempSS << "Jet collections:\n";
  TIter nextJetCont(&fJetCollArray);
  AliJetContainer * jetCont;
  while ((jetCont = static_cast<AliJetContainer *>(nextJetCont()))) {
    tempSS << "\t" << jetCont->GetName() << ": " << jetCont->GetArrayName() << "\n";
  }
  // AliEventCuts
  tempSS << "AliEventCuts\n";
  tempSS << "\tEnabled: " << !fUseBuiltinEventSelection << "\n";
  // Efficiency
  tempSS << "Efficiency\n";
  tempSS << "\tSingle track efficiency identifier: " << fEfficiencyPeriodIdentifier << "\n";
  // QA
  tempSS << "QA Hists:\n";
  tempSS << "\tEnabled: " << fCreateQAHists << "\n";
  tempSS << "\tCreate cell QA hists: " << fCreateCellQAHists << "\n";
  // Use whether the pointer as null as a proxy. It's not ideal because it's not fully initialized
  // until after UserCreateOutputObjects(). But it's good enough for these purposes.
  tempSS << "\tCalculate event plane resolution (proxy of whether it's enabled - it may not be accurate): " << (fFlowQnVectorManager != nullptr) << "\n";
  // Jet matching
  tempSS << "Jet matching:\n";
  tempSS << "\tEnabled: " << fPerformJetMatching << "\n";
  tempSS << "\tMax matching distance: " << fMaxJetMatchingDistance << "\n";
  // Response matrix
  tempSS << "Response matrix:\n";
  tempSS << "\tEnabled: " << fCreateResponseMatrix << "\n";
  tempSS << "\tConstruct response from 3 jet collections: " << fResponseFromThreeJetCollections << "\n";
  tempSS << "\tMin fraction shared pt: " << fMinFractionShared << "\n";
  tempSS << "\tJet leading hadron bias type: " << fLeadingHadronBiasType << "\n";
  tempSS << "\tResponse matrix fill map: \n";
  for (auto el : fResponseMatrixFillMap) {
    tempSS << "\t\tProperty " << el.first << " applied to jet " << el.second.first << "\n";
  }

  return tempSS.str();
}

/**
 * Print jet-hadron performance task information on an output stream using the string representation provided by
 * AliAnalysisTaskEmcalJetHPerformance::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream & AliAnalysisTaskEmcalJetHPerformance::Print(std::ostream & in) const {
  in << toString();
  return in;
}

/**
 * Print basic jet-hadron performance task information using the string representation provided by
 * AliAnalysisTaskEmcalJetHPerformance::toString
 *
 * @param[in] opt Unused
 */
void AliAnalysisTaskEmcalJetHPerformance::Print(Option_t* opt) const
{
  Printf("%s", toString().c_str());
}

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

/**
 * Implementation of the output stream operator for AliAnalysisTaskEmcalJetHPerformance. Printing
 * basic jet-hadron performance task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream & operator<<(std::ostream & in, const PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance & myTask)
{
  std::ostream & result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 *
 * NOTE: We don't swap the base class values because the base class doesn't implement swap.
 */
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance & first, PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance & second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fYAMLConfig, second.fYAMLConfig);
  swap(first.fConfigurationInitialized, second.fConfigurationInitialized);
  //swap(first.fHistManager, second.fHistManager); // Skip here, because the THistManager copy constructor is private.
  //swap(first.fEmbeddingQA, second.fEmbeddingQA); // Skip here, because the THistManager copy constructor is private.
  swap(first.fCreateQAHists, second.fCreateQAHists);
  swap(first.fCreateCellQAHists, second.fCreateCellQAHists);
  swap(first.fPerformJetMatching, second.fPerformJetMatching);
  swap(first.fCreateResponseMatrix, second.fCreateResponseMatrix);
  swap(first.fEmbeddedCellsName, second.fEmbeddedCellsName);
  swap(first.fPreviousEventTrigger, second.fPreviousEventTrigger);
  swap(first.fPreviousEmbeddedEventSelected, second.fPreviousEmbeddedEventSelected);
  swap(first.fEfficiencyPeriodIdentifier, second.fEfficiencyPeriodIdentifier);
  swap(first.fFlowQnVectorManager, second.fFlowQnVectorManager);
  swap(first.fResponseMatrixFillMap, second.fResponseMatrixFillMap);
  swap(first.fResponseFromThreeJetCollections, second.fResponseFromThreeJetCollections);
  swap(first.fMinFractionShared, second.fMinFractionShared);
  swap(first.fLeadingHadronBiasType, second.fLeadingHadronBiasType);
}
