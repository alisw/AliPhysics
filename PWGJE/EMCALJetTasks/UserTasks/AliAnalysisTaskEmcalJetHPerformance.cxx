/**
 * Implemenation for AliAnalysisTaskEmcalJetHPerformance
 */

#include "AliAnalysisTaskEmcalJetHPerformance.h"

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
#include "AliClusterContainer.h"
#include "AliEmcalContainerUtils.h"
#include "AliJetContainer.h"
#include "AliEmcalJet.h"

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
  fCreateResponseMatrix(false),
  fEmbeddedCellsName("emcalCells"),
  fPreviousEventTrigger(0),
  fPreviousEmbeddedEventSelected(false),
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
  fCreateResponseMatrix(false),
  fEmbeddedCellsName("emcalCells"),
  fPreviousEventTrigger(0),
  fPreviousEmbeddedEventSelected(false),
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
  fCreateResponseMatrix(other.fCreateResponseMatrix),
  fEmbeddedCellsName(other.fEmbeddedCellsName),
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
 * Retrieve task properties from the YAML configuration.
 */
void AliAnalysisTaskEmcalJetHPerformance::RetrieveAndSetTaskPropertiesFromYAMLConfig()
{
  // Base class options
  // Recycle unused embedded events
  fYAMLConfig.GetProperty("recycleUnusedEmbeddedEventsMode", fRecycleUnusedEmbeddedEventsMode, false);
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
      std::string collectionName = "", acceptance = "";
      double R = -1;
      fYAMLConfig.GetProperty({baseName, jetName, "collection"}, collectionName, true);
      fYAMLConfig.GetProperty({baseName, jetName, "acceptance"}, acceptance, true);
      fYAMLConfig.GetProperty({baseName, jetName, "R"}, R, true);

      // Create jet container and set the name
      AliDebugStream(1) << "Creating jet from jet collection name " << collectionName << " with acceptance " << acceptance << " and R=" << R << "\n";
      AliJetContainer * jetCont = AddJetContainer(collectionName.c_str(), acceptance.c_str(), R);
      jetCont->SetName(jetName.c_str());

      // Leading hadron type
      int leadingHadronType = -1;
      bool res = fYAMLConfig.GetProperty({baseName, jetName, "leadingHadronType"}, leadingHadronType, false);
      if (res) {
        AliDebugStream(1) << "Setting leading hadron type of " << leadingHadronType << " for jet cont " << jetName << "\n";
        jetCont->SetLeadingHadronType(leadingHadronType);
      }
    }
    else {
      AliInfoStream() << "Unable to find definition of jet container corresponding to \"" << jetName << "\"\n";
    }
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

/**
 * Setup and allocate histograms related to QA histograms.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupQAHists()
{
  // Cell level QA
  SetupCellQAHists();

  // Clusters
  AliEmcalContainer* cont = 0;
  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliClusterContainer*>(nextClusColl()))) {
    // Cluster time vs energy
    std::string name = "QA/%s/fHistClusterEnergyVsTime";
    std::string title = name + ";E_{cluster} (GeV);t_{cluster} (s)";
    fHistManager.CreateTH2(TString::Format(name.c_str(), cont->GetName()), TString::Format(title.c_str(), cont->GetName()), 1000, 0, 100, 300, -300e-9, 300e-9);
  }

  // Jets
  TIter nextJetColl(&fJetCollArray);
  while ((cont = static_cast<AliJetContainer*>(nextJetColl()))) {
    // Jet pT
    std::string name = "QA/%s/fHistJetPt";
    std::string title = name + ";p_{T} (GeV)";
    fHistManager.CreateTH1(TString::Format(name.c_str(), cont->GetName()), TString::Format(title.c_str(), cont->GetName()), 500, 0, 250);
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

  if (fCreateQAHists) {
    QAHists();
  }
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
  // No additional processing is necessary
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
  FillCellQAHists();

  // Clusters
  AliClusterContainer* clusCont = 0;
  TIter nextClusColl(&fClusterCollArray);
  AliVCluster * cluster = 0;
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusColl()))) {
    for (auto clusIter : clusCont->accepted_momentum())
    {
      cluster = clusIter.second;
      // Intentionally plotting against raw energy
      fHistManager.FillTH2(TString::Format("QA/%s/fHistClusterEnergyVsTime", clusCont->GetName()), cluster->E(), cluster->GetTOF());
    }
  }

  // Jets
  AliJetContainer * jetCont = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(nextJetColl()))) {
    for (auto jet : jetCont->accepted())
    {
      fHistManager.FillTH1(TString::Format("QA/%s/fHistJetPt", jetCont->GetName()), jet->Pt());
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
 * Process the jets according to defined cuts and fill the response matrix.
 */
void AliAnalysisTaskEmcalJetHPerformance::ResponseMatrix()
{
  AliJetContainer * jetsHybrid = GetJetContainer("hybridLevelJets");
  AliJetContainer * jetsDetLevel = GetJetContainer("detLevelJets");
  AliJetContainer * jetsPartLevel = GetJetContainer("partLevelJets");
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

    // Apply additional selection to jet 2
    // TODO: Should we apply acceptance criteria to jet 2 here?

    // Get MC level jet
    AliEmcalJet * jetToPass = 0;
    if (fResponseFromThreeJetCollections) {
      AliEmcalJet * jet3 = jet2->ClosestJet();

      // Accept jet 3
      UInt_t rejectionReason = 0;
      if (!jetsPartLevel->AcceptJet(jet3, rejectionReason)) {
        // TODO: Store rejection reasons
        //fHistRejectionReason2->Fill(jets2->GetRejectionReasonBitPosition(rejectionReason), jet2->Pt());
        continue;
      }

      AliDebugStream(4) << "jet3: " << jet3->toString() << "\n";
      AliDebugStream(4) << "jet3 address: " << jet3 << "\n";

      // Use for the response
      AliDebugStream(4) << "Using part level jet for response\n";
      jetToPass = jet3;
    }
    else {
      // Use for the response
      AliDebugStream(4) << "Using det level jet for response\n";
      jetToPass = jet2;
    }

    // Fill response
    FillResponseMatrix(jet1, jetToPass);
  }

}

/**
 * If given multiple jet collections, handle creating a reasponse matrix
 */
void AliAnalysisTaskEmcalJetHPerformance::FillResponseMatrix(AliEmcalJet * jet1, AliEmcalJet * jet2)
{
  if (!jet1 || !jet2) {
    AliErrorStream() << "Null jet passed to fill response matrix";
  }

  AliDebugStream(3) << "About to create ResponseMatrixFillWrappers\n";
  AliDebugStream(4) << "jet1: " << jet1->toString() << "\n";
  AliDebugStream(4) << "jet2: " << jet2->toString() << "\n";
  // Create map from jetNumber to jet and initialize the objects
  std::map<unsigned int, ResponseMatrixFillWrapper> jetNumberToJet = {
    std::make_pair(1, CreateResponseMatrixFillWrapper(jet1)),
    std::make_pair(2, CreateResponseMatrixFillWrapper(jet2))
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
AliAnalysisTaskEmcalJetHPerformance::ResponseMatrixFillWrapper AliAnalysisTaskEmcalJetHPerformance::CreateResponseMatrixFillWrapper(AliEmcalJet * jet) const
{
  ResponseMatrixFillWrapper wrapper;
  if (!jet) {
    AliErrorStream() << "Must pass valid jet to create object.\n";
    return wrapper;
  }
  wrapper.fPt = jet->Pt();
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
  AliAnalysisDataContainer * outputContainer = mgr->CreateContainer(task->GetName(),
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
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
  tempSS << "Recycle unused embedded events: " << fRecycleUnusedEmbeddedEventsMode << "\n";
  tempSS << "Jet collections:\n";
  TIter next(&fJetCollArray);
  AliJetContainer * jetCont;
  while ((jetCont = static_cast<AliJetContainer *>(next()))) {
    tempSS << "\t" << jetCont->GetName() << ": " << jetCont->GetArrayName() << "\n";
  }
  tempSS << "AliEventCuts\n";
  tempSS << "\tEnabled: " << !fUseBuiltinEventSelection << "\n";
  tempSS << "QA Hists:\n";
  tempSS << "\tEnabled: " << fCreateQAHists << "\n";
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
  swap(first.fCreateResponseMatrix, second.fCreateResponseMatrix);
  swap(first.fEmbeddedCellsName, second.fEmbeddedCellsName);
  swap(first.fResponseMatrixFillMap, second.fResponseMatrixFillMap);
  swap(first.fResponseFromThreeJetCollections, second.fResponseFromThreeJetCollections);
  swap(first.fMinFractionShared, second.fMinFractionShared);
  swap(first.fLeadingHadronBiasType, second.fLeadingHadronBiasType);
}
