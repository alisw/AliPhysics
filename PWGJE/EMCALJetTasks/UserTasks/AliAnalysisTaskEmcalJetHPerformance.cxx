/**
 * Implemenation for AliAnalysisTaskEmcalJetHPerformance
 */

#include "AliAnalysisTaskEmcalJetHPerformance.h"
#include "AliAnalysisTaskEmcalJetHUtils.h"

#include <map>
#include <vector>
#include <iostream>

#include <TObject.h>
#include <TCollection.h>
#include <TAxis.h>
#include <THnSparse.h>

#include <AliAnalysisManager.h>
#include <AliLog.h>

#include "yaml-cpp/yaml.h"
#include "THistManager.h"
#include "AliJetContainer.h"
#include "AliEmcalJet.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance);

namespace PWGJE {
namespace EMCALJetTasks {

/**
 * Default constructor
 */
AliAnalysisTaskEmcalJetHPerformance::AliAnalysisTaskEmcalJetHPerformance():
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetHPerformance", kFALSE),
  fYAMLConfig(),
  fConfigurationPath(""),
  fConfigurationInitialized(false),
  fHistManager(),
  fEmbeddingQA(),
  fCreateResponseMatrix(false),
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
  fConfigurationPath(""),
  fConfigurationInitialized(false),
  fHistManager(name),
  fCreateResponseMatrix(false),
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
  fConfigurationPath(other.fConfigurationPath),
  fConfigurationInitialized(other.fConfigurationInitialized),
  fHistManager(other.fHistManager.GetName()),
  fCreateResponseMatrix(other.fCreateResponseMatrix),
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
 * Retrieve task properties from the YAML configuration
 */
void AliAnalysisTaskEmcalJetHPerformance::RetrieveAndSetTaskPropertiesFromYAMLConfig()
{
  // Same ordering as in the constructor (for consistency)
  std::string baseName = "enable";
  fYAMLConfig.GetProperty({baseName, "responseMatrix"}, fCreateResponseMatrix, false);

  // Response matrix properties
  baseName = "responseMatrix";
  fYAMLConfig.GetProperty({baseName, "useThreeJetCollections"}, fResponseFromThreeJetCollections, false);
  fYAMLConfig.GetProperty({baseName, "minFractionSharedPt"}, fMinFractionShared, false);
  std::string hadronBiasStr = "";
  baseName = "jets";
  bool res = fYAMLConfig.GetProperty({baseName, "leadingHadronBiasType"}, hadronBiasStr, false);
  // Only attempt to set the property if it is retrieved successfully
  if (res) {
    fLeadingHadronBiasType = AliAnalysisTaskEmcalJetHUtils::fgkLeadingHadronBiasMap.at(hadronBiasStr);
  }

  // Base class options
  // Recycle unused embedded events
  fYAMLConfig.GetProperty("recycleUnusedEmbeddedEventsMode", fRecycleUnusedEmbeddedEventsMode, false);
  // Centrality
  std::pair<double, double> centRange;
  res = fYAMLConfig.GetProperty("centralityRange", centRange, false);
  if (res) {
    SetCentRange(centRange.first, centRange.second);
  }
}

/**
 * Create jet containers based on the jet values defined in the YAML config
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupJetContainersFromYAMLConfig()
{
  std::string baseName = "jets";
  std::vector <std::string> jetNames = {"hybridLevelJets", "detLevelJets", "partLevelJets"};
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

  // Setup and initialize the YAML config
  if (fConfigurationPath != "") {
    AliInfoStream() << "Adding YAML configuration found at \"" << fConfigurationPath << "\"\n";
    int configPosition = fYAMLConfig.AddConfiguration(fConfigurationPath, "yamlConfig");
    if (configPosition < 0) {
      // Return immediately
      return fConfigurationInitialized;
    }
  }
  else {
    AliInfoStream() << "No YAML configuration fileanme passed.\n";
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
  // Check that the task was initialized
  if (!fConfigurationInitialized) {
    AliFatal("Task was not initialized. Please ensure that Initialize() was called!");
  }

  // Reinitialize the YAML configuration
  fYAMLConfig.Reinitialize();

  // Call the base class
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Create histograms
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
 * Setup and allocate histograms related to creating a response matrix.
 */
void AliAnalysisTaskEmcalJetHPerformance::SetupResponseMatrixHists()
{
  // Main response matrix THnSparse
  std::string name = "fHistResponseMatrix";
  std::string title = "fHistResponseMatrix";
  //histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCal", jets->GetArrayName().Data());

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

  if (fCreateResponseMatrix) {
    CreateResponseMatrix();
  }

  return kTRUE;
}

/**
 * Handle the functionality necessary to create the response matrix.
 */
void AliAnalysisTaskEmcalJetHPerformance::CreateResponseMatrix()
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
    // Get jet the det level jet from the hybrid jet
    AliEmcalJet * jet2 = jet1->ClosestJet();
    if(!jet2) continue;

    AliDebugStream(4) << "jet2: " << jet2->toString() << "\n";

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
  std::string histName = "fHistResponseMatrix";
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

  // Create task and configure as desired
  AliAnalysisTaskEmcalJetHPerformance * task = new AliAnalysisTaskEmcalJetHPerformance(taskName.c_str());
  // Configure via the YAML configuration
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
  tempSS << "Jet collections:\n";
  TIter next(&fJetCollArray);
  AliJetContainer * jetCont;
  while ((jetCont = static_cast<AliJetContainer *>(next()))) {
    tempSS << "\t" << jetCont->GetName() << ": " << jetCont->GetArrayName() << "\n";
  }
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
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550
 */
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance & first, PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance & second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fYAMLConfig, second.fYAMLConfig);
  swap(first.fConfigurationPath, second.fConfigurationPath);
  swap(first.fConfigurationInitialized, second.fConfigurationInitialized);
  //swap(first.fHistManager, second.fHistManager); // Skip here, as there is no copy constructor for THistManager
  swap(first.fCreateResponseMatrix, second.fCreateResponseMatrix);
  swap(first.fResponseMatrixFillMap, second.fResponseMatrixFillMap);
  swap(first.fResponseFromThreeJetCollections, second.fResponseFromThreeJetCollections);
  swap(first.fMinFractionShared, second.fMinFractionShared);
  swap(first.fLeadingHadronBiasType, second.fLeadingHadronBiasType);
}

