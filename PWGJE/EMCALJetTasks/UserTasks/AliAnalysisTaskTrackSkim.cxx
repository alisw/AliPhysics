/**
 * Track skimming task for jet analysis
 */

#include "AliAnalysisTaskTrackSkim.h"

#include <bitset>
#include <cmath>

#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>

#include <AliAODEvent.h>
#include <AliAODMCHeader.h>
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisManager.h>
#include <AliGenPythiaEventHeader.h>
#include <AliLog.h>
#include <AliMCEvent.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliVTrack.h>

#include "AliParticleContainer.h"

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim);
/// \endcond

namespace PWGJE
{
namespace EMCALJetTasks
{

const std::map<std::string, AliMCAnalysisUtils::generator> AliAnalysisTaskTrackSkim::fgkMCUtilsGeneratorMap = {
  { "kPythia", AliMCAnalysisUtils::generator::kPythia },
  { "kHerwig",  AliMCAnalysisUtils::generator::kHerwig }
};

/**
 * Track skimming task for jet analysis
 */

/**
 * Default constructor.
 */
AliAnalysisTaskTrackSkim::AliAnalysisTaskTrackSkim()
 : AliAnalysisTaskEmcalJet("AliAnalysisTaskTrackSkim", kTRUE),
  fYAMLConfig(),
  fConfigurationInitialized{false},
  fChargeScaling{false},
  fEncodeAdditionalParticleInformation{false},
  fMCAnalysisUtilsGenerator{AliMCAnalysisUtils::generator::kPythia},
  fMCAnalysisUtils{nullptr},
  fNAcceptedFirstTrackCollection{nullptr},
  fTriggerBitINT7{false},
  fCentralityForTree{-1},
  fEventPlaneV0MForTree{-1},
  fTriggerBitCentral{false},
  fTriggerBitSemiCentral{false},
  fFirstTrackCollectionKinematics{},
  fFirstTrackCollectionLabels{},
  fSecondTrackCollectionKinematics{},
  fSecondTrackCollectionLabels{},
  fTreeSkim{nullptr}
{
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/**
 * Standard constructor.
 */
AliAnalysisTaskTrackSkim::AliAnalysisTaskTrackSkim(const char* name)
 : AliAnalysisTaskEmcalJet(name, kTRUE),
  fYAMLConfig(),
  fConfigurationInitialized{false},
  fChargeScaling{false},
  fEncodeAdditionalParticleInformation{false},
  fMCAnalysisUtilsGenerator{AliMCAnalysisUtils::generator::kPythia},
  fMCAnalysisUtils{nullptr},
  fNAcceptedFirstTrackCollection{nullptr},
  fTriggerBitINT7{false},
  fCentralityForTree{-1},
  fEventPlaneV0MForTree{-1},
  fTriggerBitCentral{false},
  fTriggerBitSemiCentral{false},
  fFirstTrackCollectionKinematics{},
  fFirstTrackCollectionLabels{},
  fSecondTrackCollectionKinematics{},
  fSecondTrackCollectionLabels{},
  fTreeSkim{nullptr}
{
  // Standard constructor.
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/**
 * Copy constructor.
 *
 * Copying this task will not copy hists or trees.
 * It is up to the user to create this objects use `UserCreateOutputObjects(...)`.
 */
AliAnalysisTaskTrackSkim::AliAnalysisTaskTrackSkim(
 const AliAnalysisTaskTrackSkim& other)
 : fYAMLConfig(other.fYAMLConfig),
  fConfigurationInitialized(other.fConfigurationInitialized),
  fChargeScaling{other.fChargeScaling},
  fEncodeAdditionalParticleInformation{fEncodeAdditionalParticleInformation},
  fMCAnalysisUtilsGenerator{other.fMCAnalysisUtilsGenerator},
  fMCAnalysisUtils{nullptr},
  fNAcceptedFirstTrackCollection{nullptr},
  fTriggerBitINT7{other.fTriggerBitINT7},
  fCentralityForTree{other.fCentralityForTree},
  fEventPlaneV0MForTree{other.fEventPlaneV0MForTree},
  fTriggerBitCentral{other.fTriggerBitCentral},
  fTriggerBitSemiCentral{other.fTriggerBitSemiCentral},
  fFirstTrackCollectionKinematics{other.fFirstTrackCollectionKinematics},
  fFirstTrackCollectionLabels{other.fFirstTrackCollectionLabels},
  fSecondTrackCollectionKinematics{other.fSecondTrackCollectionKinematics},
  fSecondTrackCollectionLabels{other.fSecondTrackCollectionLabels},
  fTreeSkim{nullptr}
{
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
AliAnalysisTaskTrackSkim& AliAnalysisTaskTrackSkim::operator=(
 AliAnalysisTaskTrackSkim other)
{
  swap(*this, other);
  return *this;
}

AliMCAnalysisUtils * AliAnalysisTaskTrackSkim::GetMCAnalysisUtils()
{
  if (!fMCAnalysisUtils){
    fMCAnalysisUtils = new AliMCAnalysisUtils();
    fMCAnalysisUtils->SetMCGenerator(fMCAnalysisUtilsGenerator);
  }
  return fMCAnalysisUtils;
}

/**
 * Retrieve task properties from the YAML configuration.
 */
void AliAnalysisTaskTrackSkim::RetrieveAndSetTaskPropertiesFromYAMLConfig()
{
  // General options
  std::string baseName = "general";
  // Centrality
  std::pair<double, double> centRange;
  bool res = fYAMLConfig.GetProperty({ baseName, "centralityRange" }, centRange, false);
  if (res) {
    AliDebugStream(3) << "Setting centrality range of (" << centRange.first << ", " << centRange.second << ").\n";
    fMinCent = centRange.first;
    fMaxCent = centRange.second;
  }
  // Scale pt by the charge to encode the charge
  res = fYAMLConfig.GetProperty({ baseName, "chargeScaling" }, fChargeScaling, false);
  if (res) {
    AliDebugStream(3) << "Charge scaling enabled.\n";
  }
  // Settings for additional encoded particle information
  res = fYAMLConfig.GetProperty({baseName, "encodeAdditionalParticleInformation"}, fEncodeAdditionalParticleInformation, false);
  if (res) {
    AliDebugStream(3) << "Will encode additional particle level information.\n";
  }
  // MC utils generator
  std::string tempStr = "";
  res = fYAMLConfig.GetProperty({baseName, "MCAnalysisUtilsGenerator"}, tempStr, false);
  if (res) {
    AliDebugStream(3) << "MC analysis utils generator set to '" << tempStr << "'.\n";
    fMCAnalysisUtilsGenerator = fgkMCUtilsGeneratorMap.at(tempStr);
  }
}

/**
 * Initialize task.
 */
bool AliAnalysisTaskTrackSkim::Initialize()
{
  fConfigurationInitialized = false;
  // Always initialize for streaming purposes
  fYAMLConfig.Initialize();

  // Setup task based on the properties defined in the YAML config
  AliDebugStream(2) << "Configuring task from the YAML configuration.\n";
  RetrieveAndSetTaskPropertiesFromYAMLConfig();
  AliDebugStream(2) << "Finished configuring via the YAML configuration.\n";

  // Print the results of the initialization
  // Print outside of the ALICE Log system to ensure that it is always available!
  std::cout << *this;

  fConfigurationInitialized = true;
  return fConfigurationInitialized;
}
/**
  * Define particle kinematics branches.
  */
void AliAnalysisTaskTrackSkim::AddParticleKinematicsToMap(std::map<std::string, std::vector<float>> & map)
{
  map["pt"] = {};
  map["eta"] = {};
  map["phi"] = {};
}
/**
  * Define particle label branches.
  */
void AliAnalysisTaskTrackSkim::AddParticleLabelsToMap(std::map<std::string, std::vector<int>> & map)
{
  map["particle_ID"] = {};
  map["label"] = {};
  if (fEncodeAdditionalParticleInformation) {
    map["encoded_information"] = {};
  }
}

/**
 * Setup output tree.
 *
 * Branches are only enabled if the output is actually needed.
 */
void AliAnalysisTaskTrackSkim::SetupTree()
{
  std::string outputName = GetOutputSlot(2)->GetContainer()->GetName();
  // Add version to tree name so we can keep track of the outputs.
  outputName += "_v" + std::to_string(fOutputVersion);
  fTreeSkim = new TTree(outputName.c_str(), outputName.c_str());
  // Ensure that the size of the trees in memory is reasonable
  int nEnabledTrees = 1;
  fTreeSkim->SetMaxVirtualSize(1.e+8/nEnabledTrees);

  // Event level properties
  fTreeSkim->Branch("run_number", &fRunNumber);
  fTreeSkim->Branch("trigger_bit_INT7", &fTriggerBitINT7);
  // pp specific
  // nothing
  // Pythia specific
  if (fIsPythia) {
    // Nothing
  }
  // PbPb specific
  // NOTE: We use the centrality selection as a proxy for PbPb because the beam type isn't reliable
  //       until after the tree is already setup.
  if ((fMinCent != -999) && (fMaxCent != -999)) {
    fTreeSkim->Branch("centrality", &fCentralityForTree);
    fTreeSkim->Branch("event_plane_V0M", &fEventPlaneV0MForTree);
    fTreeSkim->Branch("trigger_bit_central", &fTriggerBitCentral);
    fTreeSkim->Branch("trigger_bit_semi_central", &fTriggerBitSemiCentral);
  }

  // Track skim level
  // Always want the first set of kinematics variables.
  AddParticleKinematicsToMap(fFirstTrackCollectionKinematics);
  if (fIsPythia) {
    // Add further information for pythia.
    AddParticleLabelsToMap(fFirstTrackCollectionLabels);
    // NOTE: At detector level, the particle ID is apparently returned as 0.
    //       So we skip storing it here by removing it.
    fFirstTrackCollectionLabels.erase("particle_ID");
    // NOTE: We also skip the encoded information since we don't have a use for it
    //       for data/det level as of April 2023.
    // NOTE: Only remove if it's already there...
    if (fEncodeAdditionalParticleInformation) {
      fFirstTrackCollectionLabels.erase("encoded_information");
    }
    AddParticleKinematicsToMap(fSecondTrackCollectionKinematics);
    AddParticleLabelsToMap(fSecondTrackCollectionLabels);
  }

  // Add all of the track variables that we've created into the tree, defining the branches.
  // NOTE: We access the value via the map (even though we have access in the iterator) to ensure
  //       we get the right memory address (instead of the iterator memory address). We could also dereference
  //       it, but I think this is clearer.
  // First track collection
  for (auto const & p : fFirstTrackCollectionKinematics) {
    fTreeSkim->Branch(("particle_data_" + p.first).c_str(), &fFirstTrackCollectionKinematics[p.first]);
  }
  if (fIsPythia) {
    for (auto const & p : fFirstTrackCollectionLabels) {
      fTreeSkim->Branch(("particle_data_" + p.first).c_str(), &fFirstTrackCollectionLabels[p.first]);
    }
    // Second track collection (for particle level)
    for (auto const & p : fSecondTrackCollectionKinematics) {
      fTreeSkim->Branch(("particle_gen_" + p.first).c_str(), &fSecondTrackCollectionKinematics[p.first]);
    }
    for (auto const & p : fSecondTrackCollectionLabels) {
      fTreeSkim->Branch(("particle_gen_" + p.first).c_str(), &fSecondTrackCollectionLabels[p.first]);
    }
  }

  // Print out track variables that were added. This is generally for debugging, but is also good for logging purposes.
  std::cout << "Post tree creation.\n";
  std::cout << *this;
  std::cout << "Jet substructure tree:\n";
  fTreeSkim->Print();
}

/**
 * Create output objects.
 */
void AliAnalysisTaskTrackSkim::UserCreateOutputObjects()
{
  // First call the base class
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Check that the task was initialized
  if (!fConfigurationInitialized) {
    AliFatal("Task was not initialized. Please ensure that Initialize() was called!");
  }
  // Reinitialize the YAML configuration
  fYAMLConfig.Reinitialize();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fNAcceptedFirstTrackCollection = new TH1F("fNAcceptedFirstTrackCollection", "fNAcceptedFirstTrackCollection", 4000, -0.5, 3999.5);
  fOutput->Add(fNAcceptedFirstTrackCollection);

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i = 0; i < fOutput->GetEntries(); ++i) {
    TH1* h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1) {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if (hn)
      hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  // Setup the tree output.
  SetupTree();

  PostData(1, fOutput);
  PostData(2, fTreeSkim);
}

Bool_t AliAnalysisTaskTrackSkim::IsEventSelected()
{
  if (!AliAnalysisTaskEmcalJet::IsEventSelected()) {
    return false;
  }

  // Centrality selection (only if specified)
  // Performing here rather than via AliEventCuts because I don't want to deal with manual configuration of AliEventCuts.
  if ((fMinCent != -999) && (fMaxCent != -999)) {
    if ((fCent > fMaxCent) || (fCent < fMinCent)) {
      return false;
    }
  }

  return true;
}

Bool_t AliAnalysisTaskTrackSkim::Run()
{
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().

  return kTRUE;
}

void AliAnalysisTaskTrackSkim::ClearTree()
{
  // For simplicity, we reset regardless of whether we're storing them.
  fTriggerBitINT7 = false;
  fCentralityForTree = -1;
  fEventPlaneV0MForTree = -1;
  fTriggerBitCentral = false;
  fTriggerBitSemiCentral = false;
  // Reset all of the variables in the map.
  for (auto & p : fFirstTrackCollectionKinematics) {
    p.second.clear();
  }
  for (auto & p : fFirstTrackCollectionLabels) {
    p.second.clear();
  }
  // Second track collection (for particle level)
  for (auto & p : fSecondTrackCollectionKinematics) {
    p.second.clear();
  }
  for (auto & p : fSecondTrackCollectionLabels) {
    p.second.clear();
  }
}

std::string AliAnalysisTaskTrackSkim::GetMCHeaderName()
{
  AliGenEventHeader * eventHeader = MCEvent()->GenEventHeader();
  return eventHeader->ClassName();
}

/**
  * NOTE: We encode only for particle level because this won't work for detector level
  */
void AliAnalysisTaskTrackSkim::EncodeAdditionalParticleLevelInformation(AliVParticle * particle)
{
  std::bitset<32> information;

  // Check photon types
  AliMCAnalysisUtils * mcAnalysisUtils = GetMCAnalysisUtils();
  TString headerName = GetMCHeaderName();
  // Using 1 for cluster E as convention since I don't think it's relevant for particle level.
  Int_t tag = GetMCAnalysisUtils()->CheckOrigin(particle->GetLabel(), MCEvent(), headerName, 1.);
  if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPrompt)) {
    information.set(EncodedInformation_t::kPhotonPrompt);
  }
  if (GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCFragmentation)) {
    information.set(EncodedInformation_t::kPhotonFragmentation);
  }

  // Finally, store the encoded information
  fSecondTrackCollectionLabels["encoded_information"].push_back(static_cast<int>(information.to_ulong()));
}

Bool_t AliAnalysisTaskTrackSkim::FillHistograms()
{
  // Reset
  ClearTree();

  // Generally, we try to set everything within reason for simplicity, even if we're not storing it.
  // Event level properties
  // Triggers
  ULong64_t fTriggerMask = static_cast<AliVAODHeader*>(InputEvent()->GetHeader())->GetOfflineTrigger();
  fTriggerBitINT7 = static_cast<bool>(fTriggerMask & AliVEvent::kINT7);
  fTriggerBitCentral = static_cast<bool>(fTriggerMask & AliVEvent::kCentral);
  fTriggerBitSemiCentral = static_cast<bool>(fTriggerMask & AliVEvent::kSemiCentral);
  // PbPb properties
  fCentralityForTree = static_cast<float>(fCent);
  // NOTE: Uncorrected.
  fEventPlaneV0MForTree = static_cast<float>(fEPV0);

  // Get particles
  auto particles = GetParticleContainer(0);
  if (!particles) {
    AliErrorStream() << GetName() << ": Unable to retrieve particles in slot 0!\n";
    return false;
  }

  for (auto particle : particles->accepted()) {
    if(fChargeScaling && particle->Charge() < 0) {
        fFirstTrackCollectionKinematics["pt"].push_back(particle->Pt()*-1);
    }
    else {
        fFirstTrackCollectionKinematics["pt"].push_back(particle->Pt());
    }
    fFirstTrackCollectionKinematics["eta"].push_back(particle->Eta());
    fFirstTrackCollectionKinematics["phi"].push_back(particle->Phi());
    if (fIsPythia) {
      // At detector level, the particle ID is apparently returned as 0.
      // So we skip storing it here. If we want the particle ID, we have to match
      // to the particle with the label.
      // Label
      fFirstTrackCollectionLabels["label"].push_back(particle->GetLabel());
    }
  }
  fNAcceptedFirstTrackCollection->Fill(fFirstTrackCollectionKinematics["pt"].size());
  // If we have no accepted tracks, don't bother filling - it's just a (small) waste of space.
  if (fFirstTrackCollectionKinematics["pt"].size() == 0) {
    return false;
  }

  if (fIsPythia) {
    auto particlesGen = GetParticleContainer(1);
    if (!particlesGen) {
      AliErrorStream() << GetName() << ": Unable to retrieve particles in slot 1!\n";
      return false;
    }

    for (auto particle : particlesGen->accepted()) {
      fSecondTrackCollectionKinematics["pt"].push_back(particle->Pt());
      fSecondTrackCollectionKinematics["eta"].push_back(particle->Eta());
      fSecondTrackCollectionKinematics["phi"].push_back(particle->Phi());
      // PDG ID
      fSecondTrackCollectionLabels["particle_ID"].push_back(particle->PdgCode());
      // Label
      fSecondTrackCollectionLabels["label"].push_back(particle->GetLabel());
      // Encode additional particle level information
      EncodeAdditionalParticleLevelInformation(particle);
    }
  }

  // Store the output
  fTreeSkim->Fill();

  return true;
}

/**
 * Add task to an existing analysis manager.
 *
 * @param[in] suffix Suffix string to attach to the task name
 */
AliAnalysisTaskTrackSkim* AliAnalysisTaskTrackSkim::AddTaskTrackSkim(const std::string& suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliErrorClass("No analysis manager to connect to.");
    return nullptr;
  }

  // Setup task name
  std::string taskName = "AliAnalysisTaskTrackSkim";
  std::string suffixName(suffix);
  if (suffixName != "") {
    taskName += "_";
    taskName += suffixName;
  }

  // Create task and configure as desired.
  AliAnalysisTaskTrackSkim* task = new AliAnalysisTaskTrackSkim(taskName.c_str());
  // Set general defaults.
  task->SetNCentBins(5);

  // Configuration is via YAML, so nothing else to be done here.
  mgr->AddTask(task);

  // Create containers for input/output
  std::string listOutputName = taskName;
  std::string treeOutputName = taskName + "_tree";

  std::string outputFile = AliAnalysisManager::GetCommonFileName();
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  AliAnalysisDataContainer* outputContainer = mgr->CreateContainer(
   listOutputName.c_str(), TList::Class(), AliAnalysisManager::kOutputContainer, outputFile.c_str());
  mgr->ConnectOutput(task, 1, outputContainer);

  AliAnalysisDataContainer* coutput2 = mgr->CreateContainer(treeOutputName.c_str(), TTree::Class(),
                               AliAnalysisManager::kOutputContainer, outputFile.c_str());
  mgr->ConnectOutput(task, 2, coutput2);

  return task;
}

/**
 * Prints information about the task.
 *
 * @return std::string containing information about the task.
 */
std::string AliAnalysisTaskTrackSkim::toString() const
{
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  tempSS << "Track skim analysis task:\n";
  tempSS << "General:\n";
  tempSS << "\tCentrality selection: ";
  if ((fMinCent != -999) && (fMaxCent != -999)) {
    tempSS << "(" << fMinCent << ", " << fMaxCent << ")\n";
  }
  else {
    tempSS << "disabled.\n";
  }
  tempSS << "\tCharge scaling: " << fChargeScaling << "\n";
  tempSS << "\tEncode additional particle information: " << fEncodeAdditionalParticleInformation << "\n";
  if (fEncodeAdditionalParticleInformation) {
    tempSS << "\t\tMC analysis utils generator: " << fMCAnalysisUtilsGenerator << "\n";
  }
  // Particle variables:
  tempSS << "Particle variables:\n";
  tempSS << "First track collection:\n";
  for (const auto & p : fFirstTrackCollectionKinematics) {
    tempSS << "\t" << p.first << "\n";
  }
  for (const auto & p : fFirstTrackCollectionLabels) {
    tempSS << "\t" << p.first << "\n";
  }
  tempSS << "Second track collection: \n";
  if (fSecondTrackCollectionKinematics.size()) {
    for (const auto & p : fSecondTrackCollectionKinematics) {
      tempSS << "\t" << p.first << "\n";
    }
    for (const auto & p : fSecondTrackCollectionLabels) {
      tempSS << "\t" << p.first << "\n";
    }
  }
  else {
    tempSS << "\tDisabled\n";
  }
  return tempSS.str();
}

/**
 * Print task information on an output stream using the string representation provided by
 * AliAnalysisTaskTrackSkim::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream& AliAnalysisTaskTrackSkim::Print(std::ostream& in) const
{
  in << toString();
  return in;
}

/**
 * Print task information using the string representation provided by
 * AliAnalysisTaskTrackSkim::toString
 *
 * @param[in] opt Unused
 */
void AliAnalysisTaskTrackSkim::Print(Option_t* opt) const { Printf("%s", toString().c_str()); }

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

/**
 * Hardest kt analysis task.
 */

/**
 * Implementation of the output stream operator for AliAnalysisTaskTrackSkim. Printing
 * basic task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim& myTask)
{
  std::ostream& result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 *
 * NOTE: We don't swap the base class values because the base class doesn't implement swap.
 */
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim& first,
     PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim& second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fYAMLConfig, second.fYAMLConfig);
  swap(first.fConfigurationInitialized, second.fConfigurationInitialized);
  swap(first.fChargeScaling, second.fChargeScaling);
  swap(first.fEncodeAdditionalParticleInformation, second.fEncodeAdditionalParticleInformation);
  swap(first.fMCAnalysisUtilsGenerator, second.fMCAnalysisUtilsGenerator);
  swap(first.fMCAnalysisUtils, second.fMCAnalysisUtils);
  swap(first.fNAcceptedFirstTrackCollection, second.fNAcceptedFirstTrackCollection);
  swap(first.fTriggerBitINT7, second.fTriggerBitINT7),
  swap(first.fCentralityForTree, second.fCentralityForTree),
  swap(first.fEventPlaneV0MForTree, second.fEventPlaneV0MForTree),
  swap(first.fTriggerBitCentral, second.fTriggerBitCentral),
  swap(first.fTriggerBitSemiCentral, second.fTriggerBitSemiCentral),
  swap(first.fFirstTrackCollectionKinematics, second.fFirstTrackCollectionKinematics),
  swap(first.fFirstTrackCollectionLabels, second.fFirstTrackCollectionLabels),
  swap(first.fSecondTrackCollectionKinematics, second.fSecondTrackCollectionKinematics),
  swap(first.fSecondTrackCollectionLabels, second.fSecondTrackCollectionLabels),
  swap(first.fTreeSkim, second.fTreeSkim);
}
