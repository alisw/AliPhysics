/**
 * Hardest kt jet substructure task. Adapted from AliAnalysisTaskJetDynamicalGrooming
 */

#include "AliAnalysisTaskJetHardestKt.h"

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
#include <AliVTrack.h>

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliEmcalPythiaInfo.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliRhoParameter.h"

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt);
/// \endcond

namespace PWGJE
{
namespace EMCALJetTasks
{

/**
 * Hardest kt grooming analysis task.
 */
const std::map<std::string, AliAnalysisTaskJetHardestKt::JetShapeType_t> AliAnalysisTaskJetHardestKt::fgkJetShapeTypeMap = {
  { "kMCTrue", AliAnalysisTaskJetHardestKt::kMCTrue },
  { "kTrueDet", AliAnalysisTaskJetHardestKt::kTrueDet },
  { "kData", AliAnalysisTaskJetHardestKt::kData },
  { "kDetEmb", AliAnalysisTaskJetHardestKt::kDetEmb },
  { "kDetEmbPart", AliAnalysisTaskJetHardestKt::kDetEmbPart },
  { "kPythiaDef", AliAnalysisTaskJetHardestKt::kPythiaDef },
  { "kDetEmbPartPythia", AliAnalysisTaskJetHardestKt::kDetEmbPartPythia },
  { "kGenOnTheFly", AliAnalysisTaskJetHardestKt::kGenOnTheFly }
};
const std::map<std::string, AliAnalysisTaskJetHardestKt::JetShapeSub_t> AliAnalysisTaskJetHardestKt::fgkJetShapeSubMap = {
  { "kNoSub", AliAnalysisTaskJetHardestKt::kNoSub },
  { "kConstSub", AliAnalysisTaskJetHardestKt::kConstSub },
  { "kDerivSub", AliAnalysisTaskJetHardestKt::kDerivSub },
  { "kEventSub", AliAnalysisTaskJetHardestKt::kEventSub }
};
const std::map<std::string, AliAnalysisTaskJetHardestKt::JetSelectionType_t> AliAnalysisTaskJetHardestKt::fgkJetSelectionMap = {
  { "kInclusive", AliAnalysisTaskJetHardestKt::kInclusive },
  { "kRecoil", AliAnalysisTaskJetHardestKt::kRecoil }
};
const std::map<std::string, AliAnalysisTaskJetHardestKt::DerivSubtrOrder_t> AliAnalysisTaskJetHardestKt::fgkDerivSubtrOrderMap = {
  { "kSecondOrder", AliAnalysisTaskJetHardestKt::kSecondOrder },
  { "kFirstOrder", AliAnalysisTaskJetHardestKt::kFirstOrder }
};
const std::map<std::string, AliAnalysisTaskJetHardestKt::GroomingMethod_t> AliAnalysisTaskJetHardestKt::fgkGroomingMethodMap = {
  { "kLeadingKt", AliAnalysisTaskJetHardestKt::kLeadingKt },
  { "kDynamicalZ", AliAnalysisTaskJetHardestKt::kDynamicalZ },
  { "kDynamicalKt", AliAnalysisTaskJetHardestKt::kDynamicalKt },
  { "kDynamicalTime", AliAnalysisTaskJetHardestKt::kDynamicalTime },
  { "kDynamicalCore", AliAnalysisTaskJetHardestKt::kDynamicalCore }
};

/**
 * Default constructor.
 */
AliAnalysisTaskJetHardestKt::AliAnalysisTaskJetHardestKt()
 : AliAnalysisTaskEmcalJet("AliAnalysisTaskJetHardestKt", kTRUE),
  fYAMLConfig(),
  fConfigurationInitialized(false),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fCheckResolution(kFALSE),
  fSubjetCutoff(0.1),
  fMinPtConst(1),
  fGroomingMethod(GroomingMethod_t::kLeadingKt),
  fHardCutoff(0),
  fDoTwoTrack(kFALSE),
  fCutDoubleCounts(kTRUE),
  fPhiCutValue(0.02),
  fEtaCutValue(0.02),
  fMagFieldPolarity(1),
  fDerivSubtrOrder(0),
  fStoreDetLevelJets(kFALSE),
  fEnableSubjetMatching(false),
  fSubstructureVariables(),
  fPtJet(nullptr),
  fHCheckResolutionSubjets(nullptr),
  fTreeSubstructure(nullptr)
{
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/**
 * Standard constructor.
 */
AliAnalysisTaskJetHardestKt::AliAnalysisTaskJetHardestKt(const char* name)
 : AliAnalysisTaskEmcalJet(name, kTRUE),
  fYAMLConfig(),
  fConfigurationInitialized(false),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fCheckResolution(kFALSE),
  fSubjetCutoff(0.1),
  fMinPtConst(1),
  fGroomingMethod(GroomingMethod_t::kLeadingKt),
  fHardCutoff(0),
  fDoTwoTrack(kFALSE),
  fCutDoubleCounts(kTRUE),
  fPhiCutValue(0.02),
  fEtaCutValue(0.02),
  fMagFieldPolarity(1),
  fDerivSubtrOrder(0),
  fStoreDetLevelJets(kFALSE),
  fEnableSubjetMatching(false),
  fSubstructureVariables(),
  fPtJet(nullptr),
  fHCheckResolutionSubjets(nullptr),
  fTreeSubstructure(nullptr)
{
  // Standard constructor.
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/**
 * Copy constructor.
 *
 * Copying this task will not copy the splittings objects or output hists and trees.
 * It is up to the user to create this objects use `UserCreateOutputObjects(...)`.
 */
AliAnalysisTaskJetHardestKt::AliAnalysisTaskJetHardestKt(
 const AliAnalysisTaskJetHardestKt& other)
 : fYAMLConfig(other.fYAMLConfig), fConfigurationInitialized(other.fConfigurationInitialized),
  fMinFractionShared(other.fMinFractionShared),
  fJetShapeType(other.fJetShapeType),
  fJetShapeSub(other.fJetShapeSub),
  fJetSelection(other.fJetSelection),
  fPtThreshold(other.fPtThreshold),
  fRMatching(other.fRMatching),
  fCentSelectOn(other.fCentSelectOn),
  fCentMin(other.fCentMin),
  fCentMax(other.fCentMax),
  fCheckResolution(other.fCheckResolution),
  fSubjetCutoff(other.fSubjetCutoff),
  fMinPtConst(other.fMinPtConst),
  fGroomingMethod(other.fGroomingMethod),
  fHardCutoff(other.fHardCutoff),
  fDoTwoTrack(other.fDoTwoTrack),
  fCutDoubleCounts(other.fCutDoubleCounts),
  fPhiCutValue(other.fPhiCutValue),
  fEtaCutValue(other.fEtaCutValue),
  fMagFieldPolarity(other.fMagFieldPolarity),
  fDerivSubtrOrder(other.fDerivSubtrOrder),
  fStoreDetLevelJets(other.fStoreDetLevelJets),
  fEnableSubjetMatching(other.fEnableSubjetMatching),
  fSubstructureVariables(other.fSubstructureVariables),
  fPtJet(nullptr),
  fHCheckResolutionSubjets(nullptr),
  fTreeSubstructure(nullptr)
{
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
AliAnalysisTaskJetHardestKt& AliAnalysisTaskJetHardestKt::operator=(
 AliAnalysisTaskJetHardestKt other)
{
  swap(*this, other);
  return *this;
}

/**
 * Retrieve task properties from the YAML configuration.
 */
void AliAnalysisTaskJetHardestKt::RetrieveAndSetTaskPropertiesFromYAMLConfig()
{
  // Same ordering as in the constructor (for consistency)
  std::string baseName = "jetSelection";
  std::string tempStr = "";
  fYAMLConfig.GetProperty({baseName, "minFractionShared"}, fMinFractionShared, false);
  bool res = fYAMLConfig.GetProperty({baseName, "jetShapeType"}, tempStr, false);
  if (res) {
    fJetShapeType = fgkJetShapeTypeMap.at(tempStr);
  }
  res = fYAMLConfig.GetProperty({baseName, "jetShapeSub"}, tempStr, false);
  if (res) {
    fJetShapeSub = fgkJetShapeSubMap.at(tempStr);
  }
  res = fYAMLConfig.GetProperty({baseName, "jetSelection"}, tempStr, false);
  if (res) {
    fJetSelection = fgkJetSelectionMap.at(tempStr);
  }
  fYAMLConfig.GetProperty({baseName, "ptThresholdForStoringOutput"}, fPtThreshold, false);
  fYAMLConfig.GetProperty({baseName, "maxMatchingDistance"}, fRMatching, false);

  // Options
  baseName = "general";
  fYAMLConfig.GetProperty({baseName, "centralitySelection"}, fCentSelectOn, false);
  // Only both setting the centrality range if it's enabled.
  if (fCentSelectOn) {
    std::pair<double, double> centRange;
    res = fYAMLConfig.GetProperty({ baseName, "centralityRange" }, centRange, false);
    if (res) {
      AliDebugStream(3) << "Setting centrality range of (" << centRange.first << ", " << centRange.second << ").\n";
      fCentMin = centRange.first;
      fCentMax = centRange.second;
    }
  }
  // Check subjet resolution
  fYAMLConfig.GetProperty({baseName, "checkSubjetResolution"}, fCheckResolution, false);
  fYAMLConfig.GetProperty({baseName, "angularSubjetResolutionCutoff"}, fSubjetCutoff, false);
  fYAMLConfig.GetProperty({baseName, "constituentPtCutoff"}, fMinPtConst, false);
  // Grooming properties
  res = fYAMLConfig.GetProperty({baseName, "groomingMethod"}, tempStr, false);
  if (res) {
    fGroomingMethod = fgkGroomingMethodMap.at(tempStr);
  }
  fYAMLConfig.GetProperty({baseName, "hardCutoff"}, fHardCutoff, false);
  fYAMLConfig.GetProperty({baseName, "considerTwoTrackEffects"}, fDoTwoTrack, false);
  // Default is true
  fYAMLConfig.GetProperty({baseName, "cutDoubleCounts"}, fCutDoubleCounts, false);

  // Two track effects settings
  baseName = "twoTrackEffects";
  fYAMLConfig.GetProperty({baseName, "phiCutValue"}, fPhiCutValue, false);
  fYAMLConfig.GetProperty({baseName, "etaCutValue"}, fEtaCutValue, false);
  fYAMLConfig.GetProperty({baseName, "magFieldPolarity"}, fMagFieldPolarity, false);

  // Back to general options
  baseName = "general";
  res = fYAMLConfig.GetProperty({baseName, "derivSubtrOrder"}, tempStr, false);
  if (res) {
    fDerivSubtrOrder = fgkDerivSubtrOrderMap.at(tempStr);
  }
  fYAMLConfig.GetProperty({baseName, "storeDetLevelJets"}, fStoreDetLevelJets, false);
  fYAMLConfig.GetProperty({baseName, "subjetMatching"}, fEnableSubjetMatching, false);
}

/**
 * Initialize task.
 */
bool AliAnalysisTaskJetHardestKt::Initialize()
{
  fConfigurationInitialized = false;

  // Ensure that we have at least one configuration in the YAML config.
  /*if (fYAMLConfig.DoesConfigurationExist(0) == false) {
    // No configurations exist. Return immediately.
    return fConfigurationInitialized;
  }*/

  // Always initialize for streaming purposes
  fYAMLConfig.Initialize();

  // Setup task based on the properties defined in the YAML config
  AliDebugStream(2) << "Configuring task from the YAML configuration.\n";
  RetrieveAndSetTaskPropertiesFromYAMLConfig();
  /*SetupJetContainersFromYAMLConfig();
  SetupParticleContainersFromYAMLConfig();
  SetupClusterContainersFromYAMLConfig();*/
  AliDebugStream(2) << "Finished configuring via the YAML configuration.\n";

  // Print the results of the initialization
  // Print outside of the ALICE Log system to ensure that it is always available!
  std::cout << *this;

  fConfigurationInitialized = true;
  return fConfigurationInitialized;
}

std::string AliAnalysisTaskJetHardestKt::GroomingMethodName() const
{
  // Retrieve the base name from the enumeration map. In principle, this isn't the fastest thing,
  // but it should be good enough.
  std::string tempGroomingMethod = GetKeyFromMapValue(fGroomingMethod, fgkGroomingMethodMap);
  // Transform name, removing the leading "k" and adding underscores at capitals.
  // ie: "kLeadingKt" -> "leading_kt"
  tempGroomingMethod = tempGroomingMethod.substr(1);
  std::string groomingMethod = "";
  // There's probably a much better way to do this, but it's quick and easy.
  for (std::size_t i = 0; i < tempGroomingMethod.size(); i++)
  {
    // We skip for the first letter because we don't want a leading "_"
    if (isupper(tempGroomingMethod[i]) && i != 0) {
      groomingMethod += "_";
    }
    // We use tolower() for everything just in case. It simplifies the logic and shouldn't hurt anything.
    groomingMethod += tolower(tempGroomingMethod[i]);
  }

  // Now determine the rest of the name.
  if (fHardCutoff > 0) {
    // Includes leading zero
    groomingMethod += "_z_cut_0";
    // Assumes that fHardCutoff is between 0 and 1.
    groomingMethod += std::to_string(static_cast<int>(fHardCutoff * 10));
  }
  return groomingMethod;
}

void AliAnalysisTaskJetHardestKt::AddSubstructureVariablesToMap(const std::string & prefix)
{
  std::string groomingMethod = GroomingMethodName();
  fSubstructureVariables[prefix + "_jet_pt"] = -1;
  fSubstructureVariables[prefix + "_leading_track_pt"] = -1;
  if (prefix == "data" && (fJetShapeType == kDetEmbPartPythia || fJetShapeType == kData)) {
    fSubstructureVariables[prefix + "_leading_track_pt_sub"] = -1;
  }
  fSubstructureVariables[groomingMethod + "_" + prefix + "_kt"] = -1;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_z"] = -1;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_delta_R"] = -1;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_n_to_split"] = -1;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_n_groomed_to_split"] = -1;
  // Equivalent to nsd
  fSubstructureVariables[groomingMethod + "_" + prefix + "_n_passed_grooming"] = -1;
}

/**
 * Setup output tree.
 *
 * Branches are only enabled if the output is actually needed.
 *
 * NOTE: This tree uses a number of more subtle ROOT tree functions. It relies heavily on information from
 * the ROOT users guide: https://root.cern.ch/root/htmldoc/guides/users-guide/Trees.html
 */
void AliAnalysisTaskJetHardestKt::SetupTree()
{
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeSubstructure = new TTree(nameoutput, nameoutput);

  // Define substructure variables in our map, which will then be used to define tree branches.
  // Data is always included.
  AddSubstructureVariablesToMap("data");

  // Matched and det level are only relevant for pythia or embedded pythia.
  if (fJetShapeType == kDetEmbPartPythia || fJetShapeType == kPythiaDef) {
    AddSubstructureVariablesToMap("matched");
    if (fStoreDetLevelJets) {
      AddSubstructureVariablesToMap("det_level");
    }
  }
  // Add pythia info from embedding.
  if (fJetShapeType == kDetEmbPartPythia) {
    fSubstructureVariables["pt_hard_bin"] = -1;
    fSubstructureVariables["pt_hard"] = -1;
    fSubstructureVariables["cross_section"] = -1;
    fSubstructureVariables["n_trials"] = -1;
  }
  // Alternatively if in pythia, we'll fill it out directly from the variables that are already available.
  if (fIsPythia) {
    // Will be automatically filled by AliAnalysisTaskEmcal.
    fTreeSubstructure->Branch("pt_hard_bin", fPtHardInitialized ? &fPtHardBinGlobal : &fPtHardBin);
    fTreeSubstructure->Branch("pt_hard", &fPtHard);
  }

  // Add appropriate subjet matching fields.
  std::string groomingMethod = GroomingMethodName();
  if (fEnableSubjetMatching) {
    if (fJetShapeType == kDetEmbPartPythia) {
      // Hybrid-det level matching (hybrid will be labeled as "data" usually, but we specialize here because we already have to specialize
      // to get the matching arguments correct).
      std::string name = groomingMethod;
      name += "_hybrid_det_level_matching_";
      fSubstructureVariables[name + "leading"] = -1;
      fSubstructureVariables[name + "subleading"] = -1;

      // Momentum fraction of det level subjets contained in the hybrid jet.
      fSubstructureVariables[name + "leading_pt_fraction_in_hybrid_jet"] = -1;
      fSubstructureVariables[name + "subleading_pt_fraction_in_hybrid_jet"] = -1;
    }
    if (fJetShapeType == kDetEmbPartPythia || fJetShapeType == kPythiaDef) {
      // det level-true level matching (true level will be labeled as "matched" usually, but we specialize here because we already have to specialize
      // to get the matching arguments correct).
      std::string name = groomingMethod;
      name += "_det_level_true_matching_";
      fSubstructureVariables[name + "leading"] = -1;
      fSubstructureVariables[name + "subleading"] = -1;
    }
  }

  // Add all of the substructure variables that we've created into the tree, defining the branches.
  for (auto const & p : fSubstructureVariables) {
    // NOTE: We access the value via the map (even though we have access in the iterator) to ensure
    //       we get the right memory address (instead of the iterator memory address).
    fTreeSubstructure->Branch(p.first.c_str(), &fSubstructureVariables[p.first], TString::Format("%s/F", p.first.c_str()));
  }

  // Print out substructure variables that were added. This is generally for debugging, but is also good for logging purposes.
  std::cout << "Post tree creation.\n";
  std::cout << *this;
  std::cout << "Jet substructure tree:\n";
  fTreeSubstructure->Print();
}

/**
 * Create output objects.
 */
void AliAnalysisTaskJetHardestKt::UserCreateOutputObjects()
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

  fPtJet = new TH1F("fPtJet", "fPtJet", 100, 0, 200);
  fOutput->Add(fPtJet);

  //// jetpt, deltaR, pt residual 1, pt residual 2, z residual
  const Int_t dimResol = 5;
  const Int_t nBinsResol[dimResol] = { 10, 10, 80, 80, 80 };
  const Double_t lowBinResol[dimResol] = { 0, 0, -1, -1, -1 };
  const Double_t hiBinResol[dimResol] = { 200, 0.3, 1, 1, 1 };
  fHCheckResolutionSubjets = new THnSparseF("fHCheckResolutionSubjets", "Mom.Resolution of Subjets vs opening angle",
                       dimResol, nBinsResol, lowBinResol, hiBinResol);
  fOutput->Add(fHCheckResolutionSubjets);

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
  PostData(2, fTreeSubstructure);
}

Bool_t AliAnalysisTaskJetHardestKt::Run()
{
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().

  return kTRUE;
}

double AliAnalysisTaskJetHardestKt::DynamicalGrooming(const fastjet::PseudoJet& subjet1,
                               const fastjet::PseudoJet& subjet2,
                               const fastjet::PseudoJet& parent, const double R,
                               const double a) const
{
  double deltaR = subjet1.delta_R(subjet2);
  double z = 0;
  if (parent.pt() > 0) {
    z = subjet1.pt() / parent.pt();
  }
  return z * (1 - z) * parent.pt() * std::pow(deltaR / R, a);
}

double AliAnalysisTaskJetHardestKt::CalculateDynamicalZ(const fastjet::PseudoJet& subjet1,
                              const fastjet::PseudoJet& subjet2,
                              const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 0.1);
}

double AliAnalysisTaskJetHardestKt::CalculateDynamicalKt(const fastjet::PseudoJet& subjet1,
                              const fastjet::PseudoJet& subjet2,
                              const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 1);
}

double AliAnalysisTaskJetHardestKt::CalculateDynamicalTime(const fastjet::PseudoJet& subjet1,
                               const fastjet::PseudoJet& subjet2,
                               const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 2);
}

double AliAnalysisTaskJetHardestKt::CalculateDynamicalCore(const fastjet::PseudoJet& subjet1,
                               const fastjet::PseudoJet& subjet2,
                               const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 0.5);
}

Bool_t AliAnalysisTaskJetHardestKt::FillHistograms()
{
  // Container zero is always the base container: the data container, the
  // embedded subtracted in the case of embedding or the detector level in case
  // of pythia
  AliJetContainer* jetCont = GetJetContainer(0);
  // The jetR must be the same for jet collections that we consider, so we just retrieve it once.
  double jetR = jetCont->GetJetRadius();

  if (fCentSelectOn) {
    if ((fCent > fCentMax) || (fCent < fCentMin)) {
      return 0;
    }
  }

  std::string groomingMethod = GroomingMethodName();
  if (jetCont) {
    for (auto jet1 : jetCont->accepted()) {
      if (!jet1) {
        continue;
      }
      AliEmcalJet* jet2 = nullptr;
      AliEmcalJet* jet3 = nullptr;
      fPtJet->Fill(jet1->Pt());
      AliEmcalJet* jetUS = nullptr;
      Int_t ifound = 0, jfound = 0;
      Int_t ilab = -1, jlab = -1;

      // Clear out previously stored values. I suspect that the tree would clear this out, but
      // it is safer to ensure that is the case by doing it ourselves.
      for (const auto & p : fSubstructureVariables) {
        fSubstructureVariables[p.first] = -1;
      }

      // The embedding mode
      // the matching is done between unsubtracted embedded jets and detector
      // level jets unsubtracted and subtracted jets share the label. Once we
      // identify the corresponding unsubtracted jet, jetUS, then we fetch jet2,
      // which is the matched detector level jet In the case we are not
      // considering constituent subtraction, then the detector-level matched jet
      // is the one that was directly matched to the base jet1. Then, the
      // particle-level jet jet3 is obtained as the matched one to jet2 In short,
      // there are 2 consecutive matchings, between particle-level (jet3) and
      // detector-level (jet2) pythia jets and between jet2 and the embedding
      // unsubtracted jet. Note that the matching obtained via ClosestJet is
      // purely geometrical. So below we require a fraction of the probe momentum
      // to be reconstructed in the embedded jet.
      if (fJetShapeType == kDetEmbPartPythia) {
        AliJetContainer* jetContUS = GetJetContainer(2);

        if (fJetShapeSub == kConstSub) {
          for (Int_t i = 0; i < jetContUS->GetNJets(); i++) {
            jetUS = jetContUS->GetJet(i);
            if (jetUS->GetLabel() == jet1->GetLabel()) {
              ifound++;
              if (ifound == 1) {
                ilab = i;
              }
            }
          }
          if (ilab == -1) {
            continue;
          }
          jetUS = jetContUS->GetJet(ilab);
          jet2 = jetUS->ClosestJet();
        }

        if (fJetShapeSub == kEventSub) {
          jetUS = jet1->ClosestJet();
          if (!jetUS) {
            continue;
          }
          jet2 = jetUS->ClosestJet();
        }

        if (!(fJetShapeSub == kConstSub) && !(fJetShapeSub == kEventSub)) {
          jet2 = jet1->ClosestJet();
        }

        if (!jet2) {
          AliDebugStream(3) << "jet2 does not exist, returning\n";
          continue;
        }

        // AliJetContainer *jetContPart=GetJetContainer(3);
        jet3 = jet2->ClosestJet();

        if (!jet3) {
          AliDebugStream(3) << "jet3 does not exist, returning\n";
          continue;
        }
        AliDebugStream(3) << "jet 3 exists" << jet3->Pt() << "\n";

        if (fCheckResolution) {
          CheckSubjetResolution(jet2, jet3);
        }

        Double_t fraction = 0;
        if (!(fJetShapeSub == kConstSub) && !(fJetShapeSub == kEventSub)) {
          fraction = jetCont->GetFractionSharedPt(jet1);
        }
        if ((fJetShapeSub == kConstSub) || (fJetShapeSub == kEventSub)) {
          fraction = jetContUS->GetFractionSharedPt(jetUS);
        }

        if (fraction < fMinFractionShared) {
          continue;
        }

        // Lastly, add the pythia information from the embedding helper if available.
        const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
        if (embeddingHelper) {
          fSubstructureVariables["pt_hard_bin"] = embeddingHelper->GetPtHardBin();
          fSubstructureVariables["pt_hard"] = embeddingHelper->GetPythiaPtHard();
          fSubstructureVariables["cross_section"] = embeddingHelper->GetPythiaXSection();
          fSubstructureVariables["n_trials"] = embeddingHelper->GetPythiaTrials();
        }
      }

      // this is the mode to run over pythia to produce a det-part response
      // here we have also added the constituent-subtraction case, but we don't
      // use it normally in pp the matching is purely geometrical
      if (fJetShapeType == kPythiaDef) {
        AliJetContainer* jetContUS = GetJetContainer(2);
        AliJetContainer* jetContPart = GetJetContainer(3);

        if (fJetShapeSub == kConstSub) {
          for (Int_t i = 0; i < jetContUS->GetNJets(); i++) {
            jetUS = jetContUS->GetJet(i);
            if (jetUS->GetLabel() == jet1->GetLabel()) {
              ifound++;
              if (ifound == 1)
                ilab = i;
            }
          }
          if (ilab == -1)
            continue;
          jetUS = jetContUS->GetJet(ilab);
          jet2 = jetUS->ClosestJet();

          if (!jet2) {
            AliDebugStream(3) << "jet2 does not exist, returning\n";
            continue;
          }

          for (Int_t j = 0; j < jetContPart->GetNJets(); j++) {
            jet3 = jetContPart->GetJet(j);
            if (!jet3)
              continue;
            if (jet3->GetLabel() == jet2->GetLabel()) {
              jfound++;
              if (jfound == 1)
                jlab = j;
            }
          }
          if (jlab == -1)
            continue;
          jet3 = jetContPart->GetJet(jlab);
          if (!jet3) {
            AliDebugStream(3) << "jet3 does not exist, returning\n";
            continue;
          }
        }
        if (!(fJetShapeSub == kConstSub))
          jet3 = jet1->ClosestJet();
        if (!jet3) {
          AliDebugStream(3) << "jet3 does not exist, returning\n";
          continue;
        }

        if (fCheckResolution) {
          CheckSubjetResolution(jet1, jet3);
        }
      }

      Double_t ptSubtracted = 0;
      if (fJetShapeSub == kConstSub || fJetShapeSub == kEventSub) {
        ptSubtracted = jet1->Pt();
      }

      else if (fJetShapeSub == kDerivSub) {
        ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
      }

      else if (fJetShapeSub == kNoSub) {
        if ((fJetShapeType == kData) || (fJetShapeType == kDetEmbPartPythia)) {
          ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
        }
        else if ((fJetShapeType == kPythiaDef) || (fJetShapeType == kMCTrue) || (fJetShapeType == kGenOnTheFly)) {
          ptSubtracted = jet1->Pt();
        }
      }

      if (ptSubtracted < fPtThreshold) {
        continue;
      }

      if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1)) {
        continue;
      }

      fSubstructureVariables["data_jet_pt"] = ptSubtracted;

      // The double counting cut should be applied to the unsubtracted hybrid max track pt. Confusingly,
      // MaxTrackPt() called on the constituent subtracted jet returns the _unsubtracted_ jet pt. Here, we
      // store the unsubtracted max track pt. Then, during setup for declustering, we'll store the leading
      // subtracted constituent. That way, we'll be able to apply the double counting cut offline, regardless
      // of the whether we want to use the subtracted or unsubtracted case.
      //
      // We also maintain the ability to apply the double counting cut during analysis.
      fSubstructureVariables["data_leading_track_pt"] = jet1->MaxTrackPt();
      if (fCutDoubleCounts == kTRUE && fJetShapeType == kDetEmbPartPythia) {
        // Cut if the leading hybrid track pt is greater than the leading detector level track pt.
        if (jet1->MaxTrackPt() > jet2->MaxTrackPt()) {
          continue;
        }
      }

      // Determine the splittings of the main (data) jet.
      std::shared_ptr<SelectedSubjets> dataSubjets = IterativeParents(jet1, "data", true, jetR);

      // If appropriate, then fill the matched and/or detector level jets.
      std::shared_ptr<SelectedSubjets> matchedSubjets = nullptr;
      std::shared_ptr<SelectedSubjets> detLevelSubjets = nullptr;
      if (fJetShapeType == kPythiaDef) {
        fSubstructureVariables["matched_jet_pt"] = jet3->Pt();
        // NOTE: For the embedded cases, this uses the unsubtracted max track pt even though
        //       the jet pt is constituent subtracted.
        fSubstructureVariables["matched_leading_track_pt"] = jet3->MaxTrackPt();
        matchedSubjets = IterativeParents(jet3, "matched", false, jetR);
      }

      if (fJetShapeType == kDetEmbPartPythia) {
        fSubstructureVariables["matched_jet_pt"] = jet3->Pt();
        // NOTE: For the embedded cases, this uses the unsubtracted max track pt even though
        //       the jet pt is constituent subtracted.
        fSubstructureVariables["matched_leading_track_pt"] = jet3->MaxTrackPt();
        matchedSubjets = IterativeParents(jet3, "matched", false, jetR);
        if (fStoreDetLevelJets) {
          fSubstructureVariables["det_level_jet_pt"] = jet2->Pt();
          // NOTE: For the embedded cases, this uses the unsubtracted max track pt even though
          //       the jet pt is constituent subtracted.
          fSubstructureVariables["det_level_leading_track_pt"] = jet2->MaxTrackPt();
          detLevelSubjets = IterativeParents(jet2, "det_level", false, jetR);
        }
      }

      // Perform subjet matching if you're found splittings.
      if (fEnableSubjetMatching) {
        if (fJetShapeType == kDetEmbPartPythia) {
          // Hybrid-det level matching
          // Before June 2020, we've used distance based matching because the subtracted constituents had a new global index
          // (since it's a new track collection). However, I eventually realized that the track label is passed onto the subtracted
          // constituents. So we can take advantage of the label to perform hybrid-det level matching.
          StoreSubjetMatching(detLevelSubjets, dataSubjets, false, "hybrid_det_level_matching");
          // det level-part matching
          StoreSubjetMatching(matchedSubjets, detLevelSubjets, false, "det_level_true_matching");

          // Check ultimate location of the det level subjets.
          // Did they end up in the hybrid jet at all?
          SubjetsInHybridJet(detLevelSubjets, jet1);
        }

        if (fJetShapeType == kPythiaDef) {
          // det level-part matching
          StoreSubjetMatching(matchedSubjets, dataSubjets, false, "det_level_true_matching");
        }
      }

      fTreeSubstructure->Fill();
    }
  }

  return kTRUE;
}

void AliAnalysisTaskJetHardestKt::StoreSubjetMatching(const std::shared_ptr<SelectedSubjets> & generatorLikeSubjets, const std::shared_ptr<SelectedSubjets> & measuredLikeSubjets, bool matchUsingDistance, std::string matchingPrefix)
{
  int leadingSubjetStatus = -1, subleadingSubjetStatus = -1;
  // We need subjets from both det level and hybrid to perform the matching. If both are missing, it defaults to -1.
  if (generatorLikeSubjets && measuredLikeSubjets) {
    // Leading
    if (SubjetContainedInSubjet(generatorLikeSubjets->leading, generatorLikeSubjets->leadingConstituents, measuredLikeSubjets->leading, measuredLikeSubjets->leadingConstituents, matchUsingDistance)) {
      leadingSubjetStatus = 1;
    }
    else if (SubjetContainedInSubjet(generatorLikeSubjets->leading, generatorLikeSubjets->leadingConstituents, measuredLikeSubjets->subleading, measuredLikeSubjets->subleadingConstituents, matchUsingDistance)) {
      leadingSubjetStatus = 2;
    }
    else {
      leadingSubjetStatus = 3;
    }
    // Subleading
    if (SubjetContainedInSubjet(generatorLikeSubjets->subleading, generatorLikeSubjets->subleadingConstituents, measuredLikeSubjets->subleading, measuredLikeSubjets->subleadingConstituents, matchUsingDistance)) {
      subleadingSubjetStatus = 1;
    }
    else if (SubjetContainedInSubjet(generatorLikeSubjets->subleading, generatorLikeSubjets->subleadingConstituents, measuredLikeSubjets->leading, measuredLikeSubjets->leadingConstituents, matchUsingDistance)) {
      subleadingSubjetStatus = 2;
    }
    else {
      subleadingSubjetStatus = 3;
    }
  }
  else {
    if (measuredLikeSubjets) {
      leadingSubjetStatus = 0;
      subleadingSubjetStatus = 0;
    }
  }
  //std::cout << "Leading status=" << leadingSubjetStatus << "\n";

  // Store the results.
  std::string groomingMethod = GroomingMethodName();
  fSubstructureVariables[groomingMethod + "_" + matchingPrefix + "_leading"] = leadingSubjetStatus;
  fSubstructureVariables[groomingMethod + "_" + matchingPrefix + "_subleading"] = subleadingSubjetStatus;
}

/**
 * Check whether detector level (generator-like) subjets are stored in the hybrid jet.
 * This enables the study of what happens to subjets which aren't matched.
 *
 * Note: Since the detector level subjet could in principle be far into the subjet, just because
 *       the hybrid and detector level jets were matched doesn't mean that the detector level subjet
 *       will trivially be contained in the hybrid.
 *
 * @param[in] generatorLikeSubjets Detector level selected subjets.
 * @param[in] hybridJet Hybrid jet where the subjets are expected to be contained.
 *
 * @returns None. The values are stored in the tree.
 */
void AliAnalysisTaskJetHardestKt::SubjetsInHybridJet(const std::shared_ptr<SelectedSubjets> & generatorLikeSubjets, AliEmcalJet* hybridJet)
{
  std::vector<fastjet::PseudoJet> hybridConstituents;
  fastjet::PseudoJet pseudoTrack;
  for (int constituentIndex = 0; constituentIndex < hybridJet->GetNumberOfTracks(); constituentIndex++) {
    AliVParticle* part = hybridJet->Track(constituentIndex);
    if (!part) {
      continue;
    }
    // Set the PseudoJet and add it to the inputs.
    pseudoTrack.reset(part->Px(), part->Py(), part->Pz(), part->E());
    pseudoTrack.set_user_index(GetConstituentID(constituentIndex, part, hybridJet));

    hybridConstituents.push_back(pseudoTrack);
  }

  float leadingPtFractionInHybrid = 0;
  float subleadingPtFractionInHybrid = 0;
  if (generatorLikeSubjets) {
    // We can match by label because the det level labels are propagated.
    leadingPtFractionInHybrid = SubjetSharedMomentum(generatorLikeSubjets->leadingConstituents, hybridConstituents, false) / generatorLikeSubjets->leading.pt();
    subleadingPtFractionInHybrid = SubjetSharedMomentum(generatorLikeSubjets->subleadingConstituents, hybridConstituents, false) / generatorLikeSubjets->subleading.pt();
  }

  // Store the momentum fraction.
  std::string groomingMethod = GroomingMethodName();
  fSubstructureVariables[groomingMethod + "_hybrid_det_level_matching_leading_pt_fraction_in_hybrid_jet"] = leadingPtFractionInHybrid;
  fSubstructureVariables[groomingMethod + "_hybrid_det_level_matching_subleading_pt_fraction_in_hybrid_jet"] = subleadingPtFractionInHybrid;
}

/**
 * Perform subjet matching between measured-like and generator-like jets.
 *
 * Using embedding as an example, we can have two possible combinations:
 * - measured-like is hybrid, generator-like is detector level.
 * - generator-like is detector level, generator-like is part level.
 *
 * @param[in] generatorLikeSubjet The generator-like subjet.
 * @param[in] generatorLikeSubjetConstituents Generator-like subjet constituents.
 * @param[in] measuredLikeSubjet The measured-like subjet. Currently unused.
 * @param[in] measuredLikeSubjetConstituents Measured-like subjet constituents.
 * @param[in] matchUsingDistance If true, matching using distance. If false, we match using `user_index()`.
 *
 * @returns True if the subjets match (meaning more than 50% of the generator-like constituents pt is in the measured-like subjet).
 */
bool AliAnalysisTaskJetHardestKt::SubjetContainedInSubjet(const fastjet::PseudoJet & generatorLikeSubjet, const std::vector<fastjet::PseudoJet> & generatorLikeSubjetConstituents,
                             const fastjet::PseudoJet & measuredLikeSubjet, const std::vector<fastjet::PseudoJet> & measuredLikeSubjetConstituents,
                             bool matchUsingDistance)
{
  double subjetSharedMomentum = SubjetSharedMomentum(generatorLikeSubjetConstituents, measuredLikeSubjetConstituents, matchUsingDistance);

  //std::cout << "fraction=" << subjetSharedMomentum / generatorLikeSubjet.pt() << "\n";
  if ((subjetSharedMomentum / generatorLikeSubjet.pt()) > 0.5) {
    return true;
  }
  return false;
}

/**
 * Calculate the fraction of subjet momentum in another set of constituents. This fraction can be stored or used to
 * determine whether one subjet is stored in another (see `SubjetContainedInSubjet(...)`).
 *
 * Using embedding as an example, we can have two possible combinations:
 * - measured-like is hybrid, generator-like is detector level.
 * - generator-like is detector level, generator-like is part level.
 *
 * @param[in] generatorLikeSubjetConstituents Generator-like subjet constituents.
 * @param[in] measuredLikeSubjetConstituents Measured-like subjet constituents.
 * @param[in] matchUsingDistance If true, matching using distance. If false, we match using `user_index()`.
 *
 * @returns True if the subjets match (meaning more than 50% of the generator-like constituents pt is in the measured-like subjet).
 */
double AliAnalysisTaskJetHardestKt::SubjetSharedMomentum(const std::vector<fastjet::PseudoJet> & generatorLikeSubjetConstituents,
                             const std::vector<fastjet::PseudoJet> & measuredLikeSubjetConstituents,
                             bool matchUsingDistance)
{
  double sumPt = 0;
  double delta = 0.01;

  for (const auto & generatorLikeConstituent : generatorLikeSubjetConstituents) {
    double generatorLikeEta = generatorLikeConstituent.eta();
    double generatorLikePhi = generatorLikeConstituent.phi();
    int generatorLikeIndex = generatorLikeConstituent.user_index();
    for (const auto & measuredLikeConstituent : measuredLikeSubjetConstituents) {
      if (matchUsingDistance) {
        double dEta = std::abs(measuredLikeConstituent.eta() - generatorLikeEta);
        if (dEta > delta) {
          continue;
        }
        double dPhi = std::abs(measuredLikeConstituent.phi() - generatorLikePhi);
        if (dPhi > delta) {
          continue;
        }
      }
      else {
        // Match using user_index()
        if (measuredLikeConstituent.user_index() != generatorLikeIndex) {
          continue;
        }
      }
      sumPt += generatorLikeConstituent.pt();
      // We've matched once - no need to match again.
      // Otherwise, the run the risk of summing a generator-like constituent pt twice.
      break;
    }
  }

  return sumPt;
}

/**
 * Get the constituent ID to be stored in `PseudoJet::user_index()`.
 *
 * Ensure they are uniquely identified (per type of jet) using the id. We don't use the global index (accessed via
 * `TrackAt(int)`) because they won't be the same for the subtracted and unsubtracted constituents (they are different
 * TClonesArrays). Instead, we label the part and det level with the MClabel (ie. for those which have it available),
 * and for the data (where it always equals -1), we use the global index (offset sufficiently) so it won't overlap with
 * other values.
 *
 * NOTE: We don't have the same restrictions in setting this value as for the general substructure extraction because
 *       we don't use this ID for matching subjets to constituents.
 *
 * @param[in] constituentIndex Index of the constituent in the collection where it's stored (for example, in a jet).
 * @param[in] part Jet constituent.
 * @param[in] jet Jet containing the constituent.
 *
 * @returns The unique constituent ID.
 */
int AliAnalysisTaskJetHardestKt::GetConstituentID(int constituentIndex, AliVParticle * part, AliEmcalJet * jet)
{
  // NOTE: Usually, we would use the global offset defined for the general subtracter extraction task. But we don't want to
  //       depend on that task, so we just define it here locally.
  int id = part->GetLabel() != -1 ? part->GetLabel() : (jet->TrackAt(constituentIndex) + 2000000);
  return id;
}

/**
 * Determine and iterate through jet splittings for a given jet. The output is stored in the given
 * jet substructure output container.
 *
 * @param[in] jet Jet to be declustered.
 * @param[in] prefix Prefix under which the jet splitting properties will be stored.
 * @param[in] isData If True, treat the splitting as coming from data. This means that ghosts are utilized and track resolution may be considered.
 */
std::shared_ptr<SelectedSubjets> AliAnalysisTaskJetHardestKt::IterativeParents(AliEmcalJet* jet, const std::string & prefix, bool isData, double jetR)
{
  AliDebugStream(1) << "Beginning iteration through the splittings.\n";
  std::vector<fastjet::PseudoJet> inputVectors;
  fastjet::PseudoJet pseudoTrack;
  double maxLeadingTrackPt = 0;
  for (int constituentIndex = 0; constituentIndex < jet->GetNumberOfTracks(); constituentIndex++) {
    AliVParticle* part = jet->Track(constituentIndex);
    if (!part) {
      continue;
    }
    if (isData == true && fDoTwoTrack == kTRUE && CheckClosePartner(jet, part)) {
      continue;
    }
    // Set the PseudoJet and add it to the inputs.
    pseudoTrack.reset(part->Px(), part->Py(), part->Pz(), part->E());
    pseudoTrack.set_user_index(GetConstituentID(constituentIndex, part, jet));
    inputVectors.push_back(pseudoTrack);

    if (part->Pt() > maxLeadingTrackPt) {
      maxLeadingTrackPt = part->Pt();
    }
  }
  // Store the max subtracted leading track pt. This doesn't have to be equal to the leading track pt for PbPb data or hybrid in the embedded case.
  if (isData && (fJetShapeType == kDetEmbPartPythia || fJetShapeType == kData)) {
    fSubstructureVariables[prefix + "_leading_track_pt_sub"] = maxLeadingTrackPt;
  }

  // Define substructure variables.
  float hardnessMeasure = -0.005;
  float kt = -0.005;
  float zg = -0.005;
  float rg = -0.005;
  int nToSplit = 0;
  int nGroomedToSplit = 0;
  int nPassedGrooming = 0;

  std::shared_ptr<SelectedSubjets> selectedSubjets;
  try {
    fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
    fastjet::JetDefinition jetDef(jetalgo, 1., fastjet::RecombinationScheme::E_scheme, fastjet::BestFJ30);
    // For area calculation (when desired)
    fastjet::GhostedAreaSpec ghost_spec(1, 1, 0.05);
    fastjet::AreaDefinition areaDef(fastjet::passive_area, ghost_spec);
    // We use a pointer for the CS because it has to stay in scope while we explore the splitting history.
    fastjet::ClusterSequence * cs = nullptr;
    if (isData) {
      cs = new fastjet::ClusterSequenceArea(inputVectors, jetDef, areaDef);
    }
    else {
      cs = new fastjet::ClusterSequence(inputVectors, jetDef);
    }
    std::vector<fastjet::PseudoJet> outputJets = cs->inclusive_jets(0);

    fastjet::PseudoJet jj;
    jj = outputJets[0];

    // Determine the substructure properties by performing the grooming.
    int nSplit = 0;
    // Setup
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    fastjet::PseudoJet leadingSelectedSubjet;
    fastjet::PseudoJet subleadingSelectedSubjet;
    while (jj.has_parents(j1, j2)) {
      // Setup
      nSplit += 1;
      if (j1.perp() < j2.perp()) {
        swap(j1, j2);
      }

      // Calculate the splitting properties.
      double xz = j2.perp() / (j2.perp() + j1.perp());
      double xdeltaR = j1.delta_R(j2);
      double xkt = j2.perp() * sin(xdeltaR);

      // First check the hard cutoff. This lets us combine it with dynamical grooming if desired.
      if (xz > fHardCutoff) {
        // If we pass any hard cutoff that we've set, then we will always find a splitting,
        // so we already increment the count.
        nPassedGrooming += 1;

        // Determine if we've found a new hardest splitting according to our hardness measure.
        double xHardnessMeasure = -0.05;
        switch (fGroomingMethod) {
          case kLeadingKt:
            xHardnessMeasure = xkt;
            break;
          case kDynamicalZ:
            xHardnessMeasure = CalculateDynamicalZ(j1, j2, jj, jetR);
            break;
          case kDynamicalKt:
            xHardnessMeasure = CalculateDynamicalKt(j1, j2, jj, jetR);
            break;
          case kDynamicalTime:
            xHardnessMeasure = CalculateDynamicalTime(j1, j2, jj, jetR);
            break;
          case kDynamicalCore:
            xHardnessMeasure = CalculateDynamicalCore(j1, j2, jj, jetR);
            break;
          default:
            AliFatalF("Unrecognized grooming method %s", GroomingMethodName().c_str());
            break;
        }
        // Update all of the properties if selected.
        if (xHardnessMeasure > hardnessMeasure) {
          hardnessMeasure = xHardnessMeasure;
          zg = xz;
          rg = xdeltaR;
          kt = xkt;
          nToSplit = nSplit;
          nGroomedToSplit = nPassedGrooming;
          // Store the subjets so we can do matching later.
          leadingSelectedSubjet = j1;
          subleadingSelectedSubjet = j2;
        }
      }

      // Follow the iterative (hardest) splitting.
      jj = j1;
    }

    // Store the selected subjets and their constituents to return them.
    if (kt > 0) {
      selectedSubjets = std::shared_ptr<SelectedSubjets>(new SelectedSubjets{leadingSelectedSubjet, leadingSelectedSubjet.constituents(),
                         subleadingSelectedSubjet, subleadingSelectedSubjet.constituents()});
    }

    // Cleanup the allocated cluster sequence.
    delete cs;
  } catch (const fastjet::Error&) {
    AliError(" [w] FJ Exception caught.");
  }

  // Store the extracted variables into the tree.
  std::string groomingMethod = GroomingMethodName();
  fSubstructureVariables[groomingMethod + "_" + prefix + "_kt"] = kt;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_z"] = zg;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_delta_R"] = rg;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_n_to_split"] = nToSplit;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_n_groomed_to_split"] = nGroomedToSplit;
  fSubstructureVariables[groomingMethod + "_" + prefix + "_n_passed_grooming"] = nPassedGrooming;

  return selectedSubjets;
}


void AliAnalysisTaskJetHardestKt::CheckSubjetResolution(AliEmcalJet* jet, AliEmcalJet* jetM)
{
  // Setup for jet
  std::vector<fastjet::PseudoJet> inputVectors;
  fastjet::PseudoJet pseudoTrack;

  for (int constituentIndex = 0; constituentIndex < jet->GetNumberOfTracks(); constituentIndex++) {
    AliVParticle* part = jet->Track(constituentIndex);
    if (!part) {
      continue;
    }
    pseudoTrack.reset(part->Px(), part->Py(), part->Pz(), part->E());
    pseudoTrack.set_user_index(jet->TrackAt(constituentIndex) + 100);
    inputVectors.push_back(pseudoTrack);
  }

  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
  fastjet::JetDefinition jetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30);

  // Setup for jetM
  std::vector<fastjet::PseudoJet> inputVectorsM;
  fastjet::PseudoJet pseudoTrackM;

  for (int constituentIndex = 0; constituentIndex < jetM->GetNumberOfTracks(); constituentIndex++) {
    AliVParticle* part = jetM->Track(constituentIndex);
    if (!part) {
      continue;
    }
    pseudoTrack.reset(part->Px(), part->Py(), part->Pz(), part->E());
    pseudoTrack.set_user_index(jetM->TrackAt(constituentIndex) + 100);
    inputVectorsM.push_back(pseudoTrack);
  }

  fastjet::JetAlgorithm jetalgoM(fastjet::cambridge_algorithm);
  fastjet::JetDefinition jetDefM(jetalgoM, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30);

  try {
    fastjet::ClusterSequence clustSeqSA(inputVectors, jetDef);
    std::vector<fastjet::PseudoJet> outputJets;
    outputJets = clustSeqSA.inclusive_jets(0);

    fastjet::ClusterSequence clustSeqSAM(inputVectorsM, jetDefM);
    std::vector<fastjet::PseudoJet> outputJetsM;
    outputJetsM = clustSeqSAM.inclusive_jets(0);

    fastjet::PseudoJet jj, jjM;
    fastjet::PseudoJet j1, j1M;
    fastjet::PseudoJet j2, j2M;
    jj = outputJets[0];
    jjM = outputJetsM[0];

    double z1 = 0;
    double z2 = 0;
    double zcut = 0.1;
    while ((jj.has_parents(j1, j2)) && (z1 < zcut)) {
      if (j1.perp() < j2.perp())
        swap(j1, j2);

      z1 = j2.perp() / (j1.perp() + j2.perp());
      jj = j1;
    }
    if (z1 < zcut)
      return;

    while ((jjM.has_parents(j1M, j2M)) && (z2 < zcut)) {
      if (j1M.perp() < j2M.perp())
        swap(j1M, j2M);

      z2 = j2M.perp() / (j1M.perp() + j2M.perp());
      jjM = j1M;
    }
    if (z2 < zcut)
      return;

    double delta_R1 = j1.delta_R(j1M);
    double delta_R2 = j2.delta_R(j2M);
    double delta_R = j1.delta_R(j2);
    double residz = (z1 - z2) / z2;
    double resid1 = (j1.perp() - j1M.perp()) / j1M.perp();
    double resid2 = (j2.perp() - j2M.perp()) / j2M.perp();

    if ((delta_R1 < fSubjetCutoff) && (delta_R2 < fSubjetCutoff)) {
      Double_t ResolEntries[5] = { outputJets[0].perp(), delta_R, resid1, resid2, residz };
      fHCheckResolutionSubjets->Fill(ResolEntries);
    }

  } catch (const fastjet::Error&) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

/**
 * Check if other tracks in a jet are close to a given particle.
 *
 * @param[in] jet Jet whose constituents are to be checked.
 * @param[in] part1 Particle to check other particles in relation to.
 */
bool AliAnalysisTaskJetHardestKt::CheckClosePartner(const AliEmcalJet* jet, const AliVParticle * part1)
{
  for (unsigned int i = 0; i < jet->GetNumberOfTracks(); i++) {
    AliVParticle * part2 = jet->Track(i);
    if (part2 == part1) {
      continue;
    }
    Double_t phi1 = part1->Phi();
    Double_t phi2 = part2->Phi();
    Double_t chg1 = part1->Charge();
    Double_t chg2 = part2->Charge();
    Double_t ptv1 = part1->Pt();
    Double_t ptv2 = part2->Pt();
    Double_t deta = part2->Eta() - part1->Eta();
    const Float_t kLimit = fPhiCutValue * 3;

    if (TMath::Abs(part1->Eta() - part2->Eta()) < fEtaCutValue * 2.5 * 3) {
      Float_t initdpsinner = (phi2 - TMath::ASin(0.075 * chg2 * fMagFieldPolarity * 0.8 / ptv2) -
                  (phi1 - TMath::ASin(0.075 * chg1 * fMagFieldPolarity * 0.8 / ptv1)));

      Float_t initdpsouter = (phi2 - TMath::ASin(0.075 * chg2 * fMagFieldPolarity * 2.5 / ptv2) -
                  (phi1 - TMath::ASin(0.075 * chg1 * fMagFieldPolarity * 2.5 / ptv1)));

      initdpsinner = TVector2::Phi_mpi_pi(initdpsinner);
      initdpsouter = TVector2::Phi_mpi_pi(initdpsouter);

      if (TMath::Abs(initdpsinner) < kLimit || TMath::Abs(initdpsouter) < kLimit ||
        initdpsinner * initdpsouter < 0) {
        Double_t mindps = 1e5;

        for (Double_t rad = 0.8; rad < 2.51; rad += 0.01) {
          Double_t dps = (phi2 - TMath::ASin(0.075 * chg2 * fMagFieldPolarity * rad / ptv2) -
                  (phi1 - TMath::ASin(0.075 * chg1 * fMagFieldPolarity * rad / ptv1)));
          dps = TVector2::Phi_mpi_pi(dps);
          if (TMath::Abs(dps) < TMath::Abs(mindps))
            mindps = dps;
        }
        if (TMath::Abs(mindps) < fPhiCutValue && TMath::Abs(deta) < fEtaCutValue)
          return kTRUE;
      }
    }
  }
  return kFALSE;
}

Bool_t AliAnalysisTaskJetHardestKt::RetrieveEventObjects()
{
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

/**
 * Add task to an existing analysis manager.
 *
 * @param[in] suffix Suffix string to attach to the task name
 */
AliAnalysisTaskJetHardestKt* AliAnalysisTaskJetHardestKt::AddTaskJetHardestKt(
 const char* njetsBase, const char* njetsUS, const char* njetsTrue, const char* njetsPartLevel, const Double_t R,
 const char* nrhoBase, const char* ntracks, const char* ntracksUS, const char* ntracksPartLevel, const char* nclusters,
 const char* ntracksTrue, const char* type, const char* CentEst, Int_t pSel,
 AliAnalysisTaskJetHardestKt::JetShapeType_t jetShapeType,
 AliAnalysisTaskJetHardestKt::JetShapeSub_t jetShapeSub,
 AliAnalysisTaskJetHardestKt::JetSelectionType_t jetSelection,
 Float_t minpTHTrigger, Float_t maxpTHTrigger, Float_t acut,
 AliAnalysisTaskJetHardestKt::DerivSubtrOrder_t derivSubtrOrder,
 const std::string& suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliErrorClass("No analysis manager to connect to.");
    return nullptr;
  }

  // Setup task name
  std::string taskName = "AliAnalysisTaskJetHardestKt";
  taskName += "_";
  taskName += njetsBase;
  std::string suffixName(suffix);
  if (suffixName != "") {
    taskName += "_";
    taskName += suffixName;
  }

  // Create task and configure as desired.
  AliAnalysisTaskJetHardestKt* task = new AliAnalysisTaskJetHardestKt(taskName.c_str());
  // Set a few general default.
  task->SetNCentBins(5);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);

  AliParticleContainer* trackCont = nullptr; // = task->AddTrackContainer(ntracks);

  if ((jetShapeSub == AliAnalysisTaskJetHardestKt::kConstSub ||
     jetShapeSub == AliAnalysisTaskJetHardestKt::kEventSub) &&
    ((jetShapeType == AliAnalysisTaskJetHardestKt::kData) ||
     (jetShapeType == AliAnalysisTaskJetHardestKt::kDetEmbPartPythia) ||
     (jetShapeType == AliAnalysisTaskJetHardestKt::kPythiaDef))) {
    trackCont = task->AddParticleContainer(ntracks);
  } else
    trackCont = task->AddTrackContainer(ntracks);

  // Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer* trackContUS = task->AddTrackContainer(ntracksUS);
  // Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer* trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  // Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
  if (jetShapeType == AliAnalysisTaskJetHardestKt::kDetEmbPartPythia)
    trackContTrue->SetIsEmbedding(true);
  AliParticleContainer* trackContPartLevel = 0;

  if ((jetShapeSub == AliAnalysisTaskJetHardestKt::kConstSub) &&
    ((jetShapeType == AliAnalysisTaskJetHardestKt::kMCTrue) ||
     (jetShapeType == AliAnalysisTaskJetHardestKt::kPythiaDef))) {
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  } else
    trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
  if (jetShapeType == AliAnalysisTaskJetHardestKt::kDetEmbPartPythia) {
    trackContPartLevel->SetIsEmbedding(true);
  }

  AliClusterContainer* clusterCont = task->AddClusterContainer(nclusters);

  AliJetContainer* jetContBase = nullptr;
  AliJetContainer* jetContUS = nullptr;
  AliJetContainer* jetContTrue = nullptr;
  AliJetContainer* jetContPart = nullptr;
  TString strType(type);

  if ((jetShapeType == AliAnalysisTaskJetHardestKt::kMCTrue ||
     (jetShapeType == AliAnalysisTaskJetHardestKt::kGenOnTheFly))) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackContPartLevel);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
  }

  if (jetShapeType == AliAnalysisTaskJetHardestKt::kData) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
      if (jetShapeSub == AliAnalysisTaskJetHardestKt::kConstSub)
        jetContBase->SetAreaEmcCut(-2);
    }
  }

  if (jetShapeType == AliAnalysisTaskJetHardestKt::kDetEmbPartPythia) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);

      if (jetShapeSub == AliAnalysisTaskJetHardestKt::kConstSub)
        jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue, strType, R);
    if (jetContTrue) {
      // This probably doesn't matter, but leaving just in case.
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
    }

    if (jetShapeSub == AliAnalysisTaskJetHardestKt::kConstSub ||
      jetShapeSub == AliAnalysisTaskJetHardestKt::kEventSub) {
      jetContUS = task->AddJetContainer(njetsUS, strType, R);
      if (jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(acut);
      }
    }

    jetContPart = task->AddJetContainer(njetsPartLevel, strType, R);
    if (jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(acut);
    }
  }

  if (jetShapeType == AliAnalysisTaskJetHardestKt::kPythiaDef) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }

    jetContTrue = task->AddJetContainer(njetsTrue, strType, R);
    if (jetContTrue) {
      // This probably doesn't matter, but leaving just in case.
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
    }

    if (jetShapeSub == AliAnalysisTaskJetHardestKt::kConstSub) {
      jetContUS = task->AddJetContainer(njetsUS, strType, R);
      if (jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(acut);
      }
    }

    jetContPart = task->AddJetContainer(njetsPartLevel, strType, R);
    if (jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(acut);
    }
  }

  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

  // Configuration is via YAML, so nothing else to be done here.
  mgr->AddTask(task);

  // Create containers for input/output
  std::string listOutputName = task->DetermineOutputContainerName(taskName);
  std::string treeOutputName = task->DetermineOutputContainerName(taskName + "Tree");

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

std::string AliAnalysisTaskJetHardestKt::DetermineOutputContainerName(std::string containerName) const
{
  if (fJetShapeType == AliAnalysisTaskJetHardestKt::kMCTrue)
    containerName += "_MCTrue";
  if (fJetShapeType == AliAnalysisTaskJetHardestKt::kData)
    containerName += "_Data";
  if (fJetShapeType == AliAnalysisTaskJetHardestKt::kPythiaDef)
    containerName += "_PythiaDef";
  if (fJetShapeSub == AliAnalysisTaskJetHardestKt::kNoSub)
    containerName += "_NoSub";
  if (fJetShapeSub == AliAnalysisTaskJetHardestKt::kConstSub)
    containerName += "_ConstSub";
  if (fJetShapeSub == AliAnalysisTaskJetHardestKt::kEventSub)
    containerName += "_EventSub";
  if (fJetShapeSub == AliAnalysisTaskJetHardestKt::kDerivSub)
    containerName += "_DerivSub";
  if (fJetSelection == AliAnalysisTaskJetHardestKt::kInclusive)
    containerName += "_Incl";

  return containerName;
}

/**
 * Prints information about the task.
 *
 * @return std::string containing information about the task.
 */
std::string AliAnalysisTaskJetHardestKt::toString() const
{
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  tempSS << "Dynmical grooming analysis task:\n";
  tempSS << "Jet properties:\n";
  tempSS << "\tShared momentum fraction: " << fMinFractionShared << "\n";
  tempSS << "\tJet shape type: " << GetKeyFromMapValue(fJetShapeType, fgkJetShapeTypeMap) << "\n";
  tempSS << "\tJet shape sub: " << GetKeyFromMapValue(fJetShapeSub, fgkJetShapeSubMap) << "\n";
  tempSS << "\tJet selection: " << GetKeyFromMapValue(fJetSelection, fgkJetSelectionMap) << "\n";
  tempSS << "\tPt threshold: " << fPtThreshold << " GeV/c\n";
  tempSS << "\tMax matching distance: " << fRMatching << "\n";
  tempSS << "General:\n";
  tempSS << "\tCentrality selection: ";
  if (fCentSelectOn) {
    tempSS << "(" << fCentMin << ", " << fCentMax << ")\n";
  }
  else {
    tempSS << "disabled.\n";
  }
  tempSS << "\tCheck subjet resolution: " << fCheckResolution << "\n";
  tempSS << "\tSubjet cutoff: " << fSubjetCutoff << "\n";
  tempSS << "\tConstituent pt cutoff: " << fMinPtConst << " GeV/c\n";
  tempSS << "\tGrooming method: " << GetKeyFromMapValue(fGroomingMethod, fgkGroomingMethodMap) << "\n";
  tempSS << "\tHard cutoff: " << fHardCutoff << "\n";
  tempSS << "\tGrooming method output name: " << GroomingMethodName() << "\n";
  tempSS << "\tConsider two track effects: " << fDoTwoTrack << "\n";
  tempSS << "\tCut double counting: " << fCutDoubleCounts << "\n";
  tempSS << "Two track effects settings:\n";
  tempSS << "\tPhi cut value: " << fPhiCutValue << "\n";
  tempSS << "\tEta cut value: " << fEtaCutValue << "\n";
  tempSS << "\tMagnetic field polarity: " << fMagFieldPolarity << "\n";
  tempSS << "Miscellaneous:\n";
  tempSS << "\tDerivative subtracter order: " << fDerivSubtrOrder << "\n";
  tempSS << "\tStore detector level jets: " << fStoreDetLevelJets << "\n";
  tempSS << "\tEnable subjet matching: " << fEnableSubjetMatching << "\n";
  // Substructure variables:
  tempSS << "Substructure variables:\n";
  for (const auto & p : fSubstructureVariables) {
    tempSS << "\t" << p.first << "\n";
  }
  // Jet containers
  tempSS << "Attached jet containers:\n";
  for (int i = 0; i < fJetCollArray.GetEntries(); i++) {
    tempSS << "\t" << GetJetContainer(i)->GetName() << "\n";
  }
  return tempSS.str();
}

/**
 * Print task information on an output stream using the string representation provided by
 * AliAnalysisTaskJetHardestKt::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream& AliAnalysisTaskJetHardestKt::Print(std::ostream& in) const
{
  in << toString();
  return in;
}

/**
 * Print task information using the string representation provided by
 * AliAnalysisTaskJetHardestKt::toString
 *
 * @param[in] opt Unused
 */
void AliAnalysisTaskJetHardestKt::Print(Option_t* opt) const { Printf("%s", toString().c_str()); }

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

/**
 * Hardest kt analysis task.
 */

/**
 * Implementation of the output stream operator for AliAnalysisTaskJetHardestKt. Printing
 * basic task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt& myTask)
{
  std::ostream& result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 *
 * NOTE: We don't swap the base class values because the base class doesn't implement swap.
 */
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt& first,
     PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt& second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fYAMLConfig, second.fYAMLConfig);
  swap(first.fConfigurationInitialized, second.fConfigurationInitialized);
  swap(first.fMinFractionShared, second.fMinFractionShared);
  swap(first.fJetShapeType, second.fJetShapeType);
  swap(first.fJetShapeSub, second.fJetShapeSub);
  swap(first.fJetSelection, second.fJetSelection);
  swap(first.fPtThreshold, second.fPtThreshold);
  swap(first.fRMatching, second.fRMatching);
  swap(first.fCentSelectOn, second.fCentSelectOn);
  swap(first.fCentMin, second.fCentMin);
  swap(first.fCentMax, second.fCentMax);
  swap(first.fCheckResolution, second.fCheckResolution);
  swap(first.fSubjetCutoff, second.fSubjetCutoff);
  swap(first.fMinPtConst, second.fMinPtConst);
  swap(first.fGroomingMethod, second.fGroomingMethod);
  swap(first.fHardCutoff, second.fHardCutoff);
  swap(first.fDoTwoTrack, second.fDoTwoTrack);
  swap(first.fCutDoubleCounts, second.fCutDoubleCounts);
  swap(first.fPhiCutValue, second.fPhiCutValue);
  swap(first.fEtaCutValue, second.fEtaCutValue);
  swap(first.fMagFieldPolarity, second.fMagFieldPolarity);
  swap(first.fDerivSubtrOrder, second.fDerivSubtrOrder);
  swap(first.fStoreDetLevelJets, second.fStoreDetLevelJets);
  swap(first.fSubstructureVariables, second.fSubstructureVariables);
  swap(first.fPtJet, second.fPtJet);
  swap(first.fHCheckResolutionSubjets, second.fHCheckResolutionSubjets);
  swap(first.fTreeSubstructure, second.fTreeSubstructure);
}
