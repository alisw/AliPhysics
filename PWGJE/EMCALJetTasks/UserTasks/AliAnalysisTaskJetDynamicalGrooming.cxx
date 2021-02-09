/**
 * Dynamical grooming jet substructure task. Adapted from AliAnalysisTaskJetDynamicalGrooming
 */

#include "AliAnalysisTaskJetDynamicalGrooming.h"

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

#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliEmcalPythiaInfo.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliRhoParameter.h"

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::SubstructureTree::Subjets);
/// \endcond

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::SubstructureTree::JetSplittings);
/// \endcond

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::SubstructureTree::JetConstituents);
/// \endcond

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::SubstructureTree::JetSubstructureSplittings);
/// \endcond

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming);
/// \endcond

namespace PWGJE
{
namespace EMCALJetTasks
{
namespace SubstructureTree
{

/**
 * Subjets
 */

/**
 * Default constructor
 */
Subjets::Subjets():
  fSplittingNodeIndex{},
  fPartOfIterativeSplitting{},
  fConstituentIndices{}
{
  // Nothing more to be done.
}

/**
 * Copy constructor
 */
Subjets::Subjets(const Subjets& other)
 : fSplittingNodeIndex{other.fSplittingNodeIndex},
  fPartOfIterativeSplitting{other.fPartOfIterativeSplitting},
  fConstituentIndices{other.fConstituentIndices}
{
  // Nothing more to be done.
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
Subjets& Subjets::operator=(Subjets other)
{
  swap(*this, other);
  return *this;
}

bool Subjets::Clear()
{
  fSplittingNodeIndex.clear();
  fPartOfIterativeSplitting.clear();
  fConstituentIndices.clear();
  return true;
}

std::tuple<unsigned short, bool, const std::vector<unsigned short>> Subjets::GetSubjet(int i) const
{
  return std::make_tuple(fSplittingNodeIndex.at(i), fPartOfIterativeSplitting.at(i), fConstituentIndices.at(i));
}

void Subjets::AddSubjet(const unsigned short splittingNodeIndex, const bool partOfIterativeSplitting, const std::vector<unsigned short> & constituentIndices)
{
  fSplittingNodeIndex.emplace_back(splittingNodeIndex);
  // NOTE: emplace_back isn't supported for std::vector<bool> until c++14.
  fPartOfIterativeSplitting.push_back(partOfIterativeSplitting);
  // Originally, we stored the constituent indices and their jagged indices separately to try to coax ROOT
  // into storing the nested vectors in a columnar format. However, even with that design, uproot can't
  // recreate the nested jagged array without a slow python loop. So we just store the indices directly
  // and wait for uproot 4. See: https://stackoverflow.com/q/60250877/12907985
  fConstituentIndices.emplace_back(constituentIndices);
}

/**
 * Prints information about the task.
 *
 * @return std::string containing information about the task.
 */
std::string Subjets::toString() const
{
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  tempSS << "Subjets:\n";
  for (std::size_t i = 0; i < fSplittingNodeIndex.size(); i++)
  {
    tempSS << "#" << (i + 1) << ": Splitting Node: " << fSplittingNodeIndex.at(i)
        << ", part of iterative splitting = " << fPartOfIterativeSplitting.at(i)
        << ", number of jet constituents = " << fConstituentIndices.at(i).size() << "\n";
  }
  return tempSS.str();
}

/**
 * Print task information on an output stream using the string representation provided by
 * Subjets::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream& Subjets::Print(std::ostream& in) const
{
  in << toString();
  return in;
}

/**
 * Print task information using the string representation provided by
 * Subjets::toString
 *
 * @param[in] opt Unused
 */
void Subjets::Print(Option_t* opt) const { Printf("%s", toString().c_str()); }

/**
 * Jet splittings
 */

/**
 * Default constructor.
 */
JetSplittings::JetSplittings():
  fKt{},
  fDeltaR{},
  fZ{},
  fParentIndex{}
{
  // Nothing more to be done.
}

/**
 * Copy constructor
 */
JetSplittings::JetSplittings(const JetSplittings& other)
 : fKt{other.fKt},
  fDeltaR{other.fDeltaR},
  fZ{other.fZ},
  fParentIndex{other.fParentIndex}
{
  // Nothing more to be done.
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
JetSplittings& JetSplittings::operator=(JetSplittings other)
{
  swap(*this, other);
  return *this;
}

bool JetSplittings::Clear()
{
  fKt.clear();
  fDeltaR.clear();
  fZ.clear();
  fParentIndex.clear();
  return true;
}

void JetSplittings::AddSplitting(float kt, float deltaR, float z, short i)
{
  fKt.emplace_back(kt);
  fDeltaR.emplace_back(deltaR);
  fZ.emplace_back(z);
  fParentIndex.emplace_back(i);
}

std::tuple<float, float, float, short> JetSplittings::GetSplitting(int i) const
{
  return std::make_tuple(fKt.at(i), fDeltaR.at(i), fZ.at(i), fParentIndex.at(i));
}

/**
 * Prints information about the task.
 *
 * @return std::string containing information about the task.
 */
std::string JetSplittings::toString() const
{
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  tempSS << "Jet splittings:\n";
  for (std::size_t i = 0; i < fKt.size(); i++)
  {
    tempSS << "#" << (i + 1) << ": kT = " << fKt.at(i)
        << ", deltaR = " << fDeltaR.at(i) << ", z = " << fZ.at(i)
        << ", parent = " << fParentIndex.at(i) << "\n";
  }
  return tempSS.str();
}

/**
 * Print task information on an output stream using the string representation provided by
 * JetSplittings::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream& JetSplittings::Print(std::ostream& in) const
{
  in << toString();
  return in;
}

/**
 * Print task information using the string representation provided by
 * JetSplittings::toString
 *
 * @param[in] opt Unused
 */
void JetSplittings::Print(Option_t* opt) const { Printf("%s", toString().c_str()); }

/**
 * Jet constituents.
 */
// Constant shift to be applied to storing the global index. This ensures that the value is
// large enough that we're never going to overlap with the MCLabel. It also needs to be large
// enough such that we'll never overlap with values from the index map. 2000000 allows for 20
// collections, which should be more than enough.
const int JetConstituents::fgkGlobalIndexOffset = 2000000;

/**
 * Default constructor.
 */
JetConstituents::JetConstituents():
  fPt{},
  fEta{},
  fPhi{},
  fID{}
{
  // Nothing more to be done.
}

/**
 * Copy constructor
 */
JetConstituents::JetConstituents(const JetConstituents& other)
 : fPt{other.fPt},
  fEta{other.fEta},
  fPhi{other.fPhi},
  fID{other.fID}
{
  // Nothing more to be done.
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
JetConstituents& JetConstituents::operator=(JetConstituents other)
{
  swap(*this, other);
  return *this;
}

bool JetConstituents::Clear()
{
  fPt.clear();
  fEta.clear();
  fPhi.clear();
  fID.clear();
  return true;
}

void JetConstituents::AddJetConstituent(const AliVParticle * part, const int & id)
{
  fPt.emplace_back(part->Pt());
  fEta.emplace_back(part->Eta());
  fPhi.emplace_back(part->Phi());
  // NOTE: We don't use the user_index() here because we need to use that for indexing the
  //       constituents that are contained in each subjet.
  fID.emplace_back(id);
}

std::tuple<float, float, float, int> JetConstituents::GetJetConstituent(int i) const
{
  return std::make_tuple(fPt.at(i), fEta.at(i), fPhi.at(i), fID.at(i));
}

/**
 * Prints information about the task.
 *
 * @return std::string containing information about the task.
 */
std::string JetConstituents::toString() const
{
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  tempSS << "Jet constituents:\n";
  for (std::size_t i = 0; i < fPt.size(); i++)
  {
    tempSS << "#" << (i + 1) << ": pt = " << fPt.at(i)
        << ", eta = " << fEta.at(i) << ", phi = " << fPhi.at(i)
        << ", ID = " << fID.at(i) << "\n";
  }
  return tempSS.str();
}

/**
 * Print task information on an output stream using the string representation provided by
 * JetConstituents::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream& JetConstituents::Print(std::ostream& in) const
{
  in << toString();
  return in;
}

/**
 * Print task information using the string representation provided by
 * JetConstituents::toString
 *
 * @param[in] opt Unused
 */
void JetConstituents::Print(Option_t* opt) const { Printf("%s", toString().c_str()); }

/**
 * Jet substructure splittings container.
 */

/**
 * Default constructor.
 */
JetSubstructureSplittings::JetSubstructureSplittings():
  fJetPt{0},
  fJetConstituents{},
  fJetSplittings{},
  fSubjets{}
{
  // Nothing more to be done.
}

/**
 * Copy constructor
 */
JetSubstructureSplittings::JetSubstructureSplittings(
 const JetSubstructureSplittings& other)
 : fJetPt{other.fJetPt},
  fJetConstituents{other.fJetConstituents},
  fJetSplittings{other.fJetSplittings},
  fSubjets{other.fSubjets}
{
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
JetSubstructureSplittings& JetSubstructureSplittings::operator=(
 JetSubstructureSplittings other)
{
  swap(*this, other);
  return *this;
}

bool JetSubstructureSplittings::Clear()
{
  fJetPt = 0;
  fJetConstituents.Clear();
  fJetSplittings.Clear();
  fSubjets.Clear();
  return true;
}

/**
 * Add a jet constituent to the object.
 *
 * @param[in] part Constituent to be added.
 */
void JetSubstructureSplittings::AddJetConstituent(const AliVParticle * part, const int & id)
{
  fJetConstituents.AddJetConstituent(part, id);
}

/**
 * Add a jet splitting to the object.
 *
 * @param[in] kt Kt of the splitting.
 * @param[in] deltaR Delta R between the subjets.
 * @param[in] z Momentum sharing between the subjets.
 */
void JetSubstructureSplittings::AddSplitting(float kt, float deltaR, float z, short parentIndex)
{
  fJetSplittings.AddSplitting(kt, deltaR, z, parentIndex);
}

/**
 * Add a subjet to the object.
 *
 * @param[in] part Constituent to be added.
 */
void JetSubstructureSplittings::AddSubjet(const unsigned short splittingNodeIndex, const bool partOfIterativeSplitting,
                     const std::vector<unsigned short>& constituentIndices)
{
  return fSubjets.AddSubjet(splittingNodeIndex, partOfIterativeSplitting, constituentIndices);
}

std::tuple<float, float, float, int> JetSubstructureSplittings::GetJetConstituent(int i) const
{
  return fJetConstituents.GetJetConstituent(i);
}

std::tuple<float, float, float, short> JetSubstructureSplittings::GetSplitting(int i) const
{
  return fJetSplittings.GetSplitting(i);
}

std::tuple<unsigned short, bool, const std::vector<unsigned short>> JetSubstructureSplittings::GetSubjet(int i) const
{
  return fSubjets.GetSubjet(i);
}

/**
 * Prints information about the task.
 *
 * @return std::string containing information about the task.
 */
std::string JetSubstructureSplittings::toString() const
{
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  tempSS << "Splitting information: ";
  tempSS << "Jet pt = " << fJetPt << "\n";
  tempSS << fSubjets;
  tempSS << fJetSplittings;
  tempSS << fJetConstituents;
  return tempSS.str();
}

/**
 * Print task information on an output stream using the string representation provided by
 * JetSubstructureSplittings::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream& JetSubstructureSplittings::Print(std::ostream& in) const
{
  in << toString();
  return in;
}

/**
 * Print task information using the string representation provided by
 * JetSubstructureSplittings::toString
 *
 * @param[in] opt Unused
 */
void JetSubstructureSplittings::Print(Option_t* opt) const { Printf("%s", toString().c_str()); }

} /* namespace SubstructureTree */

/**
 * Dynamical grooming analysis task.
 */
const std::map<std::string, AliAnalysisTaskJetDynamicalGrooming::JetShapeType_t> AliAnalysisTaskJetDynamicalGrooming::fgkJetShapeTypeMap = {
  { "kMCTrue", AliAnalysisTaskJetDynamicalGrooming::kMCTrue },
  { "kTrueDet", AliAnalysisTaskJetDynamicalGrooming::kTrueDet },
  { "kData", AliAnalysisTaskJetDynamicalGrooming::kData },
  { "kDetEmb", AliAnalysisTaskJetDynamicalGrooming::kDetEmb },
  { "kDetEmbPart", AliAnalysisTaskJetDynamicalGrooming::kDetEmbPart },
  { "kPythiaDef", AliAnalysisTaskJetDynamicalGrooming::kPythiaDef },
  { "kDetEmbPartPythia", AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia },
  { "kGenOnTheFly", AliAnalysisTaskJetDynamicalGrooming::kGenOnTheFly }
};
const std::map<std::string, AliAnalysisTaskJetDynamicalGrooming::JetShapeSub_t> AliAnalysisTaskJetDynamicalGrooming::fgkJetShapeSubMap = {
  { "kNoSub", AliAnalysisTaskJetDynamicalGrooming::kNoSub },
  { "kConstSub", AliAnalysisTaskJetDynamicalGrooming::kConstSub },
  { "kDerivSub", AliAnalysisTaskJetDynamicalGrooming::kDerivSub },
  { "kEventSub", AliAnalysisTaskJetDynamicalGrooming::kEventSub }
};
const std::map<std::string, AliAnalysisTaskJetDynamicalGrooming::JetSelectionType_t> AliAnalysisTaskJetDynamicalGrooming::fgkJetSelectionMap = {
  { "kInclusive", AliAnalysisTaskJetDynamicalGrooming::kInclusive },
  { "kRecoil", AliAnalysisTaskJetDynamicalGrooming::kRecoil }
};
const std::map<std::string, AliAnalysisTaskJetDynamicalGrooming::DerivSubtrOrder_t> AliAnalysisTaskJetDynamicalGrooming::fgkDerivSubtrOrderMap = {
  { "kSecondOrder", AliAnalysisTaskJetDynamicalGrooming::kSecondOrder },
  { "kFirstOrder", AliAnalysisTaskJetDynamicalGrooming::kFirstOrder }
};

/**
 * Default constructor.
 */
AliAnalysisTaskJetDynamicalGrooming::AliAnalysisTaskJetDynamicalGrooming()
 : AliAnalysisTaskEmcalJet("AliAnalysisTaskJetDynamicalGrooming", kTRUE),
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
  fHardCutoff(0),
  fDoTwoTrack(kFALSE),
  fCutDoubleCounts(kTRUE),
  fPhiCutValue(0.02),
  fEtaCutValue(0.02),
  fMagFieldPolarity(1),
  fDerivSubtrOrder(0),
  fStoreDetLevelJets(kFALSE),
  fStoreRecursiveSplittings(false),
  fDataJetSplittings(),
  fMatchedJetSplittings(),
  fDetLevelJetSplittings(),
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
AliAnalysisTaskJetDynamicalGrooming::AliAnalysisTaskJetDynamicalGrooming(const char* name)
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
  fHardCutoff(0),
  fDoTwoTrack(kFALSE),
  fCutDoubleCounts(kTRUE),
  fPhiCutValue(0.02),
  fEtaCutValue(0.02),
  fMagFieldPolarity(1),
  fDerivSubtrOrder(0),
  fStoreDetLevelJets(kFALSE),
  fStoreRecursiveSplittings(false),
  fDataJetSplittings(),
  fMatchedJetSplittings(),
  fDetLevelJetSplittings(),
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
AliAnalysisTaskJetDynamicalGrooming::AliAnalysisTaskJetDynamicalGrooming(
 const AliAnalysisTaskJetDynamicalGrooming& other)
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
  fHardCutoff(other.fHardCutoff),
  fDoTwoTrack(other.fDoTwoTrack),
  fCutDoubleCounts(other.fCutDoubleCounts),
  fPhiCutValue(other.fPhiCutValue),
  fEtaCutValue(other.fEtaCutValue),
  fMagFieldPolarity(other.fMagFieldPolarity),
  fDerivSubtrOrder(other.fDerivSubtrOrder),
  fStoreDetLevelJets(other.fStoreDetLevelJets),
  fStoreRecursiveSplittings(other.fStoreRecursiveSplittings),
  fDataJetSplittings(),
  fMatchedJetSplittings(),
  fDetLevelJetSplittings(),
  fPtJet(nullptr),
  fHCheckResolutionSubjets(nullptr),
  fTreeSubstructure(nullptr)
{
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
AliAnalysisTaskJetDynamicalGrooming& AliAnalysisTaskJetDynamicalGrooming::operator=(
 AliAnalysisTaskJetDynamicalGrooming other)
{
  swap(*this, other);
  return *this;
}

/**
 * Retrieve task properties from the YAML configuration.
 */
void AliAnalysisTaskJetDynamicalGrooming::RetrieveAndSetTaskPropertiesFromYAMLConfig()
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
  fYAMLConfig.GetProperty({baseName, "storeRecursiveSplittings"}, fStoreRecursiveSplittings, false);
}

/**
 * Initialize task.
 */
bool AliAnalysisTaskJetDynamicalGrooming::Initialize()
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

/**
 * Setup output tree.
 *
 * Branches are only enabled if the output is actually needed.
 *
 * NOTE: This tree uses a number of more subtle ROOT tree functions. It relies heavily on information from
 * the ROOT users guide: https://root.cern.ch/root/htmldoc/guides/users-guide/Trees.html
 */
void AliAnalysisTaskJetDynamicalGrooming::SetupTree()
{
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeSubstructure = new TTree(nameoutput, nameoutput);
  // Ensure that splitting the object members into their own branches. It doesn't hurt to pass a higher
  // number to a object that doesn't need to split, so we pass the same to everything.
  int splitLevel = 4;
  // Buffer size for storing the branch until it is saved to file.
  // To quote ROOT, "A recommended buffer size is 32000 bytes if you have less than 50 branches."
  int bufferSize = 32000;
  // Including the "." ensures that the identically named members in the splitting functions are stored
  // with unique names. This will simplify analysis later.
  fTreeSubstructure->Branch("data.", &fDataJetSplittings, bufferSize, splitLevel);

  if (fJetShapeType == kDetEmbPartPythia || fJetShapeType == kPythiaDef) {
    fTreeSubstructure->Branch("matched.", &fMatchedJetSplittings, bufferSize, splitLevel);
    if (fStoreDetLevelJets) {
      fTreeSubstructure->Branch("detLevel.", &fDetLevelJetSplittings, bufferSize, splitLevel);
    }

    // Also add the unsubtracted leading kt (which is what we actually want to use for the double counting cut)
    if (fJetShapeType == kDetEmbPartPythia) {
      fTreeSubstructure->Branch("data_leading_track_pt", &fDataLeadingTrackPtUnsub);
    }
  }

  if (fIsPythia) {
    // Will be automatically filled by AliAnalysisTaskEmcal.
    fTreeSubstructure->Branch("ptHardBin", fPtHardInitialized ? &fPtHardBinGlobal : &fPtHardBin);
    fTreeSubstructure->Branch("ptHard", &fPtHard);
  }
}

/**
 * Create output objects.
 */
void AliAnalysisTaskJetDynamicalGrooming::UserCreateOutputObjects()
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

Bool_t AliAnalysisTaskJetDynamicalGrooming::Run()
{
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().

  return kTRUE;
}

double AliAnalysisTaskJetDynamicalGrooming::DynamicalGrooming(const fastjet::PseudoJet& subjet1,
                               const fastjet::PseudoJet& subjet2,
                               const fastjet::PseudoJet& parent, const double R,
                               const double a) const
{
  double deltaR = subjet1.delta_R(subjet2);
  double z = subjet1.pt() / parent.pt();
  return z * (1 - z) * parent.pt() * std::pow(deltaR / R, a);
}

double AliAnalysisTaskJetDynamicalGrooming::CalculateZDrop(const fastjet::PseudoJet& subjet1,
                              const fastjet::PseudoJet& subjet2,
                              const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 0.1);
}

double AliAnalysisTaskJetDynamicalGrooming::CalculateKtDrop(const fastjet::PseudoJet& subjet1,
                              const fastjet::PseudoJet& subjet2,
                              const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 1);
}

double AliAnalysisTaskJetDynamicalGrooming::CalculateTimeDrop(const fastjet::PseudoJet& subjet1,
                               const fastjet::PseudoJet& subjet2,
                               const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 2);
}

Bool_t AliAnalysisTaskJetDynamicalGrooming::FillHistograms()
{
  AliJetContainer* jetCont = GetJetContainer(0);
  // container zero is always the base container: the data container, the
  // embedded subtracted in the case of embedding or the detector level in case
  // of pythia

  if (fCentSelectOn) {
    if ((fCent > fCentMax) || (fCent < fCentMin)) {
      return 0;
    }
  }

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

      // Clear out previously stored values
      fDataJetSplittings.Clear();
      fDataLeadingTrackPtUnsub = 0;
      fMatchedJetSplittings.Clear();
      fDetLevelJetSplittings.Clear();

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

      fDataJetSplittings.SetJetPt(ptSubtracted);

      // The double counting cut should be applied to the unsubtracted hybrid max track pt. Confusingly,
      // MaxTrackPt() called on the constituent subtracted jet returns the _unsubtracted_ jet pt. Since we store
      // all jet constituents, we can already find the max subtracted constituent pt offline. However, we also need
      // the unsubtracted max track pt. So we store it explicitly alongside the rest of the information.
      // That way, we'll be able to apply the double counting cut offline, regardless of the whether we want to
      // use the subtracted or unsubtracted case.
      //
      // We also maintain the ability to apply the double counting cut during analysis.
      if (fJetShapeType == kDetEmbPartPythia || fJetShapeType == kData) {
        fDataLeadingTrackPtUnsub = jet1->MaxTrackPt();
      }
      if (fJetShapeType == kDetEmbPartPythia) {
        // Cut if the leading hybrid track pt is greater than the leading detector level track pt.
        if (fCutDoubleCounts == kTRUE && (jet1->MaxTrackPt() > jet2->MaxTrackPt())) {
          continue;
        }
      }

      // Determine the splittings of the main (data) jet.
      IterativeParents(jet1, fDataJetSplittings, true);

      // If appropriate, then fill the matched and/or detector level jets.
      if (fJetShapeType == kPythiaDef) {
        /*kMatched = 1;
        if (fJetShapeSub == kConstSub)
          kMatched = 3;*/

        fMatchedJetSplittings.SetJetPt(jet3->Pt());
        IterativeParents(jet3, fMatchedJetSplittings, false);
      }

      if (fJetShapeType == kDetEmbPartPythia) {
        /*if (fJetShapeSub == kConstSub)
          kMatched = 3;
        if (fJetShapeSub == kDerivSub)
          kMatched = 2;*/
        fMatchedJetSplittings.SetJetPt(jet3->Pt());
        IterativeParents(jet3, fMatchedJetSplittings, false);
        if (fStoreDetLevelJets) {
          fDetLevelJetSplittings.SetJetPt(jet2->Pt());
          IterativeParents(jet2, fDetLevelJetSplittings, false);
        }
      }

      fTreeSubstructure->Fill();
    }
  }

  return kTRUE;
}

/**
 * Determine and iterate through jet splittings for a given jet. The output is stored in the given
 * jet substructure output container.
 *
 * @param[in] jet Jet to be declustered.
 * @param[in, out] jetSplittings Jet substructure output object which will be used to store the splittings.
 * @param[in] isData If True, treat the splitting as coming from data. This means that ghosts are utilized and track resolution may be considered.
 */
void AliAnalysisTaskJetDynamicalGrooming::IterativeParents(AliEmcalJet* jet,
                              SubstructureTree::JetSubstructureSplittings & jetSplittings,
                              bool isData)
{
  AliDebugStream(1) << "Beginning iteration through the splittings.\n";
  std::vector<fastjet::PseudoJet> inputVectors;
  fastjet::PseudoJet pseudoTrack;
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
    // NOTE: This must be the constituent index to allow the subjets to properly determine which constituents are included
    //       in each subjet.
    pseudoTrack.set_user_index(constituentIndex);
    inputVectors.push_back(pseudoTrack);

    // Also store the jet constituents in the output
    // Ensure they are uniquely identified (per type of jet) using the id. We don't use the global index (accessed via
    // `TrackAt(int)`) because they won't be the same for the subtracted and unsubtracted constituents (they are different
    // TClonesArrays). Instead, we label the part and det level with the MClabel (ie. for those which have it available),
    // and for the data (where it always equals -1), we use the global index (offset sufficiently) so it won't overlap with
    // other values.
    int id = part->GetLabel() != -1 ? part->GetLabel() : (jet->TrackAt(constituentIndex) + SubstructureTree::JetConstituents::GetGlobalIndexOffset());
    jetSplittings.AddJetConstituent(part , id);
  }

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

    // Store the jet splittings.
    int splittingNodeIndex = -1;
    ExtractJetSplittings(jetSplittings, jj, splittingNodeIndex, true);

    // Cleanup the allocated cluster sequence.
    delete cs;
  } catch (const fastjet::Error&) {
    AliError(" [w] FJ Exception caught.");
  }

  return;
}

void AliAnalysisTaskJetDynamicalGrooming::ExtractJetSplittings(SubstructureTree::JetSubstructureSplittings & jetSplittings, fastjet::PseudoJet & inputJet, int splittingNodeIndex, bool followingIterativeSplitting)
{
  fastjet::PseudoJet j1;
  fastjet::PseudoJet j2;
  if (inputJet.has_parents(j1, j2) == false) {
    // No parents, so we're done - just return.
    return;
  }

  // j1 should always be the harder of the two subjets.
  if (j1.perp() < j2.perp()) {
    swap(j1, j2);
  }

  // We have a splitting. Record the properties.
  double z = j2.perp() / (j2.perp() + j1.perp());
  double delta_R = j1.delta_R(j2);
  double xkt = j2.perp() * sin(delta_R);
  // Add the splitting node.
  jetSplittings.AddSplitting(xkt, delta_R, z, splittingNodeIndex);
  // Determine which splitting parent the subjets will point to (ie. the one that
  // we just stored). It's stored at the end of the splittings array. (which we offset
  // by -1 to stay within the array).
  splittingNodeIndex = jetSplittings.GetNumberOfSplittings() - 1;
  // Store the subjets
  std::vector<unsigned short> j1ConstituentIndices, j2ConstituentIndices;
  for (auto constituent: j1.constituents()) {
    j1ConstituentIndices.emplace_back(constituent.user_index());
  }
  for (auto constituent: j2.constituents()) {
    j2ConstituentIndices.emplace_back(constituent.user_index());
  }
  jetSplittings.AddSubjet(splittingNodeIndex, followingIterativeSplitting, j1ConstituentIndices);
  jetSplittings.AddSubjet(splittingNodeIndex, false, j2ConstituentIndices);

  // Recurse as necessary to get the rest of the splittings.
  ExtractJetSplittings(jetSplittings, j1, splittingNodeIndex, followingIterativeSplitting);
  if (fStoreRecursiveSplittings == true) {
    ExtractJetSplittings(jetSplittings, j2, splittingNodeIndex, false);
  }
}

void AliAnalysisTaskJetDynamicalGrooming::CheckSubjetResolution(AliEmcalJet* jet, AliEmcalJet* jetM)
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
bool AliAnalysisTaskJetDynamicalGrooming::CheckClosePartner(const AliEmcalJet* jet, const AliVParticle * part1)
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

Bool_t AliAnalysisTaskJetDynamicalGrooming::RetrieveEventObjects()
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
AliAnalysisTaskJetDynamicalGrooming* AliAnalysisTaskJetDynamicalGrooming::AddTaskJetDynamicalGrooming(
 const char* njetsBase, const char* njetsUS, const char* njetsTrue, const char* njetsPartLevel, const Double_t R,
 const char* nrhoBase, const char* ntracks, const char* ntracksUS, const char* ntracksPartLevel, const char* nclusters,
 const char* ntracksTrue, const char* type, const char* CentEst, Int_t pSel,
 AliAnalysisTaskJetDynamicalGrooming::JetShapeType_t jetShapeType,
 AliAnalysisTaskJetDynamicalGrooming::JetShapeSub_t jetShapeSub,
 AliAnalysisTaskJetDynamicalGrooming::JetSelectionType_t jetSelection,
 Float_t minpTHTrigger, Float_t maxpTHTrigger, Float_t acut,
 AliAnalysisTaskJetDynamicalGrooming::DerivSubtrOrder_t derivSubtrOrder,
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
  std::string taskName = "AliAnalysisTaskJetDynamicalGrooming";
  taskName += "_";
  taskName += njetsBase;
  std::string suffixName(suffix);
  if (suffixName != "") {
    taskName += "_";
    taskName += suffixName;
  }

  // Create task and configure as desired.
  AliAnalysisTaskJetDynamicalGrooming* task = new AliAnalysisTaskJetDynamicalGrooming(taskName.c_str());
  // Set a few general default.
  task->SetNCentBins(5);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);

  AliParticleContainer* trackCont = nullptr; // = task->AddTrackContainer(ntracks);

  if ((jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub ||
     jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kEventSub) &&
    ((jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kData) ||
     (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia) ||
     (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kPythiaDef))) {
    trackCont = task->AddParticleContainer(ntracks);
  } else
    trackCont = task->AddTrackContainer(ntracks);

  // Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer* trackContUS = task->AddTrackContainer(ntracksUS);
  // Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer* trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  // Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia)
    trackContTrue->SetIsEmbedding(true);
  AliParticleContainer* trackContPartLevel = 0;

  if ((jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub) &&
    ((jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kMCTrue) ||
     (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kPythiaDef))) {
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  } else
    trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia)
    trackContPartLevel->SetIsEmbedding(true);
  // Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);

  AliClusterContainer* clusterCont = task->AddClusterContainer(nclusters);

  AliJetContainer* jetContBase = nullptr;
  AliJetContainer* jetContUS = nullptr;
  AliJetContainer* jetContTrue = nullptr;
  AliJetContainer* jetContPart = nullptr;
  TString strType(type);

  if ((jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kMCTrue ||
     (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kGenOnTheFly))) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackContPartLevel);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
  }

  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kData) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
      if (jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub)
        jetContBase->SetAreaEmcCut(-2);
    }
  }

  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);

      if (jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub)
        jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue, strType, R);
    if (jetContTrue) {
      // This probably doesn't matter, but leaving just in case.
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
    }

    if (jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub ||
      jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kEventSub) {
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

  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kPythiaDef) {
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

    if (jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub) {
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

std::string AliAnalysisTaskJetDynamicalGrooming::DetermineOutputContainerName(std::string containerName) const
{
  if (fJetShapeType == AliAnalysisTaskJetDynamicalGrooming::kMCTrue)
    containerName += "_MCTrue";
  if (fJetShapeType == AliAnalysisTaskJetDynamicalGrooming::kData)
    containerName += "_Data";
  if (fJetShapeType == AliAnalysisTaskJetDynamicalGrooming::kPythiaDef)
    containerName += "_PythiaDef";
  if (fJetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kNoSub)
    containerName += "_NoSub";
  if (fJetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub)
    containerName += "_ConstSub";
  if (fJetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kEventSub)
    containerName += "_EventSub";
  if (fJetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kDerivSub)
    containerName += "_DerivSub";
  if (fJetSelection == AliAnalysisTaskJetDynamicalGrooming::kInclusive)
    containerName += "_Incl";

  return containerName;
}

/**
 * Prints information about the task.
 *
 * @return std::string containing information about the task.
 */
std::string AliAnalysisTaskJetDynamicalGrooming::toString() const
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
  tempSS << "\tHard cutoff: " << fHardCutoff << "\n";
  tempSS << "\tConsider two track effects: " << fDoTwoTrack << "\n";
  tempSS << "\tCut double counting: " << fCutDoubleCounts << "\n";
  tempSS << "Two track effects settings:\n";
  tempSS << "\tPhi cut value: " << fPhiCutValue << "\n";
  tempSS << "\tEta cut value: " << fEtaCutValue << "\n";
  tempSS << "\tMagnetic field polarity: " << fMagFieldPolarity << "\n";
  tempSS << "Miscellaneous:\n";
  tempSS << "\tDerivative subtracter order: " << fDerivSubtrOrder << "\n";
  tempSS << "\tStore detector level jets: " << fStoreDetLevelJets << "\n";
  tempSS << "\tStore recursive jet splittings (instead of just iterative): " << fStoreRecursiveSplittings << "\n";
  // Jet containers
  tempSS << "Attached jet containers:\n";
  for (int i = 0; i < fJetCollArray.GetEntries(); i++)
  {
    tempSS << "\t" << GetJetContainer(i)->GetName() << "\n";
  }
  return tempSS.str();
}

/**
 * Print task information on an output stream using the string representation provided by
 * AliAnalysisTaskJetDynamicalGrooming::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream& AliAnalysisTaskJetDynamicalGrooming::Print(std::ostream& in) const
{
  in << toString();
  return in;
}

/**
 * Print task information using the string representation provided by
 * AliAnalysisTaskJetDynamicalGrooming::toString
 *
 * @param[in] opt Unused
 */
void AliAnalysisTaskJetDynamicalGrooming::Print(Option_t* opt) const { Printf("%s", toString().c_str()); }

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

/**
 * Subjets
 */

/**
 * Implementation of the output stream operator for SubstructureTree::Subjets. Printing
 * basic task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::SubstructureTree::Subjets& myTask)
{
  std::ostream& result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 */
void swap(PWGJE::EMCALJetTasks::SubstructureTree::Subjets& first,
     PWGJE::EMCALJetTasks::SubstructureTree::Subjets& second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fSplittingNodeIndex, second.fSplittingNodeIndex);
  swap(first.fPartOfIterativeSplitting, second.fPartOfIterativeSplitting);
  swap(first.fConstituentIndices, second.fConstituentIndices);
}

/**
 * JetSplittings
 */

/**
 * Implementation of the output stream operator for SubstructureTree::JetSplittings. Printing
 * basic task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::SubstructureTree::JetSplittings& myTask)
{
  std::ostream& result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 */
void swap(PWGJE::EMCALJetTasks::SubstructureTree::JetSplittings& first,
     PWGJE::EMCALJetTasks::SubstructureTree::JetSplittings& second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fKt, second.fKt);
  swap(first.fDeltaR, second.fDeltaR);
  swap(first.fZ, second.fZ);
  swap(first.fParentIndex, second.fParentIndex);
}

/**
 * JetConstituents
 */

/**
 * Implementation of the output stream operator for SubstructureTree::JetConstituents. Printing
 * basic task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::SubstructureTree::JetConstituents& myTask)
{
  std::ostream& result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 */
void swap(PWGJE::EMCALJetTasks::SubstructureTree::JetConstituents& first,
     PWGJE::EMCALJetTasks::SubstructureTree::JetConstituents& second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fPt, second.fPt);
  swap(first.fEta, second.fEta);
  swap(first.fPhi, second.fPhi);
  swap(first.fID, second.fID);
}

/**
 * Jet substructure splittings
 */

/**
 * Implementation of the output stream operator for JetSubstructureSplittings. Printing
 * basic task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::SubstructureTree::JetSubstructureSplittings& myTask)
{
  std::ostream& result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 */
void swap(PWGJE::EMCALJetTasks::SubstructureTree::JetSubstructureSplittings& first,
     PWGJE::EMCALJetTasks::SubstructureTree::JetSubstructureSplittings& second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fJetPt, second.fJetPt);
  swap(first.fJetConstituents, second.fJetConstituents);
  swap(first.fJetSplittings, second.fJetSplittings);
  swap(first.fSubjets, second.fSubjets);
}

/**
 * Dynamical grooming analysis task.
 */

/**
 * Implementation of the output stream operator for AliAnalysisTaskJetDynamicalGrooming. Printing
 * basic task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& myTask)
{
  std::ostream& result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 *
 * NOTE: We don't swap the base class values because the base class doesn't implement swap.
 */
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& first,
     PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& second)
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
  swap(first.fHardCutoff, second.fHardCutoff);
  swap(first.fDoTwoTrack, second.fDoTwoTrack);
  swap(first.fCutDoubleCounts, second.fCutDoubleCounts);
  swap(first.fPhiCutValue, second.fPhiCutValue);
  swap(first.fEtaCutValue, second.fEtaCutValue);
  swap(first.fMagFieldPolarity, second.fMagFieldPolarity);
  swap(first.fDerivSubtrOrder, second.fDerivSubtrOrder);
  swap(first.fStoreDetLevelJets, second.fStoreDetLevelJets);
  swap(first.fStoreRecursiveSplittings, second.fStoreRecursiveSplittings);
  swap(first.fDataJetSplittings, second.fDataJetSplittings);
  swap(first.fMatchedJetSplittings, second.fMatchedJetSplittings);
  swap(first.fDetLevelJetSplittings, second.fDetLevelJetSplittings);
  swap(first.fPtJet, second.fPtJet);
  swap(first.fHCheckResolutionSubjets, second.fHCheckResolutionSubjets);
  swap(first.fTreeSubstructure, second.fTreeSubstructure);
}
