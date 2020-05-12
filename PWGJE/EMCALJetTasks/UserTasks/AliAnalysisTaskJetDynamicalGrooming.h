#ifndef ALIANALYSISTASKJETDYNAMICALGROOMING_H
#define ALIANALYSISTASKJETDYNAMICALGROOMING_H

/**
 * @class AliAnalysisTaskJetDynamicalGrooming
 * @brief Jet substructure with dynamical grooming.
 *
 * %Analysis task for jet substructure utilizing dynamical grooming.
 * Adapted from AliAnalysisTaskNewJetSubstructure.
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
 * @date 4 Feb 2020
 */

#include <map>
#include <string>

class TH1;
class TH2;
class TH3;
class TH3F;
class TTree;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJetFinder;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliFJWrapper.h"
#include "AliEmcalParticleJetConstituent.h"
#include "AliYAMLConfiguration.h"

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWGJE {
namespace EMCALJetTasks {
namespace SubstructureTree {
  class Subjets;
  class JetSplittings;
  class JetConstituents;
  class JetSubstructureSplittings;
}
  class AliAnalysisTaskJetDynamicalGrooming;
} // namespace EMCALJetTasks
} // namespace PWGJE
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::SubstructureTree::Subjets& myTask);
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::SubstructureTree::JetSplittings& myTask);
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::SubstructureTree::JetConstituents& myTask);
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::SubstructureTree::JetSubstructureSplittings& myTask);
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& myTask);
void swap(PWGJE::EMCALJetTasks::SubstructureTree::Subjets& first,
     PWGJE::EMCALJetTasks::SubstructureTree::Subjets& second);
void swap(PWGJE::EMCALJetTasks::SubstructureTree::JetSplittings& first,
     PWGJE::EMCALJetTasks::SubstructureTree::JetSplittings& second);
void swap(PWGJE::EMCALJetTasks::SubstructureTree::JetConstituents& first,
     PWGJE::EMCALJetTasks::SubstructureTree::JetConstituents& second);
void swap(PWGJE::EMCALJetTasks::SubstructureTree::JetSubstructureSplittings& first,
     PWGJE::EMCALJetTasks::SubstructureTree::JetSubstructureSplittings& second);
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& first,
     PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& second);

namespace PWGJE {
namespace EMCALJetTasks {
namespace SubstructureTree {

class Subjets {
 public:
  // TODO: Fully update and document!
  Subjets();
  // Additional constructors
  Subjets(const Subjets & other);
  Subjets& operator=(Subjets other);
  friend void ::swap(Subjets & first, Subjets & second);
  // Avoid implementing move since c++11 is not allowed in the header
  virtual ~Subjets() {}

  /// Reset the properties for the next filling of the tree.
  bool Clear();

  // Getters and setters
  void AddSubjet(const unsigned short splittingNodeIndex, const bool partOfIterativeSplitting,
          const std::vector<unsigned short>& constituentIndices);
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  std::tuple<unsigned short, bool, const std::vector<unsigned short>> GetSubjet(int i) const;
  #endif

  // Printing
  std::string toString() const;
  friend std::ostream & ::operator<<(std::ostream &in, const Subjets &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;

 protected:
  std::vector<unsigned short> fSplittingNodeIndex;        ///<  Index of the parent splitting node.
  std::vector<bool> fPartOfIterativeSplitting;            ///<  True if the splitting is follow an iterative splitting.
  std::vector<std::vector<unsigned short>> fConstituentIndices;        ///<  Constituent jet indices (ie. index by the stored jet constituents, not the global index).

  /// \cond CLASSIMP
  ClassDef(Subjets, 2) // Subjets from splittings.
  /// \endcond
};

class JetSplittings {
 public:
  JetSplittings();
  // Additional constructors
  JetSplittings(const JetSplittings & other);
  JetSplittings& operator=(JetSplittings other);
  friend void ::swap(JetSplittings & first, JetSplittings & second);
  // Avoid implementing move since c++11 is not allowed in the header
  virtual ~JetSplittings() {}

  /// Reset the properties for the next filling of the tree.
  bool Clear();

  // Getters and setters
  void AddSplitting(float kt, float deltaR, float z, short parentIndex);
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  std::tuple<float, float, float, short> GetSplitting(int i) const;
  #endif
  unsigned int GetNumberOfSplittings() const { return fKt.size(); }

  // Printing
  std::string toString() const;
  friend std::ostream & ::operator<<(std::ostream &in, const JetSplittings &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;

 protected:
  std::vector<float> fKt;             ///<  kT between the subjets.
  std::vector<float> fDeltaR;         ///<  Delta R between the subjets.
  std::vector<float> fZ;              ///<  Momentum sharing of the splitting.
  std::vector<short> fParentIndex;    ///<  Index of the parent splitting.

  /// \cond CLASSIMP
  ClassDef(JetSplittings, 1) // Jet splittings.
  /// \endcond
};

class JetConstituents
{
 public:
  // TODO: Fully update and document!
  JetConstituents();
  // Additional constructors
  JetConstituents(const JetConstituents & other);
  JetConstituents& operator=(JetConstituents other);
  friend void ::swap(JetConstituents & first, JetConstituents & second);
  // Avoid implementing move since c++11 is not allowed in the header
  virtual ~JetConstituents() {}

  /// Reset the properties for the next filling of the tree.
  bool Clear();

  // Getters and setters
  void AddJetConstituent(const AliVParticle* part, const int & id);
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  std::tuple<float, float, float, int> GetJetConstituent(int i) const;
  #endif
  static const int GetGlobalIndexOffset() { return fgkGlobalIndexOffset; }

  // Printing
  std::string toString() const;
  friend std::ostream & ::operator<<(std::ostream &in, const JetConstituents &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;

 protected:
  static const int fgkGlobalIndexOffset;  ///<  Offset for GlobalIndex values in the ID to ensure it never conflicts with the label.

  std::vector<float> fPt;                 ///<  Jet constituent pt
  std::vector<float> fEta;                ///<  Jet constituent eta
  std::vector<float> fPhi;                ///<  Jet constituent phi
  std::vector<int> fID;                   ///<  Jet constituent identifier. MClabel (via GetLabel()) or global index (with offset defined here).

  /// \cond CLASSIMP
  ClassDef(JetConstituents, 2) // Jet constituents.
  /// \endcond
};

/**
 * @class JetSubstructureSplittings
 * @brief Jet substructure splittings.
 *
 * Jet substructure splitting properties. There is sufficient information to calculate any
 * additional splitting properties.
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
 * @date 9 Feb 2020
 */
class JetSubstructureSplittings {
 public:
  JetSubstructureSplittings();
  // Additional constructors
  JetSubstructureSplittings(const JetSubstructureSplittings & other);
  JetSubstructureSplittings& operator=(JetSubstructureSplittings other);
  friend void ::swap(JetSubstructureSplittings & first, JetSubstructureSplittings & second);
  // Avoid implementing move since c++11 is not allowed in the header
  virtual ~JetSubstructureSplittings() {}

  /// Reset the properties for the next filling of the tree.
  bool Clear();

  // Setters
  void SetJetPt(float pt) { fJetPt = pt; }
  void AddJetConstituent(const AliVParticle* part, const int & id);
  void AddSplitting(float kt, float deltaR, float z, short parentIndex);
  void AddSubjet(const unsigned short splittingNodeIndex, const bool partOfIterativeSplitting,
          const std::vector<unsigned short>& constituentIndices);
  // Getters
  float GetJetPt() { return fJetPt; }
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  std::tuple<float, float, float, int> GetJetConstituent(int i) const;
  std::tuple<float, float, float, short> GetSplitting(int i) const;
  std::tuple<unsigned short, bool, const std::vector<unsigned short>> GetSubjet(int i) const;
  #endif
  unsigned int GetNumberOfSplittings() { return fJetSplittings.GetNumberOfSplittings(); }

  // Printing
  std::string toString() const;
  friend std::ostream & ::operator<<(std::ostream &in, const JetSubstructureSplittings &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;

 private:
  // Jet properties
  float fJetPt;                                           ///<  Jet pt.
  SubstructureTree::JetConstituents fJetConstituents;     ///<  Jet constituents
  SubstructureTree::JetSplittings fJetSplittings;         ///<  Jet splittings.
  SubstructureTree::Subjets fSubjets;                     ///<  Subjets within the jet.

  /// \cond CLASSIMP
  ClassDef(JetSubstructureSplittings, 2) // Jet splitting properties.
  /// \endcond
};

} /* namespace SubstructureTree */

class AliAnalysisTaskJetDynamicalGrooming : public AliAnalysisTaskEmcalJet
{
 public:
  enum JetShapeType_t {
    kMCTrue = 0,            //!<! Generated jets only
    kTrueDet = 1,           //!<! Detector and generated jets
    kData = 2,              //!<! Raw data
    kDetEmb = 3,            //!<! Detector embedded jets
    kDetEmbPart = 4,        //!<!
    kPythiaDef = 5,         //!<!
    kDetEmbPartPythia = 6,  //!<!
    kGenOnTheFly = 7        //!<!
  };
  enum JetShapeSub_t { kNoSub = 0, kConstSub = 1, kDerivSub = 2, kEventSub = 3 };
  enum JetSelectionType_t { kInclusive = 0, kRecoil = 1 };
  enum DerivSubtrOrder_t { kSecondOrder = 0, kFirstOrder = 1 };
  static const std::map<std::string, JetShapeType_t> fgkJetShapeTypeMap; //!<! Map from name to jet shape type used with the YAML config
  static const std::map<std::string, JetShapeSub_t> fgkJetShapeSubMap; //!<! Map from name to jet shape sub used with the YAML config
  static const std::map<std::string, JetSelectionType_t> fgkJetSelectionMap; //!<! Map from name to jet selection used with the YAML config
  static const std::map<std::string, DerivSubtrOrder_t> fgkDerivSubtrOrderMap; //!<! Map from name to derivative subtracter order used with the YAML config

  AliAnalysisTaskJetDynamicalGrooming();
  AliAnalysisTaskJetDynamicalGrooming(const char* name);
  // Additional constructors
  AliAnalysisTaskJetDynamicalGrooming(const AliAnalysisTaskJetDynamicalGrooming & other);
  AliAnalysisTaskJetDynamicalGrooming& operator=(AliAnalysisTaskJetDynamicalGrooming other);
  friend void ::swap(AliAnalysisTaskJetDynamicalGrooming & first, AliAnalysisTaskJetDynamicalGrooming & second);
  // Avoid implementing move since c++11 is not allowed in the header
  virtual ~AliAnalysisTaskJetDynamicalGrooming() {}

  void UserCreateOutputObjects();

  // Setters
  void SetMinFractionShared(Double_t f) { fMinFractionShared = f; }
  void SetJetShapeType(JetShapeType_t t) { fJetShapeType = t; }
  void SetJetShapeSub(JetShapeSub_t t) { fJetShapeSub = t; }
  void SetJetSelection(JetSelectionType_t t) { fJetSelection = t; }
  void SetJetPtThreshold(Float_t f) { fPtThreshold = f; }
  void SetRMatching(Float_t f) { fRMatching = f; }

  void SetCentralitySelectionOn(Bool_t t) { fCentSelectOn = t; }
  void SetCheckResolution(Bool_t t) { fCheckResolution = t; }
  void SetSubjetCutoff(Float_t t) { fSubjetCutoff = t; }
  void SetMinPtConst(Float_t t) { fMinPtConst = t; }
  void SetHardCutoff(Float_t t) { fHardCutoff = t; }
  void SetDoTwoTrack(Bool_t t) { fDoTwoTrack = t; }
  void SetCutDoubleCounts(Bool_t t) { fCutDoubleCounts = t; }
  void SetMagFieldPol(Float_t t) { fMagFieldPolarity = t; }
  void SetMinCentrality(Float_t t) { fCentMin = t; }
  void SetMaxCentrality(Float_t t) { fCentMax = t; }
  void SetDerivativeSubtractionOrder(Int_t c) { fDerivSubtrOrder = c; }
  void SetDetLevelJetsOn(Bool_t t) { fStoreDetLevelJets = t; }
  void SetStoreRecursiveJetSplittings(bool t = true) { fStoreRecursiveSplittings = t; }

  // Initialize the task
  // Configuration is handled via the YAML configuration file
  bool Initialize();
  void AddConfigurationFile(const std::string & configurationPath, const std::string & configName = "") { fYAMLConfig.AddConfiguration(configurationPath, configName); }

  //static AliAnalysisTaskJetDynamicalGrooming * AddTaskJetDynamicalGrooming(const std::string & suffix = "");

  static AliAnalysisTaskJetDynamicalGrooming* AddTaskJetDynamicalGrooming(
   const char* njetsBase, const char* njetsUS, const char* njetsTrue, const char* njetsPartLevel, const Double_t R,
   const char* nrhoBase, const char* ntracks, const char* ntracksUS, const char* ntracksPartLevel,
   const char* nclusters, const char* ntracksTrue, const char* type, const char* CentEst, Int_t pSel,
   AliAnalysisTaskJetDynamicalGrooming::JetShapeType_t jetShapeType = AliAnalysisTaskJetDynamicalGrooming::kMCTrue,
   AliAnalysisTaskJetDynamicalGrooming::JetShapeSub_t jetShapeSub = AliAnalysisTaskJetDynamicalGrooming::kNoSub,
   AliAnalysisTaskJetDynamicalGrooming::JetSelectionType_t jetSelection =
    AliAnalysisTaskJetDynamicalGrooming::kInclusive,
   Float_t minpTHTrigger = 0., Float_t maxpTHTrigger = 0., Float_t acut = 0.6,
   AliAnalysisTaskJetDynamicalGrooming::DerivSubtrOrder_t derivSubtrOrder =
    AliAnalysisTaskJetDynamicalGrooming::kSecondOrder,
   const std::string& suffix = "");

  // Printing
  std::string toString() const;
  friend std::ostream & ::operator<<(std::ostream &in, const AliAnalysisTaskJetDynamicalGrooming &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;

  // Helpers
  std::string DetermineOutputContainerName(std::string containerName) const;

 protected:
  Bool_t RetrieveEventObjects();
  Bool_t Run();
  Bool_t FillHistograms();

  void RetrieveAndSetTaskPropertiesFromYAMLConfig();
  void SetupTree();

  template<typename T>
  std::string GetKeyFromMapValue(const T & value, const std::map<std::string, T> & valuesMap) const;

  /**
   * Calculate the dynamical grooming measure of the hardness.
   *
   * @param[in] subjet1 Leading subjet.
   * @param[in] subjet1 Subleading subjet.
   * @param[in] parent Parent of the subjets.
   * @param[in] R Jet resolution parameter.
   * @param[in] a Characteristic parameter for calculating the hardness. Different values correspond to different measures. See the paper.
   *
   */
  double DynamicalGrooming(const fastjet::PseudoJet & subjet1, const fastjet::PseudoJet & subjet2, const fastjet::PseudoJet & parent, const double R, const double a) const;
  /// Calculate zDrop (a = 0.1) for the most symmetric momentum sharing.
  double CalculateZDrop(const fastjet::PseudoJet & subjet1, const fastjet::PseudoJet & subjet2, const fastjet::PseudoJet & parent, const double R) const;
  /// Calculate KtDrop (a = 1) for the hardest kt.
  double CalculateKtDrop(const fastjet::PseudoJet & subjet1, const fastjet::PseudoJet & subjet2, const fastjet::PseudoJet & parent, const double R) const;
  /// Calculate TimeDrop (a = 2) for the earliest splitting.
  double CalculateTimeDrop(const fastjet::PseudoJet & subjet1, const fastjet::PseudoJet & subjet2, const fastjet::PseudoJet & parent, const double R) const;

  void IterativeParents(AliEmcalJet* jet, SubstructureTree::JetSubstructureSplittings& jetSplittings, bool isData);
  void ExtractJetSplittings(SubstructureTree::JetSubstructureSplittings & jetSplittings, fastjet::PseudoJet & inputJet, int splittingNodeIndex, bool followingIterativeSplitting);
  void CheckSubjetResolution(AliEmcalJet* fJet, AliEmcalJet* fJetM);
  bool CheckClosePartner(const AliEmcalJet* jet, const AliVParticle * part1);

  // Basic configuration
  PWG::Tools::AliYAMLConfiguration fYAMLConfig; ///<  YAML configuration file.
  bool fConfigurationInitialized;               ///<  True if the task configuration has been successfully initialized.

  Float_t fMinFractionShared;                   ///<  Only fill histos for jets if shared fraction larger than X
  JetShapeType_t fJetShapeType;                 ///<  Jet type to be used
  JetShapeSub_t fJetShapeSub;                   ///<  Jet subtraction to be used
  JetSelectionType_t fJetSelection;             ///<  Jet selection: inclusive/recoil jet
  Float_t fPtThreshold;                         ///<  Minimum jet pt.
  Float_t fRMatching;                           ///<  R matching distance.

  Bool_t fCentSelectOn;      ///<  switch on/off centrality selection
  Float_t fCentMin;          ///<  min centrality value
  Float_t fCentMax;          ///<  max centrality value
  Bool_t fCheckResolution;   ///<  check subjet energy resolution
  Float_t fSubjetCutoff;     ///<  angular cutoff for subjets at det/gen level
  Float_t fMinPtConst;       ///<  constituent pt cutoff
  Float_t fHardCutoff;       ///<  hard cutoff in the iterative declustering
  Bool_t fDoTwoTrack;        ///<  switch to consider 2 track effects
  Bool_t fCutDoubleCounts;   ///<  turn off to avoid true-hybrid cuts to suppress double counting
  Float_t fPhiCutValue;      ///<  cuts from HBT
  Float_t fEtaCutValue;      ///<  cuts from HBT
  Float_t fMagFieldPolarity; ///<  polarity, to calculate phimin
  Int_t fDerivSubtrOrder;    ///<  Order of the derivative subtraction.
  Bool_t fStoreDetLevelJets; ///<  If True, store the detector level jet quantities
  bool fStoreRecursiveSplittings; ///<  If true, recursive splittings will be stored.

  // Tree variables
  SubstructureTree::JetSubstructureSplittings fDataJetSplittings;       ///<  Data jet splittings.
  SubstructureTree::JetSubstructureSplittings fMatchedJetSplittings;    ///<  Matched jet splittings.
  SubstructureTree::JetSubstructureSplittings fDetLevelJetSplittings;   ///<  Det level (intermediate match) jet splittings.

  TH1F* fPtJet;                                       //!<! Jet pt

  THnSparse* fHCheckResolutionSubjets;                //!<! To evaluate energy resolution of subjets as function of aperture angle

  TTree* fTreeSubstructure; //!<! Tree with tagging variables subtracted MC or true MC or raw

 private:
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetDynamicalGrooming, 2)  // Jet dynamical grooming
  /// \endcond
};

template<typename T>
std::string AliAnalysisTaskJetDynamicalGrooming::GetKeyFromMapValue(const T & value, const std::map<std::string, T> & valuesMap) const
{
  auto findResult = std::find_if(std::begin(valuesMap), std::end(valuesMap), [&](const std::pair<std::string, T> &pair)
  {
    return pair.second == value;
  });
  std::string returnValue = "Invalid value: " + std::to_string(static_cast<int>(value));
  if (findResult != std::end(valuesMap))
  {
    returnValue = findResult->first;
  }
  return returnValue;
}

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif  /* ALIANALYSISTASKJETDYNAMICALGROOMING_H */
