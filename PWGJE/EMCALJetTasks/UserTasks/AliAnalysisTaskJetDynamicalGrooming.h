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
#include "AliYAMLConfiguration.h"

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWGJE { namespace EMCALJetTasks { class AliAnalysisTaskJetDynamicalGrooming; } }
std::ostream & operator<< (std::ostream &in, const PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming &myTask);
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming & first, PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming & second);

namespace PWGJE {
namespace EMCALJetTasks {

static const int nSubstructureVariables = 30;

class AliAnalysisTaskJetDynamicalGrooming : public AliAnalysisTaskEmcalJet
{
 public:
  enum JetShapeType_t {
    kMCTrue = 0,  // generated jets only
    kTrueDet = 1, // detector and generated jets
    kData = 2,    // raw data
    kDetEmb = 3,  // detector embedded jets
    kDetEmbPart = 4,
    kPythiaDef = 5,
    kDetEmbPartPythia = 6,
    kGenOnTheFly = 7
  };
  enum JetShapeSub_t { kNoSub = 0, kConstSub = 1, kDerivSub = 2, kEventSub = 3 };
  enum JetSelectionType_t { kInclusive = 0, kRecoil = 1 };
  enum DerivSubtrOrder_t { kSecondOrder = 0, kFirstOrder = 1 };

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
  void SetJetContainer(Int_t c) { fContainer = c; }
  void SetMinFractionShared(Double_t f) { fMinFractionShared = f; }
  void SetJetShapeType(JetShapeType_t t) { fJetShapeType = t; }
  void SetJetShapeSub(JetShapeSub_t t) { fJetShapeSub = t; }
  void SetJetSelection(JetSelectionType_t t) { fJetSelection = t; }
  void SetJetPtThreshold(Float_t f) { fPtThreshold = f; }
  void SetRMatching(Float_t f) { fRMatching = f; }

  void SetCentralitySelectionOn(Bool_t t) { fCentSelectOn = t; }
  void SetOneConstSelectionOn(Bool_t t) { fOneConstSelectOn = t; }
  void SetCheckTracksOn(Bool_t t) { fTrackCheckPlots = t; }
  void SetFillLundMC(Bool_t t) { fDoFillMCLund = t; }
  void SetCheckResolution(Bool_t t) { fCheckResolution = t; }
  void SetSubjetCutoff(Float_t t) { fSubjetCutoff = t; }
  void SetMinPtConst(Float_t t) { fMinPtConst = t; }
  void SetHardCutoff(Float_t t) { fHardCutoff = t; }
  void SetDoTwoTrack(Bool_t t) { fDoTwoTrack = t; }
  void SetCutDoubleCounts(Bool_t t) { fCutDoubleCounts = t; }
  void SetDoAreaIterative(Bool_t t) { fDoAreaIterative = t; }
  void SetPowerAlgorithm(Float_t t) { fPowerAlgo = t; }
  void SetMagFieldPol(Float_t t) { fMagFieldPolarity = t; }
  void SetMinCentrality(Float_t t) { fCentMin = t; }
  void SetMaxCentrality(Float_t t) { fCentMax = t; }
  void SetDerivativeSubtractionOrder(Int_t c) { fDerivSubtrOrder = c; }
  void SetDetLevelJetsOn(Bool_t t) { fStoreDetLevelJets = t; }

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

  void SetupTree();

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

  Float_t GetJetMass(AliEmcalJet* jet, Int_t jetContNb);
  Float_t Angularity(AliEmcalJet* jet, Int_t jetContNb);
  Float_t GetJetAngularity(AliEmcalJet* jet, Int_t jetContNb);
  Double_t RelativePhi(Double_t mphi, Double_t vphi);
  void IterativeParents(AliEmcalJet* fJet, AliJetContainer* fJetCont);
  //void IterativeParentsAreaBased(AliEmcalJet* fJet, AliJetContainer* fJetCont);
  void IterativeParentsMCAverage(AliEmcalJet* fJet, Int_t km, Double_t& aver1, Double_t& aver2, Double_t& aver3,
                  Double_t& aver4, Double_t& aver5, Double_t& aver6, Double_t& aver7, Double_t& aver8);
  void CheckSubjetResolution(AliEmcalJet* fJet, AliJetContainer* fJetCont, AliEmcalJet* fJetM,
                AliJetContainer* fJetContM);
  Bool_t CheckClosePartner(Int_t index, AliEmcalJet* fJet, AliVParticle* fTrack, AliParticleContainer* fTrackCont);

  // Basic configuration
  PWG::Tools::AliYAMLConfiguration fYAMLConfig; ///< YAML configuration file.
  bool fConfigurationInitialized;     ///<  True if the task configuration has been successfully initialized.

  Int_t fContainer;               // jets to be analyzed 0 for Base, 1 for subtracted.
  Float_t fMinFractionShared;     // only fill histos for jets if shared fraction
                  // larger than X
  JetShapeType_t fJetShapeType;     // jet type to be used
  JetShapeSub_t fJetShapeSub;       // jet subtraction to be used
  JetSelectionType_t fJetSelection; // Jet selection: inclusive/recoil jet
  Float_t fShapesVar[nSubstructureVariables];         // jet shapes used for the tagging
  Float_t fPtThreshold;
  Float_t fRMatching;

  Bool_t fCentSelectOn;      // switch on/off centrality selection
  Float_t fCentMin;          // min centrality value
  Float_t fCentMax;          // max centrality value
  Bool_t fOneConstSelectOn;  // switch on/off one constituent selection
  Bool_t fTrackCheckPlots;   // switch on qa plots
  Bool_t fDoFillMCLund;      // to fill the matched mc plane
  Bool_t fCheckResolution;   // check subjet energy resolution
  Float_t fSubjetCutoff;     // angular cutoff for subjets at det/gen level
  Float_t fMinPtConst;       // constituent pt cutoff
  Float_t fHardCutoff;       // hard cutoff in the iterative declustering
  Bool_t fDoTwoTrack;        // switch to consider 2 track effects
  Bool_t fCutDoubleCounts;   // turn off to avoid true-hybrid cuts to suppress double counting
  Bool_t fDoAreaIterative;   // subtract the area in the declustering
  Float_t fPowerAlgo;        // power of the generickt algorithm
  Float_t fPhiCutValue;      // cuts from HBT
  Float_t fEtaCutValue;      // cuts from HBT
  Float_t fMagFieldPolarity; // polarity, to calculate phimin
  Int_t fDerivSubtrOrder;
  Bool_t fStoreDetLevelJets; // store the detector level jet quantities

  TH1F* fPtJet;

  THnSparse* fHLundIterative;          //       iterative declustering
  THnSparse* fHLundIterativeMC;        //       iterative declustering
  THnSparse* fHLundIterativeMCDet;     //       iterative declustering
  THnSparse* fHCheckResolutionSubjets; //     to evaluate energy resolution of subjets
                     //     as function fo apperture angle

  TTree* fTreeSubstructure; // Tree with tagging variables subtracted MC or true
               // MC or raw

 private:
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetDynamicalGrooming, 1)  // Jet dynamical grooming
  /// \endcond
};

/*class AliJetDynamicalGroomingTree : public TNamed {

  AliJetDynamicalGroomingTree();

 private:
  TTree* fJetTree;                             //!<! tree structure
  Bool_t fInitialized;                         ///< init state of tree

  /// \cond CLASSIMP
  ClassDef(AliEmcalJetTree, 12) // Jet dynamical grooming tree
  /// \endcond
};*/

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif  /* ALIANALYSISTASKJETDYNAMICALGROOMING_H */
