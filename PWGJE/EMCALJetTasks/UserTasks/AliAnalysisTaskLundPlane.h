#ifndef ALIANALYSISTASKLUNDPLANE_H
#define ALIANALYSISTASKLUNDPLANE_H

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
class AliFJWrapper;
#include "AliAnalysisTaskEmcalJet.h"
#include "AliFJWrapper.h"
#include <vector>
class AliAnalysisTaskLundPlane : public AliAnalysisTaskEmcalJet {
public:
  enum JetShapeType {
    kMCTrue = 0,  // generated jets only
    kTrueDet = 1, // detector and generated jets
    kData = 2,    // raw data
    kDetEmb = 3,  // detector embedded jets
    kDetEmbPart = 4,
    kPythiaDef = 5,
    kDetEmbPartPythia = 6,
    kGenOnTheFly = 7
      };
  enum JetShapeSub { kNoSub = 0, kConstSub = 1, kDerivSub = 2, kEventSub=3 };
  enum JetSelectionType { kInclusive = 0, kRecoil = 1 };

  enum DerivSubtrOrder { kSecondOrder = 0, kFirstOrder = 1 };

  AliAnalysisTaskLundPlane();
  AliAnalysisTaskLundPlane(const char *name);
  virtual ~AliAnalysisTaskLundPlane();

  void UserCreateOutputObjects();
  void Terminate(Option_t *option);

  // Setters
  void SetJetContainer(Int_t c) { fContainer = c; }
  void SetMinFractionShared(Double_t f) { fMinFractionShared = f; }
  void SetJetShapeType(JetShapeType t) { fJetShapeType = t; }
  void SetJetShapeSub(JetShapeSub t) { fJetShapeSub = t; }
  void SetJetSelection(JetSelectionType t) { fJetSelection = t; }
  void SetJetPtThreshold(Float_t f) { fPtThreshold = f; }
  void SetRMatching(Float_t f) { fRMatching = f; }

  void SetCentralitySelectionOn(Bool_t t) { fCentSelectOn = t; }
  void SetOneConstSelectionOn(Bool_t t) { fOneConstSelectOn = t; }
  void SetCheckTracksOn(Bool_t t) { fTrackCheckPlots = t; }
  void SetFillLundMC(Bool_t t) { fDoFillMCLund = t; }
  void SetCheckResolution(Bool_t t) { fCheckResolution = t; }
  void SetSubjetCutoff(Float_t t) { fSubjetCutoff = t; }
  void SetMinPtConst(Float_t t) { fMinPtConst = t;}
  void SetHardCutoff(Float_t t) { fHardCutoff = t; }
  void SetDoTwoTrack(Bool_t t) { fDoTwoTrack = t; }
  void SetCutDoubleCounts(Bool_t t) {fCutDoubleCounts = t;}
  void SetDoAreaIterative(Bool_t t) { fDoAreaIterative = t; }
  void SetPowerAlgorithm(Float_t t) { fPowerAlgo = t; }
  void SetMagFieldPol(Float_t t) { fMagFieldPolarity = t; }
  void SetMinCentrality(Float_t t) { fCentMin = t; }
  void SetMaxCentrality(Float_t t) { fCentMax = t; }
  void SetDerivativeSubtractionOrder(Int_t c) { fDerivSubtrOrder = c; }
  void SetDetLevelJetsOn(Bool_t t) { fStoreDetLevelJets = t; }
  void SetDoSubJetStudy(Bool_t t) { fDoSubJet = t; }

protected:
  Bool_t RetrieveEventObjects();
  Bool_t Run();
  Bool_t FillHistograms();

  Float_t GetJetMass(AliEmcalJet *jet, Int_t jetContNb);
  Float_t Angularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb);
  Double_t RelativePhi(Double_t mphi, Double_t vphi);
  void IterativeDeclustering(AliEmcalJet *fJet, AliJetContainer *fJetCont);

  void IterativeDeclusteringMC(AliEmcalJet *fJet, Int_t km);

  Int_t fContainer; ///< jets to be analyzed 0 for Base, 1 for subtracted.
  Float_t fMinFractionShared; ///< only fill histos for jets if shared fraction
                              // larger than X
  JetShapeType fJetShapeType; ///< jet type to be used
  JetShapeSub fJetShapeSub;   ///< jet subtraction to be used
  JetSelectionType fJetSelection; ///< Jet selection: inclusive/recoil jet
  Float_t fShapesVar[20];         ///< jet shapes used for the tagging
  Float_t fPtThreshold; ///<
  Float_t fRMatching; ///<

  Bool_t fCentSelectOn;      ///< switch on/off centrality selection
  Float_t fCentMin;          ///< min centrality value
  Float_t fCentMax;          ///< max centrality value
  Bool_t fOneConstSelectOn;  ///< switch on/off one constituent selection
  Bool_t fTrackCheckPlots;   ///< switch on qa plots
  Bool_t fDoFillMCLund;      ///< to fill the matched mc plane
  Bool_t fCheckResolution;   ///< check subjet energy resolution
  Float_t fSubjetCutoff;     ///< angular cutoff for subjets at det/gen level
  Float_t fMinPtConst;       ///< constituent pt cutoff
  Float_t fHardCutoff;       ///< hard cutoff in the iterative declustering
  Bool_t fDoTwoTrack;        ///< switch to consider 2 track effects
  Bool_t fCutDoubleCounts;   ///< turn off to avoid true-hybrid cuts to suppress double counting
  Bool_t fDoAreaIterative;   ///<  subtract the area in the declustering
  Float_t fPowerAlgo;        ///< power of the generickt algorithm
  Float_t fPhiCutValue;      ///< cuts from HBT
  Float_t fEtaCutValue;      ///< cuts from HBT
  Float_t fMagFieldPolarity; ///< polarity, to calculate phimin
  Int_t fDerivSubtrOrder; ///<
  Bool_t fStoreDetLevelJets; ///< store the detector level jet quantities
  Bool_t fDoSubJet; ///< store the detector level jet quantities


 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_angle;
 
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_kt;

 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_z;
 
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_energy;
 
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_eta1;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_eta2;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_phi1;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_phi2;
 

std::vector<std::vector<Double_t>>            fShapesVar_Splittings_angle_part;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_kt_part;
std::vector<std::vector<Double_t>>            fShapesVar_Splittings_z_part; 
std::vector<std::vector<Double_t>>            fShapesVar_Splittings_energy_part;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_eta1_part;
std::vector<std::vector<Double_t>>            fShapesVar_Splittings_eta2_part;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_phi1_part;
std::vector<std::vector<Double_t>>            fShapesVar_Splittings_phi2_part;
 
 Double_t                                      fShapesVar_Splittings_ptjet;
 Double_t                                      fShapesVar_Splittings_ptjet_part; 

  TTree *fTreeSplittings; ///< Tree with tagging variables subtracted MC or true
                            // MC or raw

private:
  AliAnalysisTaskLundPlane(
      const AliAnalysisTaskLundPlane &); // not implemented
  AliAnalysisTaskLundPlane &
  operator=(const AliAnalysisTaskLundPlane &); // not implemented

  ClassDef(AliAnalysisTaskLundPlane, 9)
};
#endif
