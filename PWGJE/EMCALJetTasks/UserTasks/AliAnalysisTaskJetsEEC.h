#ifndef ALIANALYSISTASKJetsEEC_H
#define ALIANALYSISTASKJetsEEC_H

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
class AliAnalysisTaskJetsEEC : public AliAnalysisTaskEmcalJet {
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

  AliAnalysisTaskJetsEEC();
  AliAnalysisTaskJetsEEC(const char *name);
  virtual ~AliAnalysisTaskJetsEEC();

  void UserCreateOutputObjects();
  void Terminate(Option_t *option);

// Setters
  void SetJetContainer(Int_t c) { fContainer = c; }
  void SetJetShapeType(JetShapeType t) { fJetShapeType = t; }
  void SetJetShapeSub(JetShapeSub t) { fJetShapeSub = t; }
  void SetJetSelection(JetSelectionType t) { fJetSelection = t; }
  void SetJetPtThreshold(Float_t f) { fPtThreshold = f; }
//  void SetRMatching(Float_t f) { fRMatching = f; } //MIGHT NEED
  void SetCentralitySelectionOn(Bool_t t) { fCentSelectOn = t; }
  void SetOneConstSelectionOn(Bool_t t) { fOneConstSelectOn = t; }
  void SetCheckTracksOn(Bool_t t) { fTrackCheckPlots = t; }
  void SetCheckResolution(Bool_t t) { fCheckResolution = t; }
  void SetMinPtConst(Float_t t) { fMinPtConst = t;}
  void SetHardCutoff(Float_t t) { fHardCutoff = t; }
  void SetDoTwoTrack(Bool_t t) { fDoTwoTrack = t; }
  void SetCutDoubleCounts(Bool_t t) {fCutDoubleCounts = t;}
//  void SetDoAreaIterative(Bool_t t) { fDoAreaIterative = t; }
  void SetPowerAlgorithm(Float_t t) { fPowerAlgo = t; }
  void SetMinCentrality(Float_t t) { fCentMin = t; }
  void SetMaxCentrality(Float_t t) { fCentMax = t; }
  void SetDerivativeSubtractionOrder(Int_t c) { fDerivSubtrOrder = c; }
  void SetDetLevelJetsOn(Bool_t t) { fStoreDetLevelJets = t; }
//  void SetMatchRadius(Float_t t) { fMatchR = t; } //MIGHT NEED
  void SetStoreTrig(Bool_t t) {fStoreTrig = t;}
  

protected:
  Bool_t RetrieveEventObjects();
  Bool_t Run();
  Bool_t FillHistograms();

  Float_t GetJetMass(AliEmcalJet *jet, Int_t jetContNb);
  Float_t Angularity(AliEmcalJet *jet, Int_t jetContNb); //MIGHT NOT NEED
  Float_t GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb); //MIGHT NOT NEED
  
  int GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet);
  Double_t GetDownscaleWeight(string tstring);
  void ComputeEEC(AliEmcalJet *fJet, AliJetContainer *fJetCont);
  void RunChanged(Int_t nr);
  Int_t fContainer; ///< jets to be analyzed 0 for Base, 1 for subtracted.
  Float_t fMinFractionShared; ///< only fill histos for jets if shared fraction
                              // larger than X
  JetShapeType fJetShapeType; ///< jet type to be used
  JetShapeSub fJetShapeSub;   ///< jet subtraction to be used //MIGHT NOT NEED THIS
  JetSelectionType fJetSelection; ///< Jet selection: inclusive/recoil jet
  Float_t fPtThreshold; ///<


Bool_t fCentSelectOn;      ///< switch on/off centrality selection
Float_t fCentMin;          ///< min centrality value
Float_t fCentMax;          ///< max centrality value
Bool_t fOneConstSelectOn;  ///< switch on/off one constituent selection
Bool_t fTrackCheckPlots;   ///< switch on qa plots

Bool_t fCheckResolution;   ///< check subjet energy resolution

Float_t fMinPtConst;       ///< constituent pt cutoff
Float_t fHardCutoff;       ///< hard cutoff in the iterative declustering
Bool_t fDoTwoTrack;        ///< switch to consider 2 track effects
Bool_t fCutDoubleCounts;   ///< turn off to avoid true-hybrid cuts to suppress double counting
//Bool_t fDoAreaIterative;   ///<  subtract the area in the declustering
Float_t fPowerAlgo;        ///< power of the generickt algorithm
Float_t fPhiCutValue;      ///< cuts from HBT
Float_t fEtaCutValue;      ///< cuts from HBT
Int_t fDerivSubtrOrder; //What does this do?????????
Bool_t fStoreDetLevelJets; ///< store the detector level jet quantities


///MIGHT NEED MIGHT NEED MIGHT NEED
//Bool_t fDoFillMC_EEC;      ///< to fill the matched mc plane
//  Float_t fRMatching; ///<


//MIGHT NEED MIGHT NEED
//  Float_t fMatchR; ///<the matching radius
 
  Bool_t fStoreTrig; ///<storing the trigger class
//  Bool_t fMatch; ///< do the matching in the task
  
  TH1D *EEC_hist; //initializing histogram with correlators
  TH1D *jet_pt_hist; //initializing histogram with jet pt
  TH2D *EEC_pt_hist; //initializing 2D histogram with correlators and jet pt
  TH2D *EEC_pt_hist_log; //initializing 2D histogram with correlators and jet pt on a log scale
  
  TH1D *EEEC_hist; //initializing histogram with correlators
//  TH1D *jet_pt_hist; //initializing histogram with jet pt
  TH2D *EEEC_pt_hist; //initializing 2D histogram with correlators and jet pt
  TH2D *EEEC_pt_hist_log; //initializing 2D histogram with correlators and jet pt on a log scale
  
//MIGHT NEED MIGHT NEED
//  TH1D *fHtrueAll1D;      ///<  histogram of all truth level jets
//  TH1D *fHtrueMatch1D;      ///<  histogram of matched truth level jets

private:
  AliAnalysisTaskJetsEEC(
      const AliAnalysisTaskJetsEEC &); // not implemented
  AliAnalysisTaskJetsEEC &
  operator=(const AliAnalysisTaskJetsEEC &); // not implemented

  ClassDef(AliAnalysisTaskJetsEEC, 22)
};
#endif
