#ifndef ALIANALYSISTASKSUBJETFRACTION_H
#define ALIANALYSISTASKSUBJETFRACTION_H

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
#include "AliClusterContainer.h"
#include "TF1.h"

class AliAnalysisTaskSubJetFraction : public AliAnalysisTaskEmcalJet {
 public:
  
  enum JetShapeType {
    kTrue = 0,   // generated jets only 
    kTrueDet =1,  // detector and generated jets  
    kData   = 2,  // raw data 
    kDetEmbPart = 3,  //detector embedded jets
    kSim = 4,
    kGenOnTheFly=5 //fast truth level MC
  };
  enum JetShapeSub {
    kNoSub = 0, 
    kConstSub = 1,
    kDerivSub = 2 
  };
  enum JetSelectionType {
    kInclusive = 0,
    kRecoil = 1
  };
  enum DerivSubtrOrder {
    kSecondOrder = 0,
    kFirstOrder = 1
  };
  enum TreeSize {
    nVar = 28
  };

  AliAnalysisTaskSubJetFraction();
  AliAnalysisTaskSubJetFraction(const char *name);
  virtual ~AliAnalysisTaskSubJetFraction();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)                             { fContainer     = c   ;}
  void SetMinFractionShared(Double_t f)                     { fMinFractionShared = f   ;}
  void SetJetShapeType(JetShapeType t)                      { fJetShapeType       = t   ;}
  void SetJetShapeSub(JetShapeSub t)                        { fJetShapeSub     = t   ;}
  void SetJetSelection(JetSelectionType t)                  { fJetSelection    = t   ;}
  void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ;}
  void SetJetRadius(Double_t JetRadius)                     {fJetRadius=JetRadius;}
  void SetRMatching(Float_t f)                              { fRMatching = f ;}
  void SetPtTriggerSelections(Float_t PtMin, Float_t PtMax) { fPtMinTriggerHadron = PtMin; fPtMaxTriggerHadron = PtMax;}
  void SetAngularWindowRecoilJet (Float_t Window)                {fRecoilAngularWindow = Window; }
  Float_t GetMinPtTriggerSelection()                        {return fPtMinTriggerHadron;}
  Float_t GetMaxPtTriggerSelection()                        {return fPtMaxTriggerHadron;}
  void SetCentralitySelectionOn(Bool_t t)                   { fCentSelectOn = t;}
  void SetMinCentrality(Float_t t)                          { fCentMin = t ;}
  void SetMaxCentrality(Float_t t)                          { fCentMax = t ;}
  void SetSemigoodCorrect(Int_t yesno)                 {fSemigoodCorrect=yesno;}
  void SetHolePos(Float_t poshole)                        { fHolePos = poshole;}
  void SetHoleWidth(Float_t holewidth)                  { fHoleWidth = holewidth;}
  void SetSubJetAlgorithm(Int_t SubJetAlgorithm)        {fSubJetAlgorithm=SubJetAlgorithm;}
  void SetSubJetRadius(Float_t SubJetRadius)            {fSubJetRadius=SubJetRadius;}
  void SetSubJetMinPt(Float_t SubJetMinPt)              {fSubJetMinPt=SubJetMinPt;}
  void SetRMatched(Double_t RMatched)                     {fRMatched=RMatched;}
  void SetSharedFractionPtMin(Double_t SharedFractionPtMin) {fSharedFractionPtMin=SharedFractionPtMin;}
  void SetDerivativeSubtractionOrder(Int_t Order)              {fDerivSubtrOrder = Order;}
  void SetFullTree(Bool_t FullTree)                         {fFullTree = FullTree;}
  void SetBetaSD(Double_t BetaSD)                           {fBeta_SD = BetaSD;}
  void SetZCut(Double_t ZCut)                               {fZCut = ZCut;}
  void SetReclusteringAlgorithm(Int_t ReclusteringAlgorithm)     {fReclusteringAlgorithm = ReclusteringAlgorithm;}
  void SetSoftDropOn(Int_t SoftDropOn)                           {fSoftDropOn = SoftDropOn;}
  void SetMLOn(Int_t MLOn)                           {fMLOn = MLOn;}
  Int_t GetSoftDropOn()                                          {return fSoftDropOn;}
  Int_t GetMLOn()                                          {return fMLOn;}
  void SetRandomisationEqualPt(Bool_t RandmosationEqualPt)       {fRandomisationEqualPt = RandmosationEqualPt;}
  
  void SetNsubUnNormMeasure( Bool_t NsubMeasure)              {fNsubMeasure= NsubMeasure;}



AliAnalysisTaskSubJetFraction* AddTaskAliAnalysisTaskSubJetFraction(const char * njetsData, //data jets
								    const char * njetsTrue, //Pyhthia Particle Level
								    const char * njetsDet,
								    const char * njetsHybridUs,
								    const char * njetsHybridS,
								    const Double_t R,
								    const char * nrhoBase, 
								    const char * ntracksData,
                                                                    const char * ntracksTrue,
                                                                    const char * ntracksDet, 
								    const char * ntracksHybridUs,
								    const char * ntracksHybridS,
								    const char * nclusters,
								    const char *type,				      
								    const char *CentEst,
								    Double_t fSharedFractionPtMin,
								    Int_t SubJetAlgorithm,
								    Float_t SubJetRadius,
								    Float_t SubJetMinPt,
								    Int_t       pSel,
								    TString     trigClass      = "",
								    TString     kEmcalTriggers = "",
								    TString     tag            = "",
								    AliAnalysisTaskSubJetFraction::JetShapeType jetShapeType = AliAnalysisTaskSubJetFraction::kTrue, // tobefixedbyauthor
								    AliAnalysisTaskSubJetFraction::JetShapeSub jetShapeSub = AliAnalysisTaskSubJetFraction::kNoSub, // tobefixedbyauthor
								    AliAnalysisTaskSubJetFraction::JetSelectionType jetSelection =AliAnalysisTaskSubJetFraction::kInclusive, // tobefixedbyauthor
								    Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0., AliAnalysisTaskSubJetFraction::DerivSubtrOrder derivSubtrOrder = AliAnalysisTaskSubJetFraction::kSecondOrder, Int_t SoftDropOn=0, Int_t MLOn=0);





  

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Int_t                               SelectTriggerHadron(Float_t PtMin, Float_t PtMax);
  Double_t                            RelativePhiEventPlane(Double_t EventPlane, Double_t Phi);
  Double_t                            RelativePhi(Double_t Phi1, Double_t Phi2);
  Double_t                            Angularity(AliEmcalJet *Jet, Int_t JetContNb);
  Double_t                            PTD(AliEmcalJet *Jet, Int_t JetContNb);
  AliEmcalJetFinder                   *Recluster(AliEmcalJet *Jet, Int_t JetContNb, Double_t SubJetRadius, Double_t SubJetMinPt, Int_t Algorithm, const char* Name);
  Double_t                            SubJetOrdering(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Index);
  Double_t                            SubJetFraction(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Add, Bool_t Loss);
  Double_t                            NSubJettiness(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer, Int_t N, Int_t A, Int_t B);
  Double_t                            FjNSubJettiness(AliEmcalJet *Jet, Int_t JetContNb, Int_t N, Int_t Algorithm, Double_t Beta, Int_t Option=0, Double_t Beta_SD=0, Double_t ZCut=0.1);
  Double_t                            FjNSubJettinessFastJet(std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> Jet_ClusterSequence, Int_t JetContNb,Int_t N, Int_t Algorithm, Double_t Beta, Int_t Option=0, Double_t Beta_SD=0, Double_t ZCut=0.1);
  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *>                  ModifyJet(AliEmcalJet* Jet, Int_t JetContNb, TString Modification);
  std::vector<fastjet::PseudoJet>     RandomiseTracks(AliEmcalJet *Jet,std::vector<fastjet::PseudoJet> fInputVectors);
  std::vector<fastjet::PseudoJet>     AddExtraProng(std::vector<fastjet::PseudoJet> fInputVectors, Double_t Distance, Double_t PtFrac);
  std::vector<fastjet::PseudoJet>     AddkTTracks(AliEmcalJet *Jet,std::vector<fastjet::PseudoJet> fInputVectors, Double_t QHat,Double_t Xlength, Int_t NAdditionalTracks);
 
  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Float_t                             fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  JetSelectionType                    fJetSelection;               // Jet selection: inclusive/recoil jet
  TreeSize                            fTreeSize;
  Double_t                            fShapesVar[nVar];                  // jet shapes used for the tagging
  std::vector<std::vector<Double_t>>            fShapesVar_Tracks_Rec;
  std::vector<std::vector<Double_t>>            fShapesVar_Tracks_Truth;
  Float_t                             fPtThreshold;
  Float_t                             fRMatching;

  Float_t                             fPtMinTriggerHadron;                   //min - max pT for trigger particle in case of recoil jet  
  Float_t                             fPtMaxTriggerHadron;
  Float_t                             fRecoilAngularWindow;             //angular window for btb recoil analysis 
  Int_t                               fSemigoodCorrect;             //if==1 we run over semigood runs
  Float_t                             fHolePos;                          //position in radians of the bad TPC sector
  Float_t                             fHoleWidth;                       //width of the hole in radians
  TRandom3                            *fRandom;                     //! Random number generator
  Bool_t                              fCentSelectOn;                // switch on/off centrality selection
  Float_t                             fCentMin;                     // min centrality value
  Float_t                             fCentMax;                     // max centrality value
  Int_t                               fSubJetAlgorithm;
  Float_t                             fSubJetRadius;
  Float_t                             fSubJetMinPt; 
  Double_t                            fJetRadius;
  Double_t                            fRMatched; 
  Double_t                            fSharedFractionPtMin;
  Double_t                            Background_Median;
  Double_t                            Background_Fluc;
  Int_t                               fDerivSubtrOrder;
  Bool_t                              fFullTree;
  Double_t                            fBeta_SD;
  Double_t                            fZCut;
  Int_t                               fReclusteringAlgorithm;
  Int_t                               fSoftDropOn;
  Int_t                               fMLOn;
  Bool_t                              fRandomisationEqualPt;

  Bool_t                              fNsubMeasure;
  TH1F                                *fhPtTriggerHadron;
  TH1F                                *fhJetPt;
  TH1F                                *fhJetPt_1;
  TH1F                                *fhJetPt_2;
  TH1F                                *fhJetPhi;
  TH1F                                *fhJetPhi_1;
  TH1F                                *fhJetPhi_2;
  TH1F                                *fhJetEta;
  TH1F                                *fhJetEta_1;
  TH1F                                *fhJetEta_2;
  TH1F                                *fhJetMass;
  TH1F                                *fhJetMass_1;
  TH1F                                *fhJetMass_2;
  TH1F                                *fhJetRadius;
  TH1F                                *fhJetRadius_1;
  TH1F                                *fhJetRadius_2;
  TH1F                                *fhJetCounter;
  TH1F                                *fhJetCounter_1;
  TH1F                                *fhJetCounter_2;
  TH1F                                *fhNumberOfJetTracks;
  TH1F                                *fhNumberOfJetTracks_1;
  TH1F                                *fhNumberOfJetTracks_2;
  TH1F                                *fhSubJetPt;
  TH1F                                *fhSubJetPt_1;
  TH1F                                *fhSubJetPt_2;
  TH1F                                *fhSubJetMass;
  TH1F                                *fhSubJetMass_1;
  TH1F                                *fhSubJetMass_2;
  TH1F                                *fhSubJetRadius;
  TH1F                                *fhSubJetRadius_1;
  TH1F                                *fhSubJetRadius_2;
  TH1F                                *fhSubJetCounter;
  TH1F                                *fhSubJetCounter_1;
  TH1F                                *fhSubJetCounter_2;
  TH1F                                *fhNumberOfSubJetTracks;
  TH1F                                *fhNumberOfSubJetTracks_1;
  TH1F                                *fhNumberOfSubJetTracks_2;
  TH1F                                *fhEventCounter;
  TH1F                                *fhEventCounter_1;
  TH1F                                *fhEventCounter_2;
  TH1F                                *fhPhiTriggerHadronJet;
  TH1F                                *fhPhiTriggerHadronEventPlane;
  TH1F                                *fhTrackPhi;
  TH1F                                *fhTrackPhi_Cut;
  TH1F                                *fhPhiTriggerHadronEventPlaneTPC;
  TH2F                                *fh2PtTriggerHadronJet;
  TH2F                                *fh2PtRatio; 
  TH2D                                *fhSubJettiness1CheckRatio_FJ_AKT;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_KT;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_CA;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_WTA_KT;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_WTA_CA;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_OP_AKT;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_OP_KT;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_OP_CA;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_OP_WTA_KT;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_OP_WTA_CA;
  TH2D                                *fhSubJettiness1CheckRatio_FJ_MIN;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_AKT;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_KT;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_CA;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_WTA_KT;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_WTA_CA;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_OP_AKT;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_OP_KT;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_OP_CA;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_OP_WTA_KT;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_OP_WTA_CA;
  TH2D                                *fhSubJettiness2CheckRatio_FJ_MIN;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_AKT;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_KT;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_CA;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_WTA_KT;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_WTA_CA;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_OP_AKT;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_OP_KT;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_OP_CA;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_OP_WTA_KT;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_OP_WTA_CA;
  TH2D                                *fhSubJettiness2to1CheckRatio_FJ_MIN;
  TH1D                                *fhSubJettiness1_FJ_AKT;
  TH1D                                *fhSubJettiness1_FJ_KT;
  TH1D                                *fhSubJettiness1_FJ_CA;
  TH1D                                *fhSubJettiness1_FJ_WTA_KT;
  TH1D                                *fhSubJettiness1_FJ_WTA_CA;
  TH1D                                *fhSubJettiness1_FJ_OP_AKT;
  TH1D                                *fhSubJettiness1_FJ_OP_KT;
  TH1D                                *fhSubJettiness1_FJ_OP_CA;
  TH1D                                *fhSubJettiness1_FJ_OP_WTA_KT;
  TH1D                                *fhSubJettiness1_FJ_OP_WTA_CA;
  TH1D                                *fhSubJettiness1_FJ_MIN;
  TH1D                                *fhSubJettiness2_FJ_AKT;
  TH1D                                *fhSubJettiness2_FJ_KT;
  TH1D                                *fhSubJettiness2_FJ_CA;
  TH1D                                *fhSubJettiness2_FJ_WTA_KT;
  TH1D                                *fhSubJettiness2_FJ_WTA_CA;
  TH1D                                *fhSubJettiness2_FJ_OP_AKT;
  TH1D                                *fhSubJettiness2_FJ_OP_KT;
  TH1D                                *fhSubJettiness2_FJ_OP_CA;
  TH1D                                *fhSubJettiness2_FJ_OP_WTA_KT;
  TH1D                                *fhSubJettiness2_FJ_OP_WTA_CA;
  TH1D                                *fhSubJettiness2_FJ_MIN;
  TH1D                                *fhSubJettiness2to1_FJ_AKT;
  TH1D                                *fhSubJettiness2to1_FJ_KT;
  TH1D                                *fhSubJettiness2to1_FJ_CA;
  TH1D                                *fhSubJettiness2to1_FJ_WTA_KT;
  TH1D                                *fhSubJettiness2to1_FJ_WTA_CA;
  TH1D                                *fhSubJettiness2to1_FJ_OP_AKT;
  TH1D                                *fhSubJettiness2to1_FJ_OP_KT;
  TH1D                                *fhSubJettiness2to1_FJ_OP_CA;
  TH1D                                *fhSubJettiness2to1_FJ_OP_WTA_KT;
  TH1D                                *fhSubJettiness2to1_FJ_OP_WTA_CA;
  TH1D                                *fhSubJettiness2to1_FJ_MIN; 
  TTree                               *fTreeResponseMatrixAxis;  //Tree with tagging variables subtracted MC or true MC or raw
  TTree                               *fTreeTracks;

 private:
  AliAnalysisTaskSubJetFraction(const AliAnalysisTaskSubJetFraction&);            // not implemented
  AliAnalysisTaskSubJetFraction &operator=(const AliAnalysisTaskSubJetFraction&); // not implemented

  ClassDef(AliAnalysisTaskSubJetFraction, 7)
};
#endif
