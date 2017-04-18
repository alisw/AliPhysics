#ifndef ALIANALYSISTASKRECOILJETYIELD_H
#define ALIANALYSISTASKRECOILJETYIELD_H

#define LOG_NO_WARNING

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
const Int_t nBranch = 18;
class AliAnalysisTaskRecoilJetYield : public AliAnalysisTaskEmcalJet {
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

  AliAnalysisTaskRecoilJetYield();
  AliAnalysisTaskRecoilJetYield(const char *name);
  virtual ~AliAnalysisTaskRecoilJetYield();

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
  void SetDoSoftDrop(Bool_t SoftDrop)                       {fDoSoftDrop = SoftDrop;}
  
  void SetNsubUnNormMeasure( Bool_t NsubMeasure)              {fNsubMeasure= NsubMeasure;}
  void SetRhoName(const char *r)                              {fRhoName = r;}
  void SetSoftDropRecluster(Int_t n)                          {fReclusterAlgo = n;} //0 = CA, 1 = anti-kt, 2 = kt
  void SetDoSubDetMatching(Bool_t b)                          {fSubMatching = b;}
  

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Double_t                            RelativePhi(Double_t Phi1, Double_t Phi2);
  Double_t                            RelativePhiEventPlane(Double_t EventPlane, Double_t Phi);

  Double_t                            Angularity(AliEmcalJet *Jet, Int_t JetContNb);
  Double_t                            PTD(AliEmcalJet *Jet, Int_t JetContNb);
  void                                SoftDrop(AliEmcalJet *fJet,AliJetContainer *fJetCont, double zcut, double beta, Bool_t fTruthJet);
  Int_t                               SelectTriggerHadron(Float_t PtMin, Float_t PtMax);
  Double_t                            GetFractionSharedPt_SubMatching(const AliEmcalJet *jet, AliParticleContainer *cont2 = 0x0) const;
  static Bool_t                       SameParticle(const AliVParticle* part1, const AliVParticle* part2, Double_t dist = 1.e-4);


  
  
  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Float_t                             fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  JetSelectionType                    fJetSelection;               // Jet selection: inclusive/recoil jet
  Double_t                            fJetInfoVar[nBranch];                  // jet shapes used for the tagging
  Float_t                             fPtThreshold;
  Float_t                             fRMatching;

  Float_t                             fPtMinTriggerHadron;                   //min - max pT for trigger particle in case of recoil jet  
  Float_t                             fPtMaxTriggerHadron;
  Float_t                             fRecoilAngularWindow;             //angular window for btb recoil analysis 
  Int_t                               fSemigoodCorrect;             //if==1 we run over semigood runs
  Float_t                             fHolePos;                          //position in radians of the bad TPC sector
  Float_t                             fHoleWidth;                       //width of the hole in radians 
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
  Bool_t                              fDoSoftDrop;
  Double_t                            fRho;
  TString                             fRhoName;
  AliRhoParameter                     *fRhoParam;
  Bool_t                              fSubMatching;

  Bool_t                              fNsubMeasure;
  Int_t                               fReclusterAlgo;
  
  TH1F                                *fhJetPt;
  TH1F                                *fhJetPhi;
  TH1F                                *fhJetEta;
  TH1F                                *fhJetMass;
  TH1F                                *fhJetRadius;
  TH1F                                *fhJetCounter;
  TH1F                                *fhNumberOfJetTracks;
  TH1F                                *fhJetArea;
  TH1F                                *fhPtTriggerHadron;
  TH2F                                *fh2PtTriggerHadronJet;
  TH1F                                *fhPhiTriggerHadronJet;
  TH1F                                *fhPhiTriggerHadronEventPlane;
  TH1F                                *fhPhiTriggerHadronEventPlaneTPC;
  TH1F                                *fhTrackPt;
  TH2F                                *fhGroomedPtvJetPt;
  TH1F                                *fhDroppedBranches;
  TH1F                                *fhDetJetPt_Incl;
  TH1F                                *fhDetJetPt_Matched;
  TTree                               *fTreeJetInfo;  //Tree with tagging variables subtracted MC or true MC or raw
  

 private:
  AliAnalysisTaskRecoilJetYield(const AliAnalysisTaskRecoilJetYield&);            // not implemented
  AliAnalysisTaskRecoilJetYield &operator=(const AliAnalysisTaskRecoilJetYield&); // not implemented

  ClassDef(AliAnalysisTaskRecoilJetYield, 6)
};
#endif

