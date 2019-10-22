#ifndef ALIANALYSISTASKEMCALQGTAGGING_H
#define ALIANALYSISTASKEMCALQGTAGGING_H

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


class AliAnalysisTaskEmcalQGTagging : public AliAnalysisTaskEmcalJet {
 public:
  
  enum JetShapeType {
    kMCTrue = 0,   // generated jets only
    kTrueDet =1,  // detector and generated jets  
    kData   = 2,  // raw data 
    kDetEmb = 3,  //detector embedded jets
    kDetEmbPart = 4,
    kPythiaDef = 5,
    kDetEmbPartPythia=6,
    kGenOnTheFly = 7
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

  AliAnalysisTaskEmcalQGTagging();
  AliAnalysisTaskEmcalQGTagging(const char *name);
  virtual ~AliAnalysisTaskEmcalQGTagging();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)                             { fContainer     = c   ; }
  void SetMinFractionShared(Double_t f)                     { fMinFractionShared = f   ; }
  void SetJetShapeType(JetShapeType t)                      { fJetShapeType       = t   ; }
  void SetJetShapeSub(JetShapeSub t)                        { fJetShapeSub     = t   ; }
  void SetJetSelection(JetSelectionType t)                  { fJetSelection    = t   ; }
  void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ; }
  void SetRMatching(Float_t f)                              { fRMatching = f ;}
  void SetSelectShapes(Int_t c)                                {fSelectedShapes = c;}
  void SetPtTriggerSelections(Float_t minpT, Float_t maxpT) { fminpTTrig = minpT; fmaxpTTrig = maxpT; }
  void SetAngularWindowRecoilJet (Float_t t)                {fangWindowRecoil = t; }
  Float_t GetMinPtTriggerSelection()                        {return fminpTTrig;}
  Float_t GetMaxPtTriggerSelection()                        {return fmaxpTTrig;}
  void SetCentralitySelectionOn(Bool_t t)                   { fCentSelectOn = t;}
  void SetOneConstSelectionOn(Bool_t t)                     { fOneConstSelectOn =t;}
   void SetCheckTracksOn(Bool_t t)                         { fTrackCheckPlots =t;}
   void SetCheckResolution(Bool_t t)                       {fCheckResolution = t;} 
   void SetSubjetCutoff(Float_t t)                            {fSubjetCutoff = t;}
   void SetMinPtConst(Float_t t)                              {fMinPtConst = t;}
   void SetHardCutoff(Float_t t)                            {fHardCutoff = t;}
   void SetDoTwoTrack(Bool_t t)                             {fDoTwoTrack = t;}
   void SetDoAreaIterative(Bool_t t)                        {fDoAreaIterative =t;}
   void SetMagFieldPol(Float_t t)                           {fMagFieldPolarity=t;}
  void SetMinCentrality(Float_t t)                          { fCentMin = t ; }
  void SetMaxCentrality(Float_t t)                          { fCentMax = t ; }
  void SetSemigoodCorrect(Int_t yesno)                 {fSemigoodCorrect=yesno;}
  void SetHolePos(Float_t poshole)                        { fHolePos = poshole; }
  void SetHoleWidth(Float_t holewidth)                  { fHoleWidth = holewidth; }
  void SetDerivativeSubtractionOrder(Int_t c)              {fDerivSubtrOrder = c;}
 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Float_t                            GetJetMass(AliEmcalJet *jet,Int_t jetContNb);
  Float_t                            Angularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            Coronna(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetCoronna(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            PTD(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetpTD(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            Circularity(AliEmcalJet *jet, Int_t jetContNb); 
  Float_t                            GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            LeSub(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb);
  Float_t                            GetSigma2(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            Sigma2(AliEmcalJet *jet, Int_t jetContNb);

  Int_t                              SelectTrigger(Float_t minpT, Float_t maxpT);
  Double_t                           RelativePhi(Double_t mphi, Double_t vphi);
  void                                RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont);
  void                                RecursiveParentsMCAverage(AliEmcalJet *fJet,Int_t km, Double_t &aver1, Double_t &aver2, Double_t &aver3, Double_t &aver4);
  void                               CheckSubjetResolution(AliEmcalJet *fJet,AliJetContainer *fJetCont, AliEmcalJet *fJetM,AliJetContainer *fJetContM);
  Bool_t                              CheckClosePartner(Int_t index, AliEmcalJet *fJet,AliVParticle *fTrack, AliParticleContainer *fTrackCont);
  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Float_t                             fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  JetSelectionType                    fJetSelection;               // Jet selection: inclusive/recoil jet  
  Float_t                             fShapesVar[12];                  // jet shapes used for the tagging
  Float_t                             fPtThreshold;
  Float_t                             fRMatching;
  Int_t                                 fSelectedShapes;                //chose set of shapes 
  Float_t                             fminpTTrig;                   //min - max pT for trigger particle in case of recoil jet  
  Float_t                             fmaxpTTrig;
  Float_t                             fangWindowRecoil;             //angular window for btb recoil analysis 
  Int_t                                fSemigoodCorrect;             //if==1 we run over semigood runs
  Float_t                             fHolePos;                          //position in radians of the bad TPC sector
  Float_t                             fHoleWidth;                       //width of the hole in radians 
  Bool_t                              fCentSelectOn;                // switch on/off centrality selection
  Float_t                             fCentMin;                     // min centrality value
  Float_t                             fCentMax;                     // max centrality value
  Bool_t                              fOneConstSelectOn;                // switch on/off one constituent selection
  Bool_t                              fTrackCheckPlots;              //switch on qa plots
  Bool_t                              fCheckResolution;              //check subjet energy resolution 
  Float_t                             fSubjetCutoff;                 //angular cutoff for subjets at det/gen level
  Float_t                             fMinPtConst;                   //constituent pt cutoff   
  Float_t                             fHardCutoff;                   //hard cutoff in the iterative declustering 
  Bool_t                              fDoTwoTrack;                    //switch to consider 2 track effects
  Bool_t                              fDoAreaIterative;               //subtract the area in the declustering
  Float_t                             fPhiCutValue;                  //cuts from HBT
  Float_t                             fEtaCutValue;                  //cuts from HBT
  Float_t                             fMagFieldPolarity;             //polarity, to calculate phimin 
  Int_t                               fDerivSubtrOrder;

  
  TH2F                                *fh2ResponseUW;
  TH2F                                *fh2ResponseW;
  TH2F                                *fPhiJetCorr6;
  TH2F                                *fPhiJetCorr7;
  TH2F                                *fEtaJetCorr6;
  TH2F                                *fEtaJetCorr7;
  TH2F                                *fPtJetCorr;
  TH1F                                *fPtJet;
  TH2F                                *fhpTjetpT; //control p[lot fo the recoil analysis
  TH1F                                *fhPt;
  TH1F                                *fhPhi;
  TH1F                                *fhTrackPhi;   //control plot for tracks 
  THnSparse                           *fHLundIterative;//       iterative declustering
  THnSparse                           *fHCheckResolutionSubjets;//     to evaluate energy resolution of subjets as function fo apperture angle 
  TH2F                                *fNbOfConstvspT;

  TTree           *fTreeObservableTagging;  // Tree with tagging variables subtracted MC or true MC or raw 

 private:
  AliAnalysisTaskEmcalQGTagging(const AliAnalysisTaskEmcalQGTagging&);            // not implemented
  AliAnalysisTaskEmcalQGTagging &operator=(const AliAnalysisTaskEmcalQGTagging&); // not implemented

  ClassDef(AliAnalysisTaskEmcalQGTagging, 6)
};
#endif

