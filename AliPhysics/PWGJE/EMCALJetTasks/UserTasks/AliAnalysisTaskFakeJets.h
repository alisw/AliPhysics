#ifndef ALIANALYSISTASKFAKEJETS_H
#define ALIANALYSISTASKFAKEJETS_H

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



class AliAnalysisTaskFakeJets : public AliAnalysisTaskEmcalJet {
 public:
  
  enum JetShapeType {
    kTrue = 0,   // generated jets only 
    kTrueDet =1,  // detector and generated jets  
    kData   = 2,  // raw data 
    kDetEmb = 3,  //detector embedded jets
    kDetEmbPart=4,
    kPythiaDef = 5,
    kDetEmbPartPythia=6
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

  AliAnalysisTaskFakeJets();
  AliAnalysisTaskFakeJets(const char *name);
  virtual ~AliAnalysisTaskFakeJets();

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
  void SetPtTriggerSelections(Float_t minpT, Float_t maxpT) { fminpTTrig = minpT; fmaxpTTrig = maxpT; }
  void SetAngularWindowRecoilJet (Float_t t)               {fangWindowRecoil = t; }
  void SetSubjetRadius(Float_t t)                            {fSubjetRadius = t;}
 void SetJetRadius(Float_t t)                            {fJetRadius = t;}
   Float_t GetMinPtTriggerSelection()                        {return fminpTTrig;}
  Float_t GetMaxPtTriggerSelection()                        {return fmaxpTTrig;}
  void SetCentralitySelectionOn(Bool_t t)                   { fCentSelectOn = t;}
  void SetOneConstSelectionOn(Bool_t t)                     { fOneConstSelectOn =t;}
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
   Float_t                            AngularitySquared(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetAngularitySquared(AliEmcalJet *jet, Int_t jetContNb);

  Float_t                            PTD(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetpTD(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            Circularity(AliEmcalJet *jet, Int_t jetContNb); 
  Float_t                            GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            LeSub(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb);
  Float_t                            GetJetHardTrack(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            HardTrack(AliEmcalJet *jet, Int_t jetContNb);
 Float_t                            GetJetSecHardTrack(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            SecHardTrack(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetCoreFrac(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            CoreFrac(AliEmcalJet *jet, Int_t jetContNb);  
  AliEmcalJetFinder                   *Recluster(AliEmcalJet *Jet, Int_t JetContNb, Double_t SubJetRadius, Double_t SubJetMinPt, Int_t Algorithm, const char* Name);
Double_t                            SubJetOrdering(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Index);
  Double_t                            NSubJettiness(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer, Int_t N, Int_t A, Int_t B);
   Double_t                           GetSubjetFraction(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer);
  
   Double_t                           SubjetFraction(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer);
  Int_t                              SelectTrigger(Float_t minpT, Float_t maxpT);
  Double_t                           RelativePhi(Double_t mphi, Double_t vphi);

  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Float_t                             fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  JetSelectionType                    fJetSelection;               // Jet selection: inclusive/recoil jet  
  Float_t                            *fShapesVar;                  // jet shapes used for the tagging
  Float_t                             fPtThreshold;
  Float_t                             fRMatching;

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
  Int_t                                fDerivSubtrOrder;
  Float_t                             fSubjetRadius;
  Float_t                              fJetRadius;
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
  TH2F                                *fNbOfConstvspT;

  TTree           *fTreeFakeJets;  //Tree with tagging variables subtracted MC or true MC or raw 

 private:
  AliAnalysisTaskFakeJets(const AliAnalysisTaskFakeJets&);            // not implemented
  AliAnalysisTaskFakeJets &operator=(const AliAnalysisTaskFakeJets&); // not implemented

  ClassDef(AliAnalysisTaskFakeJets, 6)
};
#endif

