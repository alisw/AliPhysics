#ifndef ALIANALYSISTASKEMCALJETAXISRES_H
#define ALIANALYSISTASKEMCALJETAXISRES_H

#include "AliAnalysisTaskEmcalJet.h"
#include "AliFJWrapper.h"
#include "AliEventPoolManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"

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
class AliEventPoolManager;
class AliQnCorrectionsManager;

class AliAnalysisTaskEmcalJetAxisRes : public AliAnalysisTaskEmcalJet {
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


  AliAnalysisTaskEmcalJetAxisRes();
  AliAnalysisTaskEmcalJetAxisRes(const char *name);
  virtual ~AliAnalysisTaskEmcalJetAxisRes();

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
  void SetMinCentrality(Float_t t)                          { fCentMin = t ; }
  void SetMaxCentrality(Float_t t)                          { fCentMax = t ; }
  void SetSemigoodCorrect(Int_t yesno)                 {fSemigoodCorrect=yesno;}
  void SetHolePos(Float_t poshole)                        { fHolePos = poshole; }
  void SetHoleWidth(Float_t holewidth)                  { fHoleWidth = holewidth; }
  void SetDerivativeSubtractionOrder(Int_t c)              {fDerivSubtrOrder = c;}

  void SetMinPtForJetAxisRecalc(Double_t pt=4.)               { fMinPtForJetRecalc = pt ; }
  void SetIsUseEventSelDPG     (Bool_t isUseDPG=kTRUE)        { fIsEventSelDPG = isUseDPG ; }
  void SetEPName               (const char *epname)           { fEPName = epname ; }           
  void SetMinPartJetPt         (Double_t pt=10.)              { fMinJetPtPart = pt ; }           
  void SetMinDetJetPt          (Double_t pt=10.)              { fMinJetPtDet  = pt ; }           
  TObjArray*                   CloneAndReduceTrackList();
 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  Bool_t                              EventSelectionDPG();

  Float_t                            GetJetMass(AliEmcalJet *jet,Int_t jetContNb=0);
  Float_t                            Angularity(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            Coronna(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetCoronna(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            PTD(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetpTD(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            Circularity(AliEmcalJet *jet, Int_t jetContNb=0); 
  Float_t                            GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            LeSub(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb=0);
  Float_t                            GetSigma2(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            Sigma2(AliEmcalJet *jet, Int_t jetContNb=0);

  Int_t                              SelectTrigger(Float_t minpT, Float_t maxpT);
  Double_t                           RelativePhi(Double_t mphi, Double_t vphi);
   void                                RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont, Int_t ReclusterAlgo);
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
  Int_t                               fDerivSubtrOrder;

  AliQnCorrectionsManager             *fFlowQnVectorMgr;  
  AliEventCuts                        fEventCuts;
  Bool_t                              fIsEventSelDPG;
  Double_t                            fMinJetPtPart;
  Double_t                            fMinJetPtDet;
  Double_t                            fMinPtForJetRecalc;
  Double_t                            fEP;
  TString                             fEPName;
  Double_t                            fDivEvPlane[9];
  Double_t                            fDivPt[6];
  AliEventPoolManager                 *fPoolMgr;                 //!<! Event pool manager
  
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
  THnSparse                           *fHLundIterative;// 
  TH2F                                *fNbOfConstvspT;

  TH2D                                *fh2JetdEtadPhiEmbVsDet [9];
  TH2D                                *fh2JetdEtadPhiEmbVsPart[9];
  TH2D                                *fh2JetdEtadPhiDetVsPart;

  TH2D                                *fh2JetdEtadPhiRecalcEmbVsRecalcDet [9];
  TH2D                                *fh2JetdEtadPhiRecalcEmbVsRecalcPart[9];
  TH2D                                *fh2JetdEtadPhiRecalcDetVsRecalcPart;

  TH2D                                *fh2JetdEtadPhiRecalcEmbVsDet [9];
  TH2D                                *fh2JetdEtadPhiRecalcEmbVsPart[9];
  TH2D                                *fh2JetdEtadPhiRecalcDetVsPart;

  TH2D                                *fh2JetHadCorr[9][6];
  TH2D                                *fh2JetHadCorrNear[9][6];
  TH2D                                *fh2JetHadCorrNearSideBand[9][6];

  TH2D                                *fh2JetHadCorrMix[9][6];
  TH2D                                *fh2JetHadCorrMixNear[9][6];
  TH2D                                *fh2JetHadCorrMixNearSideBand[9][6];


  TTree           *fTreeObservableTagging;  // Tree with tagging variables subtracted MC or true MC or raw 

 private:
  AliAnalysisTaskEmcalJetAxisRes(const AliAnalysisTaskEmcalJetAxisRes&);            // not implemented
  AliAnalysisTaskEmcalJetAxisRes &operator=(const AliAnalysisTaskEmcalJetAxisRes&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetAxisRes, 6)
};
#endif

