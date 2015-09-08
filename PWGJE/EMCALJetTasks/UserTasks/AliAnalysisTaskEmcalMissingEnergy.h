#ifndef ALIANALYSISTASKEMCALMISSINGENERGY_H
#define ALIANALYSISTASKEMCALMISSINGENERGY_H

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

#include "AliAnalysisTaskEmcalJet.h"
#include "AliEmcalJetFinder.h"
#include <vector>

class AliAnalysisTaskEmcalMissingEnergy : public AliAnalysisTaskEmcalJet {
 public:
  
  enum JetShapeType {
    kTrue = 0,   // generated jets only 
    kTrueDet =1,  // detector and generated jets  
    kData   = 2,  // raw data 
    kDetEmb = 3,  //detector embedded jets
    kDetEmbPart = 4,
    kPythiaDef = 5
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

  AliAnalysisTaskEmcalMissingEnergy();
  AliAnalysisTaskEmcalMissingEnergy(const char *name);
  virtual ~AliAnalysisTaskEmcalMissingEnergy();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)                             { fContainer     = c   ; }
  void SetMinFractionShared(Double_t f)                     { fMinFractionShared = f   ; }
  void SetJetShapeType(JetShapeType t)                      { fJetShapeType       = t   ; }
  void SetJetShapeSub(JetShapeSub t)                        { fJetShapeSub     = t   ; }
  void SetJetSelection(JetSelectionType t)                  { fJetSelection    = t   ; }
  void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ; }
  void SetJetPtMin(Float_t f)                               { fJetPtMin     = f   ; }
  void SetTrackPtMin(Float_t f)                             { fTrackPtMin     = f   ; }
  void SetRMatching(Float_t f)                              { fRMatching = f ;}
  void SetJetRadius(Float_t f)                              { fJetRadius = f;}
  void SetSubJetRadius(Float_t f)                           { fSubJetRadius = f;}
  void SetPtTriggerSelections(Float_t minpT, Float_t maxpT) { fminpTTrig = minpT; fmaxpTTrig = maxpT; }
  void SetAngularWindowRecoilJet (Float_t t)                {fangWindowRecoil = t; }
  Float_t GetMinPtTriggerSelection()                        {return fminpTTrig;}
  Float_t GetMaxPtTriggerSelection()                        {return fmaxpTTrig;}
  void SetCentralitySelectionOn(Bool_t t)                   { fCentSelectOn = t;}
  void SetHTon(Bool_t t)                                    { fHTon = t;}
  void SetMinCentrality(Float_t t)                          { fCentMin = t ; }
  void SetMaxCentrality(Float_t t)                          { fCentMax = t ; }
  void SetSemigoodCorrect(Int_t yesno)                      {fSemigoodCorrect=yesno;}
  void SetHolePos(Float_t poshole)                          { fHolePos = poshole; }
  void SetHoleWidth(Float_t holewidth)                      { fHoleWidth = holewidth; }
 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Float_t                            GetJetMass(AliEmcalJet *jet,Int_t jetContNb);
  Float_t                            Angularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            PTD(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetpTD(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            Circularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            TauDen(AliEmcalJet *mainJet, Int_t jetContNb);
  Float_t                            Tau1Num(AliEmcalJet *jet, AliEmcalJet *subJet1hardest, Int_t jetContNb);
  Float_t                            Tau1Num_full(AliEmcalJet *mainJet, AliEmcalJetFinder *finder, Int_t jetContNb);
  Float_t                            Tau2Num(AliEmcalJet *jet, AliEmcalJet *subJet1hardest, AliEmcalJet *subJet2hardest, Int_t jetContNb);
  Float_t                            Tau3Num(AliEmcalJet *jet, AliEmcalJet *subJet1hardest, AliEmcalJet *subJet2hardest, AliEmcalJet *subJet3hardest, Int_t jetContNb);
  Int_t *                            JetHard(AliEmcalJetFinder *finder);
  Float_t*                           N_subjettiness(AliEmcalJet *mainJet, Int_t jetContNb);
  Float_t                            GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            LeSub(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb);
  Float_t                            GetSigma2(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            Sigma2(AliEmcalJet *jet, Int_t jetContNb);

  Int_t                              SelectTrigger(Float_t minpT, Float_t maxpT);
  Double_t                           RelativePhi(Double_t mphi, Double_t vphi);
  Double_t                           RelativePhiFancy(Double_t mphi, Double_t vphi);

  Double_t                           Minimum(Double_t x, Double_t y);
  Double_t                           Minimum(Double_t x, Double_t y, Double_t z);
  Double_t                           R_distance(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2);
  Double_t                           R_distance(AliEmcalJet *jet1, AliEmcalJet *jet2);
  Double_t                           R_distance(AliEmcalJet *jet,  AliPicoTrack *part);

  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Float_t                             fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  JetSelectionType                    fJetSelection;               // Jet selection: inclusive/recoil jet  

  Float_t                            *fSubstructureVar;                  // my observables for the substructure tree
  Float_t                            *fHadronTriggerVar;                  // my observables for the substructure tree

  Float_t                             fPtThreshold;                   
  Float_t                             fJetPtMin;                          // Jet Pt min under consideration
  Float_t                             fTrackPtMin;                        // track pt minimum as input to reclustering
  Float_t                             fRMatching;

  Float_t                             fJetRadius;                  // radius of the main jet finding
  Float_t                             fSubJetRadius;               // radius of the sub jet finding


  Float_t                             fminpTTrig;                   //min - max pT for trigger particle in case of recoil jet  
  Float_t                             fmaxpTTrig;
  Float_t                             fangWindowRecoil;             //angular window for btb recoil analysis 
  Int_t                               fSemigoodCorrect;             //if==1 we run over semigood runs
  Float_t                             fHolePos;                          //position in radians of the bad TPC sector
  Float_t                             fHoleWidth;                       //width of the hole in radians 
  Bool_t                              fCentSelectOn;                // switch on/off centrality selection
  Bool_t                              fHTon;                       // switch on/off Hadron Trigger analysis
  Float_t                             fCentMin;                     // min centrality value
  Float_t                             fCentMax;                     // max centrality value
  
  TH2F                                *fh2ResponseUW;
  TH2F                                *fh2ResponseW;
  TH2F                                *fPhiJetCorr6;
  TH2F                                *fPhiJetCorr7;
  TH2F                                *fEtaJetCorr6;
  TH2F                                *fEtaJetCorr7;
  TH2F                                *fPtJetCorr;
  TH1F                                *fPtJet;
  TH2F                                *fhpTjetpT; //control plot fo the recoil analysis
  TH1F                                *fhPt;
  TH1F                                *fhCentrality;
  TH1F                                *fhPhi;
  TH1F                                *fhJetPt;
  TH1F                                *fhJetPhi;
  TH1F                                *fhTrackPt;
  TH1F                                *fhTrackPt_JT;
  TH1F                                *fhTriggerPt;
  TH1F                                *fhJetTriggeredPt;
  TH2F                                *fhTriggerVSjetPt;
  TH2F                                *fhdPhiTrigVSjetPt;
  TH2F                                *fhdPhiTrigVSjetPtNOCUT;
  TH2F                                *fhdPhiPartVSjetPt;
  TH3F                                *fhdPhiPartVSjetPtVSPartPt;
  TH1F                                *fhTau1;
  TH1F                                *fhTau2;
  TH1F                                *fhTau3;
  TH2F                                *fhTau1vsTau2;
  TH3F                                *fhTau1vsTau2vsTau3;
  TH2F                                *fhTau1vsTau3;
  TH2F                                *fhTau2vsTau3;
  TH1F                                *fhTau2OverTau1;
  TH1F                                *fhTau3OverTau2;
  TH2F                                *fhTau2vsJetPt;
  TH2F                                *fhTau2OverTau1vsJetPt;


  TTree                               *fSubstructure;
  TTree                               *fHadronTrigger;

 private:
  AliAnalysisTaskEmcalMissingEnergy(const AliAnalysisTaskEmcalMissingEnergy&);            // not implemented
  AliAnalysisTaskEmcalMissingEnergy &operator=(const AliAnalysisTaskEmcalMissingEnergy&); // not implemented

  ClassDef(AliAnalysisTaskEmcalMissingEnergy, 4)
    };
#endif

