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

#include "AliAnalysisTaskEmcalJet.h"



class AliAnalysisTaskSubJetFraction : public AliAnalysisTaskEmcalJet {
 public:
  
  enum JetShapeType {
    kTrue = 0,   // generated jets only 
    kTrueDet =1,  // detector and generated jets  
    kData   = 2,  // raw data 
    kDetEmb = 3  //detector embedded jets 
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

  AliAnalysisTaskSubJetFraction();
  AliAnalysisTaskSubJetFraction(const char *name);
  virtual ~AliAnalysisTaskSubJetFraction();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)                             { fContainer     = c   ; }
  void SetMinFractionShared(Double_t f)                     { fMinFractionShared = f   ; }
  void SetJetShapeType(JetShapeType t)                      { fJetShapeType       = t   ; }
  void SetJetShapeSub(JetShapeSub t)                        { fJetShapeSub     = t   ; }
  void SetJetSelection(JetSelectionType t)                  { fJetSelection    = t   ; }
  void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ; }
  void SetJetRadius(Double_t JetRadius)                      {fJetRadius=JetRadius;}
  void SetRMatching(Float_t f)                              { fRMatching = f ;}
  void SetPtTriggerSelections(Float_t minpT, Float_t maxpT) { fminpTTrig = minpT; fmaxpTTrig = maxpT; }
  void SetAngularWindowRecoilJet (Float_t t)                {fangWindowRecoil = t; }
  Float_t GetMinPtTriggerSelection()                        {return fminpTTrig;}
  Float_t GetMaxPtTriggerSelection()                        {return fmaxpTTrig;}
  void SetCentralitySelectionOn(Bool_t t)                   { fCentSelectOn = t;}
  void SetMinCentrality(Float_t t)                          { fCentMin = t ; }
  void SetMaxCentrality(Float_t t)                          { fCentMax = t ; }
  void SetSemigoodCorrect(Int_t yesno)                 {fSemigoodCorrect=yesno;}
  void SetHolePos(Float_t poshole)                        { fHolePos = poshole; }
  void SetHoleWidth(Float_t holewidth)                  { fHoleWidth = holewidth; }
  void SetSubJetAlgorithm(Int_t SubJetAlgorithm)        {fSubJetAlgorithm=SubJetAlgorithm;}
  void SetSubJetRadius(Float_t SubJetRadius)            {fSubJetRadius=SubJetRadius;}
  void SetSubJetMinPt(Float_t SubJetMinPt)              {fSubJetMinPt=SubJetMinPt;}
 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Int_t                              SelectTrigger(Float_t minpT, Float_t maxpT);
  // Double_t                           RelativePhi(Double_t mphi, Double_t vphi);
  Double_t                           RelativePhi(Double_t Phi1, Double_t Phi2);

  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Float_t                             fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  JetSelectionType                    fJetSelection;               // Jet selection: inclusive/recoil jet  
  Double_t                            *fShapesVar;                  // jet shapes used for the tagging
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
  Int_t                               fSubJetAlgorithm;
  Float_t                             fSubJetRadius;
  Float_t                             fSubJetMinPt; 
  Double_t                             fJetRadius;

  Double_t                              Background_Median;
  Double_t                              Background_Fluc;
  

  TH1F                                *fhJetPt;
  TH1F                                *fhJetPhi;
  TH1F                                *fhJetEta;
  TH1F                                *fhJetMass;
  TH1F                                *fhJetRadius;
  TH1F                                *fhJetAngularity;
  TH1F                                *fhJetPTD;
  TH1F                                *fhNumberOfJetTracks;
  TH1F                                *fhSubJetPt;
  TH1F                                *fhSubJetRadius;
  TH1F                                *fhSubJetPtFrac;
  TH1F                                *fhSubJetPtFrac2;
  TH1F                                *fhSubJetPtLoss;
  TH1F                                *fhSubJetPtLoss2;
  TH1F                                *fhSubJetEnergyFrac;
  TH1F                                *fhSubJetEnergyFrac2;
  TH1F                                *fhSubJetEnergyLoss;
  TH1F                                *fhSubJetEnergyLoss2;
  TH1F                                *fhEventCounter;  
  TH1F                                *fhJetCounter;
  TH1F                                *fhSubJetCounter;
  TH1F                                *fhPtRatio;
  TH1F                                *fhParticleSubJetPtFrac;
  TH1F                                *fhDetectorSubJetPtFrac;
  TH1F                                *fhParticleSubJetPtFrac2;
  TH1F                                *fhDetectorSubJetPtFrac2;
  TH1F                                *fhSubJetPtFracRatio;
  TH1F                                *fhSubJetPtFrac2Ratio;
  TH1F                                *fhParticleSubJetiness1;
  TH1F                                *fhDetectorSubJetiness1;
  TH1F                                *fhSubJetiness1Ratio;
  TH1F                                *fhSubJetiness1;
  TH1F                                *fhParticleSubJetiness2;
  TH1F                                *fhDetectorSubJetiness2;
  TH1F                                *fhSubJetiness2Ratio;
  TH1F                                *fhSubJetiness2;
  TH1F                                *fh2to1SubJetinessRatio;
  TTree                               *fTreeResponseMatrixAxis;  //Tree with tagging variables subtracted MC or true MC or raw 

 private:
  AliAnalysisTaskSubJetFraction(const AliAnalysisTaskSubJetFraction&);            // not implemented
  AliAnalysisTaskSubJetFraction &operator=(const AliAnalysisTaskSubJetFraction&); // not implemented

  ClassDef(AliAnalysisTaskSubJetFraction, 4)
};
#endif

