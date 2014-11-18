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

#include "AliAnalysisTaskEmcalJet.h"



class AliAnalysisTaskEmcalQGTagging : public AliAnalysisTaskEmcalJet {
 public:
  enum JetShapeType {
    kTrue = 0,   // generated jets only 
    kTrueDet =1,  // detector and generated jets  
    kData   = 2,  // raw data 
    kDetEmb = 3,  //detector embedded jets 
  };

  enum JetShapeSub {
    kNoSub = 0, 
    kConstSub = 1,
    kDerivSub = 2 
  } ;

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
  void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ; }
  void SetRMatching(Float_t f)                              { fRMatching = f ;}

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
  Float_t                            GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            LeSub(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb);
  Float_t                            GetSigma2(AliEmcalJet *jet, Int_t jetContNb);
  Float_t                            Sigma2(AliEmcalJet *jet, Int_t jetContNb);


  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Float_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  Float_t                            *fShapesVar;                  // jet shapes used for the tagging
  Float_t                             fPtThreshold;
  Float_t                             fRMatching;
  
  TH2F                                *fPhiJetCorr6;
  TH2F                                *fPhiJetCorr7;
  TH2F                                *fEtaJetCorr6;
  TH2F                                *fEtaJetCorr7;
  TH2F                                *fPtJetCorr;
  TH1F                                *fPtJet;
  

  TTree           *fTreeObservableTagging;  //Tree with tagging variables subtracted MC or true MC or raw 

 private:
  AliAnalysisTaskEmcalQGTagging(const AliAnalysisTaskEmcalQGTagging&);            // not implemented
  AliAnalysisTaskEmcalQGTagging &operator=(const AliAnalysisTaskEmcalQGTagging&); // not implemented

  ClassDef(AliAnalysisTaskEmcalQGTagging, 3)
};
#endif

