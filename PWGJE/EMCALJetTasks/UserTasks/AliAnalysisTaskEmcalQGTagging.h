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
    kRaw   = 0,  //mass form anti-kt 4-vector
    kConstSub = 1,   //constituent subtracted jetshape
    kTrue  = 2, 
    kDeriv = 3 //area based subtracted jet mass
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
  void SetMCTask(Int_t f)                                   { fIsMC       = f   ; }
  void SetEmbeddingTask(Int_t f)                            { fIsEmbedding       = f   ; }
  void SetIsConstSub(Bool_t f)                              { fIsConstSub     = f   ; }
  void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ; }
  void SetRMatching(Float_t f)                              { fRMatching = f ;}

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Float_t                            GetJetMass(AliEmcalJet *jet);
  Float_t                            Angularity(AliEmcalJet *jet);
  Float_t                            GetJetAngularity(AliEmcalJet *jet);
  Float_t                            PTD(AliEmcalJet *jet);
  Float_t                            GetJetpTD(AliEmcalJet *jet);
  Float_t                            Circularity(AliEmcalJet *jet); 
  Float_t                            GetJetCircularity(AliEmcalJet *jet);
  Float_t                            LeSub(AliEmcalJet *jet);
  Float_t                            GetJetLeSub(AliEmcalJet *jet);
  Float_t                            GetJetNumberOfConstituents(AliEmcalJet *jet);
  Float_t                            GetSigma2(AliEmcalJet *jet);
  Float_t                            Sigma2(AliEmcalJet *jet);


  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Float_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;                // jet mass type to be used
  Float_t                            *fShapesVar;                  // jet shapes used for the tagging
  Int_t                               fIsMC;
  Int_t                               fIsEmbedding;
  Bool_t                              fIsConstSub;
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

  ClassDef(AliAnalysisTaskEmcalQGTagging, 1)
};
#endif

