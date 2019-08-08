#ifndef AliAnalysisTaskEmcalJetShapeExtra_H
#define AliAnalysisTaskEmcalJetShapeExtra_H

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



class AliAnalysisTaskEmcalJetShapeExtra : public AliAnalysisTaskEmcalJet {
 public:


  enum JetShapeType {

    kData   = 2,  // raw data
  
    kPythiaDef = 5,
 
  };
  enum JetShapeSub {
    kNoSub = 0, 
    kConstSub = 1,
    kDerivSub = 2
  };

 enum DerivSubtrOrder {
    kSecondOrder = 0,
    kFirstOrder = 1
  };
    AliAnalysisTaskEmcalJetShapeExtra();
  AliAnalysisTaskEmcalJetShapeExtra(const char *name);
  virtual ~AliAnalysisTaskEmcalJetShapeExtra();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)                             { fContainer     = c   ; }
  
  void SetJetShapeType(JetShapeType t)                      { fJetShapeType       = t   ; }
  void SetJetShapeSub(JetShapeSub t)                        { fJetShapeSub     = t   ; }
  
 void SetDerivativeSubtractionOrder(Int_t c)              {fDerivSubtrOrder = c;}
  
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

 
  Double_t                           RelativePhi(Double_t mphi, Double_t vphi);

  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 

  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  //JetSelectionType                    fJetSelection;               // Jet selection: inclusive/recoil jet
  Float_t                             fShapesVar[18];                  // jet shapes used for the tagging

  Int_t                               fDerivSubtrOrder;

  
  TH2F                                *fh2ResponseUW;

  TH1F                                *fPtJet;
    TH1F                                *fNumberOfJet;
 TH1F                                *fRho;
  TH2F                                *fNbOfConstvspT;

  TTree           *fTreeObservableTagging;  // Tree with tagging variables subtracted MC or true MC or raw 

 private:
  AliAnalysisTaskEmcalJetShapeExtra(const AliAnalysisTaskEmcalJetShapeExtra&);            // not implemented
  AliAnalysisTaskEmcalJetShapeExtra &operator=(const AliAnalysisTaskEmcalJetShapeExtra&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetShapeExtra, 6)
};
#endif

