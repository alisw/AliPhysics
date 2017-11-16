#ifndef ALIANALYSISTASKRECURSIVESOFTDROP_H
#define ALIANALYSISTASKRECURSIVESOFTDROP_H

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
#include "FJ_includes.h"

class AliAnalysisTaskRecursiveSoftDrop : public AliAnalysisTaskEmcalJet {
 public:
  

  enum JetShapeSub {
    kNoSub = 0, 
    kConstSub = 1
  };

  AliAnalysisTaskRecursiveSoftDrop();
  AliAnalysisTaskRecursiveSoftDrop(const char *name);
  virtual ~AliAnalysisTaskRecursiveSoftDrop();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)                             { fContainer     = c   ;}
  void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ;}
  void SetCentralitySelectionOn(Bool_t t)                   { fCentSelectOn = t;}
  void SetMinCentrality(Float_t t)                          { fCentMin = t ;}
  void SetMaxCentrality(Float_t t)                          { fCentMax = t ;}
  void SetJetShapeSub(JetShapeSub t)                        { fJetShapeSub     = t   ;}
  

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  void                                RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont, Int_t ReclusterAlgo,Bool_t bTruth);

  
  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Double_t                            fShapesVar_Det[4];                  // jet shapes used for the tagging
  Double_t                            fShapesVar_True[4];                  // jet shapes used for the tagging
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  Float_t                             fPtThreshold;
  
  Bool_t                              fCentSelectOn;                // switch on/off centrality selection
  Float_t                             fCentMin;                     // min centrality value
  Float_t                             fCentMax;                     // max centrality value
  Double_t                            fJetRadius;

  
  TH1F                                *fhJetPt;
  TH1F                                *fhJetPhi;
  TH1F                                *fhJetEta;
 
  TTree                               *fTreeRecursive_Det;  
  TTree                               *fTreeRecursive_True; 

 private:
  AliAnalysisTaskRecursiveSoftDrop(const AliAnalysisTaskRecursiveSoftDrop&);            // not implemented
  AliAnalysisTaskRecursiveSoftDrop &operator=(const AliAnalysisTaskRecursiveSoftDrop&); // not implemented

  ClassDef(AliAnalysisTaskRecursiveSoftDrop, 1)
};
#endif

