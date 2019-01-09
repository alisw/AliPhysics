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
  enum JetType {
    kData = 0, 
    kEmb = 1,
    kTrueDet = 2
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
  void SetJetType(JetType t)                                { fJetType     = t   ;}
  void SetReclusterAlgo(Int_t a)                            { fReclusteringAlgo = a;}
  void AddMedScat(Bool_t b, Float_t f,Int_t n)              { fAddMedScat = b; fAddMedScatPtFrac = f; fAddMedScatN = n;}
  void DoSubJetAreaSub(Bool_t b)                            { fDoSubJetAreaSub = b;}

  static AliAnalysisTaskRecursiveSoftDrop* AddTaskRecursiveSoftDrop(

							     const char * njetsData, //data jets
							     const char * njetsTrue, //Pythia Particle Level
							     const char * njetsDet,
							     const char * njetsHybridUs,
							     const char * njetsHybridS,
							     const Double_t R,
							     const char * nrhoBase, 
							     const char * ntracksData,
							     const char * ntracksTrue,
							     const char * ntracksDet, 
							     const char * ntracksHybridUs,
							     const char * ntracksHybridS,
							     const char *type,				      
							     const char *CentEst,
							     Int_t       pSel,
							     TString     trigClass      = "",
							     TString     kEmcalTriggers = "",
							     TString     tag            = "",
							     AliAnalysisTaskRecursiveSoftDrop::JetShapeSub jetShapeSub = JetShapeSub::kConstSub,
							     AliAnalysisTaskRecursiveSoftDrop::JetType fjetType = JetType::kData
							     );


 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  void                                RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont,Bool_t bTruth);

  
  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted. 
  Double_t                            fShapesVar_Det[5];       // jet shapes used for the tagging
  Double_t                            fShapesVar_True[5];      // jet shapes used for the tagging
  JetShapeSub                         fJetShapeSub;            // jet subtraction to be used
  JetType                             fJetType;                // jet type data/embedded
  Float_t                             fPtThreshold;            // jet pt threshold
  Float_t                             fSharedFractionPtMin;    // minimum pt shared fraction to be used to match jets
  Int_t                               fReclusteringAlgo;
  
  Bool_t                              fCentSelectOn;                // switch on/off centrality selection
  Float_t                             fCentMin;                     // min centrality value
  Float_t                             fCentMax;                     // max centrality value
  Double_t                            fJetRadius;                   // radius used in jet finding
  Float_t                             fAddMedScatPtFrac;
  Float_t                             fAddMedScatN;
  Bool_t                              fAddMedScat;
  Bool_t                              fDoSubJetAreaSub;
  
  TH1F                                *fhJetPt;
  TH1F                                *fhJetPhi;
  TH1F                                *fhJetEta;
  TH1F                                *fhDetJetPt_Matched;
 
  TTree                               *fTreeRecursive_Det;  
  TTree                               *fTreeRecursive_True; 

 private:
  AliAnalysisTaskRecursiveSoftDrop(const AliAnalysisTaskRecursiveSoftDrop&);            // not implemented
  AliAnalysisTaskRecursiveSoftDrop &operator=(const AliAnalysisTaskRecursiveSoftDrop&); // not implemented

  ClassDef(AliAnalysisTaskRecursiveSoftDrop, 1)
};
#endif

