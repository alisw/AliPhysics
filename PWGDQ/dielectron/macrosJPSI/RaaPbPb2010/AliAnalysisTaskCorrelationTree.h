#ifndef ALIANALYSISTASKCORRELATIONTREE_H
#define ALIANALYSISTASKCORRELATIONTREE_H

// Analysis task for creating a reduced tree containing event, track and resonance candidate information
// Author: Ionut-Cristian Arsene (i.c.arsene@gsi.de)

#include "TList.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisCuts;
class TTree;
class TFile;
class AliESDv0Cuts;
class AliKFVertex;
class AliCorrelationReducedEvent;
class AliCorrelationReducedEventFriend;
class AliCorrelationReducedPair;
class AliDielectron;
class AliFlowTrackCuts;

//_________________________________________________________________________
class AliAnalysisTaskCorrelationTree : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskCorrelationTree();
  AliAnalysisTaskCorrelationTree(const char *name);
  virtual ~AliAnalysisTaskCorrelationTree(){  }

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(UInt_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }
  void SetRejectPileup(Bool_t pileup=kTRUE)     { fRejectPileup=pileup;     }
  
  // Cuts for selection of event to be written to tree
  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}
  // Cuts for selecting tracks included in the tree
  void SetTrackFilter(AliAnalysisCuts * const filter) {fTrackFilter=filter;}
  // Cuts for selecting tracks to be used for Q vector calculation
  void SetFlowTrackFilter(AliAnalysisCuts * const filter) {fFlowTrackFilter = filter;}
  
  // Cuts for selecting V0s
  void SetK0sPionCuts(AliAnalysisCuts * const filter) {fK0sPionCuts=filter;}
  void SetLambdaProtonCuts(AliAnalysisCuts * const filter) {fLambdaProtonCuts=filter;}
  void SetLambdaPionCuts(AliAnalysisCuts * const filter) {fLambdaPionCuts=filter;}
  void SetK0sCuts(AliESDv0Cuts* const cuts) {fK0sCuts = cuts;}
  void SetLambdaCuts(AliESDv0Cuts* const cuts) {fLambdaCuts = cuts;}
  void SetK0sMassRange(Double_t min=0.4, Double_t max=0.6) {fK0sMassRange[0]=min; fK0sMassRange[1]=max;}
  void SetLambdaMassRange(Double_t min=1.08, Double_t max=1.15) {fLambdaMassRange[0]=min; fLambdaMassRange[1]=max;}
  void SetV0Histograms(AliDielectronHistos * const histos) {fV0Histos=histos;}
  
  // Add dielectron objects to the list. These contain cuts and histogram definitions
  void AddDielectron(AliDielectron * const die) { fListDielectron.Add(die); }
 
 private:

  TList fListDielectron;             // List of dielectron framework instances
  TList fListHistos;                 //! List of histogram managers in the dielectron framework classes
  
  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  Bool_t fRejectPileup;              // pileup rejection wanted

  AliAnalysisCuts *fEventFilter;     // event filter
  AliAnalysisCuts *fTrackFilter;     // filter for the hadrons to be correlated with the dielectrons
  AliAnalysisCuts *fFlowTrackFilter; // filter for the barrel tracks to be used for the Q-vector
  
  AliESDv0Cuts *fK0sCuts;            // v0 standard filter for K0s->pi+pi-
  AliESDv0Cuts *fLambdaCuts;         // v0 standard filter for Lambda0->p + pi
  AliAnalysisCuts *fK0sPionCuts;     // filter for pions from K0s
  AliAnalysisCuts *fLambdaProtonCuts;// filter for protons from Lambda
  AliAnalysisCuts *fLambdaPionCuts;  // filter for pions from Lambda
  Double_t fK0sMassRange[2];           // mass range for allowed K0s pairs
  Double_t fLambdaMassRange[2];        // mass range for allowed Lambda pairs
  AliDielectronHistos* fV0Histos;    // histogram manager for V0s

  TFile *fTreeFile;                  // output file containing the tree
  TTree *fTree;                      //! Reduced event tree
  TTree *fFriendTreeFile;            // output file containing the friend tree
  TTree *fFriendTree;                //! Reduced event tree with friend info (event plane, etc.)
  AliCorrelationReducedEvent *fReducedEvent;    // reduced event wise information
  AliCorrelationReducedEventFriend *fReducedEventFriend;    // friend reduced event wise information
  
  AliFlowTrackCuts* fFlowTrackCuts;   // flow track cuts object
  
  void FillEventInfo();                     // fill reduced event information
  void FillFriendEventInfo();               // fill reduced event friend information
  void FillTrackInfo();                     // fill reduced track information
  void FillDielectronPairInfo(AliDielectron* die, Short_t iDie);  // fill dielectron reduced pair information
  void FillV0PairInfo();                    // fill V0 reduced pair information
  AliCorrelationReducedPair* FillV0PairInfo(AliESDv0* v0, Int_t id, AliESDtrack* legPos, AliESDtrack* legNeg, AliKFVertex* vtxKF, Bool_t chargesAreCorrect);
  UChar_t EncodeTPCClusterMap(AliESDtrack* track);
  void FillCaloClusters();
  
  AliAnalysisTaskCorrelationTree(const AliAnalysisTaskCorrelationTree &c);
  AliAnalysisTaskCorrelationTree& operator= (const AliAnalysisTaskCorrelationTree &c);

  ClassDef(AliAnalysisTaskCorrelationTree, 1); //Analysis Task for creating a reduced event information tree 
};
#endif
