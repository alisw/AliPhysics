// Analysis task for creating a reduced tree containing event, track and resonance candidate information 
// Author: Ionut-Cristian Arsene (i.c.arsene@gsi.de,i.c.arsene@cern.ch)                                 
// 2012/06/21

#ifndef ALIANALYSISTASKREDUCEDTREE_H
#define ALIANALYSISTASKREDUCEDTREE_H 

#include "TList.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisCuts;
class TTree;
class TFile;
class AliESDv0Cuts;
class AliESDv0KineCuts;
class AliKFVertex;
class AliReducedEvent;
class AliReducedEventFriend;
class AliReducedPair;
class AliDielectron;

//_________________________________________________________________________
class AliAnalysisTaskReducedTree : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskReducedTree();
  AliAnalysisTaskReducedTree(const char *name);
  virtual ~AliAnalysisTaskReducedTree(){  }

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
  void SetGammaElectronCuts(AliAnalysisCuts* const filter) {fGammaElectronCuts=filter;}
  void SetK0sCuts(AliESDv0Cuts* const cuts) {fK0sCuts = cuts;}
  void SetLambdaCuts(AliESDv0Cuts* const cuts) {fLambdaCuts = cuts;}
  void SetGammaConvCuts(AliESDv0KineCuts* const cuts) {fGammaConvCuts = cuts;}
  void SetV0OpenCuts(AliESDv0KineCuts* const cuts) {fV0OpenCuts = cuts;}
  void SetV0StrongCuts(AliESDv0KineCuts* const cuts) {fV0StrongCuts = cuts;}
  void SetK0sMassRange(Double_t min=0.4, Double_t max=0.6) {fK0sMassRange[0]=min; fK0sMassRange[1]=max;}
  void SetLambdaMassRange(Double_t min=1.08, Double_t max=1.15) {fLambdaMassRange[0]=min; fLambdaMassRange[1]=max;}
  void SetGammaConvMassRange(Double_t min=0.0, Double_t max=0.1) {fGammaMassRange[0]=min; fGammaMassRange[1]=max;}
  void SetV0Histograms(AliDielectronHistos * const histos) {fV0Histos=histos;}
  
  // Toggle on/off information branches
  void SetFillTrackInfo(Bool_t flag=kTRUE)        {fFillTrackInfo = flag;}
  void SetFillDielectronInfo(Bool_t flag=kTRUE)   {fFillDielectronInfo = flag;}
  void SetFillV0Info(Bool_t flag=kTRUE)           {fFillV0Info = flag;}
  void SetFillGammaConversions(Bool_t flag=kTRUE) {fFillGammaConversions = flag;}
  void SetFillK0s(Bool_t flag=kTRUE)              {fFillK0s = flag;}
  void SetFillLambda(Bool_t flag=kTRUE)           {fFillLambda = flag;}
  void SetFillALambda(Bool_t flag=kTRUE)          {fFillALambda = flag;}
  void SetFillCaloClusterInfo(Bool_t flag=kTRUE)  {fFillCaloClusterInfo = flag;}
  void SetFillFMDSectorInfo(Bool_t flag=kFALSE)   {fFillFMDSectorInfo = flag;}
  void SetFillFMDChannelInfo(Bool_t flag=kFALSE)  {fFillFMDChannelInfo = flag;}
  void SetFillFriendInfo(Bool_t flag=kTRUE)       {fFillFriendInfo = flag;}
  
  // Add dielectron objects to the list. These contain cuts and histogram definitions
  void AddDielectron(AliDielectron * const die) { fListDielectron.Add(die); }
 
 private:

  TList fListDielectron;             // List of dielectron framework instances
  TList fListHistos;                 //! List of histogram managers in the dielectron framework classes
  
  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  Bool_t fRejectPileup;              // pileup rejection wanted

  Bool_t fFillTrackInfo;             // fill track information
  Bool_t fFillDielectronInfo;        // fill dielectrons
  Bool_t fFillV0Info;                // fill the V0 information
  Bool_t fFillGammaConversions;      // fill gamma conversions
  Bool_t fFillK0s;                   // fill the K0s V0s
  Bool_t fFillLambda;                // fill the lambda V0s
  Bool_t fFillALambda;               // fill the anti-lambda V0s
  Bool_t fFillCaloClusterInfo;       // fill the calorimeter clusters
  Bool_t fFillFMDSectorInfo;         // fill the FMD info for every sector
  Bool_t fFillFMDChannelInfo;        // fill the FMD info for every channel
  Bool_t fFillFriendInfo;            // fill friend tree information

  AliAnalysisCuts *fEventFilter;     // event filter
  AliAnalysisCuts *fTrackFilter;     // filter for the hadrons to be correlated with the dielectrons
  AliAnalysisCuts *fFlowTrackFilter; // filter for the barrel tracks to be used for the Q-vector
  
  AliESDv0Cuts *fK0sCuts;            // v0 standard filter for K0s->pi+pi-
  AliESDv0Cuts *fLambdaCuts;         // v0 standard filter for Lambda0->p + pi
  AliESDv0KineCuts *fGammaConvCuts;  // v0 standard filter for gamma conversions
  AliAnalysisCuts *fK0sPionCuts;     // filter for pions from K0s
  AliAnalysisCuts *fLambdaProtonCuts;   // filter for protons from Lambda
  AliAnalysisCuts *fLambdaPionCuts;     // filter for pions from Lambda
  AliAnalysisCuts *fGammaElectronCuts;  // filter for electrons from gamma conversions
  AliESDv0KineCuts *fV0OpenCuts;       // v0 strong filter for tagged V0s
  AliESDv0KineCuts *fV0StrongCuts;     // v0 strong filter for tagged V0s
  
    
  Double_t fK0sMassRange[2];         // mass range for allowed K0s pairs
  Double_t fLambdaMassRange[2];      // mass range for allowed Lambda pairs
  Double_t fGammaMassRange[2];       // mass range for allowed Gamma conversion pairs
  AliDielectronHistos* fV0Histos;    // histogram manager for V0s

  TFile *fTreeFile;                  //! output file containing the tree
  TTree *fTree;                      //! Reduced event tree
  TFile *fFriendTreeFile;            //! output file containing the friend tree
  TTree *fFriendTree;                //! Reduced event tree with friend info (event plane, etc.)
  AliReducedEvent *fReducedEvent;    //! reduced event wise information
  AliReducedEventFriend *fReducedEventFriend;    //! friend reduced event wise information
  
  void FillEventInfo();                     // fill reduced event information
  void FillFriendEventInfo();               // fill reduced event friend information
  void FillTrackInfo();                     // fill reduced track information
  void FillDielectronPairInfo(AliDielectron* die, Short_t iDie);  // fill dielectron reduced pair information
  void FillV0PairInfo();                    // fill V0 reduced pair information
  AliReducedPair* FillV0PairInfo(AliESDv0* v0, Int_t id, AliESDtrack* legPos, AliESDtrack* legNeg, AliKFVertex* vtxKF, Bool_t chargesAreCorrect);
  UChar_t EncodeTPCClusterMap(AliVParticle* track, Bool_t isAOD);
  void FillCaloClusters();
  void FillFMDInfo();
  Int_t GetSPDTrackletMultiplicity(AliVEvent* event, Float_t lowEta, Float_t highEta);
  
  AliAnalysisTaskReducedTree(const AliAnalysisTaskReducedTree &c);
  AliAnalysisTaskReducedTree& operator= (const AliAnalysisTaskReducedTree &c);

  ClassDef(AliAnalysisTaskReducedTree, 3); //Analysis Task for creating a reduced event information tree 
};
#endif
