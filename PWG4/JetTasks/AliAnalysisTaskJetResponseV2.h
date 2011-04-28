#ifndef ALIANALYSISTASKJETRESPONSEV2_H
#define ALIANALYSISTASKJETRESPONSEV2_H

class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class AliESDEvent;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"

class AliAnalysisTaskJetResponseV2 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskJetResponseV2();
  AliAnalysisTaskJetResponseV2(const char *name);
  virtual ~AliAnalysisTaskJetResponseV2();

  virtual void     LocalInit() {Init();}
  virtual void     Init();
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(const Option_t*);
  virtual Int_t    GetNInputTracks();

  virtual AliVEvent::EOfflineTriggerTypes GetOfflineTrgMask() const { return fOfflineTrgMask; }
  virtual void     GetBranchNames(TString &branch1, TString &branch2) const { branch1 = fJetBranchName[0]; branch2 = fJetBranchName[1]; }
  virtual Int_t    GetMinContribVtx() const { return fMinContribVtx; };
  virtual Float_t  GetVtxZMin() const { return fVtxZMin; }
  virtual Float_t  GetVtxZMax() const { return fVtxZMax; }
  virtual Int_t    GetEvtClassMin() const { return fEvtClassMin; }
  virtual Int_t    GetEvtClassMax() const { return fEvtClassMax; }
  virtual Float_t  GetCentMin() const { return fCentMin; }
  virtual Float_t  GetCentMax() const { return fCentMax; }
  virtual Int_t    GetNInputTracksMin() const { return fNInputTracksMin; }
  virtual Int_t    GetNInputTracksMax() const { return fNInputTracksMax; } 
  virtual Float_t  GetReactionPlaneBin() const { return fReactionPlaneBin; }
  virtual Float_t  GetJetEtaMin() const { return fJetEtaMin; }
  virtual Float_t  GetJetEtaMax() const { return fJetEtaMax; }
  virtual Float_t  GetJetPtMin() const { return fJetPtMin; }
  virtual Float_t  GetJetPtFractionMin() const { return fJetPtFractionMin; }
  virtual Int_t    GetNMatchJets() const { return fNMatchJets; }

  virtual void     SetBranchNames(const TString &branch1, const TString &branch2);
  virtual void     SetOfflineTrgMask(AliVEvent::EOfflineTriggerTypes mask) { fOfflineTrgMask = mask; }
  virtual void     SetMinContribVtx(Int_t n) { fMinContribVtx = n; }
  virtual void     SetVtxZMin(Float_t z) { fVtxZMin = z; }
  virtual void     SetVtxZMax(Float_t z) { fVtxZMax = z; }
  virtual void     SetEvtClassMin(Int_t evtClass) { fEvtClassMin = evtClass; }
  virtual void     SetEvtClassMax(Int_t evtClass) { fEvtClassMax = evtClass; }
  virtual void     SetCentMin(Float_t cent) { fCentMin = cent; }
  virtual void     SetCentMax(Float_t cent) { fCentMax = cent; }
  virtual void     SetNInputTracksMin(Int_t nTr) { fNInputTracksMin = nTr; }
  virtual void     SetNInputTracksMax(Int_t nTr) { fNInputTracksMax = nTr; }
  virtual void     SetReactionPlaneBin(Int_t rpBin) { fReactionPlaneBin = rpBin; }
  virtual void     SetJetEtaMin(Float_t eta) { fJetEtaMin = eta; }
  virtual void     SetJetEtaMax(Float_t eta) { fJetEtaMax = eta; }
  virtual void     SetJetPtMin(Float_t pt) { fJetPtMin = pt; }
  virtual void     SetJetPtFractionMin(Float_t pt) { fJetPtFractionMin = pt; }
  virtual void     SetNMatchJets(Int_t n) { fNMatchJets = n; }

 private:
  // ESD/AOD events
  AliESDEvent *fESD;    //! ESD object
  AliAODEvent *fAOD;    //! AOD event

  // jets to compare
  TString fJetBranchName[2]; //  name of jet branches to compare
  TList *fListJets[2];       //! jet lists

  // event selection
  AliVEvent::EOfflineTriggerTypes fOfflineTrgMask; // mask of offline triggers to accept
  Int_t   fMinContribVtx; // minimum number of track contributors for primary vertex
  Float_t fVtxZMin;	  // lower bound on vertex z
  Float_t fVtxZMax;	  // upper bound on vertex z
  Int_t   fEvtClassMin;	  // lower bound on event class
  Int_t   fEvtClassMax;	  // upper bound on event class
  Float_t fCentMin;	  // lower bound on centrality
  Float_t fCentMax;	  // upper bound on centrality
  Int_t   fNInputTracksMin;  // lower bound of nb. of input tracks
  Int_t   fNInputTracksMax;  // upper bound of nb. of input tracks
  Float_t fReactionPlaneBin; // reaction plane bin
  Float_t fJetEtaMin;     // lower bound on eta for found jets
  Float_t fJetEtaMax;     // upper bound on eta for found jets
  Float_t fJetPtMin;      // minimum jet pT
  Float_t fJetPtFractionMin; // minimum fraction for positiv match of jets
  Int_t   fNMatchJets;       // maximal nb. of jets taken for matching
  Double_t fMatchMaxDist;     // maximal distance of matching jets
  

  // output objects
  const Int_t fkNbranches;                   //! number of branches to be read
  const Int_t fkEvtClasses;                  //! number of event classes
  TList *fOutputList;                        //! output data container
  TH1I  *fHistEvtSelection;                  //! event selection statistic
  THnSparse *fhnEvent;                       //! variables per event
  THnSparse *fhnJetsRp;                      //! variables per jet
  THnSparse *fhnJetsDeltaPt;                 //! variables per jet
  THnSparse *fhnJetsEta;                     //! variables per jet
  THnSparse *fhnJetsArea;                    //! variables per jet
  THnSparse *fhnJetsBeforeCut1;               //! variables per jet before acceptance cut
  THnSparse *fhnJetsBeforeCut2;               //! variables per jet before acceptance cut
  
  AliAnalysisTaskJetResponseV2(const AliAnalysisTaskJetResponseV2&); // not implemented
  AliAnalysisTaskJetResponseV2& operator=(const AliAnalysisTaskJetResponseV2&); // not implemented

  ClassDef(AliAnalysisTaskJetResponseV2, 1);
};

#endif
