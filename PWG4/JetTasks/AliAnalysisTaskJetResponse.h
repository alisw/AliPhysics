#ifndef ALIANALYSISTASKJETRESPONSE_H
#define ALIANALYSISTASKJETRESPONSE_H

class TH1F;
class TH2F;
class AliESDEvent;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"

class AliAnalysisTaskJetResponse : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskJetResponse();
  AliAnalysisTaskJetResponse(const char *name);
  virtual ~AliAnalysisTaskJetResponse();

  virtual void     LocalInit() {Init();}
  virtual void     Init();
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(const Option_t*);

  virtual AliVEvent::EOfflineTriggerTypes GetOfflineTrgMask() const { return fOfflineTrgMask; }
  virtual void     GetBranchNames(TString &branch1, TString &branch2) const { branch1 = fJetBranchName[0]; branch2 = fJetBranchName[1]; }
  virtual Float_t  GetVtxZMin() const { return fVtxZMin; }
  virtual Float_t  GetVtxZMax() const { return fVtxZMax; }
  virtual Int_t    GetEvtClassMin() const { return fEvtClassMin; }
  virtual Int_t    GetEvtClassMax() const { return fEvtClassMax; }
  virtual Float_t  GetCentMin() const { return fCentMin; }
  virtual Float_t  GetCentMax() const { return fCentMax; }
  virtual Float_t  GetJetEtaMin() const { return fJetEtaMin; }
  virtual Float_t  GetJetEtaMax() const { return fJetEtaMax; }
  virtual Float_t  GetJetDeltaEta() const { return fJetDeltaEta; }
  virtual Float_t  GetJetDeltaPhi() const { return fJetDeltaPhi; }

  virtual void     SetBranchNames(const TString &branch1, const TString &branch2);
  virtual void     SetOfflineTrgMask(AliVEvent::EOfflineTriggerTypes mask) { fOfflineTrgMask = mask; }
  virtual void     SetVtxZMin(Float_t z) { fVtxZMin = z; }
  virtual void     SetVtxZMax(Float_t z) { fVtxZMax = z; }
  virtual void     SetEvtClassMin(Int_t evtClass) { fEvtClassMin = evtClass; }
  virtual void     SetEvtClassMax(Int_t evtClass) { fEvtClassMax = evtClass; }
  virtual void     SetCentMin(Float_t cent) { fCentMin = cent; }
  virtual void     SetCentMax(Float_t cent) { fCentMax = cent; }
  virtual void     SetJetEtaMin(Float_t eta) { fJetEtaMin = eta; }
  virtual void     SetJetEtaMax(Float_t eta) { fJetEtaMax = eta; }
  virtual void     SetJetDeltaEta(Float_t eta) { fJetDeltaEta = eta; }
  virtual void     SetJetDeltaPhi(Float_t phi) { fJetDeltaPhi = phi; }

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
  Float_t fJetEtaMin;     // lower bound on eta for found jets
  Float_t fJetEtaMax;     // upper bound on eta for found jets
  Float_t fJetDeltaEta;   // max difference in eta to match leading jets
  Float_t fJetDeltaPhi;   // max difference in phi to match leading jets

  // output objects
  const Int_t fkNbranches;                   //! number of branches to be read
  const Int_t fkEvtClasses;                  //! number of event classes
  TList *fOutputList;                        //! output data container
  TH1F *fHistEvtSelection;                   //! event selection statistic
  TH1F **fHistPtLeadingJet;                  //! pt of leading jet
  TH2F **fHistEtaPhiLeadingJet;              //! eta-phi of leading jet
  TH2F **fHistEtaPhiLeadingJetCut;           //! eta-phi of leading jet in eta acceptace per event class
  TH2F  *fHistDeltaEtaDeltaPhiLeadingJet;    //! delta eta vs. delta phi of leading jets
  TH2F **fHistDeltaEtaEtaLeadingJet;         //! delta eta vs. eta of leading jet per event class
  TH2F **fHistDeltaPtEtaLeadingJet;          //! delta eta vs. eta of leading jet per event class
  TH2F  *fHistPtPtExtra;                     //! jet pt response
  TH2F **fHistPtResponse;                    //! histograms per event class
  TH2F **fHistPtSmearing;                    //! emb-jet pt vs (emb+UE - emb) pt
  TH2F **fHistdR;                            //! shift dR of jets vs (emb+UE - emb) pt
  TH2F **fHistArea;                          //! area of jets vs (emb+UE - emb) pt

  AliAnalysisTaskJetResponse(const AliAnalysisTaskJetResponse&); // not implemented
  AliAnalysisTaskJetResponse& operator=(const AliAnalysisTaskJetResponse&); // not implemented

  ClassDef(AliAnalysisTaskJetResponse, 1);
};

#endif
