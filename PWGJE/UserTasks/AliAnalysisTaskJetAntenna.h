#ifndef ALIANALYSISTASKJETANTENNA_H
#define ALIANALYSISTASKJETANTENNA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// This task computes several jet observables like 
// the fraction of energy in inner and outer coronnas,
// the distance from track to jet axis and a 
// correlation strength distribution of particles inside jets.    
// Author: lcunquei@cern.ch
// *******************************************

class TH1F;
class TH1I;
class TH2F;
class TH3F;
class THnSparse;
class TRandom3;
class AliESDEvent;
class AliAODExtension;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"

class AliAnalysisTaskJetAntenna : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskJetAntenna();
   AliAnalysisTaskJetAntenna(const char *name);
   virtual ~AliAnalysisTaskJetAntenna();
   virtual void     LocalInit() {Init();}
   virtual void     Init();
   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option);
   virtual void     Terminate(const Option_t*);

   Double_t RelativePhi(Double_t angle1,Double_t angle2);     
   Int_t   GetPhiBin(Double_t phi);
   Int_t      GetPtHardBin(Double_t ptHard);
  
   virtual AliVEvent::EOfflineTriggerTypes GetOfflineTrgMask() const { return fOfflineTrgMask; }
   virtual void     GetBranchNames(TString &branch1, TString &branch2) const { branch1 = fJetBranchName[0]; branch2 = fJetBranchName[1]; }
   virtual Bool_t   GetIsPbPb() const { return fIsPbPb; }
   virtual Int_t    GetMinContribVtx() const { return fMinContribVtx; };
   virtual Float_t  GetVtxZMin() const { return fVtxZMin; }
   virtual Float_t  GetVtxZMax() const { return fVtxZMax; }
   virtual Int_t    GetEvtClassMin() const { return fEvtClassMin; }
   virtual Int_t    GetEvtClassMax() const { return fEvtClassMax; }
   virtual Float_t  GetCentMin() const { return fCentMin; }
   virtual Float_t  GetCentMax() const { return fCentMax; }
   virtual Int_t    GetNInputTracksMin() const { return fNInputTracksMin; }
   virtual Int_t    GetNInputTracksMax() const { return fNInputTracksMax; } 
   virtual Float_t  GetJetEtaMin() const { return fJetEtaMin; }
   virtual Float_t  GetJetEtaMax() const { return fJetEtaMax; }
   virtual Float_t  GetJetPtMin() const { return fJetPtMin; }
   virtual Float_t  GetJetPtFractionMin() const { return fJetPtFractionMin; }
   virtual Int_t    GetNMatchJets() const { return fNMatchJets; }
   virtual void     SetBranchNames(const TString &branch1, const TString &branch2);
   virtual void     SetBackgroundBranch(TString &branch) { fBackgroundBranch = branch;}
   virtual void     SetIsPbPb(Bool_t b=kTRUE) { fIsPbPb = b; }
   virtual void     SetOfflineTrgMask(AliVEvent::EOfflineTriggerTypes mask) { fOfflineTrgMask = mask; }
   virtual void     SetMinContribVtx(Int_t n) { fMinContribVtx = n; }
   virtual void     SetVtxZMin(Float_t z) { fVtxZMin = z; }
   virtual void     SetVtxZMax(Float_t z) { fVtxZMax = z; }
   virtual void     SetEvtClassMin(Int_t evtClass) { fEvtClassMin = evtClass; }
   virtual void     SetEvtClassMax(Int_t evtClass) { fEvtClassMax = evtClass; }
   virtual void     SetFilterMask(UInt_t i){fFilterMask = i;}
   virtual void     SetFilterMaskBestPt(UInt_t i){fFilterMaskBestPt = i;}
   virtual void     SetFilterType(Int_t iType){fFilterType=iType;}
  
   virtual void     SetCentMin(Float_t cent) { fCentMin = cent; }
   virtual void     SetCentMax(Float_t cent) { fCentMax = cent; }
   virtual void     SetNInputTracksMin(Int_t nTr) { fNInputTracksMin = nTr; }
   virtual void     SetNInputTracksMax(Int_t nTr) { fNInputTracksMax = nTr; }
   virtual void     SetRequireITSRefit(Int_t nref) {fRequireITSRefit=nref;}
   virtual void     SetSharedClusterCut(Int_t docut){fApplySharedClusterCut=docut;}
      
   virtual void     SetNRPBins(Int_t bins){fNRPBins=bins;}
   virtual void     SetSemigoodCorrect(Int_t yesno){fSemigoodCorrect=yesno;}
   virtual void     SetDoMatching(Bool_t b=kFALSE){fDoMatching=b;}
   virtual void     SetHolePos(Float_t poshole) { fHolePos = poshole; }
   virtual void     SetHoleWidth(Float_t holewidth) { fHoleWidth = holewidth; }
   virtual void     SetTrackTypeRec(Int_t i){fTrackTypeRec = i;}
   virtual void     SetTMCut(Int_t i){fCutTM=i;}
   virtual void     SetJetEtaMin(Float_t eta) { fJetEtaMin = eta; }
   virtual void     SetJetEtaMax(Float_t eta) { fJetEtaMax = eta; }
   virtual void     SetJetPtMin(Float_t pt) { fJetPtMin = pt; }
   virtual void     SetJetTriggerExclude(UChar_t i) { fJetTriggerExcludeMask = i; }
   virtual void     SetJetPtFractionMin(Float_t frac) { fJetPtFractionMin = frac; }
   virtual void     SetNMatchJets(Int_t n) { fNMatchJets = n; }
   virtual void     SetKeepJets(Bool_t b = kTRUE) { fKeepJets = b; }
   virtual void     SetNonStdFile(char* c){fNonStdFile = c;} 
    enum {kTrackUndef = 0, kTrackAOD, kTrackKineAll,kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};

private:
   // ESD/AOD events
   AliESDEvent *fESD;    //! ESD object
   AliAODEvent *fAODIn;    //! AOD event for AOD input tracks
    AliAODEvent *fAODOut;    //! AOD event 
    AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
   Int_t   GetListOfTracks(TList *list);
   Int_t   GetListOfTracksExtra(TList* list);
   Int_t   SelectTrigger(TList *list,Double_t minT,Double_t maxT,Int_t &number);
   Int_t   GetHardestTrackBackToJet(AliAODJet *jet);
   Int_t   GetListOfTracksCloseToJet(TList *list,AliAODJet *jet);
  
   // jets to compare
   TString fJetBranchName[2]; //  name of jet branches to compare
   TList *fListJets[2];       //! jet lists

   TString fBackgroundBranch;
   TString       fNonStdFile; // name of delta aod file to catch the extension
   // event selection
   Bool_t fIsPbPb;         // is Pb-Pb (fast embedding) or p-p (detector response)
   AliVEvent::EOfflineTriggerTypes fOfflineTrgMask; // mask of offline triggers to accept
   Int_t   fMinContribVtx; // minimum number of track contributors for primary vertex
   Float_t fVtxZMin;	  // lower bound on vertex z
   Float_t fVtxZMax;	  // upper bound on vertex z
   Int_t   fEvtClassMin;	  // lower bound on event class
   Int_t   fEvtClassMax;	  // upper bound on event class
   UInt_t  fFilterMask;            // filter bit for slecected tracks
   UInt_t  fFilterMaskBestPt;      // filter bit for selected hig pt tracks (best quality)
   UInt_t  fFilterType;           // type of slected tracks parrallel to filtermask
   Float_t fCentMin;	  // lower bound on centrality
   Float_t fCentMax;	  // upper bound on centrality
   Int_t   fNInputTracksMin;  // lower bound of nb. of input tracks
   Int_t   fNInputTracksMax;  // upper bound of nb. of input tracks
   Int_t   fRequireITSRefit;
   Int_t   fApplySharedClusterCut; // flag to apply shared cluster cut (needed for some AODs where this cut was not applied in the filtering)
   Int_t   fTrackTypeRec;
   Float_t   fRPAngle;
   Int_t   fNRPBins;
   Int_t   fSemigoodCorrect;
   Bool_t fDoMatching;
   Float_t fHolePos;
   Float_t fHoleWidth; 
   Float_t fCutTM;             //lower pt cut for particles entering the TM axis calculation 
   Float_t fJetEtaMin;        // lower bound on eta for found jets
   Float_t fJetEtaMax;        // upper bound on eta for found jets
   Int_t   fNevents;          // number of events
   Int_t   fTindex;           // index reference
   Float_t fJetPtMin;         // minimum jet pT
   UChar_t fJetTriggerExcludeMask; // mask for jet triggeres to exclude
   Float_t fJetPtFractionMin; // minimum fraction for positiv match of jets
   Int_t   fNMatchJets;       // maximal nb. of jets taken for matching
   Double_t fMatchMaxDist;     // maximal distance of matching jets
   Bool_t  fKeepJets;          // keep jets with negative pt after background subtraction


   // output objects
   const Int_t fkNbranches;                   //! number of branches to be read
   const Int_t fkEvtClasses;                  //! number of event classes
   
   TList *fOutputList;                        //! output data container
   
    TH1I  *fHistEvtSelection;                  //! event selection statistic
     TH2F*      fh2JetEntries;             //centrality bias of triggers 
     TH2F*      fh2Circularity;             //jet density
     TH2F*      fh2JetAxisPhi;              //smearing jet axis
     THnSparse   *fhnJetTM;               //Recoil jet spectrum


   AliAnalysisTaskJetAntenna(const AliAnalysisTaskJetAntenna&); // not implemented
   AliAnalysisTaskJetAntenna& operator=(const AliAnalysisTaskJetAntenna&); // not implemented

   ClassDef(AliAnalysisTaskJetAntenna, 6);
};

#endif

