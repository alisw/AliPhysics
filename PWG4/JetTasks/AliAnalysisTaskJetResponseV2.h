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

   virtual Int_t      GetNInputTracks();
   virtual THnSparse* NewTHnSparseF(const char* name, UInt_t entries, UInt_t opt);
   virtual void       GetDimParams(Int_t iEntry, Bool_t hr, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
   virtual Int_t      GetPtHardBin(Double_t ptHard);
   virtual Double_t   GetPt(AliAODJet* j, Int_t mode);

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
   virtual void     SetCentMin(Float_t cent) { fCentMin = cent; }
   virtual void     SetCentMax(Float_t cent) { fCentMax = cent; }
   virtual void     SetNInputTracksMin(Int_t nTr) { fNInputTracksMin = nTr; }
   virtual void     SetNInputTracksMax(Int_t nTr) { fNInputTracksMax = nTr; }
   virtual void     SetJetEtaMin(Float_t eta) { fJetEtaMin = eta; }
   virtual void     SetJetEtaMax(Float_t eta) { fJetEtaMax = eta; }
   virtual void     SetJetPtMin(Float_t pt) { fJetPtMin = pt; }
   virtual void     SetJetTriggerExclude(UChar_t i) { fJetTriggerExcludeMask = i; }
   virtual void     SetJetPtFractionMin(Float_t frac) { fJetPtFractionMin = frac; }
   virtual void     SetNMatchJets(Int_t n) { fNMatchJets = n; }
   virtual void     SetFillEvent(Bool_t b) { fbEvent = b; }
   virtual void     SetFillJetsMismatch1(Bool_t b) { fbJetsMismatch1 = b; }
   virtual void     SetFillJetsMismatch2(Bool_t b) { fbJetsMismatch2 = b; }
   virtual void     SetFillJetsRp(Bool_t b) { fbJetsRp = b; }
   virtual void     SetFillJetsDeltaPt(Bool_t b) { fbJetsDeltaPt = b; }
   virtual void     SetFillJetsEta(Bool_t b) { fbJetsEta = b; }
   virtual void     SetFillJetsPhi(Bool_t b) { fbJetsPhi = b; }
   virtual void     SetFillJetsArea(Bool_t b) { fbJetsArea = b; }
   virtual void     SetFillJetsBeforeCut1(Bool_t b) { fbJetsBeforeCut1 = b; }
   virtual void     SetFillJetsBeforeCut2(Bool_t b) { fbJetsBeforeCut2 = b; }
   virtual void     SetKeepJets(Bool_t b = kTRUE) { fKeepJets = b; }

private:
   // ESD/AOD events
   AliESDEvent *fESD;    //! ESD object
   AliAODEvent *fAOD;    //! AOD event

   // jets to compare
   TString fJetBranchName[2]; //  name of jet branches to compare
   TList *fListJets[2];       //! jet lists

   TString fBackgroundBranch;

   // event selection
   Bool_t fIsPbPb;         // is Pb-Pb (fast embedding) or p-p (detector response)
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
   Float_t fJetEtaMin;     // lower bound on eta for found jets
   Float_t fJetEtaMax;     // upper bound on eta for found jets
   Float_t fJetPtMin;      // minimum jet pT
   UChar_t fJetTriggerExcludeMask; // mask for jet triggeres to exclude
   Float_t fJetPtFractionMin; // minimum fraction for positiv match of jets
   Int_t   fNMatchJets;       // maximal nb. of jets taken for matching
   Double_t fMatchMaxDist;     // maximal distance of matching jets
   Bool_t  fKeepJets;          // keep jets with negative pt after background subtraction


   // output objects
   const Int_t fkNbranches;                   //! number of branches to be read
   const Int_t fkEvtClasses;                  //! number of event classes
   TList *fOutputList;                        //! output data container
   Bool_t fbEvent;                            // fill fhnEvent
   Bool_t fbJetsMismatch1;                    // fill fhnJetsMismatch1
   Bool_t fbJetsMismatch2;                    // fill fhnJetsMismatch2
   Bool_t fbJetsRp;                           // fill fhnJetsRp
   Bool_t fbJetsDeltaPt;                      // fill fhnJetsDeltaPt
   Bool_t fbJetsEta;                          // fill fhnJetsEta
   Bool_t fbJetsPhi;                          // fill fhnJetsEta
   Bool_t fbJetsArea;                         // fill fhnJetsArea
   Bool_t fbJetsBeforeCut1;                   // fill fhnJetsBeforeCut1
   Bool_t fbJetsBeforeCut2;                   // fill fhnJetsBeforeCut2
   TH1I  *fHistEvtSelection;                  //! event selection statistic
   TH1I  *fHistJetSelection;                  //! jet selection statistic
   TH2F  *fh2JetSelection;                    //! jet selection statistic, with probe jet pt
   THnSparse *fhnEvent;                       //! variables per event
   THnSparse *fhnJetsMismatch1;               //! variables per jet
   THnSparse *fhnJetsMismatch2;               //! variables per jet
   THnSparse *fhnJetsRp;                      //! variables per jet
   THnSparse *fhnJetsDeltaPt;                 //! variables per jet
   THnSparse *fhnJetsEta;                     //! variables per jet
   THnSparse *fhnJetsPhi;                     //! variables per jet
   THnSparse *fhnJetsArea;                    //! variables per jet
   THnSparse *fhnJetsBeforeCut1;               //! variables per jet before acceptance cut
   THnSparse *fhnJetsBeforeCut2;               //! variables per jet before acceptance cut

   AliAnalysisTaskJetResponseV2(const AliAnalysisTaskJetResponseV2&); // not implemented
   AliAnalysisTaskJetResponseV2& operator=(const AliAnalysisTaskJetResponseV2&); // not implemented

   ClassDef(AliAnalysisTaskJetResponseV2, 4);
};

#endif

