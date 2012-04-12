#ifndef ALIANALYSISTASKJETCORE_H
#define ALIANALYSISTASKJETCORE_H

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
class AliESDEvent;
class AliAODExtension;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"

class AliAnalysisTaskJetCore : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskJetCore();
   AliAnalysisTaskJetCore(const char *name);
   virtual ~AliAnalysisTaskJetCore();
   virtual void     LocalInit() {Init();}
   virtual void     Init();
   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option);
   virtual void     Terminate(const Option_t*);

   virtual Int_t      GetNInputTracks();
     
   Double_t RelativePhi(Double_t angle1,Double_t angle2);     
    virtual THnSparse* NewTHnSparseF(const char* name, UInt_t entries);
   virtual void       GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
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
   virtual void     SetRadioFrac(Float_t radiofrac) { fRadioFrac = radiofrac; }
   virtual void     SetMinDist(Float_t minDist) { fMinDist = minDist; }
   virtual void     SetCentMin(Float_t cent) { fCentMin = cent; }
   virtual void     SetCentMax(Float_t cent) { fCentMax = cent; }
   virtual void     SetNInputTracksMin(Int_t nTr) { fNInputTracksMin = nTr; }
   virtual void     SetNInputTracksMax(Int_t nTr) { fNInputTracksMax = nTr; }
   virtual void     SetAngStructCloseTracks(Int_t yesno){fAngStructCloseTracks=yesno;}
   virtual void     SetCheckMethods(Int_t yesno){fCheckMethods=yesno;}
   virtual void     SetEventMixing(Int_t yesno){fDoEventMixing=yesno;}
   virtual void     SetJetEtaMin(Float_t eta) { fJetEtaMin = eta; }
   virtual void     SetJetEtaMax(Float_t eta) { fJetEtaMax = eta; }
   virtual void     SetJetPtMin(Float_t pt) { fJetPtMin = pt; }
   virtual void     SetJetTriggerExclude(UChar_t i) { fJetTriggerExcludeMask = i; }
   virtual void     SetJetPtFractionMin(Float_t frac) { fJetPtFractionMin = frac; }
   virtual void     SetNMatchJets(Int_t n) { fNMatchJets = n; }
   virtual void     SetFillEvent(Bool_t b) { fbEvent = b; }
   virtual void     SetKeepJets(Bool_t b = kTRUE) { fKeepJets = b; }
   virtual void     SetNonStdFile(char* c){fNonStdFile = c;} 


private:
   // ESD/AOD events
   AliESDEvent *fESD;    //! ESD object
   AliAODEvent *fAODIn;    //! AOD event for AOD input tracks
    AliAODEvent *fAODOut;    //! AOD event 
    AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
   Int_t   GetListOfTracks(TList *list);
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
   Float_t fRadioFrac;                          //!size of the concentric cone
   Float_t fMinDist;   
   Float_t fCentMin;	  // lower bound on centrality
   Float_t fCentMax;	  // upper bound on centrality
   Int_t   fNInputTracksMin;  // lower bound of nb. of input tracks
   Int_t   fNInputTracksMax;  // upper bound of nb. of input tracks
   Int_t   fAngStructCloseTracks;//only constituents or all tracks with R<0.8 for the angular structure
   Int_t   fCheckMethods;     //to look into more detail into the core
   Int_t   fDoEventMixing;
   Float_t fJetEtaMin;        // lower bound on eta for found jets
   Float_t fJetEtaMax;        // upper bound on eta for found jets
   Int_t   fNevents;          // number of events
   Int_t   fTindex;           // index reference
   Int_t   fTrigBufferIndex;  //index for the buffering
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
   Bool_t fbEvent;                            // fill fhnEvent
   TH1I  *fHistEvtSelection;                  //! event selection statistic
   THnSparse *fhnDeltaR;                     //! variables per jet
   THnSparse *fhnMixedEvents;                //!mixed events matrix
 

   TH2F      *fh2JetCoreMethod1C10;          //Energy fraction in the core C10 method 1
   TH2F      *fh2JetCoreMethod2C10;          //Energy fraction in the core C10 method 2
   TH2F      *fh2JetCoreMethod1C20;          //Energy fraction in the core C20 method 1 
   TH2F      *fh2JetCoreMethod2C20;          //Energy fraction in the core C20 method 2
   TH2F      *fh2JetCoreMethod1C30;          //Energy fraction in the core C30 method 1
   TH2F      *fh2JetCoreMethod2C30;          //Energy fraction in the core C30 method 2
   TH2F      *fh2JetCoreMethod1C60;          //Energy fraction in the core C60 method 1
   TH2F      *fh2JetCoreMethod2C60;          //Energy fraction in the core C60 method 2

     TH2F*      fh2AngStructpt1C10;         //Average two particle correlation function:C10 pt1
     TH2F*      fh2AngStructpt2C10;         //C10 pt2
     TH2F*      fh2AngStructpt3C10;         //C10 pt3
     TH2F*      fh2AngStructpt4C10;         //C10 pt4
     TH2F*      fh2AngStructpt1C20;         //C20 pt1
     TH2F*      fh2AngStructpt2C20;         //C20 pt2
     TH2F*      fh2AngStructpt3C20;         //C20 pt3 
     TH2F*      fh2AngStructpt4C20;         //C20 pt4
     TH2F*      fh2AngStructpt1C30;         //C30 pt1
     TH2F*      fh2AngStructpt2C30;         //C30 pt2
     TH2F*      fh2AngStructpt3C30;         //C30 pt3
     TH2F*      fh2AngStructpt4C30;         //C30 pt4
     TH2F*      fh2AngStructpt1C60;         //C60 pt1
     TH2F*      fh2AngStructpt2C60;         //C60 pt2
     TH2F*      fh2AngStructpt3C60;         //C60 pt3
     TH2F*      fh2AngStructpt4C60;         //C60 pt4
     TH2F*      fh2JetsumHT3R2;             //jet shape 02
     TH2F*      fh2JetsumHT3R2gap;             //jet shape 02

     TH2F*      fh2JetsumHT3R4;             //jet shape 02
     TH2F*      fh2JetsumHT3R4gap;             //jet shape 02

     TH2F*      fh2JetsumHT3R6;             //jet shape 02
     TH2F*      fh2JetsumHT3R6gap;             //jet shape 02

     TH2F*      fh2JetsumHT3R8;             //jet shape 02
     TH2F*      fh2JetsumHT3R8gap;             //jet shape 02

     TH2F*      fh2JetsumHT3R10;             //jet shape 12
     TH2F*      fh2JetsumHT3R10gap;             //jet shape 12


     TH3F*      fh3spectriggered;           //triggered spectra
     TH3F*      fh3specbiased;              //biased spectra
     TH3F*      fh3spectot;                 //the two combined  
     TH3F*      fh3spectotb;                 //the two combined    
     
     Double_t            fTrigBuffer[10][7];      //!buffer for triggers   

   AliAnalysisTaskJetCore(const AliAnalysisTaskJetCore&); // not implemented
   AliAnalysisTaskJetCore& operator=(const AliAnalysisTaskJetCore&); // not implemented

   ClassDef(AliAnalysisTaskJetCore, 4);
};

#endif

