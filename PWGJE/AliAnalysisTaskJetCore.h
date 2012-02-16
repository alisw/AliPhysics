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
   virtual void     SetRadioFrac(Float_t radiofrac) { fRadioFrac = radiofrac; }
   virtual void     SetMinDist(Float_t minDist) { fMinDist = minDist; }
   virtual void     SetCentMin(Float_t cent) { fCentMin = cent; }
   virtual void     SetCentMax(Float_t cent) { fCentMax = cent; }
   virtual void     SetNInputTracksMin(Int_t nTr) { fNInputTracksMin = nTr; }
   virtual void     SetNInputTracksMax(Int_t nTr) { fNInputTracksMax = nTr; }
   virtual void     SetAngStructCloseTracks(Int_t yesno){fAngStructCloseTracks=yesno;}
   virtual void     SetJetEtaMin(Float_t eta) { fJetEtaMin = eta; }
   virtual void     SetJetEtaMax(Float_t eta) { fJetEtaMax = eta; }
   virtual void     SetJetPtMin(Float_t pt) { fJetPtMin = pt; }
   virtual void     SetJetTriggerExclude(UChar_t i) { fJetTriggerExcludeMask = i; }
   virtual void     SetJetPtFractionMin(Float_t frac) { fJetPtFractionMin = frac; }
   virtual void     SetNMatchJets(Int_t n) { fNMatchJets = n; }
   virtual void     SetFillEvent(Bool_t b) { fbEvent = b; }
   virtual void     SetKeepJets(Bool_t b = kTRUE) { fKeepJets = b; }
   virtual void SetNonStdFile(char* c){fNonStdFile = c;} 


private:
   // ESD/AOD events
   AliESDEvent *fESD;    //! ESD object
   AliAODEvent *fAOD;    //! AOD event
    AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
   Int_t   GetListOfTracks(TList *list);
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
   Float_t fRadioFrac;                          //!size of the concentric cone
   Float_t fMinDist;   
   Float_t fCentMin;	  // lower bound on centrality
   Float_t fCentMax;	  // upper bound on centrality
   Int_t   fNInputTracksMin;  // lower bound of nb. of input tracks
   Int_t   fNInputTracksMax;  // upper bound of nb. of input tracks
   Int_t   fAngStructCloseTracks;//only constituents or all tracks with R<0.8 for the angular structure
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
   TH1I  *fHistEvtSelection;                  //! event selection statistic
   TH1I  *fHistJetSelection;                  //! jet selection statistic
   TH2F  *fh2JetSelection;                    //! jet selection statistic, with 
     

   TH2F      *fh2JetCoreMethod1C10;          //Energy fraction in the core C10 method 1
   TH2F      *fh2JetCoreMethod2C10;          //Energy fraction in the core C10 method 2
   TH2F      *fh2JetCoreMethod3C10;          //Energy fraction in the core C10 method 3
   TH2F      *fh2JetCoreMethod1C20;          //Energy fraction in the core C20 method 1 
   TH2F      *fh2JetCoreMethod2C20;          //Energy fraction in the core C20 method 2
   TH2F      *fh2JetCoreMethod3C20;          //Energy fraction in the core C20 method 3
   TH2F      *fh2JetCoreMethod1C30;          //Energy fraction in the core C30 method 1
   TH2F      *fh2JetCoreMethod2C30;          //Energy fraction in the core C30 method 2
   TH2F      *fh2JetCoreMethod3C30;          //Energy fraction in the core C30 method 3
   TH2F      *fh2JetCoreMethod1C60;          //Energy fraction in the core C60 method 1
   TH2F      *fh2JetCoreMethod2C60;          //Energy fraction in the core C60 method 2
   TH2F      *fh2JetCoreMethod3C60;          //Energy fraction in the core C60 method 3
   TH2F      *fh2JetCoreMethod3C10lead;          //Energy fraction in the core C30 method 3
   TH2F      *fh2JetCoreMethod3C20lead;          //Energy fraction in the core C60 method 1
   TH2F      *fh2JetCoreMethod3C30lead;          //Energy fraction in the core C60 method 2
   TH2F      *fh2JetCoreMethod3C60lead;          //Energy fraction in the core C60 method 3
   TH2F      *fh2JetCoreMethod3C10sublead;          //Energy fraction in the core C30 method 3
   TH2F      *fh2JetCoreMethod3C20sublead;          //Energy fraction in the core C60 method 1
   TH2F      *fh2JetCoreMethod3C30sublead;          //Energy fraction in the core C60 method 2
   TH2F      *fh2JetCoreMethod3C60sublead;          //Energy fraction in the core C60 method 3

   TH2F      *fh2SumPtInC10;                  //energy fraction in inner corona C10
   TH2F      *fh2SumPtInC20;                  //energy fraction in inner corona C20 
   TH2F      *fh2SumPtInC30;                  //energy fraction in inner corona C30
   TH2F      *fh2SumPtInC60;                  //energy fraction in inner corona C60
   TH2F      *fh2SumPtInC10lead;              //energy fraction in inner corona C10 leading
   TH2F      *fh2SumPtInC20lead;              //energy fraction in inner corona C20 leading
   TH2F      *fh2SumPtInC30lead;              //energy fraction in inner corona C30 leading
   TH2F      *fh2SumPtInC60lead;              //energy fraction in inner corona C60 leading
   TH2F      *fh2SumPtInC10sublead;           //energy fraction in inner corona C10 subleading
   TH2F      *fh2SumPtInC20sublead;           //energy fraction in inner corona C20 subleading
   TH2F      *fh2SumPtInC30sublead;           //energy fraction in inner corona C30 subleading
   TH2F      *fh2SumPtInC60sublead;           //energy fraction in inner corona C60 subleading
   TH2F      *fh2SumPtOutC10;                 //energy fraction in outer corona C10
   TH2F      *fh2SumPtOutC20;                 //energy fraction in outer corona C20 
   TH2F      *fh2SumPtOutC30;                 //energy fraction in outer corona C30
   TH2F      *fh2SumPtOutC60;                 //energy fraction in outer corona C60
   TH2F      *fh2SumPtOutC10lead;              //energy fraction in outer corona C10 leading
   TH2F      *fh2SumPtOutC20lead;              //energy fraction in outer corona C20 leading
   TH2F      *fh2SumPtOutC30lead;              //energy fraction in outer corona C30 leading
   TH2F      *fh2SumPtOutC60lead;              //energy fraction in outer corona C60 leading
   TH2F      *fh2SumPtOutC10sublead;           //energy fraction in outer corona C10 subleading
   TH2F      *fh2SumPtOutC20sublead;           //energy fraction in outer corona C20 subleading
   TH2F      *fh2SumPtOutC30sublead;           //energy fraction in outer corona C30 subleading
   TH2F      *fh2SumPtOutC60sublead;           //energy fraction in outer corona C60 subleading
   TH2F      *fh2SumPtInC10bkg;        //expected from background inner C10
   TH2F      *fh2SumPtInC20bkg;        //expected from background inner C20
   TH2F      *fh2SumPtInC30bkg;        //expected from background inner C30
   TH2F      *fh2SumPtInC60bkg;        //expected from background inner C60
   TH2F      *fh2SumPtInC10bkglead;        //expected from background inner C10 lead
   TH2F      *fh2SumPtInC20bkglead;        //expected from background inner C20 lead
   TH2F      *fh2SumPtInC30bkglead;        //expected from background inner C30 lead
   TH2F      *fh2SumPtInC60bkglead;        //expected from background inner C60 lead
   TH2F      *fh2SumPtInC10bkgsublead;        //expected from background inner C10 sublead
   TH2F      *fh2SumPtInC20bkgsublead;        //expected from background inner C20 sublead  
   TH2F      *fh2SumPtInC30bkgsublead;        //expected from background inner C30  sublead
   TH2F      *fh2SumPtInC60bkgsublead;        //expected from background inner C60 sublead

   TH2F      *fh2SumPtOutC10bkg;       //expected from background outer C10
   TH2F      *fh2SumPtOutC20bkg;       //expected from background outer C10
   TH2F      *fh2SumPtOutC30bkg;       //expected from background outer C10
   TH2F      *fh2SumPtOutC60bkg;       //expected from background outer C10
   TH2F      *fh2SumPtOutC10bkglead;       //expected from background outer C10 lead
   TH2F      *fh2SumPtOutC20bkglead;       //expected from background outer C10 lead
   TH2F      *fh2SumPtOutC30bkglead;       //expected from background outer C10 lead
   TH2F      *fh2SumPtOutC60bkglead;       //expected from background outer C10 lead
   TH2F      *fh2SumPtOutC10bkgsublead;       //expected from background outer C10 sublead
   TH2F      *fh2SumPtOutC20bkgsublead;       //expected from background outer C10 sublead
   TH2F      *fh2SumPtOutC30bkgsublead;       //expected from background outer C10 sublead
   TH2F      *fh2SumPtOutC60bkgsublead;       //expected from background outer C10 sublead  
   
     TH2F*      fh2DeltaRC10pt1;            //Jet track R distance:C10 pt1
     TH2F*      fh2DeltaRC20pt1;            //C20 pt1 
     TH2F*      fh2DeltaRC30pt1;            //C30 pt1
     TH2F*      fh2DeltaRC60pt1;            //C60 pt1
     TH2F*      fh2DeltaRC10pt2;            //C10 pt2   
     TH2F*      fh2DeltaRC20pt2;            //C20 pt2
     TH2F*      fh2DeltaRC30pt2;            //C30 pt2
     TH2F*      fh2DeltaRC60pt2;            //C60 pt2 
     TH2F*      fh2DeltaRC10pt3;            //C10 pt3
     TH2F*      fh2DeltaRC20pt3;            //C20 pt3
     TH2F*      fh2DeltaRC30pt3;            //C30 pt3
     TH2F*      fh2DeltaRC60pt3;            //C60 pt3
     TH2F*      fh2DeltaRC10pt4;            //C10 pt4
     TH2F*      fh2DeltaRC20pt4;            //C20 pt4
     TH2F*      fh2DeltaRC30pt4;            //C30 pt4
     TH2F*      fh2DeltaRC60pt4;            //C60 pt4 
     TH2F*      fh2DeltaEtaC10pt1;          //The same but eta distance:C10 pt1
     TH2F*      fh2DeltaEtaC20pt1;          //C20 pt1
     TH2F*      fh2DeltaEtaC30pt1;          //C30 pt1
     TH2F*      fh2DeltaEtaC60pt1;          //C60 pt1
     TH2F*      fh2DeltaEtaC10pt2;          //C10 pt2  
     TH2F*      fh2DeltaEtaC20pt2;          //C20 pt2
     TH2F*      fh2DeltaEtaC30pt2;          //C30 pt2
     TH2F*      fh2DeltaEtaC60pt2;          //C60 pt2
     TH2F*      fh2DeltaEtaC10pt3;          //C10 pt3
     TH2F*      fh2DeltaEtaC20pt3;          //C20 pt3
     TH2F*      fh2DeltaEtaC30pt3;          //C30 pt3
     TH2F*      fh2DeltaEtaC60pt3;          //C60 pt3
     TH2F*      fh2DeltaEtaC10pt4;          //C10 pt4
     TH2F*      fh2DeltaEtaC20pt4;          //C20 pt4
     TH2F*      fh2DeltaEtaC30pt4;          //C30 pt4
     TH2F*      fh2DeltaEtaC60pt4;          //C60 pt4
     TH2F*      fh2DeltaPhiC10pt1;          //The same but phi distance:C10 pt1
     TH2F*      fh2DeltaPhiC20pt1;          //C20 pt1
     TH2F*      fh2DeltaPhiC30pt1;          //C30 pt1
     TH2F*      fh2DeltaPhiC60pt1;          //C60 pt1
     TH2F*      fh2DeltaPhiC10pt2;          //C10 pt2
     TH2F*      fh2DeltaPhiC20pt2;          //C20 pt2 
     TH2F*      fh2DeltaPhiC30pt2;          //C30 pt2
     TH2F*      fh2DeltaPhiC60pt2;          //C60 pt2 
     TH2F*      fh2DeltaPhiC10pt3;          //C10 pt3
     TH2F*      fh2DeltaPhiC20pt3;          //C20 pt3
     TH2F*      fh2DeltaPhiC30pt3;          //C30 pt3
     TH2F*      fh2DeltaPhiC60pt3;          //C60 pt3
     TH2F*      fh2DeltaPhiC10pt4;          //C10 pt4
     TH2F*      fh2DeltaPhiC20pt4;          //C20 pt4
     TH2F*      fh2DeltaPhiC30pt4;          //C30 pt4
     TH2F*      fh2DeltaPhiC60pt4;          //C60 pt4

     TH2F*      fh2DeltaRC10pt1lead;            //Jet track R distance:C10 pt1
     TH2F*      fh2DeltaRC20pt1lead;            //C20 pt1 
     TH2F*      fh2DeltaRC30pt1lead;            //C30 pt1
     TH2F*      fh2DeltaRC60pt1lead;            //C60 pt1
     TH2F*      fh2DeltaRC10pt2lead;            //C10 pt2   
     TH2F*      fh2DeltaRC20pt2lead;            //C20 pt2
     TH2F*      fh2DeltaRC30pt2lead;            //C30 pt2
     TH2F*      fh2DeltaRC60pt2lead;            //C60 pt2 
     TH2F*      fh2DeltaRC10pt3lead;            //C10 pt3
     TH2F*      fh2DeltaRC20pt3lead;            //C20 pt3
     TH2F*      fh2DeltaRC30pt3lead;            //C30 pt3
     TH2F*      fh2DeltaRC60pt3lead;            //C60 pt3
     TH2F*      fh2DeltaRC10pt4lead;            //C10 pt4
     TH2F*      fh2DeltaRC20pt4lead;            //C20 pt4
     TH2F*      fh2DeltaRC30pt4lead;            //C30 pt4
     TH2F*      fh2DeltaRC60pt4lead;            //C60 pt4 
     TH2F*      fh2DeltaEtaC10pt1lead;          //The same but eta distance:C10 pt1
     TH2F*      fh2DeltaEtaC20pt1lead;          //C20 pt1
     TH2F*      fh2DeltaEtaC30pt1lead;          //C30 pt1
     TH2F*      fh2DeltaEtaC60pt1lead;          //C60 pt1
     TH2F*      fh2DeltaEtaC10pt2lead;          //C10 pt2  
     TH2F*      fh2DeltaEtaC20pt2lead;          //C20 pt2
     TH2F*      fh2DeltaEtaC30pt2lead;          //C30 pt2
     TH2F*      fh2DeltaEtaC60pt2lead;          //C60 pt2
     TH2F*      fh2DeltaEtaC10pt3lead;          //C10 pt3
     TH2F*      fh2DeltaEtaC20pt3lead;          //C20 pt3
     TH2F*      fh2DeltaEtaC30pt3lead;          //C30 pt3
     TH2F*      fh2DeltaEtaC60pt3lead;          //C60 pt3
     TH2F*      fh2DeltaEtaC10pt4lead;          //C10 pt4
     TH2F*      fh2DeltaEtaC20pt4lead;          //C20 pt4
     TH2F*      fh2DeltaEtaC30pt4lead;          //C30 pt4
     TH2F*      fh2DeltaEtaC60pt4lead;          //C60 pt4
     TH2F*      fh2DeltaPhiC10pt1lead;          //The same but phi distance:C10 pt1
     TH2F*      fh2DeltaPhiC20pt1lead;          //C20 pt1
     TH2F*      fh2DeltaPhiC30pt1lead;          //C30 pt1
     TH2F*      fh2DeltaPhiC60pt1lead;          //C60 pt1
     TH2F*      fh2DeltaPhiC10pt2lead;          //C10 pt2
     TH2F*      fh2DeltaPhiC20pt2lead;          //C20 pt2 
     TH2F*      fh2DeltaPhiC30pt2lead;          //C30 pt2
     TH2F*      fh2DeltaPhiC60pt2lead;          //C60 pt2 
     TH2F*      fh2DeltaPhiC10pt3lead;          //C10 pt3
     TH2F*      fh2DeltaPhiC20pt3lead;          //C20 pt3
     TH2F*      fh2DeltaPhiC30pt3lead;          //C30 pt3
     TH2F*      fh2DeltaPhiC60pt3lead;          //C60 pt3
     TH2F*      fh2DeltaPhiC10pt4lead;          //C10 pt4
     TH2F*      fh2DeltaPhiC20pt4lead;          //C20 pt4
     TH2F*      fh2DeltaPhiC30pt4lead;          //C30 pt4
     TH2F*      fh2DeltaPhiC60pt4lead;          //C60 pt4

     TH2F*      fh2DeltaRC10pt1sublead;            //Jet track R distance:C10 pt1
     TH2F*      fh2DeltaRC20pt1sublead;            //C20 pt1 
     TH2F*      fh2DeltaRC30pt1sublead;            //C30 pt1
     TH2F*      fh2DeltaRC60pt1sublead;            //C60 pt1
     TH2F*      fh2DeltaRC10pt2sublead;            //C10 pt2   
     TH2F*      fh2DeltaRC20pt2sublead;            //C20 pt2
     TH2F*      fh2DeltaRC30pt2sublead;            //C30 pt2
     TH2F*      fh2DeltaRC60pt2sublead;            //C60 pt2 
     TH2F*      fh2DeltaRC10pt3sublead;            //C10 pt3
     TH2F*      fh2DeltaRC20pt3sublead;            //C20 pt3
     TH2F*      fh2DeltaRC30pt3sublead;            //C30 pt3
     TH2F*      fh2DeltaRC60pt3sublead;            //C60 pt3
     TH2F*      fh2DeltaRC10pt4sublead;            //C10 pt4
     TH2F*      fh2DeltaRC20pt4sublead;            //C20 pt4
     TH2F*      fh2DeltaRC30pt4sublead;            //C30 pt4
     TH2F*      fh2DeltaRC60pt4sublead;            //C60 pt4 
     TH2F*      fh2DeltaEtaC10pt1sublead;          //The same but eta distance:C10 pt1
     TH2F*      fh2DeltaEtaC20pt1sublead;          //C20 pt1
     TH2F*      fh2DeltaEtaC30pt1sublead;          //C30 pt1
     TH2F*      fh2DeltaEtaC60pt1sublead;          //C60 pt1
     TH2F*      fh2DeltaEtaC10pt2sublead;          //C10 pt2  
     TH2F*      fh2DeltaEtaC20pt2sublead;          //C20 pt2
     TH2F*      fh2DeltaEtaC30pt2sublead;          //C30 pt2
     TH2F*      fh2DeltaEtaC60pt2sublead;          //C60 pt2
     TH2F*      fh2DeltaEtaC10pt3sublead;          //C10 pt3
     TH2F*      fh2DeltaEtaC20pt3sublead;          //C20 pt3
     TH2F*      fh2DeltaEtaC30pt3sublead;          //C30 pt3
     TH2F*      fh2DeltaEtaC60pt3sublead;          //C60 pt3
     TH2F*      fh2DeltaEtaC10pt4sublead;          //C10 pt4
     TH2F*      fh2DeltaEtaC20pt4sublead;          //C20 pt4
     TH2F*      fh2DeltaEtaC30pt4sublead;          //C30 pt4
     TH2F*      fh2DeltaEtaC60pt4sublead;          //C60 pt4
     TH2F*      fh2DeltaPhiC10pt1sublead;          //The same but phi distance:C10 pt1
     TH2F*      fh2DeltaPhiC20pt1sublead;          //C20 pt1
     TH2F*      fh2DeltaPhiC30pt1sublead;          //C30 pt1
     TH2F*      fh2DeltaPhiC60pt1sublead;          //C60 pt1
     TH2F*      fh2DeltaPhiC10pt2sublead;          //C10 pt2
     TH2F*      fh2DeltaPhiC20pt2sublead;          //C20 pt2 
     TH2F*      fh2DeltaPhiC30pt2sublead;          //C30 pt2
     TH2F*      fh2DeltaPhiC60pt2sublead;          //C60 pt2 
     TH2F*      fh2DeltaPhiC10pt3sublead;          //C10 pt3
     TH2F*      fh2DeltaPhiC20pt3sublead;          //C20 pt3
     TH2F*      fh2DeltaPhiC30pt3sublead;          //C30 pt3
     TH2F*      fh2DeltaPhiC60pt3sublead;          //C60 pt3
     TH2F*      fh2DeltaPhiC10pt4sublead;          //C10 pt4
     TH2F*      fh2DeltaPhiC20pt4sublead;          //C20 pt4
     TH2F*      fh2DeltaPhiC30pt4sublead;          //C30 pt4
     TH2F*      fh2DeltaPhiC60pt4sublead;          //C60 pt4





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




   AliAnalysisTaskJetCore(const AliAnalysisTaskJetCore&); // not implemented
   AliAnalysisTaskJetCore& operator=(const AliAnalysisTaskJetCore&); // not implemented

   ClassDef(AliAnalysisTaskJetCore, 4);
};

#endif

