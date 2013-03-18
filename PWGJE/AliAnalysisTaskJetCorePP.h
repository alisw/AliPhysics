#ifndef ALIANALYSISTASKJETCOREPP_H
#define ALIANALYSISTASKJETCOREPP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// This task performs hadron-trigger recoil jet correlations 
// Output pT spectrum of jet given trigger pT 
// Author: filip krizek 1st March 2013
// *******************************************

class TH1F;
class TH1D;
class TH1I;
class TH2F;
class TH3F;
class THnSparse;
class AliESDEvent;
class AliAODExtension;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"

class AliAnalysisTaskJetCorePP : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskJetCorePP();
   AliAnalysisTaskJetCorePP(const char *name);
   AliAnalysisTaskJetCorePP(const AliAnalysisTaskJetCorePP& a); 
   AliAnalysisTaskJetCorePP& operator=(const AliAnalysisTaskJetCorePP& a); // not implemented
   virtual ~AliAnalysisTaskJetCorePP();
   virtual void  LocalInit() {Init();}
   virtual void  Init();
   virtual void  UserCreateOutputObjects();
   virtual void  UserExec(Option_t *option);
   virtual void  Terminate(const Option_t*);
 
   virtual void  SetBranchName(const TString &name){ fJetBranchName = name; } 
   virtual void  SetNonStdFile(char* c){fNonStdFile = c;} 
   virtual void  SetSystem(Int_t sys) { fSystem = sys; } 
   virtual void  SetJetR(Float_t jR) { fJetParamR = jR; }
   virtual void  SetOfflineTrgMask(AliVEvent::EOfflineTriggerTypes mask) { fOfflineTrgMask = mask; } 
   virtual void  SetMinContribVtx(Int_t n) { fMinContribVtx = n; } 
   virtual void  SetVtxZMin(Float_t z) { fVtxZMin = z; }
   virtual void  SetVtxZMax(Float_t z) { fVtxZMax = z; } 
   virtual void  SetFilterMask(UInt_t i){fFilterMask = i;} 
   virtual void  SetCentMin(Float_t cent) { fCentMin = cent; }
   virtual void  SetCentMax(Float_t cent) { fCentMax = cent; } 
   virtual void  SetJetEtaMin(Float_t eta) { fJetEtaMin = eta; }
   virtual void  SetJetEtaMax(Float_t eta) { fJetEtaMax = eta; } 
   virtual void  SetTriggerEtaCut(Float_t eta) { fTriggerEtaCut = eta; }
   virtual void  SetTrackEtaCut(Float_t eta) { fTrackEtaCut = eta; }
   virtual void  SetTrackLowPtCut(Float_t pt) { fTrackLowPtCut=pt; } 

   Double_t RelativePhi(Double_t angle1, Double_t angle2); 

private:
   //private member functions
   Int_t   GetListOfTracks(TList *list); //returns index of trig and track list 
   Double_t GetBackgroundInPerpCone(Float_t jetR, Double_t jetPhi, Double_t jetEta, TList* trkList); //sums pT in the cone perp in phi to jet
 
   //private member objects
   AliESDEvent *fESD;    //! ESD object
   AliAODEvent *fAODIn;  //! AOD event for AOD input tracks
   AliAODEvent *fAODOut; //! AOD event 
   AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD

   // jets to compare
   TString fJetBranchName; //  name of jet branch 
   TList  *fListJets;      //! jet lists  

   TString fNonStdFile;    // name of delta aod file to catch the extension

   // event selection
   Int_t   fSystem;        // collision system  pp=0, pPb=1  
   Float_t fJetParamR;     // jet cone resolution (radius) R 
   AliVEvent::EOfflineTriggerTypes fOfflineTrgMask; // mask of offline trigs 
   Int_t   fMinContribVtx; // min numb of trk contrib for prim vertex 
   Float_t fVtxZMin;	   // lower bound on vertex z 
   Float_t fVtxZMax;	   // upper bound on vertex z 
   UInt_t  fFilterMask;    // filter bit for slected tracks  
   Float_t fCentMin;	   // lower bound on centrality 
   Float_t fCentMax;	   // upper bound on centrality 
   Float_t fJetEtaMin;     // lower bound on eta for found jets 
   Float_t fJetEtaMax;     // upper bound on eta for found jets 
   Float_t fTriggerEtaCut; // lower bound on eta for trigger track
   Float_t fTrackEtaCut;   // upper bound on eta for trigger track 
   Float_t fTrackLowPtCut; // upper bound on eta for trigger track
   
   
   TList *fOutputList;          //! output data container 
   TH1I  *fHistEvtSelection;    //! event selection statistic 
   TH2F      *fh2Ntriggers;     //trigger pT versus centrality 
   THnSparse *fHJetSpec;      //Recoil jet spectrum  
   
   //Diagnostics
   THnSparse *fHJetDensity;       //density of jet with A>0.07  //fk
   THnSparse *fHJetDensityA4;     //density of jets with A>0.4 //fk
   TH2D *fhJetPhi;     //Azimuthal distribution of jets
   TH2D *fhTriggerPhi; //Azimuthal distribution of trigger hadron
   TH2D *fhJetEta;     //Pseudorapidity distribution of jets
   TH2D *fhTriggerEta; //Pseudorapidity distribution of trigger hadron
   TH1D *fhVertexZ;    //z vertex distribution 
   TH1D *fhVertexZAccept;    //z vertex distribution after cut
   TH1D *fhContribVtx;    //contributors to vertex 
   TH1D *fhContribVtxAccept;    //contributors to vertex after cut
   TH1D *fhDphiTriggerJet;  //Deltaphi between trigger and jet 
   TH1D *fhDphiTriggerJetAccept;  //Deltaphi between trigger and jet after cut
   TH1D *fhCentrality;  //Deltaphi between trigger and jet 
   TH1D *fhCentralityAccept;  //Deltaphi between trigger and jet after cut

   THnSparse *fHJetPtRaw;      //bg unsubtr. vs bg subtr. pT spectrum of jets vs jet area
   THnSparse *fHLeadingJetPtRaw; //bg unsubtr. vs bg. subtr. leading jet pT vs area 
   THnSparse *fHDphiVsJetPtAll;   //Dphitrigger-jet  versus jet pt for all jets given pTtrigg  
   THnSparse *fHRhoFastJetVsRhoCone; //fast jet rho vs perp cone rho given pT jet

   const Double_t fkAcceptance; //eta times phi  Alice coverage  
   Double_t fConeArea;      //cone area pi*R^2
  
   ClassDef(AliAnalysisTaskJetCorePP, 2);  //has to end with number larger than 0
};

#endif

