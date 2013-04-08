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
class TList;
class THnSparse;
class TArrayI; 
class TProfile;
class TFile;
class TKey;
class AliESDEvent;
class AliAODExtension;
class AliAODEvent;
class AliGenPythiaEventHeader;

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
   virtual Bool_t Notify();
 
   virtual void  SetBranchName(const TString &name){ fJetBranchName = name; } 
   virtual void  SetBranchNameMC(const TString &name){ fJetBranchNameMC = name; } 
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
   //Double_t GetBackgroundInPerpCone(Float_t jetR, Double_t jetPhi, Double_t jetEta, TList* trkList); //sums pT in the cone perp in phi to jet
   Bool_t SelectMCGenTracks(AliVParticle *trk, TList *trkList, Double_t &ptLeading, Int_t &index, Int_t counter);
   void FillEffHistos(TList *recList, TList *genList);

   //private member objects
   AliESDEvent *fESD;    //! ESD object
   AliAODEvent *fAODIn;  //! AOD event for AOD input tracks
   AliAODEvent *fAODOut; //! AOD event 
   AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD

   // jets to compare
   TString fJetBranchName; //  name of jet branch 
   TString fJetBranchNameMC; //  name of jet branch 
   TList  *fListJets;      //! jet list reconstructed level
   TList  *fListJetsGen;   //! jet list generator level  

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

   //MC generator level
   TH2D      *fhJetPtGenVsJetPtRec; //jet respose matrix  
   TH1D      *fhJetPtGen;           //generated pT spectrum of jets  
   TH2F      *fh2NtriggersGen; //trigger pT versus centrality in generator level
   THnSparse *fHJetSpecGen;    //Recoil jet spectrum in generator level 
   TH2D      *fhPtTrkTruePrimRec; // pt spectrum of true reconstructed primary tracks    
   TH2D      *fhPtTrkTruePrimGen; // pt spectrum of true generated primary track    
   TH2D      *fhPtTrkSecOrFakeRec; // pt spectrum of reconstructed fake or secondary tracks    
   
   Bool_t fIsMC;   //flag analysis on MC data with true and on the real data false
   TArrayI faGenIndex;   // labels of particles on MC generator level  
   TArrayI faRecIndex;   // labels of particles on reconstructed track level
   const Double_t fkAcceptance; //eta times phi  Alice coverage  
   const Double_t fkDeltaPhiCut; //Delta phi cut on  trigger-jet distance in azimuth
 
   TProfile*     fh1Xsec;   //! pythia cross section and trials
   TH1F*         fh1Trials; //! trials are added
   TH1F*         fh1AvgTrials; //! trials are added
   TH1F*         fh1PtHard;  //! Pt har of the event...      
   TH1F*         fh1PtHardNoW;  //! Pt har of the event without weigt      
   TH1F*         fh1PtHardTrials;  //! Number of trials
   Float_t       fAvgTrials;       // Average number of trials

   ClassDef(AliAnalysisTaskJetCorePP, 4);  //has to end with number larger than 0
};

#endif

