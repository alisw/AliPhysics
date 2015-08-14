#ifndef ALIANALYSISTASKJETCOREPP_H
#define ALIANALYSISTASKJETCOREPP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// This task performs hadron-trigger recoil jet correlations 
// Output pT spectrum of jet given trigger pT 
// Author: filip krizek 16th March 2013
// *******************************************

class TH1F;
class TH1D;
class TH1I;
class TH2F;
class TH3F;
class TList;
class TClonesArray;
class THnSparse;
class TRandom3;
class TArrayI; 
class TProfile;
class TFile;
class TKey;
class AliESDEvent;
class AliAODExtension;
class AliAODEvent;
class AliGenPythiaEventHeader;
class AliMCEvent;    //FK//
class AliMCEventHandler; //FK//
class AliGenEventHeader; //FK//

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
   virtual void  SetBranchNameChargMC(const TString &name){ fJetBranchNameChargMC = name; } 
   virtual void  SetBranchNameKine(const TString &name){ fJetBranchNameKine = name; } 
   virtual void  SetBranchNameFullMC(const TString &name){ fJetBranchNameFullMC = name; } 
   virtual void  SetBranchNameBg(const TString &name){ fJetBranchNameBg = name; } 
   virtual void  SetBranchNameBgChargMC(const TString &name){ fJetBranchNameBgChargMC = name; } 
   virtual void  SetBranchNameBgKine(const TString &name){ fJetBranchNameBgKine = name; } 
   virtual void  SetNonStdFile(char* c){fNonStdFile = c;} 
   virtual void  SetSystem(Int_t sys) { fSystem = sys; } 
   virtual void  SetJetR(Float_t jR) { fJetParamR = jR; }
   virtual void  SetBgJetR(Float_t bgjR) { fBgJetParamR = bgjR; }
   virtual void  SetBgMaxJetPt(Float_t mpt){ fBgMaxJetPt = mpt;}
   virtual void  SetRndTrials(Int_t nt){ fnTrials = nt;}
   virtual void  SetFreeAreaFrac(Float_t frac){ fJetFreeAreaFrac = frac;}
   virtual void  SetBgConeR(Float_t cr){ fBgConeR = cr; } 
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
   virtual void  SetTriggerType(Int_t tt){ fHardest=tt;}
   virtual void  SetEventNumberRangeLow(Int_t rl){ fEventNumberRangeLow=rl;}
   virtual void  SetEventNumberRangeHigh(Int_t rh){ fEventNumberRangeHigh=rh;}  
   virtual void  SetTriggerPtRangeLow(Float_t tl){ fTriggerPtRangeLow=tl;}   
   virtual void  SetTriggerPtRangeHigh(Float_t th){ fTriggerPtRangeHigh=th;}  
   virtual void  SetFillResponseMatrix(Bool_t brm){ fFillRespMx = brm;}
   virtual void  SetBinning(Bool_t bbb) { fDoubleBinning = bbb; } 
   virtual void  SetUseExchangeContainerInput(Bool_t b){ fUseExchContainer = b;} 

   Double_t RelativePhi(Double_t angle1, Double_t angle2); 

private:
   //private member functions
   Int_t   GetListOfTracks(TList *list); //returns index of trig and track list

   Bool_t SelectMCGenTracks(AliVParticle *trk, TList *trkList, Double_t &ptLeading, Int_t &index, Int_t counter);
   void FillEffHistos(TList *recList, TList *genList);
   
   void EstimateBgRhoMedian(TList *listJet, TList* listPart, Double_t &rhoMedian, Int_t mode);//median method to estimate bg
   void EstimateBgCone(TList *listJet, TList* listPart, Double_t &rhoPerpCone);//perp cone method to estimate bg
   void ReadTClonesArray(TString bname, TList *list); //init jets lists

   void FillDeltaPt(TList *jetList, TList *trkList, Double_t rhoMedian, Double_t rhoCones);

   //private member objects
   AliESDEvent *fESD;    //! ESD object
   AliAODEvent *fAODIn;  //! AOD event for AOD input tracks
   AliAODEvent *fAODOut; //! AOD event 
   AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
   AliMCEvent           *fMcEvent;    //! MC event                       
   AliInputEventHandler *fMcHandler;  //! MCEventHandler                 


   // jets to compare
   TString fJetBranchName; //  name of jet branch 
   TString fJetBranchNameChargMC;   //  name of jet branch output AOD
   TString fJetBranchNameKine; //  name of jet branch kine
   TString fJetBranchNameFullMC; //  name of jet branch 
   TString fJetBranchNameBg; //  name of bg (kt) jet branch 
   TString fJetBranchNameBgChargMC; //  name of bg (kT) jet branch 
   TString fJetBranchNameBgKine; //  name of bg (kT) jet branch 
   TList  *fListJets;      //! jet list reconstructed level
   TList  *fListJetsGen;   //! jet list generator level 
   TList  *fListJetsGenFull; //! jet list generator level full jets 
   TList  *fListJetsBg;      //! jet list reconstructed level to be removed from bg
   TList  *fListJetsBgGen;   //! jet list generator level to be removed from bg  


   TString fNonStdFile;    // name of delta aod file to catch the extension

   // event selection
   Int_t   fSystem;        // collision system  pp=0, pPb=1  
   Float_t fJetParamR;     // jet cone resolution (radius) R 
   Float_t fBgJetParamR;   // jet cone resolution (radius) R of jet to be removed from bg
   Float_t fBgMaxJetPt;    // max pt of jets accepted in bg 
   Float_t fBgConeR;       //perp cone R used to assess bg
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
   Bool_t  fUseExchContainer; //use exhange container
   
   TList *fOutputList;          //! output data container 
   TH1I  *fHistEvtSelection;    //! event selection statistic 
   TH2F      *fh2Ntriggers;     //trigger pT versus centrality 
   THnSparse *fHJetSpec;      //Recoil jet spectrum  
   THnSparse *fHJetSpecSubUeMedian; //Recoil jet spectrum, jet pT corrected by kT median  
   THnSparse *fHJetSpecSubUeCone;  //Recoil jet spectrum, jet pT corrected by perp cone rho 
   
   THnSparse *fHJetPhiCorr; // Dphi distribution jet-triger
   THnSparse *fHJetPhiCorrSubUeMedian; // Dphi distribution jet-triger
   THnSparse *fHJetPhiCorrSubUeCone; // Dphi distribution jet-triger

   //Diagnostics
   THnSparse *fHJetUeMedian;   //UE background from kT median
   THnSparse *fHJetUeCone;      //UE background from perp cone 
   THnSparse *fHRhoUeMedianVsCone;    //EBE UE from perp cone
   //THnSparse *fHJetDensity;       //density of jet with A>0.07  //fk
   //THnSparse *fHJetDensityA4;     //density of jets with A>0.4 //fk
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
   TH1D *fhNofMultipleTriggers; // The number of additional triggers in events with at least one trigger 
   TH1D *fhNofMultipleTriggersCone; // The number of additional triggers in events with at least one trigger 
   TH1D *fhNofMultipleTriggersConeLow; // The number of additional triggers in events with at least one trigger 
   TH1D *fhNofMultipleTriggersConeHigh; // The number of additional triggers in events with at least one trigger 
   TH1D *fhDeltaRMultTriggersLow; // Angular distributions between trigger and assoc
   TH1D *fhDeltaRMultTriggersHigh; // Angular distributions between trigger and assoc
   TH1D *fhDeltaPhiMultTriggersLow; // Delta phi between trigger and assoc single incl trigger
   TH1D *fhDeltaPhiMultTriggersHigh; // Delta phi between trigger and assoc single incl triger

   TH1D *fhDeltaPhiMultTriggersInclLow; // Delta phi between trigger and assoc  incl trigg
   TH1D *fhDeltaPhiMultTriggersInclHigh; // Delta phi between trigger and assoc incl trigg
   TH1D *fhInclTrigCounter; // count the total number of inclusive triggers

   TH1D *fhDeltaPtConeBg;  // delta pt of bg fluctuations with rho based on cone bg
   TH1D *fhDeltaPtMedianBg;// delta pt of bg fluctuations whith rho based on bg median
   //THnSparse *fHJetPtRaw;      //bg unsubtr. vs bg subtr. pT spectrum of jets vs jet area
   //THnSparse *fHLeadingJetPtRaw; //bg unsubtr. vs bg. subtr. leading jet pT vs area 
   //THnSparse *fHDphiVsJetPtAll;   //Dphitrigger-jet  versus jet pt for all jets given pTtrigg  

   //MC generator level
   TH2D      *fhJetPtGenVsJetPtRec; //jet respose matrix  
   TH2D      *fhJetPtGenVsJetPtRecSubUeMedian; //jet respose matrix both pT with subtracted kT median bg 
   TH2D      *fhJetPtGenVsJetPtRecSubUeCone; //jet respose matrix both pT with subtracted weighted kT median bg 
   TH1D      *fhJetPtGen;           //generated pT spectrum of jets  
   TH1D      *fhJetPtSubUeMedianGen; //generated pT spectrum of jets with subtracted kT median  
   TH1D      *fhJetPtSubUeConeGen;    //generated pT spectrum of jets with perp cone
   TH2D      *fhJetPtResolutionVsPtGen; // pTjet,rec-pTjet,gen/ pTjet,gen  versus pT jet,gen
   TH2D      *fhJetPtResolutionVsPtConeGen;//pTjet,rec-pTjet,gen/ pTjet,gen  versus pT jet,gen
   TH2D      *fhJetPtGenChargVsJetPtGenFull; //generated pT spectrum of full jets
   TH1D      *fhJetPtGenFull; // generated pT spectrum of full jets
   TH2F      *fh2NtriggersGen; //trigger pT versus centrality in generator level
   THnSparse *fHJetSpecGen;    //Recoil jet spectrum at generator level 
   THnSparse *fHJetSpecSubUeMedianGen;  //Recoil jet spectrum at gen level, jet pT corrected by kT median 
   THnSparse *fHJetSpecSubUeConeGen; //Recoil jet spectrum at gen level, jet pT corrected with rho from cone
   THnSparse *fHJetPhiCorrGen; // Dphi distribution jet-triger
   THnSparse *fHJetPhiCorrSubUeMedianGen; // Dphi distribution jet-triger
   THnSparse *fHJetPhiCorrSubUeConeGen; // Dphi distribution jet-triger
   THnSparse *fHJetUeMedianGen;   //UE background from kT median
   THnSparse *fHJetUeConeGen;      //UE background from Perp Cone 
   TH2D      *fhPtTrkTruePrimRec; // pt spectrum of true reconstructed primary tracks    
   TH2D      *fhPtTrkTruePrimGen; // pt spectrum of true generated primary track    
   TH2D      *fhPtTrkSecOrFakeRec; // pt spectrum of reconstructed fake or secondary tracks    
   THnSparse *fHRhoUeMedianVsConeGen; //EBE UE from Median vs Perp Cone  generator level 
  
   TH1D  *fhEntriesToMedian; //how many entries were used to calculate
   TH1D  *fhEntriesToMedianGen; //how many entries were used to calculate in MC
   TH1D  *fhCellAreaToMedian; //how many entries were used to calculate
   TH1D  *fhCellAreaToMedianGen; //how many entries were used to calculate in MC
 
   TH1D *fhNofMultipleTriggersGen; // The number of additional triggers in events with at least one trigger 
   TH1D *fhNofMultipleTriggersConeGen; // The number of additional triggers in events with at least one trigger in R<0.4 15-50
   TH1D *fhNofMultipleTriggersConeGenLow; // The number of additional triggers in events with at least one trigger in R<0.4 15-20
   TH1D *fhNofMultipleTriggersConeGenHigh; // The number of additional triggers in events with at least one trigger in R<0.4 20-50 
   TH1D *fhDeltaRMultTriggersGenLow; // Angular distributions between trigger and assoc
   TH1D *fhDeltaRMultTriggersGenHigh; // Angular distributions between trigger and assoc
   TH1D *fhDeltaPhiMultTriggersGenLow; // Angular distributions between trigger and assoc  15-20
   TH1D *fhDeltaPhiMultTriggersGenHigh; // Angular distributions between trigger and assoc 20-50

   TH1D *fhDeltaPhiMultTriggersInclGenLow; // Delta phi between trigger and assoc  incl trigg
   TH1D *fhDeltaPhiMultTriggersInclGenHigh; // Delta phi between trigger and assoc incl trigg
   TH1D *fhInclTrigCounterGen; // count the total number of inclusive triggers

   TH1D *fhNofMultipleTriggersConeGenA; // The number of additional triggers in events with at least one trigger in R<0.4 
   TH1D *fhNofMultipleTriggersConeGenALow; // The number of additional triggers in events with at least one trigger in R<0.4 15-20 
   TH1D *fhNofMultipleTriggersConeGenAHigh; // The number of additional triggers in events with at least one trigger in R<0.4  20-50
   TH1D *fhDeltaRMultTriggersGenALow; //Delta R for eloss scenation in assoc 15-20 bin 
   TH1D *fhDeltaPhiMultTriggersGenALow;//Delta phi for eloss scenation in assoc 15-20 bin 
   TH1D *fhDeltaRMultTriggersGenAHigh; //Delta R for eloss scenation in assoc 20-50 bin 
   TH1D *fhDeltaPhiMultTriggersGenAHigh;//Delta phi for eloss scenation in assoc 20-50 bin 
   TH1D *fhNofTriggersGenA; //Count triggers in eloss scenario

   TH1D *fhNofMultipleTriggersConeGenB; // The number of additional triggers in events with at least one trigger in R<0.4 15-50
   TH1D *fhNofMultipleTriggersConeGenBLow; // The number of additional triggers in events with at least one trigger in R<0.4 15-20
   TH1D *fhNofMultipleTriggersConeGenBHigh; // The number of additional triggers in events with at least one trigger in R<0.4 20-50 
   TH1D *fhDeltaRMultTriggersGenBLow; //Delta R for eloss scenation in assoc 15-20 bin 
   TH1D *fhDeltaPhiMultTriggersGenBLow;//Delta phi for eloss scenation in assoc 15-20 bin 
   TH1D *fhDeltaRMultTriggersGenBHigh; //Delta R for eloss scenation in assoc 20-50 bin 
   TH1D *fhDeltaPhiMultTriggersGenBHigh;//Delta phi for eloss scenation in assoc 20-50 bin 
   TH1D *fhNofTriggersGenB; //Count triggers in eloss scenario


   TH1D *fhTriggerCounterGenLevel; // The number of sing incl rec TT tracks that have gen level particle assigned within the same TT bin 
   TH1D *fhDeltaRMultTriggersGenLevelLow; // corresp. genereator level TT track combined with generator level Assoc tracks 15-20 
   TH1D *fhDeltaPhiMultTriggersGenLevelLow; //corresp. genereator level TT track combined with generator level Assoc tracks 15-20 
   TH1D *fhDeltaRMultTriggersGenLevelHigh;  // corresp. genereator level TT track combined with generator level Assoc tracks 20-50 
   TH1D *fhDeltaPhiMultTriggersGenLevelHigh;// corresp. genereator level TT track combined with generator level Assoc tracks 20-50

   Bool_t fIsChargedMC;   //flag analysis on MC data with true and on the real+kine data false
   Bool_t fIsKine;       //flag analysis on kine data with true and on the real+MC data false
   Bool_t fIsFullMC;   //flag analysis on MC data with true and on the real+kine data false
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

   
   Int_t   fHardest;               // trigger type 0=single incl, 1=LP 
   Int_t   fEventNumberRangeLow;   // lower range of selected event numbers  
   Int_t   fEventNumberRangeHigh;  // high range of selected event numbers  
   Float_t fTriggerPtRangeLow;   // lower range of selected trigger pt
   Float_t fTriggerPtRangeHigh;  // upper range of selected trigger pt

   Bool_t  fFillRespMx;    //fill response matrix files


   TRandom3* fRandom;           // TRandom3 
   Int_t fnTrials;  //number of random trials to measure cell area
   Float_t fJetFreeAreaFrac; //fraction of area in cell free of jets  
   const Int_t  fnPhi; //number of cells in phi
   const Int_t  fnEta; //number of cells in eta
   const Double_t fEtaSize; //cell size in eta 
   const Double_t fPhiSize; //cell size in phi
   const Double_t fCellArea; //cell area
   Double_t fSafetyMargin; //enlarge a bit the jet size to avoid contamination of UE

   Bool_t fDoubleBinning; //0=use 2 GeV/c bins  ; 1= use 1 GeV/c bins
 
   ClassDef(AliAnalysisTaskJetCorePP, 18);  //has to end with number larger than 0
};

#endif

