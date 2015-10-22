#ifndef ALIANALYSISTASKSELC2ELELAMBDAFROMAODTRACKS_H
#define ALIANALYSISTASKSELC2ELELAMBDAFROMAODTRACKS_H

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */ 

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliPID.h"
#include "AliRDHFCutsLctoeleLambdafromAODtracks.h"

/// \class AliAnalysisTaskSELc2eleLambdafromAODtracks

class THnSparse;
class TH1F;
class TClonesArray;
class AliAODRecoCascadeHF;
class AliESDVertex;
class AliAODMCParticle;
class AliEventPoolManager;
class AliNormalizationCounter;

class AliAnalysisTaskSELc2eleLambdafromAODtracks : public AliAnalysisTaskSE 
{
 public:
  AliAnalysisTaskSELc2eleLambdafromAODtracks();
  AliAnalysisTaskSELc2eleLambdafromAODtracks(const Char_t* name, AliRDHFCutsLctoeleLambdafromAODtracks* cuts, Bool_t writeVariableTree=kTRUE);
  virtual ~AliAnalysisTaskSELc2eleLambdafromAODtracks();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void FillROOTObjects(AliAODRecoCascadeHF *elobj, AliAODv0 *v0, AliAODTrack *trk, TClonesArray *mcArray);
  void FillMixROOTObjects(TLorentzVector *et, TLorentzVector *ev, Double_t *v0info, Int_t charge);
  void FillElectronROOTObjects(AliAODTrack *trk, TClonesArray *mcArray);
  void FillV0ROOTObjects(AliAODv0 *v0, TClonesArray *mcArray);
  void FillMCROOTObjects(AliAODMCParticle *part, AliAODMCParticle *mcepart, AliAODMCParticle *mcv0part, Int_t decaytype);
  void FillMCEleROOTObjects(AliAODMCParticle *mcepart, TClonesArray *mcArray);
  void FillMCV0ROOTObjects(AliAODMCParticle *mcv0part, TClonesArray *mcArray);
  Bool_t MakeMCAnalysis(TClonesArray *mcArray);
  void MakeAnalysis(AliAODEvent *aod, TClonesArray *mcArray);

  void SelectV0( const AliVEvent *event,Int_t nV0,Int_t &nSeleV0, Bool_t *seleV0Flags, TClonesArray *mcArray);
  void SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags, TClonesArray *mcArray);

  /// set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}
  void SetUseCentralityV0M(Bool_t centon) {fUseCentralityV0M = centon;}
  Bool_t GetUseCentralityV0M() const {return fUseCentralityV0M;}
  void SetWriteEachVariableTree(Bool_t a) {fWriteEachVariableTree = a;}
  Bool_t GetWriteEachVariableTree() const {return fWriteEachVariableTree;}
  void SetWriteMCVariableTree(Bool_t a) {fWriteMCVariableTree = a;}
  Bool_t GetWriteMCVariableTree() const {return fWriteMCVariableTree;}
  void SetMCEventType(Int_t theevt) {fMCEventType = theevt;}
  Bool_t GetMCEventType() const {return fMCEventType;}

  void SetReconstructPrimVert(Bool_t a) { fReconstructPrimVert=a; }

  AliAODRecoCascadeHF* MakeCascadeHF(AliAODv0 *casc, AliAODTrack *trk, AliAODEvent *aod, AliAODVertex *vert);
  AliAODVertex* ReconstructSecondaryVertex(AliAODv0 *casc, AliAODTrack *trk, AliAODEvent *aod);
	Int_t MatchToMC(AliAODRecoCascadeHF *elobj, TClonesArray *mcArray, Int_t *pdgele_array, Int_t *pdgv0_array, Int_t *labelele_array, Int_t *labelv0_array,  Int_t &ngen_ele, Int_t &ngen_v0);


  /// mixing
  void SetEventMixingWithPools(){fDoEventMixing=1;}
  void SetEventMixingOff(){fDoEventMixing=0;}
	void SetNumberOfEventsForMixing(Int_t events){fNumberOfEventsForMixing=events;}
	void SetPoolPVzBinLimits(Int_t Nzvtxbins,const Double_t *ZvtxBins){
		fNzVtxBins = Nzvtxbins;
		for(int ix = 0;ix<fNzVtxBins+1;ix++){fZvtxBins[ix] = ZvtxBins[ix];}
	}
	void SetPoolCentBinLimits(Int_t Ncentbins,const Double_t *CentBins){
		fNCentBins = Ncentbins;
		for(int ix = 0;ix<fNCentBins+1;ix++){fCentBins[ix] = CentBins[ix];}
	}
  void DoEventMixingWithPools(Int_t index);
  void ResetPool(Int_t poolIndex);
  Int_t GetPoolIndex(Double_t zvert, Double_t mult);


 private:

  AliAnalysisTaskSELc2eleLambdafromAODtracks(const AliAnalysisTaskSELc2eleLambdafromAODtracks &source);
  AliAnalysisTaskSELc2eleLambdafromAODtracks& operator=(const AliAnalysisTaskSELc2eleLambdafromAODtracks& source); 

  void DefineTreeVariables();
  void DefineEleTreeVariables();
  void DefineV0TreeVariables();
  void DefineMCTreeVariables();
  void DefineMCEleTreeVariables();
  void DefineMCV0TreeVariables();
  void DefineGeneralHistograms();
  void DefineAnalysisHistograms();

  AliAODVertex *CallPrimaryVertex(AliAODv0 *v0, AliAODTrack *trk, AliAODEvent *evt);
  AliAODVertex* PrimaryVertex(const TObjArray *trkArray,AliVEvent *event);

  Bool_t fUseMCInfo;                 /// Use MC info
  TList *fOutput;                    //!<! User output slot 1 // general histos
  TList *fOutputAll;                 //!<! User Output slot 3  //analysis histograms 
  TList *fListCuts;                  //!<! User output slot 2 // Cuts 
  TH1F *fCEvents;                    //!<! Histogram to check selected events
  TH1F *fHTrigger;                   //!<! Histogram to check Trigger
  TH1F *fHCentrality;                //!<! Histogram to check Centrality
  AliRDHFCutsLctoeleLambdafromAODtracks *fAnalCuts;// Cuts - sent to output slot 2
  Bool_t fIsEventSelected;          /// flag for event selected
  Bool_t    fWriteVariableTree;     /// flag to decide whether to write the candidate variables on a tree variables
  Bool_t    fWriteEachVariableTree;     /// flag to decide whether to write the candidate variables on a tree variables
  Bool_t    fWriteMCVariableTree;     /// flag to decide whether to write the candidate variables on a tree variables
  TTree    *fVariablesTree;         //!<! tree of the candidate variables after track selection on output slot 4
  TTree    *fEleVariablesTree;         //!<! tree of the candidate variables after track selection on output slot 4
  TTree    *fV0VariablesTree;         //!<! tree of the candidate variables after track selection on output slot 4
  TTree    *fMCVariablesTree;         //!<! tree of the candidate variables after track selection on output slot 4
  TTree    *fMCEleVariablesTree;         //!<! tree of the candidate variables after track selection on output slot 4
  TTree    *fMCV0VariablesTree;         //!<! tree of the candidate variables after track selection on output slot 4
  Bool_t fReconstructPrimVert;       ///Reconstruct primary vertex excluding candidate tracks
  Bool_t fIsMB;       /// MB trigger event
  Bool_t fIsSemi;     /// SemiCentral trigger event
  Bool_t fIsCent;     /// Central trigger event
  Bool_t fIsINT7;     /// INT7 trigger event
  Bool_t fIsEMC7;     /// EMC7 trigger event
  Float_t *fCandidateVariables;   //!<! variables to be written to the tree
  Float_t *fCandidateEleVariables;   //!<! variables to be written to the tree
  Float_t *fCandidateV0Variables;   //!<! variables to be written to the tree
  Float_t *fCandidateMCVariables;   //!<! variables to be written to the tree
  Float_t *fCandidateMCEleVariables;   //!<! variables to be written to the tree
  Float_t *fCandidateMCV0Variables;   //!<! variables to be written to the tree
  AliAODVertex *fVtx1;            /// primary vertex
  AliESDVertex *fV1;              /// primary vertex
  Float_t  fVtxZ;         /// zVertex
  Double_t fBzkG;                 /// magnetic field value [kG]
  Float_t  fCentrality;           /// Centrality
  Int_t  fRunNumber;           /// Run Number
  Float_t  fTriggerCheck;         /// Stores trigger information
  Bool_t  fUseCentralityV0M;         /// Stores trigger information
  Int_t  fEvNumberCounter;         /// EvNumber counter
	Int_t fMCEventType; /// MC eventtype to analyze 1: ccbar 2: bbbar (rest is assigned to 1 or 2  with 50% prob)

  //--------------------- My histograms ------------------
  THnSparse* fHistoEleLambdaMass;         //!<! e-Lambda mass spectra
  THnSparse* fHistoEleLambdaMassRS;         //!<! e-Lambda mass spectra right sign
  THnSparse* fHistoEleLambdaMassWS;         //!<! e-Lambda mass spectra wrong sign
  THnSparse* fHistoEleLambdaMassRSMix;         //!<! e-Lambda mass spectra right sign (mixed event)
  THnSparse* fHistoEleLambdaMassWSMix;         //!<! e-Lambda mass spectra wrong sign (mixed event)
  THnSparse* fHistoEleLambdaMassRSSide;         //!<! e-Lambda mass spectra right sign (mixed event)
  THnSparse* fHistoEleLambdaMassWSSide;         //!<! e-Lambda mass spectra wrong sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtRS;         //!<! e-Lambda mass vs elept spectra right sign
  THnSparse* fHistoEleLambdaMassvsElePtWS;         //!<! e-Lambda mass vs elept spectra wrong sign
  THnSparse* fHistoEleLambdaMassvsElePtRSMix;         //!<! e-Lambda mass vs elept spectra right sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtWSMix;         //!<! e-Lambda mass vs elept spectra wrong sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtRSSide;         //!<! e-Lambda mass vs elept spectra right sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtWSSide;         //!<! e-Lambda mass vs elept spectra wrong sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtRS1;         //!<! e-Lambda mass vs elept spectra right sign
  THnSparse* fHistoEleLambdaMassvsElePtWS1;         //!<! e-Lambda mass vs elept spectra wrong sign
  THnSparse* fHistoEleLambdaMassvsElePtRSMix1;         //!<! e-Lambda mass vs elept spectra right sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtWSMix1;         //!<! e-Lambda mass vs elept spectra wrong sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtRSSide1;         //!<! e-Lambda mass vs elept spectra right sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtWSSide1;         //!<! e-Lambda mass vs elept spectra wrong sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtRS2;         //!<! e-Lambda mass vs elept spectra right sign
  THnSparse* fHistoEleLambdaMassvsElePtWS2;         //!<! e-Lambda mass vs elept spectra wrong sign
  THnSparse* fHistoEleLambdaMassvsElePtRSMix2;         //!<! e-Lambda mass vs elept spectra right sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtWSMix2;         //!<! e-Lambda mass vs elept spectra wrong sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtRSSide2;         //!<! e-Lambda mass vs elept spectra right sign (mixed event)
  THnSparse* fHistoEleLambdaMassvsElePtWSSide2;         //!<! e-Lambda mass vs elept spectra wrong sign (mixed event)
  TH2F* fHistoElePtRS;         //!<! e spectra right sign
  TH2F* fHistoElePtWS;         //!<! e spectra wrong sign
  TH2F* fHistoElePtRSMix;         //!<! e spectra right sign (mixed event)
  TH2F* fHistoElePtWSMix;         //!<! e spectra wrong sign (mixed event)

	//Efficiency calculation
  THnSparse* fHistoEleLambdaMassMCS;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassMCGen;         //!<! EFficiency calculation denominator
  THnSparse* fHistoEleLambdaMassvsElePtMCS;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtMCGen;         //!<! EFficiency calculation denominator
  THnSparse* fHistoEleLambdaMassvsElePtMCS1;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtMCGen1;         //!<! EFficiency calculation denominator
  THnSparse* fHistoEleLambdaMassvsElePtMCS2;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtMCGen2;         //!<! EFficiency calculation denominator
  TH2F* fHistoElePtMCS;         //!<! EFficiency calculation numerator
  TH2F* fHistoElePtMCGen;         //!<! EFficiency calculation denominator

  THnSparse* fHistoElePtvsEtaRS;         //!<! e spectra (right-sign)
  THnSparse* fHistoElePtvsEtaWS;         //!<! e spectra (wrong-sign)
  THnSparse* fHistoElePtvsEtaRSMix;         //!<! e spectra (right-sign, mix)
  THnSparse* fHistoElePtvsEtaWSMix;         //!<! e spectra (wrong-sign, mix)
  THnSparse* fHistoElePtvsEtaMCS;         //!<! e spectra efficiency numerator
  THnSparse* fHistoElePtvsEtaMCGen;         //!<! e spectra efficiency denominator

  THnSparse* fHistoElePtvsLambdaPtRS;         //!<! e-Xi spectra (right-sign)
  THnSparse* fHistoElePtvsLambdaPtWS;         //!<! e-Xi spectra (wrong-sign)
  THnSparse* fHistoElePtvsLambdaPtRSMix;         //!<! e-Xi spectra (right-sign, mix)
  THnSparse* fHistoElePtvsLambdaPtWSMix;         //!<! e-Xi spectra (wrong-sign, mix)
  THnSparse* fHistoElePtvsLambdaPtMCS;         //!<! e-Xi spectra efficiency numerator
  THnSparse* fHistoElePtvsLambdaPtvsLcPtMCS;         //!<! e-Xi spectra efficiency numerator
  THnSparse* fHistoElePtvsLambdaPtMCGen;         //!<! e-Xi spectra efficiency denominator
  THnSparse* fHistoElePtvsLambdaPtvsLcPtMCGen;         //!<! e-Xi spectra efficiency denominator
  THnSparse* fHistoElePtvsLambdaPtMCLcGen;         //!<! e-Xi spectra efficiency denominator

  THnSparse* fHistoElePtvsd0RS;         //!<! e pt-d0 spectra (right-sign)
  THnSparse* fHistoElePtvsd0WS;         //!<! e pt-d0 spectra (wrong-sign)
  THnSparse* fHistoElePtvsd0RSMix;         //!<! e pt-d0 spectra (right-sign, mix)
  THnSparse* fHistoElePtvsd0WSMix;         //!<! e pt-d0 spectra (wrong-sign, mix)
  THnSparse* fHistoElePtvsd0MCS;         //!<! e pt-d0 spectra (right-sign) 
  THnSparse* fHistoElePtvsd0PromptMCS;         //!<! e pt-d0 spectra (right-sign) 
  THnSparse* fHistoElePtvsd0BFeeddownMCS;         //!<! e pt-d0 spectra (right-sign) 

	//Feeddown from Xic0
  THnSparse* fHistoEleLambdaMassFeeddownXic0MCS;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassFeeddownXic0MCGen;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXic0MCS;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXic0MCGen;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXic0MCS1;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXic0MCGen1;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXic0MCS2;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXic0MCGen2;         //!<! EFficiency calculation numerator
  TH2F* fHistoElePtFeeddownXic0MCS;         //!<! EFficiency calculation numerator
  TH2F* fHistoElePtFeeddownXic0MCGen;         //!<! EFficiency calculation numerator
  THnSparse* fHistoElePtvsEtaFeeddownXic0MCS;         //!<! e spectra efficiency numerator
  THnSparse* fHistoElePtvsEtaFeeddownXic0MCGen;         //!<! e spectra efficiency numerator
  THnSparse* fHistoElePtvsLambdaPtFeeddownXic0MCS;         //!<! e-Xi spectra efficiency numerator
  THnSparse* fHistoElePtvsLambdaPtFeeddownXic0MCGen;         //!<! e-Xi spectra efficiency numerator

	//Feeddown from XicPlus
  THnSparse* fHistoEleLambdaMassFeeddownXicPlusMCS;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassFeeddownXicPlusMCGen;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS1;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen1;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS2;         //!<! EFficiency calculation numerator
  THnSparse* fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen2;         //!<! EFficiency calculation numerator
  TH2F* fHistoElePtFeeddownXicPlusMCS;         //!<! EFficiency calculation numerator
  TH2F* fHistoElePtFeeddownXicPlusMCGen;         //!<! EFficiency calculation numerator
  THnSparse* fHistoElePtvsEtaFeeddownXicPlusMCS;         //!<! e spectra efficiency numerator
  THnSparse* fHistoElePtvsEtaFeeddownXicPlusMCGen;         //!<! e spectra efficiency numerator
  THnSparse* fHistoElePtvsLambdaPtFeeddownXicPlusMCS;         //!<! e-Xi spectra efficiency numerator
  THnSparse* fHistoElePtvsLambdaPtFeeddownXicPlusMCGen;         //!<! e-Xi spectra efficiency numerator

	//Checking histograms
  TH1F* fHistoBachPt;      //!<! Bachelor pT histogram
  TH1F* fHistoBachPtMCS;      //!<! Bachelor pT histogram (efficiency numerator)
  TH1F* fHistoBachPtMCGen;      //!<! Bachelor pT histogram (efficiency denominator)
  TH1F* fHistod0Bach;      //!<! Bachelor d0 histogram
  TH2F* fHistoLambdaMassvsPt;     //!<! Lambda mass vs pt histogram
  TH2F* fHistoLambdaMassvsPtMCS;     //!<! Lambda mass vs pt histogram
  TH2F* fHistoLambdaMassvsPtMCGen;     //!<! Lambda mass vs pt histogram
  TH2F* fHistoLambdaPtvsDl;     //!<! Lambda pt vs decay length
  TH2F* fHistoLambdaPtvsDlSide;     //!<! Lambda pt vs decay length
  TH2F* fHistoLambdaPtvsDlMCS;     //!<! Lambda pt vs decay length
  TH2F* fHistoLambdaPtvsDlFeeddownXi0MCS;    //!<! Lambda pt vs decay length
  TH2F* fHistoLambdaPtvsDlFeeddownXiMinusMCS;    //!<! Lambda pt vs decay length
  TH2F* fHistoLambdaPtvsDlFeeddownOmegaMCS;    //!<! Lambda pt vs decay length
  TH2F* fHistoK0sMassvsPt;     //!<! K0s mass vs pt histogram
  TH2F* fHistoElectronTPCPID;     //!<! TPC electron PID
  TH2F* fHistoElectronTOFPID;     //!<! TOF electron PID
  TH2F* fHistoElectronTPCSelPID;     //!<! TPC electron PID after selection
  TH2F* fHistoElectronTOFSelPID;     //!<! TOF electron PID after selection
  TH2F* fHistoElectronTPCPIDSelTOF;     //!<! TPC electron PID after TOF 3 sigma cut
  TH2F* fHistoElectronTPCPIDSelTOFSmallEta;     //!<! TPC electron PID after TOF 3 sigma cut (|eta|<0.6)
  TH2F* fHistoElectronTPCPIDSelTOFLargeEta;     //!<! TPC electron PID after TOF 3 sigma cut (0.8>|eta|>0.6)
  TH2F* fHistoElectronTPCPIDSelTOFEtaDep[8];     //!<! TPC electron PID after TOF 3 sigma cut (0.8>|eta|>0.6)
  THnSparse* fHistoElePtvsCutVarsRS[17];         //!<! e pt- cut variables (Right-sign)
  THnSparse* fHistoElePtvsCutVarsWS[17];         //!<! e pt- cut variables (Wrong-sign)
  THnSparse* fHistoElePtvsCutVarsMCS[17];         //!<! e pt- cut variables (MCS)
  TH2F* fHistoElectronQovPtvsPhi;     //!<! Electron phi distribution
  TH2F* fHistoLambdaQovPtvsPhi;     //!<! Lambda phi distribution
  THnSparse* fHistoLcMCGen;         //!<! Lc in mcArray
  THnSparse* fHistoLcMCGen1;         //!<! Lc in mcArray
  THnSparse* fHistoLcMCGen2;         //!<! Lc in mcArray
  THnSparse* fHistoLcMCS;         //!<! Lc in mcArray
  THnSparse* fHistoLcMCS1;         //!<! Lc in mcArray
  THnSparse* fHistoLcMCS2;         //!<! Lc in mcArray
  THnSparse* fHistoFeedDownXic0MCGen;     //!<! Xic0 in mcArray
  THnSparse* fHistoFeedDownXic0MCGen1;     //!<! Xic0 in mcArray
  THnSparse* fHistoFeedDownXic0MCGen2;     //!<! Xic0 in mcArray
  THnSparse* fHistoFeedDownXic0MCS;     //!<! Xic0 in mcArray
  THnSparse* fHistoFeedDownXic0MCS1;     //!<! Xic0 in mcArray
  THnSparse* fHistoFeedDownXic0MCS2;     //!<! Xic0 in mcArray
  THnSparse* fHistoFeedDownXicPlusMCGen;     //!<! XicPlus in mcArray
  THnSparse* fHistoFeedDownXicPlusMCGen1;     //!<! XicPlus in mcArray
  THnSparse* fHistoFeedDownXicPlusMCGen2;     //!<! XicPlus in mcArray
  THnSparse* fHistoFeedDownXicPlusMCS;     //!<! XicPlus in mcArray
  THnSparse* fHistoFeedDownXicPlusMCS1;     //!<! XicPlus in mcArray
  THnSparse* fHistoFeedDownXicPlusMCS2;     //!<! XicPlus in mcArray
  THnSparse* fHistoFeedDownSigmaMCGen;     //!<! Sigma in mcArray
  THnSparse* fHistoFeedDownSigmaMCGen1;     //!<! Sigma in mcArray
  THnSparse* fHistoFeedDownSigmaMCGen2;     //!<! Sigma in mcArray
  THnSparse* fHistoFeedDownSigmaMCS;     //!<! Sigma in mcArray
  THnSparse* fHistoFeedDownSigmaMCS1;     //!<! Sigma in mcArray
  THnSparse* fHistoFeedDownSigmaMCS2;     //!<! Sigma in mcArray
  THnSparse* fHistoLcElectronMCGen;         //!<! Lc in mcArray
  THnSparse* fHistoLcElectronMCGen1;         //!<! Lc in mcArray
  THnSparse* fHistoLcElectronMCGen2;         //!<! Lc in mcArray
  THnSparse* fHistoLcElectronMCS;         //!<! Lc in mcArray
  THnSparse* fHistoLcElectronMCS1;         //!<! Lc in mcArray
  THnSparse* fHistoLcElectronMCS2;         //!<! Lc in mcArray
  THnSparse* fHistoElectronFeedDownXic0MCGen;     //!<! Xic0 in mcArray
  THnSparse* fHistoElectronFeedDownXic0MCGen1;     //!<! Xic0 in mcArray
  THnSparse* fHistoElectronFeedDownXic0MCGen2;     //!<! Xic0 in mcArray
  THnSparse* fHistoElectronFeedDownXic0MCS;     //!<! Xic0 in mcArray
  THnSparse* fHistoElectronFeedDownXic0MCS1;     //!<! Xic0 in mcArray
  THnSparse* fHistoElectronFeedDownXic0MCS2;     //!<! Xic0 in mcArray
  THnSparse* fHistoElectronFeedDownXicPlusMCGen;     //!<! XicPlus in mcArray
  THnSparse* fHistoElectronFeedDownXicPlusMCGen1;     //!<! XicPlus in mcArray
  THnSparse* fHistoElectronFeedDownXicPlusMCGen2;     //!<! XicPlus in mcArray
  THnSparse* fHistoElectronFeedDownXicPlusMCS;     //!<! XicPlus in mcArray
  THnSparse* fHistoElectronFeedDownXicPlusMCS1;     //!<! XicPlus in mcArray
  THnSparse* fHistoElectronFeedDownXicPlusMCS2;     //!<! XicPlus in mcArray
  THnSparse* fHistoElectronMCGen;         //!<! electron in mcArray (only from charmed baryon)
  THnSparse* fHistoLambdaMCGen;         //!<! Lambda in mcArray (only from charmed baryon)

  THnSparse* fHistoElePtvsV0dlRS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlRS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlRS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlRSSide;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlRSSide1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlRSSide2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlRSMix;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlRSMix1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlRSMix2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWSSide;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWSSide1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWSSide2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWSMix;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWSMix1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlWSMix2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlMCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlMCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlMCS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlFeedDownXic0MCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlFeedDownXic0MCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlFeedDownXic0MCS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlFeedDownXicPlusMCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlFeedDownXicPlusMCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dlFeedDownXicPlusMCS2;         //!<! Feeddown subtraction using Lambda vertex distribution

  THnSparse* fHistoElePtvsV0dcaRS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaRS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaRS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaRSSide;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaRSSide1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaRSSide2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaRSMix;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaRSMix1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaRSMix2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWSSide;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWSSide1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWSSide2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWSMix;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWSMix1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaWSMix2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaMCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaMCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaMCS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaFeedDownXic0MCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaFeedDownXic0MCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaFeedDownXic0MCS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaFeedDownXicPlusMCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaFeedDownXicPlusMCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoElePtvsV0dcaFeedDownXicPlusMCS2;         //!<! Feeddown subtraction using Lambda vertex distribution

  THnSparse* fHistoEleLambdaPtvsV0dlRS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlRS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlRS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlRSSide;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlRSSide1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlRSSide2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlRSMix;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlRSMix1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlRSMix2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWSSide;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWSSide1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWSSide2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWSMix;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWSMix1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlWSMix2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlMCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlMCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlMCS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownXic0MCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownXic0MCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownXic0MCS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS2;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS1;         //!<! Feeddown subtraction using Lambda vertex distribution
  THnSparse* fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS2;         //!<! Feeddown subtraction using Lambda vertex distribution

	TH2D *fHistoResponseElePt; //!<! Response function electron pT <- True ept
	TH2D *fHistoResponseElePt1; //!<! Response function electron pT <- True ept
	TH2D *fHistoResponseElePt2; //!<! Response function electron pT <- True ept
	TH2D *fHistoResponseEleLambdaPt; //!<! Response function e-Lambda pT <- XicPt
	TH2D *fHistoResponseEleLambdaPt1; //!<! Response function e-Lambda pT <- XicPt
	TH2D *fHistoResponseEleLambdaPt2; //!<! Response function e-Lambda pT <- XicPt
	TH2D *fHistoResponseEleLambdaPtFeeddownXic0; //!<! Response function e-Lambda pT <- XicPt
	TH2D *fHistoResponseEleLambdaPtFeeddownXic01; //!<! Response function e-Lambda pT <- XicPt
	TH2D *fHistoResponseEleLambdaPtFeeddownXic02; //!<! Response function e-Lambda pT <- XicPt
	TH2D *fHistoResponseEleLambdaPtFeeddownXicPlus; //!<! Response function e-Lambda pT <- XicPt
	TH2D *fHistoResponseEleLambdaPtFeeddownXicPlus1; //!<! Response function e-Lambda pT <- XicPt
	TH2D *fHistoResponseEleLambdaPtFeeddownXicPlus2; //!<! Response function e-Lambda pT <- XicPt
  THnSparse* fHistoLcPtvseleLambdaPtvsElePtvsLambdaPt;         //!<! pT correlation

  AliNormalizationCounter *fCounter;//!<! Counter for normalization
	TH1F *fHistonEvtvsRunNumber;//!<! nevt vs runnumber
	TH1F *fHistonElevsRunNumber;//!<! nele vs runnumber
	TH1F *fHistonLambdavsRunNumber;//!<! nlambda vs runnumber
	TH1F *fHistoMCEventType;//!<! MC even type

  //Mixing
  Int_t fDoEventMixing; /// flag for event mixing
  Int_t  fNumberOfEventsForMixing; /// maximum number of events to be used in event mixing
	Int_t fNzVtxBins;								/// number of z vrtx bins
	Double_t fZvtxBins[100];						// [fNzVtxBinsDim]
	Int_t fNCentBins;								/// number of centrality bins
	Double_t fCentBins[100];						// [fNCentBinsDim]
  Int_t  fNOfPools; /// number of pools
  TTree** fEventBuffer;   //!<! structure for event mixing
	TObjString *fEventInfo; ///unique event id for mixed event check
  TObjArray* fElectronTracks; /// array of electron-compatible tracks
  TObjArray* fV0Tracks1; /// array of lambda-compatible tracks
  TObjArray* fV0Tracks2; /// array of antilambda-compatible tracks
	std::vector<Double_t> fV0dlArray1; /// array of lambda-compatible tracks' information
  std::vector<Double_t> fV0dlArray2; /// array of antilambda-compatible tracks' information
	std::vector<Double_t> fV0dcaArray1; /// array of lambda-compatible tracks' information
  std::vector<Double_t> fV0dcaArray2; /// array of antilambda-compatible tracks' information

  /// \cond CLASSIMP 
  ClassDef(AliAnalysisTaskSELc2eleLambdafromAODtracks,14); /// class for Lc->e Lambda
  /// \endcond 
};
#endif

