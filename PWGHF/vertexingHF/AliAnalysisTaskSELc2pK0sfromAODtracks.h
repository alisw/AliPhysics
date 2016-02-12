#ifndef ALIANALYSISTASKSELC2PK0SFROMAODTRACKS_H
#define ALIANALYSISTASKSELC2PK0SFROMAODTRACKS_H

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
#include "TVector.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliPID.h"
#include "AliRDHFCutsLctopK0sfromAODtracks.h"

/// \class AliAnalysisTaskSELc2pK0sfromAODtracks

class THnSparse;
class TH1F;
class TH2F;
class TH3F;
class TClonesArray;
class AliAODRecoCascadeHF;
class AliESDVertex;
class AliAODMCParticle;
class AliEventPoolManager;
class AliNormalizationCounter;

class AliAnalysisTaskSELc2pK0sfromAODtracks : public AliAnalysisTaskSE 
{
 public:
  AliAnalysisTaskSELc2pK0sfromAODtracks();
  AliAnalysisTaskSELc2pK0sfromAODtracks(const Char_t* name, AliRDHFCutsLctopK0sfromAODtracks* cuts, Bool_t writeVariableTree=kTRUE);
  virtual ~AliAnalysisTaskSELc2pK0sfromAODtracks();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void FillROOTObjects(AliAODRecoCascadeHF *lcobj, AliAODv0 *v0, AliAODTrack *trk, AliAODTrack *trkpid, AliAODEvent *aod, TClonesArray *mcarray);
  void FillMixROOTObjects(TLorentzVector *pt, TLorentzVector *ev, TVector *tinfo, TVector *vinfo);
  Bool_t MakeMCAnalysis(TClonesArray *mcArray);
  void MakeAnalysis(AliAODEvent *aod, TClonesArray *mcArray);
  void FillProtonROOTObjects(AliAODTrack *trk, TClonesArray *mcArray);
  void FillMCROOTObjects(AliAODMCParticle *part, AliAODMCParticle *mcepart, AliAODMCParticle *mcv0part, Int_t decaytype);
  void FillV0ROOTObjects(AliAODv0 *v0, TClonesArray *mcArray);
	void FillMCProtonROOTObjects(AliAODMCParticle *mcepart, TClonesArray *mcArray);
	void FillMCV0ROOTObjects(AliAODMCParticle *mcv0part, TClonesArray *mcArray);
  void SelectV0( const AliVEvent *event,Int_t nV0,Int_t &nSeleV0, Bool_t *seleV0Flags, TClonesArray *mcArray);
  void SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Int_t *seleFlags, TClonesArray *mcArray);
	Int_t MatchToMC(AliAODRecoCascadeHF *elobj, TClonesArray *mcArray, Int_t *pdgarray_pr, Int_t *pdgarray_v0, Int_t *labelarray_pr, Int_t *labelarray_v0,  Int_t &ngen_pr, Int_t &ngen_v0);


  /// set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}
  void SetUseCentralityV0M(Bool_t centon) {fUseCentralityV0M = centon;}
  Bool_t GetUseCentralityV0M() const {return fUseCentralityV0M;}

  void SetReconstructPrimVert(Bool_t a) { fReconstructPrimVert=a; }

  AliAODRecoCascadeHF* MakeCascadeHF(AliAODv0 *casc, AliAODTrack *trk, AliAODTrack *trkpid, AliAODEvent *aod, AliAODVertex *vert);
  AliAODVertex* ReconstructSecondaryVertex(AliAODv0 *casc, AliAODTrack *trk, AliAODEvent *aod);

  void StoreGlobalTrackReference(AliAODTrack *track, Int_t);
  void ResetGlobalTrackReference();

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
	void SetPoolRPBinLimits(Int_t Nrpbins,const Double_t *rpbins){
		fNRPBins = Nrpbins;
		for(int ix = 0;ix<fNRPBins+1;ix++){fRPBins[ix] = rpbins[ix];}
	}
  void DoEventMixingWithPools(Int_t index);
  void ResetPool(Int_t poolIndex);
  Int_t GetPoolIndex(Double_t zvert, Double_t mult, Double_t rp);

 private:

  AliAnalysisTaskSELc2pK0sfromAODtracks(const AliAnalysisTaskSELc2pK0sfromAODtracks &source);
  AliAnalysisTaskSELc2pK0sfromAODtracks& operator=(const AliAnalysisTaskSELc2pK0sfromAODtracks& source); 

  void DefineTreeVariables();
  void DefineProtonTreeVariables();
  void DefineV0TreeVariables();
  void DefineMCTreeVariables();
  void DefineMCProtonTreeVariables();
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
  TH1F *fHReactionPlane;                //!<! Histogram to check Reaction plane
  AliRDHFCutsLctopK0sfromAODtracks *fAnalCuts; /// Cuts - sent to output slot 2
  Bool_t fIsEventSelected;          /// flag for event selected
  Bool_t    fWriteVariableTree;     /// flag to decide whether to write the candidate variables on a tree variables
  Bool_t    fWriteEachVariableTree;     /// flag to decide whether to write the candidate variables on a tree variables
  Bool_t    fWriteMCVariableTree;     /// flag to decide whether to write the candidate variables on a tree variables
  TTree    *fVariablesTree;         //!<! tree of the candidate variables after track selection on output slot 4
  TTree    *fProtonVariablesTree;         //!<! tree of the candidate variables after track selection on output slot 5
  TTree    *fV0VariablesTree;         //!<! tree of the candidate variables after track selection on output slot 6
  TTree    *fMCVariablesTree;         //!<! tree of the candidate variables after track selection on output slot 7
  TTree    *fMCProtonVariablesTree;         //!<! tree of the candidate variables after track selection on output slot 9
  TTree    *fMCV0VariablesTree;         //!<! tree of the candidate variables after track selection on output slot 10
  Bool_t fReconstructPrimVert;       /// Reconstruct primary vertex excluding candidate tracks
  Bool_t fIsMB;       /// MB trigger event
  Bool_t fIsSemi;     /// SemiCentral trigger event
  Bool_t fIsCent;     /// Central trigger event
  Bool_t fIsINT7;     /// INT7 trigger event
  Bool_t fIsEMC7;     /// EMC7 trigger event
  Float_t *fCandidateVariables;   //!<! variables to be written to the tree
  Float_t *fCandidateProtonVariables;   //!<! variables to be written to the tree
  Float_t *fCandidateV0Variables;   //!<! variables to be written to the tree
  Float_t *fCandidateMCVariables;   //!<! variables to be written to the tree
  Float_t *fCandidateMCProtonVariables;   //!<! variables to be written to the tree
  Float_t *fCandidateMCV0Variables;   //!<! variables to be written to the tree
  AliAODVertex *fVtx1;            /// primary vertex
  AliESDVertex *fV1;              /// primary vertex
  Float_t  fVtxZ;         /// zVertex
  Double_t fBzkG;                 /// magnetic field value [kG]
  Float_t  fCentrality;           /// Centrality
  Float_t  fReactionPlane;           /// ReactionPlane
  Int_t  fRunNumber;           /// Run Number
  Float_t  fTriggerCheck;         /// Stores trigger information
  Bool_t  fUseCentralityV0M;         /// Stores trigger information
  Int_t  fEvNumberCounter;         /// EvNumber counter
	AliNormalizationCounter *fCounter;//!<! Counter for normalization

  //--------------------- My histograms ------------------
	TH1F *fHistonEvtvsRunNumber;     //!<! QA histogram
	TH1F *fHistonProtonvsRunNumber;     //!<! QA histogram
	TH1F *fHistonK0svsRunNumber;     //!<! QA histogram

  THnSparse* fHistoLcMCGen;         //!<! Lc Generation
  THnSparse* fHistoLcK0SpMass;         //!<! Lc mass spectra
  THnSparse* fHistoLcK0SpMassMix;         //!<! Lc mass spectra
  THnSparse* fHistoLcK0SpMassCoarse;         //!<! Lc mass spectra
  THnSparse* fHistoLcK0SpMassMixCoarse;         //!<! Lc mass spectra
  THnSparse* fHistoK0spCorrelation;         //!<! Lc mass spectra
  THnSparse* fHistoK0spCorrelationMix;         //!<! Lc mass spectra
  THnSparse* fHistoK0spCorrelationMCS;         //!<! Lc mass spectra
  THnSparse* fHistoLcK0SpMassMCS;         //!<! Lc mass spectra
  THnSparse* fHistoLcK0SpPi0MassMCS;         //!<! Lc mass spectra
  THnSparse* fHistoLcKPluspMass;         //!<! Lc mass spectra
  THnSparse* fHistoLcKMinuspMass;         //!<! Lc mass spectra
  THnSparse* fHistoLcKPluspMassMix;         //!<! Lc mass spectra
  THnSparse* fHistoLcKMinuspMassMix;         //!<! Lc mass spectra

  TH2F* fHistoBachPt;      //!<! Bachelor pT histogram
  TH2F* fHistoBachPtMCS;      //!<! Bachelor pT histogram (efficiency numerator)
  TH2F* fHistoBachPtMCGen;      //!<! Bachelor pT histogram (efficiency denominator)
  TH2F* fHistoKaonPt;      //!<! Kaon pT histogram
  TH2F* fHistoKaonPtMCS;      //!<! Kaon pT histogram (efficiency numerator)
  TH2F* fHistoKaonPtMCGen;      //!<! Kaon pT histogram (efficiency denominator)
  TH3F* fHistoK0sMassvsPt;     //!<! K0s mass vs pt histogram
  TH3F* fHistoK0sMassvsPtMCS;     //!<! K0s mass vs pt histogram
  TH3F* fHistoK0sMassvsPtMCGen;     //!<! K0s mass vs pt histogram
  TH1F* fHistod0Bach;      //!<! Bachelor d0 histogram
  TH1F* fHistod0V0;        //!<! V0 d0 histogram
  TH1F* fHistod0d0;        //!<! Bachelor d0 * V0 d0 histogram
  TH1F* fHistoV0CosPA;     //!<! V0 cosine pointing angle to primary vertex
  TH1F* fHistoProbProton;  //!<! Probability to be proton histogram
  TH1F* fHistoDecayLength; //!<! Decay length histogram
  TH1F* fHistoK0SMass;     //!<! K0s mass histogram
  TH1F* fHistoMassTagV0Min; //!<! electron-any mass
  TH1F* fHistoMassTagV0SameSignMin; //!<! electron-any mass

	TH2D *fHistoResponseLcPt; //!<! Response function Lc pT
	TH2D *fHistoResponseLcPt1; //!<! Response function Lc pT
	TH2D *fHistoResponseLcPt2; //!<! Response function Lc pT

  // Store pointers to global tracks for pid and dca
  AliAODTrack **fGTI;                //! Array of pointers, just nicely sorted according to the id
  Int_t *fGTIndex;                //! Array of integers to keep the index of tpc only track
  const UShort_t  fTrackBuffSize;          //! Size of the above array, ~12000 for PbPb

  //Mixing
  Int_t fDoEventMixing; /// flag for event mixing
  Int_t  fNumberOfEventsForMixing; /// maximum number of events to be used in event mixing
	Int_t fNzVtxBins;								/// number of z vrtx bins
	Double_t fZvtxBins[100];						// [fNzVtxBinsDim]
	Int_t fNCentBins;								/// number of centrality bins
	Double_t fCentBins[100];						// [fNCentBinsDim]
	Int_t fNRPBins;								/// number of rp bins
	Double_t fRPBins[100];						// [fNRPBinsDim]
  Int_t  fNOfPools; /// number of pools
  TTree** fEventBuffer;   //!<! structure for event mixing
	TObjString *fEventInfo; ///unique event id for mixed event check
  TObjArray* fProtonTracks; /// array of electron-compatible tracks
  TObjArray* fV0Tracks; /// array of lambda-compatible tracks
	TObjArray* fProtonCutVarsArray; /// array of RDHF cut information
	TObjArray* fV0CutVarsArray; /// array of RDHF cut information


  /// \cond CLASSIMP     
  ClassDef(AliAnalysisTaskSELc2pK0sfromAODtracks,5); /// class for Lc->p K0
  /// \endcond
};
#endif

