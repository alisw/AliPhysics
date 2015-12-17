#ifndef AliAnalysisTaskSEFT2Simulation_H
#define AliAnalysisTaskSEFT2Simulation_H

//*************************************************************************
// Class AliAnalysisTaskSEFT2Simulation
// AliAnalysisTaskSE for fast track parameterization of detector response
// Author - Framework:	J.Stiller	, jstiller@cern.ch
// Author - FT2:		R. Shahoyan	, ruben.shahoyan@cern.ch
// Started: 10-02-2015
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TH1F.h>
#include "AliAnalysisTaskSE.h"
#include <TObject.h>
#include "AliExternalTrackParam.h"
#include "FT2.h"
#include "AliPIDResponse.h"
#include "AliDecayerEvtGen.h"

class AliESDEvent;
class AliFilteredTreeEventCuts;
class AliFilteredTreeAcceptanceCuts;
class AliITSUReconstructor;
class AliITSURecoDet;
class AliITSURecoLayer;
class AliITSURecoSens;
class TParticle;
class AliVertex;
class AliESDEvent;
class AliESDfriend;
class TTreeSRedirector;

//class FT2;

class AliAnalysisTaskSEFT2Simulation : public AliAnalysisTaskSE
{
public:
	//
	AliAnalysisTaskSEFT2Simulation();
	virtual ~AliAnalysisTaskSEFT2Simulation();
	//
	// Implementation of interface methods
	virtual void UserCreateOutputObjects();
	virtual void Init();
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *option);
	// HIJING particle identification
	TString GetGenerator(Int_t label, AliMCEvent* mcEvent);
	void    GetTrackPrimaryGenerator(Int_t label,AliMCEvent *mcEvent,AliStack *stack,TString &nameGen);
	//
	// functions for esd
	Int_t AddTrack(const AliESDtrack *p){return fesd->AddTrack(p);}
	void AddPrimaryVertex(AliESDVertex	*vtx){return fesd->SetPrimaryVertexTracks(vtx);}
	void Reset(){return fesd->Reset();}
	//
	void SetDebugLevel(Int_t debug) { fDebug = debug; }
	void SetMonitorFT2Ouput(Bool_t monitor) { fUseMonitorTree = monitor ;}
	// getters and setters for selection cuts
	void SetEventCuts(AliFilteredTreeEventCuts* const cuts)								{ fFilteredTreeEventCuts = cuts; }
	void SetAcceptanceCuts(AliFilteredTreeAcceptanceCuts* const cuts)			{ fFilteredTreeAcceptanceCuts = cuts; }
	AliFilteredTreeEventCuts* GetEventCuts() const												{ return fFilteredTreeEventCuts; }
	AliFilteredTreeAcceptanceCuts* GetAcceptanceCuts() const							{ return fFilteredTreeAcceptanceCuts; }
	//
	// Setting flags for FT2 processing
	void SetTuneOnDataOrMC(Bool_t tuneOnDataOrMC) {fTuneOnDataOrMC = tuneOnDataOrMC ;}
	void SetSimMat(Bool_t simMat) { fSimMat = simMat ;}
	void SetUsePid(Bool_t usePID) { fUsePID = usePID ;}
	void SetUseKalman(Bool_t useKalman) { fUseKalman = useKalman ;}
	void SetAllowDecay(Bool_t allowDecay) { fAllowDecay = allowDecay;}
	void SetAllowAbsorbtion(Bool_t allowAbsorption) { fAllowAbsorption = allowAbsorption;}
	void SetUseConversionExtension(Bool_t allowConversions) { fUseConverisons = allowConversions;}
	void SetdNdY(Double_t dNdY) { fdNdY = dNdY;}
	void SetMaxStepTGeo(Double_t maxStepTGeo) {fMaxStepTGeo = maxStepTGeo;}
	void SetTPCParaFile( TString tpcParameterizationFile) {fTpcParameterizationFile = tpcParameterizationFile;}
	void SetXSectionFile( TString crossSectionFile) {fCrossSectionFile = crossSectionFile;}
	void SetRunNumber( Int_t runnumber) {fRunNumber = runnumber;}
	void SetStandaloneTune( Int_t stdTune) {fStandaloneTune = stdTune;}
	void SetStreamLevel(Int_t level) {fStreamLevel = level;}
	static Int_t GetMCTrueTrackMult(AliMCEvent *const mcEvent, AliFilteredTreeEventCuts *const evtCuts, AliFilteredTreeAcceptanceCuts *const accCuts);
	//
	Bool_t IsFromConversion(Int_t label, AliStack *const stack);
	Bool_t IsFromMaterial(Int_t label, AliStack *const stack);
	Bool_t IsFromStrangeness(Int_t label, AliStack *const stack);
	Bool_t IsFastMcFromWeakDecay(Int_t label, AliStack *const stack);
	TParticle *GetMother(TParticle *const particle, AliStack *const stack);
	//
private:
	//
	AliAnalysisTaskSEFT2Simulation(const AliAnalysisTaskSEFT2Simulation &source);
	AliAnalysisTaskSEFT2Simulation& operator=(const AliAnalysisTaskSEFT2Simulation& source);
	//
	Bool_t ApplyTrackCuts(FTProbe *FETtrack); // apply track cuts already on FT2 level
	AliESDtrack *PrepareESDtrack(FTProbe *FETtrack,TParticle *part,Int_t iParticle); // setting up of ESD track like object
	AliESDVertex *SmearPrimaryVertex(Float_t x, Float_t y, Float_t z); // first estimate on primary vertex
	Int_t fDebug; // debug level
	Bool_t fUseMonitorTree; // produce scaled output tree
	AliESDEvent *fesd; // smeared esd event
	TTree		*fESDTree; // esd tree to contain ESD event
	TH1F	*fNentries;	// Monitor histogram
	TH1F	*fCpuEventWatch;	// CPU time for event
	TH1F	*fRealEventWatch;	// Real time for event
	TH1F	*fCpuSingleTrackWatch;	// CPU time for 1000 tracks
	TH1F	*fRealSingleTrackWatch;	// Real time for 1000 tracks
	//
	AliFilteredTreeEventCuts      *fFilteredTreeEventCuts;      // event cuts
	AliFilteredTreeAcceptanceCuts *fFilteredTreeAcceptanceCuts; // acceptance cuts
	//
	TTreeSRedirector* fTreeSRedirector;      //! temp tree to dump output
	TTree* fMCEffTree;											//! FT2 performance monitoring tree
	AliESDtrack* fDummyTrack; //! dummy track for tree init
	FTProbe *fDummyFT2Track; //!dummy ft2 track for tree init
	
	// *** FT2 Initialization and processing flags ***
	FT2 *fFET; // FT2 pointer
	//
	Int_t fStreamLevel;
	Bool_t fTuneOnDataOrMC;	// use tune on data (or mc) option
	Bool_t fSimMat;	// simulate material?
	Bool_t fUsePID;	// use PID hypothesis in tracking?
	Bool_t fUseKalman; // use kalman?
	Bool_t fAllowDecay; // use decay ?
	Bool_t fAllowAbsorption; // use absorption?
	Bool_t fUseConverisons; // use FT2 with extended algrotihm suitable for conversions
	//
	Double_t fdNdY; // estimate on dNdy, typically 2100 for 0-10% centrality
	Double_t fMaxStepTGeo; // typically 1.0
	//
	TString fTpcParameterizationFile; // modular input from parameterization of MC or data
	TString	fCrossSectionFile; // nuclear cross-section for decay and absorption from PDG live
	Int_t	fRunNumber; // run number
	Bool_t fStandaloneTune;

	ClassDef(AliAnalysisTaskSEFT2Simulation,1);
};

#endif
