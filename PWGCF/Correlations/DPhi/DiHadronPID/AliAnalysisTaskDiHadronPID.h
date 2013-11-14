#ifndef ALIANALYSYSTASKDIHADRONPID_H
#define ALIANALYSYSTASKDIHADRONPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
* See cxx source for full Copyright notice */ 
/* $Id$ */

#include "AliAnalysisTaskSE.h"
#include "AliEventPoolManager.h"
#include "AliAODTrackCutsDiHadronPID.h"
#include "AliAODEventCutsDiHadronPID.h"
#include "TObjArray.h"
#include "THn.h"

class AliAnalysisTaskDiHadronPID : public AliAnalysisTaskSE {

public:
	// Constructors/Destructors.
	AliAnalysisTaskDiHadronPID();
	AliAnalysisTaskDiHadronPID(const char* name);
	virtual ~AliAnalysisTaskDiHadronPID();

private:
	AliAnalysisTaskDiHadronPID(const AliAnalysisTaskDiHadronPID&);
	AliAnalysisTaskDiHadronPID& operator=(const AliAnalysisTaskDiHadronPID&);

public:
	// Methods from AliAnalysisTaskSE.
	void UserCreateOutputObjects();
	void LocalInit();
	void UserExec(Option_t*);
	void Terminate(Option_t*);

	// Are all cut objects provided to the task?
	Bool_t ReadyToStart() const {return (fEventCuts && fTrackCutsTrigger && fTrackCutsAssociated);}

	// Setters.
    void SetEventCuts(AliAODEventCutsDiHadronPID* eventcuts) {fEventCuts = eventcuts;}
    void SetTrackCutsTrigger(AliAODTrackCutsDiHadronPID* trackcuts) {fTrackCutsTrigger = trackcuts;}
    void SetTrackCutsAssociated(AliAODTrackCutsDiHadronPID* trackcuts) {fTrackCutsAssociated = trackcuts;} 

    void SetNDEtaBins(Int_t nbins) {fNDEtaBins = nbins;}
    void SetNDPhiBins(Int_t nbins) {fNDPhiBins = nbins;}
    void SetMinEventsForMixing(Int_t nevents) {fMinNEventsForMixing = nevents;}
    void SetPoolTrackDepth(Int_t trackdepth) {fPoolTrackDepth = trackdepth;}
    void SetPoolSize(Int_t poolsize) {fPoolSize = poolsize;}
    void SetMixEvents(Bool_t mixevents = kTRUE) {fMixEvents = mixevents;}
    void SetMixTriggers(Bool_t mixtriggers = kTRUE) {fMixTriggers = mixtriggers;}

	void SetCalculateMismatch(Bool_t calcmismatch) {fCalculateMismatch = calcmismatch;}
	void SetMakeTOFCorrelations(Bool_t makeTOF) {fMakeTOFcorrelations = makeTOF;}
	void SetMakeTOFTPCCorrelationsPi(Bool_t makeTOFTPC) {fMakeTOFTPCcorrelationsPi = makeTOFTPC;}
	void SetMakeTOFTPCCorrelationsKa(Bool_t makeTOFTPC) {fMakeTOFTPCcorrelationsKa = makeTOFTPC;}
	void SetMakeTOFTPCCorrelationsPr(Bool_t makeTOFTPC) {fMakeTOFTPCcorrelationsPr = makeTOFTPC;}		
	void SetTOFIntervalFactorTOFTPC(Double_t factor = 1.) {fTOFIntervalFactorTOFTPC = factor;}
	void SetExtendPtAxis(Bool_t extendptaxis) {fExtendPtAxis = extendptaxis;}

	// Overrides methods from AliAnalyisTaskSE.
	void SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kMB);
	void SetDebugLevel(Int_t level);

	// Getters.
	Int_t GetNDEtaBins() const {return fNDEtaBins;}
	Int_t GetNDPhiBins() const {return fNDPhiBins;}
	Int_t GetMinEventsForMixing() const {return fMinNEventsForMixing;}
	Int_t GetPoolTrackDepth() const {return fPoolTrackDepth;}
	Int_t GetPoolSize() const {return fPoolSize;}
	Bool_t GetMixEvents() const {return fMixEvents;}
	Bool_t GetMixTriggers() const {return fMixTriggers;}

	Bool_t GetCalculateMismatch() const {return fCalculateMismatch;}
	Bool_t GetMakeTOFCorrelations() const {return fMakeTOFcorrelations;}
	Bool_t GetMakeTOFTPCCorrelationsPi() const {return fMakeTOFTPCcorrelationsPi;}
	Bool_t GetMakeTOFTPCCorrelationsKa() const {return fMakeTOFTPCcorrelationsKa;}
	Bool_t GetMakeTOFTPCCorrelationsPr() const {return fMakeTOFTPCcorrelationsPr;}		
	Double_t GetTOFIntervalFactorTOFTPC() const {return fTOFIntervalFactorTOFTPC;}
	Bool_t GetExtendPtAxis() const {return fExtendPtAxis;}

private:
	//void FillGlobalTracksArray();
	Bool_t LoadExtMismatchHistos();
	Double_t GenerateRandomHit(Double_t eta);
	void PrintPoolManagerContents();

private:

	// PID Response Object.
	AliPIDResponse*					fPIDResponse;				//! PID Response.

	// Event Cuts Object.
	AliAODEventCutsDiHadronPID*		fEventCuts;					//

	// Track Cuts Object.
	AliAODTrackCutsDiHadronPID*		fTrackCutsTrigger;			//
	AliAODTrackCutsDiHadronPID*		fTrackCutsAssociated;		//

	// Event Pool Manager.
	AliEventPoolManager*			fPoolMgr;					//! Event pool manager.

	// Track Arrays.
	TObjArray* 						fTriggerTracks;				//! 
	TObjArray* 						fAssociatedTracks;			//!
	// TObjArray* fGlobalTracksArray; 

	// Current Event.
	AliAODEvent*					fCurrentAODEvent;			//! Current AOD Event.

	// Output List.
	TList* 							fOutputList;				//! Output List.

	// Histograms.
	TH1F*							fPtSpectrumTOFbins;			//! Pt Spectrum.
	TH3F*							fCorrelationsTOFbins;		//! Correlations Histogram.
	TH3F*							fMixedEventsTOFbins;		//! Mixed Events Histogram.	
	
	TH1F*							fPtSpectrumTOFTPCbins;		//! Pt Spectrum.
	TH3F*							fCorrelationsTOFTPCbins;	//! Correlations Histogram.
	TH3F*							fMixedEventsTOFTPCbins;		//! Mixed Events Histogram.
	TObjArray*						fMixedEventsTOFTPCbinsPID;	//! Mixed Events Histograms, with PID.

	TObjArray*						fTOFhistos;					//! Array holding all correlation functions with TOF information.
	TObjArray*						fTOFmismatch;				//! Array holding mismatches, using fTOFPtAxis.
	TAxis*							fTOFPtAxis;					//! P_t axis used for the TOF correlation histograms.
	TObjArray*						fTOFTPChistos;				//! Array holding all correlation functions with TOF and TPC information.
	TObjArray*						fTOFTPCmismatch;			//! Array holding mismatches, using fTOFTPCPtAxis.
	TAxis*							fTOFTPCPtAxis;				//! P_t axis used for the TOF/ TPC correlation histograms.

	// Settings.
	Int_t							fNDEtaBins;					//
	Int_t							fNDPhiBins;					//
	Int_t							fMinNEventsForMixing;		// Pool needs at least this many events for mixing.
	Int_t							fPoolTrackDepth;			// For the pool.
	Int_t							fPoolSize;					// 
	Bool_t							fMixEvents;					// NOT YET IMPLEMENTED.
	Bool_t							fMixTriggers;				// If true, triggers are mixed, if not true, associateds are mixed.
	Bool_t							fCalculateMismatch;			//

	// TOF mismatch stuff.
	TH1F*							fT0Fill;					//
	TH2F*							fLvsEta;					//
	TObjArray*						fLvsEtaProjections;			//		

	// Flags.
	Bool_t							fMakeTOFcorrelations;		//
	Bool_t							fMakeTOFTPCcorrelationsPi;	//
	Bool_t							fMakeTOFTPCcorrelationsKa;	//
	Bool_t							fMakeTOFTPCcorrelationsPr;	//
	Double_t						fTOFIntervalFactorTOFTPC;	// Makes the TOF interval longer while keeping the resolution constant.
	Bool_t							fExtendPtAxis;				// Extends pT 

	ClassDef(AliAnalysisTaskDiHadronPID, 3);

};

#endif
