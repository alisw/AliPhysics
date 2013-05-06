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
	void SetDebugLevel(Int_t debuglevel) {fDebug = debuglevel;}

	void SetMakeTOFCorrelations(Bool_t makeTOF) {fMakeTOFcorrelations = makeTOF;}
	void SetMakeTOFTPCCorrelations(Bool_t makeTOFTPC) {fMakeTOFTPCcorrelations = makeTOFTPC;}

	// Getters.
	Int_t GetNDEtaBins() const {return fNDEtaBins;}
	Int_t GetNDPhiBins() const {return fNDPhiBins;}
	Int_t GetMinEventsForMixing() const {return fMinNEventsForMixing;}
	Int_t GetPoolTrackDepth() const {return fPoolTrackDepth;}
	Int_t GetPoolSize() const {return fPoolSize;}
	Bool_t GetMixEvents() const {return fMixEvents;}
	Bool_t GetMixTriggers() const {return fMixTriggers;}
	Int_t GetDebugLevel() const {return fDebug;}	

	Bool_t GetMakeTOFCorrelations() const {return fMakeTOFcorrelations;}
	Bool_t GetMakeTOFTPCCorrelations() const {return fMakeTOFTPCcorrelations;}

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
	TH1F*							fPtSpectrum;				//! Pt Spectrum.
	TH3F*							fCorrelations;				//! Correlations Histogram.
	TH3F*							fMixedEvents;				//! Mixed Events Histogram.	
	
	TObjArray*						fTOFhistos;					//! Array holding all correlation functions with TOF information.
	TAxis*							fTOFPtAxis;					//! P_t axis used for the TOF correlation histograms.
	TObjArray*						fTOFTPChistos;				//! Array holding all correlation functions with TOF and TPC information.
	TAxis*							fTOFTPCPtAxis;				//! P_t axis used for the TOF/ TPC correlation histograms.

	// Settings.
	Int_t							fNDEtaBins;					//
	Int_t							fNDPhiBins;					//
	Int_t							fMinNEventsForMixing;		// Pool needs at least this many events for mixing.
	Int_t							fPoolTrackDepth;			// For the pool.
	Int_t							fPoolSize;					// 
	Bool_t							fMixEvents;					// NOT YET IMPLEMENTED.
	Bool_t							fMixTriggers;				// If true, triggers are mixed, if not true, associateds are mixed.
	Bool_t							fCalculateTOFmismatch;		//

	// TOF mismatch stuff.
	TH1F*							fT0Fill;					//
	TH2F*							fLvsEta;					//
	TObjArray*						fLvsEtaProjections;			//		

	// Flags.
	Bool_t							fMakeTOFcorrelations;		//
	Bool_t							fMakeTOFTPCcorrelations;	//
	Int_t							fDebug;						// Debug flag.

	ClassDef(AliAnalysisTaskDiHadronPID,3);

};

#endif
