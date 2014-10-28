/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//
// AliCDMesonTracks
// for
// AliAnalysisTaskCDMeson
//
//  Author:
//  Felix Reidt <Felix.Reidt@cern.ch>
//
// class applies the track cuts and provides access to the tracks
//
//

// ROOT classes
#include "TObjArray.h"

// AliRoot classes
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrackCuts.h"
#include "AliVTrack.h"
#include "AliMultiplicity.h"

// header of the class itsself
#include "AliCDMesonTracks.h"
#include "assert.h"


//------------------------------------------------------------------------------
AliCDMesonTracks::AliCDMesonTracks()
	: fAODEvent(0x0)
	, fESDEvent(0x0)
	, fDoAOD(kFALSE)
	, fDoSoft(kTRUE)
	, fIsValid(kFALSE)
	, fNTrk0(-999)
	, fNch(-999)
	, fNsoft(-999)
	, fNcombined(-999)
	, fNITSpureSA(-999)
	, fNtrackletsCentralBarrel(-999)
	, fNtrackletsForward(-999)
	, fTracks(0x0)
	, fSoftTracks(0x0)
{
	//
	// constructor
	//
}


//------------------------------------------------------------------------------
AliCDMesonTracks::AliCDMesonTracks(const AliCDMesonTracks& tracks)
	: fAODEvent(tracks.fAODEvent)
	, fESDEvent(tracks.fESDEvent)
	, fDoAOD(tracks.fDoAOD)
	, fDoSoft(tracks.fDoSoft)
	, fIsValid(tracks.fIsValid)
	, fNTrk0(tracks.fNTrk0)
	, fNch(tracks.fNch)
	, fNsoft(tracks.fNsoft)
	, fNcombined(tracks.fNcombined)
	, fNITSpureSA(tracks.fNITSpureSA)
	, fNtrackletsCentralBarrel(tracks.fNtrackletsCentralBarrel)
	, fNtrackletsForward(tracks.fNtrackletsForward)
	, fTracks(tracks.fTracks)
	, fSoftTracks(tracks.fSoftTracks)
{
	//
	// copy constructor
	//
}


//------------------------------------------------------------------------------
AliCDMesonTracks& AliCDMesonTracks::operator=(const AliCDMesonTracks& tracks)
{
	//
	// assignment operator
	//

			if (this != &tracks) {
			// create deep copy ...
		}
		return *this;
}


//------------------------------------------------------------------------------
AliCDMesonTracks::~AliCDMesonTracks()
{
	//
	// destructor
	//

	if (fTracks) {
		fTracks->SetOwner(kTRUE);
		fTracks->Clear();
		delete fTracks;
		fTracks = 0x0;
	}
	if (fSoftTracks) {
		fSoftTracks->SetOwner(kTRUE);
		fSoftTracks->Clear();
		delete fSoftTracks;
		fSoftTracks = 0x0;
	}
}


//------------------------------------------------------------------------------
Bool_t AliCDMesonTracks::ProcessEvent(AliAODEvent* aodEvent,
                                      AliESDEvent* esdEvent,
                                      Bool_t doSoft /* = kTRUE */)
{
	//
	// this function controlls the processing of an event, after the processing is
	// done, the results can be used via the getters
	//

	if ((aodEvent && esdEvent) || (!aodEvent && !esdEvent)) {
		// only one event type is allowed, no event type is also not allowed
		return fIsValid = fDoAOD = kFALSE;
	}
	else if (aodEvent) {
		fAODEvent = aodEvent;
		fESDEvent = 0x0;
		fDoAOD = kTRUE;
	}
	else if (esdEvent) {
		fAODEvent = 0x0;
		fESDEvent = esdEvent;
		fDoAOD = kFALSE;
	}
	fDoSoft = doSoft;

	ApplyCuts();
	if (!fDoAOD) GetRemainingTracklets();
	fIsValid = kTRUE;
	return kTRUE;
}


//------------------------------------------------------------------------------
AliVTrack* AliCDMesonTracks::GetTrack(UInt_t index) const
{
	//
	// provides access to the selected tracks, normal tracks have lower indices
	// than soft tracks
	//

	if (!fIsValid) return 0x0; // cut selection wasn't properly done

	if ((((Int_t)index >= fNch) && !fDoSoft) || // soft not enabled
	    ((Int_t)index >= fNcombined)) return 0x0; // soft enabled
	// index out of range

	if ((Int_t)index < fNch) return (AliVTrack*)((*fTracks)[index]);
	else if ((Int_t)index < fNcombined) {
		return (AliVTrack*)((*fSoftTracks)[index - fNch]);
	}
	else return 0x0; // something went wrong
}


//------------------------------------------------------------------------------
Double_t AliCDMesonTracks::GetInvariantMass(Bool_t includeSoftTracks /*=kTRUE*/)
{
	//
	// compute the invariant mass of all accepted tracks in the event
	//

	TLorentzVector sum;
	for (Int_t iTrack = 0; iTrack < fTracks->GetEntriesFast(); ++iTrack) {
		AliVTrack* track = (AliVTrack*)fTracks->UncheckedAt(iTrack);
		TLorentzVector temp(track->Px(), track->Py(), track->Pz(), track->E());
		sum += temp;
	}
	if (includeSoftTracks) {
		for (Int_t iSoftTrack = 0; iSoftTrack < fSoftTracks->GetEntriesFast();
		     ++iSoftTrack) {
			AliVTrack* track = (AliVTrack*)fSoftTracks->UncheckedAt(iSoftTrack);
			TLorentzVector temp(track->Px(), track->Py(), track->Pz(), track->E());
			sum += temp;
		}
	}
	return sum.M();
}


//------------------------------------------------------------------------------
void AliCDMesonTracks::ApplyCuts()
{
	//
	// steers the track selection process
	//

	fNTrk0 = (fDoAOD) ?
		fAODEvent->GetNumberOfTracks() : fESDEvent->GetNumberOfTracks();

	fNch = -999;
	fNsoft = -999;
	fNcombined = -999;
	fNITSpureSA = -999;

	if (fDoAOD) {
		CutTrack(fAODEvent); // ordinary tracks
		CutTrack(fAODEvent, 2); // kITSpureSA
	}
	else {
		CutTrack(fESDEvent); // ordinary tracks
		CutTrack(fESDEvent, 2); // kITSpureSA
	}

	if (fDoSoft) { // do soft tracks
		if (fDoAOD) {
			CutTrack(fAODEvent, 1);
		}
		else {
			CutTrack(fESDEvent, 1);
		}

		fNcombined = fNch + fNsoft;
		for (Int_t iSoft = 0; iSoft < fNsoft; iSoft++) {
			Int_t iTrk = 0;
			while (iTrk < fNch) {
				// TODO find some criterion to match them properly!
				// check whether they are really complementary if not
				// an error will be raised
				// already contained in both arrays, exit loop
				++iTrk; // next track
			}
		}
	}
	else { // do not care about soft tracks
		fNcombined = fNch;
	}
}


//------------------------------------------------------------------------------
void AliCDMesonTracks::CutTrack(AliESDEvent *ESDEvent, Int_t mode /* = 0 */)
{
	//
	//CutTrack to be combined with the AOD function // TODO
	//

	const Double_t etacut = 0.9;

	AliESDtrackCuts esdTrackCuts;

	if (mode == 0) { // default mode
		// cuts for normal tracks (ITS + TPC)
		// i.e. GetStandardITSTPCTrackCuts2010(kTRUE);
		// (same, just typed in full detail...)

		if (fTracks) {
			fTracks->Clear();
			delete fTracks;
			fTracks = 0x0;
		}

		// TPC
		esdTrackCuts.SetMinNClustersTPC(70);
		//esdTrackCuts.SetMinNCrossedRowsTPC(70);
		//esdTrackCuts.SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);

		esdTrackCuts.SetMaxChi2PerClusterTPC(4);
		esdTrackCuts.SetAcceptKinkDaughters(kFALSE);
		esdTrackCuts.SetRequireTPCRefit(kTRUE);
		// ITS
		esdTrackCuts.SetRequireITSRefit(kTRUE);
		esdTrackCuts.SetClusterRequirementITS(AliESDtrackCuts::kSPD,
		                                      AliESDtrackCuts::kAny);
		// 7*(0.0026+0.0050/pt^1.01)
		esdTrackCuts.SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

		esdTrackCuts.SetMaxDCAToVertexZ(2);
		esdTrackCuts.SetDCAToVertex2D(kFALSE);
		esdTrackCuts.SetRequireSigmaToVertex(kFALSE);

		esdTrackCuts.SetEtaRange(-etacut, etacut);

		fTracks = esdTrackCuts.GetAcceptedTracks(ESDEvent);
		fNch = fTracks->GetEntriesFast();
	}
	else if (mode == 1) {
		// cuts for soft tracks (ITS only - kITSsa)

		if (fSoftTracks) {
			fSoftTracks->Clear();
			delete fSoftTracks;
			fSoftTracks = 0x0;
		}

		esdTrackCuts.SetRequireITSStandAlone(kTRUE);
		esdTrackCuts.SetRequireITSPureStandAlone(kFALSE);
		esdTrackCuts.SetRequireITSRefit(kTRUE);
		esdTrackCuts.SetMinNClustersITS(4);
		esdTrackCuts.SetClusterRequirementITS(AliESDtrackCuts::kSPD,
		                                      AliESDtrackCuts::kAny);
		esdTrackCuts.SetMaxChi2PerClusterITS(1.);
		esdTrackCuts.SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55");

		esdTrackCuts.SetEtaRange(-etacut, etacut);

		fSoftTracks = esdTrackCuts.GetAcceptedTracks(ESDEvent);
		fNsoft = fSoftTracks->GetEntriesFast();
	}
	else {
		// cuts for ITSpureSA tracks used in order to get rid of them for noise
		// studies

		// selection according to cuts
		esdTrackCuts.SetRequireITSPureStandAlone(kTRUE);

		// do selection according to status bits (never tested!!!)
		//for(Int_t itrack = 0; itrack < ESDEvent->GetNumberOfTracks(); itrack++){
		//	const AliESDtrack* esdtrack = ESDEvent->GetTrack(itrack);
		//	UInt64 status = esdtrack->GetStatus();
		//	if ((status & kITSpureSA) && !(status & kITSsa)){
		//	}
		//}
		TObjArray* arr = esdTrackCuts.GetAcceptedTracks(ESDEvent);
		fNITSpureSA = arr->GetEntriesFast();
		delete arr;
		arr = 0x0;
	}
}


//------------------------------------------------------------------------------
void AliCDMesonTracks::CutTrack(AliAODEvent *AODEvent, Int_t mode /* = 0 */)
{
	//
	// CutTrack for AODs
	//

	UInt_t bit = 0x0;
	TObjArray* trks = 0x0;
	Int_t* ntrks = 0x0;

	if (mode == 0) { // default mode
		// cuts for normal tracks (ITS + TPC)
		// i.e. GetStandardITSTPCTrackCuts2010(kTRUE);
		// (same, just typed in full detail...)
		bit = 0x1 << 14;

		// prepare storage
		if (fTracks) {
			fTracks->SetOwner(kTRUE);
			fTracks->Clear();
		}
		else {
			fTracks = new TObjArray();
			fTracks->SetOwner(kTRUE);
		}

		// store where to put selected tracks
		trks = fTracks;
		ntrks = &fNch;
	}
	else if (mode == 1) {
		// cuts for soft tracks (ITS only - kITSsa)
		bit = 0x1 << 15;

		if (fSoftTracks) {
			fSoftTracks->SetOwner(kTRUE);
			fSoftTracks->Clear();
		}
		else {
			fSoftTracks = new TObjArray();
			fSoftTracks->SetOwner(kTRUE);
		}

		// sotre where to put selected tracks
		trks = fSoftTracks;
		ntrks = &fNsoft;
	}
	else {
		// cuts for ITSpureSA tracks used in order to get rid of them for noise
		// studies
		bit = 0x1 << 16;

		// do not store tracks, just count them =>
		trks = 0x0;
		ntrks = &fNITSpureSA;
	}

	for (Int_t iTrk = 0; iTrk < AODEvent->GetNumberOfTracks(); iTrk++) {
		const AliAODTrack* trk = dynamic_cast<const AliAODTrack*>(AODEvent->GetTrack(iTrk));
                assert(trk&&"Not a standard AOD");

		if (trk->TestFilterBit(bit)) {
			// test whether track was selected by that filter

			if (trks) { // add tracks to TObjArray
				trks->Add((TObject*)trk);
			}
			else { // just count them
				++(*ntrks);
			}
		}
	}

	if (trks) {
		(*ntrks) = trks->GetEntriesFast();
	}
}


//------------------------------------------------------------------------------
void AliCDMesonTracks::GetRemainingTracklets()
{
	// determines the number of tracklets in an event, which are not assigned to
	// tracks
	// this is only possible running on ESDs, for AODs this information has to be
	// preprocessed

	if (!fESDEvent) return;
	const AliMultiplicity *mult = fESDEvent->GetMultiplicity();

	if (mult) {
		// reset values
		fNtrackletsCentralBarrel = 0;
		fNtrackletsForward = 0;

		for (Int_t iTracklet = 0; iTracklet < mult->GetNumberOfTracklets();
		     iTracklet++) {
			Int_t id1 = -1, id2 = -1;
			if (!mult->GetTrackletTrackIDs(iTracklet, 0, id1, id2)) {
				float_t eta = mult->GetEta(iTracklet);

				if ((eta < -0.9) || (eta > 0.9)) {
					++fNtrackletsForward;
				}
				else {
					++fNtrackletsCentralBarrel;
				}
			}
		}
	}
}
