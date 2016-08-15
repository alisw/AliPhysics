
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
// AliPBTracks
// for 
// AliAnalysisTaskPB
//
// see AliPBTracks.cxx for description
//
//  Author:
//  Felix Reidt <Felix.Reidt@cern.ch>

#ifndef ALIPBTRACKS_H
#define ALIPBTRACKS_H

class TObjArray;

class AliAODEvent;
class AliESDEvent;
class AliVTrack;

class AliPBTracks
{
public: 
	AliPBTracks(); // constructor
	~AliPBTracks(); // destructor

	Bool_t ProcessEvent(AliAODEvent* aodEvent, AliESDEvent* esdEvent,
	                    Bool_t doSoft = kTRUE);

	Int_t GetTracksBeforeCuts() const { return (fIsValid) ? fNTrk0 : -1; }
	Int_t GetTracks() const { return (fIsValid) ? fNch : -1; }
	Int_t GetSoftTracks() const { return (fIsValid && fDoSoft) ? fNsoft : -1; }
	Int_t GetCombinedTracks() const { return (fIsValid) ? fNcombined : -1; }
	Int_t GetITSpureSACount() const { return (fIsValid) ? fNITSpureSA : -1; }
	Int_t GetRemainingTrackletsCentralBarrel() const {
		return (fIsValid) ? fNtrackletsCentralBarrel : -1;
	}
	Int_t GetRemainingTrackletsForward() const {
		return (fIsValid) ? fNtrackletsForward : -1;
	}
	AliVTrack* GetTrack(UInt_t index) const;

	Double_t GetInvariantMass(Bool_t includeSoftTracks = kTRUE);
    
protected:
	void ApplyCuts();
	void CutTrack(AliESDEvent *ESDEvent, Int_t mode = 0);
	void CutTrack(AliAODEvent *AODEvent, Int_t mode = 0);
	void GetRemainingTracklets();

	AliAODEvent* fAODEvent; // AOD event to analyze
	AliESDEvent* fESDEvent; // ESD event to analyze
	Bool_t fDoAOD; // is active for AOD processing, inactive for ESDs
	Bool_t fDoSoft; // process soft tracks
	Bool_t fIsValid; // are the stored results valid

	Int_t fNTrk0; // number of tracks before cuts
	Int_t fNch; // number of charged ITSTPC tracks after cuts
	Int_t fNsoft; // number of soft ITS standalone tracks (complementary to fNch)
	Int_t fNcombined; // fNch+fNsoft
	Int_t fNITSpureSA; // number of ITSpureSA tracks leading double counting

	Int_t fNtrackletsCentralBarrel; // tracklets not assigned to tracks within
	                                // |eta| < 0.9
	Int_t fNtrackletsForward; // tracklets not assigned to tracks with |eta| > 0.9

	TObjArray* fTracks; // storage for the standard tracks
	TObjArray* fSoftTracks; // storage for the soft tracks
private:
	// following functions are only implemented to obey the coding conventions,
	// they are not functional and hence should not be used
	AliPBTracks(const AliPBTracks& tracks);
	// copy constructor
	AliPBTracks& operator=(const AliPBTracks& tracks);
};

#endif // ALIPBTRACKS_H
