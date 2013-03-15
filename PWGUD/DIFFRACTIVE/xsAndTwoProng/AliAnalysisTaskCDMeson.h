/*************************************************************************
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
// Select events according to gap conditions, analyze two track events in pp
// collisions
//
// Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>

#ifndef ALIANALYSISTASKCDMESON_H
#define ALIANALYSISTASKCDMESON_H

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

#include "THnSparse.h" // forward declaration is not possible

class AliESDEvent;
class AliVTrack;
class AliPIDResponse;
class AliPhysicsSelection;

class TH1I;
class TH1F;
class TH2I;
class TH2F;
class TH2D;
class TList;
class TParticle;
class THnSparse;
class TObjString;
class TObjArray;

class AliCDMesonTracks;

class AliAnalysisTaskCDMeson : public AliAnalysisTaskSE
{
public:
	AliAnalysisTaskCDMeson(const char* name, Long_t state = 0x0);
	AliAnalysisTaskCDMeson();
	virtual ~AliAnalysisTaskCDMeson();

	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t *);

private:
	enum { //
		kV0 = 0,
		kFMD,
		kSPD,
		kTPC,
		kV0FMD,
		kV0FMDSPD,
		kV0FMDSPDTPC,
		kTPCSPD,
		kTPCSPDFMD,
		kTPCSPDFMDV0,
		kSPDFMD,
		kSPDFMDV0,
		kMax
	};

	AliAnalysisTaskCDMeson(const AliAnalysisTaskCDMeson  &p);
	AliAnalysisTaskCDMeson& operator=(const AliAnalysisTaskCDMeson  &p);

	void FillEtaPhiMaps() const; // controls which event contributes to which map


	// functions called by the UserExec(...), not to be called elsewhere!
	//-------------------------------------------------------------------
	Bool_t CheckInput();
	void PostOutputs(); // cares about posting the output before exiting UserExec,
	// WARNING: PostOutputs should only be used directly in the UserExec!!
	void DoEmptyEventStudy(); // does the empty event study
	Bool_t DetermineGap(); // determines the gap of all available detectors
	void AnalyzeVtx(); // calcs the distance of the pri. vertex from tracks and SPD
	void DoMultiplicityStudy(Int_t nMCprimaries);
	// compares ITSTPC and ITSTPC+soft distribution and the multiplicity response
	void DoTrackPair(Bool_t soft = kFALSE);
	// analyses track pairs and prepares THnMother(Soft) input
	Int_t DoMCTruth(); // analyses the MCtruth for corrections, returns #Primaries
	void DoMCTrackPair(TParticle* particles[], Int_t gapCond, Int_t multiplicity);
	// analyzes track pairs in MC truth
	void DetermineMCprocessType(); // determines the MC process ID

	// functions called by the DoTrackPair(...), not to be called elsewhere!
	// ---------------------------------------------------------------------
	Int_t DoPID(Int_t combCh); //retrieves the PID information for the current two
	void FillMChists(Int_t combCh); //fills the MC histograms for two track events

	// analysis task status
	//---------------------
	Bool_t fDoAOD; // true for running on AODs
	Long_t fAnalysisStatus; // stores the status bits used to distinguish
	// which histograms should be generated
	TObjString *fAnalysisStatusString; // to be stored in the output list
	Double_t fNoGapFraction; // fraction of no-gap events to be processed only
	                         // used in case of ESD input
	Double_t fReducedGapFraction; // fraction of events stored altough they
	                              // show only the requested and not the refined
	                              // gap
	Double_t fMaxVtxDst; // maximum distance of the track and SPD vertex

	// event information
	//------------------
	AliESDEvent *fESDEvent; // esd event object
	AliAODEvent *fAODEvent; // esd event object
	AliPIDResponse *fPIDResponse; // esd pid object
	AliPhysicsSelection *fPhysicsSelection; // physics selection object
	AliCDMesonTracks *fTracks; // object taking care about the track cuts
	Double_t fVtxDst; // distance of the primary vertex from tracks and from SPD
	Double_t fVtxZ; // z-position of the primary vertex from tracks
	Int_t fResidualTracks; // tracks rejected by cuts within the event
	Int_t fResidualTracklets; // SPD tracklets not assigned to tracks
	Int_t fMCprocessType; // MC process type, 0 for data
	Int_t fMCprocess; // detailed MC sub process information

	// information about the trackpair which is currently processed
	AliVTrack* fTrkPair[2]; // track objects

	Int_t fRun; // number of the run which is about to be processed
	Int_t fPIDmode; // selects set of PID cuts, 0 for 3sigma standard cuts,
	//1 for LHC11f

	// information stored for the PWA (addresses needed for the tree)
	Float_t fTheta; // theta angle in the helicity frame
	Float_t fPhi; // phi angle in the helicity frame
	Float_t fMass; // mass of the two track system
	Float_t fMomentum[3]; // momentum of the two track system
	Int_t fCurrentGapCondition; // gap condition of the current event
	Int_t fGapInformationWCent[kMax]; // gap condition for different detectors
	                                  // individually and their combinations
	                                  // requiring central act. from SPD FastOR
	Int_t fGapInformation[kMax]; // same as above, without the requirement
	                             // on central activity

	// output objects
	//---------------
	THnSparseI *fGapRun; //contains the gap information of the different detectors
	// per run (x-axis: run number, y-axis gap-state (bit encoded)

	TList *fHist; // output list (contains all histograms)
	THnSparseD *fThnMother; // THnSparse for mother pt and mass
	THnSparseD *fThnMotherSoft; // same as above, using soft tracks as well
	THnSparseD *fThnMultiplicity; // ThnSparse for multiplicity studies
	THnSparseD *fThnMC; // THnSparse for MC mother pt and mass

	THnSparseD *fThnEmptyEvents; // THnSparse used to analyse empty events
	TTree *fPWAtree; // Tree containing information needed for PWA

	// Multiplicity distributions for the different gap conditions
	TH2D *fv0ntrk; //v0bit vs. nch
	TH2D *fv0fmdntrk; //v0fmdbit vs. nch
	TH2D *fv0fmdspdntrk; //v0fmdspdbit vs. nch
	TH2D *fv0fmdspdtpcntrk; //v0fmdspdtpcbit vs. nch

	// How many tracks are removed by the current cuts
	TH2F *fv0Rmntrk; // V0 gap events vs. of removed tracks
	TH2F *fv0fmdRmntrk; // V0-FMD gap events vs. of removed tracks
	TH2F *fv0fmdspdRmntrk; // V0-FMD-SPD gap events vs. of removed tracks
	TH2F *fv0fmdspdtpcRmntrk; // V0-FMD-SPD-TPC gap events vs. of removed tracks

	TH2F *fMultStudy; // multiplicity after standard cuts and soft cuts vs. mult
	TH2F *fMultStudyV0dg; // same for with events V0 Double-Gap
	TH2F *fMultStudyV0FMDdg; // same for events with V0-FMD DG
	TH2F *fMultStudyV0FMDSPDdg; // same for events with V0-FMD-SPD DG
	TH2F *fMultStudyV0FMDSPDTPCdg; // same for events with V0-FMD-SPD-TPC DG

	TH2F *fMultResponseMC; // multiplicity response MC to reconstruction
	TH2F *fMultResponseV0dgMC; // same for events with V0 double gap
	TH2F *fMultResponseV0FMDdgMC; // with V0-FMD DG
	TH2F *fMultResponseV0FMDSPDdgMC; // with V0-FMD-SPD DG
	TH2F *fMultResponseV0FMDSPDTPCdgMC; // with V0-FMD-SPD-TPC DG

	TH2F *fMultRegionsMC; // multiplicity in the 3 regions of interest

	TH1I *fhspdV0dg; // FastOR multiplicity in SPD FastOR with V0 double gap
	TH1I *fhspdV0FMDdg; // with V0-FMD double gap
	TH1I *fhspdV0FMDSPDdg; // with V0-FMD-SPD double gap
	TH1I *fhspdV0FMDSPDTPCdg; // with V0-FMD-SPD double gap
	TH1I *fhspdAfterCuts; // fast OR multiplicity after the event (vertex) cuts

	TH1I *fGapResponseMCv0Dg; // gap status from MC and reconstruction for V0
	TH1I *fGapResponseMCv0fmdDg; // for V0-FMD
	TH1I *fGapResponseMCv0fmdspdDg; // for V0-FMD-SPD
	TH1I *fGapResponseMCv0fmdspdtpcDg; // for V0-FMD-SPD-TPC

	// Statistics flow diagrams
	TH1I *fhspd; // distribution of fastOR-multiplicity
	//SPD fastOR offline esdEvent->GetMultiplicity->TestFastOrFiredChips()
	TH2I *fhfo; // histogram containing the SPD FastOR information ??
	TH1I *fhpriVtxDist; // primary vertex distribution
	TH1I *fhpriVtxPos; // position of the primary vertex (0 within range, 1 else)
	TH1I *fhpriv; //primary vertex cut effect
	TH1I *fhntrk; //n-trk-after-cut effect w/o requiring a gap-condition
	TH1F *fhStatsFlow; // stepwise statistics flow
	TH1F *fhFOchans; // hits in the different fastOR channels

	TObjArray* fVZEROhists; // histograms for the VZERO trigger study

	// eta-phi maps containing the ESDtracks
	TH2F *fv0Map; // V0 gap events
	TH2F *fv0fmdMap; // V0-FMD gap events
	TH2F *fv0fmdspdMap; // V0-FMD-SPD gap events
	TH2F *fv0fmdspdtpcMap; // V0-FMD-SPD-TPC gap events

	// eta-phi maps containing the ESDtracks after the track cut
	TH2F *fv0MapCutted; // V0 gap events
	TH2F *fv0fmdMapCutted; // V0-FMD gap events
	TH2F *fv0fmdspdMapCutted; // V0-FMD-SPD gap events
	TH2F *fv0fmdspdtpcMapCutted; // V0-FMD-SPD-TPC gap events

	TH2F *fHitMapSPDinner; // eta-phi map containing SPD FastOR hits (inner layer)
	TH2F *fHitMapSPDouter; // eta-phi map containing SPD FastOR hits (outer layer)
	TH2F *fHitMapSPDtrklt; // eta-phi map containing SPD hits based on tracklets
	TH2F *fHitMapFMDa; // eta-phi-multiplicity map containing FMD hits at A side
	TH2F *fHitMapFMDc; // eta-phi-multiplicity map containing FMD hits at C side

	// Summed FMD
	TH1F *fFMDsum1I; // FMD 1
	TH1F *fFMDsum2I; // FMD 2 inner ring
	TH1F *fFMDsum2O; // FMD 2 outer ring
	TH1F *fFMDsum3I; // FMD 3 inner ring
	TH1F *fFMDsum3O; // FMD 3 outer ring

	// Vertex
	TH1I *fPriVtxX; // x distribution of the primary vertex
	TH1I *fPriVtxY; // y distribution of the primary vertex
	TH1I *fPriVtxZ; // z distribution of the primary vertex
	TH1I *fPriVtxDst; // distance distribution of the pri. vtx from SPD and tracks
	TH1I *fPriVtxDstV0dg; // with V0 DG
	TH1I *fPriVtxDstV0FMDdg; // with V0 FMD DG
	TH1I *fPriVtxDstV0FMDSPDdg; // with V0 FMD SPD DG
	TH1I *fPriVtxDstV0FMDSPDTPCdg; //  with V0 FMD SPD TPC DG

	// TPC Gap
	TH2F *fTPCGapDCAaSide; // dca distribution for tracks destroying the TPC gap
	TH2F *fTPCGapDCAcSide; // for V0-FMD-SPD gap events on A and C side

	// PID
	// All events but double gap events
	TH2F *fComb2trkPIDuls; // PID matrix on unlike-sign two particles events
	TH2F *fComb2trkPIDls; // PID matrix on like-sign about two particles
	// Full-gap events (V0-FMD-SPD-TPC)
	TH2F *fComb2trkPIDulsDG; // PID matrix on unlike-sign two particles events
	TH2F *fComb2trkPIDlsDG; // PID matrix on like-sign about two particles
	// 'PID' in MC
	// All events but double gap events
	TH2F *fComb2trkPIDulsMC; // PID matrix on unlike-sign two particles events
	TH2F *fComb2trkPIDlsMC; // PID matrix on like-sign about two particles
	// Full-gap events (V0-FMD-SPD-TPC)
	TH2F *fComb2trkPIDulsDGmc; // PID matrix on unlike-sign two particles events
	TH2F *fComb2trkPIDlsDGmc; // PID matrix on like-sign about two particles

	// MC Process
	TH2F *fMCProcessUls; // Process fulfilling a gap condition with two
	// unlike-sign particles
	TH2F *fMCProcessLs; // same for like-sign particles

	TH1F *fMAllTrackMass; // histogram for the invariant mass dist. of all tracks

	ClassDef(AliAnalysisTaskCDMeson, 1);
};


#endif
