#ifndef ALITRDTRACKERDEBUG_H
#define ALITRDTRACKERDEBUG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackerDebug.h 22646 2007-11-29 18:13:40Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reader for the TRD tracker debug streamer                             // 
//                                                                        // 
//  Authors:                                                              //
//                                                                        //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        // 
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDTRACKERV1_H
#include "AliTRDtrackerV1.h"
#endif

class TTree;
class TCanvas;
class TTreeSRedirector;
class AliTRDtrackV1;
class AliTRDseedV1;

class AliTRDtrackerDebug : public AliTRDtrackerV1
{
public:
	AliTRDtrackerDebug();
	~AliTRDtrackerDebug();

	void				Draw(Option_t *);

	Bool_t      Init();
	Bool_t      Open(const char *method);
	Int_t       Process();

	void        ResidualsClustersTrack(const AliTRDseedV1 *tracklet);
	void        ResidualsClustersTracklet(const AliTRDseedV1 *tracklet) const;
	void        ResidualsClustersParametrisation(const AliTRDseedV1 *tracklet) const;
	void        ResidualsTrackletsTrack() const;
	
	void        AnalyseTiltedRiemanFit();
	void        AnalyseMinMax();
	void        AnalyseFindable(Char_t *treename);

	TCanvas*    PlotSeedingConfiguration(const Char_t *direction, Int_t event, Int_t Candidate);
	TCanvas*    PlotFullTrackFit(Int_t event, Int_t candidate, Int_t iteration = -1, const Char_t *direction = "y");
	
	static Int_t GetEventNumber(){ return fgEventNumber; }
	static Int_t GetTrackNumber(){ return fgTrackNumber; }
	static Int_t GetCandidateNumber(){ return fgCandidateNumber; }
	
	static void SetEventNumber(Int_t eventNumber){ fgEventNumber = eventNumber; }
	static void SetTrackNumber(Int_t trackNumber){ fgTrackNumber = trackNumber; }
	static void SetCandidateNumber(Int_t candidateNumber){ fgCandidateNumber = candidateNumber; }
			
private:
	AliTRDtrackerDebug(const AliTRDtrackerDebug &);
	AliTRDtrackerDebug& operator=(const AliTRDtrackerDebug &);

	Float_t     GetTrackRadius(Float_t a, Float_t b, Float_t c) const;
	Float_t     GetTrackCurvature(Float_t a, Float_t b, Float_t c) const;
	Float_t     GetDCA(Float_t a, Float_t b, Float_t c) const;

	TTreeSRedirector *fOutputStreamer;                 //!Output streamer
	TTree            *fTree;       // debug tree
	AliTRDseedV1     *fTracklet;   // current tracklet
	AliTRDtrackV1    *fTrack;      // current TRD track
	Int_t            fNClusters;   // N clusters for current track
	Float_t          fAlpha;       // sector
	
	static Int_t fgEventNumber;				//  Event Number in the tracking code
	static Int_t fgTrackNumber;       //  Track Number per Event
	static Int_t fgCandidateNumber;   //  Candidate Number per event (Set in MakeSeeds)

	ClassDef(AliTRDtrackerDebug, 1) // debug suite of the TRD tracker
};

#endif

