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
//                                                                        // 
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDTRACKERV1_H
#include "AliTRDtrackerV1.h"
#endif

class TTree;
class TTreeSRedirector;
class AliTRDtrackV1;
class AliTRDseedV1;

class AliTRDtrackerDebug : public AliTRDtrackerV1
{
public:
	AliTRDtrackerDebug();
	~AliTRDtrackerDebug();

	void				Draw(const Option_t *);

	Bool_t      Init();
	Bool_t      Open(const char *method);
	Int_t       Process();

	void        ResidualsClustersTrack(const AliTRDseedV1 *tracklet);
	void        ResidualsClustersTracklet(const AliTRDseedV1 *tracklet) const;
	void        ResidualsClustersParametrisation(const AliTRDseedV1 *tracklet) const;
	void        ResidualsTrackletsTrack() const;


private:
	AliTRDtrackerDebug(const AliTRDtrackerDebug &);
	AliTRDtrackerDebug& operator=(const AliTRDtrackerDebug &);


  TTreeSRedirector *fOutputStreamer;                 //!Output streamer
	TTree            *fTree;       // debug tree
	AliTRDseedV1     *fTracklet;   // current tracklet
	AliTRDtrackV1    *fTrack;      // current TRD track
	Int_t            fNClusters;   // N clusters for current track
	Float_t          fAlpha;       // sector

	ClassDef(AliTRDtrackerDebug, 1) // debug suite of the TRD tracker
};

#endif

