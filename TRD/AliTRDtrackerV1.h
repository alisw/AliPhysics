#ifndef ALITRDTRACKERV1_H
#define ALITRDTRACKERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD tracker                                                       //
//                                                                        //
//  Authors:                                                              //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
/////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDTRACKER_H
#include "AliTRDtracker.h"
#endif

#define DEBUG

/**************************************************************************
 * Class Status see source file                                           *
 **************************************************************************/
 
class TFile;
class TTreeSRedirector;
class TClonesArray;

class AliRieman;
class AliESDEvent;

class AliTRDseedV1;
class AliTRDstackLayer;
class AliTRDtrackerFitter;
class AliTRDrecoParam;
class AliTRDtrackV1;
class AliTrackPoint;
class AliTRDtrackerV1 : public AliTRDtracker
{

 public:
	enum{
		kNTimeBins = 35,
		kNPlanes = 6,
		kNSeedPlanes = 4,
		kMaxTracksStack = 100,
		kNConfigs = 15
	};
	AliTRDtrackerV1(AliTRDrecoParam *p = 0x0);
	AliTRDtrackerV1(const TFile *in, AliTRDrecoParam *p);
	virtual ~AliTRDtrackerV1();
  
	Int_t          Clusters2Tracks(AliESDEvent *esd);
	static void    GetExtrapolationConfig(Int_t iconfig, Int_t planes[2]);
	static void    GetSeedingConfig(Int_t iconfig, Int_t planes[4]);
	Int_t          FollowBackProlongation(AliTRDtrackV1 &t);
	Int_t          FollowProlongation(AliTRDtrackV1 &t);
	Int_t          PropagateBack(AliESDEvent *event);
	Int_t          RefitInward(AliESDEvent *event);
	void           SetRecoParam(AliTRDrecoParam *p){fRecoParam = p;}
	void           UnloadClusters();

 protected:
	Double_t       BuildSeedingConfigs(AliTRDstackLayer *layer, Int_t *configs);
	Int_t          Clusters2TracksSM( AliTRDtracker::AliTRDtrackingSector *sector, AliESDEvent *esd);
	Int_t          Clusters2TracksStack(AliTRDstackLayer *layer, TClonesArray *esdTrackList);
	void           CookLabel(AliKalmanTrack *pt, Float_t wrong) const;
	Int_t          GetSeedingLayers(AliTRDstackLayer *layers, Double_t *params);
	void           GetMeanCLStack(AliTRDstackLayer *layers, Int_t *planes, Double_t *params);
	AliTRDseedV1*  GetTracklet(AliTRDtrackV1 *trk, Int_t plane, Int_t &idx);
	virtual Bool_t GetTrackPoint(Int_t index, AliTrackPoint &p) const;
	AliTRDcluster *FindSeedingCluster(AliTRDstackLayer *layers, AliTRDseedV1/*AliRieman*/ *sfit);
	
	Double_t       MakeSeedingPlanes(AliTRDstackLayer *layer);
	AliTRDstackLayer *MakeSeedingLayer(AliTRDstackLayer *layers, Int_t Plane);
	Int_t          MakeSeeds(AliTRDstackLayer *layers, AliTRDseedV1 *sseed, Int_t *ipar);
	AliTRDtrackV1*   MakeTrack(AliTRDseedV1 *seeds, Double_t *params);
  Int_t          PropagateToX(AliTRDtrackV1 &t, Double_t xToGo, Double_t maxStep);
	Int_t          SetTracklet(AliTRDseedV1 *tracklet);

private:
	AliTRDtrackerV1(const AliTRDtrackerV1 &tracker);
	AliTRDtrackerV1 &operator=(const AliTRDtrackerV1 &tracker);
	Double_t       CookLikelihood(AliTRDseedV1 *cseed, Int_t planes[4], Double_t *chi2);
	void           ImproveSeedQuality(AliTRDstackLayer *layer, AliTRDseedV1 *cseed);

private:

	static Double_t      fgTopologicQA[kNConfigs];         //  Topologic quality
	Double_t             fTrackQuality[kMaxTracksStack];  //  Track quality 
	Int_t                fSeedLayer[kMaxTracksStack];     //  Seed layer
	Int_t                fSieveSeeding;                   //! Seeding iterator
	TClonesArray        *fTracklets;                      // List of tracklets for all sectors
	AliTRDrecoParam     *fRecoParam;                      //  Reconstruction parameters
	AliTRDtrackerFitter *fFitter;                         //! Fitter class of the tracker

	ClassDef(AliTRDtrackerV1, 1)                          //  Stand alone tracker development class

};
#endif
