//========================================================================  
// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  
// See cxx source for full Copyright notice                                
//========================================================================  
//                       
//                       Class AliEMCALTracker 
//                      -----------------------
// Implementation of the track matching method between barrel tracks and
// EMCAL clusters.
// Besides algorithm implementation, some cuts are required to be set
// in order to define, for each track, an acceptance window where clusters
// are searched to find best match (if any).
// The class accepts as input an ESD container, and works directly on it,
// simply setting, for each of its tracks, the fEMCALindex flag, for each
// track which is matched to a cluster.
// In order to use method, one must launch PropagateBack().
//
// ------------------------------------------------------------------------
// author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//=========================================================================

#ifndef ALIEMCALTRACKER_H
#define ALIEMCALTRACKER_H

#include "AliTracker.h"
#include "TMath.h"
class TList;
class TTree;
class TObjArray;
class AliESDEvent;
class AliESDCaloCluster;
class AliEMCALTrack;
class AliEMCALRecPoint;
class AliEMCALGeometry;

class AliEMCALTracker : public AliTracker 
{
public:

	AliEMCALTracker();
	AliEMCALTracker(const AliEMCALTracker &t);
	AliEMCALTracker& operator=(const AliEMCALTracker &source);
	
	virtual ~AliEMCALTracker() {Clear();}
	
	virtual void        Clear(Option_t *option="ALL");
	virtual Int_t       Clusters2Tracks(AliESDEvent*) {return -1;}
	        void        InitParameters();
	virtual Int_t       LoadClusters(TTree*);
	        Int_t       LoadClusters(AliESDEvent* esd);
	        Int_t       LoadTracks(AliESDEvent* esd);
	virtual Int_t       PropagateBack(AliESDEvent* esd);
	virtual Int_t       RefitInward(AliESDEvent*) {return -1;}
	virtual void        UnloadClusters();
	virtual AliCluster* GetCluster(Int_t) const {return NULL;};
	TTree*              SearchTrueMatches();
	void                SetCorrection(Double_t rho, Double_t x0) {fRho=rho;fX0=x0;}
	void                SetCutAlpha(Double_t min, Double_t max) {fCutAlphaMin=min;fCutAlphaMax=max;}
	void                SetCutAngle(Double_t value) {fCutAngle=value;}
	void                SetCutX(Double_t value) {fCutX=value;}
	void                SetCutY(Double_t value) {fCutY=value;}
	void                SetCutZ(Double_t value) {fCutZ=value;}
	void                SetGeometry(AliEMCALGeometry *geom) {fGeom=geom;}
	void                SetMaxDistance(Double_t value) {fMaxDist=value;}
	void                SetCutNITS(Double_t value) {fCutNITS=value;}
	void                SetCutNTPC(Double_t value) {fCutNTPC=value;}
	void                SetNumberOfSteps(Int_t n) {fNPropSteps=n;if(!n)SetTrackCorrectionMode("NONE");}
	void                SetTrackCorrectionMode(Option_t *option);
	
	enum {
		kUnmatched = -99999
	};
	
	class  AliEMCALMatchCluster : public TObject
	{
	public:
		AliEMCALMatchCluster(Int_t ID, AliEMCALRecPoint *recPoint);
		AliEMCALMatchCluster(Int_t ID, AliESDCaloCluster *caloCluster);
		virtual ~AliEMCALMatchCluster() { }
		//----------------------------------------------------------------------------
		Int_t     Index() const {return fIndex;}
		Int_t     Label() const {return fLabel;}
		Double_t  X() const {return fX;}
		Double_t  Y() const {return fY;} 
		Double_t  Z() const {return fZ;}
		Double_t   Phi() const {return TMath::ATan2(fY, fX);}
	private:
		Int_t     fIndex;  // index of cluster in its native container (ESD or TClonesArray)
		Int_t     fLabel;  // track label of assigned cluster
		Double_t  fX;      // global X position
		Double_t  fY;      // global Y position
		Double_t  fZ;      // global Z position
	};
	
	class  AliEMCALMatch : public TObject
	{
	public:
		AliEMCALMatch();
		AliEMCALMatch(const AliEMCALMatch& copy);
		virtual ~AliEMCALMatch() { }
		//----------------------------------------------------------------------------
		Bool_t&   CanBeSaved() {return fCanBeSaved;}
		Int_t     Compare(const TObject *obj) const;
		Double_t  GetDistance() const {return fDistance;}
		Double_t  GetPt() const {return fPt;}
		Int_t     GetIndexC() const {return fIndexC;}
		Int_t     GetIndexT() const {return fIndexT;}
		Bool_t    IsSortable() const {return kTRUE;}
		void      SetIndexC(Int_t icl) {fIndexC=icl;}
		void      SetIndexT(Int_t itr) {fIndexT=itr;}
		void      SetDistance(Double_t dist) {fDistance=dist;}
		void      SetPt(Double_t pt) {fPt=pt;}
	private:
		Bool_t     fCanBeSaved;  // when true, this match can be saved, otherwise it will not be
		Int_t      fIndexC;      // cluster index in 'fClusters' array
		Int_t      fIndexT;      // track index in 'fTracks' array
		Double_t   fDistance;    // track - cluster distance
		Double_t   fPt;          // track pt
	};

private:

	Double_t   AngleDiff(Double_t angle1, Double_t angle2);
	Double_t   CheckPair(AliEMCALTrack *tr, AliEMCALMatchCluster *cluster);
	Double_t   CheckPairV2(AliEMCALTrack *tr, AliEMCALMatchCluster *cluster);
	Double_t   CheckPairV3(AliEMCALTrack *tr, AliEMCALMatchCluster *cluster);
	Int_t      CreateMatches();
	Bool_t     PropagateToEMCAL(AliEMCALTrack *tr);
	Int_t      SolveCompetitions();
	
	enum ETrackCorr { 
		kTrackCorrNone  = 0, // do not correct for energy loss
		kTrackCorrMMB   = 1, // use MeanMaterialBudget() function to evaluate correction
		kTrackCorrFixed = 2  // use fixed "X0" and "rho" parameters to correct
	};
	
	Int_t       fNPropSteps;      // number of propagation steps (when correcting). If =0 no correction is done
	ETrackCorr  fTrackCorrMode;   // track correction mode
	
	Double_t    fCutX;		      // cut on X difference
	Double_t    fCutY;		      // cut on Y difference
	Double_t    fCutZ;		      // cut on Z difference
	Double_t    fCutAlphaMin;     // cut on difference between track 'alpha' and phi
	Double_t    fCutAlphaMax;     // cut on difference between track 'alpha' and phi
	Double_t    fCutAngle;        // cut on angle between track projection and cluster
	Double_t    fMaxDist;         // maximum allowed total distance between track proj an cluster
	Double_t fCutNITS;         // mimimum number of track hits in the ITS
	Double_t fCutNTPC;         // mimimum number of track hits in the TPC

	Double_t    fRho;             // energy correction: density
	Double_t    fX0;              // energy correction: radiation length

	TObjArray  *fTracks;          //! collection of tracks
	TObjArray  *fClusters;        //! collection of EMCAL clusters (ESDCaloCluster or EMCALRecPoint)
	TList      *fMatches;         //! collection of matches between tracks and clusters
	
	AliEMCALGeometry *fGeom;      //! EMCAL geometry
	
	ClassDef(AliEMCALTracker, 3)  // EMCAL "tracker"
};

#endif
