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

class TList;
class TObjArray;
class AliESD;
class AliESDCaloCluster;
class AliEMCALTrack;
class AliEMCALRecPoint;

class AliEMCALTracker : public AliTracker 
{
public:

	AliEMCALTracker();
	AliEMCALTracker(const AliEMCALTracker &t);
	AliEMCALTracker& operator=(const AliEMCALTracker &source);
	
	virtual ~AliEMCALTracker() {Clear();}
	
	virtual void        Clear(Option_t *option="ALL");
	virtual Int_t       Clusters2Tracks(AliESD*) {return -1;}
	virtual Int_t       LoadClusters(TTree*);
	        Int_t       LoadClusters(AliESD* esd);
	        Int_t       LoadTracks(AliESD* esd);
	virtual Int_t       PropagateBack(AliESD* esd);
	virtual Int_t       RefitInward(AliESD*) {return -1;}
	virtual void        UnloadClusters();
	virtual AliCluster* GetCluster(Int_t) const {return NULL;};
	void                SetCorrection(Double_t rho, Double_t x0) {fRho=rho;fX0=x0;}
	void                SetCutAlpha(Double_t min, Double_t max) {fCutAlphaMin=min;fCutAlphaMax=max;}
	void                SetCutAngle(Double_t value) {fCutAngle=value;}
	void                SetCutX(Double_t value) {fCutX=value;}
	void                SetCutY(Double_t value) {fCutY=value;}
	void                SetCutZ(Double_t value) {fCutZ=value;}
	void                SetMaxDistance(Double_t value) {fMaxDist=value;}
	void                SetNumberOfSteps(Int_t n) {fNPropSteps=n;if(!n)SetTrackCorrectionMode("NONE");}
	void                SetTrackCorrectionMode(Option_t *option);
	
	class  AliEMCALMatchCluster : public TObject
	{
	public:
		           AliEMCALMatchCluster(Int_t ID, AliEMCALRecPoint *recPoint);
				   AliEMCALMatchCluster(Int_t ID, AliESDCaloCluster *caloCluster);
		virtual   ~AliEMCALMatchCluster() { }
		Int_t&     Index() {return fIndex;}
		Int_t&     Label() {return fLabel;}
		Double_t&  X() {return fX;}
		Double_t&  Y() {return fY;}
		Double_t&  Z() {return fZ;}
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
		         ~AliEMCALMatch() { }
		Bool_t&   CanBeSaved() {return fCanBeSaved;}
		Int_t     Compare(const TObject *obj) const;
		Double_t  GetDistance() {return fDistance;}
		Int_t     GetIndexC() {return fIndexC;}
		Int_t     GetIndexT() {return fIndexT;}
		Bool_t    IsSortable() const {return kTRUE;}
		void      SetIndexC(Int_t icl) {fIndexC=icl;}
		void      SetIndexT(Int_t itr) {fIndexT=itr;}
		void      SetDistance(Double_t dist) {fDistance=dist;}
	private:
		Bool_t     fCanBeSaved;  // when true, this match can be saved, otherwise it will not be
		Int_t      fIndexC;      // cluster index in 'fClusters' array
		Int_t      fIndexT;      // track index in 'fTracks' array
		Double_t   fDistance;    // track - cluster distance
	};

private:

	Double_t   AngleDiff(Double_t angle1, Double_t angle2);
	Double_t   CheckPair(AliEMCALTrack *tr, AliEMCALMatchCluster *cluster);
	Double_t   CheckPairV2(AliEMCALTrack *tr, AliEMCALMatchCluster *cluster);
	Int_t      CreateMatches();
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
	
	Double_t    fRho;             // energy correction: density
	Double_t    fX0;              // energy correction: radiation length

	TObjArray  *fTracks;          //! collection of ESD tracks
	TObjArray  *fClusters;        //! collection of EMCAL clusters (ESDCaloCluster or EMCALRecPoint)
	TList      *fMatches;         //! collection of matches between tracks and clusters
	
	ClassDef(AliEMCALTracker, 1)  // EMCAL "tracker"
};

#endif
