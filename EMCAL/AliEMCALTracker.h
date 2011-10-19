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
#include <TMath.h>
#include <TVector3.h>
class TList;
class TTree;
class TObjArray;
class AliESDEvent;
class AliVCluster;
class AliESDCaloCluster;
class AliEMCALTrack;
class AliExternalTrackParam;
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
	virtual AliCluster* GetCluster(Int_t) const {return NULL;}
	void                SetCutEta(Double_t value) {fCutEta=value;}
	void                SetCutPhi(Double_t value) {fCutPhi=value;}
	void                SetGeometry(AliEMCALGeometry *geom) {fGeom=geom;}
	void                SetCutPt(Double_t value)   {fCutPt=value;}
	void                SetCutNITS(Double_t value) {fCutNITS=value;}
	void                SetCutNTPC(Double_t value) {fCutNTPC=value;}
	void                SetStepLength(Float_t length) {fStep=length;}
	void                SetTrackCorrectionMode(Option_t *option);

	static Bool_t ExtrapolateTrackToEMCalSurface(AliExternalTrackParam *trkParam, Double_t emcalR, Double_t mass, Double_t step, Double_t &eta, Double_t &phi);
	static Bool_t ExtrapolateTrackToPosition(AliExternalTrackParam *trkParam, Float_t *clsPos, Double_t mass, Double_t step, Double_t &tmpEta, Double_t &tmpPhi);
	static Bool_t ExtrapolateTrackToCluster(AliExternalTrackParam *trkParam, AliVCluster *cluster, Double_t mass, Double_t step, Double_t &tmpEta, Double_t &tmpPhi);

	enum {	kUnmatched = -99999 };
	
	class  AliEMCALMatchCluster : public TObject
	{
	public:
		AliEMCALMatchCluster(Int_t ID, AliEMCALRecPoint *recPoint);
		AliEMCALMatchCluster(Int_t ID, AliESDCaloCluster *caloCluster);
		virtual ~AliEMCALMatchCluster() { }
		//----------------------------------------------------------------------------
		Int_t     Index() const {return fIndex;}
		Double_t  X() const {return fX;}
		Double_t  Y() const {return fY;} 
		Double_t  Z() const {return fZ;}
	private:
		Int_t     fIndex;  // index of cluster in its native container (ESD or TClonesArray)
		Double_t  fX;      // global X position
		Double_t  fY;      // global Y position
		Double_t  fZ;      // global Z position
	};
   
private:
	Int_t  FindMatchedCluster(AliESDtrack *track);
	
	enum ETrackCorr { 
		kTrackCorrNone  = 0, // do not correct for energy loss
		kTrackCorrMMB   = 1, // use MeanMaterialBudget() function to evaluate correction
	};

	Double_t    fCutPt;           // mimimum pT cut on tracks
	Double_t    fCutNITS;         // mimimum number of track hits in the ITS
	Double_t    fCutNTPC;         // mimimum number of track hits in the TPC
	
	Float_t     fStep;            // Length of each step in propagation
	ETrackCorr  fTrackCorrMode;   // Material budget correction mode
	Double_t    fClusterWindow;   // Select clusters in the window to be matched to tracks
	Double_t    fCutEta;	      // cut on eta difference
	Double_t    fCutPhi;	      // cut on phi difference

	TObjArray  *fTracks;          //! collection of tracks
	TObjArray  *fClusters;        //! collection of EMCAL clusters (ESDCaloCluster or EMCALRecPoint)
	
	AliEMCALGeometry *fGeom;      //! EMCAL geometry
	
	ClassDef(AliEMCALTracker, 5)  // EMCAL "tracker"
};

#endif
