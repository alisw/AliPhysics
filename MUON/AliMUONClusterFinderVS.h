#ifndef ALIMUONCLUSTERFINDERVS_H
#define ALIMUONCLUSTERFINDERVS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  MUON Cluster Finder Class                 //
////////////////////////////////////////////////
#include "AliMUONHitMap.h"
#include "TF1.h"
#include "AliMUONClusterFinder.h"
#include "AliMUONSegmentation.h"

class AliMUONClusterFinderVS : 
 public AliMUONClusterFinder
{
 public:
    AliMUONClusterFinderVS
	(AliMUONSegmentation *segmentation1, AliMUONSegmentation *segmentation2,
	 AliMUONResponse *response,
	 TClonesArray *digits1, TClonesArray *digits2,
	 Int_t chamber);
    AliMUONClusterFinderVS();
    AliMUONClusterFinderVS(const AliMUONClusterFinderVS& clusterFinder);
    virtual ~AliMUONClusterFinderVS(){;}
// Set segmentation model    
    virtual void SetSegmentation(AliMUONSegmentation *seg1, AliMUONSegmentation *seg2)
	{
	fSegmentation[0]=seg1;
	fSegmentation[1]=seg2;
	}
// Set pointer to digits
    virtual void SetDigits(TClonesArray *MUONdigits1, TClonesArray *MUONdigits2);
    
// Get Segmentation
    virtual AliMUONSegmentation*  Segmentation(Int_t i);
// Get Number of Digits
    virtual Int_t NDigits(Int_t i);
// Get Digits
    virtual TClonesArray* Digits(Int_t i);
// Get HitMap
    virtual AliMUONHitMap* HitMap(Int_t i);
    
// Search for raw clusters
    virtual void FindRawClusters();
// Find cluster    
    virtual void  FindCluster(Int_t i, Int_t j, Int_t cath, AliMUONRawCluster &c);
// Decluster
    virtual void Decluster(AliMUONRawCluster *cluster);
//  Perform split by local maxima  
    virtual void   SplitByLocalMaxima(AliMUONRawCluster *cluster);
    virtual void   FindLocalMaxima(AliMUONRawCluster *cluster);
    virtual void   Split(AliMUONRawCluster * cluster);
    
//  Perform Double Mathieson Fit
    Bool_t DoubleMathiesonFit(AliMUONRawCluster *c, Int_t cath);
    Float_t CombiDoubleMathiesonFit(AliMUONRawCluster *c);
    Float_t SingleMathiesonFit(AliMUONRawCluster *c, Int_t cath);
    Float_t CombiSingleMathiesonFit(AliMUONRawCluster *c);    
//  Build up full cluster information    
    virtual void   FillCluster(AliMUONRawCluster *cluster, Int_t flag, Int_t cath);
    virtual void   FillCluster(AliMUONRawCluster *cluster, Int_t cath);
    virtual void   FillCluster(AliMUONRawCluster *cluster) {
	FillCluster(cluster,1,0);}
    // Add a new raw cluster    
    virtual void AddRawCluster(const AliMUONRawCluster cluster);
    
    virtual void SetTracks(Int_t t1, Int_t t2) 
	{
	    fTrack[0]=t1;
	    fTrack[1]=t2;
	}
    
    virtual Bool_t TestTrack(Int_t t) {
	if (fTrack[0]==-1 || fTrack[1]==-1) {
	    return kTRUE;
	} else if (t==fTrack[0] || t==fTrack[1]) {
	    return kTRUE;
	} else {
	    return kFALSE;
	}
    }
    //  Assignment operator
    AliMUONClusterFinderVS & operator = (const AliMUONClusterFinderVS& rhs);
protected:
    TClonesArray*           fDigits2;            // Digits
    Int_t                   fNdigits2;           // Number of Digits    
    AliMUONHitMapA1*        fHitMap2;            // Hit Map
    AliMUONDigit*           fDig[100][2];        // current list of digits 
    Int_t                   fIx[100][2];         // current list of x-pad-coord.
    Int_t                   fIy[100][2];         // current list of y-pad-coord.
    Float_t                 fX[100][2];          // current list of x-coord.
    Float_t                 fY[100][2];          // current list of y-coord.
    Int_t                   fIndLocal[100][2];   // indices of local maxima
    Int_t                   fNLocal[2];          // Number of local maxima
    Int_t                   fQ[100][2];          // current list of charges
    Int_t                   fMul[2];             // current multiplicity
// Current Fit
    Double_t                 fXFit[2];         // x-coordinate
    Double_t                 fYFit[2];         // y-coordinate
    Double_t                 fQrFit[2];        // charge ratio
    Float_t                  fChi2[2];         // chi2 of fit
    Float_t                  fXInit[2];        // start values
    Float_t                  fYInit[2];        // start values
    Float_t                  fQrInit[2];       // start values
    Int_t                    fFitStat;         // status of fit
    
// Selected track for debugging
    Int_t                    fTrack[2];        // Only digits with main contributions from these tracks are
                                               // considered 
//  Return pointer to raw clusters    
    ClassDef(AliMUONClusterFinderVS,1) //Class for clustering and reconstruction of space points
};
#endif










