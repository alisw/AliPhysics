#ifndef ALIMUONCLUSTERFINDER_H
#define ALIMUONCLUSTERFINDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  MUON Cluster Finder Class                 //
////////////////////////////////////////////////
#include "TObject.h"
class AliMUONHitMap;
class AliMUONHitMapA1;
class AliMUONDigit;
class TClonesArray;
class AliMUONSegmentation;
class AliMUONResponse;
class AliMUONRawCluster;

class TF1;

const Int_t kMaxNeighbours = 24; // max number of neighbours

class AliMUONClusterFinder :
 public TObject
{
 public:
    AliMUONClusterFinder
	(AliMUONSegmentation *segmentation,
	 AliMUONResponse *response, TClonesArray *digits, Int_t chamber);
    AliMUONClusterFinder(const AliMUONClusterFinder& clusterFinder);
    AliMUONClusterFinder();
    virtual  ~AliMUONClusterFinder(); 
// Set segmentation model    
    virtual void SetSegmentation(
	AliMUONSegmentation *segmentation){
	fSegmentation=segmentation;
    }
// Set response model    
    virtual void SetResponse(AliMUONResponse *response) {
	fResponse=response;
    }
// Set pointer to digits
    virtual void SetDigits(TClonesArray *MUONdigits);
    virtual void SetDigits(TClonesArray *MUONdigits1,
			   TClonesArray *MUONdigits2 ) {;}

// Set current chamber id    
    virtual void SetChamber(Int_t ich){
	fChamber=ich;
    }
// Add a new raw cluster    
    virtual void AddRawCluster(const AliMUONRawCluster cluster);
// Search for raw clusters
    virtual void FindRawClusters();
// Find cluster    
    virtual void  FindCluster(Int_t i, Int_t j, AliMUONRawCluster &c);
// Decluster
    virtual void Decluster(AliMUONRawCluster *cluster);
// Set max. number of pads per local cluster
    virtual void SetNperMax(Int_t npermax=5) {fNperMax = npermax;}
// Decluster ?
    virtual void SetDeclusterFlag(Int_t flag=1) {fDeclusterFlag =flag;}
// Set max. cluster size ; bigger clusters will deconvoluted
    virtual void SetClusterSize(Int_t clsize=5) {fClusterSize = clsize;}
// Self Calibration of COG 
    virtual void CalibrateCOG();
// Perform fit to sinoidal function    
    virtual void SinoidalFit(Float_t x, Float_t y, TF1 &func);
//
    virtual void CorrectCOG(){;}
//  True if 3-cluster is centred
    virtual Bool_t Centered(AliMUONRawCluster *cluster);
//  Perform split by local maxima  
    virtual void   SplitByLocalMaxima(AliMUONRawCluster *cluster);
//  Perform Double Mathieson Fit
    Bool_t  DoubleMathiesonFit(AliMUONRawCluster *c);
    Bool_t SingleMathiesonFit(AliMUONRawCluster *c);    
//  Build up full cluster information    
    virtual void   FillCluster(AliMUONRawCluster *cluster, Int_t n);
    virtual void   FillCluster(AliMUONRawCluster *cluster) {
	FillCluster(cluster,1);}
    virtual void SetTracks(Int_t, Int_t) {;}
    virtual Bool_t TestTrack(Int_t) {return kTRUE;}    
    
//  Return pointer to raw clusters    
    TClonesArray* RawClusters(){return fRawClusters;}
//  Assignment operator
    AliMUONClusterFinder & operator = (const AliMUONClusterFinder& rhs);
    
protected:
    TClonesArray*           fDigits;         // Digits
    Int_t                   fNdigits;        // Number of Digits
    AliMUONSegmentation*    fSegmentation;   // Chamber segmentation
    AliMUONResponse*        fResponse;       // Chamber Response
    TClonesArray*           fRawClusters;    // Raw Clusters
    Int_t                   fChamber;        // Chamber Number
    Int_t                   fNRawClusters;   // Number of Raw Clusters
    AliMUONHitMapA1*        fHitMap;         // Hit Map
    TF1*                    fCogCorr;        // Systematic correction function 
    Int_t                   fNperMax;        // Maximum number of pads per
                                             // local maximum 
    Int_t                   fDeclusterFlag;  // flaf for declusterin
    Int_t                   fClusterSize;    // cluster size 
    Int_t                   fNPeaks;         // number of local maxima
//  Current cluster    
    AliMUONDigit*           fDig[100];        // current list of digits 
    Int_t                   fIx[100];         // current list of x-pad-coord.
    Int_t                   fIy[100];         // current list of y-pad-coord.
    Float_t                 fX[100];          // current list of x-coord.
    Float_t                 fY[100];          // current list of y-coord.
    Int_t                   fIndLocal[100];   // indices of local maxima
    Int_t                   fNLocal;          // Number of local maxima
    Int_t                   fQ[100];          // current list of charges
    Int_t                   fMul;             // current multiplicity
    ClassDef(AliMUONClusterFinder,1) //Class for clustering and reconstruction of space points
};
#endif







