#ifndef AliRICHClusterFinder_H
#define AliRICHClusterFinder_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////
//  RICH Cluster Finder Class                 //
////////////////////////////////////////////////
#include "AliRICHHitMap.h"
#include "TF1.h"
class AliRICHClusterFinder :
 public TObject
{
public:
    TClonesArray*           fDigits;
    Int_t                   fNdigits;
protected:
    AliRICHSegmentation*    fSegmentation;
    AliRICHResponse*        fResponse;
    TClonesArray*           fRawClusters;
    Int_t                   fChamber;
    Int_t                   fNRawClusters;
    AliRICHHitMapA1*        fHitMap;
    TF1*                    fCogCorr;
    Int_t                   fNperMax;
    Int_t                   fDeclusterFlag;
    Int_t                   fClusterSize;
    Int_t                   fNPeaks; 
 public:
    AliRICHClusterFinder
	(AliRICHSegmentation *segmentation,
	 AliRICHResponse *response, TClonesArray *digits, Int_t chamber);
    AliRICHClusterFinder();
    ~AliRICHClusterFinder(){delete fRawClusters;}
    virtual void SetSegmentation(
	AliRICHSegmentation *segmentation){
	fSegmentation=segmentation;
    }
    virtual void SetResponse(AliRICHResponse *response) {
	fResponse=response;
    }

    virtual void SetDigits(TClonesArray *RICHdigits) {
	fDigits=RICHdigits;
	fNdigits = fDigits->GetEntriesFast();
    }
    
    virtual void SetChamber(Int_t ich){
	fChamber=ich;
    }
    
    virtual void AddRawCluster(const AliRICHRawCluster);
    // Search for raw clusters
    virtual void FindRawClusters();
    virtual void  FindCluster(Int_t i, Int_t j, AliRICHRawCluster &c);
    // Decluster
    virtual void Decluster(AliRICHRawCluster *cluster);
    // Set max. Number of pads per local cluster
    virtual void SetNperMax(Int_t npermax=5) {fNperMax = npermax;}
    // Decluster ?
    virtual void SetDeclusterFlag(Int_t flag=1) {fDeclusterFlag =flag;}
    // Set max. cluster size ; bigger clusters will be rejected
    virtual void SetClusterSize(Int_t clsize=5) {fClusterSize = clsize;}
    // Self Calibration of COG 
    virtual void CalibrateCOG();
    virtual void SinoidalFit(Float_t x, Float_t y, TF1 &func);
    //
    virtual void CorrectCOG(){;}
    
    //
    virtual Bool_t Centered(AliRICHRawCluster *cluster);
    virtual void   SplitByLocalMaxima(AliRICHRawCluster *cluster);
    virtual void   FillCluster(AliRICHRawCluster *cluster, Int_t);
    virtual void   FillCluster(AliRICHRawCluster *cluster) {
	FillCluster(cluster,1);}
    TClonesArray* RawClusters(){return fRawClusters;}
    ClassDef(AliRICHClusterFinder,1) //Class for clustering and reconstruction of space points
};
#endif







