#ifndef ALIRICHCLUSTERFINDER_H
#define ALIRICHCLUSTERFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////
//  RICH Cluster Finder Class                 //
////////////////////////////////////////////////
class AliRICHHitMapA1;

#include "TF1.h"
#include "TObject.h"
class  TClonesArray;
class AliSegmentation;
class AliRICHRawCluster;
class AliRICHResponse;
class TClonesArray;


class AliRICHClusterFinder : public TObject
{
 public:
    AliRICHClusterFinder
	(AliSegmentation *segmentation,
	 AliRICHResponse *response, TClonesArray *digits, Int_t chamber);
    AliRICHClusterFinder();
    AliRICHClusterFinder(const AliRICHClusterFinder & ClusterFinder);
    virtual ~AliRICHClusterFinder();
    virtual void SetSegmentation(
	AliSegmentation *segmentation){
	fSegmentation=segmentation;
    }
    virtual void SetResponse(AliRICHResponse *response) {
	fResponse=response;
    }

    virtual void SetDigits(TClonesArray *RICHdigits);
    
    virtual void SetChamber(Int_t ich){
	fChamber=ich;
    }
    
    virtual void AddRawCluster(const AliRICHRawCluster c);
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
    virtual void   FillCluster(AliRICHRawCluster *cluster, Int_t flag);
    virtual void   FillCluster(AliRICHRawCluster *cluster) {
	FillCluster(cluster,1);}
   
    TClonesArray* RawClusters(){return fRawClusters;}
    AliRICHClusterFinder& operator=(const AliRICHClusterFinder& rhs);
    ClassDef(AliRICHClusterFinder,1) //Class for clustering and reconstruction of space points

protected:
    AliSegmentation*        fSegmentation;                 //Segmentation model
    AliRICHResponse*        fResponse;                     //Response model
    TClonesArray*           fRawClusters;                  //Raw clusters list
    AliRICHHitMapA1*        fHitMap;                       //Hit Map with digit positions
    TF1*                    fCogCorr;                      //Correction for center of gravity
    TClonesArray*           fDigits;                       //List of digits
    Int_t                   fNdigits;                      //Number of digits
    Int_t                   fChamber;                      //Chamber number
    Int_t                   fNRawClusters;                 //Number of raw clusters
    Int_t                   fNperMax;                      //Number of pad hits per local maximum
    Int_t                   fDeclusterFlag;                //Split clusters flag
    Int_t                   fClusterSize;                  //Size of cluster 
    Int_t                   fNPeaks;                       //Number of maxima in the cluster
};
#endif







