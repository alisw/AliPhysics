#ifndef AliMUONClusterFinderv0_H
#define AliMUONClusterFinderv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  MUON Cluster Finder Class                 //
////////////////////////////////////////////////
#include "AliMUONClusterFinder.h"
//#include "TF1.h"

class AliMUONClusterFinderv0 :
public AliMUONClusterFinder {
//public:
//    TClonesArray*           fDigits;
//    Int_t                   fNdigits;
//protected:
//    AliMUONsegmentation*    fSegmentation;
//    AliMUONresponse*        fResponse;
//    TClonesArray*           fRawClusters;
//    Int_t                   fChamber;
//    Int_t                   fNRawClusters;
//    AliMUONHitMapA1*        fHitMap;
//    TF1*                   fCogCorr;
    
 public:
    AliMUONClusterFinderv0
	(AliMUONsegmentation *segmentation,
	 AliMUONresponse *response, TClonesArray *digits, Int_t chamber);
    AliMUONClusterFinderv0();
    ~AliMUONClusterFinderv0(){delete fRawClusters;}
//    virtual void SetSegmentation(
//	AliMUONsegmentation *segmentation){
//	fSegmentation=segmentation;
//   }
//    virtual void SetResponse(AliMUONresponse *response) {
//	fResponse=response;
//   }

//    virtual void SetDigits(TClonesArray *MUONdigits) {
//	fDigits=MUONdigits;
//	fNdigits = fDigits->GetEntriesFast();
//   }
    
//    virtual void SetChamber(Int_t ich){
//	fChamber=ich;
//   }
    
//    virtual void AddRawCluster(const AliMUONRawCluster);
    // Search for raw clusters
    virtual void FindRawClusters();
//    virtual void  FindCluster(Int_t i, Int_t j, AliMUONRawCluster &c);
    // Specific methods
    virtual void SetOffset(AliMUONRawCluster *cluster);
    virtual Int_t PeakOffsetAndCoordinates(Int_t DigitIndex, Float_t *X, Float_t *Y);
    // Decluster
    virtual void Decluster(AliMUONRawCluster *cluster);
    // Self Calibration of COG 
//    virtual void CalibrateCOG();
//    virtual void SinoidalFit(Float_t x, Float_t y, TF1 &func);
    //
//    virtual void CorrectCOG(){;}
    
    //
    virtual Bool_t Centered(AliMUONRawCluster *cluster);
    virtual void   SplitByLocalMaxima(AliMUONRawCluster *cluster);
//    virtual void   FillCluster(AliMUONRawCluster *cluster);
    void DumpCluster(AliMUONRawCluster *cluster);
    
    
    TClonesArray* RawClusters(){return fRawClusters;}
    ClassDef(AliMUONClusterFinderv0,1) //Class for clustering and reconstruction of space points
};
#endif







