#ifndef AliMUONClusterFinderV0_H
#define AliMUONClusterFinderV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  MUON Cluster Finder Class                 //
////////////////////////////////////////////////
#include "AliMUONClusterFinder.h"
//#include "TF1.h"

class AliMUONClusterFinderV0 :
public AliMUONClusterFinder {
 public:
    AliMUONClusterFinderV0
	(AliMUONSegmentation *segmentation,
	 AliMUONResponse *response, TClonesArray *digits, Int_t chamber);
    AliMUONClusterFinderV0();
    ~AliMUONClusterFinderV0(){delete fRawClusters;}
    virtual void FindRawClusters();
    // Specific methods
    virtual void SetOffset(AliMUONRawCluster *cluster);
    virtual Int_t PeakOffsetAndCoordinates(Int_t DigitIndex, Float_t *X, Float_t *Y);
    // Decluster
    virtual void Decluster(AliMUONRawCluster *cluster);
    //
    virtual Bool_t Centered(AliMUONRawCluster *cluster);
    virtual void   SplitByLocalMaxima(AliMUONRawCluster *cluster);
    void AliMUONClusterFinderV0::DumpCluster(class AliMUONRawCluster *);
    
    
    TClonesArray* RawClusters(){return fRawClusters;}
    ClassDef(AliMUONClusterFinderV0,1) //Class for clustering and reconstruction of space points
};
#endif







