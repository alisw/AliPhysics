/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//====================================================================================================================================================
//
//      Class for finding and building the clusters of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#ifndef AliMFTClusterFinder_H
#define AliMFTClusterFinder_H

#include "AliLog.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "AliMFTDigit.h"
#include "AliMFTCluster.h"
#include "AliMFTSegmentation.h"
#include "TTree.h"
#include "TMath.h"
#include "AliMFTConstants.h"

//====================================================================================================================================================

class AliMFTClusterFinder : public TObject {

public:

  AliMFTClusterFinder();
  ~AliMFTClusterFinder();

  void Init(const Char_t *nameGeomFile);
  
  void MakeClusterBranch(TTree *treeCluster);
  void SetClusterTreeAddress(TTree *treeCluster);
  void CreateClusters();

  void DigitsToClusters(const TObjArray *pDigitList);

  void StartEvent();

private:
 
  static const Int_t fNMaxDigitsPerCluster = AliMFTConstants::fNMaxDigitsPerCluster;
  static const Int_t fNMaxPlanes = AliMFTConstants::fNMaxPlanes;
  static const Double_t fCutForAvailableDigits;
  static const Double_t fCutForAttachingDigits;

  TClonesArray *fClustersPerPlane[fNMaxPlanes];    // ![fNPlanes] list of clusters [per plane]

  TClonesArray *fDigitsInCluster;
  AliMFTDigit *fCurrentDigit;
  AliMFTCluster *fCurrentCluster;

  AliMFTSegmentation *fSegmentation;

  Int_t fNPlanes;

  AliMFTClusterFinder(const AliMFTClusterFinder &source);
  AliMFTClusterFinder& operator=(const AliMFTClusterFinder &source);

  ClassDef(AliMFTClusterFinder,1) 

};

//====================================================================================================================================================

#endif
