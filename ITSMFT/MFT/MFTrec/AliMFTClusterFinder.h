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
#include "TStopwatch.h"

//====================================================================================================================================================

class AliMFTClusterFinder : public TObject {

public:

  AliMFTClusterFinder();
  ~AliMFTClusterFinder();
  virtual void Clear(const Option_t* /*opt*/);
  void Init(const Char_t *nameGeomFile);
  
  void MakeClusterBranch(TTree *treeCluster);
  void SetClusterTreeAddress(TTree *treeCluster);
  void CreateClusters();

  void ApplyMisalignment(Bool_t applyMisalignment) { fApplyMisalignment = applyMisalignment; }

  void DigitsToClusters(const TObjArray *pDigitList);

  void StartEvent();

private:
 
  static const Int_t fNMaxDigitsPerCluster = AliMFTConstants::fNMaxDigitsPerCluster;
  static const Int_t fNMaxPlanes = AliMFTConstants::fNMaxPlanes;
  static const Int_t fNMaxDetElemPerPlane = AliMFTConstants::fNMaxDetElemPerPlane;
  static const Double_t fCutForAvailableDigits;
  static const Double_t fCutForAttachingDigits;
  static const Double_t fMisalignmentMagnitude;

  TClonesArray *fClustersPerPlane[fNMaxPlanes];    //! [fNPlanes] list of clusters [per plane]

  TClonesArray *fDigitsInCluster;                  //!
  AliMFTDigit *fCurrentDigit;                      //!
  AliMFTCluster *fCurrentCluster;                  //!

  AliMFTSegmentation *fSegmentation;               //!
 
  Int_t fNPlanes;

  Bool_t fApplyMisalignment;                       // For MC, waiting for OCDB...

  TStopwatch *fStopWatch;                          //!

  AliMFTClusterFinder(const AliMFTClusterFinder &source);
  AliMFTClusterFinder& operator=(const AliMFTClusterFinder &source);

  ClassDef(AliMFTClusterFinder,1) 

};

//====================================================================================================================================================

#endif
