/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/////////////////////////////////////////////////////////////////////
// Author: A.Mastroserio, C.Terrevoli                              //
//         annalisa.mastroserio@cern.ch                            //
//         cristina.terrevoli@ba.infn.it                           //
// Alternative cluster finder. Usage instructions below.           //
//                                                                 //
//  For each event:                                                //
//  1. Call StartEvent()                                           //
//  2. For each pixel hit:                                         //
//     Call ProcessHit(..)                                         //
//  3. Call FinishEvent()                                          //
//  4. Access cluster information for this event by methods:       //
//     GetClusterCount(layer)                                      //
//     GetClusterMeanCol(layer,index)                              //
//     GetClusterMeanRow(layer,index)                              //
//     GetClusterSize(layer,index)                                 //
//     GetClusterType(layer,index)                                 //
//                                                                 //
/////////////////////////////////////////////////////////////////////

/* $Id$ */

#ifndef ALIITSUPGRADECLUSTERFINDER_H
#define ALIITSUPGRADECLUSTERFINDER_H

#include <Rtypes.h>
#include "AliITSUpgradeClusterList.h"
#include <TObjArray.h>

class TTree;

class AliITSUpgradeClusterFinder :public TObject{

 public:
  AliITSUpgradeClusterFinder();
  ~AliITSUpgradeClusterFinder();

  void ActivateClusterTypeSearch(){fClusterTypeFlag=kTRUE;}
  void DeActivateClusterTypeSearch(){fClusterTypeFlag=kFALSE;}

  void  StartEvent();
  Int_t ProcessHit(Int_t layer, UInt_t col, UInt_t row, UShort_t charge,Int_t label[3]);
  void  FinishEvent();

  void AddLabelIndex(UInt_t col, UInt_t row);
  void SetLabels(Int_t label[3]);
  void MakeRecPointBranch(TTree *treeR);
  void SetRecPointTreeAddress(TTree *treeR);

  void DigitsToRecPoints(const TObjArray *digList);

  UInt_t  GetClusterCount(Int_t layer) const;
  Float_t GetClusterMeanCol(Int_t layer, UInt_t index);
  Float_t GetClusterMeanRow(Int_t layer, UInt_t index);
  UInt_t  GetClusterSize(Int_t layer, UInt_t index);
  UInt_t  GetClusterWidthZ(Int_t layer, UInt_t index) ;
  UInt_t  GetClusterWidthPhi(Int_t layer, UInt_t index) ;
  UInt_t  GetClusterType(Int_t layer, UInt_t index) ;
  UShort_t GetCharge(Int_t layer, UInt_t index); 
  UInt_t GetPixelCharge(UInt_t col, UInt_t row); 
  Int_t* GetLabels(Int_t layer,UInt_t index) ;
  void  PrintClusters(Int_t layer);
  void  PrintAllClusters();

 private:
  
  void   NewEvent();
  void   NewModule();
  Int_t  DoModuleClustering(Int_t Layer, UShort_t charge);
  UInt_t FindClusterRecu(Int_t col, Int_t row, UShort_t charge);
  void   ShiftClusterTypeArea(UInt_t direction);
  UInt_t GetClusterType(UInt_t size);
  UInt_t GetClusterWidthZ();
  UInt_t GetClusterWidthPhi();

  enum {kMAXCLUSTERTYPESIDEZ=3,kMAXCLUSTERTYPESIDEY=4}; // region of interest for cluster type pattern
  enum {kSHIFTRIGHT,kSHIFTDOWN};  // used for shifting the region of interest for cluster type pattern
  
  UInt_t   fNhitsLeft;     // number of hits still left to process for this module
  Bool_t   fHits[39530][39530];// hit map for this module
  UShort_t fHitCol[999999]; // these two arrays remember which pixels are hit for this module
  UShort_t fHitRow[999999]; // these two arrays remember which pixels are hit for this module
  Short_t  fOldModule;     // remember previous module (-1 at start of event)
  Bool_t   fClusterTypeFlag;   // should we classify the clusters at all
  Bool_t   fFindClusterType;   // temporary, for classifying a cluster (pattern of pixels)
  UShort_t fClusterTypeOrigCol;// temporary, for classifying a cluster (pattern of pixels)
  UShort_t fClusterTypeOrigRow;// temporary, for classifying a cluster (pattern of pixels)
  UInt_t   fColSum; // used to determine the center of a cluster
  UInt_t   fRowSum; // used to determine the center of a cluster
  UShort_t fCharge;        // cluster charge 
  Int_t    fTmpLabel[3];   // label array to be kept temporarily during the clustering procedure
  Int_t    fLabels[10];    // label array to be attached to the cluster

  UShort_t fClusterWidthMaxCol; //max column ID of the cluster
  UShort_t fClusterWidthMinCol; //min column ID of the cluster
  UShort_t fClusterWidthMaxRow; //max row ID of the cluster
  UShort_t fClusterWidthMinRow; //min row ID of the cluster
  Bool_t   fClusterTypeArea[kMAXCLUSTERTYPESIDEZ][kMAXCLUSTERTYPESIDEY];// same as above comments
  AliITSUpgradeClusterList **fClusterList; //[fNSectors] cluster container
  TObjArray *fChargeArray;  // charge identifier
  TClonesArray *fRecPoints; // used to fill treeR
  
  Int_t fNSectors;
  AliITSUpgradeClusterFinder(const AliITSUpgradeClusterFinder &source); // copy constructor
  // assignment operator
  AliITSUpgradeClusterFinder& operator=(const AliITSUpgradeClusterFinder &source);


  ClassDef(AliITSUpgradeClusterFinder,1) 

};


#endif
