//-*- Mode: C++ -*-
// $Id: AliHLTTPCHWClusterMergerV1.h 53447 2011-12-06 21:52:47Z richterm $

#ifndef ALIHLTTPCHWCLUSTERMERGERV1_H
#define ALIHLTTPCHWCLUSTERMERGERV1_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCHWClusterMergerV1.h
//  @author Matthias Richter, Sergey Gorbunov
//  @date   2011-11-25
//  @brief  Merger class for HLT TPC Hardware clusters
//          Handles merging of branch border clusters

#include "AliHLTTPCClusterMCData.h"
#include "AliHLTLogging.h"
#include <vector>
#include "TObject.h"

class AliHLTTPCRawCluster;

/**
 * @class AliHLTTPCHWClusterMergerV1
 *
 * @ingroup alihlt_base
 */
class AliHLTTPCHWClusterMergerV1 : public AliHLTLogging
{
 public:
  /// standard constructor
  AliHLTTPCHWClusterMergerV1();
  /// destructor
  ~AliHLTTPCHWClusterMergerV1();

  Int_t Init( Bool_t processingRCU2Data );

  void SetDataPointer(  AliHLTUInt8_t *data ){ fpData=data; }

  Int_t SetDataBlock(  AliHLTComponentBlockData *block);

  /// merge clusters
  int Merge();

  /// cleanup
  void Clear();

 private:
  /// copy constructor
  AliHLTTPCHWClusterMergerV1(const AliHLTTPCHWClusterMergerV1&);
  /// assignment operator
  AliHLTTPCHWClusterMergerV1& operator=(const AliHLTTPCHWClusterMergerV1&);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// helper class to store relevant data for branch border

  struct AliBorderParam {
    AliBorderParam(float a, int b): fPadPosition(a), fPatch(b){}
    Float_t fPadPosition;
    Int_t fPatch;
  };
  
  /// helper class to store relevant data for a cluster at border
  struct AliBorderRecord {
    AliBorderRecord( AliHLTTPCRawCluster *a, AliHLTTPCClusterMCLabel *b, AliHLTUInt32_t c ):fCluster(a), fMC(b), fTimeBin(c){}
    AliHLTTPCRawCluster *fCluster;
    AliHLTTPCClusterMCLabel *fMC;
    AliHLTUInt32_t fTimeBin;    
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  static bool CompareTime( const AliBorderRecord &b1, const AliBorderRecord &b2){
    return b1.fTimeBin > b2.fTimeBin;
  }
 
  static bool CompareMCWeights(const AliHLTTPCClusterMCWeight &a, const AliHLTTPCClusterMCWeight &b){
    return a.fWeight > b.fWeight;
  }
  static bool CompareMCLabels(const AliHLTTPCClusterMCWeight &a, const AliHLTTPCClusterMCWeight &b){
    return a.fMCID < b.fMCID;
  }

  static const int fkMergeWidth = 3;
  static const int fkNSlices = 36;
  static const int fkNPatches = 6;
  static const int fkMergeTimeWindow = 3;

  int fNRows;//!
  int fNRowPads;//!
  int fNBorders;//!

  AliHLTInt16_t *fMapping;//!
  std::vector<AliBorderParam> fBorders; //!
  AliHLTUInt8_t *fpData;
  AliHLTComponentBlockData **fRawClusterBlocks; //!
  AliHLTComponentBlockData **fMCBlocks; //!
};

#endif //ALIHLTTPCHWCLUSTERMERGERV1_H
