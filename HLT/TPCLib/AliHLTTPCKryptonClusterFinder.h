// $Id$

#ifndef AliHLTTPC_KRYPTONCLUSTERFINDER
#define AliHLTTPC_KRYPTONCLUSTERFINDER
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCKryptonClusterFinder.h
    @author Kenneth Aamodt kenneth.aamodt@student.uib.no
    @date   
    @brief  Krypton Cluster Finder for the TPC
*/

//#include "AliHLTLogging.h"
//#include "AliHLTTPCPad.h"
#include "AliHLTTPCClusterFinder.h"
#include "TString.h"
#include "TH1F.h"
#include "TObjArray.h"

class AliHLTTPCSpacePointData;
class AliHLTTPCDigitReader;

/**
 * @class AliHLTTPCKryptonClusterFinder
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCKryptonClusterFinder : public AliHLTTPCClusterFinder {

 public:
  /** standard constructor */
  AliHLTTPCKryptonClusterFinder();
  /** destructor */
  virtual ~AliHLTTPCKryptonClusterFinder();

  /** Rebunches the data, use on real data which has "wrong" bunches due to keeping 0 data */
  void ReBunch(const UInt_t * bunchData,Int_t bunchSize);

  /** rads the data insorted */
  void ReadDataUnsorted(void* ptr,unsigned long size);

  /** compare one pads combining neighbouring clustercandidates to a cluster */
  Bool_t ComparePads(AliHLTTPCPad *nextPad,AliHLTTPCClusters* cluster,Int_t nextPadToRead);

  /** Find clusters on the rows */
  void FindRowClusters();

  /** combines the row clusters to a krypton cluster */
  void FindKryptonClusters();

  /** checks if there is a candidate on the previous row */
  void CheckForCandidateOnPreviousRow(AliHLTTPCClusters* tmpCluster);

  /** set the selection from minrow to maxrow, used to look at a certain interval of rows */
  void SetSelection(Int_t minRow, Int_t maxRow);

  /** returns the number of krypton clusters found */
  AliHLTUInt32_t GetNKryptonClusters(){ return fNKryptonClusters;}

  /** sets the maximum size of the output buffer */
  void SetMaxOutputSize(AliHLTUInt32_t size){fMaxOutputSize=size;}
  
  vector<AliHLTUInt16_t> fHWAddressVector;                         //! transient

 private: 
  /** copy constructor prohibited */
  AliHLTTPCKryptonClusterFinder(const AliHLTTPCKryptonClusterFinder&);
  /** assignment operator prohibited */
  AliHLTTPCKryptonClusterFinder& operator=(const AliHLTTPCKryptonClusterFinder&);

  vector<Int_t> fTimebinsInBunch;                                  //! transient

  vector<Int_t> fIndexOfBunchStart;                                //! transient
  
  Int_t fMaxQOfCluster;                                            //! transient
  
  Int_t fSelectionMinRowNumber;                                    //! transient
  Int_t fSelectionMaxRowNumber;                                    //! transient

  AliHLTUInt32_t fNKryptonClusters;                                //! transient

  AliHLTUInt32_t fMaxOutputSize;                                   //! transient
  
  ClassDef(AliHLTTPCKryptonClusterFinder,2) //TPC Krypton cluster finder
};
#endif
