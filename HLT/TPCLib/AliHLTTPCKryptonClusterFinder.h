// $Id$

#ifndef AliHLTTPC_KRYPTONCLUSTERFINDER
#define AliHLTTPC_KRYPTONCLUSTERFINDER
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCKryptonClusterFinder.h
    @author Kenneth Aamodt kenneth.aamodt@student.uib.no
    @date   
    @brief  Krypton Cluster Finder for the TPC
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


//#include "AliHLTLogging.h"
//#include "AliHLTTPCPad.h"
#include "AliHLTTPCClusterFinder.h"
class AliHLTTPCSpacePointData;
class AliHLTTPCDigitReader;

class AliHLTTPCKryptonClusterFinder : public AliHLTTPCClusterFinder {

 public:
  /** standard constructor */
  AliHLTTPCKryptonClusterFinder();
  /** destructor */
  //  virtual ~AliHLTTPCKryptonClusterFinder();

  void ReBunch(const UInt_t * bunchData,Int_t bunchSize);

  void ReadDataUnsorted(void* ptr,unsigned long size);

  void FindRowClusters();

  void FindKryptonClusters();

  void CheckForCandidateOnPreviousPad(AliHLTTPCClusters* tmpCluster);

  //  Bool_t ComparePads(AliHLTTPCPad *nextPad,AliHLTTPCClusters* candidate,Int_t nextPadToRead);

 private: 
  /** copy constructor prohibited */
  AliHLTTPCKryptonClusterFinder(const AliHLTTPCKryptonClusterFinder&);
  /** assignment operator prohibited */
  AliHLTTPCKryptonClusterFinder& operator=(const AliHLTTPCKryptonClusterFinder&);

  vector<Int_t> fTimebinsInBunch;                                  //! transient

  vector<Int_t> fIndexOfBunchStart;                                //! transient

  ClassDef(AliHLTTPCKryptonClusterFinder,0) //Fast cluster finder
};
#endif
