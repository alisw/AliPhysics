// @(#) $Id: AliHLTTPCHWClusterFinderEmulator.h 43577 2010-09-15 09:08:42Z sgorbuno $
// Original: AliHLTClustFinderNew.h,v 1.13 2004/06/18 10:55:26 loizides 

#ifndef ALIHLTTPCHWCLUSTERFINDEREMULATOR_H
#define ALIHLTTPCHWCLUSTERFINDEREMULATOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCHWClusterFinderEmulator.h
//  @author Anders Vestbo, Constantin Loizides
// 	    Kenneth Aamodt kenneth.aamodt@student.uib.no
//  @brief  HLT Cluster Finder for the TPC
//  @note

#include "Rtypes.h"
#include "AliHLTDataTypes.h"


class AliHLTTPCMapping;
class AliHLTTPCClusterMCData;

/**
 * @class AliHLTTPCHWClusterFinderEmulator
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCHWClusterFinderEmulator 
{
 public:  
  
  /** standard constructor */
   AliHLTTPCHWClusterFinderEmulator();
  
  /** destructor */
  virtual ~AliHLTTPCHWClusterFinderEmulator();
  
  /** setters */
  void SetDeconvPad(Bool_t f) {fDeconvPad=f;}
  void SetDeconvTime(Bool_t f) {fDeconvTime=f;}
  
  void Init( int slice, int patch );
  
  /** Loops over all rows finding the clusters */

  int FindClusters( const AliHLTUInt32_t *buffer,
		    AliHLTUInt64_t bufferSize32,
		    AliHLTUInt32_t *output,
		    AliHLTUInt64_t &outputSize32,
		    const Int_t *mcLabels,
		    AliHLTTPCClusterMCData *outputMC
		    );
  
 protected: 

  /** copy constructor prohibited */
  AliHLTTPCHWClusterFinderEmulator(const AliHLTTPCHWClusterFinderEmulator&);
  /** assignment operator prohibited */
  AliHLTTPCHWClusterFinderEmulator& operator=(const AliHLTTPCHWClusterFinderEmulator&);
  
  Bool_t fDeconvTime;      //! deconv in time direction
  Bool_t fDeconvPad;       //! deconv in pad direction
  
  
  public :   

    static const int kMaxNTimeBins = 1024+1; //   
    static const int kFixedPoint = 12; // bits for fixed point operations

 private:

    Int_t fSlice;
    Int_t fPatch;
    
    
    AliHLTTPCMapping *fMapping;
};
#endif
