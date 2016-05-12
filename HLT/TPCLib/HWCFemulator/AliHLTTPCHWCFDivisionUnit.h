//-*- Mode: C++ -*-
// $Id$
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef ALIHLTTPCHWCFDIVISIONUNIT_H
#define ALIHLTTPCHWCFDIVISIONUNIT_H

#include "AliHLTTPCHWCFDataTypes.h"
#include "TNtuple.h"

//  @class   AliHLTTPCHWCFDivisionUnit
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Division unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note
//
class AliHLTTPCHWCFDivisionUnit
{
 public:  

  static bool CompareMCWeights(const AliHLTTPCClusterMCWeight &a, const AliHLTTPCClusterMCWeight &b){
    return a.fWeight > b.fWeight;
  }
  static bool CompareMCLabels(const AliHLTTPCClusterMCWeight &a, const AliHLTTPCClusterMCWeight &b){
    return a.fMCID < b.fMCID;
  }
  
  /** standard constructor */
  AliHLTTPCHWCFDivisionUnit();
  
  /** destructor */
  ~AliHLTTPCHWCFDivisionUnit();

  /** set debug level */
  void SetDebugLevel( int val ){ fDebug = val; }

  /** Suppress clusters wich were not mmerged (except of clusters at branch borders)
   */
  void SetSinglePadSuppression( bool val ){ fSinglePadSuppression=val; }
  
  /** Lower charge limit for clusters 
   */
  void SetClusterLowerLimit( AliHLTUInt32_t val ){ 
    fClusterLowerLimit = val << AliHLTTPCHWCFDefinitions::kFixedPoint; 
  }
  
  /** set tagging of deconvoluted clusters
   **/
  void SetTagDeconvolutedClusters( bool b ){ fTagDeconvolutedClusters = b; }

 /** initialise */
  int Init();
  
  /** input stream of data */
  int InputStream( const AliHLTTPCHWCFClusterFragment *fragment );

  /** output stream of data */
  const AliHLTTPCHWCFCluster *OutputStream();

 private: 
  
  /** copy constructor prohibited */
  AliHLTTPCHWCFDivisionUnit(const AliHLTTPCHWCFDivisionUnit&);
  /** assignment operator prohibited */
  AliHLTTPCHWCFDivisionUnit& operator=(const AliHLTTPCHWCFDivisionUnit&);  
  
  bool fSinglePadSuppression; // suppress not merged clusters
  AliHLTUInt64_t fClusterLowerLimit; // lower charge limit for clusters 
  bool fTagDeconvolutedClusters; // tag deconvoluted clusters
  const AliHLTTPCHWCFClusterFragment *fkInput; // current input 
  AliHLTTPCHWCFCluster fOutput;  // current output
  int  fDebug; // debug level
  TNtuple *fDebugNtuple;
};

#endif
