//-*- Mode: C++ -*-
// $Id$
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef ALIHLTTPCHWCFDIVISIONUNIT_H
#define ALIHLTTPCHWCFDIVISIONUNIT_H

#include "AliHLTTPCHWCFDataTypes.h"
#include "AliHLTLogging.h"
class TNtuple;
class TFile;

//  @class   AliHLTTPCHWCFDivisionUnit
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Division unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note
//
class AliHLTTPCHWCFDivisionUnit :public AliHLTLogging
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
  void SetClusterQMaxLowerLimit( AliHLTUInt32_t val){
    fClusterQMaxLowerLimit = val << AliHLTTPCHWCFDefinitions::kFixedPoint;
  }
  
  /** set tagging of deconvoluted clusters
   **/
  void SetTagDeconvolutedClusters( AliHLTUInt32_t b ){ fTagDeconvolutedClusters = b; }
  void SetTagEdgeClusters( AliHLTUInt32_t b ){ fTagEdgeClusters = b; }

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
  AliHLTUInt64_t fClusterLowerLimit; // lower total charge limit for clusters 
  AliHLTUInt64_t fClusterQMaxLowerLimit; // lower maximum charge limit for clusters 
  AliHLTUInt32_t fTagDeconvolutedClusters; // way to tag deconvoluted clusters 
                                           // 0: no tagging 
                                           // 1: tag pad, tag time if one of the time sequences is deconvoluted 
                                           // 2: tag pad, tag time if 2 consecutive time sequences are deconvoluted
										   // 3: slightly more complicated heuristics combining option 1 and 2
  AliHLTUInt32_t fTagEdgeClusters; // tag edge clusters
  const AliHLTTPCHWCFClusterFragment *fkInput; // current input 
  AliHLTTPCHWCFCluster fOutput;  // current output
  int  fDebug; // debug level
  TNtuple *fDebugNtuple; // ntuple with some cluster parameters for debugging
  TFile * fDebugFile; // file with debug ntuple
};

#endif
