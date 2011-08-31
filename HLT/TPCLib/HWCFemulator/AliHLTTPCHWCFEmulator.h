//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCHWCFEMULATOR_H
#define ALIHLTTPCHWCFEMULATOR_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliHLTTPCHWCFExtractorUnit.h"
#include "AliHLTTPCHWCFProcessorUnit.h"
#include "AliHLTTPCHWCFMergerUnit.h"
#include "AliHLTTPCHWCFDivisionUnit.h"

class AliHLTTPCClusterMCData;
class AliHLTTPCClusterMCLabel;

//  @class   AliHLTTPCHWCFEmulator
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  FPGA ClusterFinder Emulator for TPC
//  @note
//
class AliHLTTPCHWCFEmulator 
{
 public:  

  /** standard constructor */
   AliHLTTPCHWCFEmulator();
  
  /** destructor */
  virtual ~AliHLTTPCHWCFEmulator();
   
  /** set debug level */
  void SetDebugLevel( int val ){ fDebug = val; }

  /** initialisation 
   */
  void Init( const AliHLTUInt32_t *mapping, AliHLTUInt32_t configWord1, AliHLTUInt32_t configWord2 );
  
  /** Loops over all rows finding the clusters 
   */
  int FindClusters( const AliHLTUInt32_t *rawEvent,
		    AliHLTUInt32_t rawEventSize32,
		    AliHLTUInt32_t *output,
		    AliHLTUInt32_t &outputSize32,
		    const AliHLTTPCClusterMCLabel *mcLabels,
		    AliHLTUInt32_t nMCLabels,
		    AliHLTTPCClusterMCData *outputMC
		    );

  /* useful tools */

  /** read the word written in big endian format (lowest byte first) 
   */
  static AliHLTUInt32_t ReadBigEndian( AliHLTUInt32_t word );

  /** write a word in big endian format (least byte first) 
   */
  static  AliHLTUInt32_t WriteBigEndian( AliHLTUInt32_t word );

  /** create configuration word 
   **/
  static void CreateConfiguration
    ( bool doDeconvTime, bool doDeconvPad, bool doFlowControl,  
      bool doSinglePadSuppression, bool bypassMerger, 
      AliHLTUInt32_t clusterLowerLimit,AliHLTUInt32_t singleSeqLimit, 
      AliHLTUInt32_t mergerDistance, AliHLTUInt32_t timeBinWindow, AliHLTUInt32_t chargeFluctuation,
      AliHLTUInt32_t &configWord1, AliHLTUInt32_t &configWord2  );
 
  /** create default configuration word 
   **/
  static void CreateDefaultConfiguration( AliHLTUInt32_t &configWord1, AliHLTUInt32_t &configWord2 ){
    CreateConfiguration(0,0,0,1,0,0,0, 3, 5, 0, configWord1, configWord2 );
  }
  
 private: 

  /** copy constructor prohibited */
  AliHLTTPCHWCFEmulator(const AliHLTTPCHWCFEmulator&);
  /** assignment operator prohibited */
  AliHLTTPCHWCFEmulator& operator=(const AliHLTTPCHWCFEmulator&);
 
  int  fDebug; // debug level
  const AliHLTUInt32_t *fkMapping; //! mapping array
  AliHLTTPCHWCFExtractorUnit fChannelExtractor; //! transient
  AliHLTTPCHWCFProcessorUnit fChannelProcessor; //! transient
  AliHLTTPCHWCFMergerUnit    fChannelMerger; //! transient
  AliHLTTPCHWCFDivisionUnit  fDivisionUnit;   //! transient
};
#endif
