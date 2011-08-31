//-*- Mode: C++ -*-
// $Id: AliHLTTPCHWCFProcessorUnit.h 51142 2011-08-18 13:43:40Z sgorbuno $
#ifndef ALIHLTTPCHWCFPEAKFINDERUNIT_H
#define ALIHLTTPCHWCFPEAKFINDERUNIT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *


#include "AliHLTDataTypes.h"
#include "AliHLTTPCHWCFDataTypes.h"


//  @class AliHLTTPCHWCFPeakFinderUnit
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Channel Processor unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note
//
class AliHLTTPCHWCFPeakFinderUnit
{
 public:

  /** standard constructor */
  AliHLTTPCHWCFPeakFinderUnit();
  
  /** destructor */
  ~AliHLTTPCHWCFPeakFinderUnit();

  /** set debug level */
  void SetDebugLevel( int val ){ fDebug = val; }

  /** set allowed charge fluctuation for peak finding
   */
  void SetChargeFluctuation( AliHLTUInt32_t val ){ 
    fChargeFluctuation = val;
  }


  /** initialise */
  int Init();
  
  /** input stream of data */
  int InputStream( const AliHLTTPCHWCFBunch *bunch );

  /** output stream of data */
  const AliHLTTPCHWCFBunch *OutputStream();
  
 private: 

  /** copy constructor prohibited */
  AliHLTTPCHWCFPeakFinderUnit(const AliHLTTPCHWCFPeakFinderUnit&);
  /** assignment operator prohibited */
  AliHLTTPCHWCFPeakFinderUnit& operator=(const AliHLTTPCHWCFPeakFinderUnit&);  
  

  AliHLTTPCHWCFBunch fOutput; // current output
  const AliHLTTPCHWCFBunch *fkBunch; // current input
  AliHLTUInt32_t fChargeFluctuation; // allowed charge fluctuation for peak finding 
  int fDebug; // debug level
};

#endif
