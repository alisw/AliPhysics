//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCHWCFEXTRACTORUNIT_H
#define ALIHLTTPCHWCFEXTRACTORUNIT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCHWCFExtractorUnit.h
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Channel Extractor unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note

#include "AliHLTDataTypes.h"
#include "AliHLTTPCHWCFDataTypes.h"

class AliHLTTPCHWCFExtractorUnit
{
 public:  

  /** standard constructor */
  AliHLTTPCHWCFExtractorUnit();
  
  /** destructor */
  ~AliHLTTPCHWCFExtractorUnit();

  /** initialisation **/
  int Init( const AliHLTUInt32_t *mapping, const AliHLTTPCClusterMCLabel *mcLabels, AliHLTUInt32_t nMCLabels );

  /** input stream of data **/
  int InputStream( AliHLTUInt32_t word );

  /** input "end of data" signal **/  
  int InputEndOfData();

  /** output stream of data **/
  const AliHLTTPCHWCFBunch *OutputStream();

 private: 

  enum AliHLTTPCHWCFExtractorStatus{ kReadingData, kReadingRCU, kFinishing, kStop };  
  enum AliHLTTPCHWCFExtractorInputStatus{ kEmpty, kData, kEndOfData };  
  
  /** copy constructor prohibited */
  AliHLTTPCHWCFExtractorUnit(const AliHLTTPCHWCFExtractorUnit&);
  /** assignment operator prohibited */
  AliHLTTPCHWCFExtractorUnit& operator=(const AliHLTTPCHWCFExtractorUnit&);  

  AliHLTTPCHWCFExtractorStatus fStatus; // status of the unit
  AliHLTUInt32_t fInput;  // current input
  AliHLTTPCHWCFExtractorInputStatus fInputStatus; // input status
 
  const AliHLTUInt32_t *fkMapping; // mapping array
  AliHLTTPCHWCFBunch fMemory[2];  // memory for current bunch and pending output 
  AliHLTTPCHWCFBunch *fBunch;     // current bunch
  bool fPendingOutput; // is there something in the output buffer

  AliHLTInt32_t fChannelNumWordsLeft; // n 10-bit words left in the channel
  AliHLTInt32_t fBunchNumWordsLeft;// n 10-bit words left in the bunch
  AliHLTInt32_t fBunchCurrentTime; // timebin of the curent signal

  const AliHLTTPCClusterMCLabel *fkMCLabels; // pointer to mc labels
  AliHLTUInt32_t fNMCLabels;                 // N mc labels
  AliHLTUInt32_t fCurrentMCLabel;            // mc label to read next
};

#endif
