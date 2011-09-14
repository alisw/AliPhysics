//-*- Mode: C++ -*-
// $Id$
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCHWCFSupport.h
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Input interfaces for FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note

#ifndef ALIHLTTPCHWCFSUPPORT_H
#define ALIHLTTPCHWCFSUPPORT_H

#include "AliHLTLogging.h"
#include "AliHLTTPCHWCFDataTypes.h"

class AliHLTComponentBlockData;

/**
 * @class AliHLTTPCHWCFSupport 
 * The class creates input for the FPGA cluster finder emulator 
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCHWCFSupport : public AliHLTLogging
{
 public:      
  /** constructor */
  AliHLTTPCHWCFSupport();
  /** destructor */
  virtual ~AliHLTTPCHWCFSupport();
  
  /** method to read mapping file **/
  AliHLTUInt32_t *ReadMapping( int slice, int patch, const char *mappingFileName=0 ) const;
  /** method returns default mapping for given patch **/
  const AliHLTUInt32_t *GetMapping( int slice, int patch );
 
  /** method creates raw event from the HLT data block, error flag returned **/
  int CreateRawEvent( const AliHLTComponentBlockData* block, 
		      const AliHLTUInt32_t *&rawEvent, AliHLTUInt32_t &rawEventSize32,
		      const AliHLTTPCClusterMCLabel *&mcLabels,  AliHLTUInt32_t &nMCLabels );

  /** method to check raw data */
  int CheckRawData( const AliHLTUInt32_t *buffer, unsigned long bufferSize32, int patch, int slice );

  /** clean up */
  void ReleaseEventMemory();

 
 private:

  static const int fgkNSlices = 36; // n patches in TPC
  static const int fgkNPatches = 6; // n patches in TPC

  /** copy constructor prohibited */
  AliHLTTPCHWCFSupport(const AliHLTTPCHWCFSupport&);
  /** assignment operator prohibited */
  AliHLTTPCHWCFSupport& operator=(const AliHLTTPCHWCFSupport&);

  /** add 10-bit data to the 32-bit word */
  void Add10Word( AliHLTUInt32_t &nWords32, int &seek10, UInt_t data );

  AliHLTUInt32_t *fMapping[fgkNSlices][fgkNPatches]; // mapping arrays
  AliHLTUInt32_t *fEventMemory;          // memory for created event
  AliHLTTPCClusterMCLabel *fEventMCMemory; // memory for MC labels
};

#endif
