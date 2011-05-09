// $Id: AliHLTTPCHWClusterFinderSupport.h 47875 2011-02-28 09:41:02Z richterm $

#ifndef ALIHLTTPCHWCLUSTERFINDERSUPPORT_H
#define ALIHLTTPCHWCLUSTERFINDERSUPPORT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCHWClusterFinderSupport.h
/// @author Sergey Gorbunov
/// @date   
/// @brief  Input interface for AliHLTTPCHWClusterFinderEmulator
///

#include "AliHLTLogging.h"

class AliHLTComponentBlockData;

/**
 * @class AliHLTTPCHWClusterFinderSupport 
 * The class creates input for the FPGA cluster finder emulator 
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCHWClusterFinderSupport : public AliHLTLogging
{
 public:      
  /** constructor */
  AliHLTTPCHWClusterFinderSupport();
  /** destructor */
  virtual ~AliHLTTPCHWClusterFinderSupport();

  /** method creates raw event from the HLT data block, error flag returned **/
  int CreateRawEvent( const AliHLTComponentBlockData* block, 
		      const AliHLTUInt32_t *&event, AliHLTUInt64_t &eventSize32, 
		      const AliHLTInt32_t *&mcLabels,  AliHLTUInt64_t &mcLabelsSize32 );

  int CheckRawData( const AliHLTUInt32_t *buffer, unsigned long bufferSize32, int patch, int slice );

 private:
  /** copy constructor prohibited */
  AliHLTTPCHWClusterFinderSupport(const AliHLTTPCHWClusterFinderSupport&);
  /** assignment operator prohibited */
  AliHLTTPCHWClusterFinderSupport& operator=(const AliHLTTPCHWClusterFinderSupport&);
  void Add10Word( AliHLTUInt32_t &nWords32, int &seek10, UInt_t data );

  AliHLTUInt32_t *fEventMemory;
  AliHLTInt32_t *fEventMCMemory;
};

#endif
