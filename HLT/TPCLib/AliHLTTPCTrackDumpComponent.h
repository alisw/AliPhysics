//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCTRACKDUMPCOMPONENT_H
#define ALIHLTTPCTRACKDUMPCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCTrackDumpComponent.h
    @author Gaute Ovrebekk
    @date   
    @brief  Special file writer converting TPC tracks input to ASCII.
*/

#include "AliHLTFileWriter.h"

/**
 * @class AliHLTTPCTrackDumpComponent
 * A converter for track data of the TPC to ASCII output.
 * Data blocks of type fgkTrackSegmentsDataType or fgkTracksDataType is docoded
 * written in readable ASCII format to a file.
 * 
 * Component ID: \b TPCTrackDump <br>
 * Library: \b libAliHLTTPC
 *
 * See AliHLTFileWriter for arguments, further specific options
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCTrackDumpComponent : public AliHLTFileWriter {
 public:
  /** default constructor */
  AliHLTTPCTrackDumpComponent();
  /** destructor */
  virtual ~AliHLTTPCTrackDumpComponent();

  // interface functions: property getters
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  virtual AliHLTComponent* Spawn();

 protected:
  // interface functions: processing
  int InitWriter();
  int CloseWriter();
  int DumpEvent( const AliHLTComponentEventData& evtData,
		 AliHLTComponentTriggerData& trigData );
  using AliHLTDataSink::DumpEvent;

  int ScanArgument(int argc, const char** argv);

 private:
  /** copy constructor prohibited */
  AliHLTTPCTrackDumpComponent(const AliHLTTPCTrackDumpComponent&);
  /** assignment operator prohibited */
  AliHLTTPCTrackDumpComponent& operator=(const AliHLTTPCTrackDumpComponent&);

  int PrintTrack(const AliHLTComponentEventData& evtData,const AliHLTComponentBlockData* bl,Int_t &nT);

  ClassDef(AliHLTTPCTrackDumpComponent,0);
};

#endif
