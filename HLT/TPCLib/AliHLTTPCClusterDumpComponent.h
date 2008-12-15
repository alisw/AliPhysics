//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCCLUSTERDUMPCOMPONENT_H
#define ALIHLTTPCCLUSTERDUMPCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCClusterDumpComponent.h
    @author Kenneth Aamodt, Matthias Richter
    @date   
    @brief  Special file writer converting TPC clusters input to readable ASCII format.
*/

#include "AliHLTFileWriter.h"

/**
 * @class AliHLTTPCClusterDumpComponent
 * A converter for TPC clusters to ASCII output.
 * Data blocks of type :AliHLTTPCDefinition::fgkClustersDataType and origin 'TPC ' is
 * written in readable ASCII format to a file.
 * 
 * Component ID: \b TPCClusterDump <br>
 * Library: \b libAliHLTTPC
 *
 * See AliHLTFileWriter for arguments, further specific options
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -directory    <i> directory   </i> <br>
 *      the files will be put in.
 *      default './' 
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCClusterDumpComponent : public AliHLTFileWriter {
 public:
  /** default constructor */
  AliHLTTPCClusterDumpComponent();
  /** destructor */
  virtual ~AliHLTTPCClusterDumpComponent();

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
  AliHLTTPCClusterDumpComponent(const AliHLTTPCClusterDumpComponent&);
  /** assignment operator prohibited */
  AliHLTTPCClusterDumpComponent& operator=(const AliHLTTPCClusterDumpComponent&);

  /** slice */
  //Int_t fSlice;                                                    //! transient

  ClassDef(AliHLTTPCClusterDumpComponent, 1);
};

#endif
