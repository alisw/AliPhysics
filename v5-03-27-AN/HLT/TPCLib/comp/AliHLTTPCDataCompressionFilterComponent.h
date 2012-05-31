//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCDATACOMPRESSIONFILTERCOMPONENT_H
#define ALIHLTTPCDATACOMPRESSIONFILTERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCDataCompressionFilterComponent.h
/// @author Matthias Richter
/// @date   2011-11-18
/// @brief  TPC component to filter compressed cluster data blocks from 
///         HLTOUT and in the chain.

#include "AliHLTProcessor.h"
#include <map>

/**
 * @class AliHLTTPCDataCompressionFilterComponent
 * This component can filter data compressed cluster data blocks from two
 * sources:
 * - HLTOUT
 * - parent in the chain
 *
 * It is used in an emulation chain which produces all compressed cluster
 * blocks which are missing in HLTOUT. If TPC reconstruction requires HLT
 * clusters, the emulator is automatically executed and the compressed
 * data produced if raw data is available.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCDataCompressorFilter      <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types:  <br>
 * Output Data Types: <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->

 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory consumption:</h2>
 *
 * <h2>Output size:</h2>
 *
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDataCompressionFilterComponent : public AliHLTProcessor {
public:
  /// standard constructor
  AliHLTTPCDataCompressionFilterComponent();
  /// destructor
  ~AliHLTTPCDataCompressionFilterComponent();

  /// inherited from AliHLTComponent: id of the component
  virtual const char* GetComponentID();

  /// inherited from AliHLTComponent: list of data types in the vector reference
  void GetInputDataTypes( AliHLTComponentDataTypeList& );

  /// inherited from AliHLTComponent: output data type of the component.
  AliHLTComponentDataType GetOutputDataType();

  /// inherited from AliHLTComponent: multiple output data types of the component.
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /// inherited from AliHLTComponent: output data size estimator
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /// inherited from AliHLTComponent: spawn function.
  virtual AliHLTComponent* Spawn();

protected:
  /// inherited from AliHLTProcessor: data processing
  int DoEvent( const AliHLTComponentEventData& evtData, 
	       AliHLTComponentTriggerData& trigData);
  using AliHLTProcessor::DoEvent;

  /// inherited from AliHLTComponent: component initialisation and argument scan.
  int DoInit( int argc, const char** argv );

  /// inherited from AliHLTComponent: component cleanup
  int DoDeinit();

  /// inherited from AliHLTComponent: argument scan
  int ScanConfigurationArgument(int argc, const char** argv);

  /// check the HLTOUT for availability of compressed data blocks
  int InitMapFromHLTOUT(std::map<AliHLTUInt32_t, bool>& hltoutmap);

private:
  AliHLTTPCDataCompressionFilterComponent(const AliHLTTPCDataCompressionFilterComponent&);
  AliHLTTPCDataCompressionFilterComponent& operator=(const AliHLTTPCDataCompressionFilterComponent&);

  ClassDef(AliHLTTPCDataCompressionFilterComponent, 0)
};

#endif //ALIHLTTPCDATACOMPRESSIONFILTERCOMPONENT_H
