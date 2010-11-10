// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOBJECTCOMPRESSIONCOMPONENT_H
#define ALIHLTOBJECTCOMPRESSIONCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTObjectCompressionComponent.h
/// @author Matthias Richter
/// @date   2010-11-09
/// @brief  
///

#include "AliHLTProcessor.h"

/**
 * @class AliHLTObjectCompressionComponent
 * Component for compression of TObjects with a different compression level.
 * The level is adjusted by the common AliHLTComponent @ref alihlt_component_arguments.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TObjectCompressor                                       <br>
 * Library: \b libAliHLTUtil.so						    <br>
 * Input Data Types: kAliHLTAnyDataType					    <br>
 * Output Data Types: according to input blocks	                            <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -object-compression=level                                             <br>
 *      
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * Low profile: input objects are unpacked and binary copied, no streaming
 * of obejcts.
 *
 * <h2>Memory consnumption:</h2>
 * The component allocates memory of the maximum size for every input
 * object.
 *
 * <h2>Output size:</h2>
 * 
 * @ingroup alihlt_util_components
 */
class AliHLTObjectCompressionComponent : public AliHLTProcessor
{
 public:
  /// standard constructor
  AliHLTObjectCompressionComponent();
  /// destructor
  virtual ~AliHLTObjectCompressionComponent();

  /// inherited from AliHLTComponent, get component id
  virtual const char* GetComponentID() {return "TObjectCompressor";};

  /// inherited from AliHLTComponent, get the input data type
  void GetInputDataTypes( AliHLTComponentDataTypeList& );

  /// inherited from AliHLTComponent, get the output data type
  AliHLTComponentDataType GetOutputDataType();

  /// inherited from AliHLTComponent, get the output data size estimator
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /// inherited from AliHLTComponent, create a component
  virtual AliHLTComponent* Spawn() {return new AliHLTObjectCompressionComponent;}

 protected:

  /// inherited from AliHLTProcessor, data processing
  int DoEvent( const AliHLTComponentEventData& evtData,
	       AliHLTComponentTriggerData& trigData );
  
  using AliHLTProcessor::DoEvent;

  /// inherited from AliHLTComponent, component initialisation
  int DoInit( int argc, const char** argv );

  /// inherited from AliHLTComponent, scan argument
  int ScanConfigurationArgument(int argc, const char** argv);

  /// inherited from AliHLTComponent, component cleanup.
  int DoDeinit();

 private:
  /// copy constructor prohibited
  AliHLTObjectCompressionComponent(const AliHLTObjectCompressionComponent&);
  /// assignment operator prohibited
  AliHLTObjectCompressionComponent& operator=(const AliHLTObjectCompressionComponent&);

  ClassDef(AliHLTObjectCompressionComponent, 0)
};
#endif
