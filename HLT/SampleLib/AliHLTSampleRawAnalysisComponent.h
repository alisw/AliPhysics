//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTSAMPLERAWANALYSISCOMPONENT_H
#define ALIHLTSAMPLERAWANALYSISCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTSampleRawAnalysisComponent.h
/// @author Matthias Richter
/// @date   2010-08-29
/// @brief  A sample processing component for raw data
/// @ingroup alihlt_tutorial

#include "AliHLTProcessor.h"

class AliRawReaderMemory;

/**
 * @class AliHLTSampleRawAnalysisComponent
 * An example how to implement an HLT analysis component for raw data.
 * The class features the AliHLTComponent interface for HLT processing
 * components. The interface allows to run such components in either
 * the (sequential) AliSimulation/AliReconstruction framework or the
 * parallel HLT online processing framework.
 *
 * An example to run the component can be found in macro sampleRawAnalysis.C
 * in the folder HLT/exa.
 *
 * Component fetches raw data  input objects in DDL format and prints
 * some basic properties of the CommonDataHeader. It also instantiates a
 * RawReader in order to be used with some reconstruction.
 *
 * Raw data (DDL) bloxks are sent within the HLT chain with data type
 * kAliHLTDataTypeDDLRaw (see AliHLTComponentDataType). E.g Raw data
 * from TPC has data type kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC,
 * and from SPD kAliHLTDataTypeDDLRaw|kAliHLTDataOriginITSSPD. A 32 bit
 * data specification allows to differentiate between the DDL numbers
 * for one detector. For all detectors having not more than 32 links
 * the link is indicated by the corresponding bit in the specification
 * word. E.g. for link #3 bit #3 is set which gives 8 (Please note that
 * counting starts at 0).
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b SampleRawAnalysis <br>
 * Library: \b libAliHLTSample.so     <br>
 * Input Data Types: @ref kAliHLTDataTypeDDLRaw <br>
 * Output Data Types: @ref kAliHLTDataType|kAliHLTDataOriginSample
 *                         {ROOTOBAR:SMPL} <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Argument scan is implemented in the function ScanConfigurationArgument().
 * see @ref alihltcomponent-initialization-ocdb.
 * Please provide specific descriptions and implementations.
 * \li -mandatory     <i> teststring   </i> <br>
 *      
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -oprional1      <i> teststring   </i> <br>
 *      an argument with one parameter
 * \li -verbose                              <br>
 *      print some comments
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * The component has just one default CDB entry in 
 * <tt>HLT/ConfigSample/SampleRawAnalysis</tt>.
 * It does not load any configuration from the global <tt>ConfigHLT</tt>
 * folder.
 * \li -TObjString object holding a string with the configuration parameters
 *      explained above
 *
 * <h2>Performance:</h2>
 * The component does not any event data processing.
 *
 * <h2>Memory consumption:</h2>
 * The component does not any event data processing.
 *
 * <h2>Output size:</h2>
 * The component has no output data.
 *
 * @ingroup alihlt_tutorial
 */
class AliHLTSampleRawAnalysisComponent : public AliHLTProcessor {
public:
  AliHLTSampleRawAnalysisComponent();
  virtual ~AliHLTSampleRawAnalysisComponent();

  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  void GetOCDBObjectDescription( TMap* const targetMap);

  // Spawn function, return new class instance
  AliHLTComponent* Spawn();

 protected:
  // AliHLTComponent interface functions
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int ScanConfigurationArgument(int argc, const char** argv);
  int Reconfigure(const char* cdbEntry, const char* chainId);
  int ReadPreprocessorValues(const char* modules);

  using AliHLTProcessor::DoEvent;

private:
  /** copy constructor prohibited */
  AliHLTSampleRawAnalysisComponent(const AliHLTSampleRawAnalysisComponent&);
  /** assignment operator prohibited */
  AliHLTSampleRawAnalysisComponent& operator=(const AliHLTSampleRawAnalysisComponent&);

  /// verbosity level
  int fVerbosity; //! transient

  /// rawreader instance
  AliRawReaderMemory* fRawReader; //! transient

  ClassDef(AliHLTSampleRawAnalysisComponent, 0)
};
#endif
