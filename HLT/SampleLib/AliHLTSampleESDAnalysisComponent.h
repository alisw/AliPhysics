//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTSAMPLEESDANALYSISCOMPONENT_H
#define ALIHLTSAMPLEESDANALYSISCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

//  @file   AliHLTSampleESDAnalysisComponent.h
//  @author Matthias Richter
//  @date   2010-04-17
//  @brief  A sample processing component for ESD analysis.
//  @ingroup alihlt_tutorial

#include "AliHLTProcessor.h"

/**
 * @class AliHLTSampleESDAnalysisComponent
 * An example how to implement an HLT ESD analysis component.
 * The class features the AliHLTComponent interface for HLT processing
 * components. The interface allows to run such components in either
 * the (sequential) AliSimulation/AliReconstruction framework or the
 * parallel HLT online processing framework.
 *
 * An example to run the component can be found in macro sampleEsdAnalysis.C
 * in the folder HLT/exa.
 *
 * Component fetches the ESD from the input objects, loops over tracks and
 * publishes the tracks in a TObjArray of AliESDtrack.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b Sample-ESDAnalysis <br>
 * Library: \b libAliHLTSample.so     <br>
 * Input Data Types: @ref kAliHLTDataTypeESDObject <br>
 * Output Data Types: @ref kAliHLTDataTypeTObjArray|kAliHLTDataOriginSample
 *                         {ROOTOBAR:SMPL} <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Argument scan is implemented in the function ScanConfigurationArgument().
 * see @ref alihltcomponent-initialization-ocdb.
 * Please provide specific descriptions and implementations.
 * \li -mandatory1     <i> teststring   </i> <br>
 *      an argument with one parameter
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -optional1      <i> teststring   </i> <br>
 *      an argument with one parameter
 * \li -optional2                            <br>
 *      an argument without parameters
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * The component has just one default CDB entry in 
 * <tt>HLT/ConfigSample/SampleESDAnalysis</tt>.
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
 *
 * @ingroup alihlt_tutorial
 */
class AliHLTSampleESDAnalysisComponent : public AliHLTProcessor {
public:
  AliHLTSampleESDAnalysisComponent();
  virtual ~AliHLTSampleESDAnalysisComponent();

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
  AliHLTSampleESDAnalysisComponent(const AliHLTSampleESDAnalysisComponent&);
  /** assignment operator prohibited */
  AliHLTSampleESDAnalysisComponent& operator=(const AliHLTSampleESDAnalysisComponent&);

  ClassDef(AliHLTSampleESDAnalysisComponent, 0)
};
#endif
