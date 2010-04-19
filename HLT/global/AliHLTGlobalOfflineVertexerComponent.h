//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTGLOBALOFFLINEVERTEXERCOMPONENT_H
#define ALIHLTGLOBALOFFLINEVERTEXERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

//  @file   AliHLTGlobalOfflineVertexerComponent.h
//  @author Matthias Richter
//  @date   2010-04-19
//  @brief  Component wrapping the offline vertexer
//  @ingroup alihlt_global

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"

class AliVertexerTracks;

/**
 * @class AliHLTGlobalOfflineVertexerComponent
 * Component wrapping the offline vertexer
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b GlobalOfflineVertexer <br>
 * Library: \b libAliHLTGlobal.so     <br>
 * Input Data Types: @ref kAliHLTDataTypeESDObject <br>
 * Output Data Types: none <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 *
 * <h2>Default CDB entries:</h2>
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
 * @ingroup alihlt_global
 */
class AliHLTGlobalOfflineVertexerComponent : public AliHLTProcessor {
public:
  AliHLTGlobalOfflineVertexerComponent();
  virtual ~AliHLTGlobalOfflineVertexerComponent();

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
  AliHLTGlobalOfflineVertexerComponent(const AliHLTGlobalOfflineVertexerComponent&);
  /** assignment operator prohibited */
  AliHLTGlobalOfflineVertexerComponent& operator=(const AliHLTGlobalOfflineVertexerComponent&);

  AliVertexerTracks* fVertexer; //! transient

  AliHLTComponentBenchmark fBenchmark; //! transient

  ClassDef(AliHLTGlobalOfflineVertexerComponent, 0)
};
#endif
