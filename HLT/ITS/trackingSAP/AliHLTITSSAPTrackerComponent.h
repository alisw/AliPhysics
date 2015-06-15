//-*- Mode: C++ -*-
// $Id$
// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

///  @file   AliHLTITSSAPTrackerComponent.h
///  @author Ruben Shahoyan <ruben.shahoyan@cern.ch>
///  @date   August 2014
///  @brief  An ITS standalone primaries tracker/vertexer processing component for the HLT
///  Adapted from HLT/ITS/tracking/AliHLTITSTrackerComponent component

#ifndef ALIHLTITSSAPTRACKERCOMPONENT_H
#define ALIHLTITSSAPTRACKERCOMPONENT_H

#include "AliHLTProcessor.h"
#include "AliHLTDataTypes.h"
#include "AliHLTComponentBenchmark.h"
#include "AliRecoParam.h"
class AliITSSAPTracker;
class TClonesArray;

/**
 * @class AliHLTITSSAPTrackerComponent
 * The HL ITS standalone primaries tracker/vertexer component.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ITSSAPTracker                              <br>
 * Library: \b libAliHLTITS.so                              <br>
 * Input Data Types:                                        <br> 
 *    kAliHLTDataTypeTrack|kAliHLTDataOriginTPC             <br>
 *    kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD       <br>
 *    kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD       <br>
 *    kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD       <br>
 *      
 * Output Data Types:                                       <br>
 *    kAliHLTDataTypeTrack|kAliHLTDataOriginITS             <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -config1      <i> teststring   </i> <br>
 *      a configuration argument with one parameter
 * \li -config2                            <br>
 *      a configuration argument without parameters
 *
 * <h2>Default CDB entries:</h2>
 *
 * ITS/Align/Data
 * ITS/Calib/SPDNoisy
 * ITS/Calib/SPDDead
 * ITS/Calib/PITConditions
 * ITS/Calib/CalibSDD
 * ITS/Calib/RespSDD
 * ITS/Calib/DriftSpeedSDD
 * ITS/Calib/DDLMapSDD
 * ITS/Calib/MapsTimeSDD
 * ITS/Calib/NoiseSSD
 * ITS/Calib/GainSSD
 * ITS/Calib/BadChannelsSSD
 * ITS/Calib/RecoParam
 *
 * <h2>Performance:</h2>
 * TODO
 *
 * <h2>Memory consumption:</h2>
 * TODO
 *
 * <h2>Output size:</h2>
 * TODO
 * 
 * @ingroup alihlt_its_components
 */
class AliHLTITSSAPTrackerComponent : public AliHLTProcessor
{
public:
  /** standard constructor */
  AliHLTITSSAPTrackerComponent();

  /** dummy copy constructor, defined according to effective C++ style */
  AliHLTITSSAPTrackerComponent( const AliHLTITSSAPTrackerComponent& );

  /** dummy assignment op, but defined according to effective C++ style */
  AliHLTITSSAPTrackerComponent& operator=( const AliHLTITSSAPTrackerComponent& );

  /** standard destructor */
  virtual ~AliHLTITSSAPTrackerComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** @see component interface @ref AliHLTComponent::GetComponentID */
  const char* GetComponentID() ;

  /** @see component interface @ref AliHLTComponent::GetInputDataTypes */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list )  ;

  /** @see component interface @ref AliHLTComponent::GetOutputDataType */
  AliHLTComponentDataType GetOutputDataType() ;

  /** @see component interface @ref AliHLTComponent::GetOutputDataType */
  int  GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /** @see component interface @ref AliHLTComponent::GetOutputDataSize */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) ;

  /** @see component interface @ref AliHLTComponent::Spawn */
  AliHLTComponent* Spawn() ;

protected:

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component.

  /** @see component interface @ref AliHLTComponent::DoInit */
  int DoInit( int argc, const char** argv );

  /** @see component interface @ref AliHLTComponent::DoDeinit */
  int DoDeinit();

  /** reconfigure **/
  int Reconfigure( const char* cdbEntry, const char* chainId );

  /** @see component interface @ref AliHLTProcessor::DoEvent */
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
	       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr,
	       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

private:
  /** type of imposed recoparam:   
      kDefault=1,kLowMult = 2,kHighMult=4,kCosmic=8,kCalib = 16 **/
  AliRecoParam::EventSpecie_t fRecoParamType; // event type imposed

  /** if fSkipSDD>=0, forbid or allow to use sdd, otherwise rely on recoparam **/
  int fSkipSDD;   // skip sdd layers even if data are shipped by HLT

  int fMaxMissL;  // max number of active layers track can miss

  int fMaxTrackletsToRun; // don't run tracking if Ntracklets above this value

  int fMaxVtxIter;        // max iteration for vertexer
  
  float fStopScaleChange; // min scale change in vertexer to stop iterations

  float fMaxRSPDVtx;      // max R of SPD vertex to accept

  AliHLTComponentBenchmark fBenchmark;// benchmark
  AliITSSAPTracker *fTracker; // the tracker itself

  /** set configuration parameters **/
  void SetDefaultConfiguration();
  int ReadConfigurationString(  const char* arguments );
  int ReadCDBEntry( const char* cdbEntry, const char* chainId );
  int Configure( const char* cdbEntry, const char* chainId, const char *commandLine  );

  TClonesArray* fClusters;

  ClassDef( AliHLTITSSAPTrackerComponent, 0 );

};
#endif
