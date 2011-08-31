//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCHWCFEMULATORCOMPONENT_H
#define ALIHLTTPCHWCFEMULATORCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCHWCFEmulatorComponent.h
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  HLT Component interface for for FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note



#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"
#include "AliHLTTPCHWCFSupport.h"
#include "AliHLTTPCHWCFEmulator.h"

class AliHLTTPCHWCFEmulator;
class AliHLTTPCDigitReader;
class AliTPCTransform;

/**
 * @class AliHLTTPCHWCFEmulatorComponent
 * The FPGA clusterfinder emulator for TPC
 * The component implements the interface methods of the @ref AliHLTProcessor
 * The actual cluster finding algorithm is implemented in @ref AliHLTTPCHWCFEmulator
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCHWClusterFinderEmulator <br>
 * Library: \b libAliHLTTPC
 * Input Data Types: @ref kAliHLTDataTypeDDLRaw or AliHLTTPCDefinitions::fgkUnpackedRawDataType <br>
 * Output Data Types: @ref AliHLTTPCDefinitions::fgkHWClustersDataType and AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo <br> 
 *
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -deconvolute <br>  
 *      Turns on deconvolution.
 * \li -do-mc <br>  
 *      Provide mc labels for found clusters
 *
 * <h2>Default CDB entries:</h2>
 * None
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCHWCFEmulatorComponent : public AliHLTProcessor
{
 public:      
  /**
   * constructor 
   */
  AliHLTTPCHWCFEmulatorComponent();
  /** destructor */
  virtual ~AliHLTTPCHWCFEmulatorComponent();
  
  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);

  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);

 protected:
  
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 
  
     /** @copydoc AliHLTComponent::DoInit
      */
    int DoInit( int argc, const char **argv );

    /** @copydoc AliHLTComponent::DoDeinit
     */
    int DoDeinit();

    /**  @copydoc @ref AliHLTComponent::Reconfigure
     */
    int Reconfigure( const char* cdbEntry, const char* chainId );

    /**  @copydoc @ref AliHLTComponent::ScanConfigurationArgument
     */
    int ScanConfigurationArgument(int argc, const char** argv);

    /** @copydoc @ref AliHLTProcessor::DoEvent
     */
    int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		 AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		 AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

     
    using AliHLTProcessor::DoEvent;
  
 private:

  /** copy constructor prohibited */
  AliHLTTPCHWCFEmulatorComponent(const AliHLTTPCHWCFEmulatorComponent&);

  /** assignment operator prohibited */
  AliHLTTPCHWCFEmulatorComponent& operator=(const AliHLTTPCHWCFEmulatorComponent&);

  /** Set default configuration */
  void SetDefaultConfiguration();

  /** scan configuration string */
  int ReadConfigurationString(  const char* arguments );

  /** read configuration from OCDB */
  int ReadCDBEntry( const char* cdbEntry, const char* chainId );

  /** read configuration from multiple sources */
  int Configure( const char* cdbEntry, const char* chainId, const char *commandLine );

  Bool_t fDoDeconvTime;            // flag to deconvolute in time direction
  Bool_t fDoDeconvPad;             // flag to deconvolute in pad direction
  Bool_t fDoMC;                    // flag to provide MC labels
  Bool_t fDoFlowControl;           // flag to control the data
  Bool_t fDoSinglePadSuppression;  // flag for single pad suppression
  Bool_t fBypassMerger;            // flag to bypass cluster merging between pads
  AliHLTUInt32_t fClusterLowerLimit; // cut clusters at this charge value
  AliHLTUInt32_t fSingleSeqLimit;    // cut sequences at this charge value
  AliHLTUInt32_t fMergerDistance; // max. distance in mean time between two pads to be merged
  AliHLTUInt32_t fTimeBinWindow; // timebin window
  AliHLTUInt32_t fChargeFluctuation; // allowed charge fluctuation for peak finding 
  Int_t fDebug; // debug level
  AliHLTTPCHWCFSupport fCFSupport;     // !transient
  AliHLTTPCHWCFEmulator fCFEmulator;   // !transient
  AliHLTComponentBenchmark fBenchmark; // benchmark
  
};
#endif
