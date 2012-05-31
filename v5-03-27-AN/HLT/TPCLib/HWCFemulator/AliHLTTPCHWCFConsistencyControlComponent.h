#ifndef ALIHLTTPCHWCFCONSISTENCYCONTROLCOMPONENT_H
#define ALIHLTTPCHWCFCONSISTENCYCONTROLCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCHWCFConsistencyControlComponent.h
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Comparison of TPC clusters produced by FPGA clusterfinder and by FPGA Emulator 
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note



#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"
class TH1F;

/**
 * @class AliHLTTPCHWCFConsistencyControlComponent
 * The FPGA clusterfinder emulator for TPC
 * Comparison of TPC clusters produced by FPGA clusterfinder and by FPGA Emulator
 * ( see AliHLTTPCHWCFEmulator class )
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCHWCFConsistenyControl <br>
 * Library: \b libAliHLTTPC
 * Input Data Types: @ref  AliHLTTPCDefinitions::fgkHWClustersDataType <br>
 * Output Data Types: @ref  <br> 
 *
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * None
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCHWCFConsistencyControlComponent : public AliHLTProcessor
{
 public:      
  /**
   * constructor 
   */
  AliHLTTPCHWCFConsistencyControlComponent();
  /** destructor */
  virtual ~AliHLTTPCHWCFConsistencyControlComponent();
  
  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);

  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

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
  AliHLTTPCHWCFConsistencyControlComponent(const AliHLTTPCHWCFConsistencyControlComponent&);

  /** assignment operator prohibited */
  AliHLTTPCHWCFConsistencyControlComponent& operator=(const AliHLTTPCHWCFConsistencyControlComponent&);

  /** Set default configuration */
  void SetDefaultConfiguration();

  /** scan configuration string */
  int ReadConfigurationString(  const char* arguments );

  /** read configuration from OCDB */
  int ReadCDBEntry( const char* cdbEntry, const char* chainId );

  /** read configuration from multiple sources */
  int Configure( const char* cdbEntry, const char* chainId, const char *commandLine );

  AliHLTUInt64_t fNDismatch;// N inconsistent data blocks
  AliHLTUInt64_t fNBlocks;// N of data blocks processed
  AliHLTComponentBenchmark fBenchmark; // benchmark
  TH1F *fHistHeaderAll; // checked parameters of block headers
  TH1F *fHistHeaderGood; // consistent parameters of block headers
  TH1F *fHistClusterAll; // checked parameters of clusters
  TH1F *fHistClusterGood; // consistent parameters of clusters
  TH1F *fProfHeader; // ratio of good headers
  TH1F *fProfCluster; // ratio of good clusters
};
#endif
