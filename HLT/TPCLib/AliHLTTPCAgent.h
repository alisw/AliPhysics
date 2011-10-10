//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCAGENT_H
#define ALIHLTTPCAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCAgent.h
//  @author Matthias Richter
//  @date   
//  @brief  Agent of the libAliHLTTPC library
//  @note

#include "AliHLTModuleAgent.h"
#include "AliHLTOUTHandlerEquId.h"

class AliHLTOUTHandlerChain;

/**
 * @class AliHLTTPCAgent
 * This is the agent for the AliHLTTPC library.
 *
 * @ingroup alihlt_system
 */
class AliHLTTPCAgent : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTTPCAgent();
  /** destructor */
  virtual ~AliHLTTPCAgent();

  /**
   * Register all configurations belonging to this module with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader.
   * @param handler      the configuration handler
   * @param rawReader    AliRawReader instance
   * @param runloader    AliRoot runloader
   * @return neg. error code if failed
   */
  int CreateConfigurations(AliHLTConfigurationHandler* handler,
			   AliRawReader* rawReader=NULL,
			   AliRunLoader* runloader=NULL) const;

  /**
   * Get the top configurations belonging to this module.
   * A top configuration describes a processing chain. It can simply be
   * described by the last configuration(s) in the chain. 
   * The agent can adapt the configurations to be registered to the current
   * AliRoot setup by checking the runloader.
   * @param rawReader    AliRawReader instance
   * @param runloader    AliRoot runloader
   * @return number of configurations, neg. error code if failed
   */
  const char* GetReconstructionChains(AliRawReader* rawReader=NULL,
				      AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const char* GetRequiredComponentLibraries() const;

  /**
   * Register components for the AliHLTTPC library.
   * @param [in] pHandler  instance of the component handler          
   */
  int RegisterComponents(AliHLTComponentHandler* pHandler) const;

  /**
   * Get handler decription for TPC data in the HLTOUT data stream.
   * @param [in]  dt        data type of the block
   * @param [in]  spec      specification of the block
   * @param [out] desc      handler description
   * @return 1 if the agent can provide a handler, 0 if not
   */
  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;

  /**
   * Get specific handler for TPC data in the HLTOUT data stream.
   * @param [in] dt        data type of the block
   * @param [in] spec      specification of the block
   * @return pointer to handler
   */
  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt,
				     AliHLTUInt32_t spec);

  /**
   * Delete an HLTOUT handler.
   * @param pInstance      pointer to handler
   */
  int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  /**
   * The handler for TPC RAW data in the HLTOUT stream.
   */
  class AliHLTTPCRawDataHandler : public AliHLTOUTHandlerEquId {
  public:
    /** constructor */
    AliHLTTPCRawDataHandler();
    /** destructor */
    ~AliHLTTPCRawDataHandler();

    /**
     * Process a data block.
     * Decode specification and return equipment id of the data block.
     * The data itsself i untouched.
     * @return equipment id the block should be used for.
     */
    int ProcessData(AliHLTOUT* pData);

  private:

  };

  enum EOptions {
    // indicate cluster id data blocks in the HLTOUT
    kHaveCompressedClusterIdDataBlock = BIT(15)
  };

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTTPCAgent(const AliHLTTPCAgent&);
  /** assignment operator prohibited */
  AliHLTTPCAgent& operator=(const AliHLTTPCAgent&);

  /** handler for TPC raw data in the HLTOUT stream */
  AliHLTTPCRawDataHandler* fRawDataHandler; //!transient

  /** afterburner for {'TRAKSEGS':'TPC '} in the HLTOUT stream */
  AliHLTOUTHandlerChain* fTracksegsDataHandler; //! transient

  /// handler for {'CLUSTERS':'TPC '}
  AliHLTOUTHandlerChain* fClustersDataHandler; //! transient

  /// handler for data blocks related to data compression
  AliHLTOUTHandlerChain* fCompressionMonitorHandler; //! transient

  ClassDef(AliHLTTPCAgent, 0);
};

#endif
