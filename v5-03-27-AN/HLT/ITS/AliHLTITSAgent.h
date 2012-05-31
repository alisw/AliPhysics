//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTITSAGENT_H
#define ALIHLTITSAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTITSAgent.h
//  @author Matthias Richter
//  @date   25.08.2008
//  @brief  Agent of the libAliHLTITS library
//  @note

#include "AliHLTModuleAgent.h"

// raw data handler of HLTOUT data
#include "AliHLTOUTHandlerEquId.h"

/**
 * @class AliHLTITSAgent
 * This is the agent for the AliHLTSample library.<br>
 *
 * The agent implements the HLTOUT handling of raw data blocks from the
 * ITS SDD.
 * This assumes that the data blocks are sent with data type
 * {DDL_RAW :ISDD} and the bit set in the specification corresponding.
 * to detector DDL id.
 * An HLTOUT handler is implemented to extract the equipment id from
 * the specification.
 *
 * @ingroup alihlt_its
 */
class AliHLTITSAgent : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTITSAgent();
  /** destructor */
  virtual ~AliHLTITSAgent();

  /**
   * Inherited from AliHLTModuleAgent
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
   * Inherited from AliHLTModuleAgent
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
   * Inherited from AliHLTModuleAgent
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const char* GetRequiredComponentLibraries() const;

  /**
   * Inherited from AliHLTModuleAgent
   * Register components for the AliHLTSample library.
   * @param pHandler  [in] instance of the component handler          
   */
  int RegisterComponents(AliHLTComponentHandler* pHandler) const;

  /**
   * Inherited from AliHLTModuleAgent
   * Return HLTOUT handler description for a certain data block in the
   * HLTOUT payload.
   * @return 1 if module provides a handler
   */
  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;

  /**
   * Inherited from AliHLTModuleAgent
   * Create HLTOUT handler for a certain data block in the
   * HLTOUT payload. The same handler can be returned multiple
   * times, even for different data types. The framework will
   * handle this.
   * @return pointer to handler
   */
  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt,
				     AliHLTUInt32_t spec);

  /**
   * Inherited from AliHLTModuleAgent
   * Delete the instance of the handler. This is only called once
   * even if the same handler has been issued multiple times
   */
  int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  AliHLTModulePreprocessor* GetPreprocessor();

  class AliHLTOUTSDDRawDataHandler: public AliHLTOUTHandlerEquId {
  public:
    AliHLTOUTSDDRawDataHandler() {}
    ~AliHLTOUTSDDRawDataHandler() {}
    int ProcessData(AliHLTOUT* pData);
  private:
  };

  /**
   * Create configurations for all CFs for the ITS subdetectors
   * Get the necessary information from AliDAQ using the detector id.
   * @param pHandler       Instance of the configuration handler 
   * @param detectorId     Id of the detector as specified in AliDAQ
   * @param output         Target string to receive the configurations
   * @return neg. error code f failed
   */
  int CreateCFConfigurations(AliHLTConfigurationHandler* pHandler, int detectorId, TString& output) const;

 protected:

 private:

  /** copy constructor prohibited */
  AliHLTITSAgent(const AliHLTITSAgent&);
  /** assignment operator prohibited */
  AliHLTITSAgent& operator=(const AliHLTITSAgent&);

  /** handler for ITS raw data in the HLTOUT stream */
  AliHLTOUTHandlerEquId* fRawDataHandler; //!transient

  /** ROOT specific member definition */
  ClassDef(AliHLTITSAgent, 0);
};

#endif
