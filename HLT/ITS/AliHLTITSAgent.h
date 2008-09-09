// $Id$

#ifndef ALIHLTITSAGENT_H
#define ALIHLTITSAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTITSAgent.h
    @author Matthias Richter
    @date   25.08.2008
    @brief  Agent of the libAliHLTITS library
*/

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
   * Register all configurations belonging to the sample library with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader.
   * @param handler   [in] the configuration handler
   * @param rawReader [in] AliRoot RawReader instance 
   * @param runloader [in] AliRoot runloader
   * @return neg. error code if failed
   */
  int CreateConfigurations(AliHLTConfigurationHandler* handler,
			   AliRawReader* rawReader=NULL,
			   AliRunLoader* runloader=NULL) const;

  /**
   * Get the top configurations for local event reconstruction.
   * A top configuration describes a processing chain. It can simply be
   * described by the last configuration(s) in the chain. 
   * The agent can adapt the configurations to be registered to the current
   * AliRoot setup by checking the runloader.
   * @param rawReader [in] AliRoot RawReader instance 
   * @param runloader [in] AliRoot runloader
   * @return string containing the top configurations separated by blanks
   */
  const char* GetReconstructionChains(AliRawReader* rawReader=NULL,
				      AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const char* GetRequiredComponentLibraries() const;

  /**
   * Register components for the AliHLTSample library.
   * @param pHandler  [in] instance of the component handler          
   */
  int RegisterComponents(AliHLTComponentHandler* pHandler) const;

  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;
  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt,
				     AliHLTUInt32_t spec);
  int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  AliHLTModulePreprocessor* GetPreprocessor();

  class AliHLTOUTSDDRawDataHandler: public AliHLTOUTHandlerEquId {
  public:
    AliHLTOUTSDDRawDataHandler() {}
    ~AliHLTOUTSDDRawDataHandler() {}
    int ProcessData(AliHLTOUT* pData);
  private:
  };

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
