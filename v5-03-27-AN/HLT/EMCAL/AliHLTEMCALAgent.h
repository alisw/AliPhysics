
#ifndef ALIHLTEMCALAGENT_H
#define ALIHLTEMCALAGENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTEMCALAgent.h
    @author Federico Ronchetti
    @date   
    @brief  Agent of the libAliHLTEMCAL library
*/

#include "AliHLTModuleAgent.h"
#include "AliHLTOUTHandlerEquId.h"

class AliHLTOUTHandlerChain;

/**
 * @class AliHLTEMCALAgent
 * This is the agent for the AliHLTEMCAL library.
 *
 * @ingroup alihlt_system
 */
class AliHLTEMCALAgent : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTEMCALAgent();
  /** destructor */
  virtual ~AliHLTEMCALAgent();

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
   * Register components for the AliHLTEMCAL library.
   * @param pHandler  [in] instance of the component handler          
   */
  int RegisterComponents(AliHLTComponentHandler* pHandler) const;

  /**
   * Get handler decription for EMCAL data in the HLTOUT data stream.
   * @param dt        [in] data type of the block
   * @param spec      [in] specification of the block
   * @param desc      [out] handler description
   * @return 1 if the agent can provide a handler, 0 if not
   */
  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;

  /**
   * Get specific handler for EMCAL data in the HLTOUT data stream.
   * @param dt        [in] data type of the block
   * @param spec      [in] specification of the block
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
   * The handler for EMCAL RAW data in the HLTOUT stream.
   */
  class AliHLTEMCALRawDataHandler : public AliHLTOUTHandlerEquId {
  public:
    /** constructor */
    AliHLTEMCALRawDataHandler();
    /** destructor */
    ~AliHLTEMCALRawDataHandler();

    /**
     * Process a data block.
     * Decode specification and return equipment id of the data block.
     * The data itsself i untouched.
     * @return equipment id the block should be used for.
     */
    int ProcessData(AliHLTOUT* pData);

  private:

  };

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTEMCALAgent(const AliHLTEMCALAgent&);
  /** assignment operator prohibited */
  AliHLTEMCALAgent& operator=(const AliHLTEMCALAgent&);

  /** handler for EMCAL raw data in the HLTOUT stream */
  AliHLTEMCALRawDataHandler* fRawDataHandler; //!transient

  ClassDef(AliHLTEMCALAgent, 1);
};

#endif
