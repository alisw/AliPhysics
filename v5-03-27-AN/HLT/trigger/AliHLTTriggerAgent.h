// $Id$

#ifndef ALIHLTTRIGGERAGENT_H
#define ALIHLTTRIGGERAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTTriggerAgent.h
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTTrigger library
*/

#include "AliHLTModuleAgent.h"
#include "AliHLTOUTHandler.h"

class AliESDEvent;
class TArrayC;
class TFile;
class TTree;

/**
 * @class AliHLTTriggerAgent
 * This is the agent for the AliHLTTrigger library.
 *
 * @ingroup alihlt_system
 */
class AliHLTTriggerAgent : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTTriggerAgent();
  /** destructor */
  virtual ~AliHLTTriggerAgent();

  int RegisterComponents(AliHLTComponentHandler* pHandler) const;

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

  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;

  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt, AliHLTUInt32_t spec);

  int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  /**
   * The handler for trigger decision blocks in the HLTOUT stream.
   */
  class AliHLTTriggerDecisionHandler : public AliHLTOUTHandler {
  public:
    /** constructor */
    AliHLTTriggerDecisionHandler();
    /** destructor */
    ~AliHLTTriggerDecisionHandler();

    /**
     * Process a data block.
     * Decode specification and return equipment id of the data block.
     * The data itsself i untouched.
     * @return equipment id the block should be used for.
     */
    int ProcessData(AliHLTOUT* pData);

    /** inherited from AliHLTOUTHandler */
    int GetProcessedData(const AliHLTUInt8_t* &pData);

    /** inherited from AliHLTOUTHandler */
    int ReleaseProcessedData(const AliHLTUInt8_t* pData, int size);

    /** write the temporary ESD to file */
    int WriteESD();

  private:
    /** copy constructor forbidden */
    AliHLTTriggerDecisionHandler(const AliHLTTriggerDecisionHandler&);
    /** assignment operator forbidden */
    AliHLTTriggerDecisionHandler& operator=(const AliHLTTriggerDecisionHandler&);

    AliESDEvent* fESD; //!
    TArrayC* fpData;  //!
    int fSize; //!
    TFile* fpESDfile; //!
    TTree* fpESDtree; //!
  };
 protected:

 private:
  /** copy constructor prohibited */
  AliHLTTriggerAgent(const AliHLTTriggerAgent&);
  /** assignment operator prohibited */
  AliHLTTriggerAgent& operator=(const AliHLTTriggerAgent&);

  /** handler for trigger decision blocks */
  AliHLTTriggerDecisionHandler* fTriggerDecisionHandler; //!

  ClassDef(AliHLTTriggerAgent, 1);
};

#endif
