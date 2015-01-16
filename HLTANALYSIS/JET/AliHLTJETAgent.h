//-*- Mode: C++ -*-

// $Id: AliHLTJETAgent.h 28643 2008-09-09 20:43:24Z richterm $

#ifndef ALIHLTJETAGENT_H
#define ALIHLTJETAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTJETAgent.h
    @author Jochen Thaeder
    @date   13.02.2009
    @brief  Agent of the libAliHLTJET library
*/

#include "AliHLTModuleAgent.h"

// raw data handler of HLTOUT data
#include "AliHLTOUTHandlerEquId.h"

/**
 * @class AliHLTJETAgent
 * This is the agent for the AliHLTJet library.<br>
 *
 * @ingroup alihlt_jet
 */
class AliHLTJETAgent : public AliHLTModuleAgent {
 public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTJETAgent();
  
  /** destructor */
  virtual ~AliHLTJETAgent();

  /*
   * ---------------------------------------------------------------------------------
   *                            
   * ---------------------------------------------------------------------------------
   */

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
  Int_t CreateConfigurations(AliHLTConfigurationHandler* handler,
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
  const Char_t* GetReconstructionChains(AliRawReader* rawReader=NULL,
					AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const Char_t* GetRequiredComponentLibraries() const;

  /**
   * Register components for the AliHLTSample library.
   * @param pHandler  [in] instance of the component handler          
   */
  Int_t RegisterComponents(AliHLTComponentHandler* pHandler) const;

  /*
   * ---------------------------------------------------------------------------------
   *                            
   * ---------------------------------------------------------------------------------
   */

  /**
   *
   */
  AliHLTModulePreprocessor* GetPreprocessor();
  
  /*
   * ---------------------------------------------------------------------------------
   *                            
   * ---------------------------------------------------------------------------------
   */

  /**
   *
   */
  Int_t GetHandlerDescription(AliHLTComponentDataType dt,
			      AliHLTUInt32_t spec,
			      AliHLTOUTHandlerDesc& desc) const;

  /**
   *
   */
  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt,
				     AliHLTUInt32_t spec);

  /**
   *
   */
  Int_t DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  /*
     class AliHLTOUTSDDRawDataHandler: public AliHLTOUTHandlerEquId {
     public:
     AliHLTOUTSDDRawDataHandler() {}
     ~AliHLTOUTSDDRawDataHandler() {}
     int ProcessData(AliHLTOUT* pData);
     private:
     };
  */

 protected:

 private:

  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */
  
  /** Copy constructor prohibited */
  AliHLTJETAgent(const AliHLTJETAgent&);
  
  /** Assignment operator prohibited */
  AliHLTJETAgent& operator=(const AliHLTJETAgent&);

  /*
   * ---------------------------------------------------------------------------------
   *                                     Members
   * ---------------------------------------------------------------------------------
   */

  /** Handler for JET data in the HLTOUT stream */
  AliHLTOUTHandlerEquId* fRawDataHandler; //!transient

  /** ROOT specific member definition */
  ClassDef(AliHLTJETAgent, 0);
};

#endif
