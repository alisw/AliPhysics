//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTHANDLERCHAIN_H
#define ALIHLTOUTHANDLERCHAIN_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTOUTHandlerChain.h
    @author Matthias Richter
    @date   24.06.2008
    @brief  HLTOUT handler of type kChain.
*/

#include "AliHLTOUTHandler.h"
#include "TString.h"

class AliHLTSystem;
class AliHLTConfiguration;
class AliHLTConfigurationHandler;

/**
 * @class AliHLTOUTHandlerChain
 * The default HLTOUT handler for type kChain.
 *
 * The handler implements the kChain processing of HLTOUT data.
 * The ids of the chains to be run during processing are provided
 * as parameter to the constructor. The AliHLTModuleAgent
 * can just create a new instance and specify the chains in order
 * to define the HLTOUT handling of type kChain for a certain data
 * block. The same instance can be returned for multiple data blocks.
 * The handler will run once on all data blocks.
 *
 * The AliHLTOUTPublisherComponent must be used as data source in order
 * to publish the data blocks from HLTOUT into the chain. The component
 * publishes all data blocks selected for the handler. Additional
 * filter rules can be applied.
 *
 * <h2>Chain configuration</h2>
 * The tasks in the chain to be run can be defined either by
 * the AliHLTModuleAgent in conjunction with all other configurations or
 * by an implementation of CreateConfigurations().
 *
 * The handler is controlled by arguments passed to the constructor, the
 * syntax is equal to the AliHLTSystem (see AliHLTSystem::ScanOptions).
 *
 * <h2>Usage example:</h2>
 * An agent implementation for some sample histograms. Asumes a chain to
 * be registered with name 'SAMPLE-my-histo-converter'
 * <pre>
 *  AliHLTOUTHandler* AliHLTMyAgent::GetOutputHandler(AliHLTComponentDataType dt,
 *                                                    AliHLTUInt32_t spec)
 *  {
 *   // afterburner for some histograms
 *   if (dt==kAliHLTDataTypeHistogram|kAliHLTDataOriginSample) {
 *     return new AliHLTOUTHandlerChain("chains=SAMPLE-my-histo-converter");
 *   }
 *
 *   return NULL;
 *  }
 * </pre>
 *
 * <h2>Data output</h2>
 * The chain can produce output data as usual. All produced data blocks are
 * added to the HLTOUT. This means a chain can e.g. produce ESD data blocks
 * out of the HLT output by applying a converter component as an afterburner.
 * The produced output of the chain is automatically subject to HLTOUT
 * standard processing.
 *
 * HLTOUT processing sequence:
 * - first handlers of type kChain
 * - handlers of type kEsd
 * - handlers of type kProprietary
 *
 * @ingroup alihlt_aliroot_reconstruction
 */
class AliHLTOUTHandlerChain : public AliHLTOUTHandler {
 public:
  /** constructor */
  AliHLTOUTHandlerChain(const char* arguments);
  /** standard destructor */
  virtual ~AliHLTOUTHandlerChain();

  /**
   * Process a data block.
   * The handler runs a normal HLT chain for processing of the selected blocks.
   * The input of the chain is provided by the AliHLTOUTPublisher component.
   * @return equipment id the block should be used for.
   */
  virtual int ProcessData(AliHLTOUT* pData);
 protected:
  /**
   * Create configurations.
   * The configurations of the chain to be run can be defined either by
   * the AliHLTModuleAgent in conjunction with all other configurations or
   * by an implementation of the function.
   */
  virtual int CreateConfigurations(AliHLTConfigurationHandler* handler);

 private:
  /** standard constructor prohibited */
  AliHLTOUTHandlerChain();
  /** copy constructor prohibited */
  AliHLTOUTHandlerChain(const AliHLTOUTHandlerChain&);
  /** assignment operator prohibited */
  AliHLTOUTHandlerChain& operator=(const AliHLTOUTHandlerChain&);

  /**
   * Create and init AliHLTSystem.
   * Read the arguments and create the AliHLTOUTTask as data dump.
   */
  int InitSystem();

  TString fChains; //! transient
  TString fOptions; //! transient

  AliHLTSystem* fpSystem; //!transient
  bool fbHaveOutput; //!transient

  ClassDef(AliHLTOUTHandlerChain, 1)
};
#endif
