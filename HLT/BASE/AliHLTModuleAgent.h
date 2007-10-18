//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTMODULEAGENT_H
#define ALIHLTMODULEAGENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTModuleAgent.h
    @author Matthias Richter
    @date   
    @brief  Agent helper class for component libraries.
    @note   The class is used in Offline (AliRoot) context
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include <TObject.h>
#include <TList.h>
#include "AliHLTLogging.h"
#include "AliHLTConfiguration.h"
#include "AliHLTConfigurationHandler.h"
#include "AliHLTComponentHandler.h"

class AliRunLoader;
class AliRawReader;
class AliRawStream;
class AliHLTOUTHandler;
class AliHLTOUT;

/**
 * @class AliHLTModuleAgent
 * @brief Agent helper class for HLT sub modules, e.g. PHOS, TPC, Trigger
 *
 * This class implements the agent base class for the HLT sub modules.
 * The agent of a library gives information on the features of the library/
 * components, like the configurations to run and other component libraries
 * it depends on.
 * @note There must not be more than one agent per module/library.
 *
 * If a run loader is available, reconstruction is performed on simulated
 * data as part of <tt>AliSimulation</tt>, if the raw reader is present on
 * raw data as part of <tt>AliReconstruction</tt>. The configurations
 * can adapt to the two cases.
 *
 * All HLT component libraries are loaded on demand through the HLT steering
 * instance (@ref AliHLTSystem). A library can implement an agent derived 
 * from this base class, and has to define one global object of this agent
 * in the code. The agent will be registered automatically, and the features
 * can be queried when required.
 *
 * This is usually done during running the AliRoot reconstruction (see AliRoot
 * documentation on <tt> AliSimulation </tt> and <tt> AliReconstruction </tt>).
 * The HLT implemets the @ref AliHLTSimulation and @ref
 * AliHLTReconstructor which hold the HLT steering object. Several flags can
 * be specified as options via the <tt>SetRunHLT</tt> method of
 * <tt>AliSimulation</tt> and the <tt>SetOption</tt> method of 
 * <tt>AliReconstruction</tt>, including the component libraries to be loaded.
 *
 * @section alihltmoduleagent_interface Agent interface
 * The child can implement the following functions:
 * - @ref CreateConfigurations                                              <br>
 *       Create HLT configuration forming an HLT analysis chain.            <br>
 *       Reconstruction of raw data or simulated data from digits needs
 *       usually different configurations. If a run loader is available,
 *       reconstruction is performed on simulated data, on raw data if the
 *       raw reader is present.
 *
 * - @ref GetReconstructionChains                                           <br>
 *       Configurations run during event reconstruction.                    <br>
 *       Define chains to be run during the recunstruction step,
 *       Depending on the availability of AliRoot run loader or raw reader
 *                                                                          <br>
 *
 * - @ref GetRequiredComponentLibraries                                     <br>
 *       can indicate further libraries which are required for running the
 *       chains (e.g. if components of another library are used).
 *
 * - @ref RegisterComponents                                                <br>
 *       register componens, this can be used to avoid the component
 *       registration via global objects 
 *       @see @ref alihltcomponent-handling
 *                                                                          <br>
 * - @ref GetHandlerDescription                                             <br>
 *       the agent can announce which part of the HLTOUT data can be treated
 *       by the library and through which method. Different types of handlers
 *       are defined to fit the various formats of the HLTOUT data.
 *       @see AliHLTOUTHandlerType
 *
 * - @ref GetOutputHandler                                                  <br>
 *       Return AliHLTOUTHandler for a given data type and specification.
 *       This is mainly intended to treat detector proprietary data.
 *
 * - @ref GetRawStream                                                      <br>
 *       Return an AliRawStream object which is capable of treating the
 *       specified data type and specification. Rawstream must be provided
 *       for data blocks intended to be the input for AliRoot detector
 *       reconstruction by replacing the normal input stream.
 *
 * @section alihltmoduleagent_references References
 * @see @ref AliHLTReconstructor interface to the AliRoot reconstruction
 * @see @ref AliHLTAgentSample agent for the libAliHLTSample library
 *
 * @ingroup alihlt_system
 */
class AliHLTModuleAgent : public TObject, public AliHLTLogging {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTModuleAgent();
  /** destructor */
  virtual ~AliHLTModuleAgent();

  /**
   * Print status info.
   * Short summary on registered agents. This function acts globally on the
   * list of agents if no specific agent is specified.
   */
  static void PrintStatus(const char* agent=NULL);

  /**
   * Get the first agent in the list
   * @return  pointer to first agent in the list, NULL if empty
   */
  static AliHLTModuleAgent* GetFirstAgent();

  /**
   * Get the next agent in the list
   * @return  pointer to next agent in the list, NULL if end of list
   */
  static AliHLTModuleAgent* GetNextAgent();

  /**
   * Activate a component handler for this agent.
   * The @ref RegisterComponents method will be called in order to allow
   * the agent to register components. Once activated, the function can
   * be called repeatedly with the same handler and gently ignores the
   * invocation. In the current stage of development, only one handler
   * can be activated per agent. This is sufficient for the current
   * operation, but can be extended.
   * @param pHandler  [in] the component handler instance
   */
  int ActivateComponentHandler(AliHLTComponentHandler* pHandler);

  /**
   * Register all configurations belonging to this module with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader and the raw reader. <br>
   * Run loader and raw reader are usually not present at the same time.
   * If a run loader is available, reconstruction is performed on simulated
   * data, if the raw reader is present on raw data. The configurations
   * can adapt to the two cases.
   *
   * @param handler   [in] the configuration handler
   * @param rawReader [in] AliRoot RawReader instance 
   * @param runloader [in] AliRoot runloader
   * @return neg. error code if failed
   */
  virtual int CreateConfigurations(AliHLTConfigurationHandler* handler,
				   AliRawReader* rawReader=NULL,
				   AliRunLoader* runloader=NULL) const;

  /**
   * Get the top configurations for event reconstruction.
   * A top configuration describes a processing chain. It can simply be
   * described by the last configuration(s) in the chain. 
   * The agent can adapt the configurations to be registered to the current
   * AliRoot setup by checking the run loader and the raw reader.
   * @param rawReader [in] AliRoot RawReader instance 
   * @param runloader [in] AliRoot runloader
   * @return string containing the top configurations separated by blanks
   */
  virtual const char* GetReconstructionChains(AliRawReader* rawReader=NULL,
					      AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on. <br>
   * @note This is not the right place to specify libraries which this component
   * library depends. Dependencies must be linked or loaded before.
   * @return list of component libraries as a blank-separated string.
   */
  virtual const char* GetRequiredComponentLibraries() const;

  /**
   * Register components.
   * This method can be used to register components for the module instead
   * of the 'static object approach'. Registration is done by passing a
   * sample object to the AliHLTComponentHandler via
   * - @ref AliHLTComponentHandler::RegisterComponent                      <br>
   *        The sample object is owned by the agent, make sure to delete it.
   * - @ref AliHLTComponentHandler::AddComponent                           <br>
   *        Same functionality but handler deletes the object at the end.
   *
   * @param pHandler  [in] instance of the component handler          
   */
  virtual int RegisterComponents(AliHLTComponentHandler* pHandler) const;

  /**
   * IDs for output handlers.
   * The agent can provide output handlers in order to treat the output
   * data coming from the HLTOUT nodes.
   */
  enum AliHLTOUTHandlerType {
    kUnknownOutput =0,
    /** output is in ESD format */
    kEsd,
    /** agent can create a raw stream */
    kRawstream,
    /** agent provides a chain */
    kChain,
    /** agent provides detector specific handler */
    kProprietary,
    kLastOutputHandler
  };

  /**
   * Output handler description.
   * \em fModule: module name specific for the handler type
   *              - kRawStream: class name of the Rawstream class
   *              - kChain:     blank separated list of chains
   *              - kProprietary: name of the handler class
   */
  struct AliHLTOUTHandlerDesc {
    /** type of the handler */
    AliHLTOUTHandlerType    fHType;
    /** data type treated by the handler */
    AliHLTComponentDataType fDt;
    /** class or chain name */
    const char*             fModule;
  };

  /**
   * Get handler description for a data block.
   * @param dt        [in] data type of the block
   * @param spec      [in] specification of the block
   * @param desc      [out] handler description
   * @return 1 if the agent can provide a handler, 0 if not
   */
  virtual int GetHandlerDescription(AliHLTComponentDataType dt,
				    AliHLTUInt32_t spec,
				    AliHLTOUTHandlerDesc& desc) const;
  /**
   * Get handler for a data block of the HLTOUT data.
   * @param dt        [in] data type of the block
   * @param spec      [in] specification of the block
   */
  virtual AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt,
					     AliHLTUInt32_t spec) const;

  /**
   * Get raw stream for a data block.
   * @param dt        [in] data type of the block
   * @param spec      [in] specification of the block
   * @param pData     [in] data control object
   * @return Rawstream object, NULL if no Rawstream available for data type/spec
   */
  virtual AliRawStream* GetRawStream(AliHLTComponentDataType dt,
				     AliHLTUInt32_t spec,
				     const AliHLTOUT* pData) const;

  /**
   * Old method kept for backward compatibility, redirected to @ref
   * GetReconstructionChains.
   */
  const char* GetTopConfigurations(AliRunLoader* runloader=NULL) const {
    return GetReconstructionChains(NULL,runloader);
  }

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTModuleAgent(const AliHLTModuleAgent&);
  /** assignment operator prohibited */
  AliHLTModuleAgent& operator=(const AliHLTModuleAgent&);

  /**
   * Register agent in the global list.
   * @return neg. error code if failed
   */
  static int Register(AliHLTModuleAgent* pAgent);

  /**
   * Unregister agent in the global list.
   * @return neg. error code if failed
   */
  static int Unregister(AliHLTModuleAgent* pAgent);

  /** the list of active agents */
  static AliHLTModuleAgent* fAnchor;                               //! transient

  /** next element in the list */
  AliHLTModuleAgent* fpNext;                                       //! transient

  /** the current object link (list position) */
  static AliHLTModuleAgent* fCurrent;                              //! transient

  /** number of agents */
  static int fCount;                                               //! transient

  /* instance of the active component handler */
  AliHLTComponentHandler* fpComponentHandler;                      //! transient

  ClassDef(AliHLTModuleAgent, 2);
};

#endif
