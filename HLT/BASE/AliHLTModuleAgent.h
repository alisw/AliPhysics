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

class AliRunLoader;

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
 * All HLT component libraries are loaded on demand through the HLT steering
 * instance (@ref AliHLTSystem). A library can implement an agent derived 
 * from this base class, and has to define one global object of this agent
 * in the code. The agent will be registered automatically, and the features
 * can be queried when required.
 *
 * This is usually done during running the AliRoot reconstruction (see AliRoot
 * documentation on <tt> AliReconstruction</tt>). The HLT implemets the @ref
 * AliHLTReconstructor which holds the HLT steering object. Several flags can
 * be specified as options via the <tt>SetOption</tt> method of 
 * <tt>AliReconstruction</tt>, including the component libraries to be loaded.
 *
 * @section alihltmoduleagent_interface Agent interface
 * The child can implement the following functions:
 * - @ref CreateConfigurations                                              <br>
 *       create HLT configuration forming an HLT analysis chain.
 *
 * - @ref GetLocalRecConfigurations                                         <br>
 *       configurations run during local event reconstruction. 
 *       @note Local event reconstruction is the first step of the 
 *       reconstruction chain. All events are processed at once.
 *                                                                          <br>
 * - @ref GetEventRecConfigurations                                         <br>
 *       configurations run during event reconstruction.
 *       @note This is the reconstruction on event by event basis.
 *                                                                          <br>
 * - @ref GetRequiredComponentLibraries                                     <br>
 *       can indicate further libraries which are required for running the
 *       chains (e.g. if components of another library are used).
 *
 * - @ref RegisterComponents                                                <br>
 *       register componens, this can be used to avoid the component
 *       registration via global objects 
 *       @see @ref alihltcomponent-handling
 *                                                                          <br>
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
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTModuleAgent(const AliHLTModuleAgent&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTModuleAgent& operator=(const AliHLTModuleAgent&);
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
   * Register all configurations belonging to this module with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader.
   * @param handler   [in] the configuration handler
   * @param runloader [in] AliRoot runloader
   * @return neg. error code if failed
   */
  virtual int CreateConfigurations(AliHLTConfigurationHandler* handler,
				   AliRunLoader* runloader=NULL) const;

  /**
   * Get the top configurations for local event reconstruction.
   * A top configuration describes a processing chain. It can simply be
   * described by the last configuration(s) in the chain. 
   * The agent can adapt the configurations to be registered to the current
   * AliRoot setup by checking the runloader.
   * @param runloader  [in] AliRoot runloader
   * @return string containing the top configurations separated by blanks
   */
  virtual const char* GetLocalRecConfigurations(AliRunLoader* runloader=NULL) const;

  /**
   * Get the top configurations for event reconstruction.
   * The same as for @ref GetLocalRecConfigurations, but for the reconstruction on
   * event by event basis.
   * @param runloader  [in] AliRoot runloader
   * @return string containing the top configurations separated by blanks
   */
  virtual const char* GetEventRecConfigurations(AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  virtual const char* GetRequiredComponentLibraries() const;

  /**
   * Register componets.
   * This method can be used to register components for the module instead
   * of the 'static object approach'. Registration is don by passing a
   * sample object to @ref AliHLTComponentHandler::RegisterComponent<br>
   * \em Note: The sample object is owned by the agent, make sure to delete
   * it.
   */
  virtual int RegisterComponents(AliRunLoader* runloader=NULL) const;

  /**
   * Old method kept for backward compatibility, redirected to @ref
   * GetLocalRecConfigurations.
   */
  const char* GetTopConfigurations(AliRunLoader* runloader=NULL) const {
    return GetLocalRecConfigurations(runloader);
  }

 protected:

 private:
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

  ClassDef(AliHLTModuleAgent, 1);
};

#endif
