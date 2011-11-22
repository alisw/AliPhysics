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

#include <string>
#include <TObject.h>
#include <TList.h>
#include <TString.h>
#include "AliHLTLogging.h"
#include "AliHLTConfiguration.h"
#include "AliHLTConfigurationHandler.h"
#include "AliHLTComponentHandler.h"

class AliRunLoader;
class AliRawReader;
class AliRawStream;
class AliHLTOUTHandler;
class AliHLTOUT;
class AliHLTModulePreprocessor;

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
 * data as part of <tt>AliSimulation</tt>, if only the raw reader is present,
 * on raw data as part of <tt>AliReconstruction</tt>. The configurations
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
 *       reconstruction is performed on simulated data, on raw data if Run
 *       loader is NULL and only the raw reader present.
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
  AliHLTModuleAgent(const char* id);
  /** destructor */
  virtual ~AliHLTModuleAgent();

  /**
   * Get module id.
   * The module id is a string specifying the detector, or module. The
   * library must follow the naming scheme \em libAliHLTModule.so, e.g.
   * \em libAliHLTTPC.so if the module is 'TPC'
   */
  const char* GetModuleId() const;

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
   * Get string of blank separated Module Ids
   */
  static string GetAgentIds();

  /**
   * Activate a component handler for this agent.
   * The @ref RegisterComponents method will be called in order to allow
   * the agent to register components. Once activated, the function can
   * be called repeatedly with the same handler and gently ignores the
   * invocation. In the current stage of development, only one handler
   * can be activated per agent. This is sufficient for the current
   * operation, but can be extended.
   * @param [in] pHandler  the component handler instance
   */
  int ActivateComponentHandler(AliHLTComponentHandler* pHandler);

  /**
   * Register all configurations belonging to this module with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader and the raw reader. <br>
   * The presence of Run loader and raw reader determines the mode of the
   * HLT reconstruction. If a run loader is available, reconstruction is
   * performed on simulated data, a raw reader might be available in that
   * case also. When running embedded into AliReconstruction, the Run loader
   * is always NULL and the raw gives access to data. The configurations
   * can adapt to the two cases.
   *
   * @param [in] handler   the configuration handler
   * @param [in] rawReader AliRoot RawReader instance 
   * @param [in] runloader AliRoot runloader
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
   * @param [in] rawReader AliRoot RawReader instance 
   * @param [in] runloader AliRoot runloader
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
   * @param [in] pHandler  instance of the component handler          
   */
  virtual int RegisterComponents(AliHLTComponentHandler* pHandler) const;

  /**
   * Define QA plugins
   * @return blank separated list of class names
   */
  virtual const char* GetQAPlugins() const;

  /**
   * IDs for output handlers.
   * The agent can provide output handlers in order to treat the output
   * data coming from the HLTOUT nodes.
   */
  enum AliHLTOUTHandlerType {
    kUnknownOutput =0,

    /** output is in ESD format */
    kEsd,

    /** agent provides data for a RawReader
     * From the data block one or more output blocks can be
     * created idenditcal to the ddl format. The blocks are
     * provided to subsequent analysis by a RawReader instance.
     * The data block can be processed in order to provide the
     * raw data, e.g. in case of lossless compression.
     */
    kRawReader,

    /** agent can create a raw stream
     * The agent directly generates a detector specific RawStream
     * object. This is used for pre-analyzed data which will not
     * be converted back to the raw format.
     */
    kRawStream,

    /** agent provides a chain
     * The data block is fed into an analysis chain, the treatment
     * depends on the components in the chain.
     */
    kChain,

    /** agent provides detector specific handler */
    kProprietary,
    kLastOutputHandler
  };

  /**
   * Output handler description.
   * \em fModule: module name specific for the handler type
   *              - kRawReader: DDL no printed in ascii format
   *              - kRawStream: class name of the RawStream class
   *              - kChain:     blank separated list of chains
   *              - kProprietary: name of the handler class
   */
  class AliHLTOUTHandlerDesc {
  public:
    AliHLTOUTHandlerDesc() : fHType(kUnknownOutput), fDt(kAliHLTVoidDataType), fModule() {}

    AliHLTOUTHandlerDesc(AliHLTOUTHandlerType handlerType, AliHLTComponentDataType dt, const char* module) 
      : fHType(handlerType), fDt(dt), fModule(module) {}

    AliHLTOUTHandlerDesc(const AliHLTOUTHandlerDesc& src) 
      : fHType(src.fHType), fDt(src.fDt), fModule(src.fModule) {}

    AliHLTOUTHandlerDesc& operator=(const AliHLTOUTHandlerDesc& src) {
      if (this==&src) return *this;
      fHType=src.fHType; fDt=src.fDt; fModule=src.fModule; return *this;
    }

    ~AliHLTOUTHandlerDesc() {}

    bool operator==(const AliHLTOUTHandlerType handlerType) const {
      return fHType==handlerType;
    }
    /**
     * Two descriptors are equal if all members match.
     */
    bool operator==(const AliHLTOUTHandlerDesc& desc) const {
      return fDt==desc.fDt && fHType==desc.fHType && fModule==desc.fModule;
    }

    operator AliHLTOUTHandlerType() {return fHType;}
    operator AliHLTComponentDataType() {return fDt;}

  private:
    /** type of the handler */
    AliHLTOUTHandlerType    fHType;                          //!transient
    /** data type treated by the handler */
    AliHLTComponentDataType fDt;                             //!transient
    /** class or chain name */
    TString                 fModule;                         //!transient
  };

  static const AliHLTOUTHandlerDesc fgkVoidHandlerDesc; //! initializer

  /**
   * Get handler description for a data block.
   * Depending on the data type and data specification the handler must
   * provide information
   * - if it can handle the data block, and
   * - details how it will handle it, mainly the type of the handler
   *   @ref AliHLTOUTHandlerType
   * 
   * @param [in] dt        data type of the block
   * @param [in] spec      specification of the block
   * @param [out] desc      handler description
   * @return 1 if the agent can provide a handler, 0 if not
   */
  virtual int GetHandlerDescription(AliHLTComponentDataType dt,
				    AliHLTUInt32_t spec,
				    AliHLTOUTHandlerDesc& desc) const;

  /**
   * Get handler for a data block of the HLTOUT data.
   * The agent can also provide an overloaded @ref DeleteOutputHandler
   * function to implement customized clean up. It is also possible to
   * return the same instance of a handler for different data blocks.<br>
   *
   * The framework first collects the handlers for all data blocks, and
   * calls the @ref AliHLTOUTHandler::ProcessData method afterwords for
   * each handler.
   * @param [in] dt        data type of the block
   * @param [in] spec      specification of the block
   * @return pointer to handler
   */
  virtual AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt,
					     AliHLTUInt32_t spec);

  /**
   * Delete an HLTOUT handler.
   * This is the final cleanup. The framwork makes sure that the handler is
   * not used any further outside the agent. Even if the agent returned the
   * same handler several times, cleanup is invoked only once. The default
   * implementation just deletes the object.
   * @param pInstance      pointer to handler
   */
  virtual int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  /**
   * Get raw stream for a data block.
   * @param [in] dt        data type of the block
   * @param [in] spec      specification of the block
   * @param [in] pData     data control object
   * @return Rawstream object, NULL if no Rawstream available for data type/spec
   */
  // this method is likely to be moved to a specific implementation
  // of AliHLTOUTHandler
//   virtual AliRawStream* GetRawStream(AliHLTComponentDataType dt,
// 				     AliHLTUInt32_t spec,
// 				     const AliHLTOUT* pData);

  /**
   * Get the preprocessor for this component library.
   * Create an instance of the preprocessor for this component library.
   * The caller will delete it after useage.
   * @return pointer to AliHLTModulePreprocessor object.
   */
  virtual AliHLTModulePreprocessor* GetPreprocessor();

  /**
   * Old method kept for backward compatibility, redirected to @ref
   * GetReconstructionChains.
   */
  const char* GetTopConfigurations(AliRunLoader* runloader=NULL) const {
    return GetReconstructionChains(NULL,runloader);
  }

 protected:

 private:
  /** standard constructor prohibited */
  AliHLTModuleAgent();
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
  static AliHLTModuleAgent* fgAnchor;                               //! transient

  /** next element in the list */
  AliHLTModuleAgent* fpNext;                                       //! transient

  /** the current object link (list position) */
  static AliHLTModuleAgent* fgCurrent;                              //! transient

  /** number of agents */
  static int fgCount;                                               //! transient

  /** instance of the active component handler */
  AliHLTComponentHandler* fpComponentHandler;                      //! transient

  /** id of the module */
  TString fModuleId;                                               //! transient

  ClassDef(AliHLTModuleAgent, 3);
};

#endif
