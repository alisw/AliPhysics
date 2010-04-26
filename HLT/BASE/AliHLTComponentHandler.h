//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTCOMPONENTHANDLER_H
#define ALIHLTCOMPONENTHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTComponentHandler.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Global handling of HLT processing components
    @note   The handler is part of the interface and both used in the
            Online (PubSub) and Offline (AliRoot) context.
                                                                          */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include <vector>
//#include "TObject.h"
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"

class AliHLTComponent;
class AliHLTModuleAgent;
struct AliHLTAnalysisEnvironment;
struct AliHLTComponentDataType;

/**
 * @class AliHLTComponentHandler
 * The component handler controls all the processing components available in
 * the system. It also controls the component shared libraries.
 * @ingroup alihlt_component
 */
class AliHLTComponentHandler : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTComponentHandler();
  /** constructor */
  AliHLTComponentHandler(AliHLTAnalysisEnvironment* pEnv);
  /** destructor */
  virtual ~AliHLTComponentHandler();

  /**
   * Create an instance from the global sigleton.
   * Instance has to be destroyed by the Destroy function
   */
  static AliHLTComponentHandler* CreateHandler();

  /**
   * Destroy an instance of the global singleton retrieved by
   * AliHLTComponentHandler::CreateHandler()
   */
  int Destroy();

  /**
   * Library mode.
   * - kDynamic: library can be unloaded (unload forced at termination of the
   *             handler
   * - kStatic:  library persistent, once loaded it stays
   */
  enum TLibraryMode {kDynamic, kStatic};

  /**
   * Set the environment for the HLT framework.
   * The environment mainly consists of function pointers for the integration
   * of the HLT framework into a system like the PubSub online system or
   * AliRoot offline system.
   * @param pEnv    pointer to @ref AliHLTAnalysisEnvironment structure
   * @return none
   */
  void SetEnvironment(AliHLTAnalysisEnvironment* pEnv);

  /**
   * Get the current environment.
   */
  const AliHLTAnalysisEnvironment* GetEnvironment() const {return &fEnvironment;}

  /**
   * Set library mode.
   * The mode effects all loaded libraries until another mode is set.
   * @param mode             persistent library or not
   * @return previous mode
   */
  TLibraryMode SetLibraryMode(TLibraryMode mode);

  /**
   * Load a component shared library.
   * The component library needs to be loaded from the ComponentHanler in order
   * to automatically register all components in the library.
   * Registration is done by passing a sample object of the component to the
   * handler. The object has to be valid during the whole runtime and should
   * thus be a global object which is ONLY used for the purpose of registration.
   * This also ensures automatically registration at library load time.
   * @param libraryPath      const char string containing the library name/path
   * @param bActivateAgents  activate agents after loading (@ref ActivateAgents)
   * @return 0 if succeeded, neg. error code if failed
   */
  int LoadLibrary( const char* libraryPath, int bActivateAgents=1);

  /**
   * Find a symbol in a dynamically loaded library.
   * @param library      library
   * @param symbol       the symbol to find
   * @return void pointer to function
   */
  AliHLTfctVoid FindSymbol(const char* library, const char* symbol);

  /**
   * Unload a component shared library.
   * All components will be de-registered when the last instance of the
   * library was unloaded.
   * @param libraryPath  library name as specified to @ref LoadLibrary
   * @return 0 if succeeded, neg. error code if failed
   */
  int UnloadLibrary( const char* libraryPath );

  /**
   * Schedule a component for registration.
   * Full registration will be done after successfull loading of the shared
   * library.
   * @param pSample  a sample object of the component
   * @return neg. error code if failed
   */
  int ScheduleRegister(AliHLTComponent* pSample );

  /**
   * Register a component.
   * Registration is done by passing a sample object of the component to the
   * handler. The object has to be valid during the whole runtime and should
   * thus be a global object which is ONLY used for the purpose of registration.
   * @param pSample   a sample object of the component
   * @return neg. error code if failed
   */
  int RegisterComponent(AliHLTComponent* pSample);

  /**
   * Add a component and leave control of the sample object to the handler.
   * Exactly the same functionality as @ref RegisterComponent but deletes
   * the sample object at clean-up of the handler.
   * @param pSample   a sample object of the component
   * @return neg. error code if failed
   */
  int AddComponent(AliHLTComponent* pSample);

  /**
   * Registers all scheduled components.
   */
  int RegisterScheduledComponents();

  /**
   * Deregister a component.
   * @param componentID   ID of the component
   * @return neg. error code if failed
   */
  int DeregisterComponent( const char* componentID );

  /**
   * Add standard components
   * The standard components are part of the libHLTbase library and
   * need therefore a special handling.
   */
  int AddStandardComponents();

  /**
   */
  int DeleteOwnedComponents();

  /**
   * Find the ID of a component with the given output data.
   * @param dtype     data type descriptor
   * @param prevType  can be used to iterate if there are multiple components
   *                  with the same output data type.
   * @return component id
   */
  //const char* FindComponentType( AliHLTComponentDataType dtype,
  //                               const char* prevType = NULL )
  //  { return NULL;}

  /**
   * Create a component of the given name (ID).
   * The method tries to find a registerd component of id \em componentID and
   * calls the \em Spawn method of the template component. After successful
   * creation of a new object, the Init method is called in order to initialize
   * the environment and the component arguments. <br>
   * The environment is the same for all components, but each component can
   * have an additional private parameter \em pEnvParam.<br>
   * The component arguments consist of an array of strings and the array size
   * in the usual manner of the main() function.
   * @param componentID  ID of the component to create
   * @param pEnvParam    environment parameter for the component
   * @param argc         number of arguments in argv
   * @param argv         argument array like in main()
   * @param component    reference to receive the create component instance
   * @return component pointer in component, neg. error code if failed
   */
  int CreateComponent( const char* componentID, void* pEnvParam, 
		       int argc, const char** argv, AliHLTComponent*& component);

  /**
   * Create component without initializing it.
   * @param componentID  ID of the component to create
   * @param component    reference to receive the create component instance
   */
  int CreateComponent(const char* componentID, AliHLTComponent*& component );

  /**
   * Create a component of the given name (ID).
   * Introduced for backward compatibility.
   * @param componentID  ID of the component to create
   * @param pEnvParam    environment parameter for the component
   * @param component    reference to receive the create component instance
   * @return component pointer in component, neg. error code if failed
   */
  int CreateComponent( const char* componentID, void* pEnvParam, 
		       AliHLTComponent*& component ) 
    {
    return CreateComponent( componentID, pEnvParam, 0, NULL, component );
    }

  /**
   * Set the run description.
   * The run description is set globally for all components. Each component
   * is initialized from the global run description after creation and before
   * call of AliHLTComponent::Init().
   *
   * @param desc    run descriptor, currently only the run no member is used
   * @param runType originally, run type was supposed to be a number and part
   *                of the run descriptor. But it was defined as string later
   */
  int SetRunDescription(const AliHLTRunDesc* desc, const char* runType);

  /**
   * Check if a registered component has output data, e.g. is of type
   * kSource or kProcessor (see @ref AliHLTComponent::TComponentType).
   * @param componentID  ID of the component to create
   * @return 1 if component has output data, 0 if not                 <br>
   *         -ENOENT     if component does not exist
   */
  int HasOutputData( const char* componentID);

  /**
   * Print registered components to stdout.
   * @return none
   */
  void List();

  /**
   * Announce version and compilation info of the base library.
   */
  int AnnounceVersion();

  /**
   * Find a component.
   * @param componentID  ID of the component to find
   * @return index, neg. error code if failed
   */
  int FindComponentIndex(const char* componentID);

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTComponentHandler(const AliHLTComponentHandler&);
  /** assignment operator prohibited */
  AliHLTComponentHandler& operator=(const AliHLTComponentHandler&);

  /**
   * Find a component.
   * @param componentID  ID of the component to find
   * @return descriptor
   */
  AliHLTComponent* FindComponent(const char* componentID);

  /**
   * Insert component to the list
   * @param pSample      sample object of the component
   * @return neg. error code if failed
   */
  int InsertComponent(AliHLTComponent* pSample);

  /**
   * Close all libraries.
   * @return neg. error code if failed
   */
  int UnloadLibraries();

  /**
   * Activate all module agents with this component handler.
   * The function loops over all available module agents and activates
   * each agent with this component handler. During activation, the
   * dynamic component registration is carried out by the agents version
   * of @ref AliHLTModuleAgent::RegisterComponents
   *
   * Agents are identified by an id which is a string containing the
   * module name. Libraries follow the naming scheme libAliHLT<MOD>.so
   * If the library name is provided and the specific agent found in the
   * list, only that one is activated. All pending agents otherwize.
   * @param library       library to activate the agent for
   * @param blackList     blank separated list of module ids
   */
  int ActivateAgents(const char* library=NULL, const char* blackList=NULL);

  /**
   * Compound descriptor for component libraries
   */
  struct AliHLTLibHandle {
    AliHLTLibHandle() : fHandle(NULL), fName(NULL), fMode(kDynamic) {}
    /** dlopen handle */
    void* fHandle;                                                 //! transient
    /** name of the library, casted to TString* before use */
    void* fName;                                                   //! transient
    /** library mode: kStatic means never unloaded */
    TLibraryMode fMode;                                            //! transient
  };

  /**
   * Find a specific library among the loaded libraries.
   * @param library     library name/path
   * @return pointer to AliHLTLibHandle
   */
  AliHLTLibHandle* FindLibrary(const char* library);

  /**
   * Unload a component shared library.
   * All components will be de-registered when the last instance of the
   * library was unloaded.
   * @param handle       handle to the library to unload
   * @return 0 if succeeded, neg. error code if failed
   */
  int UnloadLibrary(AliHLTComponentHandler::AliHLTLibHandle &handle);

  /** list of registered components */
  vector<AliHLTComponent*> fComponentList;                         // see above 
  /** list of scheduled components */
  vector<AliHLTComponent*> fScheduleList;                          // see above 
  /** list of libraries */
  vector<AliHLTLibHandle> fLibraryList;                            // see above 
  /** running environment for the component */
  AliHLTAnalysisEnvironment fEnvironment;                         // see above 
  /** list of owned components, deleted at termination of the handler */
  vector<AliHLTComponent*> fOwnedComponents;                       // see above 
  /** library mode effects all loaded libraries until a new mode is set */
  TLibraryMode fLibraryMode;                                       // see above 

  /** run descriptor */
  AliHLTRunDesc fRunDesc;                                          //!transient 
  /** run type string */
  char* fRunType;                                                  //!transient 

  /** the global singleton */
  static AliHLTComponentHandler* fgpInstance;                      //!transient
  /** number of used instances of the global singleton */
  static int fgNofInstances;                                       //!transient 

  ClassDef(AliHLTComponentHandler, 2);

};
#endif

