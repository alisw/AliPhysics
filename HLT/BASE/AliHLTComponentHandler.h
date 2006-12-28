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
   

#include <vector>
#include "TObject.h"
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"

class AliHLTComponent;
struct AliHLTComponentEnvironment;
struct AliHLTComponentDataType;

typedef void* AliHLTLibHandle;

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
  /** destructor */
  virtual ~AliHLTComponentHandler();

  /**
   * Set the environment for the HLT framework.
   * The environment mainly consists of function pointers for the integration
   * of the HLT framework into a system like the PubSub online system or
   * AliRoot offline system.
   * @param pEnv    pointer to @ref AliHLTComponentEnvironment structure
   * @return none
   */
  void SetEnvironment(AliHLTComponentEnvironment* pEnv);

  /**
   * Load a component shared library.
   * The component library needs to be loaded from the ComponentHanler in order
   * to automatically register all components in the library.
   * Registration is done by passing a sample object of the component to the
   * handler. The object has to be valid during the whole runtime and should
   * thus be a global object which is ONLY used for the purpose of registration.
   * This also ensures automatically registration at library load time.
   * @param libraryPath  const char string containing the library name/path
   * @return 0 if succeeded, neg. error code if failed
   */
  int LoadLibrary( const char* libraryPath );

  /**
   * Unload a component shared library.
   * All components will be de-registered.
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
   * @param pSample  a sample object of the component
   * @return neg. error code if failed
   */
  int RegisterComponent(AliHLTComponent* pSample );

  /**
   * Deregister a component.
   * @param componentID   ID of the component
   * @return neg. error code if failed
   */
  int DeregisterComponent( const char* componentID );

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
   * @param componentID  ID of the component to create
   * @param pEnv         environment for the component
   * @param argc         number of arguments in argv
   * @param argv         argument array like in main()
   * @param component    reference to receive the create component instance
   * @return component pointer in component, neg. error code if failed
   */
  int CreateComponent( const char* componentID, void* pEnv, 
		       int argc, const char** argv, AliHLTComponent*& component );

  /**
   * Create a component of the given name (ID).
   * Introduced for backward compatibility.
   * @param componentID  ID of the component to create
   * @param pEnv         environment for the component
   * @param component    reference to receive the create component instance
   * @return component pointer in component, neg. error code if failed
   */
  int CreateComponent( const char* componentID, void* pEnv, 
		       AliHLTComponent*& component ) 
    {
    return CreateComponent( componentID, pEnv, 0, NULL, component );
    }

  /**
   * Print registered components to stdout.
   * @return none
   */
  void List();

  /**
   * Announce version and compilation info of the base library.
   */
  int AnnounceVersion();

 protected:

 private:
  /**
   * Find a component.
   * @param componentID  ID of the component to find
   * @return index, neg. error code if failed
   */
  int FindComponentIndex(const char* componentID);

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

  /** list of registered components */
  vector<AliHLTComponent*> fComponentList;
  /** list of scheduled components */
  vector<AliHLTComponent*> fScheduleList;
  /** list of libraries */
  vector<AliHLTLibHandle> fLibraryList;
  /** running environment for the component */
  AliHLTComponentEnvironment fEnvironment;

  ClassDef(AliHLTComponentHandler, 0);

};
#endif

