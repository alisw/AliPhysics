// @(#) $Id$

#ifndef ALIHLTCOMPONENTHANDLER_H
#define ALIHLTCOMPONENTHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHltcomponentHandler
   handler of HLT processing components
 */

#include "TObject.h"
#include "AliHLTDataTypes.h"

class AliHLTComponent;
struct AliHLTComponentEnvironment;
struct AliHLTComponent_DataType;

typedef void* AliHLTLibHandle;

class AliHLTComponentHandler {
 public:
  AliHLTComponentHandler();
  virtual ~AliHLTComponentHandler();

  void SetEnvironment(AliHLTComponentEnvironment* pEnv);

  // Load a component shared library
  int LoadLibrary( const char* libraryPath );
  int UnloadLibrary( const char* libraryPath );

  /* Component registration funcions
   * registration is done by passing a sample object of the component to the handler
   * the object has to be valid during the whole runtime and should thus be a global object
   */
  // Schedule a component for registration, full registration will be done
  // after successfull loading of the shared library
  int ScheduleRegister(AliHLTComponent* pSample );

  // Register a component with the list of available components
  int RegisterComponent(AliHLTComponent* pSample );
  int DeregisterComponent( const char* componentID );

  // Find the ID of a component with the given output data
  // prevType can be used to iterate if there are multiple components with the same output data type.
  const char* FindComponentType( AliHLTComponent_DataType, const char* prevType = NULL ) { return NULL;}

  // Create a component of the given name
  int CreateComponent( const char* componentType, void* environ_param, int argc, const char** argv, AliHLTComponent*& component );
  int CreateComponent( const char* componentType, void* environ_param, AliHLTComponent*& component ) {
    return CreateComponent( componentType, environ_param, 0, NULL, component );
  }

  /* print registered components to stdout
   */
  void List();
 protected:
  int Logging( AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message, ... );

 private:
  /* find a component
     return index
  */
  int FindComponentIndex(const char* componentID);

  /* find a component
     return descriptor
  */
  AliHLTComponent* FindComponent(const char* componentID);

  int InsertComponent(AliHLTComponent* pSample);

  // close all libraries
  int UnloadLibraries();

  /* list of registered components
   */
  vector<AliHLTComponent*> fComponentList;

  /* list of scheduled components
   */
  vector<AliHLTComponent*> fScheduleList;

  /* list of libraries
   */
  vector<AliHLTLibHandle> fLibraryList;

  AliHLTComponentEnvironment fEnvironment;

  ClassDef(AliHLTComponentHandler, 0)
    };
#endif

