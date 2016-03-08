// $Id$

#ifndef ALIHLT_EXTERNALINTERFACE_H
#define ALIHLT_EXTERNALINTERFACE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTExternalInterface.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Pure and dynamic C interface to the AliRoot HLT analysis
    @note   Utilized by the HLT Online (PubSub) framework
*/

/** 
 * @defgroup alihlt_wrapper_interface The HLT wrapper interface
 * The wrapper interface is a pure C interface which allows to use the 
 * analysis components in external applications. The interface is utilized
 * to bind the analysis code to the PubSub framework. 
 *
 * \image html PubSub_WrapperComponent.png "Wrapper interface"
 *
 * @section alihlt_wrapper_interface_general Interface functions
 * The interface is based on function signatures. The function
 * AliHLTAnalysisGetInterfaceCall(const char* signature) of type
 * AliHLTAnalysisFctGetInterfaceCall looks for a matching interface
 * function and returns the pointer if it is found.
 * @ref ALIHLTANALYSIS_FCT_GETINTERFACECALL defines the function name. 
 *
 * @section alihlt_wrapper_interface_usage Usage
 * @subsection alihlt_wrapper_interface_usage_systemcalls Getting the interface functions
 * <pre>
 * string libraryPath=ALIHLTANALYSIS_INTERFACE_LIBRARY;
 *
 * string libraryPath=gBasePath;
 * libraryPath+="/";
 * libraryPath+=ALIHLTANALYSIS_INTERFACE_LIBRARY;
 *
 * void* libHandle=dlopen(libraryPath.c_str(), RTLD_NOW);
 * if (!libHandle) {
 *   cerr << "error: can not load library " << libraryPath.c_str() << endl;
 *   return -1;
 * }
 *
 * AliHLTAnalysisFctGetInterfaceCall fctGetSystemCall=(AliHLTAnalysisFctGetInterfaceCall)dlsym(libHandle, ALIHLTANALYSIS_FCT_GETINTERFACECALL);
 * if (!fctGetSystemCall) {
 *   cerr << "error: can not find function '" << ALIHLTANALYSIS_FCT_GETINTERFACECALL << "' in " << libraryPath.c_str() << endl;
 *   return -1;
 * }
 *
 * </pre>
 *
 * @subsection alihlt_wrapper_interface_usage_init System initialization
 * <pre>
 * AliHLTAnalysisEnvironment environment;
 * memset(&environment, 0, sizeof(environment));
 *
 * // setting function pointers
 * environment.fStructSize=sizeof(environment);
 * environment.fAllocMemoryFunc=AllocMemory;
 * environment.fLoggingFunc=Logging;
 *
 * AliHLTExtFctInitSystem fctInitSystem=(AliHLTExtFctInitSystem)fctGetSystemCall("int AliHLTAnalysisInitSystem(unsigned long,AliHLTAnalysisEnvironment*,unsigned long,const char*)");
 * if (!fctInitSystem) {
 *   cerr << "error: missing AliHLTAnalysisInitSystem call" << endl;
 *   return -1;
 * }
 *
 * if ((iResult=fctInitSystem( ALIHLT_DATA_TYPES_VERSION, &environment, 0xbeef, "dummy-run" ))<0) {
 *   cerr << "InitSystem failed with " << iResult << endl;
 *   return iResult;
 * }
 *
 * </pre>
 *
 * @subsection alihlt_wrapper_interface_usage_create Load library and create component
 * <pre>
 * AliHLTExtFctLoadLibrary fctLoadLibrary=(AliHLTExtFctLoadLibrary)fctGetSystemCall("int AliHLTAnalysisLoadLibrary(const char*)");
 * if (!fctLoadLibrary) {
 *   cerr << "error: missing LoadLibrary call" << endl;
 *   return -1;
 * }
 *
 * if ((iResult=fctLoadLibrary(moduleLibrary))<0) {
 *   cerr << "error: AliHLTAnalysisLoadLibrary failed with " << iResult << endl;
 *   return iResult;
 * }
 *
 * AliHLTExtFctCreateComponent fctCreateComponent=(AliHLTExtFctCreateComponent)fctGetSystemCall("int AliHLTAnalysisCreateComponent(const char*,void*,int,const char**,AliHLTComponentHandle*,const char*)");
 * if (!fctCreateComponent) {
 *   cerr << "error: missing CreateComponent call" << endl;
 *   return -1;
 * }
 *
 * AliHLTComponentHandle handle;
 * if ((iResult=fctCreateComponent("TestProcessor", &gDummy, 0, NULL, &handle, "-chainid=test" ))<0) {
 *   cerr << "error: AliHLTAnalysisCreateComponent failed with " << iResult << endl;
 *   return iResult;
 * }
 *
 * </pre>
 * 
 * @section alihlt_wrapper_interface_cdb CDB handling
 * The interface initializes the CDB from the path found
 * in the environment variable ALIHLT_HCDBDIR. If this is empty, path is
 * set from <tt>$ALICE_ROOT</tt>. If ALIHLT_HCDBSNAPSHOT is set, this is set as snaphsot file on top of ALIHLT_HCDBDIR as default storage.
 */

/////////////////////////////////////////////////////////////////////////////////////
//
// AliHLT external interface functions
//

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Get a system call of the interface.
   * @param function signature
   * @return pointer to system call
   * @ingroup alihlt_wrapper_interface
   */
  void* AliHLTAnalysisGetInterfaceCall(const char* function);

#ifdef __cplusplus
}
#endif


/////////////////////////////////////////////////////////////////////////////////////
//
// AliHLTSystem interface functions
//

class AliHLTSystem;
class AliHLTOUT;
class AliESDEvent;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Set options for an AliHLTSystem instance.
 * The function is introduced for the sake of backward compatibility.
 * Called from AliHLTReconstructor, which loads the function dynamically.
 * @return neg. error code if failed                                     <br>
 *         -EFAULT    type cast failed                                   <br>
 *         -EINVAL    invalid parameter
 * @ingroup alihlt_system_interface
 */
int AliHLTSystemSetOptions(AliHLTSystem*, const char*);

/**
 * Process the HLTOUT data with the specified system instance.
 * The function is introduced for the sake of backward compatibility.
 * Called from AliHLTReconstructor, which loads the function dynamically.
 * @return neg. error code if failed                                     <br>
 *         -EFAULT    type cast failed                                   <br>
 *         -EINVAL    invalid parameter
 * @ingroup alihlt_system_interface
 */
int AliHLTSystemProcessHLTOUT(AliHLTSystem* pInstance, AliHLTOUT* pHLTOUT, AliESDEvent* esd);

#ifdef __cplusplus
}
#endif

#endif //ALIHLT_EXTERNALINTERFACE_H 
