// $Id$

#ifndef ALIHLTCONFIGURATIONHANDLER_H
#define ALIHLTCONFIGURATIONHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTConfigurationHandler.h
    @author Matthias Richter
    @date   
    @brief  global handling of HLT configurations.
*/

/* #include <cerrno> */
#include <TObject.h>
#include <TList.h>
/* #include "AliHLTDataTypes.h" */
/* #include "AliHLTLogging.h" */
/* #include "AliHLTDataBuffer.h" */

class AliHLTConfigurationHandler : public AliHLTLogging {
 public:
  AliHLTConfigurationHandler();
  //AliHLTConfigurationHandler(AliHLTConfiguration* pConf);
  virtual ~AliHLTConfigurationHandler();

  /*****************************************************************************
   * registration
   */

  // register a configuration to the global list of configurations
  int RegisterConfiguration(AliHLTConfiguration* pConf);

  // create a configuration and register it
  int CreateConfiguration(const char* id, const char* component, const char* sources, const char* arguments);

  // remove a configuration from the global list
  int RemoveConfiguration(AliHLTConfiguration* pConf);
  int RemoveConfiguration(const char* id);

  // find a configuration from the global list
  AliHLTConfiguration* FindConfiguration(const char* id);

  // print the registered configurations to the logging function
  void PrintConfigurations();


 private:
  static TList fListConfigurations; // the list of registered configurations
  static TList fListDynamicConfigurations; // the list of dynamic configurations (for proper cleanup)

  ClassDef(AliHLTConfigurationHandler, 0);
};

#endif
