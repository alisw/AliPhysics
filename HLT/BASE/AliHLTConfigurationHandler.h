// $Id$

#ifndef ALIHLTCONFIGURATIONHANDLER_H
#define ALIHLTCONFIGURATIONHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTConfigurationHandler.h
    @author Matthias Richter
    @date   
    @brief  Global handling of HLT configurations.
*/

#include <TList.h>

/**
 * @class AliHLTConfigurationHandler
 * @brief Global Handling of HLT configurations.
 *
 * This class implements the global handling of @ref AliHLTConfiguration objects.
 * It is a list of all configuartion descriptor currently available in the system.
 * Each @ref AliHLTConfiguration object is registerd automatically with the
 * handler and put into the list.
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTConfigurationHandler : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTConfigurationHandler();
  
  //AliHLTConfigurationHandler(AliHLTConfiguration* pConf);

  /** destructor */
  virtual ~AliHLTConfigurationHandler();

  /*****************************************************************************
   * registration
   */

  /**
   * Register a configuration to the global list of configurations.
   * @param pConf     The configuration to register
   */
  int RegisterConfiguration(AliHLTConfiguration* pConf);

  /**
   * Create a configuration and register it.
   * @param id
   * @param component
   * @param sources
   * @param arguments
   */
  int CreateConfiguration(const char* id, const char* component, const char* sources, const char* arguments);

  /**
   * Remove a configuration from the global list.
   * @param pConf     The configuration to remove
   */
  int RemoveConfiguration(AliHLTConfiguration* pConf);

  /**
   * Remove a configuration from the global list.
   * @param id     The configuration to remove
   */
  int RemoveConfiguration(const char* id);

  /**
   * Find a configuration from the global list.
   * @param id     Id of the configuration to find
   */
  AliHLTConfiguration* FindConfiguration(const char* id);

  /**
   * Print the registered configurations to the logging function.
   */
  void PrintConfigurations();


 private:
  /** the list of registered configurations */
  static TList fgListConfigurations;                              // see above
  /** the list of dynamic configurations (for proper cleanup) */
  static TList fgListDynamicConfigurations;                       // see above

  ClassDef(AliHLTConfigurationHandler, 0);
};

#endif
