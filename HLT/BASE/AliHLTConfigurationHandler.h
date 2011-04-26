//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTCONFIGURATIONHANDLER_H
#define ALIHLTCONFIGURATIONHANDLER_H
///* This file is property of and copyright by the ALICE HLT Project        * 
///* ALICE Experiment at CERN, All rights reserved.                         *
///* See cxx source for full Copyright notice                               *

/// @file   AliHLTConfigurationHandler.h
/// @author Matthias Richter
/// @date   
/// @brief  Global handling of HLT configurations.
///

#include <TList.h>

#include "AliHLTLogging.h"
class AliHLTConfiguration;
class TMap;

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
  
  /** destructor */
  virtual ~AliHLTConfigurationHandler();

  /*****************************************************************************
   * singleton handling
   */

  /**
   * Create an instance from the global sigleton.
   * Instance has to be destroyed by the Destroy function
   */
  static AliHLTConfigurationHandler* CreateHandler();

  /**
   * Destroy an instance of the global singleton retrieved by
   * AliHLTConfigurationHandler::CreateHandler()
   */
  int Destroy();

  /**
   * Get the instance of the global singleton.
   * Does not create the global instance. Returns NULL if status of the global
   * instance is 'inactive'.
   */
  static AliHLTConfigurationHandler* Instance() {
    if (!fgpInstance || !fgpInstance->IsActive()) return NULL;
    return fgpInstance;
  }

  /*****************************************************************************
   * activation, effects the availability of the global singleton via
   * AliHLTConfigurationHandler::Instance()
   */

  /// deactivate the handler, AliHLTConfiguration objects will not register
  int Deactivate() {fFlags|=kInactive; return 0;}

  /// activate the handler, AliHLTConfiguration objects will register again
  int Activate() {fFlags&=~kInactive; return 0;}

  /// check if active
  bool IsActive() const {return (fFlags&kInactive)==0;}

  /// signal a missed registration
  static int MissedRegistration(const char* name=NULL);

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

  /**
   * Print info
   * Options:
   * -  treeroot=configuration  print the dependency tree for a configuration
   * default PrintConfigurations
   */
  void Print(const char* option="");

  /**
   * Add a component substitution by component id.
   * All components of the specified component id will be replaced by the
   * substitution component, the component arguments are replaced accordingly.
   * Component substitution is in particular useful if the input to a specific
   * component should be written to file.
   */
  static int AddSubstitution(const char* componentId, const AliHLTConfiguration& subst);

  /**
   * Add a component substitution by configuration id.
   * The component of the specified configuration will be replaced by the
   * substitution component, the component arguments are replaced accordingly.
   * Component substitution is in particular useful if the input to a specific
   * component should be written to file.
   */
  static int AddSubstitution(const AliHLTConfiguration& conf , const AliHLTConfiguration& subst);

  /**
   * Find component substitution.
   */
  static const AliHLTConfiguration* FindSubstitution(const AliHLTConfiguration& conf);

 private:
  enum {
    kInactive = 0x1
  };

  /** the list of registered configurations */
  TList fgListConfigurations;                                      // see above

  /** status of the handler */
  unsigned fFlags;                                                 //! transient

  /** the global singleton */
  static AliHLTConfigurationHandler* fgpInstance;                  //!transient
  /** number of used instances of the global singleton */
  static int fgNofInstances;                                       //!transient 

  /// component substitution map
  /// key: either TObjString with component id or AliHLTConfiguration object
  static TMap* fgpSubstitutions;                                   //!transient 

  ClassDef(AliHLTConfigurationHandler, 0);
};

#endif
