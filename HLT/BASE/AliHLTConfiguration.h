// @(#) $Id$

#ifndef ALIHLTCONFIGURATION_H
#define ALIHLTCONFIGURATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTConfiguration.h
    @author Matthias Richter
    @date   
    @brief  Base class and handling of HLT configurations.
    @note   The class is used in Offline (AliRoot) context
*/

#include <cerrno>
#include <vector>
#include <TObject.h>
#include <TList.h>
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTDataBuffer.h"

class AliHLTConfigurationHandler;

/**
 * @class AliHLTConfiguration
 * @brief Description of HLT processing chains.
 *
 * This class describes a certain configuration af an HLT processing step
 * by the following parameters:
 * - a unique id string/name
 * - the id of the component
 * - the ids of the configurations it requires input from
 * - the arguments, which are passed to the component when it is initialized
 *
 * The setup of a configuration requires simply the creation of a global object
 * of @ref AliHLTConfiguration. The Configuration is automatically registered
 * in the list of available configurations maintained by the @ref
 * AliHLTConfigurationHandler. The list is used by to resolve the dependencies
 * on other configurations. Hierarchies can be built up in an easy way.
 *
 * A configuration is interpreted by the @ref AliHLTConfigurationHandler and
 * transformed into a Task List.
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTConfiguration : public TObject, public AliHLTLogging {
 public:
  /**
   * standard constructor. The configuration is automatically registered in the
   * global configuration manager
   */
  AliHLTConfiguration();
  /**
   * constructor. The configuration is automatically registered in the
   * global configuration manager
   * @param id         unique id of the configuration
   * @param component  component id
   * @param sources    blank separated list of source configuration ids
   * @param arguments  argument string passed to the component at initialization
   */
  AliHLTConfiguration(const char* id, const char* component,
		      const char* sources, const char* arguments);
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTConfiguration(const AliHLTConfiguration&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTConfiguration& operator=(const AliHLTConfiguration&);
  /** destructor */
  /** destructor */
  virtual ~AliHLTConfiguration();

  /*****************************************************************************
   * global initialization
   */

  /**
   * Global initialization of the configuration handler.
   */
  static int GlobalInit(AliHLTConfigurationHandler* pHandler);

  /**
   * Global de-init and cleanup of the global configuration handler
   */
  static int GlobalDeinit();

  /*****************************************************************************
   * properties of the configuration
   */

  /**
   * Get configuration id, a unique name
   * This is an overridden TObject function in order to return the configuration
   * name instead of the class name. Enables use of TList standard functions.
   * @return configuration id
   */
  const char *GetName() const;

  /**
   * Get id of the component.
   * The id is a unique string.
   * @return id of the component
   */
  const char* GetComponentID() {return fComponent;}

  /**
   * Print status info.
   * Short summary on id, component, sources and unresolved sources.
   */
  void PrintStatus();

  /**
   * Get a certain source.
   * @param  id of the source configuration
   * @result pointer to the corresponding configuration descriptor
   */
  AliHLTConfiguration* GetSource(const char* id);

  /**
   * Try to find a dependency recursively in the list of sources.
   * @param id       the source to search for
   * @param pTgtList (optional) target list to receive the dependency tree
   * @return
   *   0 if not found
   *   n found in the n-th level
   *   dependency list in the target list  
   */
  int FollowDependency(const char* id, TList* pTgtList=NULL);

  /**
   * Get the number of resolved sources.
   * @return number of resolved sources
   */
  int GetNofSources() {return fListSources.size();}

  /**
   * Check resolving status.
   * @param bAuto resolve if ==1 
   * @return 1 if all sources resolved
   */
  int SourcesResolved(int bAuto=0);

  /**
   * Start iteration and get the first source.
   * @result pointer to the first configuration descriptor
   */
  AliHLTConfiguration* GetFirstSource();

  /**
   * Continue iteration and get the next source.
   * @result pointer to the next configuration descriptor in the list
   */
  AliHLTConfiguration* GetNextSource();

  /**
   * Invalidate a dependency and mark the configuration to be re-evaluted. 
   * @param pConf pointer to configuration descriptor
   */
  int InvalidateSource(AliHLTConfiguration* pConf);

  /**
   * Mark the configuration to be re-evaluted.
   */
  int InvalidateSources() {fNofSources=-1; return 0;}

  /**
   * Get the arguments array.
   * @param pArgv   pointer to receive argument array pointer
   * @return argc if succeeded, neg. error code if failed
   */
  int GetArguments(const char*** pArgv);

 protected:
  

 private:
  /* extract the source configurations from the sources string
   */
  int ExtractSources();

  /* extract arguments from the argument string
   */
  int ExtractArguments();

  /* helper function to build a vector from an argument string
   */
  int InterpreteString(const char* arg, vector<char*>& argList);

  /** id of this configuration */
  const char* fID;
  /** component id of this configuration */
  const char* fComponent;

  /** the <i>sources</i> string as passed to the constructor */
  const char* fStringSources;
  /** number of resolved sources, -1 indicates re-evaluation */
  int fNofSources;
  /** list of sources */
  vector<AliHLTConfiguration*> fListSources;
  /** iterator for the above list */
  vector<AliHLTConfiguration*>::iterator fListSrcElement;

  /**
   * The argument string as passed to the constructor.
   * Specifies the arguments for the Analysys component. The string will
   * be parsed and the separated arguments stored in the @ref fArgv array
   * and @ref fArgc member.
   */
  const char* fArguments;
  /** number of arguments */
  int fArgc;
  /** argument array */
  char** fArgv;

  static AliHLTConfigurationHandler* fConfigurationHandler;

  ClassDef(AliHLTConfiguration, 0);
};

#endif
