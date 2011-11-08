//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTCONFIGURATION_H
#define ALIHLTCONFIGURATION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTConfiguration.h
/// @author Matthias Richter
/// @date   
/// @brief  HLT configuration description for a single component.
/// @note   The class is used in Offline (AliRoot) context

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
 * This class describes a configuration for an HLT component by means of
 * the following parameters:
 * - configuration id:      a unique id string/name
 * - component id:          id returned by AliHLTComponent::GetComponentID()
 * - parent configuartions: ids of configurations it requires input from
 * - component arguments:   passed to the component when it is initialized
 *
 * The definition of a configuration requires simply the creation of an object
 * of type @ref AliHLTConfiguration. 
 * <pre>
 * AliHLTConfiguration myprocessor("MyProcessor", "Dummy", "publisher", "-output_percentage 80")
 * </pre>
 *
 * The Configuration is automatically registered in the list of available
 * configurations maintained by the @ref AliHLTConfigurationHandler.
 * The list is used to resolve the dependencies on other configurations.
 * Hierarchies can be built up by specifying the configuration id of parent
 * configurations as input in the .
 * A configuration entry is persistent and must be explicitly removed from
 * the AliHLTConfigurationHandler if desired.
 *
 * The registration mechanism requires the HLT system to be available. The
 * global instance of AliHLTSystem is created and retrieved by
 * <pre>
 *   // setup the HLT system
 *   AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
 * </pre>
 *
 * A configuration is transformed into a list of AliHLTTask objects by the
 * function AliHLTSystem::BuildTaskList().
 *
 * This class is only used in the HLT offline environment, see @ref alihlt_system
 * for more details.
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
   * @param bufsize    size of the output buffer in byte, the string can contain a
   *                   number prepended by a unit, e.g. 1M, allowed units 'k' and 'M'
   */
  AliHLTConfiguration(const char* id, const char* component,
		      const char* sources, const char* arguments,
		      const char* bufsize=NULL);
  /** copy constructor */
  AliHLTConfiguration(const AliHLTConfiguration& src);
  /** assignment op */
  AliHLTConfiguration& operator=(const AliHLTConfiguration& src);
  /** destructor */
  virtual ~AliHLTConfiguration();

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
  const char* GetComponentID() const {return fComponent;}

  /**
   * Return the source string.
   */
  const char* GetSourceSettings() const {return fStringSources;}

  /**
   * Return the argument string.
   */
  const char* GetArgumentSettings() const {return fArguments;}

  /**
   * Print status info.
   * Short summary on id, component, sources and unresolved sources.
   */
  void PrintStatus() const;

  /**
   * overloaded from TObject
   * options:
   *   status  - print status including the resolved sources
   */
  virtual void Print(const char* option="") const;

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
   * @return 1 if all sources resolved, 0 if not yet extracted or not resolved
   */
  int SourcesResolved() const;

  /**
   * extract the source configurations from the sources string
   * builds up the internal list of source configurations
   * @result 1 if sources resolved, 0 if not
   */
  int ExtractSources(AliHLTConfigurationHandler* pHandler);

  /**
   * Start iteration and get the first source.
   * @result pointer to the first configuration descriptor
   */
  AliHLTConfiguration* GetFirstSource() const;

  /**
   * Continue iteration and get the next source.
   * @result pointer to the next configuration descriptor in the list
   */
  AliHLTConfiguration* GetNextSource() const;

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
  int GetArguments(const char*** pArgv) const;

  /**
   * Get output buffer size.
   * @return size in byte or -1 if not specified
   */
  int GetOutputBufferSize() const {return fBufferSize;}

  /**
   * Two configurations are considered equal if all properties match
   */
  bool operator==(const AliHLTConfiguration& c) const {
    return (fID==c.fID) && (fComponent==c.fComponent) && (fStringSources==c.fStringSources) && (fArguments==c.fArguments);
  }
  bool operator!=(const AliHLTConfiguration& c) const {
    return !(*this==c);
  }

  /**
   * Helper function to build a vector from an argument string.
   * The function allocates memory for each token. The caller is responsible
   * for cleaning the strings recursively.
   * @param arg       pointer to argument string
   * @param argList   target to receive the argument list
   */
  static int InterpreteString(const char* arg, vector<char*>& argList);

 protected:
  

 private:
  /* extract arguments from the argument string
   */
  int ExtractArguments();

  /**
   * Convert buffer size string to number
   */
  int ConvertSizeString(const char* strSize) const;

  /** id of this configuration */
  TString fID;                                                     // see above
  /** component id of this configuration */
  TString fComponent;                                              // see above

  /** the <i>sources</i> string as passed to the constructor */
  TString fStringSources;                                          // see above
  /** number of resolved sources, -1 indicates re-evaluation */
  int fNofSources;                                                 //! transient
  /** list of sources */
  vector<AliHLTConfiguration*> fListSources;                       //! transient
  /** index of the current element in the list of sources */
  int fListSrcElementIdx;                                          //! transient

  /**
   * The argument string as passed to the constructor.
   * Specifies the arguments for the Analysys component. The string will
   * be parsed and the separated arguments stored in the @ref fArgv array
   * and @ref fArgc member.
   */
  TString fArguments;                                              // see above
  /** number of arguments */
  int fArgc;                                                       //! transient
  /** argument array */
  char** fArgv;                                                    //! transient

  /** size of the output buffer */
  int fBufferSize;                                                 // see above

  ClassDef(AliHLTConfiguration, 0);
};

#endif
