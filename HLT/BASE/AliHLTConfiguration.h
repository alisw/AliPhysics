// @(#) $Id$

#ifndef ALIHLTCONFIGURATION_H
#define ALIHLTCONFIGURATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTConfiguration.h
    @author Matthias Richter
    @date   
    @brief  base class and handling of HLT configurations.
*/

#include <cerrno>
#include <TObject.h>
#include <TList.h>
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTDataBuffer.h"

class AliHLTConfigurationHandler;

/**
 * @class AliHLTConfiguration
 * @brief Description of HLT processing chains.
 * @note Definition:
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
 * @ingroup AliHLTbase
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

struct AliHLTComponent_BlockData;
class AliHLTComponent;
class AliHLTComponentHandler;

/******************************************************************************/

 /**
  * @class AliHLTTask
  * A task collects all the information which is necessary to process a certain
  * step in the HLT data processing chain:
  * - the instance of the component
  * - the data buffer which receives the result of the component and provides
  *   the data to other tasks/components
  * - a list of all dependencies
  * - a list of consumers
  * - the task object holds the configuration object 
  */
class AliHLTTask : public TObject, public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTTask();
  /** constructor 
      @param pConf pointer to configuration descriptor
      @param pCH   the HLT component handler
   */
  AliHLTTask(AliHLTConfiguration* pConf, AliHLTComponentHandler* pCH);
  /** destructor */
  virtual ~AliHLTTask();

  /**
   * Initialize the task.
   * The task is initialized with a configuration descriptor. It needs a
   * component handler instance to create the analysis component.
   * @param pConf pointer to configuration descriptor
   * @param pCH   the HLT component handler
   */
  int Init(AliHLTConfiguration* pConf, AliHLTComponentHandler* pCH);

  /**
   * Get the name of the object.
   * This is an overridden TObject function in order to return the configuration
   * name instead of the class name. Enables use of TList standard functions.
   * @return name of the configuration
   */
  const char *GetName() const;

  /**
   * Return pointer to configuration.
   * The tasks holds internally the configuration object.
   * @return pointer to configuration
   */
  AliHLTConfiguration* GetConf() const;

  /**
   * Return pointer to component, which the task internally holds.
   * <b>Never delete this object!!!</b>
   * @return instance of the component
   */
  AliHLTComponent* GetComponent() const;

  /**
   * Find a dependency with a certain <i>name id</i>. 
   * Searches in the list of dependencies for a task.
   * @param id      the id of the <b>CONFIGURATION</b><br>
   *                <b>NOTE:</b> the id does NOT specify a COMPONENT
   * @return pointer to task
   */
  AliHLTTask* FindDependency(const char* id);

  /*
   * insert block data to the list
   * the data has to come from a task the current one depends on
   * result:
   *    -EEXIST : the block data from this task has already been inserted
   *    -ENOENT : no dependencies to the task the data is coming from
   */
  /*
   * this function is most likely depricated
  int InsertBlockData(AliHLTComponent_BlockData* pBlock, AliHLTTask* pSource);
  */

  /**
   * Add a dependency for the task.
   * The task maintains a list of other tasks it depends on.
   * @param   pDep  pointer to a task descriptor
   * @return 0 if suceeded, neg error code if failed <br>
   *    -EEXIST : the dependencie exists already
   *
   */
  int SetDependency(AliHLTTask* pDep);

  /**
   * Return number of unresolved dependencies.
   * Iterate through all the configurations the task depends on and check
   * whether a corresponding task is available in the list.
   * @return number of unresolved dependencies
   */
  int CheckDependencies();

  /**
   * Check whether the current task depends on the task pTask.
   * @param pTask pointer to Task descriptor
   * @return 1 the current task depends on pTask <br>
   *         0 no dependency <br>
   *         neg. error code if failed
   */
  int Depends(AliHLTTask* pTask);

  /**
   * Find a target with a certain id.
   * Tasks which depend on the current task are referred to be <i>targets</i>. 
   * @param id      configuration id to search for
   * @return pointer to task instance
   */
  AliHLTTask* FindTarget(const char* id);

  /**
   * Insert task into target list.
   * The target list specifies all the tasks which depend on the current task.
   * @param pDep    pointer task object
   * @return >=0 if succeeded, neg. error code if failed 
   */
  int SetTarget(AliHLTTask* pDep);

  // build a monolithic array of block data
  // @param pBlockData reference to the block data target
  // @return: array size, pointer to array in the target pTgt
  //
  /* this function is most likely depricated
  int BuildBlockDataArray(AliHLTComponent_BlockData*& pBlockData);
  */

  /**
   * Prepare the task for event processing.
   * The method initializes the Data Buffer and calls the
   * @ref AliHLTComponent::Init method of the component.<br>
   * The @ref ProcessTask methode can be called an arbitrary number of times
   * as soon as the task is in <i>running</i> mode. 
   */
  int StartRun();

  /**
   * Clean-up the task after event processing.
   * The method cleans up internal structures and calls the
   * @ref AliHLTComponent::Deinit method of the component.
   */
  int EndRun();

  /**
   * Process the task.
   * If all dependencies are resolved the tasks subscribes to the data of 
   * all source tasks, builds the block descriptor and calls the
   * @ref AliHLTComponent::ProcessEvent method of the component, after
   * processing, the data blocks are released. <br>
   * The @ref StartRun method must be called before.
   */
  int ProcessTask();

  // clear the list of source data blocks
  // the list of source data blocks has to be cleared at the beginning of 
  // a new event
  /* this function is most likely depricated
  int ClearSourceBlocks();
  */

  /**
   * Determine the number of matching data block between the component and the
   * data buffer of a consumer component. It checks which data types from the
   * list of input data types of the consumer component can be provided by data
   * blocks of the current component.
   * @param pConsumerTask   the task which subscribes to the data
   * @return number of matching data blocks
   */
  int GetNofMatchingDataBlocks(const AliHLTTask* pConsumerTask);

  /**
   * Determine the number of matching data types between the component and a
   * consumer component. It checks which data types from the list of input data
   * types of the consumer component can be provided by the current component.
   * @param pConsumerTask   the task which subscribes to the data
   * @return number of matching data types
   */
  int GetNofMatchingDataTypes(const AliHLTTask* pConsumerTask);

  /**
   * Subscribe to the data of a source task.
   * The function prepares the block descriptors for subsequent use with the
   * @ref AliHLTComponent::ProcessEvent method, the method prepares all block
   * descriptors which match the input data type of the consumer the function
   * returns the number of blocks which would be prepared in case the target
   * array is big enough.
   * @param pConsumerTask   the task which subscribes to the data
   * @param arrayBlockDesc  pointer to block descriptor to be filled
   * @param iArraySize      size of the block descriptor array
   * @return number of matching data blocks, negative error code if failed
   */
  int Subscribe(const AliHLTTask* pConsumerTask,
		AliHLTComponent_BlockData* arrayBlockDesc, int iArraySize);

  /**
   * Release a block descriptor.
   * Notification from consumer task.  
   * @param pBlockDesc      descriptor of the data segment
   * @param pConsumerTask   the task which subscribed to the data
   * @return: >0 if success, negative error code if failed
   */
  int Release(AliHLTComponent_BlockData* pBlockDesc,
	      const AliHLTTask* pConsumerTask);

  /**
   * Print the status of the task with component, dependencies and targets.
   */
  void PrintStatus();

  /**
   * Search task dependency list recursively to find a dependency.
   * @param id              id of the task to search for
   * @param pTgtList        (optional) target list to receive dependency tree
   * @return 0 if not found, >0 found in the n-th level, 
             dependency list in the target list  
   */
  int FollowDependency(const char* id, TList* pTgtList=NULL);

  /**
   * Print the tree for a certain dependency either from the task or
   * configuration list.
   * Each task posseses two "link lists": The configurations are the origin
   * of the  task list. In case of an error during the built of the task list,
   * the dependencies for the task list might be incomplete. In this case the
   * configurations can give infomation on the error reason.  
   * @param id              id of the dependency to search for
   * @param bMode           0 (default) from task dependency list, <br> 
   *                        1 from configuration list
   */
  void PrintDependencyTree(const char* id, int bMode=0);

  /**
   * Get number of source tasks.
   * @return number of source tasks
   */
  int GetNofSources() {return fListDependencies.GetSize();}

 private:
  /** the configuration descriptor */
  AliHLTConfiguration* fpConfiguration;
  /** the component described by this task */
  AliHLTComponent* fpComponent;
  /** the data buffer for the component processing */
  AliHLTDataBuffer* fpDataBuffer;
  /** the list of targets (tasks which depend upon the current one) */
  TList fListTargets;
  /** the list of sources (tasks upon which the current one depends) */ 
  TList fListDependencies;

  /**
   * block data array to be passed as argument to the 
   * @ref AliHLTComponent::ProcessEvent method. 
   * Filled through subscription to source tasks (@ref Subscribe).
   */
  AliHLTComponent_BlockData* fpBlockDataArray;
  /** size of the block data array */
  int fBlockDataArraySize;

  ClassDef(AliHLTTask, 0);
};

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
