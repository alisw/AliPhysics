// $Id$

#ifndef ALIHLTTASK_H
#define ALIHLTTASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTask.h
    @author Matthias Richter
    @date   
    @brief  base class for HLT tasks
*/

#include <cerrno>
#include <TObject.h>
#include <TList.h>
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTDataBuffer.h"

struct AliHLTComponentBlockData;
class AliHLTComponent;
class AliHLTComponentHandler;

/******************************************************************************/

/**
 * @class AliHLTTask
 * A task collects all the information which is necessary to process a certain
 * step in the HLT data processing chain.
 * - the instance of the component
 *   the task object creates and deletes the component object
 * - the data buffer which receives the result of the component and provides
 *   the data to other tasks/components
 * - a list of all dependencies
 * - a list of consumers
 * - the task object holds an external pointer to the configuration object; 
 *   \b Note: the configuration object must exist through the existence of the
 *   task object!!!
 *  
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTTask : public TObject, public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTTask();
  /** constructor 
      @param pConf pointer to configuration descriptor
   */
  AliHLTTask(AliHLTConfiguration* pConf);
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTask(const AliHLTTask&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTask& operator=(const AliHLTTask&);
  /** destructor */
  virtual ~AliHLTTask();

  /**
   * Initialize the task.
   * The task is initialized with a configuration descriptor. It needs a
   * component handler instance to create the analysis component.
   * @param pConf pointer to configuration descriptor, can be NULL if it
   *              was already provided to the constructor
   * @param pCH   the HLT component handler
   */
  int Init(AliHLTConfiguration* pConf, AliHLTComponentHandler* pCH);

  /**
   * De-Initialize the task.
   * Final cleanup after the run. The analysis component is deleted.
   */
  int Deinit();

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
  int InsertBlockData(AliHLTComponentBlockData* pBlock, AliHLTTask* pSource);
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
  int BuildBlockDataArray(AliHLTComponentBlockData*& pBlockData);
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
		AliHLTComponentBlockData* arrayBlockDesc, int iArraySize);

  /**
   * Release a block descriptor.
   * Notification from consumer task.  
   * @param pBlockDesc      descriptor of the data segment
   * @param pConsumerTask   the task which subscribed to the data
   * @return: >0 if success, negative error code if failed
   */
  int Release(AliHLTComponentBlockData* pBlockDesc,
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
  /** the configuration descriptor (external pointer) */
  AliHLTConfiguration* fpConfiguration;
  /** the component described by this task (created and deleted by the task object) */
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
  AliHLTComponentBlockData* fpBlockDataArray;
  /** size of the block data array */
  int fBlockDataArraySize;

  ClassDef(AliHLTTask, 1);
};

#endif
