// @(#) $Id$

#ifndef ALIHLTSYSTEM_H
#define ALIHLTSYSTEM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTSystem.h
    @author Matthias Richter
    @date   
    @brief  Global HLT module management and AliRoot integration.
    @note   The class is used in Offline (AliRoot) context
*/

/**
 * @defgroup alihlt_system HLT integration into AliRoot
 * This section describes the HLT integration into AliRoot.
 */

#include "AliHLTLogging.h"
#include <TList.h>

class AliHLTComponentHandler;
class AliHLTConfiguration;
class AliHLTConfigurationHandler;
class AliHLTTask;

/**
 * @class AliHLTSystem
 * Main class for the HLT integration into AliRoot.
 * The class handles a list of configurations. Configurations are translated
 * into task lists which can be executed. 
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTSystem : public AliHLTLogging {
 public:
  /** default constructor */
  AliHLTSystem();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTSystem(const AliHLTSystem&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTSystem& operator=(const AliHLTSystem&);
  /** destructor */
  virtual ~AliHLTSystem();

  /**
   * Pointer to an instance of @ref AliHLTComponentHandler.
   */
  AliHLTComponentHandler* fpComponentHandler;

  /**
   * Pointer to an instance of @ref AliHLTConfigurationHandler.
   */
  AliHLTConfigurationHandler* fpConfigurationHandler;

  /**
   * Add a configuration to the end of the list.
   * @param pConf    pointer to configuration to add
   */
  int AddConfiguration(AliHLTConfiguration* pConf);

  /**
   * Insert a configuration to the end of the list after the specified configuration.
   * @param pConf    pointer to configuration to insert
   * @param pPrec    pointer to configuration to insert the new one after
   */
  int InsertConfiguration(AliHLTConfiguration* pConf, AliHLTConfiguration* pPrec);

  /**
   * Remove a configuration from the list.
   * @param pConf    pointer to configuration to delete
   */
  int DeleteConfiguration(AliHLTConfiguration* pConf);

  /**
   * Build a task list from a configuration object.
   * This method is used to build the tasks from the 'master' configuration
   * objects which are added to the HLT system handler. This is an iterative
   * process since the task might depend upon other configurations. For each
   * configuration object which has not yet been converted into a task, the
   * method will be called iteratively. Finally, after building all tasks which
   * the current one depends on have been created, the task is inserted to the
   * list of tasks with the InsertTask method.
   * @param pConf    pointer to configuration to build the task list from
   */
  int BuildTaskList(AliHLTConfiguration* pConf);

  /**
   * Clean the list of tasks and delete all the task objects.
   */
  int CleanTaskList();

  /**
   * Insert a task to the task list.
   * The method first checks whether all dependencies are resolved (i.e. exist 
   * already in the task list). During this iteration the cross links between the 
   * tasks are set as well. If all dependencies are resolved, the task is added
   * at the end of the list.
   * @param pTask    pointer to task to add
   */
  int InsertTask(AliHLTTask* pTask);

  /**
   * Find a task with an id.
   * @param id       CONFIGURATION id (not a COMPONENT id!)
   */
  AliHLTTask* FindTask(const char* id);

  /**
   * Print the task list.
   */
  void PrintTaskList();

  /**
   * Print info on an AliHLTComponentDataType structure
   */
  void PrintComponentDataTypeInfo(const AliHLTComponentDataType& dt);

  /**
   * Run the task list.
   * All tasks of the list will be subsequently processed for each event.
   * @param iNofEvents number of events
   * @return neg error code if failed
   */
  int Run(Int_t iNofEvents=1);

  /**
   * Start task list.
   * The @ref AliHLTTask::StartRun method is called for each task, the components
   * will be prepared for event processing.
   * @return neg error code if failed
   */
  int StartTasks();

  /**
   * Process task list.
   * The @ref AliHLTTask::ProcessTask method is called for each task.
   * @return neg error code if failed
   */
  int ProcessTasks(Int_t eventNo);

  /**
   * Stop task list.
   * The @ref AliHLTTask::EndRun method is called for each task, the components
   * will be cleaned after event processing.
   * @return neg error code if failed
   */
  int StopTasks();

 protected:
  int ProcessTask();
  int StartEvent();
  int ProcessEvent();
  int StopEvent();
 
 private:
/*   TList fConfList; */
/*   int fbListChanged; */

  TList fTaskList;

 private:
  ClassDef(AliHLTSystem, 1);
};
#endif

