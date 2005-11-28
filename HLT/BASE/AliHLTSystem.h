// @(#) $Id$

#ifndef ALIHLTSYSTEM_H
#define ALIHLTSYSTEM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTComponentHandler
   global HLT module management 
 */


#include "AliL3RootTypes.h"
#include "AliHLTLogging.h"
#include <TList.h>

class AliHLTComponentHandler;
class AliHLTConfiguration;
class AliHLTConfigurationHandler;
class AliHLTTask;

class AliHLTSystem : public AliHLTLogging {
 public:
  AliHLTSystem();
  virtual ~AliHLTSystem();

  /* this will change later
   */
  AliHLTComponentHandler* fpComponentHandler;
  AliHLTConfigurationHandler* fpConfigurationHandler;

  /* add a configuration to the end of the list
   */
  int AddConfiguration(AliHLTConfiguration* pConf);

  /* add a configuration to the list after the specified configuration
   */
  int InsertConfiguration(AliHLTConfiguration* pConf, AliHLTConfiguration* pPrec);

  /* remove a configuration from the list
   */
  int DeleteConfiguration(AliHLTConfiguration* pConf);

  /* build a task list from a configuration object
   * This method is used to build the tasks from the 'master' configuration
   * objects which are added to the HLT system handler. This is an iterative
   * process since the task might depend upon other configurations. For each
   * configuration object which has not yet been converted into a task, the
   * method will be called iteratively. Finally, after building all tasks which
   * the current one depends on have been created, the task is inserted to the
   * list of tasks with the InsertTask method.
   */
  int BuildTaskList(AliHLTConfiguration* pConf);

  /* clean the list of tasks and delete all the task objects
   */
  int CleanTaskList();

  /* insert a task to the task list
   * the method first checks whether all dependencies are resolved (i.e. exist 
   * already in the task list). During this iteration the cross links between the 
   * tasks are set as well. If all dependencies are resolved, the task is added
   * at the end of the list.
   */
  int InsertTask(AliHLTTask* pTask);

  /* find a task with an id
   * NOTE: 'id' denotes a CONFIGURATION, not a COMPONENT
   */
  AliHLTTask* FindTask(const char* id);

  /* print the task list
   */
  void PrintTaskList();

  /* run the task list
   */
  int Run();

 protected:
  int ProcessTask();
  int StartEvent();
  int ProcessEvent();
  int StopEvent();
 
 private:
  TList fConfList;
  int fbListChanged;

  TList fTaskList;

 private:
  ClassDef(AliHLTSystem, 0);
};
#endif

