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
class AliHLTTask;

class AliHLTSystem : public AliHLTLogging {
 public:
  AliHLTSystem();
  virtual ~AliHLTSystem();

  AliHLTComponentHandler* fpComponentHandler;

  /* add a configuration to the end of the list
   */
  int AddConfiguration(AliHLTConfiguration* pConf);

  /* add a configuration to the list after the specified configuration
   */
  int InsertConfiguration(AliHLTConfiguration* pConf, AliHLTConfiguration* pPrec);

  /* remove a configuration from the list
   */
  int DeleteConfiguration(AliHLTConfiguration* pConf);

  /* build a task list from the configuration list
   */
  int BuildTaskList(AliHLTConfiguration* pConf);

  int CleanTaskList();

  int InsertTask(AliHLTTask* pTask);

  AliHLTTask* FindTask(const char* id);

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

