// @(#) $Id$

#ifndef ALIHLTCONFIGURATION_H
#define ALIHLTCONFIGURATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTConfiguration
   base class for HLT configurations
 */

#include <errno.h>
#include <TObject.h>
#include <TList.h>
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"

class AliHLTConfigurationHandler;
/*****************************************************************************************************
 *
 * AliHLTConfiguration
 *
 * this class describes a certain configuration af an HLT processing step by the following parameters:
 *  - a unique id string/name
 *  - the id of the component
 *  - the ids of the configurations it requires input from
 *  - the arguments, which are passed to the component when it is initialized
 *
 * The setup of a configuration requires simply the creation of a global object of class AliHLTConfiguration.
 * The Configuration is automatically registered in the list of available configurations. The list is used
 * by the handler to resolve the dependencies upon other configurations. Hierarchies can be built up in
 * an easy way.
 *
 * A configuration is interpreted by the Configuration Handler and transformed into a Task List.
 */
class AliHLTConfiguration : public TObject, public AliHLTLogging {
 public:
  AliHLTConfiguration();
  AliHLTConfiguration(const char* id, const char* component, const char* sources, const char* arguments);
  virtual ~AliHLTConfiguration();

  /****************************************************************************************************
   * global initialization
   */
  static int GlobalInit(AliHLTConfigurationHandler* pHandler);

  static int GlobalDeinit();

  /****************************************************************************************************
   * properties
   */

  // configuration id, a unique name
  // overridden TObject function in order to return the configuration name instead of the class name
  // enables use of TList standard functions
  const char *GetName() const;

  // id of the component
  const char* GetComponentID() {return fComponent;}

  // print status info
  void PrintStatus();

  // get a certain source
  AliHLTConfiguration* GetSource(const char* id);

  // try to find a dependency recursively in the list of sources
  // parameter:
  //   id - the source to search for
  //   pTgtList - (optional) target list to receive the dependency tree
  // result:
  //   0 if not found
  //   n found in the n-th level
  //   dependency list in the target list  
  int FollowDependency(const char* id, TList* pTgtList=NULL);

  // get the number of resolved sources
  int GetNofSources() {return fListSources.size();}

  // 1 if all sources are resolved
  // try to resolve if bAuto==1 
  int SourcesResolved(int bAuto=0);

  // start iteration and get the first source
  AliHLTConfiguration* GetFirstSource();

  // continue iteration and get the next source
  AliHLTConfiguration* GetNextSource();

  // invalidate a dependency and mark the configuration to be re-evaluted 
  int InvalidateSource(AliHLTConfiguration* pConf);

  // mark the configuration to be re-evaluted 
  int InvalidateSources() {fNofSources=-1; return 0;}

  // return the arguments array
  int GetArguments(int* pArgc, const char*** pArgv);

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

  const char* fID;                  // id of this configuration
  const char* fComponent;           // component id of this configuration

  const char* fStringSources;                             // the 'sources' string as passed to the constructor
  int fNofSources;
  vector<AliHLTConfiguration*> fListSources;              // list of sources
  vector<AliHLTConfiguration*>::iterator fListSrcElement; // iterator for the above list

  const char* fArguments;           // the arguments string as passed to the constructor
  int fArgc;                        // number of arguments
  char** fArgv;                     // argument array

  static AliHLTConfigurationHandler* fConfigurationHandler;

  ClassDef(AliHLTConfiguration, 0);
};

struct AliHLTComponent_BlockData;
class AliHLTComponent;
class AliHLTComponentHandler;

/*****************************************************************************************************
 *
 * AliHLTTask
 *
 * 
 */
class AliHLTTask : public TObject, public AliHLTLogging {
 public:
  AliHLTTask();
  AliHLTTask(AliHLTConfiguration* fConf, AliHLTComponentHandler* pCH);
  virtual ~AliHLTTask();

  /* initialize the task
   */
  int Init(AliHLTConfiguration* fConf, AliHLTComponentHandler* pCH);

  // overridden TObject function in order to return the configuration name instead of the class name
  // enables use of TList standard functions
  const char *GetName() const;

  /* return pointer to configuration
   */
  AliHLTConfiguration* GetConf();

  /* return pointer to configuration
   */
  AliHLTComponent* GetComponent();

  /* find a dependency with name/id 
   */
  AliHLTTask* FindDependency(const char* id);

  /* insert block data to the list
   * the data has to come from a task the current one depend on
   * result:
   *    -EEXIST : the block data from this task has already been inserted
   *    -ENOENT : no dependencies to the task the data is coming from
   */
  int InsertBlockData(AliHLTComponent_BlockData* pBlock, AliHLTTask* pSource);

  /* add a dependency for the task
   */
  int SetDependency(AliHLTTask* pDep);

  /* return number of unresolved dependencies
   */
  int CheckDependencies();

  /* 1 if the current task depends upon pTask
   * 0 if no dependency
   */
  int Depends(AliHLTTask* pTask);

  /* find a target
   */
  AliHLTTask* FindTarget(const char* id);

  /* insert task into target list
   */
  int SetTarget(AliHLTTask* pDep);

  // build a monolithic array of block data
  // result: array size, pointer to array in the target pTgt
  //
  int BuildBlockDataArray(AliHLTComponent_BlockData*& pTgt);

  /* process the task if all dependencies are resolved
   * reset the source block data list
   */
  //int ProcessTask(...);

  // clear the list of source data blocks
  // the list of source data blocks has to be cleared at the beginning of 
  // a new event
  int ClearSourceBlocks();

  // print the status of the task with component, dependencies and targets
  void PrintStatus();

  // search task dependency list recursively to find a dependency 
  // parameter:
  //   id - the task to search for
  //   pTgtList - (optional) target list to receive the dependency tree
  // result:
  //   0 if not found
  //   n found in the n-th level
  //   dependency list in the target list  
  int FollowDependency(const char* id, TList* pTgtList=NULL);

  // print the tree for a certain dependency either from the task or configuration list
  // each task posseses two "link lists": The configurations are the the origin of the 
  // task list. In case of an error during the built of the task list, the dependencies 
  // for the task list might be incomplete. In this case the configurations can give 
  // infomation on the error reason  
  // parameter:
  //   id - the dependency to search for
  //   bFromConfiguration (default=0) - if 0 the task dependency list is used, if one the configuration list is used
  void PrintDependencyTree(const char* id, int bFromConfiguration=0);

 private:
  AliHLTConfiguration* fpConfiguration;
  AliHLTComponent* fpComponent;
  vector<AliHLTComponent_BlockData*> fSources;
  TList fListTargets;
  TList fListDependencies;

  AliHLTComponent_BlockData* fpBlockDataArray;

  ClassDef(AliHLTTask, 0);
};

class AliHLTConfigurationHandler : public AliHLTLogging {
 public:
  AliHLTConfigurationHandler();
  //AliHLTConfigurationHandler(AliHLTConfiguration* pConf);
  virtual ~AliHLTConfigurationHandler();

  /****************************************************************************************************
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
