//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTSYSTEM_H
#define ALIHLTSYSTEM_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTSystem.h
//  @author Matthias Richter
//  @date   2007
//  @brief  Global HLT module management and AliRoot integration.
//  @note   The class is used in Offline (AliRoot) context
// 

/**
 * @defgroup alihlt_system HLT integration into AliRoot
 * This section describes the HLT integration into AliRoot.
 * 
 * @section alihlt_system_intro General Remarks
 * The HLT analysis is organized in components which can also be used 
 * in off-line processing. Two different types of off-line applications 
 * can be distinguished:
 * - AliRoot simulation (AliSimulation)
 * - AliRoot reconstruction (AliReconstruction)
 *
 * In either case an HLT chain described by means of AliHLTConfiguration
 * builds the core of any HLT application. Special components exist, which 
 * emulate the behavior of the components of the HLT on-line data 
 * transportation framework. Together with the analysis components, this 
 * allows the full emulation of the behavoir of HLT analysis chain off-line.
 *
 * More details how to setup up such a chain can be found:
 * - The examples page under @ref tut_hltsystem
 *
 * @section alihlt_system_simulation AliRoot simulation
 * HLT has a special role in the normal data flow  of simulation and
 * reconstruction. Since the HLT reconstruction and analysis runs on-line
 * on the HLT farm, the raw data produced by HLT as a detector contains
 * already reconstructed events. Consequently, the HLT response has to be
 * simulated as well as the data of all other detectors. 
 *
 * Since the detector data is needed, the HLT simulation is run at the 
 * end of AliSimulation.
 *
 * As a matter of fact, HLT always reconstructs data, <em><b>HLT simulation
 * </b></em> means <em><b>HLT reconstruction embedded into AliRoot</b></em>.
 *
 * More Details an be found in the module: 
 * - @ref alihlt_aliroot_simulation
 *
 * @section alihlt_system_reconstruction AliRoot reconstruction
 *
 * Like all other ALICE detectors, HLT utilizes the AliReconstruction interface
 * to implement a plugin for the AliRoot reconstruction. The reconstructor can be
 * used to
 * - run HLT analysis chains in the AliRoot reconstruction <br>
 *   This option is mainly intended for the development and debugging cycle. HLT
 *   chains can be defined by means of AliHLTConfiguration and can be run either
 *   stand-alone or embedded into the AliReconstruction cycle.
 * - run the default analysis chains <br>
 *   HLT modules can define default analysis chains to be run during AliRoot
 *   reconstruction.
 * - handle the HLTOUT data<br>
 *   The HLT output stream contains multiple data blocks produced by the various
 *   components of the HLT chain. Each block might need different and even
 *   detector specific processing, like e.g. the processing of ESD objects or the
 *   handling of compressed data.
 *
 * More Details an be found in the module: 
 * - @ref alihlt_aliroot_reconstruction
 *
 */

#include "AliHLTLogging.h"
#include <TList.h>
#include <TString.h>

class AliHLTComponentHandler;
class AliHLTConfiguration;
class AliHLTConfigurationHandler;
class AliHLTTask;
class AliHLTOUT;
class AliHLTOUTTask;
class AliHLTControlTask;
class AliRunLoader;
class AliRawReader;
class AliESDEvent;
class TObjArray;
class TStopwatch;

/**
 * @class AliHLTSystem
 * Main class for the HLT integration into AliRoot.
 * @note This class is only used for the @ref alihlt_system.
 *
 * @section alihltsystem_overview Main functionality
 * The class implements the main binding functionality to run HLT
 * embedded in AliRoot:
 * - creates handlers for HLT components (AliHLTComponentHandler)
 * and configurations (AliHLTConfigurationHandler).
 * - implements the creation of a task list from configurations
 * - runs the tasks of the task list
 * - implements the interface to AliHLTReconstructor and AliHLTSimulation
 * - manages the registered module agents (AliHLTModuleAgent) 
 *
 * @section alihltsystem_functions Main functions
 * The following list gives the most important functions, see documentation
 * of functions for details.
 * - ScanOptions(const char*) initialize the system from arguments
 * - BuildTaskList(const char*) build task list for a reconstruction chain
 * - Run() processing loop for task list
 * - Reconstruct() interface to AliHLTReconstructor
 * - ProcessHLTOUT() process an HLTOUT colection
 *
 * @section alihltsystem_initialization Initialization
 * AliHLTSystem is initialized with a list of component libraries which determine
 * the actual functionality. See ScanOptions() for description of options.
 *
 * @section alihltsystem_usage Usage
 * The class AliHLTPluginBase handles the global instance of AliHLTSystem,
 * the instance is created and fetched like 
 * <pre>
 * // setup the HLT system
 * AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
 * </pre>
 *
 * @ingroup alihlt_system
 */
class AliHLTSystem : public AliHLTLogging {
 public:
  /** default constructor */
  AliHLTSystem(AliHLTComponentLogSeverity loglevel=kHLTLogDefault, const char* name="");
  /** destructor */
  virtual ~AliHLTSystem();

  /**
   * Build a task list 
   * This method is used to build the list of tasks from the configuration
   * id of a single 'master' configuration.
   *
   * Configuration entries might depend upon other configurations. For each
   * configuration which has not yet been converted into an AliHLTTask, the
   * method will be called iteratively. Finally, after building all tasks
   * providing the input for the current one, the task is inserted to the
   * list of tasks with the InsertTask() method.
   *
   * The function can be called multiple times in order to add more than one
   * chain to the system.
   *
   * @param pConf    configuration name/id
   */
  int BuildTaskList(const char* pConf);

  /**
   * Clean the list of tasks and delete all the task objects.
   */
  int CleanTaskList();

  /**
   * Insert a task to the task list.
   * The method first checks whether all dependencies are resolved (i.e. exist 
   * already in the task list). During this iteration the cross links between
   * the tasks are set as well. If all dependencies are resolved, the task is
   * added at the end of the list.
   * @param pTask    pointer to task to add
   */
  int InsertTask(AliHLTTask* pTask);

  /**
   * Add HLTOUT task to the end of the task list.
   * If one of the specified chains has output, an AliHLTOUTTask is
   * added which controls the output. All other chains are removed from the
   * AliHLTOUTTask input.
   * @return 0 if no task has been added, 1 if task has been added
   */
  int AddHLTOUTTask(const char* chains);

  AliHLTOUTTask* GetHLTOUTTask() const {return fpHLTOUTTask;}
  AliHLTComponentHandler* GetComponentHandler() const {return fpComponentHandler;}
  AliHLTConfigurationHandler* GetConfigurationHandler() const {return fpConfigurationHandler;}

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
   * Run one or more events.
   * Core of the processing loop. The method expects the task list to be
   * already created by a previous call to BuildTaskList(const char*)
   *
   * All tasks of the list will be subsequently processed for each event.
   * The system can remain started if the \em bStop parameter is 0. In that
   * case the system just waits for the next event. A specific call with
   * nofEvents==0 is needed to execute the stop sequence.
   * @param iNofEvents number of events
   * @param bStop      stop the chain after processing
   * @return number of reconstructed events, neg error code if failed
   */
  int Run(Int_t iNofEvents, int bStop, AliHLTUInt64_t trgMask);

  /**
   * Run the tasklist
   * Somehow the 64bit variable/default value did not work out on mac.
   * Re-introducing the original function, and forwarding it.
   */
  int Run(Int_t iNofEvents=1, int bStop=1)
  {return Run(iNofEvents, bStop, 0);}

  /**
   * Init all tasks from the list.
   * The @ref AliHLTTask::Init method is called for each task, the components
   * will be created.
   * @return neg error code if failed
   */
  int InitTasks();

  /**
   * Init the stopwatches for all tasks.
   * @param pStopwatches    object array of stopwatches, for types
   *                        @see AliHLTComponent::AliHLTStopwatchType
   * @return neg error code if failed
   */
  int InitBenchmarking(TObjArray* pStopwatches);

  /**
   * Stop the stopwatches for all tasks.
   * @param pStopwatches    object array of stopwatches, for types
   *                        @see AliHLTComponent::AliHLTStopwatchType
   * @return neg error code if failed
   */
  int PauseBenchmarking(TObjArray* pStopwatches) const;

  /**
   * Continue the stopwatches for all tasks.
   * @param pStopwatches    object array of stopwatches, for types
   *                        @see AliHLTComponent::AliHLTStopwatchType
   * @return neg error code if failed
   */
  int ResumeBenchmarking(TObjArray* pStopwatches) const;

  /**
   * Print benchmarking summary.
   * Optionak: clean up stop watches.
   * @param pStopwatches    object array of stopwatches
   * @param bClean          delete stop watches if 1
   * @return neg error code if failed
   */
  int PrintBenchmarking(TObjArray* pStopwatches, int bClean=0) const;

  /**
   * Start task list.
   * The @ref AliHLTTask::StartRun method is called for each task, the
   * components will be prepared for event processing.
   * @return neg error code if failed
   */
  int StartTasks();

  /**
   * Process task list.
   * The @ref AliHLTTask::ProcessTask method is called for each task.
   * @return neg error code if failed
   */
  int ProcessTasks(Int_t eventNo, AliHLTUInt64_t trgMask=0);

  /**
   * Stop task list.
   * The @ref AliHLTTask::EndRun method is called for each task, the components
   * will be cleaned after event processing.
   * @return neg error code if failed
   */
  int StopTasks();

  /**
   * Send a control event trough the chain.
   * All data sources in the chain are switched to publish a control event like
   * SOR or EOR. The event is propagated in the same way as a normal event.
   * @param dt       type of the event
   */
  int SendControlEvent(AliHLTComponentDataType dt);

  /**
   * De-init all tasks from the list.
   * The @ref AliHLTTask::Deinit method is called for each task, the components
   * will be deleted.
   * @return neg error code if failed
   */
  int DeinitTasks();

  /**
   * Cleanup all internal objects from HLTOUT processing.
   */
  int CleanHLTOUT();

  /**
   * The memory allocation function for components.
   * This function is part of the running environment of the components, 
   * see AliHLTAnalysisEnvironment
   */
  static void* AllocMemory( void* param, unsigned long size );

  /**
   * The allocation function for component EventDoneData.
   * This function is part of the running environment of the components,
   * see AliHLTAnalysisEnvironment
   */
  static int AllocEventDoneData( void* param, AliHLTEventID_t eventID, unsigned long size, AliHLTComponentEventDoneData** edd );

  /**
   * AliRoot embedded reconstruction.
   * Main entry point to execute the HLT reconstruction from AliRoot. Called
   * either by the AliHLTReconstructor plugin during AliRoot reconstruction
   * of raw data, or AliHLTSimulation during simulation of data.
   *
   * The two cases are distinguished by the availablility of the AliRunLoader
   * and raw reader.
   * - AliRoot simulation: AliRunLoader is available and is propagated to the
   *   module agents (AliHLTModuleAgent) to generate the corresponding
   *   configurations and chains, and to the AliHLTOfflineSource components.
   *   raw reader might be available depending on whether raw data was
   *   simulated or not.
   * - AliRoot reconstruction: raw reader is available and is propagated to
   *   the agents and AliHLTOfflineSource components. AliRunLoader is always
   *   NULL.
   *
   * The system remains started after the processing and just waits for the
   * next event. A specific call with nofEvents==0 is needed to execute the
   * stop sequence.
   *
   * The 'runLoader' and 'rawReader' parameters are propagated to all active
   * AliHLTOfflineDataSource's and the HLT chain is processed for the given
   * number of events. If the rawReader is NULL, reconstruction is done on
   * simulated data, from real data if a RawReader is specified.
   *
   * After setting up the system with the reconstruction parameters the
   * Run() method is called and carries out the processing of the chain.
   *
   * @param nofEvents     number of events
   * @param runLoader     the AliRoot RunLoader
   * @param rawReader     the AliRoot RawReader
   * @return number of reconstructed events, neg. error code if failed 
   */
  int Reconstruct(int nofEvents, AliRunLoader* runLoader, 
		  AliRawReader* rawReader=NULL);

  /**
   * Fill ESD for one event.
   * To be called by the AliHLTReconstructor plugin during the event loop
   * and FillESD method of the AliRoot reconstruction.
   *
   * The method is most likely deprecated as the scheme has been slightly
   * changed. The ESD is filled by the HLTOUT handlers implemented by the
   * HLT modules rather than by components within the reconstruction chain.
   * Still, HLT reconstruction chains can be run during the AliRoot
   * reconstruction, data produced by such chains is automatically added
   * to the HLTOUT stream in order to be equivalent to the online HLT.
   * The HLTOUT is processed during AliReconstruction at the end of the
   * HLT event processing, literally during the FillESD method of the AliRoot
   * reconstruction interface. The HLT module must implement HLTOUT handlers
   * and provide those through the module agent.
   *
   * @param eventNo       current event no (Note: this event number is just a
   *                      processing counter and is not related to the nature/
   *                      origin of the event
   * @param runLoader     the AliRoot runloader
   * @param esd           an AliESDEvent instance
   * @return neg. error code if failed 
   */
  int FillESD(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd);

  /**
   * Process the HLTOUT data.
   * The provided instance of AliHLTOUT provides the access to the data.
   * AliHLTSystem queries all registered module agents (AliHLTModuleAgent)
   * for the ability to treat a specific data block. As the last step the
   * ESD object is filled. Please note that the provided ESD is the hltEsd
   * in case of AliReconstructor (switched in AliReconstruction).
   */
  int ProcessHLTOUT(AliHLTOUT* pHLTOUT, AliESDEvent* esd);

  /**
   * Process all kChain-type data blocks of the HLTOUT data.
   * Handlers of type AliHLTModuleAgent::kChain are executed as the first
   * step of the processing. The function is invoked from ProcessHLTOUT(),
   * eventual output of the chain is added to the HLTOUT collection.
   */
  int ProcessHLTOUTkChain(AliHLTOUT* pHLTOUT);

  /**
   * Load component libraries.
   * @param libs          string of blank separated library names
   * @return neg. error code if failed 
   */
  int LoadComponentLibraries(const char* libs);

  /**
   * Find a symbol in a dynamically loaded library.
   * @param library      library
   * @param symbol       the symbol to find
   * @return void pointer to function
   */
  AliHLTfctVoid FindDynamicSymbol(const char* library, const char* symbol);

  /**
   * Prepare the HLT system for running.
   * - module agents are requested to register configurations
   * - task lists are built from the reconstruction chains of the modules
   *
   * @param rawReader    instance of the raw reader or NULL
   * @param runloader    optional instance of the AliRunLoader
   * @return neg. error code if failed <br>
   *         -EBUSY      system is in kRunning state <br>
   */
  int Configure(AliRawReader* rawReader, AliRunLoader* runloader=NULL);

  /**
   * Old method kept for backward compatibilty.
   *
   * @param runloader    optional instance of the AliRunLoader
   * @return neg. error code if failed <br>
   *         -EBUSY      system is in kRunning state <br>
   */
  int Configure(AliRunLoader* runloader=NULL);

  /**
   * Scan options and load component libraries.
   * The options consist of blank separated tokens. Libraries can be just
   * specified by their name, and can be excluded by adding a '!'-mark in
   * front. <br>
   * Further options
   * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
   * \li loglevel=<i>level</i> <br>
   *     logging level for this processing
   * \li frameworklog=<i>level</i> <br>
   *     logging level for framework classes
   * \li alilog=off
   *     disable redirection of log messages to AliLog class
   * \li config=<i>macro</i>
   *     configuration macro
   * \li chains=<i>configuration</i>
   *     comma separated list of configurations to be run during local
   *     reconstruction
   * \li libmode=<i>static,dynamic(default)</i>
   *     libraries are persistent if loaded in mode <i>static</i>, i.e. they
   *     can't be unloaded
   */
  int ScanOptions(const char* options);

  /**
   * Reset the HLT system.
   * Reset is not possible while the system is in running state.
   * @param bForce       force the reset
   * @return neg. error code if failed <br>
   *         -EBUSY      system is in kRunning state <br>
   */
  int Reset(int bForce=0);

  /**
   * Load the configurations specified by the module agents.
   * The runLoader is passed to the agent and allows configuration
   * selection.
   * - AliSimulation: runloader valid, raw reader might be valid
   * - AliReconstruction: runloader always NULL, raw reader valid
   *
   * @param rawReader    instance of the raw reader or NULL
   * @param runloader    optional instance of the AliRunLoader
   * @return neg. error code if failed 
   */
  int LoadConfigurations(AliRawReader* rawReader, AliRunLoader* runloader=NULL);

  /**
   * Get the reconstruction chains from all agents and build the task lists.
   * AliHLTModuleAgent implementations can define reconstruction chains
   * depending on the availibility of AliRunLoader and AliRawReader parameters.
   *
   * @param rawReader    instance of the raw reader or NULL
   * @param runloader    optional instance of the AliRunLoader
   * @return neg. error code if failed 
   */
  int BuildTaskListsFromReconstructionChains(AliRawReader* rawReader, 
					     AliRunLoader* runloader=NULL);

  enum AliHLTSystemState {
    kUninitialized       = 0x0,
    kLibrariesLoaded     = 0x1,
    kConfigurationLoaded = 0x2,
    kTaskListCreated     = 0x4,
    kReady               = 0x7,
    kStarted             = 0x8,
    kRunning             = 0x10,
    kError               = 0x1000
  };

  /**
   * Check status of the system.
   * @param flag          AliHLTSystemState value to check for
   * @return 1 if set, 0 if not
   */
  int CheckStatus(int flag);

  /**
   * Get the current status.
   * @return status flags of @ref AliHLTSystemState
   */
  int GetStatusFlags();

  /**
   * Set logging level for framework classes.
   * This sets the local logging level of this instance and all subsequent
   * framework classes to \em level.
   * @param level     local logging level for the framework classes
   */
  void SetFrameworkLog(AliHLTComponentLogSeverity level);

  /**
   * Customized logging function.
   * The name of the system and pointer is added at the beginning of each
   * message if name was set.
   */
  int LoggingVarargs(AliHLTComponentLogSeverity severity, 
		     const char* originClass, const char* originFunc,
		     const char* file, int line, ... ) const;

 protected:
 
 private:
  /** copy constructor prohibited */
  AliHLTSystem(const AliHLTSystem&);
  /** assignment operator prohibited */
  AliHLTSystem& operator=(const AliHLTSystem&);

  /**
   * Build task list from a configuration object. 
   * This method implements the configuration parsing and transformation into
   * a list of AliHLTTask objects. It has been made private in April 2007. 
   * Use BuildTaskList(const char*) instead.
   *
   * @param pConf    pointer to configuration to build the task list from
   */
  int BuildTaskList(AliHLTConfiguration* pConf);

  /**
   * Set status flags.
   */
  int SetStatusFlags(int flags);

  /**
   * clear status flags.
   */
  int ClearStatusFlags(int flags);

  /// Pointer to an instance of @ref AliHLTComponentHandler.
  AliHLTComponentHandler* fpComponentHandler;                      //! transient

  /// Pointer to an instance of @ref AliHLTConfigurationHandler.
  AliHLTConfigurationHandler* fpConfigurationHandler;              //! transient

  /** list of tasks */
  TList fTaskList;                                                 // see above

  /** the number of instances of AliHLTSystem */
  static int fgNofInstances;                                       // see above

  /** state of the object */
  int fState;                                                      // see above

  /** chains to be run during reconstruction */
  TString fChains;                                                 //!transient

  /** array of stopwatches */
  TObjArray* fStopwatches;                                         //!transient

  /** number of events processed in total */                        
  int fEventCount;                                                 //!transient

  /** number of events processed successfully */
  int fGoodEvents;                                                 //!transient

  /** array of default libraries */
  static const char* fgkHLTDefaultLibs[];                          //!transient

  /** active kChain handlers (AliHLTOUT::AliHLTOUTHandlerListEntryVector*) */
  void* fpChainHandlers;                                           //!transient

  /** active kEsd handlers (AliHLTOUT::AliHLTOUTHandlerListEntryVector*) */
  void* fpEsdHandlers;                                             //!transient

  /** active kProprietary handlers (AliHLTOUT::AliHLTOUTHandlerListEntryVector*) */
  void* fpProprietaryHandlers;                                     //!transient

  /** active HLTOUT task for the reconstruction */
  AliHLTOUTTask* fpHLTOUTTask;                                     //!transient

  /** special task to publish the control events */
  AliHLTControlTask* fpControlTask;                                //!transient

  /** name of this system instance */
  TString fName;                                                   //!transient

  /// ECS parameter string
  TString fECSParams;                                              //!transient

  ClassDef(AliHLTSystem, 13);
};

#endif
