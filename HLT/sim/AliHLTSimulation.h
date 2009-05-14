//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTSIMULATION_H
#define ALIHLTSIMULATION_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTSimulation.h
    @author Matthias Richter
    @date   
    @brief  Binding class for HLT simulation in AliRoot

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */
/**
 * @defgroup alihlt_aliroot_simulation HLT simulation in AliRoot
 * This section describes the the simulation of the HLT in AliRoot.
 *
 * @section alihlt_aliroot_simulation_intro General Remarks
 * HLT has a special role in the normal data flow  of simulation and
 * reconstruction. Since the HLT reconstruction and analysis runs on-line
 * on the HLT farm, the raw data produced by HLT as a detector contains
 * already reconstructed events. Consequently, the HLT response has to be
 * simulated as well as the data of all other detectors. Since the detector
 * data is needed, the HLT simulation is run at the end of AliSimulation.
 * As a matter of fact, HLT always reconstructs data, <em><b>HLT simulation
 * </b></em> means <em><b>HLT reconstruction embedded into AliRoot</b></em>.
 *
 * @section alihlt_aliroot_simulation_steering Steering
 * The AliHLTSimulation class is the steering class called from AliSimulation.
 * An instance of AliHLTSystem is used to run the chains defined by the
 * available libraries or a AliHLTConfiguration configuration macro.
 *
 * The libraries to be loaded can be specified as an option to AliSimulation.
 * <pre>
 * AliSimulation sim;
 * sim.SetRunHLT("libAliHLTSample.so");
 * </pre>
 * @see AliHLTSimulation for further details
 *
 * @section alihlt_aliroot_simulation_running Running
 * The actual chains to be run depend on the HLT library modules which 
 * are loaded to the system. There is a default collection of libraries 
 * defined in AliHLTSystem::fgkHLTDefaultLibs. The default libraries are 
 * loaded if nothing else is specified. The libraries implement \em agents 
 * (AliHLTModuleAgent childs) describing the properties of a module.
 *
 * @section alihlt_aliroot_simulation_examples Examples
 * - @ref tut_simulation 
 *
 * @ingroup alihlt_system
 */

#include "TObject.h"
#include "TString.h"
class AliRunLoader;
class AliHLTPluginBase;
class AliRawReader;

/**
 * @class AliHLTSimulation
 * Plugin for HLT reconstruction embedded into <tt>AliSimulation</tt>.
 *
 * The libraries to be loaded can be specified as an option to AliSimulation.
 * <pre>
 * AliSimulation sim;
 * sim.SetRunHLT("libAliHLTSample.so");
 * </pre>
 * will only load <tt>libAliHLTSample.so</tt>
 *
 * Other available options:
 * \li loglevel=<i>level</i>                                            <br>
 *     logging level for this processing, default level is 0x79 filtering
 *     out everything below level 'warning'. 0x7c allows info messages as
 *     well, 0x3f is the highest loglevel.
 * \li alilog=off                                                       <br>
 *     disable redirection of log messages to AliLog class
 * \li config=<i>macro</i>                                              <br>
 *     configuration macro: normal ROOT macro defining HLT component
 *     configurations by means of AliHLTConfiguration.
 * \li chains=<i>configuration</i>                                      <br>
 *     comma separated list of configurations to be run during simulation
 *
 *  @ingroup alihlt_aliroot_simulation
 */
class AliHLTSimulation : public TObject {
 public:
  /** create an instance of the class */
  static AliHLTSimulation* CreateInstance();

  /** delete an instance */
  static int DeleteInstance(AliHLTSimulation* pSim);

  /** init simulation */
  int Init(AliRunLoader* pRunLoader, const char* options);

  /** run simulation with an instance of the run loader */
  int Run(AliRunLoader* pRunLoader);

 private:
  /** standard constructor */
  AliHLTSimulation();
  /** copy constructor prohibited */
  AliHLTSimulation(const AliHLTSimulation&);
  /** assignment operator prohibited */
  AliHLTSimulation& operator=(const AliHLTSimulation&);
  /** standard destructor */
  ~AliHLTSimulation();

  /* current options */
  TString fOptions;                                                   //!transient

  /** base class for AliRoot HLT plugins */
  AliHLTPluginBase* fpPluginBase;                                     //!transient

  /** RAW reader instance for chains which need RAW data as input */
  AliRawReader* fpRawReader;                                          //!transient

  ClassDef(AliHLTSimulation, 2)
};

#define ALIHLTSIMULATION_LIBRARY             "libHLTsim.so"
#define ALIHLTSIMULATION_LIBRARY_VERSION     0
#define ALIHLTSIMULATION_CREATE_INSTANCE     "AliHLTSimulationCreateInstance"
#define ALIHLTSIMULATION_DELETE_INSTANCE     "AliHLTSimulationDeleteInstance"
#define ALIHLTSIMULATION_INIT                "AliHLTSimulationInit"
#define ALIHLTSIMULATION_RUN                 "AliHLTSimulationRun"
#define ALIHLTSIMULATION_GET_LIBRARY_VERSION "AliHLTSimulationGetLibraryVersion"

#ifdef __cplusplus
extern "C" {
#endif
  typedef AliHLTSimulation* (*AliHLTSimulationCreateInstance_t)();
  typedef int (*AliHLTSimulationDeleteInstance_t)(AliHLTSimulation* pSim);
  typedef int (*AliHLTSimulationInit_t)(AliHLTSimulation* pSim, AliRunLoader* pRunLoader, const char* options);
  typedef int (*AliHLTSimulationRun_t)(AliHLTSimulation* pSim, AliRunLoader* pRunLoader);
  typedef int (*AliHLTSimulationGetLibraryVersion_t)();

  /**
   * Create an instance of the AliHLTSimulation class
   */
  AliHLTSimulation* AliHLTSimulationCreateInstance();
  /**
   * Delete an instance of the AliHLTSimulation class
   */
  int AliHLTSimulationDeleteInstance(AliHLTSimulation* pSim);
  /**
   * Set options for an instance
   */
  int AliHLTSimulationInit(AliHLTSimulation* pSim, AliRunLoader* pRunLoader, const char* options);
  /**
   * Run simulation for an instance and run loader
   */
  int AliHLTSimulationRun(AliHLTSimulation* pSim, AliRunLoader* pRunLoader);
  /**
   * Get version no of the library/class interface
   */
  int AliHLTSimulationGetLibraryVersion();
#ifdef __cplusplus
}
#endif

#endif
