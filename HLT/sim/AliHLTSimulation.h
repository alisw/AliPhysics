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
 * @defgroup alihlt_simulation HLT simulation in AliRoot
 * This section describes the the simulation of the HLT in AliRoot.
 */

#include "TObject.h"
#include "TString.h"
class AliRunLoader;
class AliHLTSystem;
class AliRawReader;

/**
 * @class AliHLTSimulation
 * Base class of HLT data processing simulations.
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

  /* HLT steering object */
  AliHLTSystem* fpSystem;                                             //!transient

  /* RAW reader instance for chains which need RAW data as input */
  AliRawReader* fpRawReader;                                            //!transient

  ClassDef(AliHLTSimulation, 1)
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
