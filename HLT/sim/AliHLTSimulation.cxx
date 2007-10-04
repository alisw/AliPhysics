// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTSimulation.cxx
    @author Matthias Richter
    @date   
    @brief  Binding class for HLT simulation in AliRoot. */

#include <cassert>
#include <cerrno>
#include "AliHLTSimulation.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliHLTSystem.h"

#if ALIHLTSIMULATION_LIBRARY_VERSION != LIBHLTSIM_VERSION
#error library version in header file and lib*.pkg do not match
#endif

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSimulation);

AliHLTSimulation::AliHLTSimulation()
  :
  fOptions(),
  fpSystem(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTSimulation::~AliHLTSimulation()
{
  // see header file for function documentation
  if (fpSystem) {
    delete fpSystem;
  }
  fpSystem=NULL;
}

AliHLTSimulation* AliHLTSimulation::CreateInstance()
{
  // see header file for function documentation
  return new AliHLTSimulation;
}

int AliHLTSimulation::DeleteInstance(AliHLTSimulation* pSim)
{
  // see header file for function documentation
  assert(pSim!=NULL);
  delete pSim;
  return 0;
}

int AliHLTSimulation::Init(AliRunLoader* pRunLoader, const char* options)
{
  // init the simulation
  fOptions=options;

  if (!fpSystem) fpSystem=new AliHLTSystem;
  if (!fpSystem) {
    AliError("can not create AliHLTSystem object");
    return -ENOMEM;
  }
  if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
    AliError("HLT system in error state");
    return -EFAULT;
  }

  if (fpSystem->ScanOptions(options)<0) {
    AliError("error setting options for HLT system");
    return -EINVAL;	
  }

  if (!fpSystem->CheckStatus(AliHLTSystem::kReady)) {
    if ((fpSystem->Configure(pRunLoader))<0) {
      AliError("error during HLT system configuration");
      return -EFAULT;
    }
  }

  return 0;
}


int AliHLTSimulation::Run(AliRunLoader* pRunLoader)
{
  // HLT reconstruction for simulated data  
  if(!pRunLoader) {
    AliError("Missing RunLoader! 0x0");
    return -EINVAL;
  }

  Int_t nEvents = pRunLoader->GetNumberOfEvents();
  int iResult=0;

  if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
    AliError("HLT system in error state");
    return -EFAULT;
  }
  if ((iResult=fpSystem->Reconstruct(nEvents, pRunLoader, NULL))>=0) {
  }
  return iResult;
}


AliHLTSimulation* AliHLTSimulationCreateInstance()
{
  // see header file for function documentation
  return AliHLTSimulation::CreateInstance();
}

int AliHLTSimulationDeleteInstance(AliHLTSimulation* pSim)
{
  // see header file for function documentation
  return AliHLTSimulation::DeleteInstance(pSim);
}

int AliHLTSimulationInit(AliHLTSimulation* pSim, AliRunLoader* pRunLoader, const char* options)
{
  assert(pSim!=NULL);
  if (pSim) {
    return pSim->Init(pRunLoader, options);
  }
  return -ENODEV;
}

int AliHLTSimulationRun(AliHLTSimulation* pSim, AliRunLoader* pRunLoader)
{
  assert(pSim!=NULL);
  if (pSim) {
    return pSim->Run(pRunLoader);
  }
  return -ENODEV;
}

int AliHLTSimulationGetLibraryVersion()
{
  // see header file for function documentation
  return LIBHLTSIM_VERSION;
}
