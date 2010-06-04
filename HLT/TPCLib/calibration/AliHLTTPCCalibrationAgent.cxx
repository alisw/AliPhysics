//-*- Mode: C++ -*-

// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCCalibrationAgent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  Agent of the libAliHLTTPCCalibration library
*/

#include "AliHLTTPCCalibrationAgent.h"

#include "AliHLTTPCDefinitions.h"
#include "AliHLTOUT.h"
#include "AliHLTOUTHandlerChain.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTPCParam.h"

// component header file

//#include "AliHLTTPCCalibCEComponent.h"
//#include "AliHLTTPCCalibPulserComponent.h"
//#include "AliHLTTPCCalibPedestalComponent.h"
//#include "AliHLTTPCCalibTracksComponent.h"

#include "AliHLTTPCCalibSeedMakerComponent.h"
#include "AliHLTTPCCalibTimeComponent.h"
#include "AliHLTTPCCalibTimeGainComponent.h"
#include "AliHLTTPCCalibrationComponent.h"


/** global instance for agent registration */
AliHLTTPCCalibrationAgent gAliHLTTPCCalibrationAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCalibrationAgent)
  
/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTTPCCalibrationAgent::AliHLTTPCCalibrationAgent()
  :
  AliHLTModuleAgent("TPCCalibration"),
  fRawDataHandler(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// #################################################################################
AliHLTTPCCalibrationAgent::~AliHLTTPCCalibrationAgent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 *                            
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTTPCCalibrationAgent::CreateConfigurations(AliHLTConfigurationHandler* /*handler*/,  AliRawReader* /*rawReader*/, AliRunLoader* /*runloader*/) const {
  // see header file for class documentation
  return 0;
}

// #################################################################################
const Char_t* AliHLTTPCCalibrationAgent::GetReconstructionChains(AliRawReader* /*rawReader*/, AliRunLoader* /*runloader*/) const {
  // see header file for class documentation

  return "";
}

// #################################################################################
const Char_t* AliHLTTPCCalibrationAgent::GetRequiredComponentLibraries() const {
  // see header file for class documentation
  return "";
}

// #################################################################################
Int_t AliHLTTPCCalibrationAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const {
  // see header file for class documentation

  if (!pHandler) return -EINVAL;

  //pHandler->AddComponent(new AliHLTTPCCalibCEComponent);
  //pHandler->AddComponent(new AliHLTTPCCalibPulserComponent);
  //pHandler->AddComponent(new AliHLTTPCCalibPedestalComponent);
  //pHandler->AddComponent(new AliHLTTPCCalibTracksComponent);
  pHandler->AddComponent(new AliHLTTPCCalibSeedMakerComponent);
  pHandler->AddComponent(new AliHLTTPCCalibTimeComponent);
  pHandler->AddComponent(new AliHLTTPCCalibTimeGainComponent);
  pHandler->AddComponent(new AliHLTTPCCalibrationComponent);
  
  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                            
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTModulePreprocessor* AliHLTTPCCalibrationAgent::GetPreprocessor(){
  // see header file for class documentation
  return NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                            
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTTPCCalibrationAgent::GetHandlerDescription(AliHLTComponentDataType /*dt*/,
					    AliHLTUInt32_t /*spec*/,
					    AliHLTOUTHandlerDesc& /*desc*/) const {
  // see header file for class documentation
  return 0;
}

// #################################################################################
AliHLTOUTHandler* AliHLTTPCCalibrationAgent::GetOutputHandler(AliHLTComponentDataType /*dt*/,  AliHLTUInt32_t /*spec*/){
  // see header file for class documentation
  return NULL;
}

// #################################################################################
Int_t AliHLTTPCCalibrationAgent::DeleteOutputHandler(AliHLTOUTHandler* /*pInstance*/){
  // see header file for class documentation
  return 0;
}
