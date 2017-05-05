// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <jochen@thaeder.de>                    *
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

/** @file   AliHLTT0Agent.cxx
    @author Alla Maevskaya <Alla.Maevskaya@cern.ch>
    @brief  Agent of the libAliHLTT0 library
*/

#include <cassert>

#include "TSystem.h"
#include "AliDAQ.h"

#include "AliHLTT0Agent.h"

#include "AliHLTErrorGuard.h"

// header files of library components
#include "AliHLTT0RecoComponent.h"

// raw data handler of HLTOUT data
#include "AliHLTOUTHandlerEquId.h"
#include "AliHLTOUTHandlerEsdBranch.h"

/** global instance for agent registration */
AliHLTT0Agent gAliHLTT0Agent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTT0Agent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTT0Agent::AliHLTT0Agent() :
  AliHLTModuleAgent("T0") {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// #################################################################################
AliHLTT0Agent::~AliHLTT0Agent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTModuleAgent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

UInt_t AliHLTT0Agent::GetDetectorMask() const
{
  return AliDAQ::kT0;
}

// #################################################################################
Int_t AliHLTT0Agent::CreateConfigurations(AliHLTConfigurationHandler* handler,
					     AliRawReader* rawReader, AliRunLoader* runloader) const {
  // see header file for class documentation

  if (!handler) 
    return -EINVAL;

  if (rawReader || !runloader) {
    // AliSimulation: use the AliRawReaderPublisher if the raw reader is available
    // Alireconstruction: indicated by runloader==NULL, run always on raw data

    // -- Define the T0 raw publisher
    // -----------------------------------
    TString arg("-equipmentid 3328 -datatype 'DDL_RAW ' 'TZRO' -dataspec 0x01");
    handler->CreateConfiguration("T0-DP", "AliRawReaderPublisher", NULL , arg.Data());

    // -- Define the T0 reconstruction components
    // -----------------------------------------------
    handler->CreateConfiguration("T0-RECO", "T0Reconstruction", "T0-DP ITS-SPD-vertexer", "");
  }
  else if (runloader && !rawReader) {
    // indicates AliSimulation with no RawReader available -> run on digits
    
    // NOT Tested/ implemented yet
    /*
    handler->CreateConfiguration("V0DigitPublisher","AliLoaderPublisher",NULL,
				 "-loader T0Loader -datatype 'ALITREED' 'VZRO'");
    //handler->CreateConfiguration("Digit","T0Reconstruction","DigitPublisher","");    
    handler->CreateConfiguration("T0-RECO", "T0Reconstruction", "V0DigitPublisher", "");
    */
  }
  
  return 0;
}

// #################################################################################
const Char_t* AliHLTT0Agent::GetReconstructionChains(AliRawReader* /*rawReader*/,
							AliRunLoader* /*runloader*/) const {
  // see header file for class documentation

  // T0 called only from the EsdConverter
  return NULL;
}

// #################################################################################
const Char_t* AliHLTT0Agent::GetRequiredComponentLibraries() const {
  // see header file for class documentation
  return "libAliHLTUtil.so libAliHLTT0.so";
}

// #################################################################################
Int_t AliHLTT0Agent::RegisterComponents(AliHLTComponentHandler* pHandler) const {
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  
  pHandler->AddComponent(new AliHLTT0RecoComponent);
  
  return 0;
}

// #################################################################################
Int_t AliHLTT0Agent::GetHandlerDescription(AliHLTComponentDataType dt, AliHLTUInt32_t spec,
					     AliHLTOUTHandlerDesc& desc) const {
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginT0)) {
    desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
      HLTInfo("module %s handles data block type %s specification %d (0x%x)", 
	      GetModuleId(), AliHLTComponent::DataType2Text(dt).c_str(), spec, spec);
      return 1;
  }

  // add TObject data blocks of type {ESD_CONT:VZRO} to ESD
  if (dt==(kAliHLTDataTypeESDContent|kAliHLTDataOriginT0)) {
      desc=AliHLTOUTHandlerDesc(kEsd, dt, GetModuleId());
      HLTInfo("module %s handles data block type %s specification %d (0x%x)", 
	      GetModuleId(), AliHLTComponent::DataType2Text(dt).c_str(), spec, spec);
      return 1;
  }

  return 0;
}

// #################################################################################
AliHLTOUTHandler* AliHLTT0Agent::GetOutputHandler(AliHLTComponentDataType dt,
						     AliHLTUInt32_t /*spec*/) {
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginT0)) {
    // use the default handler
    static AliHLTOUTHandlerEquId handler;
    return &handler;
  }

  if (dt==(kAliHLTDataTypeESDContent|kAliHLTDataOriginT0)) {
    // use AliHLTOUTHandlerEsdBranch handler to add the TObject
    // to the ESD branch
    // Note: the object should have an appropriate name returned
    // by GetName(). Use SetName() to prepare the object before streaming
    static AliHLTOUTHandlerEsdBranch handler;
    return &handler;
  }

  return NULL;
}

// #################################################################################
Int_t AliHLTT0Agent::DeleteOutputHandler(AliHLTOUTHandler* pInstance) {
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;
  
  // nothing to delete, the handler have been defined static
  return 0;
}

// #################################################################################
AliHLTModulePreprocessor* AliHLTT0Agent::GetPreprocessor() {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "GtePreProcessor not implemented for this module");
  return NULL;
}
