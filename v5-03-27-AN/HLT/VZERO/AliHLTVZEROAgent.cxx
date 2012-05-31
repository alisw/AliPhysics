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

/** @file   AliHLTVZEROAgent.cxx
    @author Jochen Thaeder <jochen@thaeder.de>
    @brief  Agent of the libAliHLTVZERO library
*/

#include <cassert>

#include "TSystem.h"

#include "AliHLTVZEROAgent.h"

#include "AliHLTErrorGuard.h"

// header files of library components
#include "AliHLTVZERORecoComponent.h"

// raw data handler of HLTOUT data
#include "AliHLTOUTHandlerEquId.h"
#include "AliHLTOUTHandlerEsdBranch.h"

/** global instance for agent registration */
AliHLTVZEROAgent gAliHLTVZEROAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTVZEROAgent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTVZEROAgent::AliHLTVZEROAgent() :
  AliHLTModuleAgent("VZERO") {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// #################################################################################
AliHLTVZEROAgent::~AliHLTVZEROAgent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTModuleAgent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTVZEROAgent::CreateConfigurations(AliHLTConfigurationHandler* handler,
					     AliRawReader* rawReader, AliRunLoader* runloader) const {
  // see header file for class documentation

  if (!handler) 
    return -EINVAL;

  if (rawReader || !runloader) {
    // AliSimulation: use the AliRawReaderPublisher if the raw reader is available
    // Alireconstruction: indicated by runloader==NULL, run always on raw data

    // -- Define the VZERO raw publisher
    // -----------------------------------
    TString arg("-equipmentid 3584 -datatype 'DDL_RAW ' 'VZRO' -dataspec 0x01");
    handler->CreateConfiguration("VZERO-DP_0", "AliRawReaderPublisher", NULL , arg.Data());

    // -- Define the VZERO reconstruction components
    // -----------------------------------------------
    handler->CreateConfiguration("VZERO-RECO", "VZEROReconstruction", "VZERO-DP_0", "");
  }
  else if (runloader && !rawReader) {
    // indicates AliSimulation with no RawReader available -> run on digits
    
    /* NOT Tested/ implemented yet
      handler->CreateConfiguration("DigitPublisher","AliLoaderPublisher",NULL,
      "-loader VZEROLoader -datatype 'ALITREED' 'VZRO'");
      handler->CreateConfiguration("Digit","VZEROReconstruction","DigitPublisher","");
    */
  }
  
  return 0;
}

// #################################################################################
const Char_t* AliHLTVZEROAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
							AliRunLoader* /*runloader*/) const {
  // see header file for class documentation

  // VZERO called only from the EsdConverter
  return NULL;
}

// #################################################################################
const Char_t* AliHLTVZEROAgent::GetRequiredComponentLibraries() const {
  // see header file for class documentation
  return "libAliHLTUtil.so libAliHLTVZERO.so";
}

// #################################################################################
Int_t AliHLTVZEROAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const {
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  
  pHandler->AddComponent(new AliHLTVZERORecoComponent);
  
  return 0;
}

// #################################################################################
Int_t AliHLTVZEROAgent::GetHandlerDescription(AliHLTComponentDataType dt, AliHLTUInt32_t spec,
					     AliHLTOUTHandlerDesc& desc) const {
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginVZERO)) {
    desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
      HLTInfo("module %s handles data block type %s specification %d (0x%x)", 
	      GetModuleId(), AliHLTComponent::DataType2Text(dt).c_str(), spec, spec);
      return 1;
  }

  // add TObject data blocks of type {ESD_CONT:VZRO} to ESD
  if (dt==(kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO)) {
      desc=AliHLTOUTHandlerDesc(kEsd, dt, GetModuleId());
      HLTInfo("module %s handles data block type %s specification %d (0x%x)", 
	      GetModuleId(), AliHLTComponent::DataType2Text(dt).c_str(), spec, spec);
      return 1;
  }

  return 0;
}

// #################################################################################
AliHLTOUTHandler* AliHLTVZEROAgent::GetOutputHandler(AliHLTComponentDataType dt,
						     AliHLTUInt32_t /*spec*/) {
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginVZERO)) {
    // use the default handler
    static AliHLTOUTHandlerEquId handler;
    return &handler;
  }

  if (dt==(kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO)) {
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
Int_t AliHLTVZEROAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance) {
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;
  
  // nothing to delete, the handler have been defined static
  return 0;
}

// #################################################################################
AliHLTModulePreprocessor* AliHLTVZEROAgent::GetPreprocessor() {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "GtePreProcessor not implemented for this module");
  return NULL;
}
