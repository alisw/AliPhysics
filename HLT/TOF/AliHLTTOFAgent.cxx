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

/** @file   AliHLTTOFAgent.cxx
    @author Jochen Thaeder <jochen@thaeder.de>
    @brief  Agent of the libAliHLTTOF library
*/

#include <cassert>

#include "TSystem.h"
#include "AliDAQ.h"

#include "AliHLTTOFAgent.h"

#include "AliHLTErrorGuard.h"

// header files of library components
#include "AliHLTTOFRecoComponent.h"

// raw data handler of HLTOUT data
#include "AliHLTOUTHandlerEquId.h"
#include "AliHLTOUTHandlerEsdBranch.h"

/** global instance for agent registration */
AliHLTTOFAgent gAliHLTTOFAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTOFAgent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTTOFAgent::AliHLTTOFAgent() :
  AliHLTModuleAgent("TOF") {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// #################################################################################
AliHLTTOFAgent::~AliHLTTOFAgent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTModuleAgent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

UInt_t AliHLTTOFAgent::GetDetectorMask() const
{
  return AliDAQ::kTOF;
}

// #################################################################################
Int_t AliHLTTOFAgent::CreateConfigurations(AliHLTConfigurationHandler* handler,
					     AliRawReader* rawReader, AliRunLoader* runloader) const {
  // see header file for class documentation

  if (!handler) 
    return -EINVAL;

  if (rawReader || !runloader) {
    // AliSimulation: use the AliRawReaderPublisher if the raw reader is available
    // Alireconstruction: indicated by runloader==NULL, run always on raw data

    // -- Define the TOF raw publisher
    // -----------------------------------
    TString recoStr;
    for (Int_t iDDL = 0; iDDL < 72; iDDL++){

      
      TString arg(Form("-equipmentid %d -datatype 'DDL_RAW ' 'TOF ' -dataspec 0x%x", 1280+iDDL, iDDL));
      handler->CreateConfiguration(Form("TOF-DP_%02d", iDDL), "AliRawReaderPublisher", NULL , arg.Data());
      recoStr += Form("TOF-DP_%02d ", iDDL);

    }
      // -- Define the TOF reconstruction components
      // -----------------------------------------------
    handler->CreateConfiguration("TOF-RECO", "TOFReconstruction", recoStr.Data(), "");
  }
  else if (runloader && !rawReader) {
    // indicates AliSimulation with no RawReader available -> run on digits
    
    // NOT Tested/ implemented yet
    /*
    handler->CreateConfiguration("TOFDigitPublisher","AliLoaderPublisher",NULL,
				 "-loader TOFLoader -datatype 'ALITREED' 'VZRO'");
    //handler->CreateConfiguration("Digit","TOFReconstruction","DigitPublisher","");    
    handler->CreateConfiguration("TOF-RECO", "TOFReconstruction", "TOFDigitPublisher", "");
    */
  }
  
  return 0;
}

// #################################################################################
const Char_t* AliHLTTOFAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
							AliRunLoader* /*runloader*/) const {
  // see header file for class documentation

  // TOF called only from the EsdConverter
  return NULL;
}

// #################################################################################
const Char_t* AliHLTTOFAgent::GetRequiredComponentLibraries() const {
  // see header file for class documentation
  return "libAliHLTUtil.so libAliHLTTOF.so";
}

// #################################################################################
Int_t AliHLTTOFAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const {
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  
  pHandler->AddComponent(new AliHLTTOFRecoComponent);
  
  return 0;
}

// #################################################################################
Int_t AliHLTTOFAgent::GetHandlerDescription(AliHLTComponentDataType dt, AliHLTUInt32_t spec,
					     AliHLTOUTHandlerDesc& desc) const {
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTOF)) {
    desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
      HLTInfo("module %s handles data block type %s specification %d (0x%x)", 
	      GetModuleId(), AliHLTComponent::DataType2Text(dt).c_str(), spec, spec);
      return 1;
  }

  /*
  // add TObject data blocks of type {ESD_CONT:TOF} to ESD
  if (dt==(kAliHLTDataTypeESDContent|kAliHLTDataOriginTOF)) {
      desc=AliHLTOUTHandlerDesc(kEsd, dt, GetModuleId());
      HLTInfo("module %s handles data block type %s specification %d (0x%x)", 
	      GetModuleId(), AliHLTComponent::DataType2Text(dt).c_str(), spec, spec);
      return 1;
  }
*/
  return 0;
}

// #################################################################################
AliHLTOUTHandler* AliHLTTOFAgent::GetOutputHandler(AliHLTComponentDataType dt,
						     AliHLTUInt32_t /*spec*/) {
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTOF)) {
    // use the default handler
    static AliHLTOUTHandlerEquId handler;
    return &handler;
  }

  /*
  if (dt==(kAliHLTDataTypeESDContent|kAliHLTDataOriginTOF)) {
    // use AliHLTOUTHandlerEsdBranch handler to add the TObject
    // to the ESD branch
    // Note: the object should have an appropriate name returned
    // by GetName(). Use SetName() to prepare the object before streaming
    static AliHLTOUTHandlerEsdBranch handler;
    return &handler;
  }
  */

  return NULL;
}

// #################################################################################
Int_t AliHLTTOFAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance) {
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;
  
  // nothing to delete, the handler have been defined static
  return 0;
}

// #################################################################################
AliHLTModulePreprocessor* AliHLTTOFAgent::GetPreprocessor() {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "GetPreProcessor not implemented for this module (TOF)");
  return NULL;
}
