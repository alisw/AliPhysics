//-*- Mode: C++ -*-

// $Id: AliHLTJETAgent.cxx $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTJETAgent.cxx
    @author Jochen Thaeder
    @date   13.02.2009
    @brief  Agent of the libAliHLTJET library
*/

#include <cassert>
#include "AliHLTJETAgent.h"
#include "AliHLTOUT.h"

// component header file
#include "AliHLTJETConeJetComponent.h"
#include "AliHLTJETAnalysisComponent.h"

#ifdef HAVE_FASTJET
#include "AliHLTJETFastJetComponent.h"
#endif



/** global instance for agent registration */
AliHLTJETAgent gAliHLTJETAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETAgent)
  
/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTJETAgent::AliHLTJETAgent()
  :
  AliHLTModuleAgent("JET"),
  fRawDataHandler(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// #################################################################################
AliHLTJETAgent::~AliHLTJETAgent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 *                            
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETAgent::CreateConfigurations(AliHLTConfigurationHandler* /*handler*/,
					   AliRawReader* /*rawReader*/,
					   AliRunLoader* /*runloader*/) const {
  // see header file for class documentation
  return 0;
}

// #################################################################################
const Char_t* AliHLTJETAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						      AliRunLoader* /*runloader*/) const {
  // see header file for class documentation

  return "";
}

// #################################################################################
const Char_t* AliHLTJETAgent::GetRequiredComponentLibraries() const {
  // see header file for class documentation
  return "";
}

// #################################################################################
Int_t AliHLTJETAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const {
  // see header file for class documentation

  assert(pHandler);
  if (!pHandler) return -EINVAL;

  pHandler->AddComponent(new AliHLTJETConeJetComponent);
  pHandler->AddComponent(new AliHLTJETAnalysisComponent);
#ifdef HAVE_FASTJET
  pHandler->AddComponent(new AliHLTJETFastJetComponent);
#endif

  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                            
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTModulePreprocessor* AliHLTJETAgent::GetPreprocessor() {
  // see header file for class documentation
  return NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                            
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETAgent::GetHandlerDescription(AliHLTComponentDataType /*dt*/,
					    AliHLTUInt32_t /*spec*/,
					    AliHLTOUTHandlerDesc& /*desc*/) const {
  // see header file for class documentation
  
  // Handlers for JET data. 
  /*
    if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginJETSDD)) {
    desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
    HLTInfo("module %s handles data block type %s specification %d (0x%x)", 
    GetModuleId(), AliHLTComponent::DataType2Text(dt).c_str(), spec, spec);
    return 1;
    }
  */
  return 0;
}

// #################################################################################
AliHLTOUTHandler* AliHLTJETAgent::GetOutputHandler(AliHLTComponentDataType /*dt*/,
						   AliHLTUInt32_t /*spec*/) {
  /*
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeJet|kAliHLTDataOriginANY)) {
  // use the default handler
  if (!fRawDataHandler) {
  fRawDataHandler=new AliHLTOUTSDDRawDataHandler;
  }
  return fRawDataHandler;
  }
  */
  return NULL;
}

// #################################################################################
Int_t AliHLTJETAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance) {
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  /**
  if (pInstance==fRawDataHandler) {
    delete fRawDataHandler;
    fRawDataHandler=NULL;
    return 0;
  }

  delete pInstance;
  */
  return 0;
}

/*
int AliHLTJETAgent::AliHLTOUTSDDRawDataHandler::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;
  static int errorCount=0;
  const int maxErrorCount=10;
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
  int iResult=pData->GetDataBlockDescription(dt, spec);
  if (iResult>=0) {
    if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginJETSDD)) {
      int ddlOffset=256;//AliDAQ::DdlIDOffset("JETSDD");
      int numberOfDDLs=24;//AliDAQ::NumberOfDdls("JETSDD");
      int ddlNo=0;
      for (;ddlNo<32 && ddlNo<numberOfDDLs; ddlNo++) {
	if (spec&(0x1<<ddlNo)) break;
      }
      if (ddlNo>=32 || ddlNo>=numberOfDDLs) {
	HLTError("invalid specification 0x%08x: can not extract DDL id for data block %s", spec, AliHLTComponent::DataType2Text(dt).c_str());
	iResult=-ENODEV;
      } else if (spec^(0x1<<ddlNo)) {
	iResult=-EEXIST;
	HLTError("multiple links set in specification 0x%08x: can not extract DDL id for data block %s", spec, AliHLTComponent::DataType2Text(dt).c_str());
      } else {
	iResult=ddlOffset+ddlNo;
      }
    } else {
      if (errorCount++<10) {
	HLTError("wrong data type: expecting %s, got %s; %s",
		 AliHLTComponent::DataType2Text(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginJETSDD).c_str(),
		 AliHLTComponent::DataType2Text(dt).c_str(),
		   errorCount==maxErrorCount?"suppressing further error messages":"");
      }
      iResult=-EFAULT;
    }
  }
  return iResult;
}
*/
