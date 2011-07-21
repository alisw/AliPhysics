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

/** @file    AliHLTVZERORecoComponent.cxx
    @author  Jochen Thaeder <jochen@thaeder.de>
    @brief   VZERO reconstruction component
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "TTree.h"
#include "TMap.h"
#include "TObjString.h"
#include "TDatime.h"

#include "AliLog.h"
#include "AliRunInfo.h"
#include "AliGRPObject.h"
#include "AliRawReaderMemory.h"
#include "AliGeomManager.h"

#include "AliVZERORecoParam.h"
#include "AliVZEROReconstructor.h"

#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTVZERORecoComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTVZERORecoComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTVZERORecoComponent::AliHLTVZERORecoComponent() :
  AliHLTProcessor(),
  fRunInfo(NULL),  
  fVZERORecoParam(NULL),
  fVZEROReconstructor(NULL),
  fRawReader(NULL) {
  // an example component which implements the ALICE HLT processor
  // interface and does some analysis on the input raw data
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

// #################################################################################
AliHLTVZERORecoComponent::~AliHLTVZERORecoComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTVZERORecoComponent::GetComponentID() { 
  // see header file for class documentation
  return "VZEROReconstruction";
}

// #################################################################################
void AliHLTVZERORecoComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginVZERO);
}

// #################################################################################
AliHLTComponentDataType AliHLTVZERORecoComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO;
}

// #################################################################################
void AliHLTVZERORecoComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 1000;
  inputMultiplier = 0.5;
}

// #################################################################################
void AliHLTVZERORecoComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  // see header file for class documentation

  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigVZERO/VZEROReconstruction"),
		 new TObjString("configuration object"));

  targetMap->Add(new TObjString("GRP/GRP/Data"),
		 new TObjString("GRP object - run information"));
  targetMap->Add(new TObjString("GRP/CTP/CTPtiming"),
		 new TObjString("GRP object - CTP information"));
  targetMap->Add(new TObjString("GRP/CTP/TimeAlign"),
		 new TObjString("GRP object - CTP information"));
  targetMap->Add(new TObjString("GRP/Calib/LHCClockPhase"),
		 new TObjString("GRP object - time calibration"));

  targetMap->Add(new TObjString("VZERO/Calib/Data"),
		 new TObjString("VZERO calibration object"));
  targetMap->Add(new TObjString("VZERO/Calib/TimeDelays"),
		 new TObjString("VZERO calibration object"));
  targetMap->Add(new TObjString("VZERO/Calib/TimeSlewing"),
		 new TObjString("VZERO calibration object"));

  return;
}

// #################################################################################
AliHLTComponent* AliHLTVZERORecoComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTVZERORecoComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTVZERORecoComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation

  Int_t iResult=0;

  // -- Load GeomManager
  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }
  
  // -- Read configuration object : HLT/ConfigVZERO/VZEROReconstruction
  TString cdbPath="HLT/ConfigVZERO/";
  cdbPath+=GetComponentID();
  iResult=ConfigureFromCDBTObjString(cdbPath);

  // -- Read the component arguments
  if (iResult>=0) {
    iResult=ConfigureFromArgumentString(argc, argv);
  }

  // -- Get AliRunInfo variables
  // -----------------------------
  TObject* pOCDBEntry=LoadAndExtractOCDBObject("GRP/GRP/Data");
  AliGRPObject* pGRP=pOCDBEntry?dynamic_cast<AliGRPObject*>(pOCDBEntry):NULL;
  
  TString beamType = "";
  TString lhcState = "";
  TString runType = "";
  Float_t beamEnergy = 0.;
  UInt_t activeDetectors = 0;
  
  if (pGRP) {
    lhcState        = pGRP->GetLHCState(); 	  	   
    beamType        = pGRP->GetBeamType(); 
    runType         = pGRP->GetRunType(); 
    beamEnergy      = pGRP->GetBeamEnergy();
    activeDetectors = pGRP->GetDetectorMask();
  }
  
  // -- Initialize members
  // -----------------------
  do {
    if (iResult<0) break;

    fRawReader = new AliRawReaderMemory;
    if (!fRawReader) {
      iResult=-ENOMEM;
      break;
    }

    // AliGRPManager grpMan;
    // Bool_t status       = grpMan.ReadGRPEntry(); // Read the corresponding OCDB entry
    // status              = grpMan.SetMagField();  // Set global field instanton
    // AliRunInfo *runInfo = grpMan.GetRunInfo();   // Get instance of run info

    fRunInfo = new AliRunInfo(lhcState.Data(), beamType.Data(),
			      beamEnergy, runType.Data(), activeDetectors);
    if (!fRunInfo) {
      iResult=-ENOMEM;
      break;
    }

    fVZERORecoParam = new AliVZERORecoParam;
    if (!fVZERORecoParam) {
      iResult=-ENOMEM;
      break;
    }  

    fVZEROReconstructor = new AliVZEROReconstructor;
    if (!fVZEROReconstructor) {
      iResult=-ENOMEM;
      break;
    }

    // implement further initialization
  } while (0);

  if (iResult<0) {
    // implement cleanup

    if (fRawReader) 
      delete fRawReader;
    fRawReader = NULL;

    if (fVZERORecoParam)
      delete fVZERORecoParam;
    fVZERORecoParam = NULL;

    if (fVZEROReconstructor)
      delete fVZEROReconstructor;
    fVZEROReconstructor = NULL;

    if (fRunInfo)
      delete fRunInfo;
    fRunInfo = NULL;
  }

  if (iResult>=0) {
    fVZEROReconstructor->SetRunInfo(fRunInfo);
    fVZEROReconstructor->Init();

    fVZEROReconstructor->SetRecoParam(fVZERORecoParam);
  }

  return iResult;
}

// #################################################################################
Int_t AliHLTVZERORecoComponent::ScanConfigurationArgument(Int_t /*argc*/, const Char_t** argv) {
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.

  Int_t ii =0;
  TString argument=argv[ii];

  if (argument.IsNull()) return 0;

  return 0;
}

// #################################################################################
Int_t AliHLTVZERORecoComponent::DoDeinit() {
  // see header file for class documentation

  if (fRawReader) 
    delete fRawReader;
  fRawReader = NULL;
  
  if (fVZERORecoParam)
    delete fVZERORecoParam;
  fVZERORecoParam = NULL;
  
  if (fVZEROReconstructor)
    delete fVZEROReconstructor;
  fVZEROReconstructor = NULL;
  
  if (fRunInfo)
    delete fRunInfo;
  fRunInfo = NULL;
  
  return 0;
}

// #################################################################################
Int_t AliHLTVZERORecoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation

  Int_t iResult=0;

  // -- Only use data event
  if (!IsDataEvent()) 
    return 0;

  // -- Get VZERO raw dat a input block and set up the rawreader
  const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginVZERO);
  if (!pBlock) {
    ALIHLTERRORGUARD(1, "No VZERO input block at event %d", GetEventCount());
    return 0;
  }
  
  // -- Add input block to raw reader
  if (!fRawReader->SetMemory((UChar_t*) pBlock->fPtr, pBlock->fSize )){
    HLTError("Could not add buffer of data block  %s, 0x%08x to rawreader",
	     DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification);
    iResult = -1;
  }
  
  TTree *digitsTree = new TTree("D", "Digits Tree");
  if (!digitsTree) {
    iResult=-ENOMEM;
  }

  if (iResult >= 0) {

    // -- Set VZERO EquipmentID
    fRawReader->SetEquipmentID(3584);
  
    // -- 1. step VZERO reconstruction
    fVZEROReconstructor->ConvertDigits(fRawReader, digitsTree);

    // -- 2. step VZERO reconstruction -- fill AliESDVZERO object
    fVZEROReconstructor->FillESD(digitsTree, NULL, NULL);

    AliESDVZERO *esdVZERO = fVZEROReconstructor->GetESDVZERO();
    
    // Send info every 10 s
    const TDatime time;    
    static UInt_t lastTime=0;
    if (time.Get()-lastTime>10) {
      lastTime=time.Get();
      HLTInfo("VZERO Multiplicity A %f - C %f", esdVZERO->GetMTotV0A(), esdVZERO->GetMTotV0A() );
    }

    // -- Send AliESDVZERO
    PushBack(esdVZERO, kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO,0);
  }
  
  // -- Clean up
  delete digitsTree;
  fRawReader->ClearBuffers();   

  return iResult;
}

// #################################################################################
Int_t AliHLTVZERORecoComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
  // see header file for class documentation

  Int_t iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigVZERO/";
    cdbPath+=GetComponentID();
  }

  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}

// #################################################################################
Int_t AliHLTVZERORecoComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  return 0;
}
