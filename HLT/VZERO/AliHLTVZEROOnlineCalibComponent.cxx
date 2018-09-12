// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: David Rohr <drohr@cern.ch>                            *
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

/** @file    AliHLTVZEROOnlineCalibComponent.cxx
    @author  David Rohr <drohr@cern.ch>
    @brief   VZERO online calib component
*/

#include "TTree.h"
#include "TMap.h"
#include "TObjString.h"
#include "TDatime.h"
#include "TH1F.h"

#include "AliLog.h"
#include "AliRunInfo.h"
#include "AliGRPObject.h"
#include "AliGeomManager.h"

#include "AliVZERORecoParam.h"

#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTVZEROOnlineCalibComponent.h"
#include "AliESDVZERO.h"

#include "AliESDVZEROfriend.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTVZEROOnlineCalibComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTVZEROOnlineCalibComponent::AliHLTVZEROOnlineCalibComponent() :
  AliHLTProcessor(),
  fRunInfo(NULL),  
  fVZERORecoParam(NULL),
  fHistMult("vzero", "vzero", 100, 0, 1000),
  fEventModulo(0)
  {
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
  fHistMult.SetDirectory(0);
}

// #################################################################################
AliHLTVZEROOnlineCalibComponent::~AliHLTVZEROOnlineCalibComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTVZEROOnlineCalibComponent::GetComponentID() { 
  // see header file for class documentation
  return "VZEROOnlineCalib";
}

// #################################################################################
void AliHLTVZEROOnlineCalibComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeESDContent | kAliHLTDataOriginVZERO);
}

// #################################################################################
AliHLTComponentDataType AliHLTVZEROOnlineCalibComponent::GetOutputDataType() 
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram|kAliHLTDataOriginOut;
}

// #################################################################################
void AliHLTVZEROOnlineCalibComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 3000;
  inputMultiplier = 1;
}

// #################################################################################
void AliHLTVZEROOnlineCalibComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  // see header file for class documentation

  if (!targetMap) return;
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
  targetMap->Add(new TObjString("VZERO/Trigger/Data"),
		 new TObjString("VZERO calibration object"));
  return;
}

// #################################################################################
AliHLTComponent* AliHLTVZEROOnlineCalibComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTVZEROOnlineCalibComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTVZEROOnlineCalibComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation
 
  Int_t iResult=0;

  // -- Load GeomManager
  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }
  
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

    // implement further initialization
  } while (0);

  if (iResult<0) {
    // implement cleanup

    if (fVZERORecoParam)
      delete fVZERORecoParam;
    fVZERORecoParam = NULL;

    if (fRunInfo)
      delete fRunInfo;
    fRunInfo = NULL;
  }

  return iResult;
}

// #################################################################################
Int_t AliHLTVZEROOnlineCalibComponent::ScanConfigurationArgument(Int_t /*argc*/, const Char_t** argv) {
  Int_t ii =0;
  TString argument=argv[ii];

  if (argument.IsNull()) return 0;

  return 0;
}

// #################################################################################
Int_t AliHLTVZEROOnlineCalibComponent::DoDeinit() {
  // see header file for class documentation

  if (fVZERORecoParam)
    delete fVZERORecoParam;
  fVZERORecoParam = NULL;
  
  if (fRunInfo)
    delete fRunInfo;
  fRunInfo = NULL;
  
  return 0;
}

// #################################################################################
Int_t AliHLTVZEROOnlineCalibComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation
  Int_t iResult=0;

  // -- Only use data event
  if (!IsDataEvent()) 
    return 0;

    
  //Get input data, and proceed if we got a VZERO ESD object
  const AliESDVZERO* esdVZERO = dynamic_cast<const AliESDVZERO*>(GetFirstInputObject(kAliHLTDataTypeESDContent | kAliHLTDataOriginVZERO, "AliESDVZERO"));
  if (esdVZERO)
  {
    //Just some example, create histogram with VZero multiplicity.
    float vZEROMultiplicity = 0.f;
    for (int i = 0;i < 64;i++) vZEROMultiplicity += esdVZERO->GetMultiplicity(i);
    fHistMult.Fill(vZEROMultiplicity);
    //vZEROTriggerChargeA = esdVZERO->GetTriggerChargeA();
    //vZEROTriggerChargeC = esdVZERO->GetTriggerChargeC();
  }
  
  //If the histogram is not empty, we send it out every 16th event (to collect some statistics).
  //Depending on the pushback period set for the component, it might not be send out every time, but only after a certain amount of time. (order of 3 minutes)
  //We check whether it was really sent out, and only if so, we reset the histogram.
  //The ZMQ merging component that sits at the end of the chain will receive all histograms from all concurrent VZEROOnlineCalib components, and merge them to the final histogram.
  if ((++fEventModulo % 16 == 0) && fHistMult.GetEntries() && PushBack(&fHistMult, kAliHLTDataTypeHistogram|kAliHLTDataOriginHLT) > 0) fHistMult.Reset();
  
  return iResult;
}

// #################################################################################
Int_t AliHLTVZEROOnlineCalibComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
  // see header file for class documentation

  Int_t iResult=0;
  return iResult;
}

