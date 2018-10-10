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
#include "AliHLTCTPData.h"
#include "AliESDVZERO.h"
#include "AliESDVertex.h"

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
fVZERORecoParam(NULL)
{
    // a component meant to extract simple V0 amplitudes for an
    // online monitoring of V0 conditions. The vast majority of
    // this monitoring can only happen when some statistics is
    // accumulated and that's not here: that's in a post-
    // processing macro
    //
    // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
    //
    // NOTE: all helper classes should be instantiated in DoInit()
    
    //Initialize
    TString lDetTypes[3] = {"V0M", "V0A", "V0C"};
    TString lReadoutTypes[2] = {"", "Online"};
    TString lTrigTypes[4] = {"", "_Central", "_SemiCentral", "_C0V0H"};
    TString lTitles[2] = {"amplitude", "trigger charge"};
    
    for(Int_t i=0; i<24; i++){
        TString lName = Form("fHistMult%s%s%s", lDetTypes[i%3].Data(), lReadoutTypes[(i/3)%2].Data(), lTrigTypes[i/6].Data());
        TString lTitle = Form("%s %s", lDetTypes[i%3].Data(), lTitles[(i/3)%2].Data());
        
        fHistMult[i].SetName(lName.Data());
        fHistMult[i].SetTitle(lTitle.Data());
        fHistMult[i].SetYTitle("Event Count");
        fHistMult[i].SetXTitle(lTitles[(i/3)%2].Data());
        
        fHistMult[i].SetDirectory(0);
    }
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
    list.push_back(kAliHLTDataTypeESDVertex | kAliHLTDataOriginITSSPD);
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
    constBase = 50000;
    inputMultiplier = 0;
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
        //Init the CTP data
        if (SetupCTPData() == -ENOMEM)
        {
            HLTError("could not SetupCTPData(); ENOMEM");
            return -ENOMEM;
        }
        //implement further initialization here
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
    
    //Setup histograms based on knowledge regarding system
    Double_t lMaxV0Signal = 10000;
    if (beamType == "Pb-Pb" ||
        beamType == "PbPb" ||
        beamType == "A-A" ||
        beamType == "AA" )
        lMaxV0Signal = 50000; //reasonable range for nuclear collisions
    
    for(Int_t i=0; i<24; i++){
        fHistMult[i].Reset();
        fHistMult[i].SetBins(10000, 0, lMaxV0Signal);
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
    
    //Get vertex for quick check
    const AliESDVertex* itsSpdVertex = dynamic_cast<const AliESDVertex*>(GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITSSPD, "AliESDVertex"));
    
    bool lIsSemiCentral = true, lIsCentral = true, lIsMinBias = true, lIsC0V0H = true, lIsV0DecisionOK = true, lIsVertexPositionGood = true;
    const AliHLTCTPData* ctp = CTPData();
    if (ctp)
    {
        lIsMinBias=ctp->MatchTriggerRE("C0TVX");
        lIsCentral=ctp->MatchTriggerRE("C0V0M");
        lIsSemiCentral=ctp->MatchTriggerRE("C0VSC");
        lIsC0V0H=ctp->MatchTriggerRE("C0V0H");
    }
    if (esdVZERO && itsSpdVertex)
    {
        //printf("Vertex position: %f %f %f\n", itsSpdVertex->GetX(), itsSpdVertex->GetY(), itsSpdVertex->GetZ());
        //Rough V0 decision check (equiv to phys. sel.)
        if (esdVZERO->GetV0ADecision()!=AliVVZERO::kV0BB) lIsV0DecisionOK = false;
        if (esdVZERO->GetV0CDecision()!=AliVVZERO::kV0BB) lIsV0DecisionOK = false;
        if ( TMath::Abs(itsSpdVertex->GetZ())>10.0 ) lIsVertexPositionGood = false;
        if ( itsSpdVertex->GetNContributors() < 1 ) lIsVertexPositionGood = false; //not okay, no contributor
        //printf("Decisions: v0 decision is %d , vertex decision is %d\n", lIsV0DecisionOK, lIsVertexPositionGood);
        
        Double_t lQuantities[6] = {
            static_cast<Double_t>(esdVZERO->GetMTotV0A()+esdVZERO->GetMTotV0C()),
            static_cast<Double_t>(esdVZERO->GetMTotV0A()),
            static_cast<Double_t>(esdVZERO->GetMTotV0C()),
            static_cast<Double_t>(esdVZERO->GetTriggerChargeA()+esdVZERO->GetTriggerChargeC()),
            static_cast<Double_t>(esdVZERO->GetTriggerChargeA()),
            static_cast<Double_t>(esdVZERO->GetTriggerChargeC())
        };
        Bool_t lDecisions[4] = {lIsMinBias, lIsCentral, lIsSemiCentral, lIsC0V0H};
        
        if( lIsV0DecisionOK && lIsVertexPositionGood ){
            for(Int_t i=0; i<24; i++){
                if(lDecisions[i/6]) fHistMult[i].Fill(lQuantities[i%6]);
            }
        }
    }
    
    //If the histogram is not empty, we send it out every 16th event (to collect some statistics).
    //Depending on the pushback period set for the component, it might not be send out every time, but only after a certain amount of time. (order of 3 minutes)
    //We check whether it was really sent out, and only if so, we reset the histogram.
    //The ZMQ merging component that sits at the end of the chain will receive all histograms from all concurrent VZEROOnlineCalib components, and merge them to the final histogram.
    
    for(Int_t i=0; i<24; i++){
        if ( fHistMult[i].GetEntries() && PushBack(&fHistMult[i], kAliHLTDataTypeHistogram|kAliHLTDataOriginHLT) > 0) fHistMult[i].Reset();
    }
    
    return iResult;
}

// #################################################################################
Int_t AliHLTVZEROOnlineCalibComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
    // see header file for class documentation
    
    Int_t iResult=0;
    return iResult;
}

