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

/** @file    AliHLTTOFRecoComponent.cxx
    @author  Jochen Thaeder <jochen@thaeder.de>
    @brief   TOF reconstruction component
*/

#include "TTree.h"
#include "TMap.h"
#include "TObjString.h"
#include "TDatime.h"
#include "TClonesArray.h"

#include "AliLog.h"
#include "AliRunInfo.h"
#include "AliGRPObject.h"
#include "AliRawReaderMemory.h"
#include "AliGeomManager.h"
#include "AliDAQ.h"

#include "AliTOFRecoParam.h"
#include "AliTOFReconstructor.h"
#include "AliTOFRawStream.h"
#include "AliTOFrawData.h"
#include "AliTOFClusterFinder.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTOFRecoComponent.h"
#include "AliHLTTOFRawCluster.h"

using namespace std;
/*
#define LogError( ... ) { HLTError(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("ERROR", __VA_ARGS__); } }
#define LogInfo( ... ) { HLTInfo(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("INFO", __VA_ARGS__); } }
#define LogInspect( ... ) { HLTDebug(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("INSPECT", __VA_ARGS__); } }
#define LogDebug( ... ) { if (fDebugLevel >= 1) { HLTInfo(__VA_ARGS__); DbgLog("DEBUG", __VA_ARGS__); } }
*/
/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTOFRecoComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTTOFRecoComponent::AliHLTTOFRecoComponent() :
  AliHLTProcessor(),
  fDebugLevel(0),
  fRunInfo(NULL),  
  fTOFRecoParam(NULL),
  fTOFReconstructor(NULL),
  fRawReaderMem(NULL),
  fTOFRawStream(NULL) {
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
AliHLTTOFRecoComponent::~AliHLTTOFRecoComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTTOFRecoComponent::GetComponentID() { 
  // see header file for class documentation
  return "TOFReconstruction";
}

// #################################################################################
void AliHLTTOFRecoComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTOF);
}

// #################################################################################
AliHLTComponentDataType AliHLTTOFRecoComponent::GetOutputDataType() 
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTTOFRecoComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{ 
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back( kAliHLTDataTypeESDContent|kAliHLTDataOriginTOF);
  return tgtList.size();
}

// #################################################################################
void AliHLTTOFRecoComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 3000;
  inputMultiplier = 0.7;
}

// #################################################################################
void AliHLTTOFRecoComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  // see header file for class documentation

  if (!targetMap) return;
  //  targetMap->Add(new TObjString("HLT/ConfigTOF/TOFReconstruction"),
  //		 new TObjString("configuration object"));

  targetMap->Add(new TObjString("GRP/GRP/Data"),
		 new TObjString("GRP object - run information"));
  targetMap->Add(new TObjString("GRP/CTP/CTPtiming"),
		 new TObjString("GRP object - CTP information"));
  targetMap->Add(new TObjString("GRP/CTP/TimeAlign"),
		 new TObjString("GRP object - CTP information"));
  targetMap->Add(new TObjString("GRP/Calib/LHCClockPhase"),
		 new TObjString("GRP object - time calibration"));
  /*
  targetMap->Add(new TObjString("TOF/Calib/Data"),
		 new TObjString("TOF calibration object"));
  targetMap->Add(new TObjString("TOF/Calib/TimeDelays"),
		 new TObjString("TOF calibration object"));
  targetMap->Add(new TObjString("TOF/Calib/TimeSlewing"),
		 new TObjString("TOF calibration object"));
  targetMap->Add(new TObjString("TOF/Trigger/Data"),
		 new TObjString("TOF calibration object"));
  */
  return;
}

// #################################################################################
AliHLTComponent* AliHLTTOFRecoComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTTOFRecoComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTTOFRecoComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation
  //cout<<"\n\n\nTOF Reconstruction Init\n\n\n"<<endl;
 
  Int_t iResult=0;

  // -- Load GeomManager
  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }
  
  // -- Read configuration object : HLT/ConfigTOF/TOFReconstruction
  TString cdbPath="HLT/ConfigTOF/";
  cdbPath+=GetComponentID();
  //iResult=ConfigureFromCDBTObjString(cdbPath);

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

    fRawReaderMem = new AliRawReaderMemory;
    if (!fRawReaderMem) {
      iResult=-ENOMEM;
      break;
    }

    fTOFRawStream = new AliTOFRawStream(fRawReaderMem);
    if (!fTOFRawStream) {
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

    AliCDBEntry* e = AliCDBManager::Instance()->Get("TOF/Calib/RecoParam");
    fTOFRecoParam = (AliTOFRecoParam*)e->GetObject();
    if (!fTOFRecoParam) {
      iResult=-ENOMEM;
      break;
    }  

    fTOFReconstructor = new AliTOFReconstructor;
    if (!fTOFReconstructor) {
      iResult=-ENOMEM;
      break;
    }

    // implement further initialization
  } while (0);

  if (iResult<0) {
    // implement cleanup

    if (fRawReaderMem) 
      delete fRawReaderMem;
    fRawReaderMem = NULL;

    if (fTOFRawStream) 
      delete fTOFRawStream;
    fTOFRawStream = NULL;

    if (fTOFRecoParam)
      delete fTOFRecoParam;
    fTOFRecoParam = NULL;

    if (fTOFReconstructor)
      delete fTOFReconstructor;
    fTOFReconstructor = NULL;

    if (fRunInfo)
      delete fRunInfo;
    fRunInfo = NULL;
  }

  if (iResult>=0) {
    fTOFReconstructor->SetRunInfo(fRunInfo);
    fTOFReconstructor->Init();

    fTOFReconstructor->SetRecoParam(fTOFRecoParam);
  }

  return iResult;
}

// #################################################################################
Int_t AliHLTTOFRecoComponent::ScanConfigurationArgument(Int_t argc, const Char_t** argv) {
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

  if (!argument.CompareTo("-debug")){
    if (++ii >= argc) return -EPROTO;
    argument = argv[ii];
    fDebugLevel = argument.Atoi();
    HLTInfo("debug level set to %d.", fDebugLevel);
    return 2;
  }
  
  
  return 0;
}

// #################################################################################
Int_t AliHLTTOFRecoComponent::DoDeinit() {
  // see header file for class documentation

  if (fRawReaderMem) 
    delete fRawReaderMem;
  fRawReaderMem = NULL;

  if (fTOFRawStream) 
    delete fTOFRawStream;
  fTOFRawStream = NULL;

  if (fTOFRecoParam)
    delete fTOFRecoParam;
  fTOFRecoParam = NULL;
  
  if (fTOFReconstructor)
    delete fTOFReconstructor;
  fTOFReconstructor = NULL;
  
  if (fRunInfo)
    delete fRunInfo;
  fRunInfo = NULL;
  
  return 0;
}

// #################################################################################
Int_t AliHLTTOFRecoComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
				      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/,
				      AliHLTUInt32_t& /*size*/, vector<AliHLTComponentBlockData>& /*outputBlocks*/ ) {
  // see header file for class documentation
  Int_t iResult=0;

  // -- Only use data event
  if (!IsDataEvent()) 
    return 0;
 
  //cout<<"\n\n\nTOF Reconstruction Do Event\n\n\n"<<endl;

  // -- Get TOF raw dat a input block and set up the rawreader
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) ) {
    return 0;
  }

  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;
  
  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {
      iter = blocks+ndx;

      if (iter->fDataType != (kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTOF)) continue;
      if(!fRawReaderMem->AddBuffer((UChar_t*) iter->fPtr, iter->fSize, 1280 + iter->fSpecification)){   
	HLTError("Could not add buffer of data block  %s, 0x%08x to rawreader",
		 DataType2Text(iter->fDataType).c_str(),
		 iter->fSpecification);
      }
      //Printf("ndx = %d, fSpecification = %d", ndx, iter->fSpecification);
      //cout<<"TOF input block found with size = "<< iter->fSize << "\n";
  }
  
  //fTOFRawStream->ReadEvent();
  Int_t indexDDL = 0;
  TClonesArray * clonesRawData;
  Int_t dummy = -1;
  Int_t detectorIndex[5]={-1,-1,-1,-1,-1};
  Int_t parTOF[7];

  // This part is taken from AliTOFClusterFinder::Digits2RecPoints
  const Int_t kDDL = AliDAQ::NumberOfDdls("TOF");

  Int_t inholes = 0;
  
  for (Int_t indexDDL = 0; indexDDL < kDDL; indexDDL++) {
    fRawReaderMem->Reset();
    fTOFRawStream->LoadRawDataBuffersV2(indexDDL,fDebugLevel);
    clonesRawData = (TClonesArray*)fTOFRawStream->GetRawData();
    for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {
      AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

      if (tofRawDatum->GetTOF()==-1) continue;

      if (1) {
	//if (fDebugLevel==2) {
	Printf("indexDDL = %d", indexDDL);
	Printf("TRD = %d", tofRawDatum->GetTRM());
	Printf("TRM chain = %d", tofRawDatum->GetTRMchain());
	Printf("TDC = %d", tofRawDatum->GetTDC());
	Printf("TDC channle = %d", tofRawDatum->GetTDCchannel());
      }

      for (Int_t aa=0; aa<5; aa++) detectorIndex[aa] = -1;

      fTOFRawStream->EquipmentId2VolumeId(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
					 tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);
      dummy = detectorIndex[3];
      detectorIndex[3] = detectorIndex[4];
      detectorIndex[4] = dummy;

      if (1) {
	//if (fVerbose==2) {
	Printf("index[0] = %d", detectorIndex[0]);
	Printf("index[1] = %d", detectorIndex[1]);
	Printf("index[2] = %d", detectorIndex[2]);
	Printf("index[3] = %d", detectorIndex[3]);
	Printf("index[4] = %d", detectorIndex[4]);
      }

      /* check valid index */
      if (detectorIndex[0]==-1||detectorIndex[1]==-1||detectorIndex[2]==-1||detectorIndex[3]==-1||detectorIndex[4]==-1) continue;
      // Do not reconstruct anything in the holes
      if (detectorIndex[0]==13 || detectorIndex[0]==14 || detectorIndex[0]==15 ) { // sectors with holes
	if (detectorIndex[1]==2) { // plate with holes
	  inholes++;
	  continue;
	}
      }
      /*
      int tdc = tofRawDatum->GetTOF(); //TDC
      int tot = tofRawDatum->GetTOT(); // TOT
      int adc = tofRawDatum->GetTOT(); //ADC==TOF
      int tdcnd = 0;//raw data: no track of undecalib sim time
      int tdcraw = tofRawDatum->GetTOF(); // RAW time
      int deltabc = tofRawDatum->GetDeltaBC(); // deltaBC
      int l0l1latency = tofRawDatum->GetL0L1Latency(); // L0-L1 latency
      */
      parTOF[0] = tofRawDatum->GetTOF(); // TDC
      parTOF[1] = tofRawDatum->GetTOT(); // TOT
      parTOF[2] = tofRawDatum->GetTOT(); // raw data have ADC=TOT
      parTOF[3] = 0; //raw data: no track of the undecalib sim time
      parTOF[4] = tofRawDatum->GetTOF(); // Raw time == TDC
      parTOF[5] = tofRawDatum->GetDeltaBC(); // deltaBC
      parTOF[6] = tofRawDatum->GetL0L1Latency(); // L0-L1 latency
      Double_t posClus[3];
      Double_t covClus[6];
      UShort_t volIdClus = AliTOFClusterFinder::GetClusterVolIndex(detectorIndex);
      Int_t lab[3]={-1,-1,-1};
      Bool_t status=kTRUE;
      AliTOFClusterFinder::GetClusterPars(detectorIndex,posClus,covClus);

      Printf("tdc = %d, tot = %d, adc = %d, tdcnd = %d, tdcraw = %d, deltabc = %d, l0l1latency = %d, idx = %d, status = %d", parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4], parTOF[5], parTOF[6], volIdClus, (Int_t)status);

      AliHLTTOFRawCluster *tofCluster = new AliHLTTOFRawCluster(volIdClus, posClus[0], posClus[1], posClus[2], covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5], lab, detectorIndex, parTOF, status, -1);

      //      AliHLTTOFRawCluster *c = new AliHLTTOFRawCluster(idx, detectorIndex, quality, r, phi, tdc, tot, adc, tdcnd, tdcraw, status, deltabc, l0l1latency);
    


    //  fTOFRawStream->DecodeV2(fDebugLevel);
  
  //cout<<"TOF input block found"<<endl;
    }
  }
  
  TTree *digitsTree = new TTree("D", "Digits Tree");
  if (!digitsTree) {
    cout<<"No digit tree created"<<endl;
    iResult=-ENOMEM;
  }

  if (iResult >= 0) {
    //cout<<"ok 1"<<endl;
    // -- Set TOF EquipmentID
    //   fRawReader->SetEquipmentID(3584);
  
    // -- 1. step TOF reconstruction
    fTOFReconstructor->ConvertDigits(fRawReaderMem, digitsTree);

    // -- 2. step TOF reconstruction -- fill AliESDTOF object
    // for TOF, this does nothing 
    //fTOFReconstructor->FillESD(digitsTree, NULL, NULL);
    /*
    AliESDTOF *esdTOF = fTOFReconstructor->GetESDTOF();
    AliESDTOFfriend *esdTOFfriend = fTOFReconstructor->GetESDTOFfriend();
    */
    
    // Send info every 10 s
    const TDatime time;    
    static UInt_t lastTime=0;
    /*
    if (time.Get()-lastTime>10) {
      lastTime=time.Get();
      HLTInfo("TOF Multiplicity A %f - C %f", esdTOF->GetMTotTOFA(), esdTOF->GetMTotTOFA() );
    }
    */
    // -- Send AliESDTOF & friend object
    //    PushBack(esdTOF, kAliHLTDataTypeESDContent|kAliHLTDataOriginTOF,0);
    // PushBack( esdTOFfriend, kAliHLTDataTypeESDFriendContent|kAliHLTDataOriginTOF,0);
  }
  
  // -- Clean up
  delete digitsTree;
  fRawReaderMem->ClearBuffers();   

  return iResult;
}

// #################################################################################
Int_t AliHLTTOFRecoComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
  // see header file for class documentation

  Int_t iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigTOF/";
    cdbPath+=GetComponentID();
  }

  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}

// #################################################################################
Int_t AliHLTTOFRecoComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  return 0;
}
