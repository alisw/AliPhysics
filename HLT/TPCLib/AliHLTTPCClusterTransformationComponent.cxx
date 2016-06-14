
//***************************************************************************
//* This file is property of and copyright by the ALICE HLT Project         *
//* ALICE Experiment at CERN, All rights reserved.                          *
//*                                                                         *
//* Primary Authors: Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de *
//*                  for The ALICE HLT Project.                             *
//*                                                                         *
//* Permission to use, copy, modify and distribute this software and its    *
//* documentation strictly for non-commercial purposes is hereby granted    *
//* without fee, provided that the above copyright notice appears in all    *
//* copies and that both the copyright notice and this permission notice    *
//* appear in the supporting documentation. The authors make no claims      *
//* about the suitability of this software for any purpose. It is           *
//* provided "as is" without express or implied warranty.                   *
//***************************************************************************

/** @file   AliHLTTPCClusterTransformationComponent.cxx
    @author Sergey Gorbunov
    @date   
    @brief 
*/

#include "AliHLTTPCClusterTransformationComponent.h"
#include "AliHLTTPCClusterTransformation.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTTPCFastTransformObject.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliDAQ.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTPCcalibDB.h"

#include "TMath.h"
#include "TObjString.h" 
#include <cstdlib>
#include <cerrno>
#include <sys/time.h>

using namespace std;

ClassImp(AliHLTTPCClusterTransformationComponent) //ROOT macro for the implementation of ROOT specific class methods

const char* AliHLTTPCClusterTransformationComponent::fgkOCDBEntryClusterTransformation="HLT/ConfigTPC/TPCClusterTransformation";

AliHLTTPCClusterTransformation AliHLTTPCClusterTransformationComponent::fgTransform;
Bool_t AliHLTTPCClusterTransformationComponent::fgTimeInitialisedFromEvent = 0;

AliHLTTPCClusterTransformationComponent::AliHLTTPCClusterTransformationComponent()
:
fOfflineMode(0),
fInitializeByObjectInDoEvent(0),
fInitialized(0),
fTPCPresent(0),
fDataId(kFALSE),
fBenchmark("ClusterTransformation")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt  

  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
}

AliHLTTPCClusterTransformationComponent::~AliHLTTPCClusterTransformationComponent()
{ 
  // destructor
}

const char* AliHLTTPCClusterTransformationComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCClusterTransformation";
}

void AliHLTTPCClusterTransformationComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
  // see header file for class documentation

  list.clear(); 
  list.push_back( AliHLTTPCDefinitions::fgkRawClustersDataType  | kAliHLTDataOriginTPC  );
  list.push_back( AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC );
  list.push_back( AliHLTTPCDefinitions::fgkTPCFastTransformDataObjectDataType | kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCClusterTransformationComponent::GetOutputDataType() { 
  // see header file for class documentation

  return kAliHLTMultipleDataType;
}

int AliHLTTPCClusterTransformationComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList) { 
  // see header file for class documentation

  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkClustersDataType| kAliHLTDataOriginTPC);
  tgtList.push_back(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC );
  return tgtList.size();
}

void AliHLTTPCClusterTransformationComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
  // see header file for class documentation
  constBase = 0;
  inputMultiplier = 2.5;
}

AliHLTComponent* AliHLTTPCClusterTransformationComponent::Spawn() { 
  // see header file for class documentation

  return new AliHLTTPCClusterTransformationComponent();
}
	
int AliHLTTPCClusterTransformationComponent::DoInit( int argc, const char** argv ) 
{ 
  // see header file for class documentation
  
  int iResult=0;
  //!! iResult = ConfigureFromCDBTObjString(fgkOCDBEntryClusterTransformation);

  AliGRPManager mgr;
  mgr.ReadGRPEntry();
  fTPCPresent = ((mgr.GetGRPData()->GetDetectorMask() & AliDAQ::kTPC) != 0);

  if (!fTPCPresent) return(iResult);

  fOfflineMode = 0;

  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);

  AliTPCcalibDB *calib=AliTPCcalibDB::Instance();  
  if(!calib){
    HLTError("AliTPCcalibDB does not exist");
    return -ENOENT;
  }
  calib->SetRun(GetRunNo());
  calib->UpdateRunInformations(GetRunNo());
  
  if( !fgTransform.IsInitialised() ){
    TStopwatch timer;
    timer.Start();
    int err = 0;
    if ( fInitializeByObjectInDoEvent == 1 ) {
          HLTInfo( "Cluster Transformation will initialize on the fly in DoEvent loop via FastTransformation Data Object, skipping initialization." );
    }
    else if( fOfflineMode ) {
      err = fgTransform.Init( GetBz(), GetTimeStamp() );
	  fInitialized = true;
    } else {
       const char* defaultNotify = "";
       const char* cdbEntry = "HLT/ConfigTPC/TPCFastTransform";
       defaultNotify = " (default)";
       const char* chainId = 0;
       
       HLTInfo( "configure from entry \"%s\"%s, chain id %s", cdbEntry, defaultNotify, ( chainId != NULL && chainId[0] != 0 ) ? chainId : "<none>" );
       AliCDBEntry *pEntry = AliCDBManager::Instance()->Get( cdbEntry );//,GetRunNo());
       if ( !pEntry ) {
	 HLTError( "cannot fetch object \"%s\" from CDB", cdbEntry );
	 return -EINVAL;
       }
       const AliHLTTPCFastTransformObject *configObj = dynamic_cast<const AliHLTTPCFastTransformObject *>( pEntry->GetObject() );

       if ( !configObj ) {
	 HLTError( "configuration object \"%s\" has wrong type, required TObjString", cdbEntry );
	 return -EINVAL;
       }

       HLTInfo( "received configuration object." );
       fgTransform.Init( *configObj );
	   fInitialized = true;
    }
    timer.Stop();
    HLTInfo("Initialization time: %f / %f", timer.CpuTime(), timer.RealTime());
    if( err!=0 ){
      HLTError(Form("Cannot retrieve offline transform from AliTPCcalibDB, AliHLTTPCClusterTransformation returns %d",err));
      return -ENOENT;
    }
  }

  fDataId = kFALSE;

  return iResult;
} // end DoInit()

int AliHLTTPCClusterTransformationComponent::DoDeinit() { 
  // see header file for class documentation   
  if (fInitialized) fgTransform.DeInit();
  fInitialized = false;
  return 0;
}

int AliHLTTPCClusterTransformationComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
  // see header file for class documentation
  fDataId = kFALSE;
  return 0;//!! ConfigureFromCDBTObjString(fgkOCDBEntryClusterTransformation);
}

int AliHLTTPCClusterTransformationComponent::ScanConfigurationArgument(int argc, const char** argv){

  // see header file for class documentation

  if (argc<=0) return 0;
  int iRet = 0;
  for( int i=0; i<argc; i++ ){
    TString argument=argv[i];  
    if (argument.CompareTo("-change-dataId")==0){
      HLTDebug("Change data ID received.");
      fDataId = kTRUE;
      iRet++;
    } else if (argument.CompareTo("-offline-mode")==0){
      fOfflineMode = 1;
      HLTDebug("Offline mode set.");
      iRet++;
    } else if (argument.CompareTo("-initialize-on-the-fly")==0){
      fInitializeByObjectInDoEvent = 1;
      HLTDebug("Initialize on the fly mode set.");
      iRet++;
    } else if (argument.CompareTo("-update-object-on-the-fly")==0){
      fInitializeByObjectInDoEvent = 2;
      HLTDebug("Initialize object at startup and update on the fly mode set.");
      iRet++;
    } else {
      iRet = -EINVAL;
      HLTError("Unknown argument %s",argv[i]);     
    }
  } 
  return iRet;
}


int AliHLTTPCClusterTransformationComponent::DoEvent(const AliHLTComponentEventData& evtData, 
					          const AliHLTComponentBlockData* blocks, 
					          AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					          AliHLTUInt32_t& size, 
					          vector<AliHLTComponentBlockData>& outputBlocks ){
  // see header file for class documentation
 
  UInt_t maxOutSize = size;
  size = 0;
  int iResult = 0;

  if (!fTPCPresent) return iResult;

  if (fInitializeByObjectInDoEvent)
  {
	//Check first whether there is a new FastTransformation object
	for ( const TObject *iter = GetFirstInputObject(AliHLTTPCDefinitions::fgkTPCFastTransformDataObjectDataType); iter != NULL; iter = GetNextInputObject() )
	{
		const AliHLTTPCFastTransformObject *configObj = dynamic_cast<const AliHLTTPCFastTransformObject *>(const_cast<TObject*>( iter ) );

		if (!configObj)
		{
			HLTError( "Error getting configuration object" );
			return -EINVAL;
		}
		if (fInitialized)
		{
			HLTInfo("Received updated cluster transformation map with new calibration");
			fgTransform.DeInit();
		}
		else
		{
			HLTInfo("Received initial cluster transformation map");
		}

		HLTInfo( "received configuration object." );
		if (fgTransform.Init( *configObj ))
		{
			HLTError("Failed on-the-fly-initialization of transformation map. Error: %s", fgTransform.GetLastError());
			return(-1);
		}
		fInitialized = true;
		break;
	}
  }

  if(!IsDataEvent()) return 0;
  
  if( !fgTransform.IsInitialised() ){
    HLTError(" TPC Transformation is not initialised ");
    return -ENOENT;    
  }

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  for( unsigned long ndx=0; ndx<evtData.fBlockCnt; ndx++ ){
    
    const AliHLTComponentBlockData *iter   = blocks+ndx;
    
    fBenchmark.AddInput(iter->fSize);
    
    HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s or $s",
	     evtData.fEventID, evtData.fEventID, 
	     DataType2Text( iter->fDataType).c_str(), 
	     DataType2Text(AliHLTTPCDefinitions::fgkRawClustersDataType).c_str(), DataType2Text(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo).c_str());                       
 
    if(iter->fDataType == (AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC) ){
      // simply forward MC labels
      
      if( size+iter->fSize > maxOutSize ){
	HLTWarning( "Output buffer (%db) is too small, required %db", maxOutSize, size+iter->fSize);
	iResult  = -ENOSPC;
	break;
      }

      memcpy( outputPtr, iter->fPtr, iter->fSize );
      
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = size;
      bd.fSize = iter->fSize;
      bd.fSpecification = iter->fSpecification;     
      bd.fDataType = iter->fDataType;
      outputBlocks.push_back( bd );     
      fBenchmark.AddOutput(bd.fSize);    
      size   += bd.fSize;
      outputPtr += bd.fSize;
      continue;
    }

    if(iter->fDataType != (AliHLTTPCDefinitions::fgkRawClustersDataType  | kAliHLTDataOriginTPC)) continue;                        
        
    UInt_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter); 
    UInt_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);

    float padpitch=1.0;
    if ((int)minPartition<AliHLTTPCGeometry::GetNRowLow())
      padpitch=AliHLTTPCGeometry::GetPadPitchWidthLow();
    else
      padpitch=AliHLTTPCGeometry::GetPadPitchWidthUp();
    float zwidth=AliHLTTPCGeometry::GetZWidth();

    fBenchmark.SetName(Form("ClusterTransform slice %d patch %d",minSlice,minPartition));

    HLTDebug("minSlice: %d, minPartition: %d", minSlice, minPartition);
    
    AliHLTTPCRawClusterData* rawClusters = (AliHLTTPCRawClusterData*)(iter->fPtr);
    if( !rawClusters ) continue;

    AliHLTTPCClusterData* outPtr  = (AliHLTTPCClusterData*)outputPtr;
    outPtr->fSpacePointCnt=0;

    long maxPoints = ((long)maxOutSize-size-sizeof(AliHLTTPCClusterData))/sizeof(AliHLTTPCSpacePointData);
  
    if( rawClusters->fCount > maxPoints ){
      HLTWarning("No more space to add clusters, exiting!");
      iResult  = -ENOSPC;
      break;
    }
  
    for( UInt_t icl=0; icl<rawClusters->fCount; icl++){
      
      const AliHLTTPCRawCluster &cl = rawClusters->fClusters[icl];

      AliHLTTPCSpacePointData& c=outPtr->fSpacePoints[outPtr->fSpacePointCnt];
      int padrow=cl.GetPadRow();
      if (padrow<0) {
	// something wrong here, padrow is stored in the cluster header
	// word which has bit pattern 0x3 in bits bit 30 and 31 which was
	// not recognized
	ALIHLTERRORGUARD(1, "can not read cluster header word");
	break;
      }
      padrow+=AliHLTTPCGeometry::GetFirstRow(minPartition);
      AliHLTUInt32_t charge=cl.GetCharge();

      float pad=cl.GetPad();
      float time=cl.GetTime();
      float sigmaY2=cl.GetSigmaPad2();
      float sigmaZ2=cl.GetSigmaTime2();
      sigmaY2*=padpitch*padpitch;
      sigmaZ2*=zwidth*zwidth;
      c.SetPad( cl.GetPad() );
      c.SetTime( cl.GetTime() );
      c.SetPadRow(padrow);
      c.SetCharge(charge);
      c.SetSigmaY2(sigmaY2);
      c.SetSigmaZ2(sigmaZ2);
      c.SetQMax(cl.GetQMax());
        Float_t xyz[3];
      int err = fgTransform.Transform( minSlice, padrow, pad, time, xyz );	 
      if( err!=0 ){
	HLTWarning(Form("Cannot transform the cluster, AliHLTTPCClusterTransformation returns error %d, %s",err, fgTransform.GetLastError()));
	continue;
      }
      c.SetX(xyz[0]);
      c.SetY(xyz[1]);
      c.SetZ(xyz[2]);

      // set the cluster ID so that the cluster dump printout is the same for FCF and SCF
      c.SetID( minSlice, minPartition, outPtr->fSpacePointCnt );
	 
      HLTDebug("Cluster number %d: %f, Y: %f, Z: %f, charge: %d \n", outPtr->fSpacePointCnt, c.GetX(), c.GetY(), c.GetZ(), (UInt_t)c.GetCharge());
	 
      outPtr->fSpacePointCnt++; 
    } // end of loop over clusters
      
    HLTDebug("Number of found clusters: %d", outPtr->fSpacePointCnt);
     
    UInt_t mysize = sizeof(AliHLTTPCClusterData) + sizeof(AliHLTTPCSpacePointData)*outPtr->fSpacePointCnt;
    
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset = size;
    bd.fSize = mysize;
    bd.fSpecification = iter->fSpecification;     
    if(fDataId==kFALSE) bd.fDataType = AliHLTTPCDefinitions::fgkClustersDataType;
    else                bd.fDataType = AliHLTTPCDefinitions::fgkAlterClustersDataType;
    
    //HLTDebug("datatype: %s", DataType2Text(bd.fDataType).c_str());
     
    outputBlocks.push_back( bd );
     
    fBenchmark.AddOutput(bd.fSize);    
    size   += mysize;
    outputPtr += mysize; 
 
  } // end of loop over data blocks  
  
  fBenchmark.Stop(0);
  HLTInfo(fBenchmark.GetStatistics());
  
  return iResult;
} // end DoEvent()



void AliHLTTPCClusterTransformationComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description needed for the particular component
  if (!targetMap) return;
  
  // OCDB entries for component arguments

  //!! targetMap->Add(new TObjString(fgkOCDBEntryClusterTransformation), new TObjString("component argument for the charge threshold"));
  
  // OCDB entries to be fetched by the TAXI (access via the AliTPCcalibDB class)
  targetMap->Add(new TObjString("TPC/Calib/Parameters"),    new TObjString("unknown content"));
  targetMap->Add(new TObjString("TPC/Calib/TimeDrift"),     new TObjString("drift velocity calibration"));
  targetMap->Add(new TObjString("TPC/Calib/TimeGain"),     new TObjString("time gain  calibration"));
  targetMap->Add(new TObjString("TPC/Calib/Temperature"),   new TObjString("temperature map"));
  targetMap->Add(new TObjString("TPC/Calib/PadGainFactor"), new TObjString("gain factor pad by pad"));
  targetMap->Add(new TObjString("TPC/Calib/ClusterParam"),  new TObjString("cluster parameters"));
  targetMap->Add(new TObjString("TPC/Calib/Correction"),  new TObjString("coreection"));
  targetMap->Add(new TObjString("TPC/Calib/RecoParam"),  new TObjString("reconstruction parameters"));
 
  // OCDB entries needed to be fetched by the Pendolino
  targetMap->Add(new TObjString("TPC/Calib/AltroConfig"), new TObjString("contains the altro config, e.g. info about the L0 trigger timing"));
  targetMap->Add(new TObjString("GRP/CTP/CTPtiming"),     new TObjString("content used in the cluster coordinate transformation in relation to the L0 trigger timing"));

  // OCDB entries necessary for replaying data on the HLT cluster
  targetMap->Add(new TObjString("GRP/GRP/Data"), new TObjString("contains magnetic field info"));  
 
  // OCDB entries needed to suppress fatals/errors/warnings during reconstruction
  targetMap->Add(new TObjString("TPC/Calib/Distortion"),  new TObjString("distortion map"));
  targetMap->Add(new TObjString("TPC/Calib/GainFactorDedx"), new TObjString("gain factor dedx"));
  targetMap->Add(new TObjString("TPC/Calib/PadTime0"),    new TObjString("time0 offset pad by pad"));
  targetMap->Add(new TObjString("TPC/Calib/PadNoise"),    new TObjString("pad noise values"));
  targetMap->Add(new TObjString("TPC/Calib/Pedestals"),   new TObjString("pedestal info"));
  targetMap->Add(new TObjString("TPC/Calib/Pulser"),      new TObjString("pulser info"));
  targetMap->Add(new TObjString("TPC/Calib/CE"),          new TObjString("CE laser calibration result"));
  targetMap->Add(new TObjString("TPC/Calib/Raw"),         new TObjString("unknown content"));
  targetMap->Add(new TObjString("TPC/Calib/QA"),          new TObjString("not important"));
  targetMap->Add(new TObjString("TPC/Calib/Mapping"),     new TObjString("unknown content"));
  targetMap->Add(new TObjString("TPC/Calib/Goofie"),      new TObjString("Goofie values, not used at the moment (05.03.2010)"));
  targetMap->Add(new TObjString("TPC/Calib/HighVoltage"), new TObjString("high voltage values, not used"));
  targetMap->Add(new TObjString("TPC/Calib/Ref"),         new TObjString("unknown content"));
}
