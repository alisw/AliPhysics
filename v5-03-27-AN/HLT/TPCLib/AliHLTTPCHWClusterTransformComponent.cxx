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

/** @file   AliHLTTPCHWClusterTransformComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief 
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCHWClusterTransformComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliRawDataHeader.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCHWCFEmulator.h"
#include "AliHLTTPCHWCFData.h"
#include "AliHLTErrorGuard.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTPCcalibDB.h"

#include "TMath.h"
#include "TObjString.h" 
#include <cstdlib>
#include <cerrno>
#include <sys/time.h>

ClassImp(AliHLTTPCHWClusterTransformComponent) //ROOT macro for the implementation of ROOT specific class methods

const char* AliHLTTPCHWClusterTransformComponent::fgkOCDBEntryHWTransform="HLT/ConfigTPC/TPCHWClusterTransform";

AliHLTTPCHWClusterTransformComponent::AliHLTTPCHWClusterTransformComponent()
:
fDataId(kFALSE),
fTransform(),
fPublishRawClusters(kFALSE),
fpDecoder(NULL),
fBenchmark("HWClusterTransform")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt  

  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
}

AliHLTTPCHWClusterTransformComponent::~AliHLTTPCHWClusterTransformComponent()
{ 
  // destructor
  if (!fpDecoder) delete fpDecoder;
  fpDecoder=NULL;
}

const char* AliHLTTPCHWClusterTransformComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCHWClusterTransform";
}

void AliHLTTPCHWClusterTransformComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
  // see header file for class documentation

  list.clear(); 
  list.push_back( AliHLTTPCDefinitions::fgkHWClustersDataType );
}

AliHLTComponentDataType AliHLTTPCHWClusterTransformComponent::GetOutputDataType() { 
  // see header file for class documentation

  return AliHLTTPCDefinitions::fgkClustersDataType;
}

int AliHLTTPCHWClusterTransformComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList) { 
  // see header file for class documentation

  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkClustersDataType| kAliHLTDataOriginTPC);
  tgtList.push_back(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC );
  tgtList.push_back(AliHLTTPCDefinitions::fgkRawClustersDataType  | kAliHLTDataOriginTPC );
  return tgtList.size();
}

void AliHLTTPCHWClusterTransformComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
  // see header file for class documentation
  constBase = 0;
  inputMultiplier = 3.0;
}

AliHLTComponent* AliHLTTPCHWClusterTransformComponent::Spawn() { 
  // see header file for class documentation

  return new AliHLTTPCHWClusterTransformComponent();
}
	
int AliHLTTPCHWClusterTransformComponent::DoInit( int argc, const char** argv ) { 
// see header file for class documentation
  
  AliTPCcalibDB *calib=AliTPCcalibDB::Instance();  
  if(!calib){
    HLTError("AliTPCcalibDB does not exist");
    return -ENOENT;
  }
  calib->SetRun(GetRunNo());
  calib->UpdateRunInformations(GetRunNo());

  int err = fTransform.Init( GetBz(), GetTimeStamp() );

  if( err!=0 ){
    HLTError("Cannot retrieve offline transform from AliTPCcalibDB");
    return -ENOENT;
  }

  int iResult=0;
  iResult = ConfigureFromCDBTObjString(fgkOCDBEntryHWTransform);

  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);

  if (iResult>=0) {
    fpDecoder=new AliHLTTPCHWCFData;
    if (!fpDecoder) iResult=-ENOMEM;
  }
  
  return iResult;
} // end DoInit()

int AliHLTTPCHWClusterTransformComponent::DoDeinit() { 
  // see header file for class documentation   
  if (!fpDecoder) delete fpDecoder;
  fpDecoder=NULL;

  return 0;
}

int AliHLTTPCHWClusterTransformComponent::DoEvent(const AliHLTComponentEventData& evtData, 
					          const AliHLTComponentBlockData* blocks, 
					          AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					          AliHLTUInt32_t& size, 
					          vector<AliHLTComponentBlockData>& outputBlocks ){
  // see header file for class documentation
   
  UInt_t maxOutSize = size;
  size = 0;
  int iResult = 0;
  if(!IsDataEvent()) return 0;

  if (!fpDecoder) return -ENODEV;

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  fTransform.SetCurrentTimeStamp( GetTimeStamp() );
  
  for( unsigned long ndx=0; ndx<evtData.fBlockCnt; ndx++ ){
     
    const AliHLTComponentBlockData *iter   = blocks+ndx;
    
    fBenchmark.AddInput(iter->fSize);
    
    HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	     evtData.fEventID, evtData.fEventID, 
	     DataType2Text( iter->fDataType).c_str(), 
	     DataType2Text(AliHLTTPCDefinitions::fgkHWClustersDataType).c_str());                       
 
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

    if(iter->fDataType != (AliHLTTPCDefinitions::fgkHWClustersDataType | kAliHLTDataOriginTPC)) continue;                        
        
    UInt_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter); 
    UInt_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);
    //UInt_t maxSlice     = AliHLTTPCDefinitions::GetMaxSliceNr(*iter); 
    //UInt_t maxPartition = AliHLTTPCDefinitions::GetMaxPatchNr(*iter);
    float padpitch=1.0;
    if ((int)minPartition<AliHLTTPCTransform::GetNRowLow())
      padpitch=AliHLTTPCTransform::GetPadPitchWidthLow();
    else
      padpitch=AliHLTTPCTransform::GetPadPitchWidthUp();
    float zwidth=AliHLTTPCTransform::GetZWidth();

    fBenchmark.SetName(Form("HWClusterTransform slice %d patch %d",minSlice,minPartition));

    HLTDebug("minSlice: %d, minPartition: %d", minSlice, minPartition);
    
    AliHLTTPCClusterData* outPtr  = (AliHLTTPCClusterData*)outputPtr;
    outPtr->fSpacePointCnt=0;

    long maxPoints = ((long)maxOutSize-size-sizeof(AliHLTTPCClusterData))/sizeof(AliHLTTPCSpacePointData);

    AliHLTUInt32_t *buffer;     
    buffer = (AliHLTUInt32_t*)iter->fPtr;  
     
     // skip the first 8 32-bit CDH words
     buffer += 8;
     UInt_t bufferSize32 = ((Int_t)iter->fSize - sizeof(AliRawDataHeader) )/sizeof(AliHLTUInt32_t);

     if (fpDecoder->Init(reinterpret_cast<AliHLTUInt8_t*>(buffer), bufferSize32*sizeof(AliHLTUInt32_t))>=0 && fpDecoder->CheckVersion()>=0) {
       for (AliHLTTPCHWCFData::iterator cl=fpDecoder->begin(); cl!=fpDecoder->end(); ++cl) {
	 if(outPtr->fSpacePointCnt>=maxPoints){
	   HLTWarning("No more space to add clusters, exiting!");
	   iResult  = -ENOSPC;
	   break;
	 }

	 AliHLTTPCSpacePointData& c=outPtr->fSpacePoints[outPtr->fSpacePointCnt];
	 int padrow=cl.GetPadRow();
	 if (padrow<0) {
	   // something wrong here, padrow is stored in the cluster header
	   // word which has bit pattern 0x3 in bits bit 30 and 31 which was
	   // not recognized
	   ALIHLTERRORGUARD(1, "can not read cluster header word");
	   break;
	 }
	 padrow+=AliHLTTPCTransform::GetFirstRow(minPartition);
	 AliHLTUInt32_t charge=cl.GetCharge();

	 float pad=cl.GetPad();
	 float time=cl.GetTime();
	 float sigmaY2=cl.GetSigmaY2();
	 float sigmaZ2=cl.GetSigmaZ2();
	 sigmaY2*=padpitch*padpitch;
	 sigmaZ2*=zwidth*zwidth;
	 c.SetPadRow(padrow);
	 c.SetCharge(charge);
	 c.SetSigmaY2(sigmaY2);
	 c.SetSigmaZ2(sigmaZ2);
	 c.SetQMax(cl.GetQMax());

	 Float_t xyz[3];
	 fTransform.Transform( minSlice, padrow, pad, time, xyz );
	 c.SetX(xyz[0]);
	 c.SetY(xyz[1]);
	 c.SetZ(xyz[2]);

	 // set the cluster ID so that the cluster dump printout is the same for FCF and SCF
	 c.SetID( minSlice, minPartition, outPtr->fSpacePointCnt );
	 
	 HLTDebug("Cluster number %d: %f, Y: %f, Z: %f, charge: %d \n", outPtr->fSpacePointCnt, cluster.fX, cluster.fY, cluster.fZ, (UInt_t)cluster.fCharge);
	 
	 outPtr->fSpacePointCnt++; 
       } // end of loop over clusters
     }     
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
  
     if (fPublishRawClusters) {

       long maxRawClusters = ((long)maxOutSize-size-sizeof(AliHLTTPCRawClusterData))/sizeof(AliHLTTPCRawCluster);
       
       if( maxRawClusters<=0 ) {
	 HLTWarning("No more space to add raw clusters, exiting!");
	 iResult  = -ENOSPC;
       } else {       

	 // copy raw cluster data from input
	 
	 AliHLTTPCRawClusterData* outputRaw= (AliHLTTPCRawClusterData*)(outputPtr);
       
	 outputRaw->fVersion = 0;
	 outputRaw->fCount = 0;

	 // check if there are clusters available, if not the format might
	 // not even been decoded at that moment 
	 if (fpDecoder->GetNumberOfClusters()>0) {
	 for (AliHLTTPCHWCFData::iterator cl=fpDecoder->begin(); cl!=fpDecoder->end(); ++cl) {
	   if(outputRaw->fCount>=maxRawClusters){
	     HLTWarning("No more space to add clusters, exiting!");
	     iResult  = -ENOSPC;
	     break;
	   }
	   AliHLTTPCRawCluster &c = outputRaw->fClusters[outputRaw->fCount];
	   int padrow=cl.GetPadRow();
	   if (padrow<0) {
	     // something wrong here, padrow is stored in the cluster header
	     // word which has bit pattern 0x3 in bits bit 30 and 31 which was
	     // not recognized
	     break;
	   }
	   padrow+=AliHLTTPCTransform::GetFirstRow(minPartition);
	   AliHLTUInt32_t charge= cl.GetCharge();

	   float pad =cl.GetPad();
	   float time =cl.GetTime();
	   float sigmaP2=cl.GetSigmaY2();
	   float sigmaT2=cl.GetSigmaZ2();
	   c.SetPadRow(padrow);
	   c.SetCharge(charge);
	   c.SetPad(pad);  
	   c.SetTime(time);
	   c.SetSigmaY2(sigmaP2);
	   c.SetSigmaZ2(sigmaT2);
	   c.SetQMax(cl.GetQMax());

	   // store cluster and continue
	   outputRaw->fCount++;
	 }
	 }

	 // fill into HLT output data
	 AliHLTComponentBlockData bdRawClusters;
	 FillBlockData( bdRawClusters );
	 bdRawClusters.fOffset = size;
	 bdRawClusters.fSize = sizeof(AliHLTTPCRawClusterData)+outputRaw->fCount*sizeof(AliHLTTPCRawCluster);
	 bdRawClusters.fSpecification = iter->fSpecification;
	 bdRawClusters.fDataType = AliHLTTPCDefinitions::fgkRawClustersDataType | kAliHLTDataOriginTPC;
	 outputBlocks.push_back( bdRawClusters );
	 fBenchmark.AddOutput(bdRawClusters.fSize);
	 size   += bdRawClusters.fSize;
	 outputPtr += bdRawClusters.fSize;
       }
     }
  } // end of loop over data blocks  
  
  fBenchmark.Stop(0);
  HLTInfo(fBenchmark.GetStatistics());
  
  return iResult;
} // end DoEvent()

int AliHLTTPCHWClusterTransformComponent::ScanConfigurationArgument(int argc, const char** argv){

  // see header file for class documentation

  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  if (argument.CompareTo("-solenoidBz")==0){
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();
    if(!calib){
      HLTError("CalibDB instance cannot be created.");
      return 0;
    }
    Float_t magneticField = argument.Atof();
    calib->SetExBField(magneticField);
    HLTInfo("SolenoidBz is set to %f in the calibDB",magneticField);
    return 2;
  }

  if (argument.CompareTo("-change-dataId")==0){
    HLTDebug("Change data ID received.");
    fDataId = kTRUE;
    return 1;
  }
  
  if (argument.CompareTo("-charge-threshold")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];    
    HLTInfo("The argument -charge-threshold is deprecated.");
    return 2;
  }    

  if (argument.CompareTo("-publish-raw")==0) {
    fPublishRawClusters=kTRUE;
    return 1;
  }  
  
  // unknown argument
  return -EINVAL;
}

int AliHLTTPCHWClusterTransformComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
  // see header file for class documentation
  return ConfigureFromCDBTObjString(fgkOCDBEntryHWTransform);
}

void AliHLTTPCHWClusterTransformComponent::PrintDebug(AliHLTUInt32_t *buffer, Int_t size){
// see header file for class documentation 

  HLTInfo("The size is: %d", size);
  for(Int_t n32bit=0; n32bit<size; n32bit++){
    
    AliHLTUInt8_t *wordPtr = reinterpret_cast<AliHLTUInt8_t*>(&buffer[n32bit]);
    //    cout << "word ptr initialized"<<endl;
    for(Int_t w=3;w>=0;w--){
      //     cout <<"accessing word"<<endl;
      AliHLTUInt8_t word = wordPtr[w];
      //     cout<< "word was accessed"<<endl; 
      for(int n=7; n>=0; n--){
	//print the byte values
	if((((word>>n)<<7)&0x80) != 0){
	  printf("1");
	}
	else{
	  printf("0");
	}
      }
      printf("  ");
    }
    printf("\n");
  }
} // end of PrintDebug

void AliHLTTPCHWClusterTransformComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description needed for the particular component
  if (!targetMap) return;
  
  // OCDB entries for component arguments

  targetMap->Add(new TObjString("HLT/ConfigTPC/TPCHWClusterTransform"), new TObjString("component argument for the charge threshold"));
  
  // OCDB entries to be fetched by the TAXI (access via the AliTPCcalibDB class)
  targetMap->Add(new TObjString("TPC/Calib/Parameters"),    new TObjString("unknown content"));
  targetMap->Add(new TObjString("TPC/Calib/TimeDrift"),     new TObjString("drift velocity calibration"));
  targetMap->Add(new TObjString("TPC/Calib/Temperature"),   new TObjString("temperature map"));
  targetMap->Add(new TObjString("TPC/Calib/PadGainFactor"), new TObjString("gain factor pad by pad"));
  targetMap->Add(new TObjString("TPC/Calib/ClusterParam"),  new TObjString("cluster parameters"));
  
  // OCDB entries needed to be fetched by the Pendolino
  targetMap->Add(new TObjString("TPC/Calib/AltroConfig"), new TObjString("contains the altro config, e.g. info about the L0 trigger timing"));
  targetMap->Add(new TObjString("GRP/CTP/CTPtiming"),     new TObjString("content used in the cluster coordinate transformation in relation to the L0 trigger timing"));

  // OCDB entries necessary for replaying data on the HLT cluster
  targetMap->Add(new TObjString("GRP/GRP/Data"), new TObjString("contains magnetic field info"));  
 
  // OCDB entries needed to suppress fatals/errors/warnings during reconstruction
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
