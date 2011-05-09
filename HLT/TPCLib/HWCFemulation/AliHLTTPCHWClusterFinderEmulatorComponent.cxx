// $Id: AliHLTTPCHWClusterFinderEmulatorComponent.cxx 48710 2011-03-24 12:14:53Z richterm $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timm Steinbeck, Matthias Richter                      *
//* Developers:      Kenneth Aamodt <kenneth.aamodt@student.uib.no>        *
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

/** @file   AliHLTTPCHWClusterFinderEmulatorComponent.cxx
    @author Kenneth Aamodt <kenneth.aamodt@student.uib.no>
    @date   
    @brief  The TPC cluster finder processing component
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCHWClusterFinderEmulatorComponent.h"
#include "AliHLTTPCHWClusterFinderSupport.h"
#include "AliHLTTPCHWClusterFinderEmulator.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusters.h"
#include "AliHLTTPCDefinitions.h"
#include "AliGRPObject.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCParam.h"
#include "AliTPCTransform.h"
#include "AliRawDataHeader.h"
#include "AliHLTTPCClusterMCLabel.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "TGeoGlobalMagField.h"
#include "AliGeomManager.h"

#include <sys/time.h>


AliHLTTPCHWClusterFinderEmulatorComponent::AliHLTTPCHWClusterFinderEmulatorComponent()
  :
  fClusterFinder(NULL),  
  fDeconvTime(kFALSE),  
  fDeconvPad(kFALSE),
  fClusterDeconv(false),  
  fPatch(0),
  fGetActivePads(0),
  fFirstTimeBin(-1),
  fLastTimeBin(-1),
  fDoMC(kFALSE),
  fReleaseMemory( kFALSE ),
  fBenchmark("TPCHWClusterFinderEmulator")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCHWClusterFinderEmulatorComponent::~AliHLTTPCHWClusterFinderEmulatorComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCHWClusterFinderEmulatorComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCHWClusterFinderEmulator";
}

void AliHLTTPCHWClusterFinderEmulatorComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTPCDefinitions::fgkUnpackedRawDataType ); 	 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCHWClusterFinderEmulatorComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTTPCHWClusterFinderEmulatorComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
  tgtList.push_back(kAliHLTDataTypeHwAddr16);
  //tgtList.push_back(AliHLTTPCCADefinitions::fgkCompressedInputDataType);
  return tgtList.size();
}

void AliHLTTPCHWClusterFinderEmulatorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.  
  constBase = 0;
  inputMultiplier = (6 * 0.4);  
}

AliHLTComponent* AliHLTTPCHWClusterFinderEmulatorComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCHWClusterFinderEmulatorComponent();
}
	
int AliHLTTPCHWClusterFinderEmulatorComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  if ( fClusterFinder )
    return -EINPROGRESS;

  fClusterFinder = new AliHLTTPCHWClusterFinderEmulator();

  TObject* pOCDBEntry=LoadAndExtractOCDBObject("GRP/GRP/Data");
  AliGRPObject* pGRP=pOCDBEntry?dynamic_cast<AliGRPObject*>(pOCDBEntry):NULL;
  TString beamType;
  if (pGRP) {
    beamType=pGRP->GetBeamType();
  }

  // first configure the default
  int iResult = 0;
  /*
  TString cdbPath="HLT/ConfigTPC/";
  cdbPath+=GetComponentID();
  iResult=ConfigureFromCDBTObjString(cdbPath, beamType.Data());
  */

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  // return iResult;

  /*
  Int_t iResult=0;
  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }
  
  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }
  */


  fClusterFinder->SetDeconvPad(fDeconvPad);
  fClusterFinder->SetDeconvTime(fDeconvPad);

  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"reco");

  return iResult;
}

int AliHLTTPCHWClusterFinderEmulatorComponent::DoDeinit()
{
  // see header file for class documentation

  if ( fClusterFinder )
    delete fClusterFinder;
  fClusterFinder = NULL;
 
  return 0;
}

struct AliHLTTPCHCWEmulationDigitData{
  UInt_t fSignal;//signal
  UInt_t fTime; // time
  Int_t fLabel[3];//MC label
};

bool CompareDigits( const AliHLTTPCHCWEmulationDigitData &a, const AliHLTTPCHCWEmulationDigitData &b){
  return a.fTime < b.fTime;
}

int AliHLTTPCHWClusterFinderEmulatorComponent::DoEvent( const AliHLTComponentEventData& evtData, 
							const AliHLTComponentBlockData* blocks, 
							AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
							AliHLTUInt32_t& size, 
							vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation

  int iResult=0;
  AliHLTUInt32_t maxSize = size;
  size = 0;
  
  //fTS->SetCurrentTimeStamp(GetTimeStamp());
  //fTS->SetCurrentTimeStamp(0);

  if(!IsDataEvent()){
    return 0;
  }

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  for ( unsigned long ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      const AliHLTComponentBlockData* iter = blocks+ndx;

      if (  iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC) 
	    &&  iter->fDataType != AliHLTTPCDefinitions::fgkUnpackedRawDataType ) continue;

      int slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      int patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
 
      const char *str=Form("slice %d patch %d:", slice, patch);

      fBenchmark.AddInput(iter->fSize);
 
      if (!iter->fPtr) continue;
 
      // create input block for the HW cluster finder

      const AliHLTUInt32_t *rawEvent=0;
      const AliHLTInt32_t *mcLabels = 0;
      AliHLTUInt64_t rawEventSize32 = 0;
      AliHLTUInt64_t mcLabelsSize32 = 0;

      AliHLTTPCHWClusterFinderSupport supp;  

      if( supp.CreateRawEvent( iter, rawEvent, rawEventSize32, mcLabels, mcLabelsSize32 )<0 ) continue; 

      // book memory for the output
      
      AliHLTUInt64_t maxNClusters = rawEventSize32 + 1; // N 32-bit words in input
      AliHLTUInt64_t clustersSize32 = maxNClusters*5;
      AliHLTUInt64_t nOutputMC = maxNClusters;

      AliHLTUInt64_t headerSize = sizeof(AliRawDataHeader);                   
      AliHLTUInt8_t *outBlock = new AliHLTUInt8_t[ headerSize+clustersSize32*sizeof(AliHLTUInt32_t) ];
      
      AliHLTTPCClusterMCData *outMC = reinterpret_cast<AliHLTTPCClusterMCData *>(new AliHLTTPCClusterMCLabel[nOutputMC+1]);
      
      if( !outBlock || !outMC ){
	HLTWarning("%s Not enouth memory!!!", str);
	delete[] outBlock;
	delete[] outMC;
	continue;	
      }
      
      // fill CDH header here, since the HW clusterfinder does not receive it
      
      AliRawDataHeader *cdhHeader = reinterpret_cast<AliRawDataHeader*>(iter->fPtr);
      AliRawDataHeader *outCDHHeader = reinterpret_cast<AliRawDataHeader*>(outBlock);      
      *outCDHHeader = *cdhHeader;
      outCDHHeader->fSize = 0xFFFFFFFF;

      AliHLTUInt32_t *outClusters = reinterpret_cast<AliHLTUInt32_t*> (outBlock + headerSize);
     
      fBenchmark.Start(1);
      
      fClusterFinder->Init(slice, patch );
      
      int err = fClusterFinder->FindClusters( rawEvent, rawEventSize32, 
					      outClusters, clustersSize32, 
					      mcLabels, 
					      outMC );
      //cout<<"HCF finished with error code "<<err<<endl;
      fBenchmark.Stop(1);
      /*
      for( AliHLTUInt64_t i=0; i<clustersSize32; i+=5 ){
	AliHLTUInt32_t *c = outClusters+i;
	AliHLTUInt32_t flag = (c[0]>>30);  
	cout<<flag<<endl;
	if( flag == 0x3){ //beginning of a cluster	   	  
	  int padRow  = (c[0]>>24)&0x3f;
	  int q  = (c[0]&0xFFFFFF)>>6; 
	  Float_t tmpPad   = *((Float_t*)&c[1]);
	  Float_t tmpTime  = *((Float_t*)&c[2]);
	  //Float_t sigmaY2 = *((Float_t*)&c[3]);
	  //Float_t sigmaZ2 = *((Float_t*)&c[4]);
	  cout<<i<<" r "<<padRow<<" p,t ["<<tmpPad<<","<<tmpTime<<"] q "<<q<<endl;
	  if( outMC && outMC->fCount>0 ){
	    cout<<"mc"<<
	      " ("<<outMC->fLabels[i/5].fClusterID[0].fMCID<<" "<<outMC->fLabels[i/5].fClusterID[0].fWeight<<")"
	      " ("<<outMC->fLabels[i/5].fClusterID[1].fMCID<<" "<<outMC->fLabels[i/5].fClusterID[1].fWeight<<")"
	      " ("<<outMC->fLabels[i/5].fClusterID[2].fMCID<<" "<<outMC->fLabels[i/5].fClusterID[2].fWeight<<")"
	      
		<<endl;
	  }
	}
      }
      */      
      
      AliHLTUInt64_t outSize = headerSize + clustersSize32*sizeof(AliHLTUInt32_t);
      
      if( size + outSize <= maxSize ){
	
	memcpy( outputPtr, outBlock, outSize );
	
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = size;
	bd.fSize = outSize;
	bd.fSpecification = iter->fSpecification;
	bd.fDataType = AliHLTTPCDefinitions::fgkHWClustersDataType | kAliHLTDataOriginTPC;
	outputBlocks.push_back( bd );
	fBenchmark.AddOutput(bd.fSize);
	size+= outSize;
	outputPtr+=outSize;     
      } else {
	HLTWarning( "Output buffer (%db) is too small, required %db", maxSize, size+outSize);
	iResult=-ENOSPC;
      }

      
      delete[] outBlock;
      delete[] outMC;
      //delete[] rawEvent;
      //delete[] mcLabels;      
    }
  
  fBenchmark.Stop(0);  
  HLTInfo(fBenchmark.GetStatistics());
  return iResult;
}


int AliHLTTPCHWClusterFinderEmulatorComponent::ScanConfigurationArgument(int argc, const char** argv){

  // see header file for class documentation

  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  if (argument.CompareTo("-solenoidBz")==0){
    if (++i>=argc) return -EPROTO;
    HLTWarning("argument -solenoidBz is deprecated, magnetic field set up globally (%f)", GetBz());
    return 2;
  }

  if (argument.CompareTo("-deconvolute-time")==0){
    HLTDebug("Switching on deconvolution in time direction.");
    fDeconvTime = kTRUE;
    fClusterFinder->SetDeconvTime(fDeconvTime);
    return 1;
  }

  if (argument.CompareTo("-deconvolute-pad")==0){
    HLTDebug("Switching on deconvolution in pad direction.");
    fDeconvPad = kTRUE;
    fClusterFinder->SetDeconvPad(fDeconvPad);
    return 1;
  }

  if (argument.CompareTo("-first-timebin")==0){
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fFirstTimeBin = argument.Atoi();
    if(fFirstTimeBin>=0){
      HLTDebug("fFirstTimeBin set to %d",fFirstTimeBin);
      //fClusterFinder->SetFirstTimeBin(fFirstTimeBin);
    }
    else{
      HLTError("-first-timebin specifier is negative: %d",fFirstTimeBin);
    }
    return 2;
  }

  if (argument.CompareTo("-last-timebin")==0){
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fLastTimeBin = argument.Atoi();
    if(fLastTimeBin<AliHLTTPCTransform::GetNTimeBins()){
      HLTDebug("fLastTimeBin set to %d",fLastTimeBin);
    }
    else{
      HLTError("fLastTimeBins is too big: %d. Maximum: %d",fLastTimeBin,AliHLTTPCTransform::GetNTimeBins());
    }
    return 2;
  }
 
  if (argument.CompareTo("-do-mc")==0) {
    fDoMC=kTRUE;
    //fClusterFinder->SetDoMC(fDoMC);
    HLTDebug("Setting fDoMC to true.");
    return 1;
  }

  if (argument.CompareTo("-active-pads")==0 || argument.CompareTo("activepads")==0){
    if(argument.CompareTo("activepads" )==0){
      HLTWarning("Please change to new component argument naming scheme and use '-active-pads' instead of 'activepads'");
    }
    HLTDebug("Switching on ActivePads");
    fGetActivePads = 1;
    //fClusterFinder->SetDoPadSelection(kTRUE);
    return 1;
  }

  // unknown argument
  return -EINVAL;
}

int AliHLTTPCHWClusterFinderEmulatorComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{  
  // see header file for class documentation

  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigTPC/";
    cdbPath+=GetComponentID();
  }

  return ConfigureFromCDBTObjString(cdbPath.Data());

  /*
  int iResult=0;
  
  const char* path="HLT/ConfigTPC/ClusterFinderComponent";
  if (cdbEntry) path=cdbEntry;
  if (path) {
    HLTInfo("reconfigure from entry %s, chain id %s", path, (chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path);//,GetRunNo());
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object: %s", pString->GetString().Data());
	iResult = Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  return iResult;
  */
}

int AliHLTTPCHWClusterFinderEmulatorComponent::Configure(const char* arguments){
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();

      if (argument.IsNull()) continue;
      
      // -- deconvolute-time option
      if (argument.CompareTo("-deconvolute-time")==0){
	HLTDebug("Switching on deconvolution in time direction.");
	fDeconvTime = kTRUE;
	fClusterFinder->SetDeconvTime(fDeconvTime);
      }
      else if (argument.CompareTo("-deconvolute-pad")==0){
	HLTDebug("Switching on deconvolution in pad direction.");
	fDeconvPad = kTRUE;
	fClusterFinder->SetDeconvPad(fDeconvPad);
      }
      else if (argument.CompareTo("-timebins")==0 || argument.CompareTo("timebins" )==0){
	HLTWarning("Argument %s is depreciated after moving to the offline AliTPCTransform class for xyz calculations.",argument.Data());
	/*
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	AliHLTTPCTransform::SetNTimeBins(((TObjString*)pTokens->At(i))->GetString().Atoi());
	fClusterFinder->UpdateLastTimeBin();
	HLTInfo("number of timebins set to %d, zbin=%f", AliHLTTPCTransform::GetNTimeBins(), AliHLTTPCTransform::GetZWidth());
	*/
	if(argument.CompareTo("timebins")==0){
	  HLTWarning("Argument 'timebins' is old, please switch to new argument naming convention (-timebins). The timebins argument will still work, but please change anyway.");
	}
	
      }
      else if (argument.CompareTo("-first-timebin")==0){
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	fFirstTimeBin = ((TObjString*)pTokens->At(i))->GetString().Atoi();
	if(fFirstTimeBin>=0){
	  HLTDebug("fFirstTimeBin set to %d",fFirstTimeBin);
	  //fClusterFinder->SetFirstTimeBin(fFirstTimeBin);
	}
	else{
	  HLTError("-first-timebin specifier is negative: %d",fFirstTimeBin);
	}
      }
      else if (argument.CompareTo("-last-timebin")==0){
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	fLastTimeBin = ((TObjString*)pTokens->At(i))->GetString().Atoi();
	if(fLastTimeBin<AliHLTTPCTransform::GetNTimeBins()){
	  HLTDebug("fLastTimeBin set to %d",fLastTimeBin);
	}
	else{
	  HLTError("fLastTimeBins is too big: %d. Maximum: %d",fLastTimeBin,AliHLTTPCTransform::GetNTimeBins());
	}
      }
      else if (argument.CompareTo("-do-mc")==0) {
	fDoMC=kTRUE;
	//fClusterFinder->SetDoMC(fDoMC);
	HLTInfo("Setting fDoMC to true.");
      }
      else if (argument.CompareTo("-release-memory")==0) {
	fReleaseMemory = kTRUE;
	//fClusterFinder->SetReleaseMemory( kTRUE );
	HLTInfo("Setting fReleaseMemory to true.");
      }
      else if (argument.CompareTo("-active-pads")==0 || argument.CompareTo("activepads")==0){
	if(argument.CompareTo("activepads" )==0){
	  HLTWarning("Please change to new component argument naming scheme and use '-active-pads' instead of 'activepads'");
	}
	HLTDebug("Switching on ActivePads");
	fGetActivePads = 1;
	//fClusterFinder->SetDoPadSelection(kTRUE);
      }
      else if (argument.CompareTo("-update-calibdb")==0){
	//fClusterFinder->UpdateCalibDB();
      }
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTTPCHWClusterFinderEmulatorComponent::GetOCDBObjectDescription( TMap* const targetMap){
// Get a list of OCDB object description needed for the particular component
  
  if (!targetMap) return;
  
  // OCDB entries for component arguments
  targetMap->Add(new TObjString("HLT/ConfigTPC/TPCClusterFinderUnpacked"), new TObjString("component arguments, empty at the moment"));
    
  
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
