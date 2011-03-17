// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk <st05886@alf.uib.no>                   *
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

/// @file   AliHLTITSClusterFinderComponent.cxx
/// @author Gaute Ovrebekk <st05886@alf.uib.no>
/// @date   
/// @brief  Component to run offline clusterfinders
///

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTITSClusterFinderComponent.h" 

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliHLTITSSpacePointData.h"
#include "AliHLTITSClusterDataFormat.h"
#include <AliHLTDAQ.h>
#include "AliGeomManager.h"
#include "AliITSRecoParam.h"
#include "AliITSReconstructor.h"
#include "AliHLTITSClusterFinderSPD.h"
#include "AliHLTITSClusterFinderSSD.h"
#include "TMap.h"
#include "AliITSRecPointContainer.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include <cstdlib>
#include <cerrno>
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include <sys/time.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSClusterFinderComponent);

AliHLTITSClusterFinderComponent::AliHLTITSClusterFinderComponent(int mode)
  :
  fModeSwitch(mode),
  fInputDataType(kAliHLTVoidDataType),
  fOutputDataType(kAliHLTVoidDataType),
  fUseOfflineFinder(0),
  fNModules(0),
  fId(0),
  fNddl(0),
  fRawReader(NULL),
  fDettype(NULL),
  fgeom(NULL),
  fgeomInit(NULL),
  fSPD(NULL),
  fSSD(NULL),
  tD(NULL),
  tR(NULL),
  fSPDNModules(0),
  fSDDNModules(0),
  fSSDNModules(0),
  fFirstModule(0),
  fLastModule(0),
  fnClusters(0),
  fclusters(),
  fBenchmark(GetComponentID()),
  fOutputSizeOffset(0),
  fInputMultiplierDigits(20),
  fpLoader(NULL)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  switch(fModeSwitch){
  case kClusterFinderSPD:
    fInputDataType  = kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD;
    fOutputDataType = kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD;
    break;
  case kClusterFinderSDD: 	 
    fInputDataType  = kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSDD;
    fOutputDataType = kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD;
    break;
  case kClusterFinderSSD:
    fInputDataType  = kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD;
    fOutputDataType = kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD;
    break;
  case kClusterFinderDigits:
    fInputDataType  = kAliHLTDataTypeAliTreeD|kAliHLTDataOriginITS;
    fOutputDataType = kAliHLTDataTypeClusters|kAliHLTDataOriginITS;
    break;
  default:
    HLTFatal("unknown cluster finder");
  }
}

AliHLTITSClusterFinderComponent::~AliHLTITSClusterFinderComponent() 
{
  // see header file for class documentation
  delete fRawReader;
  delete fDettype;
  delete fgeomInit;  
  delete fSPD;
  delete fSSD;
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTITSClusterFinderComponent::GetComponentID()
{
  // see header file for class documentation
  switch(fModeSwitch){
  case kClusterFinderSPD:
    return "ITSClusterFinderSPD";
    break;
  case kClusterFinderSDD: 	 
    return "ITSClusterFinderSDD";
    break;
  case kClusterFinderSSD:
    return "ITSClusterFinderSSD";
    break;
  case kClusterFinderDigits:
    return "ITSClusterFinderDigits";
    break;
  }
  return "";
}

void AliHLTITSClusterFinderComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) 
{
  // see header file for class documentation
  list.clear(); 
  list.push_back( fInputDataType );
}

AliHLTComponentDataType AliHLTITSClusterFinderComponent::GetOutputDataType() 
{
  // see header file for class documentation
  return fOutputDataType;
}

void AliHLTITSClusterFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation
  constBase = fOutputSizeOffset;
  switch(fModeSwitch){
  case kClusterFinderDigits:
    inputMultiplier = fInputMultiplierDigits;
    break;
  case kClusterFinderSPD:
  case kClusterFinderSDD: 	 
  case kClusterFinderSSD:
  default:
    inputMultiplier = 20;
  }
}

AliHLTComponent* AliHLTITSClusterFinderComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTITSClusterFinderComponent(fModeSwitch);
}
	
Int_t AliHLTITSClusterFinderComponent::DoInit( int argc, const char** argv ) {
  // see header file for class documentation
  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"reco");

  /*
  fStatTime = 0;
  fStatTimeAll = 0;
  fStatTimeC = 0;
  fStatTimeAllC = 0;
  fStatNEv = 0;
  */
  
  Int_t runNo = GetRunNo();
  AliCDBStorage* store = AliCDBManager::Instance()->GetDefaultStorage();
  if (!store) {
    return -ENOENT;
  }

  bool cdbOK = true;
  //OCDB for SPD
  if(store->GetLatestVersion("ITS/Calib/SPDNoisy", runNo)<0){
    HLTError("SPDNoisy is not found in SPD/Calib");
    cdbOK = false;
  }
  if(store->GetLatestVersion("ITS/Calib/SPDDead", runNo)<0){
    HLTError("SPDDead is not found in SPD/Calib");
    cdbOK = false;
  }
  if(store->GetLatestVersion("TRIGGER/SPD/PITConditions", runNo)<0){
    HLTError("PITConditions is not found in TRIGGER/SPD");
    cdbOK = false;
  }
  
  //OCDB for SDD
  if(store->GetLatestVersion("ITS/Calib/CalibSDD", runNo)<0){
    HLTError("CalibSDD is not found in ITS/Calib");
    cdbOK = false;
  }
  if(store->GetLatestVersion("ITS/Calib/RespSDD", runNo)<0){
    HLTError("RespSDD is not found in ITS/Calib");
    cdbOK = false;
  }
  if(store->GetLatestVersion("ITS/Calib/DriftSpeedSDD", runNo)<0){
    HLTError("DriftSpeedSDD is not found in ITS/Calib");
    cdbOK = false;
  }
  if(store->GetLatestVersion("ITS/Calib/DDLMapSDD", runNo)<0){
    HLTError("DDLMapSDD is not found in ITS/Calib");
    cdbOK = false;
  }
  if(store->GetLatestVersion("ITS/Calib/MapsTimeSDD", runNo)<0){
    HLTError("MapsTimeSDD is not found in ITS/Calib");
    cdbOK = false;
  }

  //OCDB for SSD
  if(store->GetLatestVersion("ITS/Calib/NoiseSSD", runNo)<0){
    HLTError("NoiseSSD is not found in ITS/Calib");
    cdbOK = false;
  }
  if(store->GetLatestVersion("ITS/Calib/GainSSD", runNo)<0){
    HLTError("GainSSD is not found in ITS/Calib");
    cdbOK = false;
  }
  if(store->GetLatestVersion("ITS/Calib/BadChannelsSSD", runNo)<0){
    HLTError("BadChannelsSSD is not found in ITS/Calib");
    cdbOK = false;
  }
  
  //General reconstruction
  if(store->GetLatestVersion("GRP/CTP/Scalers", runNo)<0){
    HLTError("Scalers is not found in GRP/CTP/");
    cdbOK = false;
  }
  if(!cdbOK){return -EIO;}

  if(fModeSwitch==kClusterFinderSPD) {
    HLTDebug("using ClusterFinder for SPD");
    //fNModules=AliITSgeomTGeo::GetNDetectors(1)*AliITSgeomTGeo::GetNLadders(1) + AliITSgeomTGeo::GetNDetectors(2)*AliITSgeomTGeo::GetNLadders(2);
    fId=AliHLTDAQ::DdlIDOffset("ITSSPD");
    fNddl=AliHLTDAQ::NumberOfDdls("ITSSPD");
  }
  else if(fModeSwitch==kClusterFinderSDD) {
    HLTDebug("using ClusterFinder for SDD");
    //fNModules=AliITSgeomTGeo::GetNDetectors(3)*AliITSgeomTGeo::GetNLadders(3) + AliITSgeomTGeo::GetNDetectors(4)*AliITSgeomTGeo::GetNLadders(4);
    fId=AliHLTDAQ::DdlIDOffset("ITSSDD");
    fNddl=AliHLTDAQ::NumberOfDdls("ITSSDD");
  }
  else if(fModeSwitch==kClusterFinderSSD) {
    HLTDebug("using ClusterFinder for SSD");
    //fNModules=AliITSgeomTGeo::GetNDetectors(5)*AliITSgeomTGeo::GetNLadders(5) + AliITSgeomTGeo::GetNDetectors(6)*AliITSgeomTGeo::GetNLadders(6);
    fId=AliHLTDAQ::DdlIDOffset("ITSSSD");
    fNddl=AliHLTDAQ::NumberOfDdls("ITSSSD");
  }
  else if(fModeSwitch==kClusterFinderDigits) {
    tR = new TTree();
  }
  else{
     HLTFatal("No mode set for clusterfindercomponent");
  }

  //Removed the warning for loading default RecoParam in HLT
  AliITSRecoParam *par = AliITSRecoParam::GetLowFluxParam();
  AliITSReconstructor *rec = new AliITSReconstructor();
  rec->SetRecoParam(par);
  
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()){
    HLTError("Default CDB storage has not been set !");
    return -ENOENT;
  }

  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }

  fgeomInit = new AliITSInitGeometry();
  //fgeomInit = new AliITSInitGeometry(kvPPRasymmFMD,2);
  //fgeomInit->InitAliITSgeom(fgeom);
  fgeom = fgeomInit->CreateAliITSgeom();
 
  fNModules = fgeom->GetIndexMax();
  Int_t modperlay[6];
  for(Int_t i=0;i<6;i++)modperlay[i]=AliITSgeomTGeo::GetNDetectors(1+i)*AliITSgeomTGeo::GetNLadders(1+i);
  fSPDNModules=modperlay[0]+modperlay[1];
  fSDDNModules=modperlay[2]+modperlay[3];
  fSSDNModules=modperlay[4]+modperlay[5];
  
  if(fModeSwitch==kClusterFinderSPD) {
    fFirstModule=0;
    fLastModule=fSPDNModules;
  }
  else if(fModeSwitch==kClusterFinderSDD) {
     fFirstModule=fSPDNModules;
     fLastModule=fFirstModule + fSDDNModules;
  }
  else if(fModeSwitch==kClusterFinderSSD) {
    fFirstModule=fSPDNModules + fSDDNModules;
    fLastModule=fFirstModule + fSSDNModules;
  }
 
  //set dettype
  fDettype = new AliITSDetTypeRec();
  fDettype->SetITSgeom(fgeom); 
  fDettype->SetDefaults();
  fDettype->SetDefaultClusterFindersV2(kTRUE); 

  if ( fRawReader )
    return -EINPROGRESS;

  fRawReader = new AliRawReaderMemory();
  fSPD = new AliHLTITSClusterFinderSPD( fDettype );
  fSSD = new AliHLTITSClusterFinderSSD( fDettype, fRawReader );

  TString arguments = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !arguments.IsNull() ) arguments += " ";
    arguments += argv[i];
  }

  tD = NULL;
  //tR = NULL;

  return Configure( arguments.Data() );
}

Int_t AliHLTITSClusterFinderComponent::DoDeinit() {
  // see header file for class documentation

  if ( fRawReader )
    delete fRawReader;
  fRawReader = NULL;

  if ( fDettype )
    delete fDettype;
  fDettype = NULL;

  if ( fgeomInit )
    delete fgeomInit;
  fgeomInit = NULL;
  
  delete fSPD;
  fSPD = 0;

  delete fSSD;
  fSSD = 0;

  fUseOfflineFinder = 0;

  if (fpLoader) {
    fpLoader->UnloadDigits();
  }
  fpLoader=NULL;

  return 0;
}

int AliHLTITSClusterFinderComponent::DoEvent
(
 const AliHLTComponentEventData& evtData,
 const AliHLTComponentBlockData* /*blocks*/,
 AliHLTComponentTriggerData& /*trigData*/,
 AliHLTUInt8_t* outputPtr,
 AliHLTUInt32_t& size,
  vector<AliHLTComponentBlockData>& outputBlocks )
{
 // see header file for class documentation
  
  AliHLTUInt32_t maxBufferSize = size;
  size = 0; // output size

  if (!IsDataEvent()) return 0;
  
  if ( evtData.fBlockCnt<=0 )
    {
      HLTDebug("no blocks in event" );
      return 0;
    }

  AliHLTUInt32_t totalInputSize=0;
  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);
  for( const AliHLTComponentBlockData *i= GetFirstInputBlock(fInputDataType); i!=NULL; i=GetNextInputBlock() ){
    totalInputSize+=i->fSize;
    fBenchmark.AddInput(i->fSize);
  }

  Int_t ret = 0;

  if(fModeSwitch==kClusterFinderDigits) {

    for ( const TObject *iter = GetFirstInputObject(fInputDataType); iter != NULL; iter = GetNextInputObject() ) {  
      tD = dynamic_cast<TTree*>(const_cast<TObject*>( iter ) );
      if(!tD){
	HLTFatal("No Digit Tree found");
	return -1;
      }
      // 2010-04-17 very crude workaround: TTree objects are difficult to send
      // The actual case: Running ITS and TPC reconstruction fails at the second event
      // to read the ITS digits from the TreeD
      //
      // Reason: reading fails in TBranch::GetBasket, there a new basket is opened from
      // a TFile object. The function TBranch::GetFile returns the file object from
      // an internal fDirectory (TDirectory) object. This file is at the second event
      // set to the TPC.Digits.root. The internal mismatch creates a seg fault
      //
      // Investigation: TBranch::Streamer uses a crude assignment after creating the
      // TBranch object
      //    fDirectory = gDirectory;
      // gDirectory is obviously not set correctly. Setting the directory to a TFile
      // object for the ITS digits helps to fix the internal mess. Tried also to set
      // the Directory for the TreeD to NULL (This has only effect if ones sets it 
      // to something not NULL first, and then to NULL). But then no content, i.e.
      // ITS clusters could be retrieved.
      //
      // Conclusion: TTree objects are hardly to be sent via TMessage, there are direct
      // links to the file required anyhow.
      //
      // 2011-01-28 hotfix reloaded: accessing the files like that fails if there are
      // multiple digit files because of a large number of events. New ugly fix is to
      // use the global runloader instance to get hold on the digits tree.
      fnClusters = 0;
      AliRunLoader* pRunLoader=AliRunLoader::Instance();
      if (!pRunLoader) {
	HLTError("failed to get global runloader instance");
	return -ENOSYS;
      }
      // get the specific loader for the module
      if (!fpLoader) {
	const char* loaderType="ITSLoader";
	fpLoader=pRunLoader->GetLoader(loaderType);
	if (!fpLoader) {
	  HLTError("can not get loader \"%s\" from runloader", loaderType);
	  return -ENOSYS;
	}
	// prepare the loader
	fpLoader->LoadDigits("read");
      }
      pRunLoader->GetEvent(GetEventCount());

      tD=fpLoader->TreeD();
      tR->Reset();
      tR->SetDirectory(0);
      fDettype->SetTreeAddressD(tD);
      fDettype->MakeBranch(tR,"R");
      fDettype->SetTreeAddressR(tR);
      Option_t *opt="All";
      fBenchmark.Start(1);
      fDettype->DigitsToRecPoints(tD,tR,0,opt,0);
      fBenchmark.Stop(1);
      TClonesArray * fRecPoints = NULL;
      tR->SetBranchAddress("ITSRecPoints",&fRecPoints);
      for(Int_t treeEntry=0;treeEntry<tR->GetEntries();treeEntry++){
	tR->GetEntry(treeEntry);
	fnClusters += fRecPoints->GetEntries();
      }

      UInt_t bufferSize = fnClusters * sizeof(AliHLTITSSpacePointData) + sizeof(AliHLTITSClusterData);
      if( size + bufferSize > maxBufferSize ){
	//HLTWarning( "Output buffer size exceed (buffer size %d, required size %d)", maxBufferSize, size+bufferSize);
	if (totalInputSize>0) {
	  fInputMultiplierDigits=(float)(size+bufferSize)/totalInputSize;
	  fInputMultiplierDigits+=1.;
	} else {
	  fOutputSizeOffset=totalInputSize;
	  fInputMultiplierDigits=1.;
	}
	ret = -ENOSPC;      
	break;		
      }
      if( fnClusters>0 ){
	fBenchmark.Start(1);
	RecPointToSpacePoint(outputPtr,size);
	fBenchmark.Stop(1);
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = size;
	bd.fSize = bufferSize;
	bd.fSpecification = 0x00000000;
	bd.fDataType = GetOutputDataType();
	outputBlocks.push_back( bd );
	size += bufferSize;
	fBenchmark.AddOutput(bd.fSize);
      }
    }
  }
  else{

    AliITSRecPointContainer* rpc = AliITSRecPointContainer::Instance();

    // -- Loop over blocks
    for( const AliHLTComponentBlockData* iter = GetFirstInputBlock(fInputDataType); iter != NULL; iter = GetNextInputBlock() ) {

      if(fUseOfflineFinder){
	if(fModeSwitch==kClusterFinderSPD){rpc->ResetSPD();}
	if(fModeSwitch==kClusterFinderSSD){rpc->ResetSSD();}
      }
      if(fModeSwitch==kClusterFinderSDD){rpc->ResetSDD();}
      
      // -- Debug output of datatype --
      HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	       evtData.fEventID, evtData.fEventID, 
	       DataType2Text(iter->fDataType).c_str(), 
	       DataType2Text(fInputDataType).c_str());
      
      // -- Check for the correct data type
      if ( iter->fDataType != (fInputDataType) )  
	continue;
      
      // -- Get equipment ID out of specification
      AliHLTUInt32_t spec = iter->fSpecification;
      
      Int_t id = fId;
      for ( Int_t ii = 0; ii < fNddl ; ii++ ) {   //number of ddl's
	if ( spec & 0x00000001 ) {
	  id += ii;
	  break;
	}
	spec = spec >> 1 ;
      }
      
      // -- Set equipment ID to the raw reader
      
      if(!fRawReader){
	HLTWarning("The fRawReader pointer is NULL");
	continue;
      }

      if(!fRawReader->AddBuffer((UChar_t*) iter->fPtr, iter->fSize, id)){
	HLTWarning("Could not add buffer");
      }

      fBenchmark.Start(1);

      fnClusters = 0;

      if(fModeSwitch==kClusterFinderSPD && !fUseOfflineFinder){ fSPD->RawdataToClusters( fRawReader, fclusters ); }
      else if(fModeSwitch==kClusterFinderSSD && !fUseOfflineFinder){ fSSD->RawdataToClusters( fclusters ); }
      else{
	if(fModeSwitch==kClusterFinderSPD && fUseOfflineFinder) {fDettype->DigitsToRecPoints(fRawReader,"SPD");}
	if(fModeSwitch==kClusterFinderSSD && fUseOfflineFinder) {fDettype->DigitsToRecPoints(fRawReader,"SSD");}
	if(fModeSwitch==kClusterFinderSDD) {fDettype->DigitsToRecPoints(fRawReader,"SDD");}
	TClonesArray* clusters = NULL;
	for(int i=fFirstModule;i<fLastModule;i++){
	  clusters = rpc->UncheckedGetClusters(i);
	  if(clusters != NULL){
	    fnClusters += clusters->GetEntriesFast();
	  }
	}     
      }
      
      /* Do not work
	 if(fModeSwitch==kClusterfinderSPD){
	 fnClusters = rpc->GetNClustersInLayerFast(1) + rpc->GetNClustersInLayerFast(2);
	 }
	 if(fModeSwitch==kClusterfinderSDD){
	 fnClusters = rpc->GetNClustersInLayerFast(3) + rpc->GetNClustersInLayerFast(4);
	 }
	 if(fModeSwitch==kClusterfinderSSD){
	 fnClusters = rpc->GetNClustersInLayerFast(5) + rpc->GetNClustersInLayerFast(6);
	 }
      */

      fBenchmark.Stop(1);
      
      fRawReader->ClearBuffers();    
      
      UInt_t nClusters=fclusters.size();
      if(nClusters>0){fnClusters = nClusters;}

      UInt_t bufferSize = fnClusters * sizeof(AliHLTITSSpacePointData) + sizeof(AliHLTITSClusterData);
      if( size + bufferSize > maxBufferSize ){
	HLTWarning( "Output buffer size exceed (buffer size %d, current size %d)", maxBufferSize, size+bufferSize);
	ret = -ENOSPC;      
	break;		
      }
      if( fnClusters>0 ){

	RecPointToSpacePoint(outputPtr,size);

	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = size;
	bd.fSize = bufferSize;
	bd.fSpecification = iter->fSpecification;
	bd.fDataType = GetOutputDataType();
	outputBlocks.push_back( bd );
	size += bufferSize;
	fBenchmark.AddOutput(bd.fSize);
	if(nClusters>0){fclusters.clear();}	
      }
      
    } // input blocks
  }

  fBenchmark.Stop(0);
  HLTInfo(fBenchmark.GetStatistics());
  return ret;
}

int AliHLTITSClusterFinderComponent::Configure(const char* arguments)
{
  
  int iResult=0;
  
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  
  TObjArray* pTokens=allArgs.Tokenize(" ");
  
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;      
      if (argument.CompareTo("-use-offline-finder")==0) {
	fUseOfflineFinder = 1;
	HLTInfo("Off-line ClusterFinder will be used");
	continue;
      }
      /*
      else if (argument.CompareTo("")==0) {
	HLTInfo("");
	continue;
      }
      */
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  
  return iResult;
}

int AliHLTITSClusterFinderComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  
  const char* path="";

  switch(fModeSwitch){
  case kClusterFinderSPD:
    path = "HLT/ConfigITS/ITSClusterFinderSPD";
    break;
  case kClusterFinderSDD: 	 
    path = "HLT/ConfigITS/ITSClusterFinderSDD";
    break;
  case kClusterFinderSSD:
    path = "HLT/ConfigITS/ITSClusterFinderSSD";
    break;
  case kClusterFinderDigits:
    path = "";
    break;
  default:
    HLTFatal("unknown cluster finder");
  }

  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  
  return iResult;
}
void AliHLTITSClusterFinderComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description.
  if (!targetMap) return;
  //SPD
  targetMap->Add(new TObjString("ITS/Calib/SPDNoisy"),new TObjString("Calibration object for SPD" ));
  targetMap->Add(new TObjString("ITS/Calib/SPDDead"),new TObjString("Calibration object for SPD" ));
  targetMap->Add(new TObjString("TRIGGER/SPD/PITConditions"),new TObjString("Calibration object for SPD" ));
  //SDD
  targetMap->Add(new TObjString("ITS/Calib/CalibSDD"),new TObjString("Calibration object for SDD" ));
  targetMap->Add(new TObjString("ITS/Calib/RespSDD"),new TObjString("Calibration object for SDD" ));
  targetMap->Add(new TObjString("ITS/Calib/DriftSpeedSDD"),new TObjString("Calibration object for SDD" ));
  targetMap->Add(new TObjString("ITS/Calib/DDLMapSDD"),new TObjString("Calibration object for SDD" ));
  targetMap->Add(new TObjString("ITS/Calib/MapsTimeSDD"),new TObjString("Calibration object for SDD" ));
  //SSD
  targetMap->Add(new TObjString("ITS/Calib/NoiseSSD"),new TObjString("Calibration object for SSD" ));
  targetMap->Add(new TObjString("ITS/Calib/GainSSD"),new TObjString("Calibration object for SSD" ));
  targetMap->Add(new TObjString("ITS/Calib/BadChannelsSSD"),new TObjString("Calibration object for SSD" ));
  //General reconstruction
  targetMap->Add(new TObjString("GRP/CTP/Scalers"),new TObjString("General reconstruction object" ));
}


void AliHLTITSClusterFinderComponent::RecPointToSpacePoint(AliHLTUInt8_t* outputPtr,AliHLTUInt32_t& size){
  AliHLTITSClusterData *outputClusters = reinterpret_cast<AliHLTITSClusterData*>(outputPtr + size);
  outputClusters->fSpacePointCnt=fnClusters;    
  int clustIdx=0;
  if(fModeSwitch==kClusterFinderDigits) {
    TClonesArray * fRecPoints = NULL;
    tR->SetBranchAddress("ITSRecPoints",&fRecPoints);
    for(Int_t treeEntry=0;treeEntry<tR->GetEntries();treeEntry++){
      tR->GetEntry(treeEntry);
      fnClusters += fRecPoints->GetEntries();
      for(Int_t tCloneEntry=0;tCloneEntry<fRecPoints->GetEntries();tCloneEntry++){
	AliITSRecPoint *recpoint=(AliITSRecPoint*)fRecPoints->At(tCloneEntry);
	RecpointToOutput(outputClusters,recpoint,clustIdx);
      }
    } 
  }
  else if(fclusters.size()>0){
    for(UInt_t i=0;i<fclusters.size();i++){
      AliITSRecPoint *recpoint = (AliITSRecPoint*) &(fclusters[i]);
      RecpointToOutput(outputClusters,recpoint,clustIdx);
    }
  }
  else{
    AliITSRecPointContainer* rpc = AliITSRecPointContainer::Instance();
    TClonesArray* clusters = NULL;
    for(Int_t i=fFirstModule;i<fLastModule;i++){
      clusters = rpc->UncheckedGetClusters(i);
      for(Int_t j=0;j<clusters->GetEntriesFast();j++){
	AliITSRecPoint *recpoint = (AliITSRecPoint*) clusters->At(j);
	RecpointToOutput(outputClusters,recpoint,clustIdx);
      }
    } 
  }
}

void AliHLTITSClusterFinderComponent::RecpointToOutput(AliHLTITSClusterData *outputClusters, AliITSRecPoint *recpoint, int &clustIdx){
  outputClusters->fSpacePoints[clustIdx].fY=recpoint->GetY();
  outputClusters->fSpacePoints[clustIdx].fZ=recpoint->GetZ();
  outputClusters->fSpacePoints[clustIdx].fSigmaY2=recpoint->GetSigmaY2();
  outputClusters->fSpacePoints[clustIdx].fSigmaZ2=recpoint->GetSigmaZ2();
  outputClusters->fSpacePoints[clustIdx].fSigmaYZ=recpoint->GetSigmaYZ();
  outputClusters->fSpacePoints[clustIdx].fQ=recpoint->GetQ();
  outputClusters->fSpacePoints[clustIdx].fNy=recpoint->GetNy();
  outputClusters->fSpacePoints[clustIdx].fNz=recpoint->GetNz();
  outputClusters->fSpacePoints[clustIdx].fLayer=recpoint->GetLayer();
  outputClusters->fSpacePoints[clustIdx].fIndex=recpoint->GetDetectorIndex() | recpoint->GetPindex() | recpoint->GetNindex();
  outputClusters->fSpacePoints[clustIdx].fTracks[0]=recpoint->GetLabel(0);
  outputClusters->fSpacePoints[clustIdx].fTracks[1]=recpoint->GetLabel(1);
  outputClusters->fSpacePoints[clustIdx].fTracks[2]=recpoint->GetLabel(2);      
  clustIdx++;
}
