// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Øvrebekk <st05886@alf.uib.no>                   *
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

/** @file   AliHLTITSClusterFinderComponent.cxx
    @author Gaute Øvrebekk <st05886@alf.uib.no>
    @date   
    @brief  Component to run offline clusterfinders
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTITSClusterFinderComponent.h" 

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliHLTDataTypes.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliHLTITSSpacePointData.h"
#include "AliHLTITSClusterDataFormat.h"
#include <AliHLTDAQ.h>
#include "AliGeomManager.h"
#include "AliITSRecoParam.h"
#include "AliITSReconstructor.h"
#include "AliHLTITSClusterFinderSPD.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include <sys/time.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSClusterFinderComponent);

AliHLTITSClusterFinderComponent::AliHLTITSClusterFinderComponent(int mode)
  :
  fModeSwitch(mode),
  fNModules(0),
  fId(0),
  fNddl(0),
  fClusters(NULL),
  fRawReader(NULL),
  fDettype(NULL),
  fgeom(NULL),
  fgeomInit(NULL),
  fSPD(NULL)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (fModeSwitch!=kClusterFinderSPD &&
      fModeSwitch!=kClusterFinderSDD &&
      fModeSwitch!=kClusterFinderSSD) {
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
  }
  return "";
}

void AliHLTITSClusterFinderComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  switch(fModeSwitch){
  case kClusterFinderSPD:
    list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD );
    break;
  case kClusterFinderSDD: 	 
    list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSDD );
    break;
  case kClusterFinderSSD:
    list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD );
    break;
  }
}

AliHLTComponentDataType AliHLTITSClusterFinderComponent::GetOutputDataType() {
  // see header file for class documentation
  switch(fModeSwitch){
  case kClusterFinderSPD:
    return kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD;
    break;
  case kClusterFinderSDD: 	 
    return kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD;
    break;
  case kClusterFinderSSD:
    return kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD;
    break;
  }
  return kAliHLTVoidDataType;
}

void AliHLTITSClusterFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation
  constBase = 0;
  inputMultiplier = 100;
}

AliHLTComponent* AliHLTITSClusterFinderComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTITSClusterFinderComponent(fModeSwitch);
}
	
Int_t AliHLTITSClusterFinderComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation

  fStatTime = 0;
  fStatTimeAll = 0;
  fStatTimeC = 0;
  fStatTimeAllC = 0;
  fStatNEv = 0;

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

 //fgeomInit = new AliITSInitGeometry(kvSPD02,2);
  fgeomInit = new AliITSInitGeometry(kvPPRasymmFMD,2);
  //fgeomInit->InitAliITSgeom(fgeom);
  fgeom = fgeomInit->CreateAliITSgeom();
 
  fNModules = fgeom->GetIndexMax();

  fClusters = new TClonesArray*[fNModules]; 
  for (Int_t iModule = 0; iModule < fNModules; iModule++) {
    fClusters[iModule] = NULL;
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

  return 0;
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

  for (Int_t iModule = 0; iModule < fNModules; iModule++) {
    if(fClusters[iModule] != NULL){
      fClusters[iModule]->Delete();
      delete fClusters[iModule];
    }
    fClusters[iModule] = NULL;
  } 
  
  return 0;
}

#include "TStopwatch.h"

Int_t AliHLTITSClusterFinderComponent::DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/)
{  // see header file for class documentation

  // -- Iterator over Data Blocks --
  const AliHLTComponentBlockData* iter = NULL;

  //if (!IsDataEvent()) return 0;

  if ( evtData.fBlockCnt<=0 )
      {
	HLTDebug("no blocks in event" );
	return 0;
      }


  AliHLTComponentDataType datatype;
  if(fModeSwitch==kClusterFinderSPD) { datatype = kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD;}
  else if(fModeSwitch==kClusterFinderSDD) {datatype = kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSDD;}
  else if(fModeSwitch==kClusterFinderSSD) {datatype = kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD;}
  
  TStopwatch timer;

  // -- Loop over blocks
  for ( iter = GetFirstInputBlock(datatype); iter != NULL; iter = GetNextInputBlock() ) {
  
    // -- Debug output of datatype --
    HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	       evtData.fEventID, evtData.fEventID, 
	       DataType2Text(iter->fDataType).c_str(), 
	       DataType2Text(datatype).c_str());
    
    // -- Check for the correct data type
    if ( iter->fDataType != (datatype) )  
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

    if(!fRawReader->AddBuffer((UChar_t*) iter->fPtr, iter->fSize, id)){
      HLTWarning("Could not add buffer");
    }
    TStopwatch timer1;

    std::vector<AliITSRecPoint> vclusters;
    
    if(fModeSwitch==kClusterFinderSPD) {fSPD->RawdataToClusters( fRawReader, vclusters );}
    else{
      if(fModeSwitch==kClusterFinderSDD) {fDettype->DigitsToRecPoints(fRawReader,fClusters,"SDD");}
      if(fModeSwitch==kClusterFinderSSD) {fDettype->DigitsToRecPoints(fRawReader,fClusters,"SSD");}
      for(int i=0;i<fNModules;i++){
	if(fClusters[i] != NULL){
	  for(int j=0;j<fClusters[i]->GetEntriesFast();j++){
	    AliITSRecPoint *recpoint = (AliITSRecPoint*) (fClusters[i]->At(j));
	    vclusters.push_back(*recpoint);
	  }
	  fClusters[i]->Delete();
	  delete fClusters[i];
	}
	fClusters[i] = NULL;
      }     
    }
    
    timer1.Stop();
    fStatTime+=timer1.RealTime();
    fStatTimeC+=timer1.CpuTime();
    
    UInt_t nClusters=vclusters.size();
    
    UInt_t bufferSize = nClusters * sizeof(AliHLTITSSpacePointData) + sizeof(AliHLTITSClusterData);
    AliHLTUInt8_t *buffer = new AliHLTUInt8_t[bufferSize];
    AliHLTITSClusterData *outputClusters = reinterpret_cast<AliHLTITSClusterData*>(buffer);
    outputClusters->fSpacePointCnt=nClusters;
    
    int clustIdx=0;
    for(int i=0;i<vclusters.size();i++){
      AliITSRecPoint *recpoint = (AliITSRecPoint*) &(vclusters[i]);
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

    if(fModeSwitch==kClusterFinderSPD) {PushBack(buffer,bufferSize,kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD,iter->fSpecification);}
    else if(fModeSwitch==kClusterFinderSDD) {PushBack(buffer,bufferSize,kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD,iter->fSpecification);}
    else if(fModeSwitch==kClusterFinderSSD) {PushBack(buffer,bufferSize,kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD,iter->fSpecification);}
  
    delete[] buffer; 
    fRawReader->ClearBuffers();
    
  } //  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {    

  timer.Stop();

  fStatTimeAll+=timer.RealTime();
  fStatTimeAllC+=timer.CpuTime();
  fStatNEv++;
  //if( fStatNEv%10000==0 && fStatTimeAll>0.0 && fStatTime>0.0 && fStatTimeAllC>0.0 && fStatTimeC>0.0)
  //cout<<fStatTimeAll/fStatNEv*1.e3<<" "<<fStatTime/fStatNEv*1.e3<<" "
  //<<fStatTimeAllC/fStatNEv*1.e3<<" "<<fStatTimeC/fStatNEv*1.e3<<" ms"<<endl;

  return 0;
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
      /*
      if (argument.CompareTo("")==0) {
	HLTInfo("");
	continue;
      }

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
  
  const char* path="HLT/ConfigITS/ClusterFinderComponent";
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

