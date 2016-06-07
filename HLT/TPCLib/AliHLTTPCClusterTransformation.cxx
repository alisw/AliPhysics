// $Id: AliHLTTPCClusterTransformation.cxx 41244 2010-05-14 08:13:35Z kkanaki $

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

/** @file   AliHLTTPCClusterTransformation.cxx
    @author Kalliopi Kanaki, Sergey Gorbubnov
    @date   
    @brief 
*/


#include "AliHLTTPCClusterTransformation.h"
#include "AliHLTTPCFastTransform.h"
#include "AliHLTTPCDefinitions.h"

#include "AliCDBPath.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliTPCcalibDB.h"
#include "AliTPCTransform.h"
#include "AliTPCParam.h"
#include "AliTPCRecoParam.h"
#include "AliGeomManager.h"
#include "AliRunInfo.h"
#include "AliEventInfo.h"
#include "AliRawEventHeaderBase.h"
#include <iostream>
#include <iomanip>
#include <TGeoGlobalMagField.h>

using namespace std;

ClassImp(AliHLTTPCClusterTransformation) //ROOT macro for the implementation of ROOT specific class methods

AliRecoParam AliHLTTPCClusterTransformation::fOfflineRecoParam;

AliHLTTPCClusterTransformation::AliHLTTPCClusterTransformation()
:
  fError(),
  fFastTransform()  
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt  
}

AliHLTTPCClusterTransformation::~AliHLTTPCClusterTransformation() 
{ 
  // see header file for class documentation
}


int  AliHLTTPCClusterTransformation::Init( double FieldBz, Long_t TimeStamp )
{
  // Initialisation
 
  if(!AliGeomManager::GetGeometry()){
    AliGeomManager::LoadGeometry();
  }

  if(!AliGeomManager::GetGeometry()) return Error(-1,"AliHLTTPCClusterTransformation::Init: Can not initialise geometry");
  
  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();
 
  if(!pCalib ) return Error(-2,"AliHLTTPCClusterTransformation::Init: Calibration not found");
  
  const AliMagF * field = (AliMagF*) TGeoGlobalMagField::Instance()->GetField();
  pCalib->SetExBField(field);

 
  if( !pCalib->GetTransform() ) return Error(-3,"AliHLTTPCClusterTransformation::Init: No TPC transformation found");
  
  // -- Get AliRunInfo variables  

  AliGRPObject tmpGRP, *pGRP=0;

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  
  if(!entry) return Error(-4,"AliHLTTPCClusterTransformation::Init: No GRP object found in data base");

  {
    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
      //cout<<"Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject"<<endl;
      m->Print();
      pGRP = &tmpGRP;
      pGRP->ReadValuesFromMap(m);
    }
    else {
      //cout<<"Found an AliGRPObject in GRP/GRP/Data, reading it"<<endl;
      pGRP = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
    }
  }
    
  if( !pGRP ){
    return Error(-5,"AliHLTTPCClusterTransformation::Init: Unknown format of the GRP object in data base");
  }

  AliRunInfo runInfo(pGRP->GetLHCState(),pGRP->GetBeamType(),pGRP->GetBeamEnergy(),pGRP->GetRunType(),pGRP->GetDetectorMask());
  AliEventInfo evInfo;
  evInfo.SetEventType(AliRawEventHeaderBase::kPhysicsEvent);

  entry=AliCDBManager::Instance()->Get("TPC/Calib/RecoParam");

  if(!entry) return Error(-6,"AliHLTTPCClusterTransformation::Init: No TPC reco param entry found in data base");

  TObject *recoParamObj = entry->GetObject();
  if(!recoParamObj) return Error(-7,"AliHLTTPCClusterTransformation::Init: Empty TPC reco param entry in data base");

  if (dynamic_cast<TObjArray*>(recoParamObj)) {
    //cout<<"\n\nSet reco param from AliHLTTPCClusterTransformation: TObjArray found \n"<<endl;
    TObjArray *copy = (TObjArray*)( static_cast<TObjArray*>(recoParamObj)->Clone() );
    fOfflineRecoParam.AddDetRecoParamArray(1,copy);
  }
  else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
    //cout<<"\n\nSet reco param from AliHLTTPCClusterTransformation: AliDetectorRecoParam found \n"<<endl;
    AliDetectorRecoParam *copy = (AliDetectorRecoParam*)static_cast<AliDetectorRecoParam*>(recoParamObj)->Clone();
    fOfflineRecoParam.AddDetRecoParam(1,copy);
  } else {    
    return Error(-8,"AliHLTTPCClusterTransformation::Init: Unknown format of the TPC Reco Param entry in the data base");
  }
  
  
  fOfflineRecoParam.SetEventSpecie(&runInfo, evInfo, 0);    
 
  // 

  AliTPCRecoParam* recParam = (AliTPCRecoParam*)fOfflineRecoParam.GetDetRecoParam(1);

  if( !recParam ) return Error(-9,"AliHLTTPCClusterTransformation::Init: No TPC Reco Param entry found for the given event specification");

 
  pCalib->GetTransform()->SetCurrentRecoParam(recParam);

  // set current time stamp and initialize the fast transformation
  int err = fFastTransform.Init( pCalib->GetTransform(), TimeStamp );

  if( err!=0 ){
    return Error(-10,Form( "AliHLTTPCClusterTransformation::Init: Initialisation of Fast Transformation failed with error %d :%s",err,fFastTransform.GetLastError()) );
  }

  return 0;
}


Int_t  AliHLTTPCClusterTransformation::Init( const AliHLTTPCFastTransformObject &obj )
{
  // Initialisation

  if(!AliGeomManager::GetGeometry()){
    AliGeomManager::LoadGeometry();
  }

  if(!AliGeomManager::GetGeometry()) return Error(-1,"AliHLTTPCClusterTransformation::Init: Can not initialise geometry");
  
  // set current time stamp and initialize the fast transformation
  int err = fFastTransform.ReadFromObject( obj );

  if( err!=0 ){
    return Error(-10,Form( "AliHLTTPCClusterTransformation::Init: Initialisation of Fast Transformation failed with error %d :%s",err,fFastTransform.GetLastError()) );
  }
  return(0);
}


Bool_t AliHLTTPCClusterTransformation::IsInitialised() const 
{
  // Is the transformation initialised
  return fFastTransform.IsInitialised();
}

void AliHLTTPCClusterTransformation::DeInit()
{
  // Deinitialisation
  fFastTransform.DeInit();
}

Int_t AliHLTTPCClusterTransformation::SetCurrentTimeStamp( Long_t TimeStamp )
{
  // Set the current time stamp  

  AliTPCRecoParam* recParam = (AliTPCRecoParam*)fOfflineRecoParam.GetDetRecoParam(1);
  if( !recParam )  return Error(-1,"AliHLTTPCClusterTransformation::SetCurrentTimeStamp: No TPC Reco Param entry found");

  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();
  if(!pCalib ) return Error(-2,"AliHLTTPCClusterTransformation::Init: Calibration not found");
   
  if( !pCalib->GetTransform() ) return Error(-3,"AliHLTTPCClusterTransformation::SetCurrentTimeStamp: No TPC transformation found");
  
  pCalib->GetTransform()->SetCurrentRecoParam(recParam);

  int err = fFastTransform.SetCurrentTimeStamp( TimeStamp );
  if( err!=0 ){
    return Error(-4,Form( "AliHLTTPCClusterTransformation::SetCurrentTimeStamp: SetCurrentTimeStamp to the Fast Transformation failed with error %d :%s",err,fFastTransform.GetLastError()) );
  }
  return 0;
}

void AliHLTTPCClusterTransformation::Print(const char* /*option*/) const
{
  // print info
  fFastTransform.Print();
}


Int_t AliHLTTPCClusterTransformation::GetSize() const
{
  // total size of the object
  int size = sizeof(AliHLTTPCClusterTransformation) - sizeof(AliHLTTPCFastTransform) + fFastTransform.GetSize();
  return size;
}
