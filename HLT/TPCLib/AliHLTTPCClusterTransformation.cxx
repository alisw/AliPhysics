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

#ifdef HAVE_ALIGPU
#include "TPCFastTransform.h"
#include "TPCFastTransformManager.h"
#endif

using namespace std;

ClassImp(AliHLTTPCClusterTransformation) //ROOT macro for the implementation of ROOT specific class methods

AliRecoParam AliHLTTPCClusterTransformation::fOfflineRecoParam;

AliHLTTPCClusterTransformation::AliHLTTPCClusterTransformation()
:
  fError(),
  fTransformKind( TransformOldFastTransform ),
  fOrigTransform(NULL),
  fFastTransform(),
#ifdef HAVE_ALIGPU
  fFastTransformIRS(new GPUCA_NAMESPACE::gpu::TPCFastTransform),
  fFastTransformManager( new GPUCA_NAMESPACE::gpu::TPCFastTransformManager ),
#else
  fFastTransformIRS(NULL),
  fFastTransformManager(NULL),
#endif
  fIsMC(0)
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
#ifdef HAVE_ALIGPU
  delete fFastTransformManager;
  delete fFastTransformIRS;
#endif
}

int  AliHLTTPCClusterTransformation::Init( double /*FieldBz*/, Long_t TimeStamp, bool isMC, int useOrigTransform )
{
  TransformationKind kind = TransformOldFastTransform;
  if( useOrigTransform ) kind = TransformOriginal;
  return Init( TimeStamp, isMC, kind );
}


int  AliHLTTPCClusterTransformation::Init( Long_t TimeStamp, bool isMC, TransformationKind transformKind )
{
  // Initialisation

  fTransformKind = transformKind;

  fIsMC = isMC;
 
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
  
  if (fIsMC && !recParam->GetUseCorrectionMap()) TimeStamp = 0;
 
  pCalib->GetTransform()->SetCurrentRecoParam(recParam);
  

  // set current time stamp and initialize the fast transformation
  bool useOrigTransform = ( fTransformKind==TransformOriginal);
  int err = fFastTransform.Init( pCalib->GetTransform(), TimeStamp, useOrigTransform );
  if( err!=0 ){
    return Error(-10,Form( "AliHLTTPCClusterTransformation::Init: Initialisation of Fast Transformation failed with error %d :%s",err,fFastTransform.GetLastError()) );
  }

#ifdef HAVE_ALIGPU
  if( fTransformKind == TransformFastIRS ){
    err = fFastTransformManager->create( *fFastTransformIRS, pCalib->GetTransform(), TimeStamp );
    if( err!=0 ){
      return Error(-10,Form( "AliHLTTPCClusterTransformation::Init: Initialisation of Fast Transformation failed with error %d :%s",err,fFastTransformManager->getLastError()) );
    }
  }
#endif
  
  // offline transformation is already initialised in fFastTransform or in fFastTransformManager
  fOrigTransform = pCalib->GetTransform();

  return 0;
}


Int_t  AliHLTTPCClusterTransformation::Init( const AliHLTTPCFastTransformObject &obj )
{
  // Initialisation
 
  fTransformKind = TransformOldFastTransform;
 
  if(!AliGeomManager::GetGeometry()){
    AliGeomManager::LoadGeometry();
  }

  if(!AliGeomManager::GetGeometry()) return Error(-1,"AliHLTTPCClusterTransformation::Init: Can not initialise geometry");
  
  // Load the fast transform object
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
  
  if (fIsMC && !recParam->GetUseCorrectionMap()) TimeStamp = 0;

  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();
  if(!pCalib ) return Error(-2,"AliHLTTPCClusterTransformation::Init: Calibration not found");
   
  if( !pCalib->GetTransform() ) return Error(-3,"AliHLTTPCClusterTransformation::SetCurrentTimeStamp: No TPC transformation found");
  
  pCalib->GetTransform()->SetCurrentRecoParam(recParam);

  int err = fFastTransform.SetCurrentTimeStamp( TimeStamp );
  if( err!=0 ){
    return Error(-4,Form( "AliHLTTPCClusterTransformation::SetCurrentTimeStamp: SetCurrentTimeStamp to the Fast Transformation failed with error %d :%s",err,fFastTransform.GetLastError()) );
  }

#ifdef HAVE_ALIGPU
  if( fTransformKind == TransformFastIRS ){
    err = fFastTransformManager->updateCalibration( *fFastTransformIRS, TimeStamp );
    if( err!=0 ){
      return Error(-10,Form( "AliHLTTPCClusterTransformation::SetCurrentTimeStamp: Initialisation of Fast Transformation failed with error %d :%s",err,fFastTransformManager->getLastError()) );
    }
  }
#endif

  // time stamp for the offline transformation is already set in fFastTransform or in fFastTransformManager
  fOrigTransform = pCalib->GetTransform();

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


Int_t  AliHLTTPCClusterTransformation::Transform( int Slice, int Row, float Pad, float Time, float XYZ[] )
{
  //static int st=0;
  //st++;
#ifdef HAVE_ALIGPU
  if( fTransformKind == TransformFastIRS ){
    //if( st<100) cout<<"IRS"<<endl;
    fFastTransformIRS->Transform( Slice, Row, Pad, Time, XYZ[0], XYZ[1], XYZ[2]);
  } else
#endif
  if( fTransformKind == TransformOriginal ){
    //if( st<100) cout<<"orig"<<endl;
    Int_t sector=-99, thisrow=-99;
    AliHLTTPCGeometry::Slice2Sector( Slice, Row, sector, thisrow);
    Int_t is[]={sector};
    Double_t xx[]={static_cast<Double_t>(thisrow),Pad,Time};
    fOrigTransform->Transform(xx,is,0,1);
    for (int i = 0;i < 3;i++) XYZ[i] = xx[i];
  } else {
    //if( st<100) cout<<"old"<<endl;
    // Convert row, pad, time to X Y Z
    Int_t sector=-99, thisrow=-99;
    AliHLTTPCGeometry::Slice2Sector( Slice, Row, sector, thisrow);
    int err = fFastTransform.Transform(sector, thisrow, Pad, Time, XYZ);
    if( err!=0 ) return Error(-1,Form( "AliHLTTPCClusterTransformation::Transform: Fast Transformation failed with error %d :%s",err,fFastTransform.GetLastError()) );
  }
  
  return 0;
}
