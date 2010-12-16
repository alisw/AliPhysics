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
#include "AliHLTTPCTransform.h"

#include "AliTPCcalibDB.h"
#include "AliTPCTransform.h"
#include "AliTPCParam.h"
#include "AliTPCRecoParam.h"
#include "AliGeomManager.h"

ClassImp(AliHLTTPCClusterTransformation) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCClusterTransformation::AliHLTTPCClusterTransformation()
:
  fOfflineTransform(NULL),
  fOfflineTPCParam( NULL )
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt  
}

AliHLTTPCClusterTransformation::AliHLTTPCClusterTransformation(const AliHLTTPCClusterTransformation&)
:
  fOfflineTransform(NULL),
  fOfflineTPCParam( NULL )
{
  // copy constructor prohibited
}

AliHLTTPCClusterTransformation& AliHLTTPCClusterTransformation::operator=(const AliHLTTPCClusterTransformation&)
{
  // assignment operator prohibited
  return *this;
}

AliHLTTPCClusterTransformation::~AliHLTTPCClusterTransformation() 
{ 
  // see header file for class documentation
  delete fOfflineTransform;
}


int  AliHLTTPCClusterTransformation::Init( double FieldBz, UInt_t TimeStamp )
{
  // Initialisation

  delete fOfflineTransform;
  fOfflineTPCParam = 0;

  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();

  if(!pCalib ) return -1;

  pCalib->SetExBField(FieldBz);
  
  if(!AliGeomManager::GetGeometry()){
     AliGeomManager::LoadGeometry();
  }

  if( !pCalib->GetTransform() ) return -2; 

  fOfflineTransform = new AliTPCTransform (*pCalib->GetTransform());
  fOfflineTransform->SetCurrentRecoParam( AliTPCRecoParam::GetHLTParam() );
  fOfflineTransform->SetCurrentTimeStamp( TimeStamp );
  fOfflineTPCParam = pCalib->GetParameters(); 
  if( !fOfflineTPCParam ) return -3;

  fOfflineTPCParam->Update();
  fOfflineTPCParam->ReadGeoMatrices();

  return 0;
}


void AliHLTTPCClusterTransformation::SetCurrentTimeStamp( UInt_t TimeStamp )
{
  // Set the current time stamp
  if( fOfflineTransform ) fOfflineTransform->SetCurrentTimeStamp( TimeStamp );
}


int  AliHLTTPCClusterTransformation::Transform( int Slice, int Row, float Pad, float Time, float XYZ[] )
{
  // Convert row, pad, time to X Y Z
   	   
  Int_t sector=-99, thisrow=-99;
  AliHLTTPCTransform::Slice2Sector( Slice, Row, sector, thisrow);

  if( !fOfflineTransform ){   	   	   
    AliHLTTPCTransform::Raw2Local( XYZ, sector, thisrow, Pad, Time); 
    if(Slice>17) XYZ[1]= - XYZ[1];	   
    return 0;
  }

  Int_t iSector[1]= {sector};   
  Double_t x[3] = { thisrow, Pad, Time }; 
  fOfflineTransform->Transform(x,iSector,0,1);
  Double_t y[3]= {x[0],x[1],x[2]};

	  
  if( fOfflineTPCParam && sector<fOfflineTPCParam->GetNSector() ){
    TGeoHMatrix  *alignment = fOfflineTPCParam->GetClusterMatrix( sector );
    if ( alignment ) alignment->LocalToMaster( x, y);
  }	 

  XYZ[0] = y[0];
  XYZ[1] = y[1];
  XYZ[2] = y[2];

  return 0; 
}
