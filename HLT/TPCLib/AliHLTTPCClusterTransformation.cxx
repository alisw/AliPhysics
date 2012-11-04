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
#include "AliHLTTPCFastTransform.h"

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

using namespace std;

ClassImp(AliHLTTPCClusterTransformation) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCClusterTransformation::AliHLTTPCClusterTransformation()
:
  fOfflineTPCParam( NULL ),
  fOfflineRecoParam( NULL), 
  fLastSector(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt  
  fAliT[0] = 0.;
  fAliT[1] = 0.;
  fAliT[2] = 0.;
  SetRotationMatrix();
}
AliHLTTPCClusterTransformation::~AliHLTTPCClusterTransformation() 
{ 
  // see header file for class documentation
}


int  AliHLTTPCClusterTransformation::Init( double FieldBz, UInt_t TimeStamp )
{
  // Initialisation
  fOfflineTPCParam = 0;
  fLastSector = -1;
  fAliT[0] = 0.;
  fAliT[1] = 0.;
  fAliT[2] = 0.;
  
  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();

  if(!pCalib ) return -1;

  pCalib->SetExBField(FieldBz);
  if(!AliGeomManager::GetGeometry()){
    AliGeomManager::LoadGeometry();
  }

  if( !pCalib->GetTransform() ) return -2; 
  pCalib->GetTransform()->SetCurrentRecoParam(NULL);
  
  delete fOfflineRecoParam;
  fOfflineRecoParam = new AliRecoParam;
  if( !fOfflineRecoParam ) return -3;
  
  
  AliCDBEntry *entry=AliCDBManager::Instance()->Get("TPC/Calib/RecoParam");
  if(!entry) return -4;
  TObject *recoParamObj = entry->GetObject();
  if(!recoParamObj) return -5;
  if (dynamic_cast<TObjArray*>(recoParamObj)) {
    //cout<<"\n\nSet reco param from AliHLTTPCClusterTransformation: TObjArray found \n"<<endl;
    fOfflineRecoParam->AddDetRecoParamArray(1,dynamic_cast<TObjArray*>(recoParamObj));
  }
  else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
    //cout<<"\n\nSet reco param from AliHLTTPCClusterTransformation: AliDetectorRecoParam found \n"<<endl;
    fOfflineRecoParam->AddDetRecoParam(1,dynamic_cast<AliDetectorRecoParam*>(recoParamObj));
  }
  
  // -- Get AliRunInfo variables  

  AliGRPObject tmpGRP, *pGRP=0;

  entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  if(!entry) return -6;

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
    return -7;
  }

  AliRunInfo runInfo(pGRP->GetLHCState(),pGRP->GetBeamType(),pGRP->GetBeamEnergy(),pGRP->GetRunType(),pGRP->GetDetectorMask());
  AliEventInfo evInfo;
  evInfo.SetEventType(AliRawEventHeaderBase::kPhysicsEvent);
  
  fOfflineRecoParam->SetEventSpecie(&runInfo, evInfo, 0);    
 
  AliTPCRecoParam* recParam = (AliTPCRecoParam*)fOfflineRecoParam->GetDetRecoParam(1);

  // 
 
  pCalib->GetTransform()->SetCurrentRecoParam(recParam);

  fOfflineTPCParam = pCalib->GetParameters();
  if( !fOfflineTPCParam ) return -8;

  fOfflineTPCParam->Update();
  fOfflineTPCParam->ReadGeoMatrices();  

  SetRotationMatrix();

  // set current time stamp and initialize the fast transformation instance, if necessary

  if( !AliHLTTPCFastTransform::Instance()->IsInitialised() ) AliHLTTPCFastTransform::Instance()->Init( pCalib->GetTransform(), TimeStamp );
  else SetCurrentTimeStamp( TimeStamp );
  return 0;
}


void AliHLTTPCClusterTransformation::SetCurrentTimeStamp( UInt_t TimeStamp )
{
  // Set the current time stamp  
  AliHLTTPCFastTransform::Instance()->SetCurrentTimeStamp( TimeStamp );
}


int  AliHLTTPCClusterTransformation::Transform( int Slice, int Row, float Pad, float Time, float XYZ[] )
{
  // Convert row, pad, time to X Y Z
   	   
  Int_t sector=-99, thisrow=-99;  
  AliHLTTPCTransform::Slice2Sector( Slice, Row, sector, thisrow);
  Double_t x[3] = {0,0,0}; 
  AliHLTTPCFastTransform::Instance()->Transform(sector, thisrow, Pad, Time, x);

  if( sector!= fLastSector ){
    if( fOfflineTPCParam && sector<fOfflineTPCParam->GetNSector() ){
      TGeoHMatrix  *alignment = fOfflineTPCParam->GetClusterMatrix( sector );
      if ( alignment ){
	const Double_t *tr = alignment->GetTranslation();
	const Double_t *rot = alignment->GetRotationMatrix();
	if( tr && rot ){
	  for( int i=0; i<3; i++ ) fAliT[i] = tr[i];
	  for( int i=0; i<9; i++ ) fAliR[i] = rot[i];
	}
      }
    } else {
      fAliT[0] = 0.;
      fAliT[1] = 0.;
      fAliT[2] = 0.;
      for( int i=0; i<9; i++ ) fAliR[i] = 0;
      fAliR[0] = 1.;
      fAliR[4] = 1.;
      fAliR[8] = 1.;
    }
    fLastSector = sector;
  }
  // alignment->LocalToMaster( x, y);

  XYZ[0] = fAliT[0] + x[0]*fAliR[0] + x[1]*fAliR[1] + x[2]*fAliR[2];
  XYZ[1] = fAliT[1] + x[0]*fAliR[3] + x[1]*fAliR[4] + x[2]*fAliR[5];
  XYZ[2] = fAliT[2] + x[0]*fAliR[6] + x[1]*fAliR[7] + x[2]*fAliR[8];

  return 0; 
}

int  AliHLTTPCClusterTransformation::ReverseAlignment( float XYZ[], int slice, int padrow)
{
  // reverse the alignment correction
  Int_t sector=-99, thisrow=-99;
  AliHLTTPCTransform::Slice2Sector( slice, padrow, sector, thisrow);
  if( sector!= fLastSector ){
    if( fOfflineTPCParam && sector<fOfflineTPCParam->GetNSector() ){
      TGeoHMatrix  *alignment = fOfflineTPCParam->GetClusterMatrix( sector );
      if ( alignment ){
	const Double_t *tr = alignment->GetTranslation();
	const Double_t *rot = alignment->GetRotationMatrix();
	if(tr){
	  for( int i=0; i<3; i++ ) fAliT[i] = tr[i];
	}
	SetRotationMatrix(rot, true);
      }
    } else {
      fAliT[0] = 0.;
      fAliT[1] = 0.;
      fAliT[2] = 0.;
      SetRotationMatrix(NULL, true);
    }
    fLastSector = sector;
  }

  // correct for alignment: translation
  float xyz[3];
  xyz[0] = XYZ[0] - fAliT[0];
  xyz[1] = XYZ[1] - fAliT[1];
  xyz[2] = XYZ[2] - fAliT[2];

  // correct for alignment: rotation
  XYZ[0]=xyz[0]*fAdjR[0] + xyz[1]*fAdjR[1] + xyz[2]*fAdjR[2];
  XYZ[1]=xyz[0]*fAdjR[3] + xyz[1]*fAdjR[4] + xyz[2]*fAdjR[5];
  XYZ[2]=xyz[0]*fAdjR[6] + xyz[1]*fAdjR[7] + xyz[2]*fAdjR[8];

  return 0;
}

void AliHLTTPCClusterTransformation::SetRotationMatrix(const Double_t *rot, bool bCalcAdjugate)
{
  // set the rotation matrix and calculate the adjugate if requested
  if (rot) {
    for( int i=0; i<9; i++ ) fAliR[i] = rot[i];
    if (bCalcAdjugate) {
      CalcAdjugateRotation();
    }
    return;
  }
  for( int i=0; i<9; i++ ) {fAliR[i] = 0; fAdjR[i] = 0;}
  fAliR[0] = 1.;
  fAliR[4] = 1.;
  fAliR[8] = 1.; 
  fAdjR[0] = 1.;
  fAdjR[4] = 1.;
  fAdjR[8] = 1.; 
}

bool AliHLTTPCClusterTransformation::CalcAdjugateRotation(bool bCheck)
{
  // check rotation matrix and adjugate for consistency
  fAdjR[0]= fAliR[4]*fAliR[8]-fAliR[5]*fAliR[7];
  fAdjR[1]= fAliR[5]*fAliR[6]-fAliR[3]*fAliR[8];
  fAdjR[2]= fAliR[3]*fAliR[7]-fAliR[4]*fAliR[6];

  fAdjR[3]= fAliR[2]*fAliR[7]-fAliR[1]*fAliR[8];
  fAdjR[4]= fAliR[0]*fAliR[8]-fAliR[2]*fAliR[6];
  fAdjR[5]= fAliR[2]*fAliR[6]-fAliR[0]*fAliR[7];

  fAdjR[6]= fAliR[1]*fAliR[5]-fAliR[2]*fAliR[4];
  fAdjR[7]= fAliR[2]*fAliR[3]-fAliR[0]*fAliR[5];
  fAdjR[8]= fAliR[0]*fAliR[4]-fAliR[1]*fAliR[3];

  if (bCheck) {
    for (int r=0; r<3; r++) {
      for (int c=0; c<3; c++) {
	float a=0.;
	float expected=0.;
	if (r==c) expected=1.;
	for (int i=0; i<3; i++) {
	  a+=fAliR[3*r+i]*fAdjR[c+(3*i)];
	}
	if (TMath::Abs(a-expected)>0.00001) {
	  std::cout << "inconsistent adjugate at " << r << c << ": " << a << " " << expected << std::endl;
	  return false;
	}
      }
    }
  }
  return true;
}

void AliHLTTPCClusterTransformation::Print(const char* /*option*/) const
{
  // print info
  ios::fmtflags coutflags=std::cout.flags(); // backup cout status flags
  std::cout << "AliHLTTPCClusterTransformation for sector " << fLastSector << std::endl;

  std::cout.setf(ios_base::showpos|ios_base::showpos|ios::right);
  std::cout << "  translation: " << std::endl;
  int r=0;
  for (r=0; r<3; r++) {
    std::cout << setw(7) << fixed << setprecision(2);
    cout << "  " << fAliT[r] << std::endl;
  }
  std::cout << "  rotation and adjugated rotation: " << std::endl;
  for (r=0; r<3; r++) {
    int c=0;
    std::cout << setw(7) << fixed << setprecision(2);
    for (c=0; c<3; c++) std::cout << "  " << fAliR[3*r+c];
    std::cout << "      ";
    for (c=0; c<3; c++) std::cout << "  " << fAdjR[3*r+c];
    std::cout << endl;
  }
  std::cout.flags(coutflags); // restore the original flags
}
