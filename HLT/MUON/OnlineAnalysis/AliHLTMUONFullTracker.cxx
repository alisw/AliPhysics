/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////
///Author : Indranil Das, SINP, INDIA
///
///Email :  indra.das@saha.ac.in || indra.ehep@gmail.com
///
/// This class is created to peform the a full track reconstroction
/// for online HLT for Dimuon Detector. It is based on the method of
/// straight line tracking for slat detectors and cellular automation
/// mehtods for finding the tracks in the quadrant detectos. The track
/// segments of the front and back side of the spectromrter is either
/// compared using KalmanChi2 method (a slightly modified version of
/// AliMUONTrackReconstructorK) or alternated method of simple 
/// extrapolation of tracks through dipole magnet. Finally the Track is
/// passed through Absorber to take the MCS effect and Branson Effect
/// into account for correction of track parameters.
/////////////////////////////////////////////////////////////////////////

//ROOT Includes
#include "TVector3.h"
#include "TGeoGlobalMagField.h"

//STEER Includes
#include "AliMagF.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"

//MUON Includes
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTrackParam.h"

//MUON mappinf includes
#include "AliMpDEIterator.h"

//HLT Includes
#include "AliHLTLogging.h"

//HLT/MUON Includes
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONDataTypes.h"
#include "AliHLTMUONUtils.h"

//HLT/MUON/OnlineAnalysis includes
#include "AliHLTMUONFullTracker.h"

using namespace std;

class	AliRunInfo;
class	AliLog;
class	AliCDBEntry;
class	AliMpCDB;
class	AliMpSegmentation;
class	AliMpDDLStore;
class	AliGRPObject;
class	TString;
class	TMap;

//#define PRINT_FULL 1

#ifdef PRINT_FULL
#define PRINT_POINTS 1
#define PRINT_BACK 1
#define PRINT_FRONT 1
#define PRINT_KALMAN 1
#define PRINT_SELECT 1
#define PRINT_OUTPUT 1
#define PRINT_TRACKSEG 1
///#define PRINT_DETAIL_KALMAN 1
#endif

class AliHLTMUONFullTracker;

const Float_t AliHLTMUONFullTracker::fgkTrackDetCoordinate[3] = {
  155.179+20.0,  166.234+20.0, 
  (AliMUONConstants::DefaultChamberZ(4)+ AliMUONConstants::DefaultChamberZ(5))/2.0
};

const Double_t AliHLTMUONFullTracker::fgkAbsoedge[4] = {92.0, 317.0,443.0,499.0};
const Double_t AliHLTMUONFullTracker::fgkRadLen[3] = {13.875635, 11.273801, 1.765947};
const Double_t AliHLTMUONFullTracker::fgkRho[3] = {1.750000, 2.350000, 7.880000};
const Double_t AliHLTMUONFullTracker::fgkAtomicZ[3] = {6.000000, 11.098486,25.780000};
const Double_t AliHLTMUONFullTracker::fgkAtomicA[3] = {12.010000, 22.334156,55.299670 };


const Int_t AliHLTMUONFullTracker::fgkMaxNofCellsPerCh = 1500;
const Int_t AliHLTMUONFullTracker::fgkMaxNofPointsPerCh = 600; 
const Int_t AliHLTMUONFullTracker::fgkMaxNofCh = 11 ; /// 10tracking + 1 trigger
const Int_t AliHLTMUONFullTracker::fgkMaxNofTracks = 200;
const Int_t AliHLTMUONFullTracker::fgkMaxNofConnectedTracks = 20;



AliHLTMUONFullTracker::AliHLTMUONFullTracker() : 
  AliHLTLogging(),
  fChamberGeometryTransformer(0x0),
  fChPoint(0x0),
  fChPoint11(0x0),
  fBackTrackSeg(0x0),
  fFrontTrackSeg(0x0),
  fExtrapSt3X(0x0),
  fExtrapSt3Y(0x0),
  fInclinationBack(0x0),
  fNofConnectedfrontTrackSeg(0x0),
  fBackToFront(0x0),
  fNofPoints(0x0),
  fTrackParam(0x0),
  fHalfTrack(0x0),
  fTotNofPoints(0),
  fTotTrackSeg(0),
  fNofbackTrackSeg(0),
  fNoffrontTrackSeg(0),
  fNofConnected(0),
  fNofTracks(0),
  fDetElemList(),
  fFastTracking(kFALSE)
{

  /// Constructor of the class

  /// fChPoint = (AliHLTMUONRecHitStruct***)(malloc(10 * sizeof(AliHLTMUONRecHitStruct)));
  /// for( Int_t ich=0;ich<10;ich++)
  ///   fChPoint[ich] = (AliHLTMUONRecHitStruct**)(malloc(300 * sizeof(AliHLTMUONRecHitStruct)));

  try{
    fChPoint = new AliHLTMUONRecHitStruct**[fgkMaxNofCh-1];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fChPoint");
    throw;
  }

  for( Int_t ich=0;ich<(fgkMaxNofCh-1);ich++)
    try{
      fChPoint[ich] = new AliHLTMUONRecHitStruct*[fgkMaxNofPointsPerCh];
    }catch (const std::bad_alloc&){
      HLTError("Dynamic memory allocation failed in constructor : fChPoint[%d]",ich);
      throw;
    }
  
  try{
    fChPoint11 = new AliHLTMUONTriggerRecordStruct*[fgkMaxNofPointsPerCh];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fChPoint11");
    throw;
  }

  try{
    fBackTrackSeg = new TrackSeg[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fBackTrackSeg");
    throw;
  }

  try{
    fFrontTrackSeg = new TrackSeg[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fFrontTrackSeg");
    throw;
  }

  try{
    fExtrapSt3X = new float[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fExtrapSt3X");
    throw;
  }

  try{
    fExtrapSt3Y = new float[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fExtrapSt3Y");
    throw;
  }

  try{
    fInclinationBack = new float[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fInclinationBack");
    throw;
  }
  
  try{
    fNofConnectedfrontTrackSeg = new int[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fNofConnectedfrontTrackSeg");
    throw;
  }

  try{
    fBackToFront = new int*[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fBackToFront");
    throw;
  }

  for( Int_t itr=0;itr<fgkMaxNofTracks;itr++)
    try{
      fBackToFront[itr] = new int[fgkMaxNofConnectedTracks];
    }catch (const std::bad_alloc&){
      HLTError("Dynamic memory allocation failed in constructor : fBackToFront[%d]",itr);
      throw;
    }

  
  try{
    fNofPoints = new int[fgkMaxNofCh];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fNofPoints");
    throw;
  }
  
  try{
    fTrackParam = new AliMUONTrackParam[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fTrackParam");
    throw;
  }

  try{
    fHalfTrack  = new HalfTrack[fgkMaxNofTracks];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fHalfTrack");
    throw;
  }
  
  Clear();

}

///__________________________________________________________________________

AliHLTMUONFullTracker::~AliHLTMUONFullTracker()
{ 
  /// Destructor of the class

  ///delete fChamberGeometryTransformer;
  ///free(fChPoint);
  delete []fChPoint;
  delete []fChPoint11;
  delete []fBackTrackSeg;
  delete []fFrontTrackSeg;
  delete []fExtrapSt3X;
  delete []fExtrapSt3Y;
  delete []fInclinationBack;
  delete []fNofConnectedfrontTrackSeg;
  delete []fBackToFront;
  delete []fNofPoints;
  delete []fTrackParam;
  delete []fHalfTrack;
  
  fChamberGeometryTransformer->Delete();

  fDetElemList.clear();  

}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::Clear(){
  
  /// Clear method to be called after each event
  
#ifdef PRINT_TRACKSEG
  for(int iFrontTrackSeg=0;iFrontTrackSeg<fNoffrontTrackSeg;iFrontTrackSeg++)
    for(int ipoint=0;ipoint<4;ipoint++)
      if(fFrontTrackSeg[iFrontTrackSeg].fIndex[ipoint]!=-1)
	printf("FrontTrackSeg : %d\t%f\t%f\t%f\n",iFrontTrackSeg,
	       fChPoint[ipoint][fFrontTrackSeg[iFrontTrackSeg].fIndex[ipoint]]->fX,
	       fChPoint[ipoint][fFrontTrackSeg[iFrontTrackSeg].fIndex[ipoint]]->fY,
	       fChPoint[ipoint][fFrontTrackSeg[iFrontTrackSeg].fIndex[ipoint]]->fZ); 
  printf("\n\n");
  for(int iBackTrackSeg=0;iBackTrackSeg<fNofbackTrackSeg;iBackTrackSeg++)
    for(int ipoint=0;ipoint<4;ipoint++)
      if(fBackTrackSeg[iBackTrackSeg].fIndex[ipoint]!=-1)
	printf("BackTrackSeg : %d\t%f\t%f\t%f, nofFront : %d\n",iBackTrackSeg,
	       fChPoint[ipoint+6][fBackTrackSeg[iBackTrackSeg].fIndex[ipoint]]->fX,
	       fChPoint[ipoint+6][fBackTrackSeg[iBackTrackSeg].fIndex[ipoint]]->fY,
	       fChPoint[ipoint+6][fBackTrackSeg[iBackTrackSeg].fIndex[ipoint]]->fZ,fNofConnectedfrontTrackSeg[iBackTrackSeg]); 
#endif
  
  for( Int_t ich=0;ich<fgkMaxNofCh;ich++)
    fNofPoints[ich] = 0;
  
  fNofbackTrackSeg = 0;
  fNoffrontTrackSeg = 0;
  fNofConnected = 0;
  fNofTracks = 0;
  
  return true;
}

///__________________________________________________________________________

void AliHLTMUONFullTracker::Print()
{
  /// Print anything mainly for debugging the codes, will be removed later
  HLTInfo("\nPrint This--->\n");
}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::Init()
{
  /// Initilation to be called once, later can be used to set/load the CDB path/entries
  HLTInfo("path : %s, run number : %d",
	  (AliCDBManager::Instance())->GetDefaultStorage()->GetURI().Data(),(AliCDBManager::Instance())->GetRun());
  if (AliGeomManager::GetGeometry() == NULL){
    AliGeomManager::LoadGeometry();
    AliGeomManager::ApplyAlignObjsFromCDB("GRP MUON");
  }
  
  Double_t b[3], x[3];
  x[0] = 0.0 ; x[1] = 0.0 ; x[2] = -950.0;

  TGeoGlobalMagField::Instance()->Field(x,b);
  //The following condition is based on the fact that at the middle of the dipole the field cannot be zero
  if(TMath::AreEqualAbs(b[0],0.0,1.0e-5) and TMath::AreEqualAbs(b[1],0.0,1.0e-5) and TMath::AreEqualAbs(b[2],0.0,1.0e-5)){
    HLTWarning("Magnetic field is not set by GRP");
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));
    TGeoGlobalMagField::Instance()->Field(x,b);
  }
  
  
  HLTWarning("At (X,Y,Z) : (%6.2lf,%6.2lf,%6.2lf) Field (Bx,By,Bz) is (%6.2lf,%6.2lf,%6.2lf)",
	     x[0],x[1],x[2],b[0],b[1],b[2]);

  AliMUONTrackExtrap::SetField();
  /// see header file for class documentation
  fChamberGeometryTransformer = new AliMUONGeometryTransformer();
  if(! fChamberGeometryTransformer->LoadGeometryData()){
    HLTError("Failed to Load Geomerty Data ");
  }
  
  fDetElemList.clear();  
  for(int ich=0;ich<AliMUONConstants::NCh();ich++){
    AliMpDEIterator it;
    for ( it.First(ich); ! it.IsDone(); it.Next() )
      {
	Int_t detElemId = it.CurrentDEId();
	fDetElemList[detElemId] = detElemId;
      }
  }//chamber loop
  return true;
}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::SetInput(AliHLTInt32_t /*ddl*/, const AliHLTMUONRecHitStruct  *data, AliHLTInt32_t size)
{
  /// Set the input for rechit points
  if(int(size)>=fgkMaxNofPointsPerCh/2)
    size = 0;
  
  AliHLTUInt16_t detElemID;
  AliHLTUInt8_t chamber;
  
#ifdef PRINT_POINTS  
  printf("Received from DDL : %d, nofHits : %d, data : %p\n",ddl,size,data);
#endif
  for( Int_t ipoint=0;ipoint<int(size);ipoint++){
    if(!data){
      HLTError("Null Data pointer from HitRec");
      Clear();
      return false;
    }
    
    AliHLTMUONUtils::UnpackRecHitFlags(data->fFlags,chamber,detElemID);

    if((not fDetElemList[detElemID]) or (chamber<0 ) or (chamber>=AliMUONConstants::NTrackingCh())){
      HLTDebug("Invalid tracking detelem : %d or chamber : %d",detElemID,chamber);
      continue;
    }
    
#ifdef PRINT_POINTS  
    printf("ch : %02d, detelem : %04d, (X,Y,Z) : (%8.3f,%8.3f,%8.3f)\n",chamber,detElemID,data->fX,data->fY,data->fZ);
#endif
    fChPoint[detElemID/100-1][fNofPoints[detElemID/100-1]++]  = (AliHLTMUONRecHitStruct  *)data;
    data++;
  }

  return true;
}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::SetInput(AliHLTInt32_t /*ddl*/, const AliHLTMUONTriggerRecordStruct  *data, AliHLTInt32_t size)
{
  /// Set the input for trigrecs

  if(int(size)>=fgkMaxNofPointsPerCh/2) 
    size = 0;

  AliHLTUInt16_t detElemID;
  AliHLTUInt8_t chamber;

  for( Int_t ipoint=0;ipoint<int(size);ipoint++){
    if(!data){
      HLTError("Null Data pointer from TrigRec");
      Clear();
      return false;
    }
    fChPoint11[fNofPoints[10]++] = (AliHLTMUONTriggerRecordStruct *)data;
    for( Int_t ich=0;ich<4;ich++){
      AliHLTMUONUtils::UnpackRecHitFlags((data->fHit[ich]).fFlags,chamber,detElemID);
      if((not fDetElemList[detElemID]) or (chamber<AliMUONConstants::NTrackingCh()) or (chamber > AliMUONConstants::NCh()) ){
	HLTDebug("Invalid trigger detelem : %d or chamber : %d",detElemID,chamber);
	continue;
      }
#ifdef PRINT_POINTS  
      printf("size : %d, itrig  : %04d, ch : %02d, detelem : %04d, (X,Y,Z) : (%8.3f,%8.3f,%8.3f)\n",
	     size,ipoint,chamber,detElemID,(data->fHit[ich]).fX,(data->fHit[ich]).fY,(data->fHit[ich]).fZ);
#endif
    }///ich loop
    data++;
  }///ipoint
  return true;
}

 
///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::Run( Int_t iEvent,AliHLTMUONTrackStruct *data, AliHLTUInt32_t& size)
{
  /// Main Run call of the class
  if(iEvent<0)
    {
      size = 0;
      return true;
    }
  bool resultOk = true;

  resultOk = SlatTrackSeg();
  if(not resultOk){
    HLTDebug("Error happened in tracking through slat chambers, this event will be skipped");
  }
  HLTDebug("Finishing SlatTrackSeg");

  if(resultOk){
    resultOk = PrelimMomCalc();
    if(not resultOk){
      HLTDebug("Error happened in calculating preliminary momentum, this event will be skipped");
    }
  }
  HLTDebug("Finishing PrelimMomCalc");
  
  if(resultOk){
    resultOk = QuadTrackSeg();
    if(not resultOk){
      HLTDebug("Error happened in tracking through quadrant chambers, this event will be skipped");
    }
  }
  HLTDebug("Finishing QuadTrackSeg");


  if(resultOk){
    if(fFastTracking)
      resultOk = SelectFront();
    else
      resultOk = KalmanChi2Test();

    if(not resultOk){
      HLTDebug("Error happened in tracking through in Kalman Chi2 checking, this event will be skipped");
    }
  }
  HLTDebug("Finishing KalmanChi2Test");

  if(resultOk){
    resultOk = ExtrapolateToOrigin();
    if(not resultOk){
      HLTDebug("Error happened in tracking extrapolation, this event will be skipped");
    }
  }
  HLTDebug("Finishing ExtrapolateToOrigin");

  if(resultOk){
    resultOk = FillOutData(data,size);
    if(not resultOk){
      HLTDebug("Error happened in filling the output data, this event will be skipped");
    }
  }
  HLTDebug("Finishing FillOutData");
  
  if(!resultOk)
    size = 0;

  HLTDebug("\niEvent: %d, has tracks : %d, triggers : %d, nof slat tracks : %d, quad tracks : %d, connected : %d\n",
	   iEvent,size,fNofPoints[10],fNofbackTrackSeg,fNoffrontTrackSeg,fNofConnected);
  Clear();
  
  return resultOk;

}
 
///__________________________________________________________________________

void AliHLTMUONFullTracker::Sub(const AliHLTMUONRecHitStruct *v1, const AliHLTMUONRecHitStruct *v2, AliHLTMUONRecHitStruct *v3) const
{
  /// Subtraction of position co-odinate of two space points

  v3->fX = v1->fX - v2->fX; 
  v3->fY = v1->fY - v2->fY; 
  v3->fZ = v1->fZ - v2->fZ;
};

///__________________________________________________________________________

Double_t AliHLTMUONFullTracker::Angle(const AliHLTMUONRecHitStruct *v1, const AliHLTMUONRecHitStruct *v2)
{
  ///Angle of a straight line formed using v1 and v2

  Float_t ptot2 = ((v1->fX * v1->fX) + (v1->fY * v1->fY) + (v1->fZ * v1->fZ))*
    ((v2->fX * v2->fX) + (v2->fY * v2->fY) + (v2->fZ * v2->fZ));
  if(ptot2 <= 0) {
    return 0.0;
  } else {
    Float_t arg = ((v1->fX * v2->fX) + (v1->fY * v2->fY) + (v1->fZ * v2->fZ))/sqrt(ptot2);
    if(arg >  1.0) arg =  1.0;
    if(arg < -1.0) arg = -1.0;
    return TMath::ACos(arg);
    ///return acos(arg);
  }
  
}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::FillOutData(AliHLTMUONTrackStruct *track, AliHLTUInt32_t& size)
{
  ///Fill the output data pointers

  size = (AliHLTUInt32_t(fNofbackTrackSeg)<size) ? AliHLTUInt32_t(fNofbackTrackSeg) : size;

  Bool_t hitset[16];
  for( Int_t ibackTrackSeg=0;ibackTrackSeg<int(size);ibackTrackSeg++){

    if(fNofConnectedfrontTrackSeg[ibackTrackSeg]>0){

    
      if(not TMath::Finite(fTrackParam[ibackTrackSeg].Px()) 
	 || not TMath::Finite(fTrackParam[ibackTrackSeg].Py()) 
	 || not TMath::Finite(fTrackParam[ibackTrackSeg].Pz())) continue; 

#ifdef PRINT_OUTPUT
      printf("\nsize : %d, itrack  : %04d, sign : %2d, Pt : %8.3f, (Px,Py,Pz) : (%8.3f,%8.3f,%8.3f)\n",
	     size,ibackTrackSeg,Int_t(TMath::Sign(1.,fTrackParam[ibackTrackSeg].GetInverseBendingMomentum())),
	     TMath::Sqrt(fTrackParam[ibackTrackSeg].Px()*fTrackParam[ibackTrackSeg].Px() +
			 fTrackParam[ibackTrackSeg].Py()*fTrackParam[ibackTrackSeg].Py()),
	     fTrackParam[ibackTrackSeg].Px(),fTrackParam[ibackTrackSeg].Py(),
	     fTrackParam[ibackTrackSeg].Pz());
#endif

   
      // Bits 8 and 9 must be kept zero to prevent the track ID from conflicting with other tracker components.
      track->fId = (ibackTrackSeg << 10) | (fBackTrackSeg[ibackTrackSeg].fTrigRec & 0xFF);
      track->fTrigRec = fBackTrackSeg[ibackTrackSeg].fTrigRec;
      track->fPx = fTrackParam[ibackTrackSeg].Px();
      track->fPy = fTrackParam[ibackTrackSeg].Py();
      track->fPz = fTrackParam[ibackTrackSeg].Pz();
      
      track->fChi2 = 0;
      
      track->fInverseBendingMomentum = fTrackParam[ibackTrackSeg].GetInverseBendingMomentum();
      track->fThetaY = TMath::Tan(fTrackParam[ibackTrackSeg].GetBendingSlope());
      track->fThetaX = TMath::Tan(fTrackParam[ibackTrackSeg].GetNonBendingSlope());
      
      track->fZ = fTrackParam[ibackTrackSeg].GetZ();
      track->fY = fTrackParam[ibackTrackSeg].GetBendingCoor();
      track->fX = fTrackParam[ibackTrackSeg].GetNonBendingCoor();
      

      for( Int_t ipoint=15;ipoint>=0;ipoint--){
	track->fHit[ipoint] = AliHLTMUONConstants::NilRecHitStruct();
	hitset[ipoint] = false;
	if(ipoint >= 6 and ipoint <= 9 and fBackTrackSeg[ibackTrackSeg].fIndex[ipoint-6]!=-1 ){
	    track->fHit[ipoint] = *(fChPoint[ipoint][fBackTrackSeg[ibackTrackSeg].fIndex[ipoint-6]]);
	    hitset[ipoint] = true;
	}else if(ipoint <= 3 and fFrontTrackSeg[fBackToFront[ibackTrackSeg][0]].fIndex[ipoint]!=-1 ){
	    track->fHit[ipoint] = *(fChPoint[ipoint][fFrontTrackSeg[fBackToFront[ibackTrackSeg][0]].fIndex[ipoint]]);
	    hitset[ipoint] = true;
	}
      }
      AliHLTMUONParticleSign sign = AliHLTMUONParticleSign(Int_t(TMath::Sign(1.,fTrackParam[ibackTrackSeg].GetInverseBendingMomentum())));
      track->fFlags = AliHLTMUONUtils::PackTrackFlags(sign,hitset);
      
      track++;
      fNofTracks++;
    
    }else{
      
      
      if(not TMath::Finite(fHalfTrack[ibackTrackSeg].fPx)
	 || not TMath::Finite(fHalfTrack[ibackTrackSeg].fPy)
	 || not TMath::Finite(fHalfTrack[ibackTrackSeg].fPz)) continue;
      
#ifdef PRINT_OUTPUT
      printf("\nsize : %d, itrack  : %04d, sign : %2d, Pt : %8.3f, (Px,Py,Pz) : (%8.3f,%8.3f,%8.3f)\n",
	     size,ibackTrackSeg,Int_t(TMath::Sign(1.,fTrackParam[ibackTrackSeg].GetInverseBendingMomentum())),
	     TMath::Sqrt(fHalfTrack[ibackTrackSeg].fPx*fHalfTrack[ibackTrackSeg].fPx +
			 fHalfTrack[ibackTrackSeg].fPy*fHalfTrack[ibackTrackSeg].fPy),
	     fHalfTrack[ibackTrackSeg].fPx,fHalfTrack[ibackTrackSeg].fPy,
	     fHalfTrack[ibackTrackSeg].fPz);
#endif
      
      
      // Bits 8 and 9 must be kept zero to prevent the track ID from conflicting with other tracker components.
      track->fId = (ibackTrackSeg << 10) | (fBackTrackSeg[ibackTrackSeg].fTrigRec & 0xFF);
      track->fTrigRec = fBackTrackSeg[ibackTrackSeg].fTrigRec;
      track->fPx = fHalfTrack[ibackTrackSeg].fPx;
      track->fPy = fHalfTrack[ibackTrackSeg].fPy;
      track->fPz = fHalfTrack[ibackTrackSeg].fPz;
      
      track->fChi2 = 0;
      
      
      
      for( Int_t ipoint=15;ipoint>=0;ipoint--){
	track->fHit[ipoint] = AliHLTMUONConstants::NilRecHitStruct();
	hitset[ipoint] = false;

	if(ipoint<=9 and ipoint>=6 and fBackTrackSeg[ibackTrackSeg].fIndex[ipoint-6]!=-1 ){
	  
	  track->fHit[ipoint] = *(fChPoint[ipoint][fBackTrackSeg[ibackTrackSeg].fIndex[ipoint-6]]);
	  hitset[ipoint] = true;
	}
      }
      AliHLTMUONParticleSign sign = AliHLTMUONParticleSign(fHalfTrack[ibackTrackSeg].fCharge);
      track->fFlags = AliHLTMUONUtils::PackMansoTrackFlags(sign,hitset);
      
      track++;
      fNofTracks++;
      
    }//if nof connected is more than zero or not
  }//back track seg for loop
  
  size = fNofTracks;
  
  return true;
}

///__________________________________________________________________________


Bool_t AliHLTMUONFullTracker::SlatTrackSeg()
{

  ///Find the Slat Track Segments
  if(fNofPoints[6]==0 and fNofPoints[7]==0){
    HLTDebug("No Hits found in Stn4");
    return false;
  }else if(fNofPoints[8]==0 and fNofPoints[9]==0){
    HLTDebug("No Hits found in Stn5");
    return false;
  }

  Float_t trigX1,trigY1,trigZ1=AliMUONConstants::DefaultChamberZ(10);
  Float_t trigX2,trigY2,trigZ2=AliMUONConstants::DefaultChamberZ(12);
  Float_t extrapCh9X,extrapCh9Y,extrapCh9Z=AliMUONConstants::DefaultChamberZ(9);
  Float_t extrapCh8X,extrapCh8Y,extrapCh8Z=AliMUONConstants::DefaultChamberZ(8);
  Float_t extrapCh7X,extrapCh7Y,extrapCh7Z=AliMUONConstants::DefaultChamberZ(7);
  Float_t extrapCh6X,extrapCh6Y,extrapCh6Z=AliMUONConstants::DefaultChamberZ(6);

  Double_t distChFront,distChBack;
  Int_t nofFrontChPoints,nofBackChPoints;
  Int_t frontIndex[fgkMaxNofTracks], backIndex[fgkMaxNofTracks];
  
  AliHLTMUONRecHitStruct p2,p3,pSeg1,pSeg2,pSeg3;
  Double_t anglediff,anglediff1,anglediff2;
  Double_t minAngle = 2.0;
  
  Bool_t st5TrackletFound = false;
  Bool_t ch9PointFound = false;
  Bool_t ch8PointFound = false;
  Bool_t st4TrackletFound = false;
  Bool_t ch7PointFound = false;
  Bool_t ch6PointFound = false;

  Int_t index1,index2,index3,index4;
  IntPair cells[2][fgkMaxNofTracks]; ///cell array  for 5 stn for given trigger

  Float_t maxXDeflectionExtrap = 10.0 + 4.0;  ///simulation result 10.0
  Float_t extrapRatio = 0.2;            ///simulation result 0.2  
  Float_t circularWindow = 20.0 + 5.0 + 25.0;        ///simulatiuon result 20
  Float_t minAngleWindow = 2.0 + 1.0 + 2.0;        ///simulation result 2.0
  
  if(fFastTracking){
    maxXDeflectionExtrap = 10.0;  ///simulation result 10.0
    extrapRatio = 0.2;            ///simulation result 0.2  
    circularWindow = 20.0 ;        ///simulatiuon result 20
    minAngleWindow = 2.0;        ///simulation result 2.0
  }

  Float_t tx=0.0,ty=0.0;

  AliHLTUInt16_t detElemID,prevDetElemID=0xFFFF;
  AliHLTUInt8_t chamber;
  Int_t minTrgCh,maxTrgCh;

#ifdef PRINT_BACK
  printf("\nAliHLTMUONFullTracker::SlatTrackSeg()--Begin\n\n");
#endif
  

  for( Int_t itrig=0;itrig<fNofPoints[10];itrig++){
    
    st5TrackletFound = false;
    ch9PointFound = false;
    ch8PointFound = false;

    st4TrackletFound = false;
    ch7PointFound = false;
    ch6PointFound = false;

    minTrgCh = -1;
    maxTrgCh = -1;

    fNofCells[0] = 0;
    fNofCells[1] = 0;

    for( Int_t ihit=0;ihit<4;ihit++){

      AliHLTMUONUtils::UnpackRecHitFlags((fChPoint11[itrig]->fHit[ihit]).fFlags,chamber,detElemID);

      if(ihit==0 and detElemID!=0)
	minTrgCh = ihit;
      if(ihit==1 and detElemID!=0 and prevDetElemID==0)
	minTrgCh = ihit;
      if(ihit==2 and detElemID!=0)
	maxTrgCh = ihit;
      if(ihit==3 and detElemID!=0 and prevDetElemID==0)
	maxTrgCh = ihit;

      prevDetElemID = detElemID;
    }
    
    if(minTrgCh == -1 or maxTrgCh == -1){
      HLTDebug("Trigger hits not found in both trigger station minTrgCh : %d, maxTrgCh : %d, not harmful to HLT chain",minTrgCh,maxTrgCh);
      continue;
    }

    if( not fFastTracking){
      AliHLTMUONUtils::UnpackRecHitFlags((fChPoint11[itrig]->fHit[minTrgCh]).fFlags,chamber,detElemID);
      if(not fDetElemList[detElemID]){
	HLTDebug("Invalid detection element : %d, not harmful to HLT chain",detElemID);
	continue;
      }
      fChamberGeometryTransformer->Local2Global(detElemID,0.0,0.0,0.0,tx,ty,trigZ1);
      AliHLTMUONUtils::UnpackRecHitFlags((fChPoint11[itrig]->fHit[maxTrgCh]).fFlags,chamber,detElemID);
      if(not fDetElemList[detElemID]){
	HLTDebug("Invalid detection element : %d, not harmful to HLT chain",detElemID);
	continue;
      }
      fChamberGeometryTransformer->Local2Global(detElemID,0.0,0.0,0.0,tx,ty,trigZ2);
    }

    
    trigX1 = (fChPoint11[itrig]->fHit[minTrgCh]).fX;
    trigY1 = (fChPoint11[itrig]->fHit[minTrgCh]).fY;
    trigZ1 = (fChPoint11[itrig]->fHit[minTrgCh]).fZ;

    trigX2 = (fChPoint11[itrig]->fHit[maxTrgCh]).fX;
    trigY2 = (fChPoint11[itrig]->fHit[maxTrgCh]).fY;
    trigZ2 = (fChPoint11[itrig]->fHit[maxTrgCh]).fZ;
    
#ifdef PRINT_BACK
    printf("itrig : %d  trig 1 : (%f,%f,%f) \n",itrig,trigX1,trigY1,trigZ1);
    printf("itrig : %d  trig 2 : (%f,%f,%f) \n",itrig,trigX2,trigY2,trigZ2);
#endif

    /////////////////////////////////////////////////// Stn 5///////////////////////////////////////////////////////////////


    // #ifdef PRINT_BACK
    //     printf("\textrap9 : (%f,%f,%f)\n",extrapCh9X,extrapCh9Y,extrapCh9Z);
    //     printf("\textrap8 : (%f,%f,%f)\n",extrapCh8X,extrapCh8Y,extrapCh8Z);
    // #endif    
      
    nofFrontChPoints = 0; nofBackChPoints = 0;

    extrapCh9X = trigX1 * extrapCh9Z/trigZ1 ;
    extrapCh9Y = trigY1 + (trigY2-trigY1) * (extrapCh9Z-trigZ1)/(trigZ2 - trigZ1) ;
    for( Int_t ipointch9=0;ipointch9<fNofPoints[9];ipointch9++){

      if(not fFastTracking){
	AliHLTMUONUtils::UnpackRecHitFlags(fChPoint[9][ipointch9]->fFlags,chamber,detElemID);
	if(not fDetElemList[detElemID]){
	  HLTDebug("Invalid detection element : %d, not harmful to HLT chain",detElemID);
	  continue;
	}
	fChamberGeometryTransformer->Local2Global(detElemID,0.0,0.0,0.0,tx,ty,extrapCh9Z);
      
	extrapCh9X = trigX1 * extrapCh9Z/trigZ1 ;
	extrapCh9Y = trigY1 + (trigY2-trigY1) * (extrapCh9Z-trigZ1)/(trigZ2 - trigZ1) ;
      }
      
      if(nofBackChPoints < (fgkMaxNofTracks-1) && 
	 TMath::Abs(extrapCh9X-fChPoint[9][ipointch9]->fX) < maxXDeflectionExtrap && 
	 TMath::Abs(extrapCh9Y-fChPoint[9][ipointch9]->fY)/
	 ((fChPoint[9][ipointch9]->fX * fChPoint[9][ipointch9]->fX) + 
	  (fChPoint[9][ipointch9]->fY * fChPoint[9][ipointch9]->fY)) <= extrapRatio ){
	
	distChBack = sqrt((extrapCh9X-fChPoint[9][ipointch9]->fX)*(extrapCh9X-fChPoint[9][ipointch9]->fX) 
			  + (extrapCh9Y-fChPoint[9][ipointch9]->fY)*(extrapCh9Y-fChPoint[9][ipointch9]->fY));
	if(distChBack>circularWindow) continue;

#ifdef PRINT_BACK
	printf("\t\tpoints selected in Ch9  : (%f,%f,%f)\n",
	       distChBack,fChPoint[9][ipointch9]->fX,fChPoint[9][ipointch9]->fY,fChPoint[9][ipointch9]->fZ);
#endif

	backIndex[nofBackChPoints++] = ipointch9;
      }///if point found
    }/// ch10 loop


    extrapCh8X = trigX1 * extrapCh8Z/trigZ1 ;
    extrapCh8Y = trigY1 + (trigY2-trigY1) * (extrapCh8Z-trigZ1)/(trigZ2 - trigZ1) ;
    for( Int_t ipointch8=0;ipointch8<fNofPoints[8];ipointch8++){

      if(not fFastTracking){      
	AliHLTMUONUtils::UnpackRecHitFlags(fChPoint[8][ipointch8]->fFlags,chamber,detElemID);
	if(not fDetElemList[detElemID]){
	  HLTDebug("Invalid detection element : %d, not harmful to HLT chain",detElemID);
	  continue;
	}
	fChamberGeometryTransformer->Local2Global(detElemID,0.0,0.0,0.0,tx,ty,extrapCh8Z);
      
	extrapCh8X = trigX1 * extrapCh8Z/trigZ1 ;
	extrapCh8Y = trigY1 + (trigY2-trigY1) * (extrapCh8Z-trigZ1)/(trigZ2 - trigZ1) ;
      }

      if( nofFrontChPoints < (fgkMaxNofTracks-1) &&
	  TMath::Abs(extrapCh8X-fChPoint[8][ipointch8]->fX) < maxXDeflectionExtrap && 
	  TMath::Abs(extrapCh8Y-fChPoint[8][ipointch8]->fY)/
	  ((fChPoint[8][ipointch8]->fX * fChPoint[8][ipointch8]->fX) + 
	   (fChPoint[8][ipointch8]->fY * fChPoint[8][ipointch8]->fY)) <= extrapRatio ){
	
	distChFront = sqrt((extrapCh8X-fChPoint[8][ipointch8]->fX)*(extrapCh8X-fChPoint[8][ipointch8]->fX) 
			   + (extrapCh8Y-fChPoint[8][ipointch8]->fY)*(extrapCh8Y-fChPoint[8][ipointch8]->fY));
	
	if(distChFront>circularWindow) continue;
	
#ifdef PRINT_BACK
	printf("\t\tpoints selected in Ch8  : (%f,%f,%f)\n",
	       fChPoint[8][ipointch8]->fX,fChPoint[8][ipointch8]->fY,fChPoint[8][ipointch8]->fZ,distChFront);
#endif
	
	frontIndex[nofFrontChPoints++] = ipointch8;
      }///if point found
    }/// ch9 loop

    if(nofBackChPoints==0 and nofFrontChPoints==0) continue;
    
    minAngle = minAngleWindow;
    p3.fX = trigX1 ; p3.fY = trigY1 ; p3.fZ = trigZ1 ;
    for( Int_t ibackpoint=0;ibackpoint<nofBackChPoints;ibackpoint++){
      Sub(&p3,fChPoint[9][backIndex[ibackpoint]],&pSeg2);
      for( Int_t ifrontpoint=0;ifrontpoint<nofFrontChPoints;ifrontpoint++){
	Sub(fChPoint[9][backIndex[ibackpoint]],fChPoint[8][frontIndex[ifrontpoint]],&pSeg1);
	anglediff = TMath::RadToDeg()* Angle(&pSeg1,&pSeg2);
	// #ifdef PRINT_BACK
	// 	printf("\t\ttracklet-check-St5 : anglediff : %lf, minAngle : %lf\n",anglediff,minAngle);
	// #endif	
	if(anglediff<minAngle && fNofCells[1]<(fgkMaxNofTracks-1)){
	  st5TrackletFound = true;
	  cells[1][fNofCells[1]].fFirst =  frontIndex[ifrontpoint];
	  cells[1][fNofCells[1]].fSecond =  backIndex[ibackpoint];
	  fNofCells[1]++ ;
#ifdef PRINT_BACK
	  printf("\t\ttracklet-St5 : anglediff : %lf\n",anglediff);
	  printf("\t\t\tCh9  : (%f,%f,%f)\n",fChPoint[9][backIndex[ibackpoint]]->fX,
		 fChPoint[9][backIndex[ibackpoint]]->fY,fChPoint[9][backIndex[ibackpoint]]->fZ);
	  printf("\t\t\tCh8  : (%f,%f,%f)\n",fChPoint[8][frontIndex[ifrontpoint]]->fX,
		 fChPoint[8][frontIndex[ifrontpoint]]->fY,fChPoint[8][frontIndex[ifrontpoint]]->fZ);
#endif
	}///anglediff condition
      }///front
    }///back

    

    
    /// If tracklet not found, search for the single space point in Ch9 or in Ch8
    if(!st5TrackletFound){
      
      minAngle = minAngleWindow; 
      p3.fX = trigX2 ; p3.fY = trigY2 ; p3.fZ = trigZ2 ;
      p2.fX = trigX1 ; p2.fY = trigY1 ; p2.fZ = trigZ1 ;
      Sub(&p3,&p2,&pSeg2);
      
      ///Search in Ch9
      for( Int_t ibackpoint=0;ibackpoint<nofBackChPoints;ibackpoint++){
	Sub(&p2,fChPoint[9][backIndex[ibackpoint]],&pSeg1);
	anglediff = TMath::RadToDeg()*Angle(&pSeg1,&pSeg2);
	if(anglediff<minAngle && fNofCells[1]<(fgkMaxNofTracks-1)){
	  ch9PointFound = true;
	  cells[1][fNofCells[1]].fFirst =  -1;
	  cells[1][fNofCells[1]].fSecond =  backIndex[ibackpoint];
	  fNofCells[1]++ ;
#ifdef PRINT_BACK
	  printf("\t\tno st tracklet and single point-Ch9 : anglediff : %lf\n",anglediff);
#endif
	}
      }///back
      
      ///Search in Ch8
      for( Int_t ifrontpoint=0;ifrontpoint<nofFrontChPoints;ifrontpoint++){
	Sub(&p2,fChPoint[8][frontIndex[ifrontpoint]],&pSeg1);
	anglediff = TMath::RadToDeg()*Angle(&pSeg1,&pSeg2);
	if(anglediff<minAngle && fNofCells[1]<(fgkMaxNofTracks-1)){
	  ch8PointFound = true;
	  cells[1][fNofCells[1]].fFirst = frontIndex[ifrontpoint];
	  cells[1][fNofCells[1]].fSecond =  -1;
	  fNofCells[1]++ ;
#ifdef PRINT_BACK
	  printf("\t\tno st tracklet and single point-Ch8 : anglediff : %lf\n",anglediff);
#endif
	}
      }///front
    }///if no tracklets found condition
    
#ifdef PRINT_BACK
    printf("\tnofTracks found after stn 5 : %d\n",fNofCells[1]);
#endif
    
    if(!st5TrackletFound && !ch9PointFound && !ch8PointFound) continue;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    

    /////////////////////////////////////////////////// Stn 4///////////////////////////////////////////////////////////////
    
    // #ifdef PRINT_BACK
    //       printf("\textrap7 :  (%f,%f,%f)\n",extrapCh7X,extrapCh7Y,extrapCh7Z);
    //       printf("\textrap6 :  (%f,%f,%f)\n",extrapCh6X,extrapCh6Y,extrapCh6Z);
    // #endif    
    
    nofFrontChPoints = 0; nofBackChPoints = 0;

    extrapCh7X = trigX1 * extrapCh7Z/trigZ1 ;
    extrapCh7Y = trigY1 + (trigY2-trigY1) * (extrapCh7Z-trigZ1)/(trigZ2 - trigZ1) ;
    for( Int_t ipointch7=0;ipointch7<fNofPoints[7];ipointch7++){
     
      if(not fFastTracking){
	AliHLTMUONUtils::UnpackRecHitFlags(fChPoint[7][ipointch7]->fFlags,chamber,detElemID);
	if(not fDetElemList[detElemID]){
	  HLTDebug("Invalid detection element : %d, not harmful to HLT chain",detElemID);
	  continue;
	}
	fChamberGeometryTransformer->Local2Global(detElemID,0.0,0.0,0.0,tx,ty,extrapCh7Z);

	extrapCh7X = trigX1 * extrapCh7Z/trigZ1 ;
	extrapCh7Y = trigY1 + (trigY2-trigY1) * (extrapCh7Z-trigZ1)/(trigZ2 - trigZ1) ;
      }

      if( nofBackChPoints < (fgkMaxNofTracks-1) &&
	  TMath::Abs(extrapCh7X-fChPoint[7][ipointch7]->fX) < maxXDeflectionExtrap && 
	  TMath::Abs(extrapCh7Y-fChPoint[7][ipointch7]->fY)/
	  ((fChPoint[7][ipointch7]->fX * fChPoint[7][ipointch7]->fX) + 
	   (fChPoint[7][ipointch7]->fY * fChPoint[7][ipointch7]->fY)) <= extrapRatio ){
	
	distChBack = sqrt((extrapCh7X-fChPoint[7][ipointch7]->fX)*(extrapCh7X-fChPoint[7][ipointch7]->fX) 
			  + (extrapCh7Y-fChPoint[7][ipointch7]->fY)*(extrapCh7Y-fChPoint[7][ipointch7]->fY));
	
	if(distChBack>circularWindow) continue;
#ifdef PRINT_BACK
	printf("\t\tpoints selected in Ch7  : (%f,%f,%f)\n",
	       fChPoint[7][ipointch7]->fX,fChPoint[7][ipointch7]->fY,fChPoint[7][ipointch7]->fZ);
#endif
	
	backIndex[nofBackChPoints++] = ipointch7;
      }///if point found
    }///ch8 loop

    extrapCh6X = trigX1 * extrapCh6Z/trigZ1 ;
    extrapCh6Y = trigY1 + (trigY2-trigY1) * (extrapCh6Z-trigZ1)/(trigZ2 - trigZ1) ;
    for( Int_t ipointch6=0;ipointch6<fNofPoints[6];ipointch6++){
      
      if(not fFastTracking){
	AliHLTMUONUtils::UnpackRecHitFlags(fChPoint[6][ipointch6]->fFlags,chamber,detElemID);
	if(not fDetElemList[detElemID]){
	  HLTDebug("Invalid detection element : %d, not harmful to HLT chain",detElemID);
	  continue;
	}
	fChamberGeometryTransformer->Local2Global(detElemID,0.0,0.0,0.0,tx,ty,extrapCh6Z);	
      
      
	extrapCh6X = trigX1 * extrapCh6Z/trigZ1 ;
	extrapCh6Y = trigY1 + (trigY2-trigY1) * (extrapCh6Z-trigZ1)/(trigZ2 - trigZ1) ;
      }
      
      if(nofFrontChPoints < (fgkMaxNofTracks-1) && 
	 TMath::Abs(extrapCh6X-fChPoint[6][ipointch6]->fX) < maxXDeflectionExtrap && 
	 TMath::Abs(extrapCh6Y-fChPoint[6][ipointch6]->fY)/
	 ((fChPoint[6][ipointch6]->fX * fChPoint[6][ipointch6]->fX) + 
	  (fChPoint[6][ipointch6]->fY * fChPoint[6][ipointch6]->fY)) <= extrapRatio ){
	
	distChFront = sqrt((extrapCh6X-fChPoint[6][ipointch6]->fX)*(extrapCh6X-fChPoint[6][ipointch6]->fX) 
			   + (extrapCh6Y-fChPoint[6][ipointch6]->fY)*(extrapCh6Y-fChPoint[6][ipointch6]->fY));
	if(distChFront>circularWindow) continue;
	
#ifdef PRINT_BACK
	printf("\t\tpoints selected in Ch6  : (%f,%f,%f)\n",
	       fChPoint[6][ipointch6]->fX,fChPoint[6][ipointch6]->fY,fChPoint[6][ipointch6]->fZ);
#endif
	frontIndex[nofFrontChPoints++] = ipointch6;
      }///if point found
    }/// ch7 loop

    if(nofBackChPoints==0 and nofFrontChPoints==0) continue;
    
    minAngle = minAngleWindow;
    p3.fX = trigX1 ; p3.fY = trigY1 ; p3.fZ = trigZ1 ;
    for( Int_t ibackpoint=0;ibackpoint<nofBackChPoints;ibackpoint++){
      Sub(&p3,fChPoint[7][backIndex[ibackpoint]],&pSeg2);
      for( Int_t ifrontpoint=0;ifrontpoint<nofFrontChPoints;ifrontpoint++){
	Sub(fChPoint[7][backIndex[ibackpoint]],fChPoint[6][frontIndex[ifrontpoint]],&pSeg1);
	anglediff = TMath::RadToDeg() * Angle(&pSeg1,&pSeg2);
	if(anglediff<minAngle && fNofCells[0]<(fgkMaxNofTracks-1)){
	  st4TrackletFound = true;
	  cells[0][fNofCells[0]].fFirst =  frontIndex[ifrontpoint];
	  cells[0][fNofCells[0]].fSecond =  backIndex[ibackpoint];
	  fNofCells[0]++ ;
#ifdef PRINT_BACK
	  printf("\t\ttracklet-St4 : anglediff : %lf\n",anglediff);
	  printf("\t\t\tCh7  : (%f,%f,%f)\n",fChPoint[7][backIndex[ibackpoint]]->fX,
		 fChPoint[7][backIndex[ibackpoint]]->fY,fChPoint[7][backIndex[ibackpoint]]->fZ);
	  printf("\t\t\tCh6  : (%f,%f,%f)\n",fChPoint[6][frontIndex[ifrontpoint]]->fX,
		 fChPoint[6][frontIndex[ifrontpoint]]->fY,fChPoint[6][frontIndex[ifrontpoint]]->fZ);
#endif
	}///anglediff condn
      }///front
    }///back

    

    
    /// If tracklet not found search for the single space point in Ch7 or in Ch6
    if(!st4TrackletFound){

      minAngle = minAngleWindow; 
      p3.fX = trigX2 ; p3.fY = trigY2 ; p3.fZ = trigZ2 ;
      p2.fX = trigX1 ; p2.fY = trigY1 ; p2.fZ = trigZ1 ;
      Sub(&p3,&p2,&pSeg2);

      ///Search in Ch7
      for( Int_t ibackpoint=0;ibackpoint<nofBackChPoints;ibackpoint++){
	Sub(&p2,fChPoint[7][backIndex[ibackpoint]],&pSeg1);
	anglediff = TMath::RadToDeg()*Angle(&pSeg1,&pSeg2);
	if(anglediff<minAngle && fNofCells[0]<(fgkMaxNofTracks-1)){
	  ch7PointFound = true;
	  cells[0][fNofCells[0]].fFirst =  -1;
	  cells[0][fNofCells[0]].fSecond =  backIndex[ibackpoint];
	  fNofCells[0]++ ;
#ifdef PRINT_BACK
	  printf("\t\tno st tracklet and single point-Ch7 : anglediff : %lf\n",anglediff);
#endif
	}
      }///back
      
      ///Search in Ch6
      for( Int_t ifrontpoint=0;ifrontpoint<nofFrontChPoints;ifrontpoint++){
	Sub(&p2,fChPoint[6][frontIndex[ifrontpoint]],&pSeg1);
	anglediff = TMath::RadToDeg()*Angle(&pSeg1,&pSeg2);
	if(anglediff<minAngle && fNofCells[0]<(fgkMaxNofTracks-1)){
	  ch6PointFound = true;
	  cells[0][fNofCells[0]].fFirst = frontIndex[ifrontpoint];
	  cells[0][fNofCells[0]].fSecond =  -1;
	  fNofCells[0]++ ;
#ifdef PRINT_BACK
	  printf("\t\tno st tracklet and single point-Ch6 : anglediff : %lf\n",anglediff);
#endif
	}
      }///front
    }///if no tracklets found condition
    
#ifdef PRINT_BACK
    printf("\tnofTracks found after stn 4 : %d\n",fNofCells[0]);
#endif
    
    if(!st4TrackletFound && !ch7PointFound && !ch6PointFound) continue;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ;
      
    ////////////////////////////////////////////// Analyse and fill trackseg array////////////////////////////////////////
    ;
#ifdef PRINT_BACK
    printf("\tfNofbackTrackSeg : %d, st5TrackletFound : %d, st4TrackletFound : %d\n",fNofbackTrackSeg,st5TrackletFound,st4TrackletFound);
#endif

    if(st5TrackletFound && st4TrackletFound){
      
      minAngle = minAngleWindow;
      
      for( Int_t itrackletfront=0;itrackletfront<fNofCells[0];itrackletfront++){
	index1 = cells[0][itrackletfront].fFirst ;
	index2 = cells[0][itrackletfront].fSecond ;
	Sub(fChPoint[7][index2],fChPoint[6][index1],&pSeg1);
	for( Int_t itrackletback=0;itrackletback<fNofCells[1];itrackletback++){
	  index3 = cells[1][itrackletback].fFirst ;
	  index4 = cells[1][itrackletback].fSecond ;
	  Sub(fChPoint[8][index3],fChPoint[7][index2],&pSeg2);
	  Sub(fChPoint[9][index4],fChPoint[8][index3],&pSeg3);
	  anglediff = Angle(&pSeg1,&pSeg2) + Angle(&pSeg2,&pSeg3);
	  if(anglediff<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
	    fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = index1;
	    fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = index2;
	    fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = index3;
	    fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = index4;
	    fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
	    minAngle = anglediff;
#ifdef PRINT_BACK
	    printf("\t\ttracklet-St4 and St5 : anglediff : %lf\n",anglediff);
	    printf("\t\t\tCh9  : (%f,%f,%f)\n",fChPoint[9][index4]->fX,
		   fChPoint[9][index4]->fY,fChPoint[9][index4]->fZ);
	    printf("\t\t\tCh8  : (%f,%f,%f)\n",fChPoint[8][index3]->fX,
		   fChPoint[8][index3]->fY,fChPoint[8][index3]->fZ);
	    printf("\t\t\tCh7  : (%f,%f,%f)\n",fChPoint[7][index2]->fX,
		   fChPoint[7][index2]->fY,fChPoint[7][index2]->fZ);
	    printf("\t\t\tCh6  : (%f,%f,%f)\n",fChPoint[6][index1]->fX,
		   fChPoint[6][index1]->fY,fChPoint[6][index1]->fZ);
#endif
	  }///if minangle
	}///for of front ch
      }///for loop of back ch
      
      if(minAngle<minAngleWindow)
	fNofbackTrackSeg++;
      
    }else if(st5TrackletFound && (ch7PointFound || ch6PointFound)){
      
      
      nofFrontChPoints = 0; nofBackChPoints = 0;
      for( Int_t ifrontpoint=0;ifrontpoint<fNofCells[0];ifrontpoint++){
	if(cells[0][ifrontpoint].fFirst==-1 && nofBackChPoints<(fgkMaxNofTracks-1))
	  backIndex[nofBackChPoints++] = cells[0][ifrontpoint].fSecond;
	else if(cells[0][ifrontpoint].fSecond==-1 && nofFrontChPoints<(fgkMaxNofTracks-1))
	  frontIndex[nofFrontChPoints++] = cells[0][ifrontpoint].fFirst; 
      }
      
      minAngle = minAngleWindow;
      if(nofFrontChPoints>0 && nofBackChPoints>0){
	for( Int_t itrackletback=0;itrackletback<fNofCells[1];itrackletback++){
	  index3 = cells[1][itrackletback].fFirst ;
	  index4 = cells[1][itrackletback].fSecond ;
	  Sub(fChPoint[9][index4],fChPoint[8][index3],&pSeg3);
	  for( Int_t ibackpoint=0;ibackpoint<nofBackChPoints;ibackpoint++){
	    Sub(fChPoint[8][index3],fChPoint[7][backIndex[ibackpoint]],&pSeg2);
	    for( Int_t ifrontpoint=0;ifrontpoint<nofFrontChPoints;ifrontpoint++){
	      Sub(fChPoint[7][backIndex[ibackpoint]],fChPoint[6][frontIndex[ifrontpoint]],&pSeg1);
	      anglediff = TMath::RadToDeg()*(Angle(&pSeg1,&pSeg2) + Angle(&pSeg2,&pSeg3))/2.0;

	      if(anglediff<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
		fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = frontIndex[ifrontpoint];
		fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = backIndex[ibackpoint] ;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = index3;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = index4;
		fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
		minAngle = anglediff;
		continue;
	      }

	      Sub(fChPoint[8][index3],fChPoint[6][frontIndex[ifrontpoint]],&pSeg1);
	      anglediff1 = TMath::RadToDeg() * Angle(&pSeg1,&pSeg3);
	      anglediff2 = TMath::RadToDeg() * Angle(&pSeg2,&pSeg3);

	      if( anglediff1 < anglediff2 && anglediff1<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
		fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = frontIndex[ifrontpoint];
		fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = -1;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = index3;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = index4;
		fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
		minAngle = anglediff1;
		continue;
	      }

	      if( anglediff2 < anglediff1 && anglediff2<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
		fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = -1;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = backIndex[ibackpoint] ;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = index3;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = index4;
		fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
		minAngle = anglediff2;
	      }
	      
	    }///loop of ifrontpoint
	  }///loop of ibackpoint
	}/// for loop of St5 cells
      }else if(nofFrontChPoints>0){
	
	for( Int_t itrackletback=0;itrackletback<fNofCells[1];itrackletback++){
	  index3 = cells[1][itrackletback].fFirst ;
	  index4 = cells[1][itrackletback].fSecond ;
	  Sub(fChPoint[9][index4],fChPoint[8][index3],&pSeg3);

	  for( Int_t ifrontpoint=0;ifrontpoint<nofFrontChPoints;ifrontpoint++){
	    Sub(fChPoint[8][index3],fChPoint[6][frontIndex[ifrontpoint]],&pSeg2);
	    
	    anglediff = TMath::RadToDeg() * Angle(&pSeg2,&pSeg3);
	    if( anglediff<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = frontIndex[ifrontpoint];
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = -1;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = index3;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = index4;
	      fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
	      minAngle = anglediff;
	    }///if anglediff
	  }///backch loop
	}///tracklet loop

      }else{ /// if(nofBackChPoints>0){
	for( Int_t itrackletback=0;itrackletback<fNofCells[1];itrackletback++){
	  index3 = cells[1][itrackletback].fFirst ;
	  index4 = cells[1][itrackletback].fSecond ;
	  Sub(fChPoint[9][index4],fChPoint[8][index3],&pSeg3);
	  
	  for( Int_t ibackpoint=0;ibackpoint<nofBackChPoints;ibackpoint++){
	    Sub(fChPoint[8][index3],fChPoint[7][backIndex[ibackpoint]],&pSeg2);
	    
	    anglediff = TMath::RadToDeg() * Angle(&pSeg2,&pSeg3);
	    
	    if( anglediff<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = -1;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = backIndex[ibackpoint] ;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = index3;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = index4;
	      fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
	      minAngle = anglediff;
	    }///if anglediff
	  }///backch loop
	}///tracklet loop

      }///condn for if(nofFrontChPoints>0)

      if(minAngle<minAngleWindow)
	fNofbackTrackSeg++;

    }else if((ch9PointFound || ch8PointFound) && st4TrackletFound){

      nofFrontChPoints = 0; nofBackChPoints = 0;
      for( Int_t ibackpoint=0;ibackpoint<fNofCells[1];ibackpoint++){
	if(cells[1][ibackpoint].fFirst==-1 && nofBackChPoints<(fgkMaxNofTracks-1))
	  backIndex[nofBackChPoints++] = cells[1][ibackpoint].fSecond;
	else if(nofFrontChPoints<(fgkMaxNofTracks-1))
	  frontIndex[nofFrontChPoints++] = cells[1][ibackpoint].fFirst; 
      }
      
      minAngle = minAngleWindow;
      if(nofFrontChPoints>0 && nofBackChPoints>0){

	for( Int_t itrackletfront=0;itrackletfront<fNofCells[0];itrackletfront++){
	  index1 = cells[0][itrackletfront].fFirst ;
	  index2 = cells[0][itrackletfront].fSecond ;
	  
	  Sub(fChPoint[7][index2],fChPoint[6][index1],&pSeg1);
	  
	  for( Int_t ifrontpoint=0;ifrontpoint<nofFrontChPoints;ifrontpoint++){
	    
	    Sub(fChPoint[8][frontIndex[ifrontpoint]],fChPoint[7][index2],&pSeg2);
	    
	    for( Int_t ibackpoint=0;ibackpoint<nofBackChPoints;ibackpoint++){
	      
	      Sub(fChPoint[9][backIndex[ibackpoint]],fChPoint[8][frontIndex[ifrontpoint]],&pSeg3);
	      
	      anglediff = TMath::RadToDeg()*(Angle(&pSeg1,&pSeg2) + Angle(&pSeg2,&pSeg3))/2.0;
	      if(anglediff<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
		fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = index1;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = index2;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = frontIndex[ifrontpoint];
		fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = backIndex[ibackpoint] ;
		fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
		minAngle = anglediff;
		continue;
	      }
	      
	      Sub(fChPoint[9][backIndex[ibackpoint]],fChPoint[7][index2],&pSeg3);
	      
	      anglediff1 = TMath::RadToDeg() * Angle(&pSeg1,&pSeg2);
	      anglediff2 = TMath::RadToDeg() * Angle(&pSeg1,&pSeg3);
	      
	      if( anglediff1 < anglediff2 && anglediff1<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
		fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = index1;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = index2;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = frontIndex[ifrontpoint];
		fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = -1;
		fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
		minAngle = anglediff1;
		continue;
	      }
	      
	      if( anglediff2 < anglediff1 && anglediff2<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
		fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = index1;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = index2;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = -1;
		fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = backIndex[ibackpoint] ;
		fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
		minAngle = anglediff2;
	      }
	      
	    }///loop of ifrontpoint
	  }///loop of ibackpoint
	}/// for loop of St5 cells
      }else if(nofFrontChPoints>0){

	for( Int_t itrackletfront=0;itrackletfront<fNofCells[0];itrackletfront++){
	  index1 = cells[0][itrackletfront].fFirst ;
	  index2 = cells[0][itrackletfront].fSecond ;
	  
	  Sub(fChPoint[7][index2],fChPoint[6][index1],&pSeg1);
	  
	  for( Int_t ifrontpoint=0;ifrontpoint<nofFrontChPoints;ifrontpoint++){
	    
	    Sub(fChPoint[8][frontIndex[ifrontpoint]],fChPoint[7][index2],&pSeg2);
	    
	    anglediff = TMath::RadToDeg() * Angle(&pSeg1,&pSeg2);
	    if( anglediff<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = index1;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = index2;;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = frontIndex[ifrontpoint];
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = -1;
	      fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
	      minAngle = anglediff;
	    }///if anglediff
	  }///point loop
	}///tracklet loop

      }else{ /// if(nofBackChPoints>0){

	for( Int_t itrackletfront=0;itrackletfront<fNofCells[0];itrackletfront++){
	  index1 = cells[0][itrackletfront].fFirst ;
	  index2 = cells[0][itrackletfront].fSecond ;

	  Sub(fChPoint[7][index2],fChPoint[6][index1],&pSeg1);

	  for( Int_t ibackpoint=0;ibackpoint<nofBackChPoints;ibackpoint++){
	    Sub(fChPoint[9][backIndex[ibackpoint]],fChPoint[6][index1],&pSeg3);
	    anglediff = TMath::RadToDeg()* Angle(&pSeg1,&pSeg3);
	    if(anglediff<minAngle && fNofbackTrackSeg<(fgkMaxNofTracks-1)){
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[0] = index1;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[1] = index2;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[2] = -1;
	      fBackTrackSeg[fNofbackTrackSeg].fIndex[3] = backIndex[ibackpoint] ;
	      fBackTrackSeg[fNofbackTrackSeg].fTrigRec = fChPoint11[itrig]->fId;
	      minAngle = anglediff;
	    }
	  }///backch loop
	}///tracklet loop

      }///condn for if(nofFrontChPoints>0)

      if(minAngle<minAngleWindow)
	fNofbackTrackSeg++;
      
    }else if((ch9PointFound || ch8PointFound) && (ch7PointFound || ch6PointFound)){
      ///To Do :  To be analysed for two points out of four slat chambers
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef PRINT_BACK
    printf("\n");
#endif

  }///trigger loop

  Float_t meanX1,meanX2,meanY1,meanY2,meanZ1,meanZ2;

  for( Int_t ibacktrackseg=0;ibacktrackseg<fNofbackTrackSeg;ibacktrackseg++){
    
#ifdef PRINT_BACK
    printf("Index : (%d,%d,%d,%d) \n",fBackTrackSeg[ibacktrackseg].fIndex[0],
	   fBackTrackSeg[ibacktrackseg].fIndex[1],fBackTrackSeg[ibacktrackseg].fIndex[2],
	   fBackTrackSeg[ibacktrackseg].fIndex[3]);
#endif

    if(fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1 && fBackTrackSeg[ibacktrackseg].fIndex[1]!=-1 ){
      meanX1 = (fChPoint[6][fBackTrackSeg[ibacktrackseg].fIndex[0]]->fX 
		+ fChPoint[7][fBackTrackSeg[ibacktrackseg].fIndex[1]]->fX)/2.0 ;
      meanY1 = (fChPoint[6][fBackTrackSeg[ibacktrackseg].fIndex[0]]->fY 
		+ fChPoint[7][fBackTrackSeg[ibacktrackseg].fIndex[1]]->fY)/2.0 ;
      meanZ1 = (fChPoint[6][fBackTrackSeg[ibacktrackseg].fIndex[0]]->fZ 
		+ fChPoint[7][fBackTrackSeg[ibacktrackseg].fIndex[1]]->fZ)/2.0 ;
    }else if(fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1 && fBackTrackSeg[ibacktrackseg].fIndex[1]==-1 ){
      meanX1 = fChPoint[6][fBackTrackSeg[ibacktrackseg].fIndex[0]]->fX ;
      meanY1 = fChPoint[6][fBackTrackSeg[ibacktrackseg].fIndex[0]]->fY ;
      meanZ1 = fChPoint[6][fBackTrackSeg[ibacktrackseg].fIndex[0]]->fZ ;
    }else{
      meanX1 = fChPoint[7][fBackTrackSeg[ibacktrackseg].fIndex[1]]->fX ;
      meanY1 = fChPoint[7][fBackTrackSeg[ibacktrackseg].fIndex[1]]->fY ;
      meanZ1 = fChPoint[7][fBackTrackSeg[ibacktrackseg].fIndex[1]]->fZ ;
    }
    
    if(fBackTrackSeg[ibacktrackseg].fIndex[2]!=-1 && fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1 ){
      meanX2 = (fChPoint[8][fBackTrackSeg[ibacktrackseg].fIndex[2]]->fX 
		+ fChPoint[9][fBackTrackSeg[ibacktrackseg].fIndex[3]]->fX)/2.0 ;
      meanY2 = (fChPoint[8][fBackTrackSeg[ibacktrackseg].fIndex[2]]->fY 
		+ fChPoint[9][fBackTrackSeg[ibacktrackseg].fIndex[3]]->fY)/2.0 ;
      meanZ2 = (fChPoint[8][fBackTrackSeg[ibacktrackseg].fIndex[2]]->fZ 
		+ fChPoint[9][fBackTrackSeg[ibacktrackseg].fIndex[3]]->fZ)/2.0 ;
    }else if(fBackTrackSeg[ibacktrackseg].fIndex[2]!=-1 && fBackTrackSeg[ibacktrackseg].fIndex[3]==-1 ){
      meanX2 = fChPoint[8][fBackTrackSeg[ibacktrackseg].fIndex[2]]->fX ;
      meanY2 = fChPoint[8][fBackTrackSeg[ibacktrackseg].fIndex[2]]->fY ;
      meanZ2 = fChPoint[8][fBackTrackSeg[ibacktrackseg].fIndex[2]]->fZ ;
    }else{
      meanX2 = fChPoint[9][fBackTrackSeg[ibacktrackseg].fIndex[3]]->fX ;
      meanY2 = fChPoint[9][fBackTrackSeg[ibacktrackseg].fIndex[3]]->fY ;
      meanZ2 = fChPoint[9][fBackTrackSeg[ibacktrackseg].fIndex[3]]->fZ ;
    }
    fExtrapSt3X[ibacktrackseg] = meanX1 + (fgkTrackDetCoordinate[2]-meanZ1)*(meanX2-meanX1)/(meanZ2-meanZ1);
    fExtrapSt3Y[ibacktrackseg] = meanY1 + (fgkTrackDetCoordinate[2]-meanZ1)*(meanY2-meanY1)/(meanZ2-meanZ1);  
    fInclinationBack[ibacktrackseg] = (meanX2-meanX1)/(meanZ2-meanZ1) ;
    fNofConnectedfrontTrackSeg[ibacktrackseg] = 0;
  }///backtrigseg loop

#ifdef PRINT_BACK
  printf("AliHLTMUONFullTracker::SlatTrackSeg()--End\n");
  printf("\n\n");
#endif
  
  return true;
}
///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::PrelimMomCalc()
{
  /// momentum calculation for half tracks

  Cluster clus1,clus2;
  Float_t xf,yf,thetaDev,zf = 0.5*(AliMUONConstants::DefaultChamberZ(4) + AliMUONConstants::DefaultChamberZ(5));
  Float_t p,px,py,pz;

  for( Int_t ibacktrackseg=0;ibacktrackseg<fNofbackTrackSeg;ibacktrackseg++){
    
    
    Int_t maxIndex = (fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1)?3:2;
    Int_t maxCh = (fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1)?9:8;

    Int_t minIndex = (fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1)?0:1;
    Int_t minCh = (fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1)?6:7;
    
    clus2.fX = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fX ;
    clus2.fY = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fY ;
    clus2.fZ = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fZ ;
    
    clus1.fX = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fX ;
    clus1.fY = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fY ;
    clus1.fZ = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fZ ;

    thetaDev= (1/zf)*(clus1.fY*clus2.fZ - clus2.fY*clus1.fZ)/(clus2.fZ - clus1.fZ);
    xf = clus1.fX*zf/clus1.fZ;
    yf = clus2.fY - ((clus2.fY - clus1.fY)*(clus2.fZ-zf))/(clus2.fZ - clus1.fZ);
    p = 3.0*0.3/thetaDev;
    px = p*xf/zf;
    py = p*yf/zf;;
    pz = sqrt((p*p-(px*px + py*py)));
    if(zf<0)pz = -pz;
    
    fHalfTrack[ibacktrackseg].fCharge = (p>0)?-1:+1;
    fHalfTrack[ibacktrackseg].fPx = p*xf/zf;
    fHalfTrack[ibacktrackseg].fPy = p*yf/zf;
    fHalfTrack[ibacktrackseg].fPz = pz;
    
  }
  
  return true;
}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::QuadTrackSeg()
{
  ///Find the Track Segments in the Quadrant chambers
  if(fNofPoints[0]==0 and fNofPoints[1]==0){
    HLTDebug("No Hits found in Stn4");
    return false;
  }else if(fNofPoints[2]==0 and fNofPoints[3]==0){
    HLTDebug("No Hits found in Stn5");
    return false;
  }else if(fNofbackTrackSeg==0){
    HLTDebug("No Hits found in Stn5 and Stn4 so no tracking is done in quadrants");
    return false;
  }

  
  Float_t meanX1,meanX2,meanY1,meanY2,meanZ1,meanZ2;
  Float_t expectSt3X,expectSt3Y,inclinationFront;

  AliHLTMUONRecHitStruct pSeg1,pSeg2,pSeg3;
  Double_t anglediff;///,anglediff1,anglediff2;
  Double_t minAngle = -1.0;
  Int_t  indexMinAngleFront = -1;
  Int_t  indexMinAngleBack = -1;
  Int_t backIndex = -1;

  Int_t ch3CellPoint[fgkMaxNofCellsPerCh],ch2CellPoint[fgkMaxNofCellsPerCh],nofSt2Cells=0;
  Int_t ch1CellPoint[fgkMaxNofCellsPerCh],ch0CellPoint[fgkMaxNofCellsPerCh],nofSt1Cells=0;
  Bool_t isConnected[fgkMaxNofTracks];

  Float_t distDiffX = 4.0;                    ///simulation result 4.0
  Float_t distDiffY = 10.0 ;                   ///simulation result 4.0
  ///float closeToY0 = 10.0;                    ///simulation result 10.0 
  Float_t minAngleWindow = 2.0 ;               ///simulation result 2.0
  Float_t diffDistStX = 25.0;                  ///simulation result 25.0
  Float_t diffDistStY = 75.0;                  ///simulation result 25.0
  Float_t st3WindowX = 40.0   ;                ///simulation result 40.0
  Float_t st3WindowY = 10.0;                   ///simulation result 10.0

  if(fFastTracking){
    distDiffX = 4.0;                    ///simulation result 4.0
    distDiffY = 4.0 ;                   ///simulation result 4.0
    ///float closeToY0 = 10.0;                    ///simulation result 10.0 
    minAngleWindow = 2.0 ;               ///simulation result 2.0
    diffDistStX = 25.0;                  ///simulation result 25.0
    diffDistStY = 25.0;                  ///simulation result 25.0
    st3WindowX = 40.0   ;                ///simulation result 40.0
    st3WindowY = 10.0;                   ///simulation result 10.0
  }
  
  ///   Float_t inclinationWindow = 0.04;            ///inclination window   
  ///   Float_t st3WindowXOp2 = 40.0 ;                 ///simulation result 40.0
  ///   Float_t st3WindowYOp2 = 10.0;                ///simulation result 10.0
  
  
#ifdef PRINT_FRONT
  printf("\nAliHLTMUONFullTracker::QuadTrackSeg()--Begin\n\n");
#endif

  for( Int_t ibackpoint=0;ibackpoint<fNofPoints[3];ibackpoint++){
    for( Int_t ifrontpoint=0;ifrontpoint<fNofPoints[2];ifrontpoint++){

      if(TMath::Abs(fChPoint[3][ibackpoint]->fX - fChPoint[2][ifrontpoint]->fX)<distDiffX 
	 && TMath::Abs(fChPoint[3][ibackpoint]->fY - fChPoint[2][ifrontpoint]->fY)<distDiffY){
	 

	/// if((TMath::Abs(fChPoint[3][ibackpoint]->fY) > closeToY0 && 
	///     TMath::Abs(fChPoint[3][ibackpoint]->fY) < TMath::Abs(fChPoint[2][ifrontpoint]->fY)) ||
	///    nofSt2Cells >= (fgkMaxNofCellsPerCh-1)) continue;

	if(nofSt2Cells >= (fgkMaxNofCellsPerCh-1)) continue;
	
	// #ifdef PRINT_FRONT
	// 	printf("\t\t\tCh3  : %d, (%f,%f,%f)\n",
	// 	       nofSt2Cells,fChPoint[3][ibackpoint]->fX,fChPoint[3][ibackpoint]->fY,fChPoint[3][ibackpoint]->fZ);
	// 	printf("\t\t\tCh2  :(%f,%f,%f)\n\n",
	// 	       fChPoint[2][ifrontpoint]->fX,fChPoint[2][ifrontpoint]->fY,fChPoint[2][ifrontpoint]->fZ);
	// #endif
	
	ch3CellPoint[nofSt2Cells] = ibackpoint;
	ch2CellPoint[nofSt2Cells] = ifrontpoint;
	nofSt2Cells++; 
	
      }///if point found
    }///frontch
  }///backch

  if(nofSt2Cells==0){
    HLTDebug("No Hits found in Stn2");
    return false;
  }
  
  for( Int_t ibackpoint=0;ibackpoint<fNofPoints[1];ibackpoint++){
    for( Int_t ifrontpoint=0;ifrontpoint<fNofPoints[0];ifrontpoint++){
      
      if(TMath::Abs(fChPoint[1][ibackpoint]->fX - fChPoint[0][ifrontpoint]->fX)< distDiffX
	 && TMath::Abs(fChPoint[1][ibackpoint]->fY - fChPoint[0][ifrontpoint]->fY)<distDiffY){
	 
	
	/// if((TMath::Abs(fChPoint[1][ibackpoint]->fY) > closeToY0 && 
	///     TMath::Abs(fChPoint[1][ibackpoint]->fY) < TMath::Abs(fChPoint[0][ifrontpoint]->fY)) ||
	///    nofSt1Cells >= (fgkMaxNofCellsPerCh-1)) continue;

	if(nofSt1Cells >= (fgkMaxNofCellsPerCh-1)) continue;
	   
	
	// #ifdef PRINT_FRONT
	// 	printf("\t\t\tCh1  : %d, (%f,%f,%f)\n",
	// 	       nofSt1Cells,fChPoint[1][ibackpoint]->fX,fChPoint[1][ibackpoint]->fY,fChPoint[1][ibackpoint]->fZ);
	// 	printf("\t\t\tCh0  :(%f,%f,%f)\n\n",
	// 	       fChPoint[0][ifrontpoint]->fX,fChPoint[0][ifrontpoint]->fY,fChPoint[0][ifrontpoint]->fZ);
	// #endif
	ch1CellPoint[nofSt1Cells] = ibackpoint;
	ch0CellPoint[nofSt1Cells] = ifrontpoint;
	nofSt1Cells++;
      }///if point found
    }///frontch
  }///backch
  if(nofSt1Cells==0){
    HLTDebug("No Hits found in Stn1");
    return false;
  }
  
#ifdef PRINT_FRONT
  printf("\tnofSt1Cells : %d, nofSt2Cells : %d\n",nofSt1Cells,nofSt2Cells);
#endif
  
  for( Int_t ibacktrackseg=0;ibacktrackseg<fNofbackTrackSeg;ibacktrackseg++)
    isConnected[ibacktrackseg] = false;
  

  ///First Check for Tracklets in two front stations
  for( Int_t itrackletfront=0;itrackletfront<nofSt1Cells;itrackletfront++){
    
    Sub(fChPoint[1][ch1CellPoint[itrackletfront]],fChPoint[0][ch0CellPoint[itrackletfront]],&pSeg1);
    
    minAngle = minAngleWindow;
    indexMinAngleBack = -1;
    indexMinAngleFront = -1;

    meanX1 = (fChPoint[0][ch0CellPoint[itrackletfront]]->fX 
	      + fChPoint[1][ch1CellPoint[itrackletfront]]->fX)/2.0 ;
    meanY1 = (fChPoint[0][ch0CellPoint[itrackletfront]]->fY 
	      + fChPoint[1][ch1CellPoint[itrackletfront]]->fY)/2.0 ;
    meanZ1 = (fChPoint[0][ch0CellPoint[itrackletfront]]->fZ 
	      + fChPoint[1][ch1CellPoint[itrackletfront]]->fZ)/2.0 ;
    
    for( Int_t itrackletback=0;itrackletback<nofSt2Cells;itrackletback++){
      // #ifdef PRINT_FRONT
      //       cout<<"\tBefore "
      // 	  <<TMath::Abs(fChPoint[2][ch2CellPoint[itrackletback]]->fX - fChPoint[1][ch1CellPoint[itrackletfront]]->fX)
      // 	  <<"\t"
      // 	  <<TMath::Abs(fChPoint[2][ch2CellPoint[itrackletback]]->fY - fChPoint[1][ch1CellPoint[itrackletfront]]->fY)
      // 	  <<endl;
      // #endif
      if(TMath::Abs(fChPoint[2][ch2CellPoint[itrackletback]]->fX - 
		    fChPoint[1][ch1CellPoint[itrackletfront]]->fX) > diffDistStX || 
	 TMath::Abs(fChPoint[2][ch2CellPoint[itrackletback]]->fY - 
		    fChPoint[1][ch1CellPoint[itrackletfront]]->fY) > diffDistStY ) continue;

      meanX2 = (fChPoint[2][ch2CellPoint[itrackletback]]->fX 
		+ fChPoint[3][ch3CellPoint[itrackletback]]->fX)/2.0 ;
      meanY2 = (fChPoint[2][ch2CellPoint[itrackletback]]->fY 
		+ fChPoint[3][ch3CellPoint[itrackletback]]->fY)/2.0 ;
      meanZ2 = (fChPoint[2][ch2CellPoint[itrackletback]]->fZ 
		+ fChPoint[3][ch3CellPoint[itrackletback]]->fZ)/2.0 ;
      
      expectSt3X = meanX2 + (fgkTrackDetCoordinate[2]-meanZ2)*(meanX2-meanX1)/(meanZ2-meanZ1);
      expectSt3Y = meanY2 + (fgkTrackDetCoordinate[2]-meanZ2)*(meanY2-meanY1)/(meanZ2-meanZ1);
      inclinationFront = (meanX2-meanX1)/(meanZ2-meanZ1) ;
      
      for( Int_t ibacktrackseg=0;ibacktrackseg<fNofbackTrackSeg;ibacktrackseg++){  
	
	if( /// TMath::Abs(inclinationBack[ibacktrackseg]-inclinationFront)<0.04 && 
	   TMath::Abs((expectSt3X-fExtrapSt3X[ibacktrackseg])) < st3WindowX
	   && TMath::Abs((expectSt3Y-fExtrapSt3Y[ibacktrackseg])) < st3WindowY){
	  
	  Sub(fChPoint[2][ch2CellPoint[itrackletback]],fChPoint[1][ch1CellPoint[itrackletfront]],&pSeg2);
	  Sub(fChPoint[3][ch3CellPoint[itrackletback]],fChPoint[2][ch2CellPoint[itrackletback]],&pSeg3);
	  
	  anglediff = TMath::RadToDeg()* (Angle(&pSeg1,&pSeg2) + Angle(&pSeg2,&pSeg3));

	  // #ifdef PRINT_FRONT
	  // 	  printf("\t\t\tanglediff : %lf\n",anglediff);
	  // #endif	  
	  if(anglediff<minAngle){
	    minAngle = anglediff;
	    indexMinAngleBack = itrackletback;
	    indexMinAngleFront = itrackletfront;
	    backIndex = ibacktrackseg;
	    isConnected[ibacktrackseg] = true;
	  }
	}///matching tracklet
      }///for loop of backtrackseg
      
    }///for of back ch
    
    if(minAngle < minAngleWindow && indexMinAngleFront!=-1 
       && indexMinAngleBack!=-1 && fNoffrontTrackSeg<(fgkMaxNofTracks-1)){
      
      fFrontTrackSeg[fNoffrontTrackSeg].fIndex[0] = ch0CellPoint[indexMinAngleFront];
      fFrontTrackSeg[fNoffrontTrackSeg].fIndex[1] = ch1CellPoint[indexMinAngleFront];
      fFrontTrackSeg[fNoffrontTrackSeg].fIndex[2] = ch2CellPoint[indexMinAngleBack];
      fFrontTrackSeg[fNoffrontTrackSeg].fIndex[3] = ch3CellPoint[indexMinAngleBack];
      
      fBackToFront[backIndex][fNofConnectedfrontTrackSeg[backIndex]++] = fNoffrontTrackSeg;
      fNoffrontTrackSeg++;
    }///condition to find valid tracklet
    
  }///for loop of front ch
  
  Int_t nofNCfBackTrackSeg = 0;
  Int_t ncfBackTrackSeg[fgkMaxNofTracks];
  
  for( Int_t ibacktrackseg=0;ibacktrackseg<fNofbackTrackSeg;ibacktrackseg++)
    if(!isConnected[ibacktrackseg]) 
      ncfBackTrackSeg[nofNCfBackTrackSeg++] = ibacktrackseg;
    else
      fNofConnected++;
  

#ifdef PRINT_FRONT
  printf("\tfNofConnected : %d, nofNCfBackTrackSeg : %d\n",fNofConnected,nofNCfBackTrackSeg);
  printf("\tfNofPoints[3] : %d, fNofPoints[2] : %d\n",fNofPoints[3],fNofPoints[2]);
  if(nofNCfBackTrackSeg==0){
    printf("All fBackTrackSegs are connected with fFrontTrackSegs, no need to search further\n");
    printf("AliHLTMUONFullTracker::QuadTrackSeg()--End\n\n");
  }
#endif
  
  if(nofNCfBackTrackSeg==0) return true;


  ///Next Check for tracklet in Stn1 and space point in Stn2
  Bool_t isbackpoint=false,isfrontpoint=false;
  for( Int_t itrackletfront=0;itrackletfront<nofSt1Cells;itrackletfront++){
    Sub(fChPoint[1][ch1CellPoint[itrackletfront]],fChPoint[0][ch0CellPoint[itrackletfront]],&pSeg1);
    minAngle = minAngleWindow;
    indexMinAngleBack = -1;
    indexMinAngleFront = -1;

    for( Int_t ibackpoint=0;ibackpoint<fNofPoints[3];ibackpoint++){
      if(/// hasCh3Cells[ibackpoint] == true && 
	 TMath::Abs(fChPoint[3][ibackpoint]->fX - 
		    fChPoint[1][ch1CellPoint[itrackletfront]]->fX) > diffDistStX || 
	 TMath::Abs(fChPoint[3][ibackpoint]->fY - 
		    fChPoint[1][ch1CellPoint[itrackletfront]]->fY) > diffDistStY ) continue;
      
      expectSt3X = fChPoint[3][ibackpoint]->fX + (fgkTrackDetCoordinate[2] - fChPoint[3][ibackpoint]->fZ)*
	(fChPoint[3][ibackpoint]->fX - fChPoint[1][ch1CellPoint[itrackletfront]]->fX)/
	(fChPoint[3][ibackpoint]->fZ - fChPoint[1][ch1CellPoint[itrackletfront]]->fZ);
      expectSt3Y = fChPoint[3][ibackpoint]->fY + (fgkTrackDetCoordinate[2] - fChPoint[3][ibackpoint]->fZ)*
	(fChPoint[3][ibackpoint]->fY - fChPoint[1][ch1CellPoint[itrackletfront]]->fY)/
	(fChPoint[3][ibackpoint]->fZ - fChPoint[1][ch1CellPoint[itrackletfront]]->fZ);
      inclinationFront = (fChPoint[3][ibackpoint]->fX - fChPoint[1][ch1CellPoint[itrackletfront]]->fX)/
	(fChPoint[3][ibackpoint]->fZ - fChPoint[1][ch1CellPoint[itrackletfront]]->fZ) ;
      
      for( Int_t ibacktrackseg=0;ibacktrackseg<nofNCfBackTrackSeg;ibacktrackseg++){  
	
	if(/// TMath::Abs(inclinationBack[ncfBackTrackSeg[ibacktrackseg]]-inclinationFront)< inclinationWindow && 
	   TMath::Abs((expectSt3X-fExtrapSt3X[ncfBackTrackSeg[ibacktrackseg]])) < st3WindowX
	   && TMath::Abs((expectSt3Y-fExtrapSt3Y[ncfBackTrackSeg[ibacktrackseg]])) < st3WindowY){
	  
	  Sub(fChPoint[3][ibackpoint],fChPoint[1][ch1CellPoint[itrackletfront]],&pSeg2);
	  
	  anglediff = TMath::RadToDeg()* Angle(&pSeg1,&pSeg2) ;
	  // #ifdef PRINT_FRONT
	  // 	  printf("\t\t annglediff(Ch4) : %lf\n",anglediff);
	  // #endif	  
	  if(anglediff<minAngle){
	    minAngle = anglediff;
	    indexMinAngleBack = ibackpoint;
	    indexMinAngleFront = itrackletfront;
	    backIndex = ncfBackTrackSeg[ibacktrackseg];
	    isConnected[ncfBackTrackSeg[ibacktrackseg]] = true;
	    isbackpoint = true;
	    isfrontpoint = false;
	  }///angle min condn
	}///matching tracklet
      }///loop on not found back trackseg
    }///backpoint loop
    
    for( Int_t ifrontpoint=0;ifrontpoint<fNofPoints[2];ifrontpoint++){
      if(/// hasCh2Cells[ifrontpoint] == true && 
	 TMath::Abs(fChPoint[2][ifrontpoint]->fX - 
		    fChPoint[1][ch1CellPoint[itrackletfront]]->fX) > diffDistStX || 
	 TMath::Abs(fChPoint[2][ifrontpoint]->fY - 
		    fChPoint[1][ch1CellPoint[itrackletfront]]->fY) > diffDistStY ) continue;
      
      expectSt3X = fChPoint[2][ifrontpoint]->fX + (fgkTrackDetCoordinate[2] - fChPoint[2][ifrontpoint]->fZ)*
	(fChPoint[2][ifrontpoint]->fX - fChPoint[1][ch1CellPoint[itrackletfront]]->fX)/
	(fChPoint[2][ifrontpoint]->fZ - fChPoint[1][ch1CellPoint[itrackletfront]]->fZ);
      expectSt3Y = fChPoint[2][ifrontpoint]->fY + (fgkTrackDetCoordinate[2] - fChPoint[2][ifrontpoint]->fZ)*
	(fChPoint[2][ifrontpoint]->fY - fChPoint[1][ch1CellPoint[itrackletfront]]->fY)/
	(fChPoint[2][ifrontpoint]->fZ - fChPoint[1][ch1CellPoint[itrackletfront]]->fZ);
      inclinationFront = (fChPoint[2][ifrontpoint]->fX - fChPoint[1][ch1CellPoint[itrackletfront]]->fX)/
	(fChPoint[2][ifrontpoint]->fZ - fChPoint[1][ch1CellPoint[itrackletfront]]->fZ) ;
      // #ifdef PRINT_FRONT
      //       printf("\t\texpectSt3X : %f, expectSt3Y : %f, inclinationFront : %f\n",expectSt3X,expectSt3Y,inclinationFront);
      // #endif      
      
      for( Int_t ibacktrackseg=0;ibacktrackseg<nofNCfBackTrackSeg;ibacktrackseg++){  
	
	if( /// TMath::Abs(inclinationBack[ncfBackTrackSeg[ibacktrackseg]]-inclinationFront)< inclinationWindow && 
	   TMath::Abs((expectSt3X-fExtrapSt3X[ncfBackTrackSeg[ibacktrackseg]])) < st3WindowX
	   && TMath::Abs((expectSt3Y-fExtrapSt3Y[ncfBackTrackSeg[ibacktrackseg]])) < st3WindowY){
	  
	  Sub(fChPoint[2][ifrontpoint],fChPoint[1][ch1CellPoint[itrackletfront]],&pSeg2);

	  anglediff = TMath::RadToDeg()* Angle(&pSeg1,&pSeg2) ;
	  // #ifdef PRINT_FRONT
	  // 	  printf("\t\t annglediff(Ch3) : %lf\n",anglediff);
	  // #endif	  
	  if(anglediff<minAngle){
	    minAngle = anglediff;
	    indexMinAngleBack = ifrontpoint;
	    indexMinAngleFront = itrackletfront;
	    backIndex = ncfBackTrackSeg[ibacktrackseg];
	    isConnected[ncfBackTrackSeg[ibacktrackseg]] = true;
	    isbackpoint = false;
	    isfrontpoint = true;

	  }///angle min condn
	}///matching tracklet
      }///loop on not found back trackseg
    }///backpoint loop
    
    if(minAngle < minAngleWindow && indexMinAngleFront!=-1 && 
       indexMinAngleBack!=-1 && fNoffrontTrackSeg<(fgkMaxNofTracks-1)){
      
      fFrontTrackSeg[fNoffrontTrackSeg].fIndex[0] = ch0CellPoint[indexMinAngleFront];
      fFrontTrackSeg[fNoffrontTrackSeg].fIndex[1] = ch1CellPoint[indexMinAngleFront];
      if(isfrontpoint && !isbackpoint){
	fFrontTrackSeg[fNoffrontTrackSeg].fIndex[2] = indexMinAngleBack;
	fFrontTrackSeg[fNoffrontTrackSeg].fIndex[3] = -1;
      }else{
	fFrontTrackSeg[fNoffrontTrackSeg].fIndex[2] = -1;
	fFrontTrackSeg[fNoffrontTrackSeg].fIndex[3] = indexMinAngleBack;
      }
      fBackToFront[backIndex][fNofConnectedfrontTrackSeg[backIndex]++] = fNoffrontTrackSeg;      
      fNoffrontTrackSeg++;
    }///condition to find valid tracklet
    
  }///front ch


  Int_t nofSNCfBackTrackSeg = 0;
  Int_t sncfBackTrackSeg[fgkMaxNofTracks];
  
  for( Int_t ibacktrackseg=0;ibacktrackseg<nofNCfBackTrackSeg;ibacktrackseg++)
    if(!isConnected[ncfBackTrackSeg[ibacktrackseg]]) 
      sncfBackTrackSeg[nofSNCfBackTrackSeg++] = ncfBackTrackSeg[ibacktrackseg];
    else
      fNofConnected++;

#ifdef PRINT_FRONT
  printf("\tfNofConnected : %d, nofSNCfBackTrackSeg : %d\n",fNofConnected,nofSNCfBackTrackSeg);
  if(nofSNCfBackTrackSeg==0){
    printf("All fBackTrackSegs are connected with fFrontTrackSegs, no need to search further\n");
    printf("AliHLTMUONFullTracker::QuadTrackSeg()--End\n\n");
  }
#endif

  if(nofSNCfBackTrackSeg==0) return true;

  ///Last Check for tracklet in Stn2 and space point in Stn1
  for( Int_t itrackletback=0;itrackletback<nofSt2Cells;itrackletback++){
    Sub(fChPoint[3][ch3CellPoint[itrackletback]],fChPoint[2][ch2CellPoint[itrackletback]],&pSeg1);
    minAngle = minAngleWindow ;
    indexMinAngleBack = -1;
    indexMinAngleFront = -1;

    for( Int_t ibackpoint=0;ibackpoint<fNofPoints[1];ibackpoint++){
      if(/// hasCh1Cells[ibackpoint] == true && 
	 TMath::Abs(fChPoint[2][ch2CellPoint[itrackletback]]->fX - 
		    fChPoint[1][ibackpoint]->fX) > diffDistStX || 
	 TMath::Abs(fChPoint[2][ch2CellPoint[itrackletback]]->fY - 
		    fChPoint[1][ibackpoint]->fY) > diffDistStY) continue;
      
      expectSt3X = fChPoint[2][ch2CellPoint[itrackletback]]->fX + (fgkTrackDetCoordinate[2] - fChPoint[2][ch2CellPoint[itrackletback]]->fZ)*
	(fChPoint[2][ch2CellPoint[itrackletback]]->fX - fChPoint[1][ibackpoint]->fX)/
	(fChPoint[2][ch2CellPoint[itrackletback]]->fZ - fChPoint[1][ibackpoint]->fZ);
      expectSt3Y = fChPoint[2][ch2CellPoint[itrackletback]]->fY + (fgkTrackDetCoordinate[2] - fChPoint[2][ch2CellPoint[itrackletback]]->fZ)*
	(fChPoint[2][ch2CellPoint[itrackletback]]->fY - fChPoint[1][ibackpoint]->fY)/
	(fChPoint[2][ch2CellPoint[itrackletback]]->fZ - fChPoint[1][ibackpoint]->fZ);
      inclinationFront = (fChPoint[2][ch2CellPoint[itrackletback]]->fX - fChPoint[1][ibackpoint]->fX)/
	(fChPoint[2][ch2CellPoint[itrackletback]]->fZ - fChPoint[1][ibackpoint]->fZ) ;
      
      for( Int_t ibacktrackseg=0;ibacktrackseg<nofSNCfBackTrackSeg;ibacktrackseg++){  
	
	if(///  TMath::Abs(inclinationBack[sncfBackTrackSeg[ibacktrackseg]]-inclinationFront)<inclinationWindow && 
	   TMath::Abs((expectSt3X-fExtrapSt3X[sncfBackTrackSeg[ibacktrackseg]])) < st3WindowX
	   && TMath::Abs((expectSt3Y-fExtrapSt3Y[sncfBackTrackSeg[ibacktrackseg]])) < st3WindowY ){
	  
	  Sub(fChPoint[2][ch2CellPoint[itrackletback]],fChPoint[1][ibackpoint],&pSeg2);
	  
	  anglediff = TMath::RadToDeg()* Angle(&pSeg1,&pSeg2) ;
	  if(anglediff<minAngle){
	    minAngle = anglediff;
	    indexMinAngleBack = itrackletback;
	    indexMinAngleFront = ibackpoint;
	    backIndex = sncfBackTrackSeg[ibacktrackseg];
	    isConnected[sncfBackTrackSeg[ibacktrackseg]] = true;
	    isbackpoint = true;
	    isfrontpoint = false;
	  }///angle min condn
	}///matching tracklet
      }///loop on not found back trackseg
    }///backpoint loop
    
    for( Int_t ifrontpoint=0;ifrontpoint<fNofPoints[0];ifrontpoint++){
      if(/// hasCh0Cells[ifrontpoint] == true &&
	 TMath::Abs(fChPoint[2][ch2CellPoint[itrackletback]]->fX - 
		    fChPoint[0][ifrontpoint]->fX) > diffDistStX || 
	 TMath::Abs(fChPoint[2][ch2CellPoint[itrackletback]]->fY - 
		    fChPoint[0][ifrontpoint]->fY) > diffDistStY ) continue;
      
      expectSt3X = fChPoint[2][ch2CellPoint[itrackletback]]->fX + (fgkTrackDetCoordinate[2] - fChPoint[2][ch2CellPoint[itrackletback]]->fZ)*
	(fChPoint[2][ch2CellPoint[itrackletback]]->fX - fChPoint[0][ifrontpoint]->fX)/
	(fChPoint[2][ch2CellPoint[itrackletback]]->fZ - fChPoint[0][ifrontpoint]->fZ);
      expectSt3Y = fChPoint[2][ch2CellPoint[itrackletback]]->fY + (fgkTrackDetCoordinate[2] - fChPoint[2][ch2CellPoint[itrackletback]]->fZ)*
	(fChPoint[2][ch2CellPoint[itrackletback]]->fY - fChPoint[0][ifrontpoint]->fY)/
	(fChPoint[2][ch2CellPoint[itrackletback]]->fZ - fChPoint[0][ifrontpoint]->fZ);
      inclinationFront = (fChPoint[2][ch2CellPoint[itrackletback]]->fX - fChPoint[0][ifrontpoint]->fX)/
	(fChPoint[2][ch2CellPoint[itrackletback]]->fZ - fChPoint[0][ifrontpoint]->fZ) ;
      
      for( Int_t ibacktrackseg=0;ibacktrackseg<nofSNCfBackTrackSeg;ibacktrackseg++){  

	if(///  TMath::Abs(inclinationBack[sncfBackTrackSeg[ibacktrackseg]]-inclinationFront)<inclinationWindow && 
	   TMath::Abs((expectSt3X-fExtrapSt3X[sncfBackTrackSeg[ibacktrackseg]])) < st3WindowX
	   && TMath::Abs((expectSt3Y-fExtrapSt3Y[sncfBackTrackSeg[ibacktrackseg]])) < st3WindowY ){
	  
	  Sub(fChPoint[2][ch2CellPoint[itrackletback]],fChPoint[0][ifrontpoint],&pSeg2);
	  
	  anglediff = TMath::RadToDeg() * Angle(&pSeg1,&pSeg2) ;
	  if(anglediff<minAngle){
	    minAngle = anglediff;
	    indexMinAngleBack = itrackletback;
	    indexMinAngleFront = ifrontpoint;
	    backIndex = sncfBackTrackSeg[ibacktrackseg];
	    isConnected[sncfBackTrackSeg[ibacktrackseg]] = true;
	    isbackpoint = false;
	    isfrontpoint = true;
	    
	  }///angle min condn
	}///matching tracklet
      }///loop on not found back trackseg
    }///backpoint loop
    
    if(minAngle < minAngleWindow && indexMinAngleFront!=-1 && 
       indexMinAngleBack!=-1 && fNoffrontTrackSeg<(fgkMaxNofTracks-1)){

      if(isfrontpoint){      
	fFrontTrackSeg[fNoffrontTrackSeg].fIndex[0] = indexMinAngleFront;
	fFrontTrackSeg[fNoffrontTrackSeg].fIndex[1] = -1;
      }else{
	fFrontTrackSeg[fNoffrontTrackSeg].fIndex[0] = -1;
	fFrontTrackSeg[fNoffrontTrackSeg].fIndex[1] = indexMinAngleFront;
      }

      fFrontTrackSeg[fNoffrontTrackSeg].fIndex[2] = ch2CellPoint[indexMinAngleBack];
      fFrontTrackSeg[fNoffrontTrackSeg].fIndex[3] = ch3CellPoint[indexMinAngleBack];

      fBackToFront[backIndex][fNofConnectedfrontTrackSeg[backIndex]++] = fNoffrontTrackSeg;
      fNoffrontTrackSeg++;
    }///condition to find valid tracklet
    
  }///front ch
  
  for( Int_t ibacktrackseg=0;ibacktrackseg<nofSNCfBackTrackSeg;ibacktrackseg++)
    if(isConnected[sncfBackTrackSeg[ibacktrackseg]]) 
      fNofConnected++;
  
#ifdef PRINT_FRONT
  printf("\tfNofConnected : %d\n",fNofConnected);
  printf("Three spacepoints are found in fFrontTrackSegs\n");
  printf("AliHLTMUONFullTracker::QuadTrackSeg()--End\n\n");
#endif


  return true;
}

///__________________________________________________________________________

Double_t AliHLTMUONFullTracker::KalmanFilter(AliMUONTrackParam &trackParamAtCluster, Cluster *cluster)
{
  //// Compute new track parameters and their covariances including new cluster using kalman filter
  //// return the additional track chi2

#ifdef PRINT_DETAIL_KALMAN
  printf("AliHLTMUONFullTracker::KalmanFilter()--Begin\n\n");
#endif


  /// Get actual track parameters (p)
  TMatrixD param(trackParamAtCluster.GetParameters());
#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","param.Print() [p]");
  param.Print();
  printf("GetZ : %lf\n",trackParamAtCluster.GetZ());
#endif



  /// Get new cluster parameters (m)
  ///AliMUONVCluster *cluster = trackParamAtCluster.GetClusterPtr();
  TMatrixD clusterParam(5,1);
  clusterParam.Zero();
  ///   clusterParam(0,0) = cluster->GetX();
  ///   clusterParam(2,0) = cluster->GetY();
  clusterParam(0,0) = cluster->fX;
  clusterParam(2,0) = cluster->fY;
#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","clusterParam.Print() [m]");
  clusterParam.Print();
#endif



  /// Compute the actual parameter weight (W)
  TMatrixD paramWeight(trackParamAtCluster.GetCovariances());
#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","Covariance : [C]");
  paramWeight.Print();
#endif

  if (paramWeight.Determinant() != 0) {
    paramWeight.Invert();
  } else {
    Warning("KalmanFilter"," Determinant = 0");
    return 1.0e10;
  }

#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","Weight Matrix inverse of Covariance [W = C^-1]");
  paramWeight.Print();
#endif


  /// Compute the new cluster weight (U)
  TMatrixD clusterWeight(5,5);
  clusterWeight.Zero();
  clusterWeight(0,0) = 1. / cluster->fErrX2;
  clusterWeight(2,2) = 1. / cluster->fErrY2;
#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","clusterWeight.Print() [U]");
  printf("\tErrX2 : %lf, ErrY2 : %lf\n",cluster->fErrX2,cluster->fErrY2);
  clusterWeight.Print();
#endif




  /// Compute the new parameters covariance matrix ( (W+U)^-1 )
  TMatrixD newParamCov(paramWeight,TMatrixD::kPlus,clusterWeight);
#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","newParamCov.Print() [(W+U)]");
  newParamCov.Print();
#endif
  if (newParamCov.Determinant() != 0) {
    newParamCov.Invert();
  } else {
    Warning("RunKalmanFilter"," Determinant = 0");
    return 1.0e10;
  }
#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","newParamCov.Print() [(W+U)^-1] (new covariances[W] for trackParamAtCluster)");
  newParamCov.Print();
#endif

  /// Save the new parameters covariance matrix
  trackParamAtCluster.SetCovariances(newParamCov);
  
  /// Compute the new parameters (p' = ((W+U)^-1)U(m-p) + p)
  TMatrixD tmp(clusterParam,TMatrixD::kMinus,param);
  TMatrixD tmp2(clusterWeight,TMatrixD::kMult,tmp); /// U(m-p)
  TMatrixD newParam(newParamCov,TMatrixD::kMult,tmp2); /// ((W+U)^-1)U(m-p)
  newParam += param; /// ((W+U)^-1)U(m-p) + p
#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","newParam.Print() [p' = ((W+U)^-1)U(m-p) + p] (new parameters[p] for trackParamAtCluster)");
  newParam.Print();
#endif

  /// Save the new parameters
  trackParamAtCluster.SetParameters(newParam);
  ///   printf(Form("Pt : %lf\n",TMath::Sqrt(trackParamAtCluster.Px()*trackParamAtCluster.Px() + 
  /// 				      trackParamAtCluster.Py()*trackParamAtCluster.Py())));
  
  
  /// Compute the additional chi2 (= ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m))
  tmp = newParam; /// p'
  tmp -= param; /// (p'-p)
  TMatrixD tmp3(paramWeight,TMatrixD::kMult,tmp); /// W(p'-p)
  TMatrixD addChi2Track(tmp,TMatrixD::kTransposeMult,tmp3); /// ((p'-p)^-1)W(p'-p)
  tmp = newParam; /// p'
  tmp -= clusterParam; /// (p'-m)
  TMatrixD tmp4(clusterWeight,TMatrixD::kMult,tmp); /// U(p'-m)
  addChi2Track += TMatrixD(tmp,TMatrixD::kTransposeMult,tmp4); /// ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)

#ifdef PRINT_DETAIL_KALMAN
  Info("\tKalmanFilter","addChi2Track.Print() [additional chi2 = ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)))]");
  addChi2Track.Print();
  printf("AliHLTMUONFullTracker::KalmanFilter()--End\n\n");
#endif


  return addChi2Track(0,0);
}

///__________________________________________________________________________

Double_t AliHLTMUONFullTracker::TryOneCluster(const AliMUONTrackParam &trackParam, Cluster* cluster,
					      AliMUONTrackParam &trackParamAtCluster, Bool_t updatePropagator)
{
  //// Test the compatibility between the track and the cluster (using trackParam's covariance matrix):
  //// return the corresponding Chi2
  //// return trackParamAtCluster
  
  /// extrapolate track parameters and covariances at the z position of the tested cluster
  trackParamAtCluster = trackParam;
  AliMUONTrackExtrap::ExtrapToZCov(&trackParamAtCluster, cluster->fZ, updatePropagator);
  
  /// set pointer to cluster into trackParamAtCluster
  ///trackParamAtCluster.SetClusterPtr(cluster);
  
  /// Set differences between trackParam and cluster in the bending and non bending directions
  Double_t dX = cluster->fX - trackParamAtCluster.GetNonBendingCoor();
  Double_t dY = cluster->fY - trackParamAtCluster.GetBendingCoor();
  
  /// Calculate errors and covariances
  const TMatrixD& kParamCov = trackParamAtCluster.GetCovariances();
  Double_t sigmaX2 = kParamCov(0,0) + cluster->fErrX2;
  Double_t sigmaY2 = kParamCov(2,2) + cluster->fErrY2;
  
  /// Compute chi2
  return dX * dX / sigmaX2 + dY * dY / sigmaY2;
  
}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::TryOneClusterFast(const AliMUONTrackParam &trackParam, const Cluster* cluster)
{
  //// Test the compatibility between the track and the cluster
  //// given the track resolution + the maximum-distance-to-track value
  //// and assuming linear propagation of the track:
  //// return kTRUE if they are compatibles
  
  Float_t sigmaCutForTracking = 6.0;
  Float_t maxNonBendingDistanceToTrack = 1.0;
  Float_t maxBendingDistanceToTrack = 1.0;
  
  Double_t dZ = cluster->fZ - trackParam.GetZ();
  Double_t dX = cluster->fX - (trackParam.GetNonBendingCoor() + trackParam.GetNonBendingSlope() * dZ);
  Double_t dY = cluster->fY - (trackParam.GetBendingCoor() + trackParam.GetBendingSlope() * dZ);
  const TMatrixD& kParamCov = trackParam.GetCovariances();
  Double_t errX2 = kParamCov(0,0) + dZ * dZ * kParamCov(1,1) + 2. * dZ * kParamCov(0,1);
  Double_t errY2 = kParamCov(2,2) + dZ * dZ * kParamCov(3,3) + 2. * dZ * kParamCov(2,3);

  Double_t dXmax = sigmaCutForTracking * TMath::Sqrt(errX2) +
    maxNonBendingDistanceToTrack;
  Double_t dYmax = sigmaCutForTracking * TMath::Sqrt(errY2) +
    maxBendingDistanceToTrack;
  
  if (TMath::Abs(dX) > dXmax || TMath::Abs(dY) > dYmax) return kFALSE;
  
  return kTRUE;
  
}
///__________________________________________________________________________

void AliHLTMUONFullTracker::PropagateTracks(Double_t charge, Float_t& px, Float_t& py, Float_t& pz, 
					    Float_t& xr, Float_t& yr, Float_t& zr, Float_t zprop)
{
  ///
  /// propagate in magnetic field between hits of indices i1 and i2
  ///

  Double_t vect[7], vout[7];
  Double_t step = -5.0;
  Double_t zMax = zprop+.5;
  
  vect[0] = xr;
  vect[1] = yr;
  vect[2] = zr;
  vect[6] = sqrt((px)*(px) + (py)*(py) + (pz)*(pz));
  vect[3] = px/vect[6];
  vect[4] = py/vect[6];
  vect[5] = pz/vect[6];
  ///   cout<<"vec[2] : "<<vect[2]<<", zMax : "<<zMax<<endl;

  Int_t nSteps = 0;
  while ((vect[2] < zMax) && (nSteps < 10000)) {
    nSteps++;
    ///     OneStepRungekutta(charge, step, vect, vout);
    OneStepHelix3(charge,step,vect,vout);
    ///SetPoint(fCount,vout[0],vout[1],vout[2]);
    ///fCount++;
    ///     printf("(x,y,z) : (%f,%f,%f)\n",vout[0],vout[1],vout[2]);
    for (Int_t i = 0; i < 7; i++) {
      vect[i] = vout[i];
    }
  }
  
  xr = vect[0];
  yr = vect[1];
  zr = vect[2];
  
  px = vect[3]*vect[6];
  py = vect[4]*vect[6];
  pz = vect[5]*vect[6];
  return;
}

///__________________________________________________________________________

void AliHLTMUONFullTracker::OneStepHelix3(Double_t field, Double_t step, const Double_t *vect, Double_t *vout) const
{
  //// <pre>
  ////	******************************************************************
  ////	*								 *
  ////	*	Tracking routine in a constant field oriented		 *
  ////	*	along axis 3						 *
  ////	*	Tracking is performed with a conventional		 *
  ////	*	helix step method					 *
  ////	*								 *
  ////	*    ==>Called by : USER, GUSWIM				 *
  ////	*	Authors    R.Brun, M.Hansroul  *********		 *
  ////	*	Rewritten  V.Perevoztchikov
  ////	*								 *
  ////	******************************************************************
  //// </pre>

  Double_t hxp[3];
  Double_t h4, hp, rho, tet;
  Double_t sint, sintt, tsint, cos1t;
  Double_t f1, f2, f3, f4, f5, f6;

  const Int_t kix  = 0;
  const Int_t kiy  = 1;
  const Int_t kiz  = 2;
  const Int_t kipx = 3;
  const Int_t kipy = 4;
  const Int_t kipz = 5;
  const Int_t kipp = 6;

  const Double_t kec = 2.9979251e-4;

  /// 
  ///     ------------------------------------------------------------------
  /// 
  ///       units are kgauss,centimeters,gev/c
  /// 
  vout[kipp] = vect[kipp];
  h4 = field * kec;

  hxp[0] = - vect[kipy];
  hxp[1] = + vect[kipx];
 
  hp = vect[kipz];

  rho = -h4/vect[kipp];
  tet = rho * step;
  if (TMath::Abs(tet) > 0.15) {
    sint = TMath::Sin(tet);
    sintt = (sint/tet);
    tsint = (tet-sint)/tet;
    cos1t = 2.* TMath::Sin(0.5*tet) * TMath::Sin(0.5*tet)/tet;
  } else {
    tsint = tet*tet/36.;
    sintt = (1. - tsint);
    sint = tet*sintt;
    cos1t = 0.5*tet;
  }

  f1 = step * sintt;
  f2 = step * cos1t;
  f3 = step * tsint * hp;
  f4 = -tet*cos1t;
  f5 = sint;
  f6 = tet * cos1t * hp;
 
  vout[kix] = vect[kix] + f1*vect[kipx] + f2*hxp[0];
  vout[kiy] = vect[kiy] + f1*vect[kipy] + f2*hxp[1];
  vout[kiz] = vect[kiz] + f1*vect[kipz] + f3;
 
  vout[kipx] = vect[kipx] + f4*vect[kipx] + f5*hxp[0];
  vout[kipy] = vect[kipy] + f4*vect[kipy] + f5*hxp[1];
  vout[kipz] = vect[kipz] + f4*vect[kipz] + f6;

  return;
}

///__________________________________________________________________________

///______________________________________________________________________________

void AliHLTMUONFullTracker::OneStepRungekutta(Double_t charge, Double_t step,
					      const Double_t* vect, Double_t* vout)
{
  ////	******************************************************************
  ////	*								 *
  ////	*  Runge-Kutta method for tracking a particle through a magnetic *
  ////	*  field. Uses Nystroem algorithm (See Handbook Nat. Bur. of	 *
  ////	*  Standards, procedure 25.5.20)				 *
  ////	*								 *
  ////	*  Input parameters						 *
  ////	*	CHARGE    Particle charge				 *
  ////	*	STEP	  Step size					 *
  ////	*	VECT	  Initial co-ords,direction cosines,momentum	 *
  ////	*  Output parameters						 *
  ////	*	VOUT	  Output co-ords,direction cosines,momentum	 *
  ////	*  User routine called  					 *
  ////	*	CALL GUFLD(X,F) 					 *
  ////	*								 *
  ////	*    ==>Called by : <USER>, GUSWIM				 *
  ////	*	Authors    R.Brun, M.Hansroul  *********		 *
  ////	*		   V.Perevoztchikov (CUT STEP implementation)	 *
  ////	*								 *
  ////	*								 *
  ////	******************************************************************

  Double_t h2, h4, f[4];
  Double_t xyzt[3], a, b, c, ph,ph2;
  Double_t secxs[4],secys[4],seczs[4],hxp[3];
  Double_t g1, g2, g3, g4, g5, g6, ang2, dxt, dyt, dzt;
  Double_t est, at, bt, ct, cba;
  Double_t f1, f2, f3, f4, rho, tet, hnorm, hp, rho1, sint, cost;

  Double_t x;
  Double_t y;
  Double_t z;

  Double_t xt;
  Double_t yt;
  Double_t zt;

  Double_t maxit = 1992;
  Double_t maxcut = 11;

  const Double_t kdlt   = 1e-4;
  const Double_t kdlt32 = kdlt/32.;
  const Double_t kthird = 1./3.;
  const Double_t khalf  = 0.5;
  const Double_t kec = 2.9979251e-4;

  const Double_t kpisqua = 9.86960440109;
  const Int_t kix  = 0;
  const Int_t kiy  = 1;
  const Int_t kiz  = 2;
  const Int_t kipx = 3;
  const Int_t kipy = 4;
  const Int_t kipz = 5;

  /// *.
  /// *.    ------------------------------------------------------------------
  /// *.
  /// *             this constant is for units cm,gev/c and kgauss
  /// *
  Int_t iter = 0;
  Int_t ncut = 0;
  for(Int_t j = 0; j < 7; j++)
    vout[j] = vect[j];

  Double_t  pinv   = kec * charge / vect[6];
  Double_t tl = 0.;
  Double_t h = step;
  Double_t rest;


  do {
    rest  = step - tl;
    if (TMath::Abs(h) > TMath::Abs(rest)) h = rest;
    ///cmodif: call gufld(vout,f) changed into:
    TGeoGlobalMagField::Instance()->Field(vout,f);

    /// *
    /// *             start of integration
    /// *
    x      = vout[0];
    y      = vout[1];
    z      = vout[2];
    a      = vout[3];
    b      = vout[4];
    c      = vout[5];

    h2     = khalf * h;
    h4     = khalf * h2;
    ph     = pinv * h;
    ph2    = khalf * ph;
    secxs[0] = (b * f[2] - c * f[1]) * ph2;
    secys[0] = (c * f[0] - a * f[2]) * ph2;
    seczs[0] = (a * f[1] - b * f[0]) * ph2;
    ang2 = (secxs[0]*secxs[0] + secys[0]*secys[0] + seczs[0]*seczs[0]);
    if (ang2 > kpisqua) break;

    dxt    = h2 * a + h4 * secxs[0];
    dyt    = h2 * b + h4 * secys[0];
    dzt    = h2 * c + h4 * seczs[0];
    xt     = x + dxt;
    yt     = y + dyt;
    zt     = z + dzt;
    /// *
    /// *              second intermediate point
    /// *

    est = TMath::Abs(dxt) + TMath::Abs(dyt) + TMath::Abs(dzt);
    if (est > h) {
      if (ncut++ > maxcut) break;
      h *= khalf;
      continue;
    }

    xyzt[0] = xt;
    xyzt[1] = yt;
    xyzt[2] = zt;

    ///cmodif: call gufld(xyzt,f) changed into:
    TGeoGlobalMagField::Instance()->Field(xyzt,f);

    at     = a + secxs[0];
    bt     = b + secys[0];
    ct     = c + seczs[0];

    secxs[1] = (bt * f[2] - ct * f[1]) * ph2;
    secys[1] = (ct * f[0] - at * f[2]) * ph2;
    seczs[1] = (at * f[1] - bt * f[0]) * ph2;
    at     = a + secxs[1];
    bt     = b + secys[1];
    ct     = c + seczs[1];
    secxs[2] = (bt * f[2] - ct * f[1]) * ph2;
    secys[2] = (ct * f[0] - at * f[2]) * ph2;
    seczs[2] = (at * f[1] - bt * f[0]) * ph2;
    dxt    = h * (a + secxs[2]);
    dyt    = h * (b + secys[2]);
    dzt    = h * (c + seczs[2]);
    xt     = x + dxt;
    yt     = y + dyt;
    zt     = z + dzt;
    at     = a + 2.*secxs[2];
    bt     = b + 2.*secys[2];
    ct     = c + 2.*seczs[2];

    est = TMath::Abs(dxt)+TMath::Abs(dyt)+TMath::Abs(dzt);
    if (est > 2.*TMath::Abs(h)) {
      if (ncut++ > maxcut) break;
      h *= khalf;
      continue;
    }

    xyzt[0] = xt;
    xyzt[1] = yt;
    xyzt[2] = zt;

    ///cmodif: call gufld(xyzt,f) changed into:
    TGeoGlobalMagField::Instance()->Field(xyzt,f);

    z      = z + (c + (seczs[0] + seczs[1] + seczs[2]) * kthird) * h;
    y      = y + (b + (secys[0] + secys[1] + secys[2]) * kthird) * h;
    x      = x + (a + (secxs[0] + secxs[1] + secxs[2]) * kthird) * h;

    secxs[3] = (bt*f[2] - ct*f[1])* ph2;
    secys[3] = (ct*f[0] - at*f[2])* ph2;
    seczs[3] = (at*f[1] - bt*f[0])* ph2;
    a      = a+(secxs[0]+secxs[3]+2. * (secxs[1]+secxs[2])) * kthird;
    b      = b+(secys[0]+secys[3]+2. * (secys[1]+secys[2])) * kthird;
    c      = c+(seczs[0]+seczs[3]+2. * (seczs[1]+seczs[2])) * kthird;

    est    = TMath::Abs(secxs[0]+secxs[3] - (secxs[1]+secxs[2]))
      + TMath::Abs(secys[0]+secys[3] - (secys[1]+secys[2]))
      + TMath::Abs(seczs[0]+seczs[3] - (seczs[1]+seczs[2]));

    if (est > kdlt && TMath::Abs(h) > 1.e-4) {
      if (ncut++ > maxcut) break;
      h *= khalf;
      continue;
    }

    ncut = 0;
    /// *               if too many iterations, go to helix
    if (iter++ > maxit) break;

    tl += h;
    if (est < kdlt32)
      h *= 2.;
    cba    = 1./ TMath::Sqrt(a*a + b*b + c*c);
    vout[0] = x;
    vout[1] = y;
    vout[2] = z;
    vout[3] = cba*a;
    vout[4] = cba*b;
    vout[5] = cba*c;
    rest = step - tl;
    if (step < 0.) rest = -rest;
    if (rest < 1.e-5*TMath::Abs(step)) return;

  } while(1);

  /// angle too big, use helix

  f1  = f[0];
  f2  = f[1];
  f3  = f[2];
  f4  = TMath::Sqrt(f1*f1+f2*f2+f3*f3);
  rho = -f4*pinv;
  tet = rho * step;

  hnorm = 1./f4;
  f1 = f1*hnorm;
  f2 = f2*hnorm;
  f3 = f3*hnorm;

  hxp[0] = f2*vect[kipz] - f3*vect[kipy];
  hxp[1] = f3*vect[kipx] - f1*vect[kipz];
  hxp[2] = f1*vect[kipy] - f2*vect[kipx];

  hp = f1*vect[kipx] + f2*vect[kipy] + f3*vect[kipz];

  rho1 = 1./rho;
  sint = TMath::Sin(tet);
  cost = 2.*TMath::Sin(khalf*tet)*TMath::Sin(khalf*tet);

  g1 = sint*rho1;
  g2 = cost*rho1;
  g3 = (tet-sint) * hp*rho1;
  g4 = -cost;
  g5 = sint;
  g6 = cost * hp;

  vout[kix] = vect[kix] + g1*vect[kipx] + g2*hxp[0] + g3*f1;
  vout[kiy] = vect[kiy] + g1*vect[kipy] + g2*hxp[1] + g3*f2;
  vout[kiz] = vect[kiz] + g1*vect[kipz] + g2*hxp[2] + g3*f3;

  vout[kipx] = vect[kipx] + g4*vect[kipx] + g5*hxp[0] + g6*f1;
  vout[kipy] = vect[kipy] + g4*vect[kipy] + g5*hxp[1] + g6*f2;
  vout[kipz] = vect[kipz] + g4*vect[kipz] + g5*hxp[2] + g6*f3;

  return;
}

///______________________________________________________________________________


Bool_t AliHLTMUONFullTracker::SelectFront()
{
  /// Select the front trackseg connected with back trackseg

  Cluster clus1,clus2;
  Int_t minIndex=0,maxIndex=0;
  Int_t minCh=0,maxCh=0;
  Int_t ifronttrackseg=0;
  
  Float_t xf,yf,zf,thetaDev;
  Float_t p,spx,spy,spz,px,py,pz,sx,sy,sz,x,y,z;
  Double_t charge;
  Float_t dist,mindist;
  Int_t frontsegIndex ;

  for( Int_t ibacktrackseg=0;ibacktrackseg<fNofbackTrackSeg;ibacktrackseg++){

    if(fNofConnectedfrontTrackSeg[ibacktrackseg]<=0) continue;
    
    ///     if(fBackTrackSeg[ibacktrackseg].fIndex[2]==-1 || fBackTrackSeg[ibacktrackseg].fIndex[3]==-1) continue;

    ifronttrackseg = fBackToFront[ibacktrackseg][fNofConnectedfrontTrackSeg[ibacktrackseg]-1];
    
#ifdef PRINT_SELECT
    printf("AliHLTMUONFullTracker::SelectFront()--Begin\n\n");
    printf("\tbacktrack : %d is connected with : %d front tracks\n",
	   ibacktrackseg,fNofConnectedfrontTrackSeg[ibacktrackseg]);
#endif

    maxIndex = (fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1)?3:2;
    maxCh = (fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1)?9:8;

    minIndex = (fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1)?0:1;
    minCh = (fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1)?6:7;


    
    clus2.fX = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fX ;
    clus2.fY = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fY ;
    clus2.fZ = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fZ ;
    
    clus1.fX = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fX ;
    clus1.fY = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fY ;
    clus1.fZ = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fZ ;

    zf = 0.5*(AliMUONConstants::DefaultChamberZ(4) + AliMUONConstants::DefaultChamberZ(5)); 
    thetaDev= (1/zf)*(clus1.fY*clus2.fZ - clus2.fY*clus1.fZ)/(clus2.fZ - clus1.fZ);
    xf = clus1.fX*zf/clus1.fZ;
    yf = clus2.fY - ((clus2.fY - clus1.fY)*(clus2.fZ-zf))/(clus2.fZ - clus1.fZ);
    p = 3.0*0.3/thetaDev;
    charge = (p>0)?-1:+1;

    spx = p*xf/zf;
    spy = p*yf/zf;
    spz = sqrt((p*p-(spx*spx + spy*spy)));

    if(zf<0)spz = -spz;
    
    sx = clus1.fX; sy = clus1.fY; sz = clus1.fZ;
    mindist = 200000.0;
    frontsegIndex = -1;
    for( Int_t iconnected=0;iconnected<fNofConnectedfrontTrackSeg[ibacktrackseg];iconnected++){
      
      ifronttrackseg = fBackToFront[ibacktrackseg][iconnected];
      
      minIndex = (fFrontTrackSeg[ifronttrackseg].fIndex[3]!=-1)?3:2;
      minCh = (fFrontTrackSeg[ifronttrackseg].fIndex[3]!=-1)?3:2;
      clus1.fX = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fX ;
      clus1.fY = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fY ;
      clus1.fZ = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fZ ;
      
      x = sx; y = sy; z = sz;
      px = spx; py = spy; pz = spz;
      PropagateTracks(charge,px,py,pz,x,y,z,clus1.fZ);
      
      dist = sqrt((clus1.fX-x)*(clus1.fX-x) + 
		  (clus1.fY-y)*(clus1.fY-y));
      if(dist<mindist){
	mindist = dist;
	frontsegIndex = ifronttrackseg;
      }
    }///for loop on all connected front track segs
    
    fNofConnectedfrontTrackSeg[ibacktrackseg] = 0;
    ///have to check later
    if(frontsegIndex==-1) continue;

    fBackToFront[ibacktrackseg][fNofConnectedfrontTrackSeg[ibacktrackseg]++] = frontsegIndex; 
    
    ///     fTrackParam[ibacktrackseg] = trackParam;
    
#ifdef PRINT_SELECT
    printf("\tbacktrack : %d is connected with : %d front tracks\n",
	   ibacktrackseg,fNofConnectedfrontTrackSeg[ibacktrackseg]);
    printf("AliHLTMUONFullTracker::SelectFront()--End\n\n");

#endif
    
  }///backtrackSeg loop
  
  
  return true;
}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::KalmanChi2Test()
{
  
  /// Kalman Chi2 test for trach segments selection

  Cluster clus1,clus2;
  Int_t minIndex=0,maxIndex=0;
  Int_t minCh=0,maxCh=0;
  Int_t ifronttrackseg=0;
  for( Int_t ibacktrackseg=0;ibacktrackseg<fNofbackTrackSeg;ibacktrackseg++){

    if(fNofConnectedfrontTrackSeg[ibacktrackseg]<=0) continue;
    
    ///     if(fBackTrackSeg[ibacktrackseg].fIndex[2]==-1 || fBackTrackSeg[ibacktrackseg].fIndex[3]==-1) continue;

    ifronttrackseg = fBackToFront[ibacktrackseg][fNofConnectedfrontTrackSeg[ibacktrackseg]-1];
    
#ifdef PRINT_KALMAN
    printf("AliHLTMUONFullTracker::KalmanChi2Test()--Begin\n\n");
    printf("\tbacktrack : %d is connected with : %d front tracks, front track index %d\n",
	   ibacktrackseg,fNofConnectedfrontTrackSeg[ibacktrackseg],ifronttrackseg);
#endif
    maxIndex = (fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1)?3:2;
    maxCh = (fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1)?9:8;

    minIndex = (fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1)?0:1;
    minCh = (fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1)?6:7;
      
    
    clus2.fX = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fX ;
    clus2.fY = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fY ;
    clus2.fZ = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fZ ;
    clus2.fErrX2 =  0.020736;
    clus2.fErrY2 =  0.000100;
    
    clus1.fX = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fX ;
    clus1.fY = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fY ;
    clus1.fZ = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fZ ;
    clus1.fErrX2 =  0.020736;
    clus1.fErrY2 =  0.000100;


    AliMUONTrackParam trackParam;
    Double_t dZ =  clus2.fZ - clus1.fZ;
    trackParam.SetNonBendingCoor(clus2.fX);
    trackParam.SetBendingCoor(clus2.fY);
    trackParam.SetZ(clus2.fZ);
    trackParam.SetNonBendingSlope((clus2.fX - clus1.fX) / dZ);
    trackParam.SetBendingSlope((clus2.fY - clus1.fY) / dZ);
    Double_t bendingImpact = clus1.fY - clus1.fZ * trackParam.GetBendingSlope();
    Double_t inverseBendingMomentum = 1. 
      / AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(bendingImpact);
    trackParam.SetInverseBendingMomentum(inverseBendingMomentum);
    
#ifdef PRINT_KALMAN
    printf("\t\tCh%d  : (%f,%f,%f)\n",maxCh,clus2.fX,clus2.fY,clus2.fZ);
    printf("\t\tCh%d  : (%f,%f,%f)\n",minCh,clus1.fX,clus1.fY,clus1.fZ);
    printf("\t\tFor minCh : %d, GetBeMom : %lf\n",minCh,trackParam.GetInverseBendingMomentum());
#endif      

    ///       trackParam->SetClusterPtr(clus[8]);
      
    /// => Reset track parameter covariances at last clus (as if the other cluss did not exist)
    TMatrixD lastParamCov(5,5);
    lastParamCov.Zero();
    lastParamCov(0,0) = clus1.fErrX2;
    lastParamCov(1,1) = 100. * ( clus2.fErrX2 + clus1.fErrX2 ) / dZ / dZ;
    lastParamCov(2,2) = clus1.fErrY2;
    lastParamCov(3,3) = 100. * ( clus2.fErrY2 + clus1.fErrY2 ) / dZ / dZ;
    lastParamCov(4,4) = ((10.0*10.0 + (clus2.fZ * clus2.fZ * clus1.fErrY2 + 
				       clus1.fZ * clus1.fZ * clus2.fErrY2) 
			  / dZ / dZ) /bendingImpact / bendingImpact +
			 0.1 * 0.1) * inverseBendingMomentum * inverseBendingMomentum;
    trackParam.SetCovariances(lastParamCov);

    AliMUONTrackExtrap::AddMCSEffect(&trackParam,AliMUONConstants::ChamberThicknessInX0(minCh),1.);
    ///     AliMUONTrackExtrap::ExtrapToZCov(trackParam, AliMUONConstants::DefaultChamberZ(minCh),kTRUE);
    ///AliMUONTrackExtrap::ExtrapToZCov(&trackParam, clus1.fZ,kTRUE);

    LinearExtrapToZ(&trackParam, clus1.fZ);

    trackParam.SetTrackChi2(0.);
    
    Double_t chi2 = 0.0;

    chi2 = KalmanFilter(trackParam,&clus1);
    
    if( chi2 > 1.0e9 /* is order to check TMath::AreEqualAbs(chi2,,1.0e-5)*/  ) {
      HLTWarning("Kalman Chi2 calculation cannot be completed...skipping slat track %d (c.p = 1)",ibacktrackseg);
      fNofConnectedfrontTrackSeg[ibacktrackseg] = 0;
      continue;
    }

#ifdef PRINT_KALMAN
    printf("\t\tFor minCh : %d, Chi2 = %lf, GetBeMom : %lf\n",minCh,chi2,trackParam.GetInverseBendingMomentum());
    ///       TMatrixD para2(trackParam->GetParameters());
    ///       para2.Print();
#endif      

    AliMUONTrackParam extrapTrackParamAtCluster2,minChi2Param;
    Double_t chi2OfCluster;
    Bool_t tryonefast;
    Double_t minChi2=1000000.0;
    Int_t frontsegIndex = -1;
    extrapTrackParamAtCluster2 = trackParam ;
    minChi2Param = trackParam ;

    for( Int_t iconnected=0;iconnected<fNofConnectedfrontTrackSeg[ibacktrackseg];iconnected++){
      
      ifronttrackseg = fBackToFront[ibacktrackseg][iconnected];
      AliMUONTrackParam extrapTrackParam;
      trackParam = extrapTrackParamAtCluster2;
      
      minIndex = (fFrontTrackSeg[ifronttrackseg].fIndex[3]!=-1)?3:2;
      minCh = (fFrontTrackSeg[ifronttrackseg].fIndex[3]!=-1)?3:2;
      clus1.fX = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fX ;
      clus1.fY = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fY ;
      clus1.fZ = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fZ ;
      clus1.fErrX2 =  0.020736;
      clus1.fErrY2 =  0.000100;

    
      AliMUONTrackExtrap::AddMCSEffect(&trackParam,AliMUONConstants::ChamberThicknessInX0(minCh),1.);
      trackParam.ResetPropagator();
      AliMUONTrackExtrap::ExtrapToZCov(&trackParam,clus1.fZ, kTRUE);
      
      tryonefast = TryOneClusterFast(trackParam, &clus1);
      ///if (!tryonefast) continue;
      
      /// try to add the current cluster accuratly
      chi2OfCluster = TryOneCluster(trackParam, &clus1 , extrapTrackParam,kTRUE);
      
      extrapTrackParam.SetExtrapParameters(extrapTrackParam.GetParameters());
      extrapTrackParam.SetExtrapCovariances(extrapTrackParam.GetCovariances());
      
      chi2 = KalmanFilter(extrapTrackParam,&clus1);
      if( chi2 > 1.0e9 /* is order to check TMath::AreEqualAbs(chi2,,1.0e-5)*/  ) {
	HLTWarning("Kalman Chi2 calculation cannot be completed...skipping slat track %d (c.p = 2)",ibacktrackseg);
	fNofConnectedfrontTrackSeg[ibacktrackseg] = 0;
	continue;
      }
      if(chi2<minChi2){
	minChi2 = chi2;
	minChi2Param = extrapTrackParam;
	frontsegIndex = ifronttrackseg;
      }
    }///for loop on all connected front track segs
    
    fNofConnectedfrontTrackSeg[ibacktrackseg] = 0;
    ///have to check later
    if(frontsegIndex==-1) continue;

    fBackToFront[ibacktrackseg][fNofConnectedfrontTrackSeg[ibacktrackseg]++] = frontsegIndex; 
    ifronttrackseg = frontsegIndex;
    
    trackParam = minChi2Param ;
    
#ifdef PRINT_KALMAN
    printf("\t\t\tCh%d  : (%f,%f,%f)\n",minCh,clus1.fX,clus1.fY,clus1.fZ);
    printf("\t\t\tFor minCh : %d, Chi2 = %lf, GetBeMom : %lf\t",minCh,chi2,trackParam.GetInverseBendingMomentum());
    printf(Form("Pt : %lf\n",TMath::Sqrt(trackParam.Px()*trackParam.Px() + 
					 trackParam.Py()*trackParam.Py())));
#endif      
    
    
    minIndex = (fFrontTrackSeg[ifronttrackseg].fIndex[0]!=-1)?0:1;
    minCh = (fFrontTrackSeg[ifronttrackseg].fIndex[0]!=-1)?0:1;
    clus1.fX = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fX ;
    clus1.fY = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fY ;
    clus1.fZ = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fZ ;
    clus1.fErrX2 =  0.020736;
    clus1.fErrY2 =  0.000100;
    
    AliMUONTrackExtrap::AddMCSEffect(&trackParam,AliMUONConstants::ChamberThicknessInX0(minCh),1.);
    trackParam.ResetPropagator();
    ///AliMUONTrackExtrap::ExtrapToZCov(&trackParam,clus1.fZ, kTRUE);
    LinearExtrapToZ(&trackParam, clus1.fZ);

    tryonefast = TryOneClusterFast(trackParam, &clus1);
    ///if (!tryonefast) continue;
	  
    chi2OfCluster = TryOneCluster(trackParam, &clus1 , extrapTrackParamAtCluster2,kTRUE);
	  
    extrapTrackParamAtCluster2.SetExtrapParameters(extrapTrackParamAtCluster2.GetParameters());
    extrapTrackParamAtCluster2.SetExtrapCovariances(extrapTrackParamAtCluster2.GetCovariances());
    
    chi2 = KalmanFilter(extrapTrackParamAtCluster2,&clus1);
    if(chi2 > 1.0e9 /* is order to check TMath::AreEqualAbs(chi2,,1.0e-5)*/ ) {
      HLTWarning("Kalman Chi2 calculation cannot be completed...skipping slat track %d (c.p = 3)",ibacktrackseg);
      fNofConnectedfrontTrackSeg[ibacktrackseg] = 0;
      continue;
    }

    trackParam = extrapTrackParamAtCluster2;
    
#ifdef PRINT_KALMAN
    printf("\t\tCh%d  : (%f,%f,%f)\n",minCh,clus1.fX,clus1.fY,clus1.fZ);
    printf("\t\tFor minCh : %d, Chi2 = %lf, GetBeMom : %lf\t",minCh,chi2,trackParam.GetInverseBendingMomentum());
    printf(Form("Pt : %lf\n",TMath::Sqrt(extrapTrackParamAtCluster2.Px()*extrapTrackParamAtCluster2.Px() + 
					 extrapTrackParamAtCluster2.Py()*extrapTrackParamAtCluster2.Py())));
    trackParam.Print();
    printf("AliHLTMUONFullTracker::KalmanChi2Test()--End\n\n");
#endif      
    ///AliMUONTrackExtrap::ExtrapToVertex(&trackParam, 0., 0., 0., 0., 0.);
    //trackParam.SetInverseBendingMomentum(trackParam.GetCharge());
    fTrackParam[ibacktrackseg] = trackParam;
    

    
  }///trig loop
  

  
  return true;
}

///__________________________________________________________________________


Double_t AliHLTMUONFullTracker::BetheBloch(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicA, Double_t atomicZ)
{
  //// Returns the mean total momentum energy loss of muon with total momentum='pTotal'
  //// in the absorber layer of lenght='pathLength', density='rho', A='atomicA' and Z='atomicZ'
  Double_t muMass = 0.105658369; /// GeV
  Double_t eMass = 0.510998918e-3; /// GeV
  Double_t k = 0.307075e-3; /// GeV.g^-1.cm^2
  Double_t i = 9.5e-9; /// mean exitation energy per atomic Z (GeV)
  Double_t p2=pTotal*pTotal;
  Double_t beta2=p2/(p2 + muMass*muMass);
  
  Double_t w = k * rho * pathLength * atomicZ / atomicA / beta2;
  
  if (beta2/(1-beta2)>3.5*3.5)
    return w * (log(2.*eMass*3.5/(i*atomicZ)) + 0.5*log(beta2/(1-beta2)) - beta2);
  
  return w * (log(2.*eMass*beta2/(1-beta2)/(i*atomicZ)) - beta2);
}


///__________________________________________________________________________

Double_t AliHLTMUONFullTracker::EnergyLossFluctuation2(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicA, Double_t atomicZ)
{
  //// Returns the total momentum energy loss fluctuation of muon with total momentum='pTotal'
  //// in the absorber layer of lenght='pathLength', density='rho', A='atomicA' and Z='atomicZ'
  Double_t muMass = 0.105658369; /// GeV
  ///Double_t eMass = 0.510998918e-3; /// GeV
  Double_t k = 0.307075e-3; /// GeV.g^-1.cm^2
  Double_t p2=pTotal*pTotal;
  Double_t beta2=p2/(p2 + muMass*muMass);
  
  Double_t fwhm = 2. * k * rho * pathLength * atomicZ / atomicA / beta2; /// FWHM of the energy loss Landau distribution
  Double_t sigma2 = fwhm * fwhm / (8.*log(2.)); /// gaussian: fwmh = 2 * srqt(2*ln(2)) * sigma (i.e. fwmh = 2.35 * sigma)
  
  ///sigma2 = k * rho * pathLength * atomicZ / atomicA * eMass; /// sigma2 of the energy loss gaussian distribution
  
  return sigma2;
}


///__________________________________________________________________________

void AliHLTMUONFullTracker::LinearExtrapToZ(AliMUONTrackParam* trackParam, Double_t zEnd)
{
  //// Track parameters (and their covariances if any) linearly extrapolated to the plane at "zEnd".
  //// On return, results from the extrapolation are updated in trackParam.
  
  if ( TMath::AreEqualAbs(trackParam->GetZ(),zEnd,1.0e-5) ) return; /// nothing to be done if same z

  Double_t dZ = zEnd - trackParam->GetZ();
  trackParam->SetNonBendingCoor(trackParam->GetNonBendingCoor() + trackParam->GetNonBendingSlope() * dZ);
  trackParam->SetBendingCoor(trackParam->GetBendingCoor() + trackParam->GetBendingSlope() * dZ);
  trackParam->SetZ(zEnd);
}

void AliHLTMUONFullTracker::CorrectELossEffectInAbsorber(AliMUONTrackParam* param, Double_t eLoss)
{
  /// Energy loss correction
  Double_t nonBendingSlope = param->GetNonBendingSlope();
  Double_t bendingSlope = param->GetBendingSlope();
  param->SetInverseBendingMomentum(param->GetCharge() / (param->P() + eLoss) *
				   TMath::Sqrt(1.0 + nonBendingSlope*nonBendingSlope + bendingSlope*bendingSlope) /
				   TMath::Sqrt(1.0 + bendingSlope*bendingSlope));
}

///__________________________________________________________________________

void AliHLTMUONFullTracker::Cov2CovP(const TMatrixD &param, TMatrixD &cov)
{
  //// change coordinate system: (X, SlopeX, Y, SlopeY, q/Pyz) -> (X, SlopeX, Y, SlopeY, q*PTot)
  //// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
  
  /// charge * total momentum
  Double_t qPTot = TMath::Sqrt(1. + param(1,0)*param(1,0) + param(3,0)*param(3,0)) /
    TMath::Sqrt(1. + param(3,0)*param(3,0)) / param(4,0);
  
  /// Jacobian of the opposite transformation
  TMatrixD jacob(5,5);
  jacob.UnitMatrix();
  jacob(4,1) = qPTot * param(1,0) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,3) = - qPTot * param(1,0) * param(1,0) * param(3,0) /
    (1. + param(3,0)*param(3,0)) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,4) = - qPTot / param(4,0);
  
  /// compute covariances in new coordinate system
  TMatrixD tmp(5,5);
  tmp.MultT(cov,jacob);
  cov.Mult(jacob,tmp);

}

///__________________________________________________________________________

void AliHLTMUONFullTracker::CovP2Cov(const TMatrixD &param, TMatrixD &covP)
{
  //// change coordinate system: (X, SlopeX, Y, SlopeY, q*PTot) -> (X, SlopeX, Y, SlopeY, q/Pyz)
  //// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
  
  /// charge * total momentum
  Double_t qPTot = TMath::Sqrt(1. + param(1,0)*param(1,0) + param(3,0)*param(3,0)) /
    TMath::Sqrt(1. + param(3,0)*param(3,0)) / param(4,0);
  
  /// Jacobian of the transformation
  TMatrixD jacob(5,5);
  jacob.UnitMatrix();
  jacob(4,1) = param(4,0) * param(1,0) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,3) = - param(4,0) * param(1,0) * param(1,0) * param(3,0) /
    (1. + param(3,0)*param(3,0)) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,4) = - param(4,0) / qPTot;
  
  /// compute covariances in new coordinate system
  TMatrixD tmp(5,5);
  tmp.MultT(covP,jacob);
  covP.Mult(jacob,tmp);

}

///__________________________________________________________________________


void AliHLTMUONFullTracker::CorrectMCSEffectInAbsorber(AliMUONTrackParam* param,
						       Double_t xVtx, Double_t yVtx, Double_t zVtx,
						       Double_t absZBeg, Double_t f1, Double_t f2)
{
  //// Correct parameters and corresponding covariances using Branson correction
  //// - input param are parameters and covariances at the end of absorber
  //// - output param are parameters and covariances at vertex
  //// Absorber correction parameters are supposed to be calculated at the current track z-position
  
  /// Position of the Branson plane (spectro. (z<0))
  Double_t zB = (f1>0.) ? absZBeg - f2/f1 : 0.;
  
  LinearExtrapToZ(param,zB);
  
  /// compute track parameters at vertex
  TMatrixD newParam(5,1);
  newParam.Zero();
  newParam(0,0) = xVtx;
  newParam(1,0) = (param->GetNonBendingCoor() - xVtx) / (zB - zVtx);
  newParam(2,0) = yVtx;
  newParam(3,0) = (param->GetBendingCoor() - yVtx) / (zB - zVtx);
  newParam(4,0) = param->GetCharge() / param->P() *
    TMath::Sqrt(1.0 + newParam(1,0)*newParam(1,0) + newParam(3,0)*newParam(3,0)) /
    TMath::Sqrt(1.0 + newParam(3,0)*newParam(3,0));
  
  TMatrixD paramCovP(5,5);
  
  /// Get covariances in (X, SlopeX, Y, SlopeY, q*PTot) coordinate system
  paramCovP.Use(param->GetCovariances());
  
  Cov2CovP(param->GetParameters(),paramCovP);
  
  /// Get the covariance matrix in the (XVtx, X, YVtx, Y, q*PTot) coordinate system
  ///   TMatrixD paramCovVtx(5,5);
  Double_t element44 = paramCovP(4,4);
  paramCovP.Zero();
  paramCovP(4,4) = element44;
  
  CovP2Cov(newParam,paramCovP);
  
  /// Set parameters and covariances at vertex
  param->SetParameters(newParam);
  param->SetZ(zVtx);
  param->SetCovariances(paramCovP);
}

///__________________________________________________________________________

Bool_t AliHLTMUONFullTracker::ExtrapolateToOrigin()
{
  /// Extrapolation to origin through absorber

  Int_t minIndex=0,maxIndex=0;
  Int_t minCh=0,maxCh=0;
  Int_t ifronttrackseg = -1;
  AliMUONTrackParam trackP;
  Double_t bSlope, nbSlope;
  AliHLTMUONRecHitStruct p1,p2,pSeg1,pSeg2;
  Double_t pyz = -1.0;
  TVector3 v1,v2,v3,v4;
  Double_t eLoss1,eLoss2,eLoss3;
  Double_t b;
  Double_t zE,zB,dzE,dzB; 
  Double_t f0,f1,f2;
  Double_t f0Sum,f1Sum,f2Sum;
  Double_t fXVertex=0.0,fYVertex=0.0,fZVertex=0.0;
  Double_t thetaDev;

  for( Int_t ibacktrackseg=0;ibacktrackseg<fNofbackTrackSeg;ibacktrackseg++){
    
    if(fNofConnectedfrontTrackSeg[ibacktrackseg]<=0) continue;
    
    maxIndex = (fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1)?3:2;
    maxCh = (fBackTrackSeg[ibacktrackseg].fIndex[3]!=-1)?9:8;
    
    minIndex = (fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1)?0:1;
    minCh = (fBackTrackSeg[ibacktrackseg].fIndex[0]!=-1)?6:7;
    
    p2.fX = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fX ; 
    p2.fY = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fY ;
    p2.fZ = fChPoint[maxCh][fBackTrackSeg[ibacktrackseg].fIndex[maxIndex]]->fZ ;

    p1.fX = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fX ; 
    p1.fY = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fY ;
    p1.fZ = fChPoint[minCh][fBackTrackSeg[ibacktrackseg].fIndex[minIndex]]->fZ ;
    
    thetaDev= (p1.fY*p2.fZ - p2.fY*p1.fZ)/(p2.fZ - p1.fZ);

    Sub(&p2,&p1,&pSeg1);
    
    ifronttrackseg = fBackToFront[ibacktrackseg][0];

    maxIndex = (fFrontTrackSeg[ifronttrackseg].fIndex[3]!=-1)?3:2;
    maxCh = (fFrontTrackSeg[ifronttrackseg].fIndex[3]!=-1)?3:2;
    
    minIndex = (fFrontTrackSeg[ifronttrackseg].fIndex[0]!=-1)?0:1;
    minCh = (fFrontTrackSeg[ifronttrackseg].fIndex[0]!=-1)?0:1;

    p2.fX = fChPoint[maxCh][fFrontTrackSeg[ifronttrackseg].fIndex[maxIndex]]->fX ;
    p2.fY = fChPoint[maxCh][fFrontTrackSeg[ifronttrackseg].fIndex[maxIndex]]->fY ;
    p2.fZ = fChPoint[maxCh][fFrontTrackSeg[ifronttrackseg].fIndex[maxIndex]]->fZ ;

    p1.fX = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fX ;
    p1.fY = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fY ;
    p1.fZ = fChPoint[minCh][fFrontTrackSeg[ifronttrackseg].fIndex[minIndex]]->fZ ;

    Sub(&p2,&p1,&pSeg2);

    if(thetaDev > 0)
      pyz = (3.0*0.3/sin(Angle(&pSeg1,&pSeg2)));/// *  sqrt(x3*x3 + y3*y3)/z3 ;
    else
      pyz = -(3.0*0.3/sin(Angle(&pSeg1,&pSeg2)));/// *  sqrt(x3*x3 + y3*y3)/z3 ;
    
    nbSlope = (p2.fX - p1.fX)/(p2.fZ - p1.fZ); 
    bSlope = (p2.fY - p1.fY)/(p2.fZ - p1.fZ); 
    
    trackP.SetZ(p1.fZ);
    trackP.SetNonBendingCoor(p1.fX);
    trackP.SetNonBendingSlope(nbSlope);
    trackP.SetBendingCoor(p1.fY);
    trackP.SetBendingSlope(bSlope);
    trackP.SetInverseBendingMomentum(1.0/pyz) ;
    
    


    if(not fFastTracking){
      trackP = fTrackParam[ibacktrackseg]    ;
    }
    
    LinearExtrapToZ(&trackP,fgkAbsoedge[3]);
    v4.SetXYZ(trackP.GetNonBendingCoor(),trackP.GetBendingCoor(),trackP.GetZ());
    LinearExtrapToZ(&trackP,fgkAbsoedge[2]);
    v3.SetXYZ(trackP.GetNonBendingCoor(),trackP.GetBendingCoor(),trackP.GetZ());
    LinearExtrapToZ(&trackP,fgkAbsoedge[1]);
    v2.SetXYZ(trackP.GetNonBendingCoor(),trackP.GetBendingCoor(),trackP.GetZ());
    LinearExtrapToZ(&trackP,fgkAbsoedge[0]);
    v1.SetXYZ(trackP.GetNonBendingCoor(),trackP.GetBendingCoor(),trackP.GetZ());
    
    eLoss1 = BetheBloch(trackP.P(), (v4-v3).Mag(), fgkRho[2], fgkAtomicA[2], fgkAtomicZ[2]);
    eLoss2 = BetheBloch(trackP.P(), (v3-v2).Mag(), fgkRho[1], fgkAtomicA[1], fgkAtomicZ[1]);
    eLoss3 = BetheBloch(trackP.P(), (v2-v1).Mag(), fgkRho[0], fgkAtomicA[0], fgkAtomicZ[0]);
    
    ///       sigmaELoss1 = EnergyLossFluctuation2(trackP.P(), (v4-v3).Mag(), rho[2], atomicA[2], atomicZ[2]);
    ///       sigmaELoss2 = EnergyLossFluctuation2(trackP.P(), (v3-v2).Mag(), rho[1], atomicA[1], atomicZ[1]);
    ///       sigmaELoss3 = EnergyLossFluctuation2(trackP.P(), (v2-v1).Mag(), rho[0], atomicA[0], atomicZ[0]);

    ///     eDiff = totELoss-(eLoss1+eLoss2+eLoss3);
    ///     sigmaELossDiff = sigmaTotELoss ;///- (sigmaELoss1+sigmaELoss2+sigmaELoss3);

    ///       CorrectELossEffectInAbsorber(&trackP, 0.5*(eLoss1+eLoss2+eLoss3), 0.5*(sigmaELoss1+sigmaELoss2+sigmaELoss3));
    

    ///CorrectELossEffectInAbsorber(&trackP, totELoss,sigmaTotELoss);

    CorrectELossEffectInAbsorber(&trackP, 0.7*(eLoss1+eLoss2+eLoss3));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    f0Sum = 0.0;      f1Sum = 0.0;      f2Sum = 0.0;

    b = (v4.Z()-v1.Z())/((v4-v1).Mag());
    
    zB = v1.Z();
    zE = b*((v2-v1).Mag()) + zB;
    dzB = zB - v1.Z();
    dzE = zE - v1.Z();
    
    f0 = ((v2-v1).Mag())/fgkRadLen[0];
    f1  = (dzE*dzE - dzB*dzB) / b / b / fgkRadLen[0] /2.;
    f2 = (dzE*dzE*dzE - dzB*dzB*dzB) / b / b / b / fgkRadLen[0] / 3.;

    f0Sum += f0;
    f1Sum += f1;
    f2Sum += f2;
    
    zB = zE;
    zE = b*((v3-v2).Mag()) + zB;
    dzB = zB - v1.Z();
    dzE = zE - v1.Z();
    
    f0 = ((v3-v2).Mag())/fgkRadLen[1];
    f1  = (dzE*dzE - dzB*dzB) / b / b / fgkRadLen[1] /2.;
    f2 = (dzE*dzE*dzE - dzB*dzB*dzB) / b / b / b / fgkRadLen[1] / 3.;

    f0Sum += f0;
    f1Sum += f1;
    f2Sum += f2;
    
    zB = zE;
    zE = b*((v4-v3).Mag()) + zB;
    dzB = zB - v1.Z();
    dzE = zE - v1.Z();
    
    f0 = ((v4-v3).Mag())/fgkRadLen[2];
    f1  = (dzE*dzE - dzB*dzB) / b / b / fgkRadLen[2] /2.;
    f2 = (dzE*dzE*dzE - dzB*dzB*dzB) / b / b / b / fgkRadLen[2] / 3.;

    f0Sum += f0;
    f1Sum += f1;
    f2Sum += f2;
    
    
    ///AddMCSEffectInAbsorber(&trackP,(v4-v1).Mag(),f0Sum,f1Sum,f2Sum);

    CorrectMCSEffectInAbsorber(&trackP,fXVertex,fYVertex, fZVertex,AliMUONConstants::AbsZBeg(),f1Sum,f2Sum);
    CorrectELossEffectInAbsorber(&trackP, 0.5*(eLoss1+eLoss2+eLoss3));


    ///AliMUONTrackExtrap::ExtrapToVertex(&trackP, 0., 0., 0., 0., 0.);
    trackP.SetZ(p1.fZ);
    trackP.SetNonBendingCoor(p1.fX);
    trackP.SetBendingCoor(p1.fY);
    LinearExtrapToZ(&trackP, 0.0);

    fTrackParam[ibacktrackseg] = trackP;
    
  }///backtrackseg for loop
  
  return true;
}

