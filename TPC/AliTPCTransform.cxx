/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-------------------------------------------------------
//          Implementation of the TPC transformation class
//
//   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
//           Magnus Mager
//
//   Class for tranformation of the coordinate frame
//   Transformation  
//    local coordinate frame (sector, padrow, pad, timebine) ==>
//    rotated global (tracking) cooridnate frame (sector, lx,ly,lz)
//
//    Unisochronity  - (substract time0 - pad by pad)
//    Drift velocity - Currently common drift velocity - functionality of AliTPCParam
//    ExB effect     - 
//
//    Time of flight correction -
//                   - Depends on the vertex position
//                   - by default 
//                           
//    Usage:
//          AliTPCclustererMI::AddCluster
//          AliTPCtrackerMI::Transform
//    
//-------------------------------------------------------

/* To test it:
   cdb=AliCDBManager::Instance()
   cdb->SetDefaultStorage("local:///u/mmager/mycalib1")
   c=AliTPCcalibDB::Instance()
   c->SetRun(0)
   Double_t x[]={1.0,2.0,3.0}
   Int_t i[]={4}
   AliTPCTransform trafo
   trafo.Transform(x,i,0,1)
 */

/* $Id$ */

#include "AliTPCROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "TMath.h"
#include "AliLog.h"
#include "AliTPCExB.h"
#include "AliTPCCorrection.h"
#include "TGeoMatrix.h"
#include "AliTPCRecoParam.h"
#include "AliTPCCalibVdrift.h"
#include "AliTPCTransform.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTracker.h"
#include <AliCTPTimeParams.h>

ClassImp(AliTPCTransform)


AliTPCTransform::AliTPCTransform():
  AliTransform(),
  fCurrentRecoParam(0),       //! current reconstruction parameters
  fCurrentRun(0),             //! current run
  fCurrentTimeStamp(0)        //! current time stamp   
{
  //
  // Speed it up a bit!
  //
  for (Int_t i=0;i<18;++i) {
    Double_t alpha=TMath::DegToRad()*(10.+20.*(i%18));
    fSins[i]=TMath::Sin(alpha);
    fCoss[i]=TMath::Cos(alpha);
  }
  fPrimVtx[0]=0;
  fPrimVtx[1]=0;
  fPrimVtx[2]=0;
}
AliTPCTransform::AliTPCTransform(const AliTPCTransform& transform):
  AliTransform(transform),
  fCurrentRecoParam(transform.fCurrentRecoParam),       //! current reconstruction parameters
  fCurrentRun(transform.fCurrentRun),             //! current run
  fCurrentTimeStamp(transform.fCurrentTimeStamp)        //! current time stamp   
{
  //
  // Speed it up a bit!
  //
  for (Int_t i=0;i<18;++i) {
    Double_t alpha=TMath::DegToRad()*(10.+20.*(i%18));
    fSins[i]=TMath::Sin(alpha);
    fCoss[i]=TMath::Cos(alpha);
  }
  fPrimVtx[0]=0;
  fPrimVtx[1]=0;
  fPrimVtx[2]=0;
}

AliTPCTransform::~AliTPCTransform() {
  //
  // Destructor
  //
}

void AliTPCTransform::SetPrimVertex(Double_t *vtx){
  //
  //
  //
  fPrimVtx[0]=vtx[0];
  fPrimVtx[1]=vtx[1];
  fPrimVtx[2]=vtx[2];
}


void AliTPCTransform::Transform(Double_t *x,Int_t *i,UInt_t /*time*/,
				Int_t /*coordinateType*/) {
  // input: x[0] - pad row
  //        x[1] - pad 
  //        x[2] - time in us
  //        i[0] - sector
  // output: x[0] - x (all in the rotated global coordinate frame)
  //         x[1] - y
  //         x[2] - z
  //
  //  primvtx     - position of the primary vertex
  //                used for the TOF correction
  //                TOF of particle calculated assuming the speed-of-light and 
  //                line approximation  
  //
  if (!fCurrentRecoParam) {
    return;
  }
  Int_t row=TMath::Nint(x[0]);
  Int_t pad=TMath::Nint(x[1]);
  Int_t sector=i[0];
  AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();  
  //
  AliTPCCalPad * time0TPC = calib->GetPadTime0(); 
  AliTPCCalPad * distortionMapY = calib->GetDistortionMap(0); 
  AliTPCCalPad * distortionMapZ = calib->GetDistortionMap(1); 
  AliTPCCalPad * distortionMapR = calib->GetDistortionMap(2); 
  AliTPCParam  * param    = calib->GetParameters(); 
  AliTPCCorrection * correction = calib->GetTPCComposedCorrection();   // first user defined correction  // if does not exist  try to get it from calibDB array
  if (!correction) correction = calib->GetTPCComposedCorrection(AliTracker::GetBz());
  if (!time0TPC){
    AliFatal("Time unisochronity missing");
    return ; // make coverity happy
  }
  AliTPCCorrection * correctionDelta = calib->GetTPCComposedCorrectionDelta(); 

  if (!param){
    AliFatal("Parameters missing");
    return; // make coverity happy
  }

  Double_t xx[3];
  //  Apply Time0 correction - Pad by pad fluctuation
  //  
  if (!calib->HasAlignmentOCDB()) x[2]-=time0TPC->GetCalROC(sector)->GetValue(row,pad);
  //
  // Tranform from pad - time coordinate system to the rotated global (tracking) system
  //
  Local2RotatedGlobal(sector,x);
  //
  //
  //
  // Alignment
  //TODO:  calib->GetParameters()->GetClusterMatrix(sector)->LocalToMaster(x,xx);
  RotatedGlobal2Global(sector,x);
  
  //
  // old ExB correction 
  //
  if(fCurrentRecoParam->GetUseExBCorrection()) {

    calib->GetExB()->Correct(x,xx);

  } else {

    xx[0] = x[0];
    xx[1] = x[1];
    xx[2] = x[2];
  }

  //
  // new composed  correction  - will replace soon ExB correction
  //
  if(fCurrentRecoParam->GetUseComposedCorrection()&&correction) {
    Float_t distPoint[3]={xx[0],xx[1],xx[2]};
    correction->CorrectPoint(distPoint, sector);
    xx[0]=distPoint[0];
    xx[1]=distPoint[1];
    xx[2]=distPoint[2];
    if (correctionDelta&&fCurrentRecoParam->GetUseAlignmentTime()){  // appply time dependent correction if available and enabled
      Float_t distPointDelta[3]={xx[0],xx[1],xx[2]};
      correctionDelta->CorrectPoint(distPointDelta, sector);
      xx[0]=distPointDelta[0];
      xx[1]=distPointDelta[1];
      xx[2]=distPointDelta[2];
    }
  } 


  //
  // Time of flight correction
  // 
  if (fCurrentRecoParam->GetUseTOFCorrection()){
    const Int_t kNIS=param->GetNInnerSector(), kNOS=param->GetNOuterSector(); 
    Float_t sign=1;
    if (sector < kNIS) {
      sign = (sector < kNIS/2) ? 1 : -1;
    } else {
      sign = ((sector-kNIS) < kNOS/2) ? 1 : -1;
    }
    Float_t deltaDr =0;
    Float_t dist=0;
    dist+=(fPrimVtx[0]-x[0])*(fPrimVtx[0]-x[0]);
    dist+=(fPrimVtx[1]-x[1])*(fPrimVtx[1]-x[1]);
    dist+=(fPrimVtx[2]-x[2])*(fPrimVtx[2]-x[2]);
    dist = TMath::Sqrt(dist);
    // drift length correction because of TOF
    // the drift velocity is in cm/s therefore multiplication by 0.01
    deltaDr = (dist*(0.01*param->GetDriftV()))/TMath::C(); 
    xx[2]+=sign*deltaDr;
  }
  //
  //
  //

  //
  Global2RotatedGlobal(sector,xx);

  //
  // Apply non linear distortion correction  
  //
  if (distortionMapY ){
    // wt - to get it form the OCDB
    // ignore T1 and T2
    AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
    Double_t bzField = magF->SolenoidField()/10.; //field in T
    Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
    Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
    if (sector%36<18) ezField*=-1;
    Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ;
    Double_t c0=1./(1.+wt*wt);
    Double_t c1=wt/c0;
    
    //can be switch on for each dimension separatelly
    if (fCurrentRecoParam->GetUseFieldCorrection()&0x2)
      if (distortionMapY){
	xx[1]-= c0*distortionMapY->GetCalROC(sector)->GetValue(row,pad);
	xx[0]-= c1*distortionMapY->GetCalROC(sector)->GetValue(row,pad);
      }
    if (fCurrentRecoParam->GetUseFieldCorrection()&0x4) 
      if (distortionMapZ)
	xx[2]-=distortionMapZ->GetCalROC(sector)->GetValue(row,pad);
    if (fCurrentRecoParam->GetUseFieldCorrection()&0x8) 
      if (distortionMapR){
	xx[0]-= c0*distortionMapR->GetCalROC(sector)->GetValue(row,pad);
	xx[1]-=-c1*distortionMapR->GetCalROC(sector)->GetValue(row,pad)*wt;
      }
    
  }
  //

  //
  x[0]=xx[0];x[1]=xx[1];x[2]=xx[2];
}

void AliTPCTransform::Local2RotatedGlobal(Int_t sector, Double_t *x) const {
  //
  //  
  // Tranform coordinate from  
  // row, pad, time to x,y,z
  //
  // Drift Velocity 
  // Current implementation - common drift velocity - for full chamber
  // TODO: use a map or parametrisation!
  //
  //  
  //
  if (!fCurrentRecoParam) return;
  const  Int_t kMax =60;  // cache for 60 seconds
  static Int_t lastStamp=-1;  //cached values
  static Double_t lastCorr = 1;
  //
  AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();
  AliTPCParam  * param    = calib->GetParameters(); 
  AliTPCCalibVdrift *driftCalib = AliTPCcalibDB::Instance()->GetVdrift(fCurrentRun);
  Double_t driftCorr = 1.;
  if (driftCalib){
    //
    // caching drift correction - temp. fix
    // Extremally slow procedure
    if ( TMath::Abs((lastStamp)-Int_t(fCurrentTimeStamp))<kMax){
      driftCorr = lastCorr;
    }else{
      driftCorr = 1.+(driftCalib->GetPTRelative(fCurrentTimeStamp,0)+ driftCalib->GetPTRelative(fCurrentTimeStamp,1))*0.5;
      lastCorr=driftCorr;
      lastStamp=fCurrentTimeStamp;
      
    }
  }
  //
  // simple caching non thread save
  static Double_t vdcorrectionTime=1;
  static Double_t vdcorrectionTimeGY=0;
  static Double_t time0corrTime=0;
  static Int_t    lastStampT=-1;
  //
  if (lastStampT!=(Int_t)fCurrentTimeStamp){
    lastStampT=fCurrentTimeStamp;
    if(fCurrentRecoParam->GetUseDriftCorrectionTime()>0) {
      vdcorrectionTime = (1+AliTPCcalibDB::Instance()->
			  GetVDriftCorrectionTime(fCurrentTimeStamp, 
						  fCurrentRun,
						  sector%36>=18,
						  fCurrentRecoParam->GetUseDriftCorrectionTime()));			  
      time0corrTime= AliTPCcalibDB::Instance()->
	GetTime0CorrectionTime(fCurrentTimeStamp, 
			       fCurrentRun,
			       sector%36>=18,
			       fCurrentRecoParam->GetUseDriftCorrectionTime());	
    }
    //
    if(fCurrentRecoParam->GetUseDriftCorrectionGY()>0) {
      
      Double_t corrGy= AliTPCcalibDB::Instance()->
			GetVDriftCorrectionGy(fCurrentTimeStamp, 
					      AliTPCcalibDB::Instance()->GetRun(),
					      sector%36>=18,
					      fCurrentRecoParam->GetUseDriftCorrectionGY());
      vdcorrectionTimeGY = corrGy;
    }
  }


  if (!param){
    AliFatal("Parameters missing");
    return; // make coverity happy
  }
  Int_t row=TMath::Nint(x[0]);
  //  Int_t pad=TMath::Nint(x[1]);
  //
  const Int_t kNIS=param->GetNInnerSector(), kNOS=param->GetNOuterSector();
  Double_t sign = 1.;
  Double_t zwidth    = param->GetZWidth()*driftCorr;
  Float_t xyzPad[3];
  AliTPCROC::Instance()->GetPositionGlobal(sector, TMath::Nint(x[0]) ,TMath::Nint(x[1]), xyzPad);
  if (AliTPCRecoParam:: GetUseTimeCalibration()) zwidth*=vdcorrectionTime*(1+xyzPad[1]*vdcorrectionTimeGY);
  Double_t padWidth  = 0;
  Double_t padLength = 0;
  Double_t    maxPad    = 0;
  //
  if (sector < kNIS) {
    maxPad = param->GetNPadsLow(row);
    sign = (sector < kNIS/2) ? 1 : -1;
    padLength = param->GetPadPitchLength(sector,row);
    padWidth = param->GetPadPitchWidth(sector);
  } else {
    maxPad = param->GetNPadsUp(row);
    sign = ((sector-kNIS) < kNOS/2) ? 1 : -1;
    padLength = param->GetPadPitchLength(sector,row);
    padWidth  = param->GetPadPitchWidth(sector);
  }
  //
  // X coordinate
  x[0] = param->GetPadRowRadii(sector,row);  // padrow X position - ideal
  //
  // Y coordinate
  //
  x[1]=(x[1]-0.5*maxPad)*padWidth;
  // pads are mirrorred on C-side
  if (sector%36>17){
    x[1]*=-1;
  }
  
  //
  
  //
  // Z coordinate
  //
  if (AliTPCcalibDB::Instance()->IsTrgL0()){
    // by defualt we assume L1 trigger is used - make a correction in case of  L0
    AliCTPTimeParams* ctp = AliTPCcalibDB::Instance()->GetCTPTimeParams();
    if (ctp){
      //for TPC standalone runs no ctp info
      Double_t delay = ctp->GetDelayL1L0()*0.000000025;
      x[2]-=delay/param->GetTSample();
    }
  }
  x[2]-= param->GetNTBinsL1();
  x[2]*= zwidth;  // tranform time bin to the distance to the ROC
  x[2]-= 3.*param->GetZSigma() + time0corrTime;
  // subtract the time offsets
  x[2] = sign*( param->GetZLength(sector) - x[2]);
}

void AliTPCTransform::RotatedGlobal2Global(Int_t sector,Double_t *x) const {
  //
  // transform possition rotated global to the global
  //
  Double_t cos,sin;
  GetCosAndSin(sector,cos,sin);
  Double_t tmp=x[0];
  x[0]= cos*tmp-sin*x[1];
  x[1]=+sin*tmp+cos*x[1];
}

void AliTPCTransform::Global2RotatedGlobal(Int_t sector,Double_t *x) const {
  //
  // tranform possition Global2RotatedGlobal
  //
  Double_t cos,sin;
  GetCosAndSin(sector,cos,sin);
  Double_t tmp=x[0];
  x[0]= cos*tmp+sin*x[1];
  x[1]= -sin*tmp+cos*x[1];
}

void AliTPCTransform::GetCosAndSin(Int_t sector,Double_t &cos,
					  Double_t &sin) const {
  cos=fCoss[sector%18];
  sin=fSins[sector%18];
}


void AliTPCTransform::ApplyTransformations(Double_t */*xyz*/, Int_t /*volID*/){
  //
  // Modify global position
  // xyz    - global xyz position
  // volID  - volID of detector (sector number)
  //
  //
  
}
