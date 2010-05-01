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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Component for redoing the reconstruction from the clusters and tracks
//      
//   The new calibration data used 
//
//   In reality it overwrites the content of the ESD 
//    

/*

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chainCl = tool.MakeChain("calib.txt","Clusters",0,1);
  chainCl->Lookup();
  TChain * chainTr = tool.MakeChain("calib.txt","Tracks",0,1);
  chainTr->Lookup();



*/



//  marian.ivanov@cern.ch
// 
#include "AliTPCcalibCalib.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include "TTimeStamp.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliTracker.h"
#include "AliTPCClusterParam.h"
#include "AliTPCParam.h"

#include "AliTPCcalibDB.h"
#include "AliTPCTransform.h"
#include "AliTPCRecoParam.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTPCPointCorrection.h"
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>

ClassImp(AliTPCcalibCalib)

AliTPCcalibCalib::AliTPCcalibCalib():
AliTPCcalibBase(),
  fApplyExBCorrection(1),      // apply ExB correction
  fApplyTOFCorrection(1),      // apply TOF correction
  fApplyPositionCorrection(1), // apply position correction
  fApplySectorAlignment(1),    // apply sector alignment
  fApplyRPhiCorrection(1),     // apply R-Phi correction
  fApplyRCorrection(1)         // apply Radial correction

{
  //
  // Constructor
  //
}


AliTPCcalibCalib::AliTPCcalibCalib(const Text_t *name, const Text_t *title) 
  :AliTPCcalibBase(),
  fApplyExBCorrection(1),      // apply ExB correction
  fApplyTOFCorrection(1),      // apply TOF correction
  fApplyPositionCorrection(1), // apply position correction
  fApplySectorAlignment(1),    // apply sector alignment
  fApplyRPhiCorrection(1),     // apply R-Phi correction
  fApplyRCorrection(1)         // apply Radial correction
{  
  SetName(name);
  SetTitle(title);
}


AliTPCcalibCalib::AliTPCcalibCalib(const AliTPCcalibCalib&calib):
  AliTPCcalibBase(calib),
  fApplyExBCorrection(calib.GetApplyExBCorrection()),
  fApplyTOFCorrection(calib.GetApplyTOFCorrection()),
  fApplyPositionCorrection(calib.GetApplyPositionCorrection()),
  fApplySectorAlignment(calib.GetApplySectorAlignment()),
  fApplyRPhiCorrection(calib.GetApplyRPhiCorrection()),
  fApplyRCorrection(calib.GetApplyRCorrection())

{
  //
  // copy constructor
  //
}

AliTPCcalibCalib &AliTPCcalibCalib::operator=(const AliTPCcalibCalib&calib){
  //
  //
  //
  ((AliTPCcalibBase *)this)->operator=(calib);
  return *this;
}


AliTPCcalibCalib::~AliTPCcalibCalib() {
  //
  // destructor
  //
}


void     AliTPCcalibCalib::Process(AliESDEvent *event){
  //
  // 
  //
  if (!event) {
    return;
  }  
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!ESDfriend) {
   return;
  }
  
  if (GetDebugLevel()>20) printf("Hallo world: Im here\n");
  Int_t ntracks=event->GetNumberOfTracks();   
  //AliTPCcalibDB::Instance()->SetExBField(fMagF);

  //
  //
  //

  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track = event->GetTrack(i);     
    AliESDfriendTrack *friendTrack = (AliESDfriendTrack*) ESDfriend->GetTrack(i);
    if (!friendTrack) continue;
    track->SetFriendTrack(friendTrack);
    const AliExternalTrackParam * trackIn  = track->GetInnerParam();
    const AliExternalTrackParam * trackOut = track->GetOuterParam();
    AliExternalTrackParam * tpcOut   = (AliExternalTrackParam *)friendTrack->GetTPCOut();
    if (!trackIn) continue;
    if (!trackOut) continue;
    if (!tpcOut) continue;   
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }
    if (!seed) continue;
    RefitTrack(track, seed, event->GetMagneticField());
    (*tpcOut)=*(track->GetOuterParam());  
  }
  return;
}

Bool_t  AliTPCcalibCalib::RefitTrack(AliESDtrack * track, AliTPCseed *seed, Float_t magesd){
  //
  // Refit track
  // if magesd==0 forget the curvature

  //
  // 0 - Setup transform object
  //
  AliESDfriendTrack *friendTrack = (AliESDfriendTrack *)track->GetFriendTrack();

  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
  AliTPCParam     *param     = AliTPCcalibDB::Instance()->GetParameters();
  transform->SetCurrentRun(fRun);
  transform->SetCurrentTimeStamp((UInt_t)fTime);
  if(!fApplyExBCorrection) { // disable ExB correction in transform
    if(transform->GetCurrentRecoParam())
      transform->GetCurrentRecoParamNonConst()->SetUseExBCorrection(0);
  }
  if(!fApplyTOFCorrection) { // disable TOF correction in transform
    if(transform->GetCurrentRecoParam())
      transform->GetCurrentRecoParamNonConst()->SetUseTOFCorrection(kFALSE);
  }

  //
  // First apply calibration
  //
  //  AliTPCPointCorrection * corr =  AliTPCPointCorrection::Instance();
  for (Int_t irow=0;irow<159;irow++) {
    AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);
    if (!cluster) continue; 
    AliTPCclusterMI cl0(*cluster);
    Double_t x[3]={cluster->GetRow(),cluster->GetPad(),cluster->GetTimeBin()};
    Int_t i[1]={cluster->GetDetector()};

    transform->Transform(x,i,0,1);
    //
    // get position correction
    //
    Int_t ipad=0;
    if (cluster->GetDetector()>35) ipad=1;
    Float_t dy =0;//AliTPCClusterParam::SPosCorrection(0,ipad,cluster->GetPad(),cluster->GetTimeBin(),cluster->GetZ(),cluster->GetSigmaY2(),cluster->GetSigmaZ2(),cluster->GetMax());
    Float_t dz =0;//AliTPCClusterParam::SPosCorrection(1,ipad,cluster->GetPad(),cluster->GetTimeBin(),cluster->GetZ(),cluster->GetSigmaY2(),cluster->GetSigmaZ2(),cluster->GetMax());
    // if(fApplyPositionCorrection) {
    //x[1]-=dy;
    //x[2]-=dz;
    // }
    //
    // Apply sector alignment
    //
    Double_t dxq = AliTPCPointCorrection::SGetCorrectionSector(0,cluster->GetDetector()%36,cluster->GetX(),
							       cluster->GetY(),cluster->GetZ());
    Double_t dyq = AliTPCPointCorrection::SGetCorrectionSector(1,cluster->GetDetector()%36,cluster->GetX(),
							       cluster->GetY(),cluster->GetZ());
    Double_t dzq = AliTPCPointCorrection::SGetCorrectionSector(2,cluster->GetDetector()%36,cluster->GetX(),
							       cluster->GetY(),cluster->GetZ());
    if (0&fApplySectorAlignment){
      x[0]-=dxq;
      x[1]-=dyq;
      x[2]-=dzq;
    }
//     //
//     // Apply r-phi correction  - To be done on track level- knowing the track angle !!!
//     //
//     Double_t corrclY =  
//       corr->RPhiCOGCorrection(cluster->GetDetector(),cluster->GetRow(), cluster->GetPad(),
// 				  cluster->GetY(),cluster->GetY(), cluster->GetZ(), 0., cluster->GetMax(),2.5);
//     // R correction
//     Double_t corrR   = corr->CorrectionOutR0(kFALSE,kFALSE,cluster->GetX(),cluster->GetY(),cluster->GetZ(),cluster->GetDetector());

//     if (0&fApplyRPhiCorrection){
//       if (cluster->GetY()>0) x[1]+=corrclY;  // rphi correction
//       if (cluster->GetY()<0) x[1]-=corrclY;  // rphi correction
//     }
//     if (0&fApplyRCorrection){      
//       x[0]+=corrR;                           // radial correction
//     }

    //
    //
    //
    cluster->SetX(x[0]);
    cluster->SetY(x[1]);
    cluster->SetZ(x[2]);
    
    //
    // Apply alignemnt
    //
    if (1){
      if (!param->IsGeoRead()) param->ReadGeoMatrices();
      TGeoHMatrix  *mat = param->GetClusterMatrix(cluster->GetDetector());
      //TGeoHMatrix  mat;
      Double_t pos[3]= {cluster->GetX(),cluster->GetY(),cluster->GetZ()};
      Double_t posC[3]={cluster->GetX(),cluster->GetY(),cluster->GetZ()};
      if (mat) mat->LocalToMaster(pos,posC);
      else{
	// chack Loading of Geo matrices from GeoManager - TEMPORARY FIX
      }
      cluster->SetX(posC[0]);
      cluster->SetY(posC[1]);
      cluster->SetZ(posC[2]);
    }



    if (fStreamLevel>2){
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	(*cstream)<<"Clusters"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "triggerClass="<<&fTriggerClass<<      //  trigger	  
	  "mag="<<fMagF<<             //  magnetic field
	  "cl0.="<<&cl0<<
	  "cl.="<<cluster<<
	  "cy="<<dy<<
	  "cz="<<dz<<
	  //	  "cY="<<corrclY<<
	  //	  "cR="<<corrR<<
	  "dxq="<<dxq<<
	  "dyq="<<dyq<<
	  "dzq="<<dzq<<
	  "\n";
      }
    }
  }
  //
  //
  //
  Int_t ncl = seed->GetNumberOfClusters();
  const Double_t kResetCov=4.;
  const Double_t kSigma=5.;
  Double_t covar[15];
  for (Int_t i=0;i<15;i++) covar[i]=0;
  covar[0]=kSigma*kSigma;
  covar[2]=kSigma*kSigma;
  covar[5]=kSigma*kSigma/Float_t(ncl*ncl);
  covar[9]=kSigma*kSigma/Float_t(ncl*ncl);
  covar[14]=0.2*0.2;
  if (TMath::Abs(magesd)<0.05) {
     covar[14]=0.025*0.025;
  }
  // 
  // And now do refit
  //
  AliExternalTrackParam * trackInOld  = (AliExternalTrackParam*)track->GetInnerParam();
  AliExternalTrackParam * trackOuter = (AliExternalTrackParam*)track->GetOuterParam();
  AliExternalTrackParam * trackOutOld   = (AliExternalTrackParam *)friendTrack->GetTPCOut();
    

  AliExternalTrackParam trackIn  = *trackOutOld;
  trackIn.ResetCovariance(kResetCov);
  trackIn.AddCovariance(covar);
  if (TMath::Abs(magesd)<0.05) {
    ((Double_t&)(trackIn.GetParameter()[4]))=0.000000001;
    ((Double_t&)(trackIn.GetCovariance()[14]))=covar[14];  // fix the line
  }  

  Double_t xyz[3];
  Int_t nclIn=0,nclOut=0;
  //
  // Refit in
  //

  for (Int_t irow=159; irow>0; irow--){
    AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
    if (!cl) continue;
    if (cl->GetX()<80) continue;
    Int_t sector = cl->GetDetector();
    Float_t dalpha = TMath::DegToRad()*(sector%18*20.+10.)-trackIn.GetAlpha();
    if (TMath::Abs(dalpha)>0.01){
      if (!trackIn.Rotate(TMath::DegToRad()*(sector%18*20.+10.))) break;
    }
    Double_t r[3]={cl->GetX(),cl->GetY(),cl->GetZ()};
    Double_t cov[3]={0.01,0.,0.01}; //TODO: correct error parametrisation
    AliTPCseed::GetError(cl, &trackIn,cov[0],cov[2]);
    cov[0]*=cov[0];
    cov[2]*=cov[2];
    trackIn.GetXYZ(xyz);
    Double_t bz = AliTracker::GetBz(xyz);

    if (!trackIn.PropagateTo(r[0],bz)) continue;
    if (RejectCluster(cl,&trackIn)) continue;
    nclIn++;
    trackIn.Update(&r[1],cov);    
  }
  //
  AliExternalTrackParam trackOut = trackIn;
  trackOut.ResetCovariance(kResetCov);
  trackOut.AddCovariance(covar);
  if (TMath::Abs(magesd)<0.05) {
    ((Double_t&)(trackOut.GetParameter()[4]))=0.000000001;
    ((Double_t&)(trackOut.GetCovariance()[14]))=covar[14];  // fix the line
  }

  //
  // Refit out
  //
  //Bool_t lastEdge=kFALSE;
  for (Int_t irow=0; irow<160; irow++){
    AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
    if (!cl) continue;
    if (cl->GetX()<80) continue;
    Int_t sector = cl->GetDetector();
    Float_t dalpha = TMath::DegToRad()*(sector%18*20.+10.)-trackOut.GetAlpha();

    if (TMath::Abs(dalpha)>0.01){
      if (!trackOut.Rotate(TMath::DegToRad()*(sector%18*20.+10.))) break;
    }
    Double_t r[3]={cl->GetX(),cl->GetY(),cl->GetZ()};

    Double_t cov[3]={0.01,0.,0.01}; //TODO: correct error parametrisation    
    AliTPCseed::GetError(cl, &trackOut,cov[0],cov[2]);
    cov[0]*=cov[0];
    cov[2]*=cov[2];
    trackOut.GetXYZ(xyz);
    Double_t bz = AliTracker::GetBz(xyz);
    if (!trackOut.PropagateTo(r[0],bz)) continue;
    if (RejectCluster(cl,&trackOut)) continue;
    nclOut++;
    trackOut.Update(&r[1],cov);    	
    //if (cl->GetType()<0) lastEdge=kTRUE;
    //if (cl->GetType()>=0) lastEdge=kFALSE;    
  }
  //
  //
  //
  nclIn=0;
  trackIn  = trackOut;
  trackIn.ResetCovariance(kResetCov);
  if (TMath::Abs(magesd)<0.05) {
    ((Double_t&)(trackIn.GetParameter()[4]))=0.000000001;
    ((Double_t&)(trackIn.GetCovariance()[14]))=covar[14];  // fix the line
  }
  //
  // Refit in one more time
  //
  for (Int_t irow=159; irow>0; irow--){
    AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
    if (!cl) continue;
    if (cl->GetX()<80) continue;
    Int_t sector = cl->GetDetector();
    Float_t dalpha = TMath::DegToRad()*(sector%18*20.+10.)-trackIn.GetAlpha();
    if (TMath::Abs(dalpha)>0.01){
      if (!trackIn.Rotate(TMath::DegToRad()*(sector%18*20.+10.))) break;
    }
    Double_t r[3]={cl->GetX(),cl->GetY(),cl->GetZ()};
    Double_t cov[3]={0.01,0.,0.01}; //TODO: correct error parametrisation
    AliTPCseed::GetError(cl, &trackIn,cov[0],cov[2]);
    cov[0]*=cov[0];
    cov[2]*=cov[2];
    trackIn.GetXYZ(xyz);
    Double_t bz = AliTracker::GetBz(xyz);

    if (!trackIn.PropagateTo(r[0],bz)) continue;
    if (RejectCluster(cl,&trackIn)) continue;
    nclIn++;
    trackIn.Update(&r[1],cov);    
  }


  trackIn.Rotate(trackInOld->GetAlpha());
  trackOut.Rotate(trackOutOld->GetAlpha());
  //
  trackInOld->GetXYZ(xyz);
  Double_t bz = AliTracker::GetBz(xyz);
  trackIn.PropagateTo(trackInOld->GetX(),bz);
  //
  trackOutOld->GetXYZ(xyz);
  bz = AliTracker::GetBz(xyz);  
  trackOut.PropagateTo(trackOutOld->GetX(),bz);
  

  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      (*cstream)<<"Tracks"<<
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"triggerClass="<<&fTriggerClass<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
	"nclIn="<<nclIn<<
	"nclOut="<<nclOut<<
	"ncl="<<ncl<<
	"TrIn0.="<<trackInOld<<
	"TrOut0.="<<trackOutOld<<
	"TrIn1.="<<&trackIn<<
	"TrOut1.="<<&trackOut<<
	"\n";
    }
  }
  //
  // And now rewrite ESDtrack and TPC seed
  //
 
  (*trackInOld)  = trackIn;
  (*trackOutOld) = trackOut;
  (*trackOuter) = trackOut;
  AliExternalTrackParam *t = &trackIn;
  track->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
  seed->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
  seed->SetNumberOfClusters((nclIn+nclOut)/2);
  return kTRUE;
}



Bool_t AliTPCcalibCalib::RejectCluster(AliTPCclusterMI* cl, AliExternalTrackParam * param){
  //
  // check the acceptance of cluster
  // Cut on edge effects
  //
  Float_t kEdgeCut=2.5;
  Float_t kSigmaCut=6;

  Bool_t isReject = kFALSE;
  Float_t edgeY = cl->GetX()*TMath::Tan(TMath::Pi()/18);
  Float_t dist  = edgeY - TMath::Abs(cl->GetY());
  if (param)  dist  = TMath::Abs(edgeY - TMath::Abs(param->GetY()));
  if (dist<kEdgeCut) isReject=kTRUE;

  Double_t cov[3]={0.01,0.,0.01}; //TODO: correct error parametrisation    
  AliTPCseed::GetError(cl, param,cov[0],cov[2]);
  if (param->GetSigmaY2()<0 || param->GetSigmaZ2()<0){
    AliError("Wrong parameters");
    return kFALSE;
  }
  Double_t py = (cl->GetY()-param->GetY())/TMath::Sqrt(cov[0]*cov[0]+param->GetSigmaY2());
  Double_t pz = (cl->GetZ()-param->GetZ())/TMath::Sqrt(cov[2]*cov[2]+param->GetSigmaZ2());
  //
  if ((py*py+pz*pz)>kSigmaCut*kSigmaCut) isReject=kTRUE;
  
  return isReject;
}



