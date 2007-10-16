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
//                                                                           //
//  Comparison class for V0 information                                      //
//  responsible: 
//  marian.ivanov@cern.ch                                                    //
//
//

 



#include <stdio.h>
#include <string.h>
//ROOT includes
#include "Rtypes.h"
//
//ALIROOT includes
//
#include "AliESDtrack.h"
#include "AliTPCParam.h"
#include "AliTrackReference.h"
#include "AliTPCParamSR.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliTPCseed.h"
#include "AliITStrackMI.h"
#include "AliTRDtrack.h"
#include "AliHelix.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliESDkink.h"
#include "AliESDv0.h"
#include "AliV0.h"
//
#include "AliTreeDraw.h"
#include "AliMCInfo.h"
#include "AliGenKinkInfo.h"
#include "AliGenV0Info.h"


#include "AliESDRecV0Info.h"



ClassImp(AliESDRecV0Info)





void  AliESDRecV0Info::Update(Float_t vertex[3])
{ 

  if ( (fT1.fStatus[1]>0)&& (fT2.fStatus[1]>0)){
    Float_t distance1,distance2;
    Double_t xx[3],pp[3];
    //
    Double_t xd[3],pd[3],signd;
    Double_t xm[3],pm[3],signm;
    //
    //
    if (fT1.fITSOn&&fT2.fITSOn){
      for (Int_t i=0;i<3;i++){
	xd[i] = fT2.fITSinR1[i];
	pd[i] = fT2.fITSinP1[i];
	xm[i] = fT1.fITSinR1[i];
	pm[i] = fT1.fITSinP1[i];
      }
    }
    else{
      
      for (Int_t i=0;i<3;i++){
	xd[i] = fT2.fTPCinR1[i];
	pd[i] = fT2.fTPCinP1[i];
	xm[i] = fT1.fTPCinR1[i];
	pm[i] = fT1.fTPCinP1[i];
      }
    }
    //
    //
    signd =  fT2.fSign<0 ? -1:1;
    signm =  fT1.fSign<0 ? -1:1;

    AliHelix dhelix1(xd,pd,signd);
    dhelix1.GetMomentum(0,pp,0);
    dhelix1.Evaluate(0,xx);      
    // 
    //  Double_t x2[3],p2[3];
    //            
    AliHelix mhelix(xm,pm,signm);    
    //
    //find intersection linear
    //
    Double_t phase[2][2],radius[2];
    Int_t  points = dhelix1.GetRPHIintersections(mhelix, phase, radius,200);
    Double_t delta1=10000,delta2=10000;  

    if (points==1){
      fRs[0] = TMath::Sqrt(radius[0]);
      fRs[1] = TMath::Sqrt(radius[0]);
    }
    if (points==2){
      fRs[0] =TMath::Min(TMath::Sqrt(radius[0]),TMath::Sqrt(radius[1]));
      fRs[1] =TMath::Max(TMath::Sqrt(radius[0]),TMath::Sqrt(radius[1]));
    }
    
    if (points>0){
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    }
    if (points==2){    
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    }
    if (points==1){
      fRs[0] = TMath::Sqrt(radius[0]);
      fRs[1] = TMath::Sqrt(radius[0]);
      fDistMinR = delta1;
    }
    if (points==2){
      if (radius[0]<radius[1]){
	fRs[0] = TMath::Sqrt(radius[0]);
	fRs[1] = TMath::Sqrt(radius[1]);
	fDistMinR = delta1;
      }
      else{
	fRs[0] = TMath::Sqrt(radius[1]);
	fRs[1] = TMath::Sqrt(radius[0]);
	fDistMinR = delta2;
      }
    }
    //
    //
    distance1 = TMath::Min(delta1,delta2);
    //
    //find intersection parabolic
    //
    points = dhelix1.GetRPHIintersections(mhelix, phase, radius);
    delta1=10000,delta2=10000;  
    
    if (points>0){
      dhelix1.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    }
    if (points==2){    
      dhelix1.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    }
    
    distance2 = TMath::Min(delta1,delta2);
    if (distance2>100) fDist2 =100;
    return;
    if (delta1<delta2){
      //get V0 info
      dhelix1.Evaluate(phase[0][0],fXr);
      dhelix1.GetMomentum(phase[0][0],fPdr);
      mhelix.GetMomentum(phase[0][1],fPm);
      dhelix1.GetAngle(phase[0][0],mhelix,phase[0][1],fAngle);
      fRr = TMath::Sqrt(radius[0]);
    }
    else{
      dhelix1.Evaluate(phase[1][0],fXr);
      dhelix1.GetMomentum(phase[1][0], fPdr);
      mhelix.GetMomentum(phase[1][1], fPm);
      dhelix1.GetAngle(phase[1][0],mhelix,phase[1][1],fAngle);
      fRr = TMath::Sqrt(radius[1]);
    }
    fDist1 = TMath::Sqrt(distance1);
    fDist2 = TMath::Sqrt(distance2);      
    
    if (fDist2<10.5){
      Double_t x,alpha,param[5],cov[15];
      //
      fT1.GetESDtrack()->GetInnerExternalParameters(alpha,x,param);
      fT1.GetESDtrack()->GetInnerExternalCovariance(cov);
      AliExternalTrackParam paramm(x,alpha,param,cov);
      //
      fT2.GetESDtrack()->GetInnerExternalParameters(alpha,x,param);
      fT2.GetESDtrack()->GetInnerExternalCovariance(cov);
      AliExternalTrackParam paramd(x,alpha,param,cov);
    }    
    //            
    //   
    
    Float_t v[3] = {fXr[0]-vertex[0],fXr[1]-vertex[1],fXr[2]-vertex[2]};
    Float_t p[3] = {fPdr[0]+fPm[0], fPdr[1]+fPm[1],fPdr[2]+fPm[2]};
    
    Float_t vnorm2 = v[0]*v[0]+v[1]*v[1];
    Float_t vnorm3 = TMath::Sqrt(v[2]*v[2]+vnorm2);
    vnorm2 = TMath::Sqrt(vnorm2);
    Float_t pnorm2 = p[0]*p[0]+p[1]*p[1];
    Float_t pnorm3 = TMath::Sqrt(p[2]*p[2]+pnorm2);
    pnorm2 = TMath::Sqrt(pnorm2);
    
    fPointAngleFi = (v[0]*p[0]+v[1]*p[1])/(vnorm2*pnorm2);
    fPointAngleTh = (v[2]*p[2]+vnorm2*pnorm2)/(vnorm3*pnorm3);  
    fPointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3);
  }
}

