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
//  Time Projection Chamber                                                  //
//  Comparison macro for reconstructed tracks - ESDs V0s                                     //
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
#include "AliESDRecKinkInfo.h"



ClassImp(AliESDRecKinkInfo)



AliESDRecKinkInfo::AliESDRecKinkInfo():
  fT1(),      //track1
  fT2(),      //track2  
  fKink(),    //kink
  fDist1(0),    //info about closest distance according closest MC - linear DCA
  fDist2(0),    //info about closest distance parabolic DCA
  fInvMass(0),  //reconstructed invariant mass -
  //
  fRr(0),       // rec position of the vertex 
  fMinR(0),     // minimum radius in rphi intersection
  fDistMinR(0), // distance at minimal radius
  fPointAngleFi(0), //point angle fi
  fPointAngleTh(0), //point angle theta
  fPointAngle(0),   //point angle full
  fStatus(0),       //status -tracks 
  fRecStatus(0),    //kink -status- 0 - not found  1-good -  fake
  fMultiple(0),     // how many times was kink reconstructed
  fKinkMultiple(0) // how many times was kink reconstructed
{
  //
  // Default constructor
  //
  for (Int_t i=0; i<3; i++) {
    fPdr[i] = 0;
    fXr[i] = 0;
    fPm[i] = 0;
    fAngle[i] = 0;
  }
  for (Int_t i=0; i<2; i++) { fLab[i] = 0; } 
}

////
void  AliESDRecKinkInfo::Update()
{

  if ( (fT1.fTPCOn)&& (fT2.fTPCOn)){
    //
    // IF BOTH RECONSTRUCTED
    Float_t distance1,distance2;
    Double_t xx[3],pp[3];
    //
    Double_t xd[3],pd[3],signd;
    Double_t xm[3],pm[3],signm;
    for (Int_t i=0;i<3;i++){
      xd[i] = fT2.fTPCinR1[i];
      pd[i] = fT2.fTPCinP1[i];
      xm[i] = fT1.fTPCinR1[i];
      pm[i] = fT1.fTPCinP1[i];
    }
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
    Double_t phase[2][2] = { {0,0}, {0,0}};
    Double_t radius[2] = {0} ;
    Int_t  points = dhelix1.GetRPHIintersections(mhelix, phase, radius,200);
    Double_t delta1=10000,delta2=10000;  

    if (points==1){
      fMinR = TMath::Sqrt(radius[0]);
    }
    if (points==2){
      fMinR =TMath::Min(TMath::Sqrt(radius[0]),TMath::Sqrt(radius[1]));
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
      fMinR = TMath::Sqrt(radius[0]);
      fDistMinR = delta1;
    }
    if (points==2){
      if (radius[0]<radius[1]){
	fMinR = TMath::Sqrt(radius[0]);
	fDistMinR = delta1;
      }
      else{
	fMinR = TMath::Sqrt(radius[1]);
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
      /*
      AliESDkink kink;
      kink.Update(&paramm,&paramd);
      //      kink.Dump();
      Double_t diff  = kink.fRr-fRr;
      Double_t diff2 = kink.fDist2-fDist2;
      printf("Diff\t%f\t%f\n",diff,diff2);
      */
    }
    
    //            
    //
  }

}

