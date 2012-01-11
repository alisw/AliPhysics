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
#include "AliHelix.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliESDkink.h"
#include "AliESDv0.h"
#include "AliV0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
//
#include "AliTreeDraw.h"
#include "AliMCInfo.h"
#include "AliGenKinkInfo.h"
#include "AliGenV0Info.h"


#include "AliESDRecV0Info.h"



ClassImp(AliESDRecV0Info)


AliESDRecV0Info:: AliESDRecV0Info():
  TObject(),
  fT1(),               //track1
  fT2(),               //track2  
  fDist1(0),           //info about closest distance according closest MC - linear DCA
  fDist2(0),           //info about closest distance parabolic DCA
  fInvMass(0),         //reconstructed invariant mass -
  //
  fDistMinR(0),        // distance at minimal radius
  fRr(0),              // rec position of the vertex 
  fPointAngleFi(0),    //point angle fi
  fPointAngleTh(0),    //point angle theta
  fPointAngle(0),      //point angle full
  fV0Status(0),        // status of the kink
  fV0tpc(0),           // Vo information from reconsturction according TPC
  fV0its(0),           // Vo information from reconsturction according ITS
  fV0rec(0),           // V0 information form the reconstruction
  fV0recOff(0),        // V0 information form the reconstruction - OFFLINE
  fMultiple(0),        // how man times V0 was recostructed 
  fRecStatus(0),       // status form the reconstuction
  fV0MultipleOn(0),    // how man times was V0 reconstucted
  fV0MultipleOff(0),   // how man times was V0 reconstucted
  //
  fKFrecChi2NC(0),     //  ONLINE V0 finder non constrained chi2  
  fKFrecChi2C(0),      //  ONLINE V0 finder   constrained chi2 - prim vertex  
  fKFrecChi2CM(0),     //  ONLINE V0 finder   constrained chi2 - prim vertex+mass 
  fKFRecNC(0),         //  non constrained  
  fKFRecC(0),          //  constrained vertex
  fKFRecCM(0),         //  constrained vertex+mass
  fKFrecOffChi2NC(0),  // OFFLINE V0 finder - non constrained chi2  
  fKFrecOffChi2C(0),   // OFFLINE V0 finder -     constrained chi2 - prim vertex  
  fKFrecOffChi2CM(0),  // OFFLINE V0 finder -     constrained chi2 - prim vertex+mass
  fKFOffRecNC(0),      //  non constrained  
  fKFOffRecC(0),       //  constrained vertex
  fKFOffRecCM(0)       //  constrained vertex+mass
{
  //
  // default constructor
  //
  for (Int_t i=0; i<3; i++) {
    fPdr[i] = 0;
    fXr[i] = 0;
    fPm[i] = 0;
    fAngle[i] = 0;
  }
  for (Int_t i=0; i<2; i++) {
    fRs[i] = 0;
    fLab[i] = 0;
  }

  fV0tpc    = new AliV0();
  fV0its    = new AliV0();
  fV0rec    = new AliV0();
  fV0recOff = new AliV0();
}


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
    Double_t phase[2][2] = { {0,0},{0,0} }; 
    Double_t radius[2] = {0};
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
    Bool_t checkAll = 0; // to be checked
    if(checkAll) { 
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
}


void  AliESDRecV0Info::Reset(){
  //
  // Reset status and all counters
  //
  fDist1=-1;    //info about closest distance according closest MC - linear DCA
  fDist2=-1;    //info about closest distance parabolic DCA
  fInvMass=-1;  //reconstructed invariant mass -
  //
  fDistMinR=-1; // distance at minimal radius
  fRr=-1;       // rec position of the vertex 
  fLab[0]=-20;   //MC label of the partecle
  fLab[1]=-10;   //MC label of the partecle
  fPointAngleFi=0; //point angle fi
  fPointAngleTh=0; //point angle theta
  fPointAngle=0;   //point angle full
  //
  fV0Status= -100;       // status of the V0
  fMultiple=0;           // how man times V0 was recostructed 
  fRecStatus=0;          // status form the reconstuction - 1 reconstructed - -1 fake
  fV0MultipleOn=0;       // how man times was V0 reconstucted - onfly
  fV0MultipleOff=0;      // how man times was V0 reconstucted - offline
  //
  // AliKF variables - variables to make a selection + resoluton study
  //
  fKFrecChi2NC=0;     //  ONLINE V0 finder non constrained chi2  
  fKFrecChi2C=0;      //  ONLINE V0 finder   constrained chi2 - prim vertex  
  fKFrecChi2CM=0;     //  ONLINE V0 finder   constrained chi2 - prim vertex+mass 
  //
  fKFrecOffChi2NC=0;  // OFFLINE V0 finder - non constrained chi2  
  fKFrecOffChi2C=0;   // OFFLINE V0 finder -     constrained chi2 - prim vertex  
  fKFrecOffChi2CM=0;  // OFFLINE V0 finder -     constrained chi2 - prim vertex+mass
}



void  AliESDRecV0Info::UpdateKF(const AliESDVertex &vertex, Int_t pdg0, Int_t pdg1, Float_t mass){
  //
  // Calculate properties of V0 vertex using different type of constraints 
  //
  fKFrecChi2NC=0;     //  ONLINE V0 finder non constrained chi2  
  fKFrecChi2C=0;      //  ONLINE V0 finder   constrained chi2 - prim vertex  
  fKFrecChi2CM=0;     //  ONLINE V0 finder   constrained chi2 - prim vertex+mass 
  if (fKFRecNC) {delete fKFRecNC; fKFRecNC=0;}
  if (fKFRecC)  {delete fKFRecC;  fKFRecC=0;}
  if (fKFRecCM) {delete fKFRecCM; fKFRecCM=0;}
  //
  fKFrecOffChi2NC=0;  // OFFLINE V0 finder - non constrained chi2  
  fKFrecOffChi2C=0;   // OFFLINE V0 finder -     constrained chi2 - prim vertex  
  fKFrecOffChi2CM=0;  // OFFLINE V0 finder -     constrained chi2 - prim vertex+mass
  if (fKFOffRecNC) {delete fKFOffRecNC; fKFOffRecNC=0;}
  if (fKFOffRecC)  {delete fKFOffRecC;  fKFOffRecC=0;}
  if (fKFOffRecCM) {delete fKFOffRecCM; fKFOffRecCM=0;}
  if (fV0Status==0) return; //
  //
  AliKFVertex primVtx(vertex);
  //
  if (fV0rec && 
      TMath::Abs(fV0rec->GetParamN()->GetSigmaY2())>0.000000001&&
      TMath::Abs(fV0rec->GetParamP()->GetSigmaY2())>0.000000001
      ){
    //
    Double_t x, y, z;
    AliKFParticle p1( *(fV0rec->GetParamN()), pdg0 );
    AliKFParticle p2( *(fV0rec->GetParamP()), pdg1 );
    //
    fKFRecNC  = new AliKFParticle;
    fV0rec->GetXYZ(x,y,z);
    fKFRecNC->SetVtxGuess(x,y,z);
    *(fKFRecNC)+=p1;
    *(fKFRecNC)+=p2;
    fKFrecChi2NC =fKFRecNC->GetChi2() ;
    //
    fKFRecC  = new AliKFParticle;
    fV0rec->GetXYZ(x,y,z);
    fKFRecC->SetVtxGuess(x,y,z);
    *(fKFRecC)+=p1;
    *(fKFRecC)+=p2;
    fKFRecC->SetProductionVertex(primVtx);
    fKFrecChi2C =fKFRecC->GetChi2();
    //
    fKFRecCM  = new AliKFParticle;
    fV0rec->GetXYZ(x,y,z);
    fKFRecCM->SetVtxGuess(x,y,z);
    *(fKFRecCM)+=p1;
    *(fKFRecCM)+=p2;
    fKFRecCM->SetProductionVertex(primVtx);
    fKFRecCM->SetMassConstraint(mass);
    fKFrecChi2CM =fKFRecCM->GetChi2();    
  }

  if (fV0recOff && 
      TMath::Abs(fV0recOff->GetParamN()->GetSigmaY2())>0.000000001&&
      TMath::Abs(fV0recOff->GetParamP()->GetSigmaY2())>0.000000001
      ){
    //
    Double_t x, y, z;
    AliKFParticle p1( *(fV0recOff->GetParamN()), pdg0 );
    AliKFParticle p2( *(fV0recOff->GetParamP()), pdg1 );
    //
    fKFOffRecNC  = new AliKFParticle;
    fV0recOff->GetXYZ(x,y,z);
    fKFOffRecNC->SetVtxGuess(x,y,z);
    *(fKFOffRecNC)+=p1;
    *(fKFOffRecNC)+=p2;
    fKFrecOffChi2NC =fKFOffRecNC->GetChi2() ;
    //
    fKFOffRecC  = new AliKFParticle;
    fV0recOff->GetXYZ(x,y,z);
    fKFOffRecC->SetVtxGuess(x,y,z);
    *(fKFOffRecC)+=p1;
    *(fKFOffRecC)+=p2;
    fKFOffRecC->SetProductionVertex(primVtx);
    fKFrecOffChi2C =fKFOffRecC->GetChi2();
    //
    fKFOffRecCM  = new AliKFParticle;
    fV0recOff->GetXYZ(x,y,z);
    fKFOffRecCM->SetVtxGuess(x,y,z);
    *(fKFOffRecCM)+=p1;
    *(fKFOffRecCM)+=p2;
    fKFOffRecCM->SetProductionVertex(primVtx);
    fKFOffRecCM->SetMassConstraint(mass);
    fKFrecOffChi2CM =fKFOffRecCM->GetChi2();    
  }


}


