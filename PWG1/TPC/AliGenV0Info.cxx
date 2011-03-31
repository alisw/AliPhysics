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


///////////////////////////////////////////////////////////////////////////
/*

Origin: marian.ivanov@cern.ch
Container classes with MC infomation for V0 


*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string.h>
//ROOT includes
#include "TROOT.h"
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"

//ALIROOT includes
#include "AliRun.h"
#include "AliStack.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliTPC.h"
#include "AliTPCLoader.h"
#include "AliDetector.h"
#include "AliTrackReference.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliMagF.h"
#include "AliHelix.h"
#include "AliTrackPointArray.h"

#endif
#include "AliGenV0Info.h" 
//
// 

ClassImp(AliGenV0Info)





/////////////////////////////////////////////////////////////////////////////////
AliGenV0Info::AliGenV0Info():
  fMCd(),      //info about daughter particle - 
  fMCm(),      //info about mother particle   - first particle for V0
  fMotherP(),   //particle info about mother particle
  fMCDist1(0), //info about closest distance according closest MC - linear DCA
  fMCDist2(0),    //info about closest distance parabolic DCA
  fMCRr(0),       // rec position of the vertex 
  fMCR(0),        //exact r position of the vertex
  fInvMass(0),  //reconstructed invariant mass -
  fPointAngleFi(0), //point angle fi
  fPointAngleTh(0), //point angle theta
  fPointAngle(0)   //point angle full
{
  for (Int_t i=0;i<3; i++){
   fMCPdr[i]=0;    
   fMCX[i]=0;     
   fMCXr[i]=0;    
   fMCPm[i]=0;    
   fMCAngle[i]=0; 
   fMCPd[i]=0;    
  }
  fMCPd[3]=0;     
  for (Int_t i=0; i<2;i++){
    fPdg[i]=0;   
    fLab[i]=0;  
  }
}

void AliGenV0Info::Update(Float_t vertex[3])
{
  //
  // Update information - calculates derived variables
  //
  
  fMCPd[0] = fMCd.GetParticle().Px();
  fMCPd[1] = fMCd.GetParticle().Py();
  fMCPd[2] = fMCd.GetParticle().Pz();
  fMCPd[3] = fMCd.GetParticle().P();
  //
  fMCX[0]  = fMCd.GetParticle().Vx();
  fMCX[1]  = fMCd.GetParticle().Vy();
  fMCX[2]  = fMCd.GetParticle().Vz();
  fMCR       = TMath::Sqrt( fMCX[0]*fMCX[0]+fMCX[1]*fMCX[1]);
  //
  fPdg[0]    = fMCd.GetParticle().GetPdgCode();
  fPdg[1]    = fMCm.GetParticle().GetPdgCode();
  //
  fLab[0]    = fMCd.GetParticle().GetUniqueID();
  fLab[1]    = fMCm.GetParticle().GetUniqueID();
  //
  //
  //
  Double_t x1[3],p1[3];
  TParticle& part0 = fMCd.GetParticle();
  x1[0] = part0.Vx();      
  x1[1] = part0.Vy();
  x1[2] = part0.Vz();
  p1[0] = part0.Px();
  p1[1] = part0.Py();
  p1[2] = part0.Pz();
  if (part0.GetPDG()==0) return;

  Double_t sign = (part0.GetPDG()->Charge()>0)? -1:1;
  AliHelix dhelix1(x1,p1,sign);
  //
  //
  Double_t x2[3],p2[3];            
  //
  TParticle& part1 = fMCm.GetParticle();
  if (part1.GetPDG()==0) return;
  x2[0] = part1.Vx();      
  x2[1] = part1.Vy();      
  x2[2] = part1.Vz();      
  p2[0] = part1.Px();
  p2[1] = part1.Py();
  p2[2] = part1.Pz();
  //
  //
  if (part1.GetPDG()==0) return;
  sign = (part1.GetPDG()->Charge()>0) ? -1:1;
  AliHelix mhelix(x2,p2,sign);
  
  //
  //
  //
  //find intersection linear
  //
  Double_t distance1, distance2;
  Double_t phase[2][2] = {{0,0},{0,0}};
  Double_t radius[2] = {0};
  Int_t  points = dhelix1.GetRPHIintersections(mhelix, phase, radius);
  Double_t delta1=10000,delta2=10000;    
  if (points>0){
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
  }
  else{
    fInvMass=-1;
    return;
  }
  if (points==2){    
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
  }
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
  //
  if (delta1<delta2){
    //get V0 info
    dhelix1.Evaluate(phase[0][0],fMCXr);
    dhelix1.GetMomentum(phase[0][0],fMCPdr);
    mhelix.GetMomentum(phase[0][1],fMCPm);
    dhelix1.GetAngle(phase[0][0],mhelix,phase[0][1],fMCAngle);
    fMCRr = TMath::Sqrt(radius[0]);
  }
  else{
    dhelix1.Evaluate(phase[1][0],fMCXr);
    dhelix1.GetMomentum(phase[1][0],fMCPdr);
    mhelix.GetMomentum(phase[1][1],fMCPm);
    dhelix1.GetAngle(phase[1][0],mhelix,phase[1][1],fMCAngle);
    fMCRr = TMath::Sqrt(radius[1]);
  }     
  //            
  //
  fMCDist1 = TMath::Sqrt(distance1);
  fMCDist2 = TMath::Sqrt(distance2);      
  //
  //
  // 
  Float_t v[3] = {fMCXr[0]-vertex[0],fMCXr[1]-vertex[1],fMCXr[2]-vertex[2]};
  //Float_t v[3] = {fMCXr[0],fMCXr[1],fMCXr[2]};
  Float_t p[3] = {fMCPdr[0]+fMCPm[0], fMCPdr[1]+fMCPm[1],fMCPdr[2]+fMCPm[2]};
  Float_t vnorm2 = v[0]*v[0]+v[1]*v[1];
  Float_t vnorm3 = TMath::Sqrt(v[2]*v[2]+vnorm2);
  vnorm2 = TMath::Sqrt(vnorm2);
  Float_t pnorm2 = p[0]*p[0]+p[1]*p[1];
  Float_t pnorm3 = TMath::Sqrt(p[2]*p[2]+pnorm2);
  pnorm2 = TMath::Sqrt(pnorm2);
  //
  if (vnorm2>0){
    fPointAngleFi = (v[0]*p[0]+v[1]*p[1])/(vnorm2*pnorm2);
    fPointAngleTh = (v[2]*p[2]+vnorm2*pnorm2)/(vnorm3*pnorm3);  
    fPointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3);
  }else{
    fPointAngleFi = 1;
    fPointAngleTh = 1;  
    fPointAngle   = 1; 
  }
  Double_t mass1 = fMCd.GetMass();
  Double_t mass2 = fMCm.GetMass();

  
  Double_t e1    = TMath::Sqrt(mass1*mass1+
			      fMCPd[0]*fMCPd[0]+
			      fMCPd[1]*fMCPd[1]+
			      fMCPd[2]*fMCPd[2]);
  Double_t e2    = TMath::Sqrt(mass2*mass2+
			      fMCPm[0]*fMCPm[0]+
			      fMCPm[1]*fMCPm[1]+
			      fMCPm[2]*fMCPm[2]);
  
  fInvMass =  
    (fMCPm[0]+fMCPd[0])*(fMCPm[0]+fMCPd[0])+
    (fMCPm[1]+fMCPd[1])*(fMCPm[1]+fMCPd[1])+
    (fMCPm[2]+fMCPd[2])*(fMCPm[2]+fMCPd[2]);
  
  //  fInvMass = TMath::Sqrt((e1+e2)*(e1+e2)-fInvMass);
  fInvMass = (e1+e2)*(e1+e2)-fInvMass;
  if (fInvMass>0) fInvMass = TMath::Sqrt(fInvMass);    
}








