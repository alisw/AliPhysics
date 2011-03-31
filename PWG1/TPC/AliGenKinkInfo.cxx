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
Container classes with MC infomation

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
#include "AliGenKinkInfo.h"

//
// 

ClassImp(AliGenKinkInfo)



/////////////////////////////////////////////////////////////////////////////////

AliGenKinkInfo::AliGenKinkInfo():
  fMCd(),       //info about daughter particle - second particle for V0
  fMCm(),       //info about mother particle   - first particle for V0
  fMCDist1(0),  //info about closest distance according closest MC - linear DCA
  fMCDist2(0),  //info about closest distance parabolic DCA
  fMCRr(0),     // rec position of the vertex 
  fMCR(0)       //exact r position of the vertex
{
  //
  // default constructor
  //
  for (Int_t i=0;i<3;i++){
    fMCPdr[i]=0;
    fMCPd[i]=0;
    fMCX[i]=0;
    fMCPm[i]=0;
    fMCAngle[i]=0;
    fMCXr[i] = 0;
  }
  for (Int_t i=0; i<2; i++) {
    fPdg[i]= 0;
    fLab[i]=0;
  }
}

void AliGenKinkInfo::Update()
{
  //
  // Update information
  //  some redundancy - faster acces to the values in analysis code
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
  TParticle& pdaughter = fMCd.GetParticle();
  x1[0] = pdaughter.Vx();      
  x1[1] = pdaughter.Vy();
  x1[2] = pdaughter.Vz();
  p1[0] = pdaughter.Px();
  p1[1] = pdaughter.Py();
  p1[2] = pdaughter.Pz();
  Double_t sign = (pdaughter.GetPDG()->Charge()>0)? -1:1;
  AliHelix dhelix1(x1,p1,sign);
  //
  //
  Double_t x2[3],p2[3];            
  //
  TParticle& pmother = fMCm.GetParticle();
  x2[0] = pmother.Vx();      
  x2[1] = pmother.Vy();      
  x2[2] = pmother.Vz();      
  p2[0] = pmother.Px();
  p2[1] = pmother.Py();
  p2[2] = pmother.Pz();
  //
  const AliTrackReference & pdecay = fMCm.GetTRdecay();
  x2[0] = pdecay.X();      
  x2[1] = pdecay.Y();      
  x2[2] = pdecay.Z();      
  p2[0] = pdecay.Px();
  p2[1] = pdecay.Py();
  p2[2] = pdecay.Pz();
  //
  sign = (pmother.GetPDG()->Charge()>0) ? -1:1;
  AliHelix mhelix(x2,p2,sign);
  
  //
  //
  //
  //find intersection linear
  //
  Double_t distance1, distance2;
  Double_t phase[2][2] = { {0,0},{0,0} };
  Double_t radius[2] = {0};
  Int_t  points = dhelix1.GetRPHIintersections(mhelix, phase, radius);
  Double_t delta1=10000,delta2=10000;    
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
    
}


Float_t AliGenKinkInfo::GetQt(){
  //
  //
  Float_t momentumd = TMath::Sqrt(fMCPd[0]*fMCPd[0]+fMCPd[1]*fMCPd[1]+fMCPd[2]*fMCPd[2]);
  return TMath::Sin(fMCAngle[2])*momentumd;
}


