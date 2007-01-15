/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// A straight line is coded as a point (3 Double_t) and           //
// 3 direction cosines                                           //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TTree.h>
#include <TMath.h>

#include "AliStrLine.h"

ClassImp(AliStrLine)

//________________________________________________________
AliStrLine::AliStrLine() :
  TObject(),
  fTpar(0),
  fDebug(0)
 {
  // Default constructor
  for(Int_t i=0;i<3;i++) {
    fP0[i] = 0.;
    fCd[i] = 0.;
  }
}

//________________________________________________________
AliStrLine::AliStrLine(Double_t *point, Double_t *cd,Bool_t twopoints) :
  TObject(),
  fTpar(0),
  fDebug(0)
{
  // Standard constructor
  // if twopoints is true:  point and cd are the 3D coordinates of
  //                        two points defininig the straight line
  // if twopoint is false: point represents the 3D coordinates of a point
  //                       belonging to the straight line and cd is the
  //                       direction in space
  if(twopoints){
    InitTwoPoints(point,cd);
  }
  else {
    InitDirection(point,cd);
  }
}

//________________________________________________________
AliStrLine::AliStrLine(Float_t *pointf, Float_t *cdf,Bool_t twopoints) :
  TObject(),
  fTpar(0),
  fDebug(0)
{
  // Standard constructor - with float arguments
  // if twopoints is true:  point and cd are the 3D coordinates of
  //                        two points defininig the straight line
  // if twopoint is false: point represents the 3D coordinates of a point
  //                       belonging to the straight line and cd is the
  //                       direction in space
  Double_t point[3];
  Double_t cd[3];
  for(Int_t i=0;i<3;i++){
    point[i] = pointf[i];
    cd[i] = cdf[i];
  }
  if(twopoints){
    InitTwoPoints(point,cd);
  }
  else {
    InitDirection(point,cd);
  }
}

//________________________________________________________
void AliStrLine::InitDirection(Double_t *point, Double_t *cd){
  // Initialization from a point and a direction
  Double_t norm = 0.;
  for(Int_t i=0;i<3;i++)norm+=cd[i]*cd[i];
  if(norm) {
    norm = TMath::Sqrt(norm);
    for(Int_t i=0;i<3;i++) cd[i]/=norm;
  }
  else {
    Error("AliStrLine","Null direction cosines!!!");
  }
  SetP0(point);
  SetCd(cd);
  fTpar = 0.;
  SetDebug();
}

//________________________________________________________
void AliStrLine::InitTwoPoints(Double_t *pA, Double_t *pB){
  // Initialization from the coordinates of two
  // points in the space
  Double_t cd[3];
  for(Int_t i=0;i<3;i++)cd[i] = pB[i]-pA[i];
  InitDirection(pA,cd);
}

//________________________________________________________
AliStrLine::~AliStrLine() {
  // destructor
}

//________________________________________________________
void AliStrLine::PrintStatus() const {
  // Print current status
  cout <<"=======================================================\n";
  cout <<"Direction cosines: ";
  for(Int_t i=0;i<3;i++)cout <<fCd[i]<<"; ";
  cout <<endl;
  cout <<"Known point: ";
  for(Int_t i=0;i<3;i++)cout <<fP0[i]<<"; ";
  cout <<endl;
  cout <<"Current value for the parameter: "<<fTpar<<endl;
  cout <<" Debug flag: "<<fDebug<<endl;
}

//________________________________________________________
Int_t AliStrLine::IsParallelTo(AliStrLine *line) const {
  // returns 1 if lines are parallel, 0 if not paralel
  Double_t cd2[3];
  line->GetCd(cd2);
  Double_t vecpx=fCd[1]*cd2[2]-fCd[2]*cd2[1];
  if(vecpx!=0) return 0;
  Double_t vecpy=-fCd[0]*cd2[2]+fCd[2]*cd2[0];
  if(vecpy!=0) return 0;
  Double_t vecpz=fCd[0]*cd2[1]-fCd[1]*cd2[0];
  if(vecpz!=0) return 0;
  return 1;
}
//________________________________________________________
Int_t AliStrLine::Crossrphi(AliStrLine *line){
  // Cross 2 lines in the X-Y plane
  Double_t p2[3];
  Double_t cd2[3];
  line->GetP0(p2);
  line->GetCd(cd2);
  Double_t a=fCd[0];
  Double_t b=-cd2[0];
  Double_t c=p2[0]-fP0[0];
  Double_t d=fCd[1];
  Double_t e=-cd2[1];
  Double_t f=p2[1]-fP0[1];
  Double_t deno = a*e-b*d;
  Int_t retcode = 0;
  if(deno != 0.) {
    fTpar = (c*e-b*f)/deno;
  }
  else {
    retcode = -1;
  }
  return retcode;
}

//________________________________________________________
Int_t AliStrLine::CrossPoints(AliStrLine *line, Double_t *point1, Double_t *point2){
  // Looks for the crossing point estimated starting from the
  // DCA segment
  Double_t p2[3];
  Double_t cd2[3];
  line->GetP0(p2);
  line->GetCd(cd2);
  Int_t i;
  Double_t k1 = 0;
  for(i=0;i<3;i++)k1+=(fP0[i]-p2[i])*fCd[i];
  Double_t k2 = 0;
  for(i=0;i<3;i++)k2+=(fP0[i]-p2[i])*cd2[i];
  Double_t a11 = 0;
  for(i=0;i<3;i++)a11+=fCd[i]*cd2[i];
  Double_t a22 = -a11;
  Double_t a21 = 0;
  for(i=0;i<3;i++)a21+=cd2[i]*cd2[i];
  Double_t a12 = 0;
  for(i=0;i<3;i++)a12-=fCd[i]*fCd[i];
  Double_t deno = a11*a22-a21*a12;
  if(deno == 0.) return -1;
  fTpar = (a11*k2-a21*k1) / deno;
  Double_t par2 = (k1*a22-k2*a12) / deno;
  line->SetPar(par2);
  GetCurrentPoint(point1);
  line->GetCurrentPoint(point2);
  return 0;
}
//________________________________________________________________
Int_t AliStrLine::Cross(AliStrLine *line, Double_t *point){

  //Finds intersection between lines
  Double_t point1[3];
  Double_t point2[3];
  Int_t retcod=CrossPoints(line,point1,point2);
  if(retcod==0){
    for(Int_t i=0;i<3;i++)point[i]=(point1[i]+point2[i])/2.;
    return 0;
  }else{
    return -1;
  }
}

//___________________________________________________________
Double_t AliStrLine::GetDCA(AliStrLine *line) const{
  //Returns the distance of closest approach between two lines
  Double_t p2[3];
  Double_t cd2[3];
  line->GetP0(p2);
  line->GetCd(cd2);
  Int_t i;
  Int_t ispar=IsParallelTo(line);
  if(ispar){
    Double_t dist1q=0,dist2=0,mod=0;
    for(i=0;i<3;i++){
      dist1q+=(fP0[i]-p2[i])*(fP0[i]-p2[i]);
      dist2+=(fP0[i]-p2[i])*fCd[i];
      mod+=fCd[i]*fCd[i];
    }
    if(mod!=0){
      dist2/=mod;
      return TMath::Sqrt(dist1q-dist2*dist2);
    }else{return -1;}
  }else{
     Double_t perp[3];
     perp[0]=fCd[1]*cd2[2]-fCd[2]*cd2[1];
     perp[1]=-fCd[0]*cd2[2]+fCd[2]*cd2[0];
     perp[2]=fCd[0]*cd2[1]-fCd[1]*cd2[0];
     Double_t mod=0,dist=0;
     for(i=0;i<3;i++){
       mod+=perp[i]*perp[i];
       dist+=(fP0[i]-p2[i])*perp[i];
     }
     mod=sqrt(mod);
     if(mod!=0){
       dist/=mod;
       return TMath::Abs(dist);
     }else{return -1;}
  }
}
//________________________________________________________
void AliStrLine::GetCurrentPoint(Double_t *point) const {
  // Fills the array point with the current value on the line
  for(Int_t i=0;i<3;i++)point[i]=fP0[i]+fCd[i]*fTpar;
}

//________________________________________________________
Double_t AliStrLine::GetDistFromPoint(Double_t *point) const {
  // computes distance from point 
  AliStrLine tmp(point,(Double_t *)fCd,kFALSE);
  return this->GetDCA(&tmp);
}
