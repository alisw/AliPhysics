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

/* $Id$ */

/*
  AliTPC parallel tracker - 
  How to use?  - 
  run AliTPCFindClusters.C macro - clusters neccessary for tracker are founded
  run AliTPCFindTracksMI.C macro - to find tracks
  tracks are written to AliTPCtracks.root file
  for comparison also seeds are written to the same file - to special branch
*/

//-------------------------------------------------------
//          Implementation of the TPC tracker
//
//   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
// 
//-------------------------------------------------------

#include <TObjArray.h>
#include <TFile.h>
#include <TTree.h>
#include "Riostream.h"

#include "AliTPCtrackerMI.h"
#include "AliTPCclusterMI.h"
#include "AliTPCParam.h"
#include "AliTPCClustersRow.h"
#include "AliComplexCluster.h"
#include "AliTPCpolyTrack.h"
#include "AliRunLoader.h"
#include "TStopwatch.h"


ClassImp(AliTPCseed)
ClassImp(AliTPCKalmanSegment)



//_____________________________________________________________________________

AliTPCKalmanSegment::AliTPCKalmanSegment(){
  //
  //
  fX=fAlpha=fChi2=0;
  for (Int_t i=0;i<5;i++) fState[i] = 0.;
  for (Int_t i=0;i<15;i++) fCovariance[i] = 0.;
  fNCFoundable = 0;
  fNC          = 0;
  //  fN           = 0;
}

//_____________________________________________________________________________

void AliTPCKalmanSegment::Init(AliTPCseed* seed)
{
  // in initialization  
  // initial entrance integral chi2, fNCFoundable and fNC stored 
  fNCFoundable = seed->fNFoundable;
  fNC          = seed->GetNumberOfClusters();
  fChi2        = seed->GetChi2();
}


void AliTPCKalmanSegment::Finish(AliTPCseed* seed)
{
  //
  // in finish state vector stored and chi2 and fNC... calculated
  Double_t x;  
  Double_t state[5];  
  Double_t cov[15];  
  seed->GetExternalParameters(x,state);
  seed->GetExternalCovariance(cov);
  //float precision for tracklet 
  for (Int_t i=0;i<5;i++) fState[i] = state[i];
  for (Int_t i=0;i<15;i++) fCovariance[i] = cov[i];
  //
  // in current   seed integral track characteristic 
  // for tracklet differenciation between beginning (Init state) and Finish state 
  fNCFoundable =  seed->fNFoundable - fNCFoundable; 
  fNC          =  seed->GetNumberOfClusters() - fNC;
  fChi2        =  seed->GetChi2()-fChi2;
  //
}

void AliTPCKalmanSegment::GetState(Double_t &x, Double_t & alpha, Double_t state[5])
{
  //
  x = fX;
  alpha = fAlpha;
  for (Int_t i=0;i<5;i++) state[i] = fState[i];
} 
       
void AliTPCKalmanSegment::GetCovariance(Double_t covariance[15])
{
  //
  for (Int_t i=0;i<5;i++) covariance[i] = fCovariance[i];
}

void AliTPCKalmanSegment::GetStatistic(Int_t & nclusters, Int_t & nfoundable, Float_t & chi2)
{
  //
  //
  nclusters  = fNC;
  nfoundable = fNCFoundable;
  chi2       = fChi2;  
}



AliTPCclusterTracks::AliTPCclusterTracks(){
  // class for storing overlaping info
  fTrackIndex[0]=-1;
  fTrackIndex[1]=-1;
  fTrackIndex[2]=-1;
  fDistance[0]=1000;
  fDistance[1]=1000;
  fDistance[2]=1000;
}









Int_t AliTPCtrackerMI::UpdateTrack(AliTPCseed * track, AliTPCclusterMI* c, Double_t chi2, UInt_t i){

  Int_t sec=(i&0xff000000)>>24; 
  Int_t row = (i&0x00ff0000)>>16;
  track->fRow=(i&0x00ff0000)>>16;
  track->fSector = sec;
  //  Int_t index = i&0xFFFF;
  if (sec>=fParam->GetNInnerSector()) track->fRow += fParam->GetNRowLow(); 
  track->fClusterIndex[track->fRow] = i;
  track->fFirstPoint = row;
  if ( track->fLastPoint<row) track->fLastPoint =row;
  //

  AliTPCTrackPoint   *trpoint =track->GetTrackPoint(track->fRow);
  Float_t angle2 = track->GetSnp()*track->GetSnp();
  angle2 = TMath::Sqrt(angle2/(1-angle2)); 
  //
  //SET NEW Track Point
  //
  if (c!=0){
    //if we have a cluster
    trpoint->GetCPoint().SetY(c->GetY());
    trpoint->GetCPoint().SetZ(c->GetZ());    
    //
    trpoint->GetCPoint().SetSigmaY(c->GetSigmaY2()/(track->fCurrentSigmaY*track->fCurrentSigmaY));
    trpoint->GetCPoint().SetSigmaZ(c->GetSigmaZ2()/(track->fCurrentSigmaZ*track->fCurrentSigmaZ));
    //
    trpoint->GetCPoint().SetType(c->GetType());
    trpoint->GetCPoint().SetQ(c->GetQ());
    trpoint->GetCPoint().SetMax(c->GetMax());
    //  
    trpoint->GetCPoint().SetErrY(TMath::Sqrt(track->fErrorY2));
    trpoint->GetCPoint().SetErrZ(TMath::Sqrt(track->fErrorZ2));
    //
  }
  trpoint->GetTPoint().SetX(track->GetX());
  trpoint->GetTPoint().SetY(track->GetY());
  trpoint->GetTPoint().SetZ(track->GetZ());
  //
  trpoint->GetTPoint().SetAngleY(angle2);
  trpoint->GetTPoint().SetAngleZ(track->GetTgl());
  

  
  if (chi2>10){
    //    printf("suspicious chi2 %f\n",chi2);
  }
  //  if (track->fIsSeeding){ 
    track->fErrorY2 *= 1.2;
    track->fErrorY2 += 0.0064;    
    track->fErrorZ2 *= 1.2;   
    track->fErrorY2 += 0.005;    
 
    //}

  return track->Update(c,chi2,i);

}
//_____________________________________________________________________________
AliTPCtrackerMI::AliTPCtrackerMI(const AliTPCParam *par): 
AliTracker(), fkNIS(par->GetNInnerSector()/2), fkNOS(par->GetNOuterSector()/2)
{
  //---------------------------------------------------------------------
  // The main TPC tracker constructor
  //---------------------------------------------------------------------
  fInnerSec=new AliTPCSector[fkNIS];         
  fOuterSec=new AliTPCSector[fkNOS];

  Int_t i;
  for (i=0; i<fkNIS; i++) fInnerSec[i].Setup(par,0);
  for (i=0; i<fkNOS; i++) fOuterSec[i].Setup(par,1);

  fN=0;  fSectors=0;

  fClustersArray.Setup(par);
  fClustersArray.SetClusterType("AliTPCclusterMI");

  fSeeds=0;
  fNtracks = 0;
  fParam = par;
}

//_____________________________________________________________________________
AliTPCtrackerMI::~AliTPCtrackerMI() {
  //------------------------------------------------------------------
  // TPC tracker destructor
  //------------------------------------------------------------------
  delete[] fInnerSec;
  delete[] fOuterSec;
  if (fSeeds) {
    fSeeds->Delete(); 
    delete fSeeds;
  }
}


Double_t AliTPCtrackerMI::ErrY2(AliTPCseed* seed, AliTPCclusterMI * cl){
  //
  //
  Float_t snoise2;
  Float_t z = TMath::Abs(fParam->GetZLength()-TMath::Abs(seed->GetZ()));

  //cluster "quality"
  Float_t rsigmay = 1;
  Int_t ctype = 0;

  //standard if we don't have cluster - take MIP
  const Float_t chmip = 50.; 
  Float_t amp = chmip/0.3;  
  Float_t nel;
  Float_t nprim;
  if (cl){
    amp = cl->GetQ();
    rsigmay = cl->GetSigmaY2()/(seed->fCurrentSigmaY*seed->fCurrentSigmaY);
    ctype = cl->GetType();
  }
  

  Float_t landau=2 ;    //landau fluctuation part
  Float_t gg=2;         // gg fluctuation part
  Float_t padlength= fSectors->GetPadPitchLength(seed->GetX());


  if (fSectors==fInnerSec){
    snoise2 = 0.0004;
    nel     = 0.268*amp;
    nprim   = 0.155*amp;
    gg      = (2+0.0002*amp)/nel;
    landau  = (2.+0.12*nprim)*0.5*(amp*amp/40000.+2)/nprim;
    if (landau>1) landau=1;
  }
  else {
    snoise2 = 0.0004;
    nel     = 0.3*amp;
    nprim   = 0.133*amp;
    gg      = (2+0.0002*amp)/nel;
    landau  = (2.+0.12*nprim)*0.5*(amp*amp/40000.+2)/nprim;
    if (landau>1) landau=1;
  }


  Float_t sdiff = gg*fParam->GetDiffT()*fParam->GetDiffT()*z;
  Float_t angle2 = seed->GetSnp()*seed->GetSnp();
  angle2 = angle2/(1-angle2); 
  Float_t angular = landau*angle2*padlength*padlength/12.;
  Float_t res = sdiff + angular;
  
  
  if ((ctype==0) && (fSectors ==fOuterSec))
    res *= 0.78 +TMath::Exp(7.4*(rsigmay-1.2));

  if ((ctype==0) && (fSectors ==fInnerSec))
    res *= 0.72 +TMath::Exp(3.36*(rsigmay-1.2));


  if ((ctype>0))
    res*= TMath::Power((rsigmay+0.5),1.5)+0.0064;
  
  if (ctype<0)
    res*=2.4;  // overestimate error 2 times
  
  res+= snoise2;
 
  if (res<2*snoise2)
    res = 2*snoise2;

  seed->SetErrorY2(res);
  return res;
}




Double_t AliTPCtrackerMI::ErrZ2(AliTPCseed* seed, AliTPCclusterMI * cl){
  //
  //
  Float_t snoise2;
  Float_t z = TMath::Abs(fParam->GetZLength()-TMath::Abs(seed->GetZ()));
  //signal quality
  Float_t rsigmaz=1;
  Int_t ctype =0;

  const Float_t chmip = 50.;
  Float_t amp = chmip/0.3;  
  Float_t nel;
  Float_t nprim;
  if (cl){
    amp = cl->GetQ();
    rsigmaz = cl->GetSigmaZ2()/(seed->fCurrentSigmaZ*seed->fCurrentSigmaZ);
    ctype = cl->GetType();
  }

  //
  Float_t landau=2 ;    //landau fluctuation part
  Float_t gg=2;         // gg fluctuation part
  Float_t padlength= fSectors->GetPadPitchLength(seed->GetX());

  if (fSectors==fInnerSec){
    snoise2 = 0.0004;
    nel     = 0.268*amp;
    nprim   = 0.155*amp;
    gg      = (2+0.0002*amp)/nel;
    landau  = (2.+0.12*nprim)*0.5*(amp*amp/40000.+2)/nprim;
    if (landau>1) landau=1;
  }
  else {
    snoise2 = 0.0004;
    nel     = 0.3*amp;
    nprim   = 0.133*amp;
    gg      = (2+0.0002*amp)/nel;
    landau  = (2.+0.12*nprim)*0.5*(amp*amp/40000.+2)/nprim;
    if (landau>1) landau=1;
  }
  Float_t sdiff = gg*fParam->GetDiffT()*fParam->GetDiffT()*z;
  
  Float_t angle = seed->GetTgl();
  Float_t angular = landau*angle*angle*padlength*padlength/12.;
  Float_t res = sdiff + angular;

  if ((ctype==0) && (fSectors ==fOuterSec))
    res *= 0.81 +TMath::Exp(6.8*(rsigmaz-1.2));

  if ((ctype==0) && (fSectors ==fInnerSec))
    res *= 0.72 +TMath::Exp(2.04*(rsigmaz-1.2));
  if ((ctype>0))
    res*= TMath::Power(rsigmaz+0.5,1.5)+0.0064;  //0.31+0.147*ctype;
  if (ctype<0)
    res*=1.3;
  if ((ctype<0) &&amp<70)
    res*=1.3;  

  res += snoise2;
  if (res<2*snoise2)
     res = 2*snoise2;

  seed->SetErrorZ2(res);
  return res;
}





void AliTPCseed::Reset()
{
  //
  //PH  SetN(0);
  fNFoundable = 0;
  ResetCovariance();
  SetChi2(0);
  for (Int_t i=0;i<200;i++) fClusterIndex[i]=-1;
}


Int_t  AliTPCseed::GetProlongation(Double_t xk, Double_t &y, Double_t & z) const
{
  //-----------------------------------------------------------------
  // This function find proloncation of a track to a reference plane x=xk.
  // doesn't change internal state of the track
  //-----------------------------------------------------------------
  
  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1;
  //  Double_t y1=fP0, z1=fP1;
  Double_t c1=fP4*x1 - fP2, r1=sqrt(1.- c1*c1);
  Double_t c2=fP4*x2 - fP2, r2=sqrt(1.- c2*c2);
  
  y = fP0;
  z = fP1;
  y += dx*(c1+c2)/(r1+r2);
  z += dx*(c1+c2)/(c1*r2 + c2*r1)*fP3;
  return 0;  
}


//_____________________________________________________________________________
Double_t AliTPCseed::GetPredictedChi2(const AliTPCclusterMI *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  //Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  Double_t r00=fErrorY2, r01=0., r11=fErrorZ2;
  r00+=fC00; r01+=fC10; r11+=fC11;

  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1.e-10) {
    Int_t n=GetNumberOfClusters();
    if (n>4) cerr<<n<<" AliKalmanTrack warning: Singular matrix !\n";
    return 1e10;
  }
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  
  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  
  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}


//_________________________________________________________________________________________


Int_t AliTPCseed::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the sector - for given sector according z
  //-----------------------------------------------------------------
  AliTPCseed *t=(AliTPCseed*)o;
  if (t->fRelativeSector>fRelativeSector) return -1;
  if (t->fRelativeSector<fRelativeSector) return 1;

  Double_t z2 = t->GetZ();
  Double_t z1 = GetZ();
  if (z2>z1) return 1;
  if (z2<z1) return -1;
  return 0;
}

void AliTPCtrackerMI::RotateToLocal(AliTPCseed *seed)
{
  //rotate to track "local coordinata
  Float_t x = seed->GetX();
  Float_t y = seed->GetY();
  Float_t ymax = x*TMath::Tan(0.5*fSectors->GetAlpha());
  if (y > ymax) {
    seed->fRelativeSector= (seed->fRelativeSector+1) % fN;
    if (!seed->Rotate(fSectors->GetAlpha())) 
      return;
  } else if (y <-ymax) {
    seed->fRelativeSector= (seed->fRelativeSector-1+fN) % fN;
    if (!seed->Rotate(-fSectors->GetAlpha())) 
      return;
  }   

}




//_____________________________________________________________________________
Int_t AliTPCseed::Update(const AliTPCclusterMI *c, Double_t chisq, UInt_t index) {
  //-----------------------------------------------------------------
  // This function associates a cluster with this track.
  //-----------------------------------------------------------------
  //  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  //Double_t r00=sigmay2, r01=0., r11=sigmaz2;
  Double_t r00=fErrorY2, r01=0., r11=fErrorZ2;

  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

  Double_t k00=fC00*r00+fC10*r01, k01=fC00*r01+fC10*r11;
  Double_t k10=fC10*r00+fC11*r01, k11=fC10*r01+fC11*r11;
  Double_t k20=fC20*r00+fC21*r01, k21=fC20*r01+fC21*r11;
  Double_t k30=fC30*r00+fC31*r01, k31=fC30*r01+fC31*r11;
  Double_t k40=fC40*r00+fC41*r01, k41=fC40*r01+fC41*r11;

  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  Double_t cur=fP4 + k40*dy + k41*dz, eta=fP2 + k20*dy + k21*dz;
  if (TMath::Abs(cur*fX-eta) >= 0.9) {
    //    Int_t n=GetNumberOfClusters();
    //if (n>4) cerr<<n<<" AliTPCtrack warning: Filtering failed !\n";
    return 0;
  }

  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2  = eta;
  fP3 += k30*dy + k31*dz;
  fP4  = cur;

  Double_t c01=fC10, c02=fC20, c03=fC30, c04=fC40;
  Double_t c12=fC21, c13=fC31, c14=fC41;

  fC00-=k00*fC00+k01*fC10; fC10-=k00*c01+k01*fC11;
  fC20-=k00*c02+k01*c12;   fC30-=k00*c03+k01*c13;
  fC40-=k00*c04+k01*c14; 

  fC11-=k10*c01+k11*fC11;
  fC21-=k10*c02+k11*c12;   fC31-=k10*c03+k11*c13;
  fC41-=k10*c04+k11*c14; 

  fC22-=k20*c02+k21*c12;   fC32-=k20*c03+k21*c13;
  fC42-=k20*c04+k21*c14; 

  fC33-=k30*c03+k31*c13;
  fC43-=k40*c03+k41*c13; 

  fC44-=k40*c04+k41*c14; 

  Int_t n=GetNumberOfClusters();
  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chisq);

  return 1;
}



//_____________________________________________________________________________
Double_t AliTPCtrackerMI::f1(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -xr*yr/sqrt(xr*xr+yr*yr); 
}


//_____________________________________________________________________________
Double_t AliTPCtrackerMI::f2(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature times center of curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  
  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

//_____________________________________________________________________________
Double_t AliTPCtrackerMI::f3(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}


Int_t AliTPCtrackerMI::LoadClusters()
{
  //
  // load clusters to the memory
  Int_t j=Int_t(fClustersArray.GetTree()->GetEntries());
  for (Int_t i=0; i<j; i++) {
    fClustersArray.LoadEntry(i);
  }

  LoadOuterSectors();
  LoadInnerSectors();

  return 0;
}

void AliTPCtrackerMI::UnloadClusters()
{
  //
  // load clusters to the memory
  Int_t j=Int_t(fClustersArray.GetTree()->GetEntries());
  for (Int_t i=0; i<j; i++) {
    fClustersArray.ClearSegment(i);
  }
}



//_____________________________________________________________________________
void AliTPCtrackerMI::LoadOuterSectors() {
  //-----------------------------------------------------------------
  // This function fills outer TPC sectors with clusters.
  //-----------------------------------------------------------------
  UInt_t index;
  //Int_t j=Int_t(fClustersArray.GetTree()->GetEntries());
  Int_t j = ((AliTPCParam*)fParam)->GetNRowsTotal();
  for (Int_t i=0; i<j; i++) {
    //  AliSegmentID *s=fClustersArray.LoadEntry(i);
    AliSegmentID *s= const_cast<AliSegmentID*>(fClustersArray.At(i));
    if (!s) continue;
    Int_t sec,row;
    AliTPCParam *par=(AliTPCParam*)fClustersArray.GetParam();
    par->AdjustSectorRow(s->GetID(),sec,row);
    if (sec<fkNIS*2) continue;
    AliTPCClustersRow *clrow=fClustersArray.GetRow(sec,row);
    Int_t ncl=clrow->GetArray()->GetEntriesFast();
    while (ncl--) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(*clrow)[ncl];
      index=(((sec<<8)+row)<<16)+ncl;
      fOuterSec[(sec-fkNIS*2)%fkNOS][row].InsertCluster(c,index);
    }
  }  
  fN=fkNOS;
  fSectors=fOuterSec;
}


//_____________________________________________________________________________
void AliTPCtrackerMI::LoadInnerSectors() {
  //-----------------------------------------------------------------
  // This function fills inner TPC sectors with clusters.
  //-----------------------------------------------------------------
  UInt_t index;
  //Int_t j=Int_t(fClustersArray.GetTree()->GetEntries());
  Int_t j = ((AliTPCParam*)fParam)->GetNRowsTotal();
  for (Int_t i=0; i<j; i++) {
    //   AliSegmentID *s=fClustersArray.LoadEntry(i);
    AliSegmentID *s=const_cast<AliSegmentID*>(fClustersArray.At(i));
    if (!s) continue;
    Int_t sec,row;
    AliTPCParam *par=(AliTPCParam*)fClustersArray.GetParam();
    par->AdjustSectorRow(s->GetID(),sec,row);
    if (sec>=fkNIS*2) continue;
    AliTPCClustersRow *clrow=fClustersArray.GetRow(sec,row);
    Int_t ncl=clrow->GetArray()->GetEntriesFast();
    while (ncl--) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(*clrow)[ncl];
      index=(((sec<<8)+row)<<16)+ncl;
      fInnerSec[sec%fkNIS][row].InsertCluster(c,index);
    }
  }

  fN=fkNIS;
  fSectors=fInnerSec;
}

Int_t AliTPCtrackerMI::FollowToNext(AliTPCseed& t, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------
  //  Double_t xt=t.GetX();
  //  Int_t row = fSectors->GetRowNumber(xt)-1;
  //  if (row < nr) return 1; // don't prolongate if not information until now -
  //
  Double_t  x=fSectors->GetX(nr), ymax=fSectors->GetMaxY(nr);
  //  if (t.GetRadius()>x+10 ) return 0;

  if (!t.PropagateTo(x)) {
    t.fStopped = kTRUE;
    return 0;
  }
  // update current
  t.fCurrentSigmaY = GetSigmaY(&t);
  t.fCurrentSigmaZ = GetSigmaZ(&t);
  //  
  AliTPCclusterMI *cl=0;
  UInt_t index=0;
  const AliTPCRow &krow=fSectors[t.fRelativeSector][nr];
  Double_t sy2=ErrY2(&t)*2;
  Double_t sz2=ErrZ2(&t)*2;


  Double_t  roady  =3.*sqrt(t.GetSigmaY2() + sy2);
  Double_t  roadz = 3 *sqrt(t.GetSigmaZ2() + sz2);
  Double_t  y=t.GetY(), z=t.GetZ();

  if (TMath::Abs(TMath::Abs(y)-ymax)<krow.fDeadZone){
    t.fInDead = kTRUE;
    Int_t row = nr;
    if (fSectors==fOuterSec) row += fParam->GetNRowLow();
    t.fClusterIndex[row] = -1; 
    return 0;
  } 
  else
    {
      if (TMath::Abs(z)<(1.05*x+10)) t.fNFoundable++;
      else
	return 0;
    }   
  //calculate 
  Float_t maxdistance = roady*roady + roadz*roadz;
  if (krow) {
    for (Int_t i=krow.Find(z-roadz); i<krow; i++) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(krow[i]);
      if (c->GetZ() > z+roadz) break;
      if ( (c->GetY()-y) >  roady ) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)+(c->GetY()-y)*(c->GetY()-y);
      if (maxdistance>distance) {
	maxdistance = distance;
	cl=c;
	//	index=krow.GetIndex(i);       
	index =i;
      }
    }
  }      
  if (cl) {
    //    Double_t sy2= ErrY2(&t,cl);
    //    Double_t sz2= ErrZ2(&t,cl);
    //    Double_t chi2= t.GetPredictedChi2(cl);    
    //    UpdateTrack(&t,cl,chi2,index);   
   
    t.fCurrentCluster = cl; 
    t.fCurrentClusterIndex1 = krow.GetIndex(index);   
    t.fCurrentClusterIndex2 = index;   
    Double_t sy2=ErrY2(&t,t.fCurrentCluster);
    Double_t sz2=ErrZ2(&t,t.fCurrentCluster);

    Double_t sdistancey = TMath::Sqrt(sy2+t.GetSigmaY2());
    Double_t sdistancez = TMath::Sqrt(sz2+t.GetSigmaZ2());

    Double_t rdistancey = TMath::Abs(t.fCurrentCluster->GetY()-t.GetY());
    Double_t rdistancez = TMath::Abs(t.fCurrentCluster->GetZ()-t.GetZ());
    
    Double_t rdistance  = TMath::Sqrt(TMath::Power(rdistancey/sdistancey,2)+TMath::Power(rdistancez/sdistancez,2));


    //    printf("\t%f\t%f\t%f\n",rdistancey/sdistancey,rdistancez/sdistancez,rdistance);
    if ( (rdistancey>1) || (rdistancez>1)) return 0;
    if (rdistance>4) return 0;

    if ((rdistancey/sdistancey>2.5 || rdistancez/sdistancez>2.5) && t.fCurrentCluster->GetType()==0)  
	return 0;  //suspisiouce - will be changed

    if ((rdistancey/sdistancey>2. || rdistancez/sdistancez>2.0) && t.fCurrentCluster->GetType()>0)  
	// strict cut on overlaped cluster
	return 0;  //suspisiouce - will be changed

    if ( (rdistancey/sdistancey>1. || rdistancez/sdistancez>2.5 ||t.fCurrentCluster->GetQ()<70 ) 
	 && t.fCurrentCluster->GetType()<0)
      return 0;

    //    t.SetSampledEdx(0.3*t.fCurrentCluster->GetQ()/l,t.GetNumberOfClusters(), GetSigmaY(&t), GetSigmaZ(&t));
    UpdateTrack(&t,t.fCurrentCluster,t.GetPredictedChi2(t.fCurrentCluster),t.fCurrentClusterIndex1);

  } else {    
    if (y > ymax) {
      t.fRelativeSector= (t.fRelativeSector+1) % fN;
      if (!t.Rotate(fSectors->GetAlpha())) 
	return 0;
    } else if (y <-ymax) {
      t.fRelativeSector= (t.fRelativeSector-1+fN) % fN;
      if (!t.Rotate(-fSectors->GetAlpha())) 
	return 0;
    }   
  }
  return 1;
}


Int_t AliTPCtrackerMI::UpdateClusters(AliTPCseed& t,Int_t trindex,  Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------
  t.fCurrentCluster  = 0;
  t.fCurrentClusterIndex1 = 0;   
  t.fCurrentClusterIndex2 = 0;
   
  Double_t xt=t.GetX();
  Int_t row = fSectors->GetRowNumber(xt)-1;
  if (row < nr) return 1; // don't prolongate if not information until now -
  Double_t x=fSectors->GetX(nr);
  //  if (t.fStopped) return 0;
  //  if (t.GetRadius()>x+10 ) return 0;
  if (!t.PropagateTo(x)){
    t.fStopped =kTRUE;
    return 0;
  }
  // update current
  t.fCurrentSigmaY = GetSigmaY(&t);
  t.fCurrentSigmaZ = GetSigmaZ(&t);
    
  AliTPCclusterMI *cl=0;
  UInt_t index=0;
  AliTPCRow &krow=fSectors[t.fRelativeSector][nr];
  //
  Double_t  y=t.GetY(), z=t.GetZ();
  Double_t roady = 3.* TMath::Sqrt(t.GetSigmaY2() + t.fCurrentSigmaY*t.fCurrentSigmaY);
  Double_t roadz = 3.* TMath::Sqrt(t.GetSigmaZ2() + t.fCurrentSigmaZ*t.fCurrentSigmaZ);
  //

  Float_t maxdistance = 1000000;
  if (krow) {    
    for (Int_t i=krow.Find(z-roadz); i<krow; i++) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(krow[i]);
      if (c->GetZ() > z+roadz) break;
      if (TMath::Abs(c->GetY()-y)>roady) continue;
            
      //krow.UpdateClusterTrack(i,trindex,&t);

      Float_t dy2 = (c->GetY()- t.GetY());
      dy2*=dy2;
      Float_t dz2 = (c->GetZ()- t.GetZ());
      dz2*=dz2;
      //
      Float_t distance = dy2+dz2;
      //
      if (distance > maxdistance) continue;
      maxdistance = distance;
      cl=c;
      index=i;       
    }
  }
  t.fCurrentCluster  = cl;
  t.fCurrentClusterIndex1 = krow.GetIndex(index);   
  t.fCurrentClusterIndex2 = index;   
  return 1;
}


Int_t AliTPCtrackerMI::FollowToNextCluster(Int_t trindex, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------
  AliTPCseed & t  = *((AliTPCseed*)(fSeeds->At(trindex)));
  AliTPCRow &krow=fSectors[t.fRelativeSector][nr];
  //  Double_t pt=t.GetConvConst()/(100/0.299792458/0.2)/t.Get1Pt();
  Double_t y=t.GetY();
  Double_t ymax=fSectors->GetMaxY(nr);

  if (TMath::Abs(TMath::Abs(y)-ymax)<krow.fDeadZone){
    t.fInDead = kTRUE;
    Int_t row = nr;
    if (fSectors==fOuterSec) row += fParam->GetNRowLow();
    t.fClusterIndex[row] = -1; 
    return 0;
  } 
  else
    {
      if (TMath::Abs(t.GetZ())<(1.05*t.GetX()+10)) t.fNFoundable++;
      else
	return 0;      
    }
  
  if (t.fCurrentCluster) {
    //    Float_t l=fSectors->GetPadPitchLength();
    //    AliTPCclusterTracks * cltrack = krow.GetClusterTracks(t.fCurrentClusterIndex1);

    Double_t sy2=ErrY2(&t,t.fCurrentCluster);
    Double_t sz2=ErrZ2(&t,t.fCurrentCluster);


    Double_t sdistancey = TMath::Sqrt(sy2+t.GetSigmaY2());
    Double_t sdistancez = TMath::Sqrt(sz2+t.GetSigmaZ2());

    Double_t rdistancey = TMath::Abs(t.fCurrentCluster->GetY()-t.GetY());
    Double_t rdistancez = TMath::Abs(t.fCurrentCluster->GetZ()-t.GetZ());
    
    Double_t rdistance  = TMath::Sqrt(TMath::Power(rdistancey/sdistancey,2)+TMath::Power(rdistancez/sdistancez,2));


    //    printf("\t%f\t%f\t%f\n",rdistancey/sdistancey,rdistancez/sdistancez,rdistance);
    if ( (rdistancey>1) || (rdistancez>1)) return 0;
    if (rdistance>4) return 0;

    if ((rdistancey/sdistancey>2.5 || rdistancez/sdistancez>2.5) && t.fCurrentCluster->GetType()==0)  
	return 0;  //suspisiouce - will be changed

    if ((rdistancey/sdistancey>2. || rdistancez/sdistancez>2.0) && t.fCurrentCluster->GetType()>0)  
	// strict cut on overlaped cluster
	return 0;  //suspisiouce - will be changed

    if ( (rdistancey/sdistancey>1. || rdistancez/sdistancez>2.5 ||t.fCurrentCluster->GetQ()<70 ) 
	 && t.fCurrentCluster->GetType()<0)
      return 0;

    //    t.SetSampledEdx(0.3*t.fCurrentCluster->GetQ()/l,t.GetNumberOfClusters(), GetSigmaY(&t), GetSigmaZ(&t));
    UpdateTrack(&t,t.fCurrentCluster,t.GetPredictedChi2(t.fCurrentCluster),t.fCurrentClusterIndex1);
   
  } else {
    if (y > ymax) {
      t.fRelativeSector= (t.fRelativeSector+1) % fN;
      if (!t.Rotate(fSectors->GetAlpha())) 
	return 0;
    } else if (y <-ymax) {
      t.fRelativeSector= (t.fRelativeSector-1+fN) % fN;
      if (!t.Rotate(-fSectors->GetAlpha())) 
	return 0;
    }
  }
  return 1;
}




/*
Int_t AliTPCtrackerMI::FollowProlongationFast(AliTPCseed& t, Int_t step)
{
  //-----------------------------------------------------------------
  // fast prolongation mathod -
  // don't update track only after step clusters
  //-----------------------------------------------------------------
  Double_t xt=t.GetX();
  //
  Double_t alpha=t.GetAlpha(); 
  alpha =- fSectors->GetAlphaShift();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  t.fRelativeSector = Int_t(alpha/fSectors->GetAlpha())%fN;
  Int_t row0 = fSectors->GetRowNumber(xt); 
  Double_t x    = fSectors->GetX(row0);
  Double_t ymax = fSectors->GetMaxY(row0);
  //
  Double_t sy2=ErrY2(&t)*2;
  Double_t sz2=ErrZ2(&t)*2;
  Double_t  roady  =3.*sqrt(t.GetSigmaY2() + sy2);
  Double_t  roadz = 3 *sqrt(t.GetSigmaZ2() + sz2);
  Float_t maxdistance = roady*roady + roadz*roadz; 
  t.fCurrentSigmaY = GetSigmaY(&t);
  t.fCurrentSigmaZ = GetSigmaZ(&t);
  //
  Int_t nclusters = 0;
  Double_t y;
  Double_t z;
  Double_t yy[200];   //track prolongation
  Double_t zz[200];
  Double_t cy[200];  // founded cluster position
  Double_t cz[200];
  Double_t sy[200];  // founded cluster error
  Double_t sz[200];
  Bool_t   hitted[200]; // indication of cluster presence
  //
  
  //
  for (Int_t drow = step; drow>=0; drow--) {
    Int_t row = row0-drow;
    if (row<0) break;
    Double_t x    = fSectors->GetX(row);
    Double_t ymax = fSectors->GetMaxY(row);
    t.GetProlongation(x,y,z);
    yy[drow] =y;
    zz[drow] =z;    
    const AliTPCRow &krow=fSectors[t.fRelativeSector][row];
    if (TMath::Abs(TMath::Abs(y)-ymax)<krow.fDeadZone){
      t.fInDead = kTRUE;
      break;
    } 
    else
      {
	t.fNFoundable++;
      }  
    
    //find nearest  cluster 
    AliTPCclusterMI *cl= 0;
    if (krow) {
      for (Int_t i=krow.Find(z-roadz); i<krow; i++) {
	AliTPCclusterMI *c=(AliTPCclusterMI*)(krow[i]);
	if (c->GetZ() > z+roadz) break;
	if ( (c->GetY()-y) >  roady ) continue;
	Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)+(c->GetY()-y)*(c->GetY()-y);
	if (maxdistance>distance) {
	  maxdistance = distance;
	  cl=c;
	  //	  index=krow.GetIndex(i);       
	}	
      }              
    }  // end of seearch
    //update cluster information
    if (cl){ 
      cy[drow] = cl->GetY();
      cz[drow] = cl->GetZ();
      sy[drow] = ErrY2(&t,cl);
      sz[drow] = ErrZ2(&t,cl);
      hitted[drow] = kTRUE;
      nclusters++;
    }
    else
      hitted[drow] = kFALSE;
  }
  //if we have information - update track
  if (nclusters>0){
    Float_t sumyw0   = 0;
    Float_t sumzw0   = 0;
    Float_t sumyw   = 0;
    Float_t sumzw   = 0;
    Float_t sumrow  = 0;
    for (Int_t i=0;i<step;i++)
      if (hitted[i]){
	sumyw0+= 1/sy[i];
	sumzw0+= 1/sz[i];
	//
	sumyw+= (cy[i]-yy[i])/sy[i];
	sumzw+= (cz[i]-zz[i])/sz[i];
	sumrow+= i;
      }
    Float_t dy = sumyw/sumyw0;
    Float_t dz = sumzw/sumzw0;
    Float_t mrow = sumrow/nclusters+row0;
    Float_t x = fSectors->GetX(mrow);
    t.PropagateTo(x);
    AliTPCclusterMI cvirtual;
    cvirtual.SetZ(dz+t.GetZ());
    cvirtual.SetY(dy+t.GetY());
    t.SetErrorY2(1.2*t.fErrorY2/TMath::Sqrt(Float_t(nclusters)));
    t.SetErrorZ2(1.2*t.fErrorZ2/TMath::Sqrt(Float_t(nclusters)));
    Float_t chi2 = t.GetPredictedChi2(&cvirtual);
    t.Update(&cvirtual,chi2,0);
    Int_t ncl = t.GetNumberOfClusters();
    ncl = ncl-1+nclusters;
    t.SetN(ncl);
  }     
  return  1;
}   
*/

//_____________________________________________________________________________
Int_t AliTPCtrackerMI::FollowProlongation(AliTPCseed& t, Int_t rf) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  Double_t xt=t.GetX();
  //
  Double_t alpha=t.GetAlpha() - fSectors->GetAlphaShift();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  t.fRelativeSector = Int_t(alpha/fSectors->GetAlpha())%fN;
    
  for (Int_t nr=fSectors->GetRowNumber(xt)-1; nr>=rf; nr--) {
   
    if (FollowToNext(t,nr)==0) {
    }
  }   
  return 1;
}


//_____________________________________________________________________________
Int_t AliTPCtrackerMI::FollowBackProlongation(AliTPCseed& t, Int_t rf) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  Double_t xt=t.GetX();
  //
  Double_t alpha=t.GetAlpha() - fSectors->GetAlphaShift();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  t.fRelativeSector = Int_t(alpha/fSectors->GetAlpha())%fN;
    
  for (Int_t nr=fSectors->GetRowNumber(xt)+1; nr<=rf; nr++) {
    if (t.GetSnp()<0.8)
      FollowToNext(t,nr);							      
  }   
  return 1;
}




   
Float_t AliTPCtrackerMI::OverlapFactor(AliTPCseed * s1, AliTPCseed * s2, Int_t &sum1, Int_t & sum2)
{
  //
  //
  sum1=0;
  sum2=0;
  Int_t sum=0;
  //
  Float_t dz2 =(s1->GetZ() - s2->GetZ());
  dz2*=dz2;  
  /*
    Float_t x = s1->GetX();
    Float_t x2 = s2->GetX();
    
    Float_t ymax = x*TMath::Tan(0.5*fSectors->GetAlpha());
  */
  Float_t dy2 =TMath::Abs((s1->GetY() - s2->GetY()));
  //if (TMath::Abs(dy2)>2*ymax-3) 
  //  dy2-=2*ymax;
  dy2*=dy2;
  Float_t distance = TMath::Sqrt(dz2+dy2);
  if (distance>4.) return 0; // if there are far away  - not overlap - to reduce combinatorics
 
  Int_t offset =0;
  if (fSectors==fOuterSec) offset = fParam->GetNRowLow();
  Int_t firstpoint = TMath::Min(s1->fFirstPoint,s2->fFirstPoint);
  Int_t lastpoint = TMath::Max(s1->fLastPoint,s2->fLastPoint);
  lastpoint +=offset;
  firstpoint+=offset;
  if (lastpoint>160) 
    lastpoint =160;
  if (firstpoint<0) 
    firstpoint = 0;
  if (firstpoint<lastpoint-15) {
    firstpoint =0;
    lastpoint  =160;
  }
    
  
  for (Int_t i=firstpoint;i<lastpoint;i++){
    if (s1->fClusterIndex[i]>0) sum1++;
    if (s2->fClusterIndex[i]>0) sum2++;
    if (s1->fClusterIndex[i]==s2->fClusterIndex[i] && s1->fClusterIndex[i]>0) {
      sum++;
    }
  }
 
  Float_t summin = TMath::Min(sum1+1,sum2+1);
  Float_t ratio = (sum+1)/Float_t(summin);
  return ratio;
}

void  AliTPCtrackerMI::SignShared(AliTPCseed * s1, AliTPCseed * s2)
{
  //
  //
  if (s1->fSector!=s2->fSector) return;
  //
  Float_t dz2 =(s1->GetZ() - s2->GetZ());
  dz2*=dz2;
  Float_t dy2 =(s1->GetY() - s2->GetY());

  dy2*=dy2;
  Float_t distance = TMath::Sqrt(dz2+dy2);
  if (distance>15.) return ; // if there are far away  - not overlap - to reduce combinatorics
  //trpoint = new (pointarray[track->fRow]) AliTPCTrackPoint;
  //  TClonesArray &pointarray1 = *(s1->fPoints);
  //TClonesArray &pointarray2 = *(s2->fPoints);
  //
  for (Int_t i=0;i<160;i++){
    if (s1->fClusterIndex[i]==s2->fClusterIndex[i] && s1->fClusterIndex[i]>0) {
      //  AliTPCTrackPoint *p1  = (AliTPCTrackPoint *)(pointarray1.UncheckedAt(i));
      //AliTPCTrackPoint *p2  = (AliTPCTrackPoint *)(pointarray2.UncheckedAt(i)); 
      AliTPCTrackPoint *p1  = s1->GetTrackPoint(i);
      AliTPCTrackPoint *p2  = s2->GetTrackPoint(i);; 
      p1->fIsShared = kTRUE;
      p2->fIsShared = kTRUE;
    }
  } 
}




void  AliTPCtrackerMI::RemoveOverlap(TObjArray * arr, Float_t factor, Int_t removalindex , Bool_t shared){

  

  //
  // remove overlap - used removal factor - removal index stored in the track
  for (Int_t i=0; i<arr->GetEntriesFast(); i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (pt) RotateToLocal(pt);
  }
  arr->Sort();  // sorting according z
  arr->Expand(arr->GetEntries());
  Int_t nseed=arr->GetEntriesFast();
  //  printf("seeds \t%p \t%d\n",arr, nseed);
  //  arr->Expand(arr->GetEntries());  //remove 0 pointers
  nseed = arr->GetEntriesFast();
  Int_t removed = 0;
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    if (!(pt->IsActive())) continue;
    for (Int_t j=i+1; j<nseed; j++){
      AliTPCseed *pt2=(AliTPCseed*)arr->UncheckedAt(j);
      //
      if (!pt2) continue; 
      if (!(pt2->IsActive())) continue;
      if (TMath::Abs(pt->fRelativeSector-pt2->fRelativeSector)>0) break;
      if (TMath::Abs(pt2->GetZ()-pt->GetZ())<4){
	Int_t sum1,sum2;
	Float_t ratio = OverlapFactor(pt,pt2,sum1,sum2);
	//if (sum1==0) {
	//  pt->Desactivate(removalindex); // arr->RemoveAt(i); 
	//  break;
	//}
	if (ratio>factor){
	  //	  if (pt->GetChi2()<pt2->GetChi2()) pt2->Desactivate(removalindex);  // arr->RemoveAt(j);	      
	  Float_t ratio2 = (pt->GetChi2()*sum2)/(pt2->GetChi2()*sum1);
	  Float_t ratio3 = Float_t(sum1-sum2)/Float_t(sum1+sum2);
	  removed++;
	  if (TMath::Abs(ratio3)>0.025){  // if much more points  
	    if (sum1>sum2) pt2->Desactivate(removalindex);
	    else {
	      pt->Desactivate(removalindex); // arr->RemoveAt(i); 
	      break;
	    }
	  }
	  else{  //decide on mean chi2
	    if (ratio2<1)  
	      pt2->Desactivate(removalindex);
	    else {
	      pt->Desactivate(removalindex); // arr->RemoveAt(i); 
	      break;
	    }	    
	  }  
	  
	}  // if suspicious ratio
      }
      else
	break;
    }
  }
  //  printf("removed\t%d\n",removed);
  Int_t good =0; 
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) break;
    if (pt->GetNumberOfClusters() < pt->fNFoundable*0.5) {
      //desactivate tracks with small number of points
      //      printf("%d\t%d\t%f\n", pt->GetNumberOfClusters(), pt->fNFoundable,pt->GetNumberOfClusters()/Float_t(pt->fNFoundable));
      pt->Desactivate(10);  //desactivate  - small muber of points
    }
    if (!(pt->IsActive())) continue;
    good++;
  }
  
  
  if (shared)
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
      if (!pt) continue;
      if (!(pt->IsActive())) continue;
      for (Int_t j=i+1; j<nseed; j++){
	AliTPCseed *pt2=(AliTPCseed*)arr->UncheckedAt(j);
	if ((pt2) && pt2->IsActive()) {
	  if ( TMath::Abs(pt->fSector-pt2->fSector)>1) break;
	  SignShared(pt,pt2);
	}
      }
    }
  fNtracks = good;
  printf("\n*****\nNumber of good tracks after overlap removal\t%d\n",fNtracks);


}

void AliTPCtrackerMI::RemoveUsed(TObjArray * arr, Float_t factor, Int_t removalindex)
{

  //Loop over all tracks and remove "overlaps"
  //
  //
  Int_t nseed = arr->GetEntriesFast();  
  Int_t good =0;
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    if (!(pt->IsActive())) continue;
    Int_t noc=pt->GetNumberOfClusters();
    Int_t shared =0;
    for (Int_t i=0; i<noc; i++) {
      Int_t index=pt->GetClusterIndex(i);
      AliTPCclusterMI *c=(AliTPCclusterMI*)GetClusterMI(index); 
      if (!c) continue;
      if (c->IsUsed()) shared++;
    }
    if ((Float_t(shared)/Float_t(noc))>factor)
      pt->Desactivate(removalindex);
    else{
      good++;
      for (Int_t i=0; i<noc; i++) {
	Int_t index=pt->GetClusterIndex(i);
	AliTPCclusterMI *c=(AliTPCclusterMI*)GetClusterMI(index);  
	if (!c) continue;
	c->Use();  
      }
    }
  }
  fNtracks = good;
  printf("\n*****\nNumber of good tracks after shared removal\t%d\n",fNtracks);

}


void AliTPCtrackerMI::MakeSeedsAll()
{
  if (fSeeds == 0) fSeeds = new TObjArray;
  TObjArray * arr;
  for (Int_t sec=0;sec<fkNOS;sec+=3){
     arr = MakeSeedsSectors(sec,sec+3);
     Int_t nseed = arr->GetEntriesFast();
     for (Int_t i=0;i<nseed;i++)
       fSeeds->AddLast(arr->RemoveAt(i));
  }  
  //  fSeeds = MakeSeedsSectors(0,fkNOS);
}

TObjArray *  AliTPCtrackerMI::MakeSeedsSectors(Int_t sec1, Int_t sec2)
{
  //
  // loop over all  sectors and make seed
  //find track seeds
  TStopwatch timer;
  Int_t nup=fOuterSec->GetNRows(), nlow=fInnerSec->GetNRows();
  Int_t nrows=nlow+nup;  
  Int_t gap=Int_t(0.125*nrows), shift=Int_t(0.5*gap);
  //  if (fSeeds==0) fSeeds = new TObjArray;
  TObjArray * arr = new TObjArray;

  for (Int_t sec=sec1; sec<sec2;sec++){
    MakeSeeds(arr, sec, nup-1, nup-1-gap);
    MakeSeeds(arr, sec, nup-2-shift, nup-2-shift-gap);
  }
  gap = Int_t(0.2* nrows);
  for (Int_t sec=sec1; sec<sec2;sec++){    
    //find secondaries
    //MakeSeeds2(arr, sec, nup-1, nup-1-gap);
    MakeSeeds2(arr, sec, nup-1-shift, nup-1-shift-gap);
    //MakeSeeds2(arr, sec, nup-1-2*shift, nup-1-2*shift-gap);
    MakeSeeds2(arr, sec, nup-1-3*shift, nup-1-3*shift-gap);
    //MakeSeeds2(arr, sec, nup-1-4*shift, nup-1-4*shift-gap);
    MakeSeeds2(arr, sec, nup-1-5*shift, nup-1-5*shift-gap);    
    MakeSeeds2(arr, sec, gap, 1);    
  }
  
  Int_t nseed=arr->GetEntriesFast();
  Int_t i;    
  
    gap=Int_t(0.3*nrows);
    // continue seeds 
    for (i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i), &t=*pt; 
    if (!pt) continue;
    if (FollowProlongation(t,nup-gap)) {
    pt->fIsSeeding =kFALSE;
    continue;
    }
    delete arr->RemoveAt(i);
    }
  
     
  //
  //remove seeds which overlaps  
  RemoveOverlap(arr,0.4,1);	  
  //delete seeds - which were sign  
  nseed=arr->GetEntriesFast();
  for (i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);
    //, &t=*pt;
    if (!pt) continue;
    if ((pt->IsActive())  &&  pt->GetNumberOfClusters() > pt->fNFoundable*0.5 ) {
      //pt->Reset();
      //FollowBackProlongation(*pt,nup-1);
      //if ( pt->GetNumberOfClusters() < pt->fNFoundable*0.5 || pt->GetNumberOfClusters()<10 )
      //delete arr->RemoveAt(i);
      //else
      //	pt->Reset();
      continue;
    }
    delete arr->RemoveAt(i);    
  } 
  //RemoveOverlap(arr,0.6,1);
  return arr;
}



//_____________________________________________________________________________
void AliTPCtrackerMI::MakeSeeds(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2) {
  //-----------------------------------------------------------------
  // This function creates track seeds.
  //-----------------------------------------------------------------
  //  if (fSeeds==0) fSeeds=new TObjArray(15000);

  Double_t x[5], c[15];

  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  Double_t cs=cos(alpha), sn=sin(alpha);

  Double_t x1 =fOuterSec->GetX(i1);
  Double_t xx2=fOuterSec->GetX(i2);
  //
  //  for (Int_t ns=0; ns<fkNOS; ns++) 
  Int_t ns =sec;
    {
    Int_t nl=fOuterSec[(ns-1+fkNOS)%fkNOS][i2];
    Int_t nm=fOuterSec[ns][i2];
    Int_t nu=fOuterSec[(ns+1)%fkNOS][i2];
    const AliTPCRow& kr1=fOuterSec[ns][i1];
    AliTPCRow&  kr21 = fOuterSec[(ns-1+fkNOS)%fkNOS][i2];
    AliTPCRow&  kr22 = fOuterSec[(ns)%fkNOS][i2];
    AliTPCRow&  kr23 = fOuterSec[(ns+1)%fkNOS][i2];

    for (Int_t is=0; is < kr1; is++) {
      Double_t y1=kr1[is]->GetY(), z1=kr1[is]->GetZ();
      Double_t x3=GetX(), y3=GetY(), z3=GetZ();

      Float_t anglez = (z1-z3)/(x1-x3); 
      Float_t extraz = z1 - anglez*(x1-xx2);  // extrapolated z

      for (Int_t js=0; js < nl+nm+nu; js++) {
	const AliTPCclusterMI *kcl;
        Double_t x2,   y2,   z2;
	if (js<nl) {	 
	  if (js==0) {
	    js = kr21.Find(extraz-15.);
	    if (js>=nl) continue;
	  }	  
	  kcl=kr21[js];
	  z2=kcl->GetZ();
	  if ((extraz-z2)>10) continue;	  
	  if ((extraz-z2)<-10) {
	    js = nl-1;
	    continue;
	  }
          y2=kcl->GetY(); 
          x2= xx2*cs+y2*sn;
          y2=-xx2*sn+y2*cs;
	} else 
	  if (js<nl+nm) {
	    if (js==nl) {
	      js = nl+kr22.Find(extraz-15.);
	      if (js>=nl+nm) continue;
	    }	  	  
	    kcl=kr22[js-nl];
	    z2=kcl->GetZ();
	    if ((extraz-z2)>10) continue;
	    if ((extraz-z2)<-10) {
	      js = nl+nm-1;
	      continue;
	    }
            x2=xx2; y2=kcl->GetY(); 
	  } else {
	    //const AliTPCRow& kr2=fOuterSec[(ns+1)%fkNOS][i2];	  
	    if (js==nl+nm) {
	      js = nl+nm+kr23.Find(extraz-15.);
	      if (js>=nl+nm+nu) break;
	    }	  
	    kcl=kr23[js-nl-nm];
	    z2=kcl->GetZ(); 
	    if ((extraz-z2)>10) continue;
	    if ((extraz-z2)<-10) {
	      break;	    
	    }
            y2=kcl->GetY();
            x2=xx2*cs-y2*sn;
            y2=xx2*sn+y2*cs;
	  }

        Double_t zz=z1 - anglez*(x1-x2); 
        if (TMath::Abs(zz-z2)>10.) continue;

        Double_t d=(x2-x1)*(0.-y2)-(0.-x2)*(y2-y1);
        if (d==0.) {cerr<<"MakeSeeds warning: Straight seed !\n"; continue;}

	x[0]=y1;
	x[1]=z1;
	x[4]=f1(x1,y1,x2,y2,x3,y3);
	if (TMath::Abs(x[4]) >= 0.0066) continue;
	x[2]=f2(x1,y1,x2,y2,x3,y3);
	//if (TMath::Abs(x[4]*x1-x[2]) >= 0.99999) continue;
	x[3]=f3(x1,y1,x2,y2,z1,z2);
	if (TMath::Abs(x[3]) > 1.2) continue;
	Double_t a=asin(x[2]);
	Double_t zv=z1 - x[3]/x[4]*(a+asin(x[4]*x1-x[2]));
	if (TMath::Abs(zv-z3)>10.) continue; 

        Double_t sy1=kr1[is]->GetSigmaY2()*2, sz1=kr1[is]->GetSigmaZ2()*4;
        Double_t sy2=kcl->GetSigmaY2()*2,     sz2=kcl->GetSigmaZ2()*4;
	//Double_t sy3=400*3./12., sy=0.1, sz=0.1;
	Double_t sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;
	//Double_t sy3=25000*x[4]*x[4]*60+0.5, sy=0.1, sz=0.1;

	Double_t f40=(f1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
	Double_t f42=(f1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
	Double_t f43=(f1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
	Double_t f20=(f2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(f2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f23=(f2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
	Double_t f30=(f3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
	Double_t f31=(f3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
	Double_t f32=(f3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
	Double_t f34=(f3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;

        c[0]=sy1;
        c[1]=0.;       c[2]=sz1;
        c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
        c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
                       c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
        c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
        c[13]=f30*sy1*f40+f32*sy2*f42;
        c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;

        UInt_t index=kr1.GetIndex(is);
	AliTPCseed *track=new AliTPCseed(index, x, c, x1, ns*alpha+shift);
	track->fIsSeeding = kTRUE;
	//Int_t rc=FollowProlongation(*track, i2);
	Int_t delta4 = Int_t((i2-i1)/4.);

	FollowProlongation(*track, i1-delta4);
	if (track->GetNumberOfClusters() < track->fNFoundable/2.) {
	  delete track;
	  continue;
	}
	FollowProlongation(*track, i1-2*delta4);
	if (track->GetNumberOfClusters() < track->fNFoundable/2.) {
	  delete track;
	  continue;
	}
	FollowProlongation(*track, i1-3*delta4);
	if (track->GetNumberOfClusters() < track->fNFoundable/2.) {
	  delete track;
	  continue;
	}
	FollowProlongation(*track, i2);
	//Int_t rc = 1;
	
	track->fLastPoint = i1;  // first cluster in track position
	if (track->GetNumberOfClusters()<(i1-i2)/4 || track->GetNumberOfClusters() < track->fNFoundable/2. ) delete track;
        else arr->AddLast(track); 
      }
    }
  }
}


//_____________________________________________________________________________
void AliTPCtrackerMI::MakeSeeds2(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2) {
  //-----------------------------------------------------------------
  // This function creates track seeds - without vertex constraint
  //-----------------------------------------------------------------

  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  //  Double_t cs=cos(alpha), sn=sin(alpha);
  Int_t row0 = (i1+i2)/2;
  Int_t drow = (i1-i2)/2;
  const AliTPCRow& kr0=fSectors[sec][row0];
  const AliTPCRow& krm=fSectors[sec][row0-1];
  const AliTPCRow& krp=fSectors[sec][row0+1];
  AliTPCRow * kr=0;

  AliTPCpolyTrack polytrack;
  Int_t nclusters=fSectors[sec][row0];

  for (Int_t is=0; is < nclusters; is++) {
    const AliTPCclusterMI * cl= kr0[is];
    Double_t x = kr0.GetX();

    // Initialization of the polytrack
    polytrack.Reset();

    Double_t y0= cl->GetY();
    Double_t z0= cl->GetZ();
    polytrack.AddPoint(x,y0,z0);
    Float_t roady = 5*TMath::Sqrt(cl->GetSigmaY2()+0.2);
    Float_t roadz = 5*TMath::Sqrt(cl->GetSigmaZ2()+0.2);
    //
    x = krm.GetX();
    cl = krm.FindNearest(y0,z0,roady,roadz);
    if (cl) polytrack.AddPoint(x,cl->GetY(),cl->GetZ(),roady,roadz);
    //
    x = krp.GetX();
    cl = krp.FindNearest(y0,z0,roady,roadz);
    if (cl) polytrack.AddPoint(x,cl->GetY(),cl->GetZ(),cl->GetSigmaY2()+0.05,cl->GetSigmaZ2()+0.05);
    //
    polytrack.UpdateParameters();
    // follow polytrack
    roadz = 0.6;
    roady = 0.6;
    //
    Double_t yn,zn;
    Int_t nfoundable = polytrack.GetN();
    Int_t nfound     = nfoundable; 
    for (Int_t ddrow = 2; ddrow<drow;ddrow++){
      for (Int_t delta = -1;delta<=1;delta+=2){
	Int_t row = row0+ddrow*delta;
	kr = &(fSectors[sec][row]);
	Double_t xn = kr->GetX();
	Double_t ymax = fSectors->GetMaxY(row)-kr->fDeadZone;
	polytrack.GetFitPoint(xn,yn,zn);
	if (TMath::Abs(yn)>ymax) continue;
	nfoundable++;
	AliTPCclusterMI * cln = kr->FindNearest(yn,zn,roady,roadz);
	if (cln) {
	  polytrack.AddPoint(xn,cln->GetY(),cln->GetZ(),cln->GetSigmaY2()+0.05,cln->GetSigmaZ2()+0.05);
	  nfound++;
	}
      }
      polytrack.UpdateParameters();
      if (nfound<0.45*nfoundable) break;
    }
    if ((nfound>0.5*nfoundable) &&( nfoundable>0.4*(i1-i2))) {
      // add polytrack candidate
      Double_t x[5], c[15];
      Double_t x1,x2,x3,y1,y2,y3,z1,z2,z3;
      polytrack.GetBoundaries(x3,x1);
      x2 = (x1+x3)/2.;
      polytrack.GetFitPoint(x1,y1,z1);
      polytrack.GetFitPoint(x2,y2,z2);
      polytrack.GetFitPoint(x3,y3,z3);
      //
      //is track pointing to the vertex ?
      Double_t x0,y0,z0;
      x0=0;
      polytrack.GetFitPoint(x0,y0,z0);
      if ( (TMath::Abs(z0-GetZ())<10) && (TMath::Abs(y0-GetY())<5)){ //if yes apply vertex constraint
	//	x3 = 0;
	//y3 = GetY();
	//z3 = GetZ();
      }

      x[0]=y1;
      x[1]=z1;
      x[4]=f1(x1,y1,x2,y2,x3,y3);
      if (TMath::Abs(x[4]) >= 0.05) continue;  //MI change
      x[2]=f2(x1,y1,x2,y2,x3,y3);
      //if (TMath::Abs(x[4]*x1-x[2]) >= 0.99999) continue;
      x[3]=f3(x1,y1,x2,y2,z1,z2);
      if (TMath::Abs(x[3]) > 1.2) continue;
      if (TMath::Abs(x[2]) > 0.99) continue;
      //      Double_t a=asin(x[2]);
      

      Double_t sy=1.5, sz=1.5;
      Double_t sy1=1.5, sz1=1.5;
      Double_t sy2=1.3, sz2=1.3;
      Double_t sy3=1.5;
      //sz3=1.5;
      //Double_t sy3=400*3./12., sy=0.1, sz=0.1;
      //      Double_t sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;
      //Double_t sy3=25000*x[4]*x[4]*60+0.5, sy=0.1, sz=0.1;
      
      Double_t f40=(f1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
      Double_t f42=(f1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
      Double_t f43=(f1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
      Double_t f20=(f2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
      Double_t f22=(f2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
      Double_t f23=(f2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
      Double_t f30=(f3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
      Double_t f31=(f3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
      Double_t f32=(f3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
      Double_t f34=(f3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;
      
      c[0]=sy1;
      c[1]=0.;       c[2]=sz1;
      c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
      c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
      c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
      c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
      c[13]=f30*sy1*f40+f32*sy2*f42;
      c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
      
      UInt_t index=0;
      //kr0.GetIndex(is);
      AliTPCseed *track=new AliTPCseed(index, x, c, x1, sec*alpha+shift);
      track->fStopped =kFALSE;
      track->fIsSeeding = kTRUE;
      Int_t rc=FollowProlongation(*track, i2);	
	track->fLastPoint = i1;  // first cluster in track position
	if (rc==0 || track->GetNumberOfClusters()<(i1-i2)/4 || track->GetNumberOfClusters() < track->fNFoundable/2. ) delete track;
        else arr->AddLast(track); 
    }
  }
  
  
  
}





//_____________________________________________________________________________
Int_t AliTPCtrackerMI::ReadSeeds(const TFile *inp) {
  //-----------------------------------------------------------------
  // This function reades track seeds.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 

  TFile *in=(TFile*)inp;
  if (!in->IsOpen()) {
     cerr<<"AliTPCtrackerMI::ReadSeeds(): input file is not open !\n";
     return 1;
  }

  in->cd();
  TTree *seedTree=(TTree*)in->Get("Seeds");
  if (!seedTree) {
     cerr<<"AliTPCtrackerMI::ReadSeeds(): ";
     cerr<<"can't get a tree with track seeds !\n";
     return 2;
  }
  AliTPCtrack *seed=new AliTPCtrack; 
  seedTree->SetBranchAddress("tracks",&seed);
  
  if (fSeeds==0) fSeeds=new TObjArray(15000);

  Int_t n=(Int_t)seedTree->GetEntries();
  for (Int_t i=0; i<n; i++) {
     seedTree->GetEvent(i);
     fSeeds->AddLast(new AliTPCseed(*seed,seed->GetAlpha()));
  }
  
  delete seed;
  delete seedTree; 
  savedir->cd();
  return 0;
}

//_____________________________________________________________________________
Int_t AliTPCtrackerMI::Clusters2Tracks() {
  //-----------------------------------------------------------------
  // This is a track finder.
  //-----------------------------------------------------------------
  TTree* clustersTree = AliRunLoader::GetTreeR("TPC", kFALSE,AliConfig::fgkDefaultEventFolderName);
  if (!clustersTree) {
    Error("Clusters2Tracks", "no clusters found");
    return 1;
  }
  fClustersArray.ConnectTree(clustersTree);

  TTree* tracksTree = AliRunLoader::GetTreeT("TPC", kTRUE,AliConfig::fgkDefaultEventFolderName);
  TTree& tracktree = *tracksTree;
//  TTree seedtree("Seeds","Seeds");
  AliTPCtrack *iotrack=0;
  AliTPCseed  *ioseed=0;
  tracktree.Branch("tracks","AliTPCtrack",&iotrack,32000,0);
  TStopwatch timer;
  
  printf("Loading clusters \n");
  LoadClusters();
  printf("Time for loading clusters: \t");timer.Print();timer.Start();

  fSectors = fOuterSec;
  fN=fkNOS;
  
  
  //find track seeds
  MakeSeedsAll(); 
  printf("Time for seeding: \t"); timer.Print();timer.Start();
  Int_t nup=fOuterSec->GetNRows(), nlow=fInnerSec->GetNRows();
  Int_t nrows=nlow+nup;
  
  Int_t gap=Int_t(0.3*nrows);
  Int_t i;
  //RemoveOverlap(fSeeds,0.6,2);
  Int_t nseed=fSeeds->GetEntriesFast();
  // outer sectors parallel tracking
  ParallelTracking(fSectors->GetNRows()-gap-1,0); 
  printf("Time for parralel tracking outer sectors: \t"); timer.Print();timer.Start();

  RemoveOverlap(fSeeds, 0.4,3);	 
  printf("Time for removal overlap- outer sectors: \t");timer.Print();timer.Start();
  //parallel tracking 
  fSectors = fInnerSec;
  fN=fkNIS;  

  ParallelTracking(fSectors->GetNRows()-1,0);
  /*
    ParallelTracking(fSectors->GetNRows()-1,2*fSectors->GetNRows()/3);
    RemoveOverlap(fSeeds,0.4,5,kTRUE);
    ParallelTracking(2*fSectors->GetNRows()/3-1,fSectors->GetNRows()/3);
    RemoveOverlap(fSeeds,0.4,5,kTRUE);
    ParallelTracking(fSectors->GetNRows()/3-1,0);
  */
  printf("Number of tracks after  inner tracking  %d\n",fNtracks); 
  printf("Time for parralel tracking inner sectors: \t"); timer.Print();timer.Start();
  //
  for (Int_t i=0;i<fSeeds->GetEntriesFast();i++){
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i);    
    if (!pt) continue;   
    if (!pt->IsActive()) continue;
    pt->PropagateTo(90.);
  } 
  RemoveOverlap(fSeeds,0.4,5,kTRUE);  // remove overlap -  shared points signed 
  RemoveUsed(fSeeds,0.4,6);
  printf("Time for removal overlap- inner sectors: \t"); timer.Print();timer.Start();
  //
  // 
  
  ioseed  = (AliTPCseed*)(fSeeds->UncheckedAt(0));
  AliTPCseed * vseed = new AliTPCseed;
  vseed->fPoints = new TClonesArray("AliTPCTrackPoint",1);
  vseed->fEPoints = new TClonesArray("AliTPCExactPoint",1);
  vseed->fPoints->ExpandCreateFast(2);
  
  //TBranch * seedbranch =   
//  seedtree.Branch("seeds","AliTPCseed",&vseed,32000,99);
  //delete vseed;
  nseed=fSeeds->GetEntriesFast();

  Int_t found = 0;
  for (i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<20) continue;
    t.CookdEdx(0.02,0.6);
    CookLabel(pt,0.1); //For comparison only
    //   if ((pt->IsActive()) && (nc>Int_t(0.4*nrows))){
    if ((pt->IsActive()) && (nc>Int_t(0.5*t.fNFoundable) && (t.fNFoundable>Int_t(0.3*nrows)))){
      iotrack=pt;
      tracktree.Fill();
//     cerr<<found++<<'\r';      
    }   
    /*
      pt->RebuildSeed();
      seedbranch->SetAddress(&pt);
      
      seedtree.Fill();        
      for (Int_t j=0;j<160;j++){
      delete pt->fPoints->RemoveAt(j);
      }
      delete pt->fPoints;
      pt->fPoints =0;
    */
    delete fSeeds->RemoveAt(i);
  }
  //  fNTracks = found;
  printf("Time for track writing and dedx cooking: \t"); timer.Print();timer.Start();

  UnloadClusters();
  printf("Time for unloading cluster: \t"); timer.Print();timer.Start();

//  seedtree.Write();
  cerr<<"Number of found tracks : "<<"\t"<<found<<endl;
  
  AliRunLoader::GetDetectorLoader("TPC",AliConfig::fgkDefaultEventFolderName)->WriteTracks("OVERWRITE");

  return 0;
}


void  AliTPCtrackerMI::ParallelTracking(Int_t rfirst, Int_t rlast)
{
  //
  // try to track in parralel

  Int_t nseed=fSeeds->GetEntriesFast();
  //prepare seeds for tracking
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt; 
    if (!pt) continue;
    if (!t.IsActive()) continue;
    // follow prolongation to the first layer
    if ( (fSectors ==fInnerSec) || (t.fFirstPoint>rfirst+1))  
      FollowProlongation(t, rfirst+1);
  }


  //
  for (Int_t nr=rfirst; nr>=rlast; nr--){      
    // make indexes with the cluster tracks for given       
    //    for (Int_t i = 0;i<fN;i++)
    //  fSectors[i][nr].MakeClusterTracks();

    // find nearest cluster
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;       
      if (!pt) continue;
      if (!pt->IsActive()) continue;
      if ( (fSectors ==fOuterSec) && pt->fFirstPoint<nr) continue;
      if (pt->fRelativeSector>17) {
	continue;
      }
      UpdateClusters(t,i,nr);
    }
    // prolonagate to the nearest cluster - if founded
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i); 
      if (!pt) continue;
      if (!pt->IsActive()) continue; 
      if ((fSectors ==fOuterSec) &&pt->fFirstPoint<nr) continue;
      if (pt->fRelativeSector>17) {
	continue;
      }
      FollowToNextCluster(i,nr);
    }
    //    for (Int_t i= 0;i<fN;i++)
    //  fSectors[i][nr].ClearClusterTracks();
  }  
  
}

Float_t  AliTPCtrackerMI::GetSigmaY(AliTPCseed * seed)
{
  //
  //  
  Float_t sd2 = TMath::Abs((fParam->GetZLength()-TMath::Abs(seed->GetZ())))*fParam->GetDiffL()*fParam->GetDiffL();
  Float_t padlength =  fParam->GetPadPitchLength(seed->fSector);
  Float_t sres = (seed->fSector < fParam->GetNSector()/2) ? 0.2 :0.3;
  Float_t angular  = seed->GetSnp();
  angular = angular*angular/(1-angular*angular);
  //  angular*=angular;
  //angular  = TMath::Sqrt(angular/(1-angular));
  Float_t res = TMath::Sqrt(sd2+padlength*padlength*angular/12.+sres*sres);
  return res;
}
Float_t  AliTPCtrackerMI::GetSigmaZ(AliTPCseed * seed)
{
  //
  //
  Float_t sd2 = TMath::Abs((fParam->GetZLength()-TMath::Abs(seed->GetZ())))*fParam->GetDiffL()*fParam->GetDiffL();
  Float_t padlength =  fParam->GetPadPitchLength(seed->fSector);
  Float_t sres = fParam->GetZSigma();
  Float_t angular  = seed->GetTgl();
  Float_t res = TMath::Sqrt(sd2+padlength*padlength*angular*angular/12.+sres*sres);
  return res;
}



//_________________________________________________________________________
AliTPCclusterMI *AliTPCtrackerMI::GetClusterMI(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t sec=(index&0xff000000)>>24; 
  Int_t row=(index&0x00ff0000)>>16; 
  Int_t ncl=(index&0x0000ffff)>>00;

  AliTPCClustersRow *clrow=((AliTPCtrackerMI *) this)->fClustersArray.GetRow(sec,row);
  if (!clrow) return 0;
  return (AliTPCclusterMI*)(*clrow)[ncl];      
}

//__________________________________________________________________________
void AliTPCtrackerMI::CookLabel(AliKalmanTrack *t, Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  Int_t noc=t->GetNumberOfClusters();
  Int_t *lb=new Int_t[noc];
  Int_t *mx=new Int_t[noc];
  AliTPCclusterMI **clusters=new AliTPCclusterMI*[noc];

  Int_t i;
  for (i=0; i<noc; i++) {
     lb[i]=mx[i]=0;
     Int_t index=t->GetClusterIndex(i);
     clusters[i]=GetClusterMI(index);
  }

  Int_t lab=123456789;
  for (i=0; i<noc; i++) {
    AliTPCclusterMI *c=clusters[i];
    if (!clusters[i]) continue;
    lab=TMath::Abs(c->GetLabel(0));
    Int_t j;
    for (j=0; j<noc; j++) if (lb[j]==lab || mx[j]==0) break;
    lb[j]=lab;
    (mx[j])++;
  }

  Int_t max=0;
  for (i=0; i<noc; i++) if (mx[i]>max) {max=mx[i]; lab=lb[i];}
    
  for (i=0; i<noc; i++) {
    AliTPCclusterMI *c=clusters[i]; 
    if (!clusters[i]) continue;
    if (TMath::Abs(c->GetLabel(1)) == lab ||
        TMath::Abs(c->GetLabel(2)) == lab ) max++;
  }

  if ((1.- Float_t(max)/noc) > wrong) lab=-lab;

  else {
     Int_t tail=Int_t(0.10*noc);
     max=0;
     for (i=1; i<=tail; i++) {
       AliTPCclusterMI *c=clusters[noc-i];
       if (!clusters[i]) continue;
       if (lab == TMath::Abs(c->GetLabel(0)) ||
           lab == TMath::Abs(c->GetLabel(1)) ||
           lab == TMath::Abs(c->GetLabel(2))) max++;
     }
     if (max < Int_t(0.5*tail)) lab=-lab;
  }

  t->SetLabel(lab);

  delete[] lb;
  delete[] mx;
  delete[] clusters;
}

//_________________________________________________________________________
void AliTPCtrackerMI::AliTPCSector::Setup(const AliTPCParam *par, Int_t f) {
  //-----------------------------------------------------------------------
  // Setup inner sector
  //-----------------------------------------------------------------------
  if (f==0) {
     fAlpha=par->GetInnerAngle();
     fAlphaShift=par->GetInnerAngleShift();
     fPadPitchWidth=par->GetInnerPadPitchWidth();
     fPadPitchLength=par->GetInnerPadPitchLength();
     fN=par->GetNRowLow();
     fRow=new AliTPCRow[fN];
     for (Int_t i=0; i<fN; i++) {
       fRow[i].SetX(par->GetPadRowRadiiLow(i));
       fRow[i].fDeadZone =1.5;  //1.5 cm of dead zone
     }
  } else {
     fAlpha=par->GetOuterAngle();
     fAlphaShift=par->GetOuterAngleShift();
     fPadPitchWidth  = par->GetOuterPadPitchWidth();
     fPadPitchLength = par->GetOuter1PadPitchLength();
     f1PadPitchLength = par->GetOuter1PadPitchLength();
     f2PadPitchLength = par->GetOuter2PadPitchLength();

     fN=par->GetNRowUp();
     fRow=new AliTPCRow[fN];
     for (Int_t i=0; i<fN; i++) {
       fRow[i].SetX(par->GetPadRowRadiiUp(i)); 
       fRow[i].fDeadZone =1.5;  // 1.5 cm of dead zone
     }
  } 
}


AliTPCtrackerMI::AliTPCRow::~AliTPCRow(){
  //
  if (fClusterTracks) delete [] fClusterTracks;
  fClusterTracks = 0;
}

void AliTPCtrackerMI::AliTPCRow::MakeClusterTracks(){
  //create cluster tracks
  if (fN>0) 
    fClusterTracks = new AliTPCclusterTracks[fN];
}

void AliTPCtrackerMI::AliTPCRow::ClearClusterTracks(){
  if (fClusterTracks) delete[] fClusterTracks;
  fClusterTracks =0;
}



void AliTPCtrackerMI::AliTPCRow::UpdateClusterTrack(Int_t clindex, Int_t trindex, AliTPCseed * seed){
  //
  //
  // update information of the cluster tracks - if track is nearer then other tracks to the 
  // given track
  const AliTPCclusterMI * cl = (*this)[clindex];
  AliTPCclusterTracks * cltracks = GetClusterTracks(clindex);
  // find the distance of the cluster to the track
  Float_t dy2 = (cl->GetY()- seed->GetY());
  dy2*=dy2;
  Float_t dz2 = (cl->GetZ()- seed->GetZ());
  dz2*=dz2;
  //
  Float_t distance = TMath::Sqrt(dy2+dz2);
  if (distance> 3) 
    return;  // MI - to be changed - AliTPCtrackerParam
  
  if ( distance < cltracks->fDistance[0]){
    cltracks->fDistance[2] =cltracks->fDistance[1];
    cltracks->fDistance[1] =cltracks->fDistance[0];
    cltracks->fDistance[0] =distance;
    cltracks->fTrackIndex[2] =cltracks->fTrackIndex[1];
    cltracks->fTrackIndex[1] =cltracks->fTrackIndex[0];
    cltracks->fTrackIndex[0] =trindex; 
  }
  else
    if ( distance < cltracks->fDistance[1]){
      cltracks->fDistance[2] =cltracks->fDistance[1];  
      cltracks->fDistance[1] =distance;
      cltracks->fTrackIndex[2] =cltracks->fTrackIndex[1];
      cltracks->fTrackIndex[1] =trindex; 
    } else
      if (distance < cltracks->fDistance[2]){
	cltracks->fDistance[2] =distance;
	cltracks->fTrackIndex[2] =trindex;
      }  
}


//_________________________________________________________________________
void 
AliTPCtrackerMI::AliTPCRow::InsertCluster(const AliTPCclusterMI* c, UInt_t index) {
  //-----------------------------------------------------------------------
  // Insert a cluster into this pad row in accordence with its y-coordinate
  //-----------------------------------------------------------------------
  if (fN==kMaxClusterPerRow) {
    cerr<<"AliTPCRow::InsertCluster(): Too many clusters !\n"; return;
  }
  if (fN==0) {fIndex[0]=index; fClusters[fN++]=c; return;}
  Int_t i=Find(c->GetZ());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliTPCclusterMI*));
  memmove(fIndex   +i+1 ,fIndex   +i,(fN-i)*sizeof(UInt_t));
  fIndex[i]=index; fClusters[i]=c; fN++;
}

//___________________________________________________________________
Int_t AliTPCtrackerMI::AliTPCRow::Find(Double_t z) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster 
  //-----------------------------------------------------------------------
  if (fN==0) return 0;
  if (z <= fClusters[0]->GetZ()) return 0;
  if (z > fClusters[fN-1]->GetZ()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fClusters[m]->GetZ()) b=m+1;
    else e=m; 
  }
  return m;
}



//___________________________________________________________________
AliTPCclusterMI * AliTPCtrackerMI::AliTPCRow::FindNearest(Double_t y, Double_t z, Double_t roady, Double_t roadz) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance = roady*roady + roadz*roadz;

  AliTPCclusterMI *cl =0;
  for (Int_t i=Find(z-roadz); i<fN; i++) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;
      if ( (c->GetY()-y) >  roady ) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)+(c->GetY()-y)*(c->GetY()-y);
      if (maxdistance>distance) {
	maxdistance = distance;
	cl=c;       
      }
  }
  return cl;      
}





AliTPCseed::AliTPCseed():AliTPCtrack(){
  //
  fRow=0; 
  fRemoval =0; 
  memset(fClusterIndex,0,sizeof(Int_t)*200);
  fPoints = 0;
  fEPoints = 0;
  fNFoundable =0;
  fNShared  =0;
  fTrackPoints =0;
  fRemoval = 0;
}

AliTPCseed::AliTPCseed(const AliTPCtrack &t):AliTPCtrack(t){
  fPoints = 0;
  fEPoints = 0;
  fNShared  =0; 
  fTrackPoints =0;
  fRemoval =0;
}

AliTPCseed::AliTPCseed(const AliKalmanTrack &t, Double_t a):AliTPCtrack(t,a){
  fRow=0;
  memset(fClusterIndex,0,sizeof(Int_t)*200); 
  fPoints = 0;
  fEPoints = 0;
  fNFoundable =0; 
  fNShared  =0; 
  fTrackPoints =0;
  fRemoval =0;
}

AliTPCseed::AliTPCseed(UInt_t index, const Double_t xx[5], const Double_t cc[15], 
					Double_t xr, Double_t alpha):      
  AliTPCtrack(index, xx, cc, xr, alpha) {
  //
  //
  fRow =0;
  memset(fClusterIndex,0,sizeof(Int_t)*200); 
  fPoints = 0;
  fEPoints = 0;
  fNFoundable =0;
  fNShared  = 0;
  fTrackPoints =0;
  fRemoval =0;
}

AliTPCseed::~AliTPCseed(){
  if (fPoints) delete fPoints;
  fPoints =0;
  fEPoints = 0;
  if (fTrackPoints){
    for (Int_t i=0;i<8;i++){
      delete [] fTrackPoints[i];
    }
    delete fTrackPoints;
    fTrackPoints =0;
  }

}

AliTPCTrackPoint * AliTPCseed::GetTrackPoint(Int_t i)
{
  //
  // 
  if (!fTrackPoints) {
    fTrackPoints = new AliTPCTrackPoint*[8];
    for ( Int_t i=0;i<8;i++)
      fTrackPoints[i]=0;
  }
  Int_t index1 = i/20;
  if (!fTrackPoints[index1]) fTrackPoints[index1] = new AliTPCTrackPoint[20];
  return &(fTrackPoints[index1][i%20]);
}

void AliTPCseed::RebuildSeed()
{
  //
  // rebuild seed to be ready for storing
  fPoints = new TClonesArray("AliTPCTrackPoint",160);
  fPoints->ExpandCreateFast(160);
  fEPoints = new TClonesArray("AliTPCExactPoint",1);
  for (Int_t i=0;i<160;i++){
    AliTPCTrackPoint *trpoint = (AliTPCTrackPoint*)fPoints->UncheckedAt(i);
    *trpoint = *(GetTrackPoint(i));
  }

}

//_____________________________________________________________________________
void AliTPCseed::CookdEdx(Double_t low, Double_t up) {
  //-----------------------------------------------------------------
  // This funtion calculates dE/dX within the "low" and "up" cuts.
  //-----------------------------------------------------------------

  Float_t amp[200];
  Float_t angular[200];
  Float_t weight[200];
  Int_t index[200];
  //Int_t nc = 0;
  //  TClonesArray & arr = *fPoints; 
  Float_t meanlog = 100.;
  
  Float_t mean[4]  = {0,0,0,0};
  Float_t sigma[4] = {1000,1000,1000,1000};
  Int_t nc[4]      = {0,0,0,0};
  Float_t norm[4]    = {1000,1000,1000,1000};
  //
  //
  fNShared =0;

  for (Int_t of =0; of<4; of++){    
    for (Int_t i=of;i<160;i+=4)
      {
	//AliTPCTrackPoint * point = (AliTPCTrackPoint *) arr.At(i);
	AliTPCTrackPoint * point = GetTrackPoint(i);
	if (point==0) continue;
	if (point->fIsShared){
	  fNShared++;
	  continue;
	}
	if (point->GetCPoint().GetMax()<5) continue;
	Float_t angley = point->GetTPoint().GetAngleY();
	Float_t anglez = point->GetTPoint().GetAngleZ();
	Int_t   type   = point->GetCPoint().GetType();
	Float_t rsigmay =  point->GetCPoint().GetSigmaY();
	Float_t rsigmaz =  point->GetCPoint().GetSigmaZ();
	Float_t rsigma = TMath::Sqrt(rsigmay*rsigmaz);

	Float_t ampc   = 0;     // normalization to the number of electrons
	if (i>64){
	  ampc = 1.*point->GetCPoint().GetMax();
	  //ampc = 1.*point->GetCPoint().GetQ();	  
	  //	  AliTPCClusterPoint & p = point->GetCPoint();
	  //	  Float_t dy = TMath::Abs(Int_t( TMath::Abs(p.GetY()/0.6)) - TMath::Abs(p.GetY()/0.6)+0.5);
	  // Float_t iz =  (250.0-TMath::Abs(p.GetZ())+0.11)/0.566;
	  //Float_t dz = 
	  //  TMath::Abs( Int_t(iz) - iz + 0.5);
	  //ampc *= 1.15*(1-0.3*dy);
	  //ampc *= 1.15*(1-0.3*dz);
	  //	  Float_t zfactor = (1.05-0.0004*TMath::Abs(point->GetCPoint().GetZ()));
	  //ampc               *=zfactor; 
	}
	else{ 
	  ampc = 1.0*point->GetCPoint().GetMax(); 
	  //ampc = 1.0*point->GetCPoint().GetQ(); 
	  //AliTPCClusterPoint & p = point->GetCPoint();
	  // Float_t dy = TMath::Abs(Int_t( TMath::Abs(p.GetY()/0.4)) - TMath::Abs(p.GetY()/0.4)+0.5);
	  //Float_t iz =  (250.0-TMath::Abs(p.GetZ())+0.11)/0.566;
	  //Float_t dz = 
	  //  TMath::Abs( Int_t(iz) - iz + 0.5);

	  //ampc *= 1.15*(1-0.3*dy);
	  //ampc *= 1.15*(1-0.3*dz);
	  //	Float_t zfactor = (1.02-0.000*TMath::Abs(point->GetCPoint().GetZ()));
	  //ampc               *=zfactor; 

	}
	ampc *= 2.0;     // put mean value to channel 50
	//ampc *= 0.58;     // put mean value to channel 50
	Float_t w      =  1.;
	//	if (type>0)  w =  1./(type/2.-0.5); 
	//	Float_t z = TMath::Abs(point->GetCPoint().GetZ());
	if (i<64) {
	  ampc /= 0.6;
	  //ampc /= (1+0.0008*z);
	} else
	  if (i>128){
	    ampc /=1.5;
	    //ampc /= (1+0.0008*z);
	  }else{
	    //ampc /= (1+0.0008*z);
	  }
	
	if (type<0) {  //amp at the border - lower weight
	  // w*= 2.;
	  
	  continue;
	}
	if (rsigma>1.5) ampc/=1.3;  // if big backround
	amp[nc[of]]        = ampc;
	angular[nc[of]]    = TMath::Sqrt(1.+angley*angley+anglez*anglez);
	weight[nc[of]]     = w;
	nc[of]++;
      }
    
    TMath::Sort(nc[of],amp,index,kFALSE);
    Float_t sumamp=0;
    Float_t sumamp2=0;
    Float_t sumw=0;
    //meanlog = amp[index[Int_t(nc[of]*0.33)]];
    meanlog = 200;
    for (Int_t i=int(nc[of]*low+0.5);i<int(nc[of]*up+0.5);i++){
      Float_t ampl      = amp[index[i]]/angular[index[i]];
      ampl              = meanlog*TMath::Log(1.+ampl/meanlog);
      //
      sumw    += weight[index[i]]; 
      sumamp  += weight[index[i]]*ampl;
      sumamp2 += weight[index[i]]*ampl*ampl;
      norm[of]    += angular[index[i]]*weight[index[i]];
    }
    if (sumw<1){ 
      SetdEdx(0);  
    }
    else {
      norm[of] /= sumw;
      mean[of]  = sumamp/sumw;
      sigma[of] = sumamp2/sumw-mean[of]*mean[of];
      if (sigma[of]>0.1) 
	sigma[of] = TMath::Sqrt(sigma[of]);
      else
	sigma[of] = 1000;
      
    mean[of] = (TMath::Exp(mean[of]/meanlog)-1)*meanlog;
    //mean  *=(1-0.02*(sigma/(mean*0.17)-1.));
    //mean *=(1-0.1*(norm-1.));
    }
  }

  Float_t dedx =0;
  fSdEdx =0;
  fMAngular =0;
  //  mean[0]*= (1-0.05*(sigma[0]/(0.01+mean[1]*0.18)-1));
  //  mean[1]*= (1-0.05*(sigma[1]/(0.01+mean[0]*0.18)-1));

  
  //  dedx = (mean[0]* TMath::Sqrt((1.+nc[0]))+ mean[1]* TMath::Sqrt((1.+nc[1])) )/ 
  //  (  TMath::Sqrt((1.+nc[0]))+TMath::Sqrt((1.+nc[1])));

  Int_t norm2 = 0;
  Int_t norm3 = 0;
  for (Int_t i =0;i<4;i++){
    if (nc[i]>2&&nc[i]<1000){
      dedx      += mean[i] *nc[i];
      fSdEdx    += sigma[i]*(nc[i]-2);
      fMAngular += norm[i] *nc[i];    
      norm2     += nc[i];
      norm3     += nc[i]-2;
    }
    fDEDX[i]  = mean[i];             
    fSDEDX[i] = sigma[i];            
    fNCDEDX[i]= nc[i]; 
  }

  if (norm3>0){
    dedx   /=norm2;
    fSdEdx /=norm3;
    fMAngular/=norm2;
  }
  else{
    SetdEdx(0);
    return;
  }
  //  Float_t dedx1 =dedx;
  
  dedx =0;
  for (Int_t i =0;i<4;i++){
    if (nc[i]>2&&nc[i]<1000){
      mean[i]   = mean[i]*(1-0.12*(sigma[i]/(fSdEdx)-1.));
      dedx      += mean[i] *nc[i];
    }
    fDEDX[i]  = mean[i];                
  }
  dedx /= norm2;
  

  
  SetdEdx(dedx);
    
  //mi deDX



  //Very rough PID
  Double_t p=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt()));

  if (p<0.6) {
    if (dedx < 39.+ 12./(p+0.25)/(p+0.25)) { SetMass(0.13957); return;}
    if (dedx < 39.+ 12./p/p) { SetMass(0.49368); return;}
    SetMass(0.93827); return;
  }

  if (p<1.2) {
    if (dedx < 39.+ 12./(p+0.25)/(p+0.25)) { SetMass(0.13957); return;}
    SetMass(0.93827); return;
  }

  SetMass(0.13957); return;

}


