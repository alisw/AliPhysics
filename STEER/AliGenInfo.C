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

Macro to generate comples MC information - used for Comparison later on
How to use it?

.L $ALICE_ROOT/STEER/AliGenInfo.C+
AliGenInfoMaker *t = new AliGenInfoMaker("galice.root","genTracks.root",0,0)
t->Exec();

*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string.h>
//ROOT includes
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TString.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TView.h"
#include "TGeometry.h"
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
#include "AliMathBase.h"

#endif
#include "AliGenInfo.h" 
//
// 

AliTPCParam * GetTPCParam(){
  AliTPCParamSR * par = new AliTPCParamSR;
  par->Update();
  return par;
}

AliPointsMI::AliPointsMI(){
  fN=0;
  fX=0;
  fY=0;
  fZ=0;  
  fCapacity = 0;
  fLabel0=0;
  fLabel1=0;
}


AliPointsMI::AliPointsMI(Int_t n, Float_t *x,Float_t *y, Float_t *z){
  //
  //
  fN=n;
  fCapacity = 1000+2*n;
  fX= new Float_t[n];
  fY= new Float_t[n];
  fZ= new Float_t[n];
  memcpy(fX,x,n*sizeof(Float_t));
  memcpy(fY,y,n*sizeof(Float_t));
  memcpy(fZ,z,n*sizeof(Float_t));
  fLabel0=0;
  fLabel1=0;
}

void AliPointsMI::Reset()
{
  fN=0;
}

void AliPointsMI::Reset(AliDetector * det, Int_t particle){
  //
  // write points from detector points
  //
  Reset();
  if (!det) return;
  TObjArray *array = det->Points();
  if (!array) return;
  for (Int_t i=0;i<array->GetEntriesFast();i++){
    AliPoints * points = (AliPoints*) array->At(i);
    if (!points) continue;
    if (points->GetIndex()!=particle) continue;
    Int_t npoints = points->GetN();
    if (npoints<2) continue;
    Int_t delta = npoints/100;
    if (delta<1) delta=1;
    if (delta>10) delta=10;
    Int_t mypoints = npoints/delta;
    //
    fN = mypoints;
    if (fN>fCapacity){
      fCapacity = 1000+2*fN;
      delete []fX;
      delete []fY;
      delete []fZ;
      fX = new Float_t[fCapacity];
      fY = new Float_t[fCapacity];
      fZ = new Float_t[fCapacity]; 
    }   
    Float_t *xyz = points->GetP();
    for (Int_t ipoint=0;ipoint<mypoints;ipoint++){
      Int_t index = 3*ipoint*delta;
      fX[ipoint]=0;
      fY[ipoint]=0;
      fZ[ipoint]=0;
      if (index+2<npoints*3){
	fX[ipoint] = xyz[index];
	fY[ipoint] = xyz[index+1];
	fZ[ipoint] = xyz[index+2];
      }
    }
  }
  fLabel0 = particle;
}


AliPointsMI::~AliPointsMI(){
  fN=0;
  fCapacity =0;
  delete[] fX;
  delete[]fY;
  delete []fZ;
  fX=0;fY=0;fZ=0;  
}





////////////////////////////////////////////////////////////////////////
AliMCInfo::AliMCInfo()
{
  fTPCReferences  = new TClonesArray("AliTrackReference",10);
  fITSReferences  = new TClonesArray("AliTrackReference",10);
  fTRDReferences  = new TClonesArray("AliTrackReference",10);
  fTOFReferences  = new TClonesArray("AliTrackReference",10);
  fTRdecay.SetTrack(-1);
  fCharge = 0;
}

AliMCInfo::~AliMCInfo()
{
  if (fTPCReferences) {
    delete fTPCReferences;
  }
  if (fITSReferences){
    delete fITSReferences;
  }
  if (fTRDReferences){    
    delete fTRDReferences;  
  }
  if (fTOFReferences){    
    delete fTOFReferences;  
  }

}



void AliMCInfo::Update()
{
  //
  //
  fMCtracks =1;
  if (!fTPCReferences) {
    fNTPCRef =0;
    return;
  }
  Float_t direction=1;
  //Float_t rlast=0;
  fNTPCRef = fTPCReferences->GetEntriesFast();
  fNITSRef = fITSReferences->GetEntriesFast();
  fNTRDRef = fTRDReferences->GetEntriesFast();
  fNTOFRef = fTOFReferences->GetEntriesFast();
  
  for (Int_t iref =0;iref<fTPCReferences->GetEntriesFast();iref++){
    AliTrackReference * ref = (AliTrackReference *) fTPCReferences->At(iref);
    //Float_t r = (ref->X()*ref->X()+ref->Y()*ref->Y());
    Float_t newdirection = ref->X()*ref->Px()+ref->Y()*ref->Py(); //inside or outside
    if (iref==0) direction = newdirection;
    if ( newdirection*direction<0){
      //changed direction
      direction = newdirection;
      fMCtracks+=1;
    }
    //rlast=r;			    
  }
  //
  // decay info
  fTPCdecay=kFALSE;
  if (fTRdecay.GetTrack()>0){
    fDecayCoord[0] = fTRdecay.X();
    fDecayCoord[1] = fTRdecay.Y();
    fDecayCoord[2] = fTRdecay.Z();
    if ( (fTRdecay.R()<250)&&(fTRdecay.R()>85) && (TMath::Abs(fTRdecay.Z())<250) ){
      fTPCdecay=kTRUE;     
    }
    else{
      fDecayCoord[0] = 0;
      fDecayCoord[1] = 0;
      fDecayCoord[2] = 0;
    }
  }
  //
  //
  //digits information update
  fRowsWithDigits    = fTPCRow.RowsOn();    
  fRowsWithDigitsInn = fTPCRow.RowsOn(63); // 63 = number of inner rows
  fRowsTrackLength   = fTPCRow.Last() - fTPCRow.First();
  //
  //
  // calculate primary ionization per cm  
  if (fParticle.GetPDG()){
    fMass = fParticle.GetMass();  
    fCharge = fParticle.GetPDG()->Charge();
    if (fTPCReferences->GetEntriesFast()>0){
      fTrackRef = *((AliTrackReference*)fTPCReferences->At(0));
    }
    if (fMass>0){
      Float_t p = TMath::Sqrt(fTrackRef.Px()*fTrackRef.Px()+
			      fTrackRef.Py()*fTrackRef.Py()+
			      fTrackRef.Pz()*fTrackRef.Pz());    
      if (p>0.001){
	Float_t betagama = p /fMass;
	fPrim = AliMathBase::BetheBlochAleph(betagama);
      }else fPrim=0;
    }
  }else{
    fMass =0;
    fPrim =0;
  }  
}

/////////////////////////////////////////////////////////////////////////////////
/*
void AliGenV0Info::Update()
{
  fMCPd[0] = fMCd.fParticle.Px();
  fMCPd[1] = fMCd.fParticle.Py();
  fMCPd[2] = fMCd.fParticle.Pz();
  fMCPd[3] = fMCd.fParticle.P();
  //
  fMCX[0]  = fMCd.fParticle.Vx();
  fMCX[1]  = fMCd.fParticle.Vy();
  fMCX[2]  = fMCd.fParticle.Vz();
  fMCR       = TMath::Sqrt( fMCX[0]*fMCX[0]+fMCX[1]*fMCX[1]);
  //
  fPdg[0]    = fMCd.fParticle.GetPdgCode();
  fPdg[1]    = fMCm.fParticle.GetPdgCode();
  //
  fLab[0]    = fMCd.fParticle.GetUniqueID();
  fLab[1]    = fMCm.fParticle.GetUniqueID();
  //  
}
*/

void AliGenV0Info::Update(Float_t vertex[3])
{
  fMCPd[0] = fMCd.fParticle.Px();
  fMCPd[1] = fMCd.fParticle.Py();
  fMCPd[2] = fMCd.fParticle.Pz();
  fMCPd[3] = fMCd.fParticle.P();
  //
  fMCX[0]  = fMCd.fParticle.Vx();
  fMCX[1]  = fMCd.fParticle.Vy();
  fMCX[2]  = fMCd.fParticle.Vz();
  fMCR       = TMath::Sqrt( fMCX[0]*fMCX[0]+fMCX[1]*fMCX[1]);
  //
  fPdg[0]    = fMCd.fParticle.GetPdgCode();
  fPdg[1]    = fMCm.fParticle.GetPdgCode();
  //
  fLab[0]    = fMCd.fParticle.GetUniqueID();
  fLab[1]    = fMCm.fParticle.GetUniqueID();
  //
  //
  //
  Double_t x1[3],p1[3];
  TParticle & pdaughter = fMCd.fParticle;
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
  TParticle & pmother = fMCm.fParticle;
  x2[0] = pmother.Vx();      
  x2[1] = pmother.Vy();      
  x2[2] = pmother.Vz();      
  p2[0] = pmother.Px();
  p2[1] = pmother.Py();
  p2[2] = pmother.Pz();
  //
  //
  sign = (pmother.GetPDG()->Charge()>0) ? -1:1;
  AliHelix mhelix(x2,p2,sign);
  
  //
  //
  //
  //find intersection linear
  //
  Double_t distance1, distance2;
  Double_t phase[2][2],radius[2];
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
  fPointAngleFi = (v[0]*p[0]+v[1]*p[1])/(vnorm2*pnorm2);
  fPointAngleTh = (v[2]*p[2]+vnorm2*pnorm2)/(vnorm3*pnorm3);  
  fPointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3);
  Double_t mass1 = fMCd.fMass;
  Double_t mass2 = fMCm.fMass;

  
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



/////////////////////////////////////////////////////////////////////////////////
void AliGenKinkInfo::Update()
{
  fMCPd[0] = fMCd.fParticle.Px();
  fMCPd[1] = fMCd.fParticle.Py();
  fMCPd[2] = fMCd.fParticle.Pz();
  fMCPd[3] = fMCd.fParticle.P();
  //
  fMCX[0]  = fMCd.fParticle.Vx();
  fMCX[1]  = fMCd.fParticle.Vy();
  fMCX[2]  = fMCd.fParticle.Vz();
  fMCR       = TMath::Sqrt( fMCX[0]*fMCX[0]+fMCX[1]*fMCX[1]);
  //
  fPdg[0]    = fMCd.fParticle.GetPdgCode();
  fPdg[1]    = fMCm.fParticle.GetPdgCode();
  //
  fLab[0]    = fMCd.fParticle.GetUniqueID();
  fLab[1]    = fMCm.fParticle.GetUniqueID();
  //
  //
  //
  Double_t x1[3],p1[3];
  TParticle & pdaughter = fMCd.fParticle;
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
  TParticle & pmother = fMCm.fParticle;
  x2[0] = pmother.Vx();      
  x2[1] = pmother.Vy();      
  x2[2] = pmother.Vz();      
  p2[0] = pmother.Px();
  p2[1] = pmother.Py();
  p2[2] = pmother.Pz();
  //
  AliTrackReference & pdecay = fMCm.fTRdecay;
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
  Double_t phase[2][2],radius[2];
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



  
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// End of implementation of the class AliMCInfo
//
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
digitRow::digitRow()
{
  Reset();
}
////////////////////////////////////////////////////////////////////////
digitRow & digitRow::operator=(const digitRow &digOld)
{
  for (Int_t i = 0; i<kgRowBytes; i++) fDig[i] = digOld.fDig[i];
  return (*this);
}
////////////////////////////////////////////////////////////////////////
void digitRow::SetRow(Int_t row) 
{
  if (row >= 8*kgRowBytes) {
    cerr<<"digitRow::SetRow: index "<<row<<" out of bounds."<<endl;
    return;
  }
  Int_t iC = row/8;
  Int_t iB = row%8;
  SETBIT(fDig[iC],iB);
}

////////////////////////////////////////////////////////////////////////
Bool_t digitRow::TestRow(Int_t row)
{
//
// return kTRUE if row is on
//
  Int_t iC = row/8;
  Int_t iB = row%8;
  return TESTBIT(fDig[iC],iB);
}
////////////////////////////////////////////////////////////////////////
Int_t digitRow::RowsOn(Int_t upto)
{
//
// returns number of rows with a digit  
// count only rows less equal row number upto
//
  Int_t total = 0;
  for (Int_t i = 0; i<kgRowBytes; i++) {
    for (Int_t j = 0; j < 8; j++) {
      if (i*8+j > upto) return total;
      if (TESTBIT(fDig[i],j))  total++;
    }
  }
  return total;
}
////////////////////////////////////////////////////////////////////////
void digitRow::Reset()
{
//
// resets all rows to zero
//
  for (Int_t i = 0; i<kgRowBytes; i++) {
    fDig[i] <<= 8;
  }
}
////////////////////////////////////////////////////////////////////////
Int_t digitRow::Last()
{
//
// returns the last row number with a digit
// returns -1 if now digits 
//
  for (Int_t i = kgRowBytes-1; i>=0; i--) {
    for (Int_t j = 7; j >= 0; j--) {
      if TESTBIT(fDig[i],j) return i*8+j;
    }
  }
  return -1;
}
////////////////////////////////////////////////////////////////////////
Int_t digitRow::First()
{
//
// returns the first row number with a digit
// returns -1 if now digits 
//
  for (Int_t i = 0; i<kgRowBytes; i++) {
    for (Int_t j = 0; j < 8; j++) {
      if (TESTBIT(fDig[i],j)) return i*8+j;
    }
  }
  return -1;
}

////////////////////////////////////////////////////////////////////////
//
// end of implementation of a class digitRow
//
////////////////////////////////////////////////////////////////////////
  
////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::AliGenInfoMaker()
{
  Reset();
}

////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::AliGenInfoMaker(const char * fnGalice, const char* fnRes,
				   Int_t nEvents, Int_t firstEvent)
{
  Reset();
  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  //   fFnRes = fnRes;
  sprintf(fFnRes,"%s",fnRes);
  //
  fLoader = AliRunLoader::Open(fnGalice);
  if (gAlice){
    delete AliRunLoader::GetRunLoader();
    delete gAlice;
    gAlice = 0x0;
  }
  if (fLoader->LoadgAlice()){
    cerr<<"Error occured while l"<<endl;
  }
  Int_t nall = fLoader->GetNumberOfEvents();
  if (nEvents==0) {
    nEvents =nall;
    fNEvents=nall;
    fFirstEventNr=0;
  }    

  if (nall<=0){
    cerr<<"no events available"<<endl;
    fEventNr = 0;
    return;
  }
  if (firstEvent+nEvents>nall) {
    fEventNr = nall-firstEvent;
    cerr<<"restricted number of events availaible"<<endl;
  }
  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf,0);
}


AliMCInfo * AliGenInfoMaker::MakeInfo(UInt_t i)
{
  // 
  if (i<fNParticles) {
    if (fGenInfo[i]) return  fGenInfo[i];
    fGenInfo[i] = new AliMCInfo;  
    fNInfos++;
    return fGenInfo[i];
  }
  else 
    return 0;  
}

////////////////////////////////////////////////////////////////////////
void AliGenInfoMaker::Reset()
{
  fEventNr = 0;
  fNEvents = 0;
  fTreeGenTracks = 0;
  fFileGenTracks = 0;
  fGenInfo = 0;
  fNInfos  = 0;
  //
  //
  fDebug = 0;
  fVPrim[0] = -1000.;
  fVPrim[1] = -1000.;
  fVPrim[2] = -1000.;
  fParamTPC = 0;
}
////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::~AliGenInfoMaker()
{
  
  if (fLoader){
    fLoader->UnloadgAlice();
    gAlice = 0;
    delete fLoader;
  }
}

Int_t  AliGenInfoMaker::SetIO()
{
  //
  // 
  CreateTreeGenTracks();
  if (!fTreeGenTracks) return 1;
  //  AliTracker::SetFieldFactor(); 
 
  fParamTPC = GetTPCParam();
  //
  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::SetIO(Int_t eventNr)
{
  //
  // 
  // SET INPUT
  fLoader->SetEventNumber(eventNr);
  //
  fLoader->LoadHeader();
  fLoader->LoadKinematics();  
  fStack = fLoader->Stack();
  //
  fLoader->LoadTrackRefs();
  fLoader->LoadHits();
  fTreeTR = fLoader->TreeTR();
  //
  AliTPCLoader * tpcl = (AliTPCLoader*)fLoader->GetLoader("TPCLoader");
  tpcl->LoadDigits();
  fTreeD = tpcl->TreeD();  
  return 0;
}

Int_t AliGenInfoMaker::CloseIOEvent()
{
  fLoader->UnloadHeader();
  fLoader->UnloadKinematics();
  fLoader->UnloadTrackRefs();
  AliTPCLoader * tpcl = (AliTPCLoader*)fLoader->GetLoader("TPCLoader");
  tpcl->UnloadDigits();
  return 0;
}

Int_t AliGenInfoMaker::CloseIO()
{
  fLoader->UnloadgAlice();
  return 0;
}



////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::Exec(Int_t nEvents, Int_t firstEventNr)
{
  fNEvents = nEvents;
  fFirstEventNr = firstEventNr;
  return Exec();
}

////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::Exec()  
{
  TStopwatch timer;
  timer.Start();
  Int_t status =SetIO();
  if (status>0) return status;
  //

  for (fEventNr = fFirstEventNr; fEventNr < fFirstEventNr+fNEvents;
       fEventNr++) {
    SetIO(fEventNr);
    fNParticles = fStack->GetNtrack();
    //
    fGenInfo = new AliMCInfo*[fNParticles];
    for (UInt_t i = 0; i<fNParticles; i++) {
      fGenInfo[i]=0; 
    }
    //
    cout<<"Start to process event "<<fEventNr<<endl;
    cout<<"\tfNParticles = "<<fNParticles<<endl;
    if (fDebug>2) cout<<"\n\n\n\tStart loop over TreeTR"<<endl;
    if (TreeTRLoop()>0) return 1;
    //
    if (fDebug>2) cout<<"\n\n\n\tStart loop over TreeD"<<endl;
    if (TreeDLoop()>0) return 1;
    //
    if (fDebug>2) cout<<"\n\n\n\tStart loop over TreeK"<<endl;
    if (TreeKLoop()>0) return 1;
    if (fDebug>2) cout<<"\tEnd loop over TreeK"<<endl;
    //
    if (BuildKinkInfo()>0) return 1;
    if (BuildV0Info()>0) return 1;
    //if (BuildHitLines()>0) return 1;
    if (fDebug>2) cout<<"\tEnd loop over TreeK"<<endl;
    //
    for (UInt_t i = 0; i<fNParticles; i++) {
      if (fGenInfo[i]) delete fGenInfo[i]; 
    }
    delete []fGenInfo;
    CloseIOEvent();
  }
  //
  CloseIO();
  CloseOutputFile();

  cerr<<"Exec finished"<<endl;

  timer.Stop();
  timer.Print();
  return 0;
}
////////////////////////////////////////////////////////////////////////
void AliGenInfoMaker::CreateTreeGenTracks() 
{
  fFileGenTracks = TFile::Open(fFnRes,"RECREATE");
  if (!fFileGenTracks) {
    cerr<<"Error in CreateTreeGenTracks: cannot open file "<<fFnRes<<endl;
    return;
  }
  fTreeGenTracks = new TTree("genTracksTree","genTracksTree");  
  AliMCInfo * info = new AliMCInfo;
  fTreeGenTracks->Branch("MC","AliMCInfo",&info);
  delete info; 
  //
  AliGenKinkInfo *kinkinfo = new AliGenKinkInfo;
  fTreeKinks = new TTree("genKinksTree","genKinksTree");  
  fTreeKinks->Branch("MC","AliGenKinkInfo",&kinkinfo);
  delete kinkinfo;
  //
  AliGenV0Info *v0info = new AliGenV0Info;
  fTreeV0 = new TTree("genV0Tree","genV0Tree");  
  fTreeV0->Branch("MC","AliGenV0Info",&v0info);
  delete v0info;
  //
  //
  AliPointsMI * points0 = new AliPointsMI;  
  AliPointsMI * points1 = new AliPointsMI;  
  AliPointsMI * points2 = new AliPointsMI;  
  fTreeHitLines = new TTree("HitLines","HitLines");
  fTreeHitLines->Branch("TPC.","AliPointsMI",&points0,32000,10);
  fTreeHitLines->Branch("TRD.","AliPointsMI",&points1,32000,10);
  fTreeHitLines->Branch("ITS.","AliPointsMI",&points2,32000,10);
  Int_t label=0;
  fTreeHitLines->Branch("Label",&label,"label/I");
  //
  fTreeGenTracks->AutoSave();
  fTreeKinks->AutoSave();
  fTreeV0->AutoSave();
  fTreeHitLines->AutoSave();
}
////////////////////////////////////////////////////////////////////////
void AliGenInfoMaker::CloseOutputFile() 
{
  if (!fFileGenTracks) {
    cerr<<"File "<<fFnRes<<" not found as an open file."<<endl;
    return;
  }
  fFileGenTracks->cd();
  fTreeGenTracks->Write();  
  delete fTreeGenTracks;
  fTreeKinks->Write();
  delete fTreeKinks;
  fTreeV0->Write();
  delete fTreeV0;

  fFileGenTracks->Close();
  delete fFileGenTracks;
  return;
}

////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::TreeKLoop()
{
//
// open the file with treeK
// loop over all entries there and save information about some tracks
//

  AliStack * stack = fStack;
  if (!stack) {cerr<<"Stack was not found!\n"; return 1;}
  
  if (fDebug > 0) {
    cout<<"There are "<<fNParticles<<" primary and secondary particles in event "
	<<fEventNr<<endl;
  }  
  Int_t  ipdg = 0;
  TParticlePDG *ppdg = 0;
  // not all generators give primary vertex position. Take the vertex
  // of the particle 0 as primary vertex.
  TDatabasePDG  pdg; //get pdg table  
  //thank you very much root for this
  TBranch * br = fTreeGenTracks->GetBranch("MC");
  TParticle *particle = stack->ParticleFromTreeK(0);
  fVPrim[0] = particle->Vx();
  fVPrim[1] = particle->Vy();
  fVPrim[2] = particle->Vz();
  for (UInt_t iParticle = 0; iParticle < fNParticles; iParticle++) {
    // load only particles with TR
    AliMCInfo * info = GetInfo(iParticle);
    if (!info) continue;
    //////////////////////////////////////////////////////////////////////
    info->fLabel = iParticle;
    //
    info->fParticle = *(stack->Particle(iParticle));
    info->fVDist[0] = info->fParticle.Vx()-fVPrim[0];
    info->fVDist[1] = info->fParticle.Vy()-fVPrim[1];
    info->fVDist[2] = info->fParticle.Vz()-fVPrim[2]; 
    info->fVDist[3] = TMath::Sqrt(info->fVDist[0]*info->fVDist[0]+
				  info->fVDist[1]*info->fVDist[1]+info->fVDist[2]*info->fVDist[2]);
    //
    //
    ipdg = info->fParticle.GetPdgCode();
    info->fPdg = ipdg;
    ppdg = pdg.GetParticle(ipdg);   	   
    info->fEventNr = fEventNr;
    info->Update();
    //////////////////////////////////////////////////////////////////////    
    br->SetAddress(&info);    
    fTreeGenTracks->Fill();    
  }
  fTreeGenTracks->AutoSave();
  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////
Int_t  AliGenInfoMaker::BuildKinkInfo()
{
  //
  // Build tree with Kink Information
  //
  AliStack * stack = fStack;
  if (!stack) {cerr<<"Stack was not found!\n"; return 1;}
  
  if (fDebug > 0) {
    cout<<"There are "<<fNParticles<<" primary and secondary particles in event "
	<<fEventNr<<endl;
  }  
  //  Int_t  ipdg = 0;
  //TParticlePDG *ppdg = 0;
  // not all generators give primary vertex position. Take the vertex
  // of the particle 0 as primary vertex.
  TDatabasePDG  pdg; //get pdg table  
  //thank you very much root for this
  TBranch * br = fTreeKinks->GetBranch("MC");
  //
  AliGenKinkInfo * kinkinfo = new AliGenKinkInfo;
  //
  for (UInt_t iParticle = 0; iParticle < fNParticles; iParticle++) {
    // load only particles with TR
    AliMCInfo * info = GetInfo(iParticle);
    if (!info) continue;
    if (info->fCharge==0) continue;  
    if (info->fTRdecay.P()<0.13) continue;  //momenta cut 
    if (info->fTRdecay.R()>500)  continue;  //R cut - decay outside barrel
    TParticle & particle = info->fParticle;
    Int_t first = particle.GetDaughter(0);
    if (first ==0) continue;

    Int_t last = particle.GetDaughter(1);
    if (last ==0) last=first;
    AliMCInfo * info2 =  0;
    AliMCInfo * dinfo =  0;
    
    for (Int_t id2=first;id2<=last;id2++){
      info2 = GetInfo(id2);
      if (!info2) continue;
      TParticle &p2 = info2->fParticle;
      Double_t r = TMath::Sqrt(p2.Vx()*p2.Vx()+p2.Vy()*p2.Vy());
      if ( TMath::Abs(info->fTRdecay.R()-r)>0.01) continue;
      if (!(p2.GetPDG())) continue;
      dinfo =info2;
     
      kinkinfo->fMCm = (*info);
      kinkinfo->fMCd = (*dinfo);
      if (kinkinfo->fMCm.fParticle.GetPDG()==0 ||  kinkinfo->fMCd.fParticle.GetPDG()==0) continue;
      kinkinfo->Update();
      br->SetAddress(&kinkinfo);    
      fTreeKinks->Fill();
    }
    /*
    if (dinfo){
      kinkinfo->fMCm = (*info);
      kinkinfo->fMCd = (*dinfo);
      kinkinfo->Update();
      br->SetAddress(&kinkinfo);    
      fTreeKinks->Fill();
    }
    */
  }
  fTreeGenTracks->AutoSave();
  if (fDebug > 2) cerr<<"end of Kink Loop"<<endl;
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////
Int_t  AliGenInfoMaker::BuildV0Info()
{
  //
  // Build tree with V0 Information
  //
  AliStack * stack = fStack;
  if (!stack) {cerr<<"Stack was not found!\n"; return 1;}
  
  if (fDebug > 0) {
    cout<<"There are "<<fNParticles<<" primary and secondary particles in event "
	<<fEventNr<<endl;
  }  
  //  Int_t  ipdg = 0;
  //TParticlePDG *ppdg = 0;
  // not all generators give primary vertex position. Take the vertex
  // of the particle 0 as primary vertex.
  TDatabasePDG  pdg; //get pdg table  
  //thank you very much root for this
  TBranch * br = fTreeV0->GetBranch("MC");
  //
  AliGenV0Info * v0info = new AliGenV0Info;
  //
  //
  for (UInt_t iParticle = 0; iParticle < fNParticles; iParticle++) {
    // load only particles with TR
    AliMCInfo * info = GetInfo(iParticle);
    if (!info) continue;
    if (info->fCharge==0) continue;  
    //
    //
    TParticle & particle = info->fParticle;
    Int_t mother = particle.GetMother(0);
    if (mother <=0) continue;
    //
    TParticle * motherparticle = stack->Particle(mother);
    if (!motherparticle) continue;
    //
    Int_t last = motherparticle->GetDaughter(1);
    if (last==(int)iParticle) continue;
    AliMCInfo * info2 =  info;
    AliMCInfo * dinfo =  GetInfo(last);
    if (!dinfo) continue;
    if (!dinfo->fParticle.GetPDG()) continue;
    if (!info2->fParticle.GetPDG()) continue;
    if (dinfo){
      v0info->fMCm = (*info);
      v0info->fMCd = (*dinfo);
      v0info->fMotherP = (*motherparticle);
      v0info->Update(fVPrim);
      br->SetAddress(&v0info);    
      fTreeV0->Fill();
    }
  }
  fTreeV0->AutoSave();
  if (fDebug > 2) cerr<<"end of V0 Loop"<<endl;
  return 0;
}




////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::BuildHitLines()
{

//
// open the file with treeK
// loop over all entries there and save information about some tracks
//

  AliStack * stack = fStack;
  if (!stack) {cerr<<"Stack was not found!\n"; return 1;}
  
  if (fDebug > 0) {
    cout<<"There are "<<fNParticles<<" primary and secondary particles in event "
	<<fEventNr<<endl;
  }  
  Int_t  ipdg = 0;
  // TParticlePDG *ppdg = 0;
  // not all generators give primary vertex position. Take the vertex
  // of the particle 0 as primary vertex.
  TDatabasePDG  pdg; //get pdg table  
  //thank you very much root for this
  AliPointsMI *tpcp = new AliPointsMI;
  AliPointsMI *trdp = new AliPointsMI;
  AliPointsMI *itsp = new AliPointsMI;
  Int_t label =0;
  //
  TBranch * brtpc = fTreeHitLines->GetBranch("TPC.");
  TBranch * brtrd = fTreeHitLines->GetBranch("TRD.");  
  TBranch * brits = fTreeHitLines->GetBranch("ITS.");
  TBranch * brlabel = fTreeHitLines->GetBranch("Label");
  brlabel->SetAddress(&label);
  brtpc->SetAddress(&tpcp);
  brtrd->SetAddress(&trdp);
  brits->SetAddress(&itsp);
  //
  AliDetector *dtpc = gAlice->GetDetector("TPC");
  AliDetector *dtrd = gAlice->GetDetector("TRD");
  AliDetector *dits = gAlice->GetDetector("ITS");
 
  for (UInt_t iParticle = 0; iParticle < fNParticles; iParticle++) {
    // load only particles with TR
    AliMCInfo * info = GetInfo(iParticle);
    if (!info) continue;
    Int_t primpart = info->fPrimPart;
    ipdg = info->fParticle.GetPdgCode();
    label = iParticle;
    //
    gAlice->ResetHits();
    tpcp->Reset();
    itsp->Reset();
    trdp->Reset();
    tpcp->fLabel1 = ipdg;
    trdp->fLabel1 = ipdg;
    itsp->fLabel1 = ipdg;
    if (dtpc->TreeH()->GetEvent(primpart)){
      dtpc->LoadPoints(primpart);
      tpcp->Reset(dtpc,iParticle);
    }
    if (dtrd->TreeH()->GetEvent(primpart)){
      dtrd->LoadPoints(primpart);
      trdp->Reset(dtrd,iParticle);
    }
    if (dits->TreeH()->GetEvent(primpart)){
      dits->LoadPoints(primpart);
      itsp->Reset(dits,iParticle);
    }    
    //    
    fTreeHitLines->Fill();
    dtpc->ResetPoints();
    dtrd->ResetPoints();
    dits->ResetPoints();
  }
  delete tpcp;
  delete trdp;
  delete itsp;
  fTreeHitLines->AutoSave();
  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;
  return 0;
}


////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::TreeDLoop()
{
  //
  // open the file with treeD
  // loop over all entries there and save information about some tracks
  //
  
  Int_t nInnerSector = fParamTPC->GetNInnerSector();
  Int_t rowShift = 0;
  Int_t zero=fParamTPC->GetZeroSup()+6;  
  //
  //
  AliSimDigits digitsAddress, *digits=&digitsAddress;
  fTreeD->GetBranch("Segment")->SetAddress(&digits);
  
  Int_t sectorsByRows=(Int_t)fTreeD->GetEntries();
  if (fDebug > 1) cout<<"\tsectorsByRows = "<<sectorsByRows<<endl;
  for (Int_t i=0; i<sectorsByRows; i++) {
    if (!fTreeD->GetEvent(i)) continue;
    Int_t sec,row;
    fParamTPC->AdjustSectorRow(digits->GetID(),sec,row);
    if (fDebug > 1) cout<<sec<<' '<<row<<"                          \r";
    // here I expect that upper sectors follow lower sectors
    if (sec > nInnerSector) rowShift = fParamTPC->GetNRowLow();
    //
    digits->ExpandTrackBuffer();
    digits->First();        
    do {
      Int_t iRow=digits->CurrentRow();
      Int_t iColumn=digits->CurrentColumn();
      Short_t digitValue = digits->CurrentDigit();
      if (digitValue >= zero) {
	Int_t label;
	for (Int_t j = 0; j<3; j++) {
	  //	  label = digits->GetTrackID(iRow,iColumn,j);
	  label = digits->GetTrackIDFast(iRow,iColumn,j)-2; 	  
	  if (label >= (Int_t)fNParticles) { //don't label from bakground event
	    continue;
	  }
	  if (label >= 0 && label <= (Int_t)fNParticles) {
	    if (fDebug > 6 ) {
	      cout<<"Inner loop: sector, iRow, iColumn, label, value, row "
		  <<sec<<" "
		  <<iRow<<" "<<iColumn<<" "<<label<<" "<<digitValue
		  <<" "<<row<<endl;
	    }	
	    AliMCInfo * info = GetInfo(label);
	    if (info){
	      info->fTPCRow.SetRow(row+rowShift);
	    }
	  }
	}
      }
    } while (digits->Next());
  }
  
  if (fDebug > 2) cerr<<"end of TreeDLoop"<<endl;  
  return 0;
}


////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::TreeTRLoop()
{
  //
  // loop over TrackReferences and store the first one for each track
  //  
  TTree * treeTR = fTreeTR;
  Int_t nPrimaries = (Int_t) treeTR->GetEntries();
  if (fDebug > 1) cout<<"There are "<<nPrimaries<<" entries in TreeTR"<<endl;
  //
  //
  //track references for TPC
  TClonesArray* TPCArrayTR = new TClonesArray("AliTrackReference");
  TClonesArray* ITSArrayTR = new TClonesArray("AliTrackReference");
  TClonesArray* TRDArrayTR = new TClonesArray("AliTrackReference");
  TClonesArray* TOFArrayTR = new TClonesArray("AliTrackReference");
  TClonesArray* RunArrayTR = new TClonesArray("AliTrackReference");
  //
  if (treeTR->GetBranch("TPC"))    treeTR->GetBranch("TPC")->SetAddress(&TPCArrayTR);
  if (treeTR->GetBranch("ITS"))    treeTR->GetBranch("ITS")->SetAddress(&ITSArrayTR);
  if (treeTR->GetBranch("TRD"))    treeTR->GetBranch("TRD")->SetAddress(&TRDArrayTR);
  if (treeTR->GetBranch("TOF"))    treeTR->GetBranch("TOF")->SetAddress(&TOFArrayTR);
  if (treeTR->GetBranch("AliRun")) treeTR->GetBranch("AliRun")->SetAddress(&RunArrayTR);
  //
  //
  //
  for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) {
    treeTR->GetEntry(iPrimPart);
    //
    // Loop over TPC references
    //
    for (Int_t iTrackRef = 0; iTrackRef < TPCArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)TPCArrayTR->At(iTrackRef);            
      //
      if (trackRef->TestBit(BIT(2))){  
	//if decay 
	if (trackRef->P()<fgTPCPtCut) continue;
	Int_t label = trackRef->GetTrack(); 
	AliMCInfo * info = GetInfo(label);
	if (!info) info = MakeInfo(label);
	info->fTRdecay = *trackRef;      
      }
      //
      if (trackRef->P()<fgTPCPtCut) continue;
      Int_t label = trackRef->GetTrack();      
      AliMCInfo * info = GetInfo(label);
      if (!info) info = MakeInfo(label);
      if (!info) continue;
      info->fPrimPart =  iPrimPart;
      TClonesArray & arr = *(info->fTPCReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    //
    // Loop over ITS references
    //
    for (Int_t iTrackRef = 0; iTrackRef < ITSArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)ITSArrayTR->At(iTrackRef);            
      // 
      //
      if (trackRef->P()<fgTPCPtCut) continue;
      Int_t label = trackRef->GetTrack();      
      AliMCInfo * info = GetInfo(label);
      if ( (!info) && trackRef->Pt()<fgITSPtCut) continue;
      if (!info) info = MakeInfo(label);
      if (!info) continue;
      info->fPrimPart =  iPrimPart;
      TClonesArray & arr = *(info->fITSReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    //
    // Loop over TRD references
    //
    for (Int_t iTrackRef = 0; iTrackRef < TRDArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)TRDArrayTR->At(iTrackRef);            
      //
      if (trackRef->P()<fgTPCPtCut) continue;
      Int_t label = trackRef->GetTrack();      
      AliMCInfo * info = GetInfo(label);
      if ( (!info) && trackRef->Pt()<fgTRDPtCut) continue;
      if (!info) info = MakeInfo(label);
      if (!info) continue;
      info->fPrimPart =  iPrimPart;
      TClonesArray & arr = *(info->fTRDReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    //
    // Loop over TOF references
    //
    for (Int_t iTrackRef = 0; iTrackRef < TOFArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)TOFArrayTR->At(iTrackRef);            
      Int_t label = trackRef->GetTrack();      
      AliMCInfo * info = GetInfo(label);
      if (!info){
	if (trackRef->Pt()<fgTPCPtCut) continue;
	if ( (!info) && trackRef->Pt()<fgTOFPtCut) continue;
      }
      if (!info) info = MakeInfo(label);
      if (!info) continue;
      info->fPrimPart =  iPrimPart;
      TClonesArray & arr = *(info->fTOFReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    //
    // get dacay position
    for (Int_t iTrackRef = 0; iTrackRef < RunArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)RunArrayTR->At(iTrackRef);      
      //
      Int_t label = trackRef->GetTrack();
      AliMCInfo * info = GetInfo(label);
      if (!info) continue;
      if (!trackRef->TestBit(BIT(2))) continue;  //if not decay
      //      if (TMath::Abs(trackRef.X());
      info->fTRdecay = *trackRef;      
    }
  }
  //
  TPCArrayTR->Delete();
  delete  TPCArrayTR;
  TRDArrayTR->Delete();
  delete  TRDArrayTR;
  TOFArrayTR->Delete();
  delete  TOFArrayTR;

  ITSArrayTR->Delete();
  delete  ITSArrayTR;
  RunArrayTR->Delete();
  delete  RunArrayTR;
  //
  return 0;
}

////////////////////////////////////////////////////////////////////////
Float_t AliGenInfoMaker::TR2LocalX(AliTrackReference *trackRef,
				    AliTPCParam *paramTPC) {

  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2(x,index);
  return x[0];
}
////////////////////////////////////////////////////////////////////////



TH1F * AliComparisonDraw::DrawXY(const char * chx, const char *chy, const char* selection, 
		const char * quality, Int_t nbins, Float_t minx, Float_t maxx, Float_t miny, Float_t maxy, Int_t nBinsRes)
{
  //
  Double_t* bins = CreateLogBins(nbins, minx, maxx);
  TH2F* hRes2 = new TH2F("hRes2", "residuals", nbins, minx, maxx, nBinsRes, miny, maxy);
  char cut[1000];
  sprintf(cut,"%s&&%s",selection,quality);
  char expression[1000];
  sprintf(expression,"%s:%s>>hRes2",chy,chx);
  fTree->Draw(expression, cut, "groff");
  TH1F* hMean=0;
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, chx, chy);
  //
  delete hRes2;
  delete[] bins;
  fRes  = hRes;
  fMean = hMean;
  return hRes;
}

TH1F * AliComparisonDraw::DrawLogXY(const char * chx, const char *chy, const char* selection, 
				    const char * quality, Int_t nbins, Float_t minx, Float_t maxx, Float_t miny, Float_t maxy, Int_t nBinsRes)
{
  //
  Double_t* bins = CreateLogBins(nbins, minx, maxx);
  TH2F* hRes2 = new TH2F("hRes2", "residuals", nbins, bins, nBinsRes, miny, maxy);
  char cut[1000];
  sprintf(cut,"%s&&%s",selection,quality);
  char expression[1000];
  sprintf(expression,"%s:%s>>hRes2",chy,chx);
  fTree->Draw(expression, cut, "groff");
  TH1F* hMean=0;  
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, chx, chy);
  //
  delete hRes2;
  delete[] bins;
  fRes  = hRes;
  fMean = hMean;
  return hRes;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
TH1F * AliComparisonDraw::Eff(const char *variable, const char* selection, const char * quality, 
			      Int_t nbins, Float_t min, Float_t max)
{
  //
  //
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nbins, min, max);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nbins, min, max);
  char inputGen[1000];  
  sprintf(inputGen,"%s>>hGen", variable);
  fTree->Draw(inputGen, selection, "groff");
  char selectionRec[256];
  sprintf(selectionRec, "%s && %s", selection, quality);
  char inputRec[1000];  
  sprintf(inputRec,"%s>>hRec", variable);
  fTree->Draw(inputRec, selectionRec, "groff");
  //
  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, variable, "#epsilon [%]");
  fRes = hEff;
  delete hRec;
  delete hGen;
  return hEff;
}



///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
TH1F * AliComparisonDraw::EffLog(const char *variable, const char* selection, const char * quality, 
			      Int_t nbins, Float_t min, Float_t max)
{
  //
  //
  Double_t* bins = CreateLogBins(nbins, min, max);
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nbins, bins);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nbins, bins);
  char inputGen[1000];  
  sprintf(inputGen,"%s>>hGen", variable);
  fTree->Draw(inputGen, selection, "groff");
  char selectionRec[256];
  sprintf(selectionRec, "%s && %s", selection, quality);
  char inputRec[1000];  
  sprintf(inputRec,"%s>>hRec", variable);
  fTree->Draw(inputRec, selectionRec, "groff");
  //
  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, variable, "#epsilon [%]");
  fRes = hEff;
  delete hRec;
  delete hGen;
  delete[] bins;
  return hEff;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

Double_t* AliComparisonDraw::CreateLogBins(Int_t nBins, Double_t xMin, Double_t xMax)
{
  Double_t* bins = new Double_t[nBins+1];
  bins[0] = xMin;
  Double_t factor = pow(xMax/xMin, 1./nBins);
  for (Int_t i = 1; i <= nBins; i++)
    bins[i] = factor * bins[i-1];
  return bins;
}




void AliComparisonDraw::AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle)
{
  //
  histo->GetXaxis()->SetTitle(xAxisTitle);
  histo->GetYaxis()->SetTitle(yAxisTitle);
}


TH1F* AliComparisonDraw::CreateEffHisto(TH1F* hGen, TH1F* hRec)
{
  //
  Int_t nBins = hGen->GetNbinsX();
  TH1F* hEff = (TH1F*) hGen->Clone("hEff");
  hEff->SetTitle("");
  hEff->SetStats(kFALSE);
  hEff->SetMinimum(0.);
  hEff->SetMaximum(110.);
  //
  for (Int_t iBin = 0; iBin <= nBins; iBin++) {
    Double_t nGen = hGen->GetBinContent(iBin);
    Double_t nRec = hRec->GetBinContent(iBin);
    if (nGen > 0) {
      Double_t eff = nRec/nGen;
      hEff->SetBinContent(iBin, 100. * eff);
      Double_t error = sqrt((eff*(1.-eff)+0.01) / nGen);      
      //      if (error == 0) error = sqrt(0.1/nGen);
      //
      if (error == 0) error = 0.0001;
      hEff->SetBinError(iBin, 100. * error);
    } else {
      hEff->SetBinContent(iBin, 100. * 0.5);
      hEff->SetBinError(iBin, 100. * 0.5);
    }
  }
  return hEff;
}



TH1F* AliComparisonDraw::CreateResHisto(TH2F* hRes2, TH1F **phMean,  Bool_t drawBinFits, 
		     Bool_t overflowBinFits)
{
  TVirtualPad* currentPad = gPad;
  TAxis* axis = hRes2->GetXaxis();
  Int_t nBins = axis->GetNbins();
  TH1F* hRes, *hMean;
  if (axis->GetXbins()->GetSize()){
    hRes = new TH1F("hRes", "", nBins, axis->GetXbins()->GetArray());
    hMean = new TH1F("hMean", "", nBins, axis->GetXbins()->GetArray());
  }
  else{
    hRes = new TH1F("hRes", "", nBins, axis->GetXmin(), axis->GetXmax());
    hMean = new TH1F("hMean", "", nBins, axis->GetXmin(), axis->GetXmax());

  }
  hRes->SetStats(false);
  hRes->SetOption("E");
  hRes->SetMinimum(0.);
  //
  hMean->SetStats(false);
  hMean->SetOption("E");
 
  // create the fit function
  //TKFitGaus* fitFunc = new TKFitGaus("resFunc");
  //   TF1 * fitFunc = new TF1("G","[3]+[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  TF1 * fitFunc = new TF1("G","[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  // create canvas for fits
  TCanvas* canBinFits = NULL;
  Int_t nPads = (overflowBinFits) ? nBins+2 : nBins;
  Int_t nx = Int_t(sqrt(nPads-1.));// + 1;
  Int_t ny = (nPads-1) / nx + 1;
  if (drawBinFits) {
    canBinFits = (TCanvas*)gROOT->FindObject("canBinFits");
    if (canBinFits) delete canBinFits;
    canBinFits = new TCanvas("canBinFits", "fits of bins", 200, 100, 500, 700);
    canBinFits->Divide(nx, ny);
  }

  // loop over x bins and fit projection
  Int_t dBin = ((overflowBinFits) ? 1 : 0);
  for (Int_t bin = 1-dBin; bin <= nBins+dBin; bin++) {
    if (drawBinFits) canBinFits->cd(bin + dBin);
    TH1D* hBin = hRes2->ProjectionY("hBin", bin, bin);
    //    
    if (hBin->GetEntries() > 5) {
      fitFunc->SetParameters(hBin->GetMaximum(),hBin->GetMean(),hBin->GetRMS());
      hBin->Fit(fitFunc,"s");
      Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));

      if (sigma > 0.){
	hRes->SetBinContent(bin, TMath::Abs(fitFunc->GetParameter(2)));
	hMean->SetBinContent(bin, fitFunc->GetParameter(1));	
      }
      else{
	hRes->SetBinContent(bin, 0.);
	hMean->SetBinContent(bin,0);
      }
      hRes->SetBinError(bin, fitFunc->GetParError(2));
      hMean->SetBinError(bin, fitFunc->GetParError(1));
      
      //
      //

    } else {
      hRes->SetBinContent(bin, 0.);
      hRes->SetBinError(bin, 0.);
      hMean->SetBinContent(bin, 0.);
      hMean->SetBinError(bin, 0.);
    }
    

    if (drawBinFits) {
      char name[256];
      if (bin == 0) {
	sprintf(name, "%s < %.4g", axis->GetTitle(), axis->GetBinUpEdge(bin));
      } else if (bin == nBins+1) {
	sprintf(name, "%.4g < %s", axis->GetBinLowEdge(bin), axis->GetTitle());
      } else {
	sprintf(name, "%.4g < %s < %.4g", axis->GetBinLowEdge(bin),
		axis->GetTitle(), axis->GetBinUpEdge(bin));
      }
      canBinFits->cd(bin + dBin);
      hBin->SetTitle(name);
      hBin->SetStats(kTRUE);
      hBin->DrawCopy("E");
      canBinFits->Update();
      canBinFits->Modified();
      canBinFits->Update();
    }
    
    delete hBin;
  }

  delete fitFunc;
  currentPad->cd();
  *phMean = hMean;
  return hRes;
}


void   AliComparisonDraw::DrawFriend2D(const char * chx, const char *chy, const char* selection, TTree * tfriend)
{
  
}


void   AliComparisonDraw::GetPoints3D(const char * label, const char * chpoints, const char* selection, TTree * tpoints, Int_t color,Float_t rmin){
  //
  //
  if (!fPoints) fPoints = new TObjArray;
  if (tpoints->GetIndex()==0) tpoints->BuildIndex("Label","Label");
  AliPointsMI * points = new AliPointsMI;
  TBranch * br = tpoints->GetBranch(chpoints);
  br->SetAddress(&points);
  Int_t npoints = fTree->Draw(label,selection);
  Float_t xyz[30000];
  rmin*=rmin;
  for (Int_t i=0;i<npoints;i++){
    Int_t index = (Int_t)fTree->GetV1()[i];
    tpoints->GetEntryWithIndex(index,index);
    if (points->fN<2) continue;
    Int_t goodpoints=0;
    for (Int_t i=0;i<points->fN; i++){
      xyz[goodpoints*3]   = points->fX[i];
      xyz[goodpoints*3+1] = points->fY[i];
      xyz[goodpoints*3+2] = points->fZ[i];
      if ( points->fX[i]*points->fX[i]+points->fY[i]*points->fY[i]>rmin) goodpoints++;
    } 
    TPolyMarker3D * marker = new TPolyMarker3D(goodpoints,xyz); 
    //    marker->SeWidth(1); 
    marker->SetMarkerColor(color);
    marker->SetMarkerStyle(1);
    fPoints->AddLast(marker);
  }
}

void AliComparisonDraw::Draw3D(Int_t min, Int_t max){
  if (!fPoints) return;
  Int_t npoints = fPoints->GetEntriesFast();
  max = TMath::Min(npoints,max);
  min = TMath::Min(npoints,min);

  for (Int_t i =min; i<max;i++){
    TObject * obj = fPoints->At(i);
    if (obj) obj->Draw();
  }
}

void AliComparisonDraw:: SavePoints(const char* name)
{
  char fname[25];
  sprintf(fname,"TH_%s.root",name);
  TFile f(fname,"new");
  fPoints->Write("Tracks",TObject::kSingleKey);
  f.Close();
  sprintf(fname,"TH%s.gif",name);
  fCanvas->SaveAs(fname);
}

void   AliComparisonDraw::InitView(){
  //
  //
  AliRunLoader* rl = AliRunLoader::Open();
  rl->GetEvent(0); 
  rl->CdGAFile(); 
  //
  fCanvas = new TCanvas("cdisplay", "Cluster display",0,0,700,730);
  fView   = new TView(1);
  fView->SetRange(-330,-360,-330, 360,360, 330);
  fCanvas->Clear();
  fCanvas->SetFillColor(1);
  fCanvas->SetTheta(90.);
  fCanvas->SetPhi(0.);
  //
  TGeometry *geom =(TGeometry*)gDirectory->Get("AliceGeom");
  geom->Draw("same");

}
