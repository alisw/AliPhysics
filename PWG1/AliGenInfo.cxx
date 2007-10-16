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
#include "TSystem.h"
#include "TCanvas.h"
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
#include "AliTrackPointArray.h"

#endif
#include "AliGenInfo.h" 
//
// 

ClassImp(AliTPCdigitRow)
ClassImp(AliMCInfo)
ClassImp(AliGenV0Info)
ClassImp(AliGenKinkInfo)



////////////////////////////////////////////////////////////////////////
AliMCInfo::AliMCInfo():
  fTrackRef(),
  fTrackRefOut(),
  fTRdecay(),
  fPrimPart(0),
  fParticle(),
  fMass(0),
  fCharge(0),
  fLabel(0),
  fEventNr(),
  fMCtracks(),
  fPdg(0),
  fTPCdecay(0),
  fRowsWithDigitsInn(0),
  fRowsWithDigits(0),
  fRowsTrackLength(0),
  fPrim(0),
  fTPCRow(), 
  fNTPCRef(0),                    // tpc references counter
  fNITSRef(0),                    // ITS references counter
  fNTRDRef(0),                    // TRD references counter
  fNTOFRef(0),                    // TOF references counter
  fTPCReferences(0),
  fITSReferences(0),
  fTRDReferences(0),
  fTOFReferences(0)
{
  fTPCReferences  = new TClonesArray("AliTrackReference",10);
  fITSReferences  = new TClonesArray("AliTrackReference",10);
  fTRDReferences  = new TClonesArray("AliTrackReference",10);
  fTOFReferences  = new TClonesArray("AliTrackReference",10);
  fTRdecay.SetTrack(-1);
  fCharge = 0;
}

AliMCInfo::AliMCInfo(const AliMCInfo& info):
  TObject(),
  fTrackRef(info.fTrackRef),
  fTrackRefOut(info.fTrackRefOut),
  fTRdecay(info.GetTRdecay()),
  fPrimPart(info.fPrimPart),
  fParticle(info.fParticle),
  fMass(info.GetMass()),
  fCharge(info.fCharge),
  fLabel(info.fLabel),
  fEventNr(info.fEventNr),
  fMCtracks(info.fMCtracks),
  fPdg(info.fPdg),
  fTPCdecay(info.fTPCdecay),
  fRowsWithDigitsInn(info.fRowsWithDigitsInn),
  fRowsWithDigits(info.fRowsWithDigits),
  fRowsTrackLength(info.fRowsTrackLength),
  fPrim(info.fPrim),
  fTPCRow(info.fTPCRow), 
  fNTPCRef(info.fNTPCRef),                    // tpc references counter
  fNITSRef(info.fNITSRef),                    // ITS references counter
  fNTRDRef(info.fNTRDRef),                    // TRD references counter
  fNTOFRef(info.fNTOFRef),                    // TOF references counter
  fTPCReferences(0),
  fITSReferences(0),
  fTRDReferences(0),
  fTOFReferences(0)
{
  fTPCReferences = (TClonesArray*)info.fTPCReferences->Clone();
  fITSReferences = (TClonesArray*)info.fITSReferences->Clone();
  fTRDReferences = (TClonesArray*)info.fTRDReferences->Clone();
  fTOFReferences = (TClonesArray*)info.fTOFReferences->Clone();
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
	fPrim = TPCBetheBloch(betagama);
      }else fPrim=0;
    }
  }else{
    fMass =0;
    fPrim =0;
  }  
}

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
AliTPCdigitRow::AliTPCdigitRow()
{
  Reset();
}
////////////////////////////////////////////////////////////////////////
AliTPCdigitRow & AliTPCdigitRow::operator=(const AliTPCdigitRow &digOld)
{
  for (Int_t i = 0; i<kgRowBytes; i++) fDig[i] = digOld.fDig[i];
  return (*this);
}
////////////////////////////////////////////////////////////////////////
void AliTPCdigitRow::SetRow(Int_t row) 
{
  if (row >= 8*kgRowBytes) {
    cerr<<"AliTPCdigitRow::SetRow: index "<<row<<" out of bounds."<<endl;
    return;
  }
  Int_t iC = row/8;
  Int_t iB = row%8;
  SETBIT(fDig[iC],iB);
}

////////////////////////////////////////////////////////////////////////
Bool_t AliTPCdigitRow::TestRow(Int_t row) const
{
//
// return kTRUE if row is on
//
  Int_t iC = row/8;
  Int_t iB = row%8;
  return TESTBIT(fDig[iC],iB);
}
////////////////////////////////////////////////////////////////////////
Int_t AliTPCdigitRow::RowsOn(Int_t upto) const
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
void AliTPCdigitRow::Reset()
{
//
// resets all rows to zero
//
  for (Int_t i = 0; i<kgRowBytes; i++) {
    fDig[i] <<= 8;
  }
}
////////////////////////////////////////////////////////////////////////
Int_t AliTPCdigitRow::Last() const
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
Int_t AliTPCdigitRow::First() const
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


//_____________________________________________________________________________
Float_t AliMCInfo::TPCBetheBloch(Float_t bg)
{
  //
  // Bethe-Bloch energy loss formula
  //
  const Double_t kp1=0.76176e-1;
  const Double_t kp2=10.632;
  const Double_t kp3=0.13279e-4;
  const Double_t kp4=1.8631;
  const Double_t kp5=1.9479;

  Double_t dbg = (Double_t) bg;

  Double_t beta = dbg/TMath::Sqrt(1.+dbg*dbg);

  Double_t aa = TMath::Power(beta,kp4);
  Double_t bb = TMath::Power(1./dbg,kp5);

  bb=TMath::Log(kp3+bb);
  
  return ((Float_t)((kp2-aa-bb)*kp1/aa));
}

