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
//     Class to make a internal alignemnt of TPC chambers                    //
//
//     Different linear tranformation investigated
//     12 parameters - arbitrary linear transformation 
//      9 parameters - scaling fixed to 1
//      6 parameters - x-y rotation x-z, y-z tiliting
////
//// Only the 12 parameter version works so far...
////
////
//// 

#include "TLinearFitter.h"
#include "AliTPCcalibAlign.h"
#include "AliExternalTrackParam.h"
#include "AliTPCTracklet.h"
#include "TH1D.h"
#include "TVectorD.h"
#include "TTreeStream.h"
#include "TFile.h"

#include <iostream>
#include <sstream>
using namespace std;

ClassImp(AliTPCcalibAlign)

AliTPCcalibAlign::AliTPCcalibAlign()
  :  fDphiHistArray(72*72),
     fDthetaHistArray(72*72),
     fDyHistArray(72*72),
     fDzHistArray(72*72),
     fFitterArray12(72*72),
     fFitterArray9(72*72),
     fFitterArray6(72*72)
{
  //
  // Constructor
  //
  for (Int_t i=0;i<72*72;++i) {
    fPoints[i]=0;
  }
}

AliTPCcalibAlign::AliTPCcalibAlign(const Text_t *name, const Text_t *title)
  :AliTPCcalibBase(),  
   fDphiHistArray(72*72),
   fDthetaHistArray(72*72),
   fDyHistArray(72*72),
   fDzHistArray(72*72),
   fFitterArray12(72*72),
   fFitterArray9(72*72),
   fFitterArray6(72*72)
{
  //
  // Constructor
  //
  SetName(name);
  SetTitle(title);
  for (Int_t i=0;i<72*72;++i) {
    fPoints[i]=0;
  }
}

AliTPCcalibAlign::~AliTPCcalibAlign() {
  //
  // destructor
  //
}

void AliTPCcalibAlign::Process(AliTPCseed *seed) {
  TObjArray tracklets=
    AliTPCTracklet::CreateTracklets(seed,AliTPCTracklet::kKalman,
				    kFALSE,20,2);
 //  TObjArray trackletsL=
//     AliTPCTracklet::CreateTracklets(seed,AliTPCTracklet::kLinear,
// 				    kFALSE,20,2);
//   TObjArray trackletsQ=
//     AliTPCTracklet::CreateTracklets(seed,AliTPCTracklet::kQuadratic,
// 				    kFALSE,20,2);
//   TObjArray trackletsR=
//     AliTPCTracklet::CreateTracklets(seed,AliTPCTracklet::kRiemann,
// 				    kFALSE,20,2);
  tracklets.SetOwner();
 //  trackletsL.SetOwner();
//   trackletsQ.SetOwner();
//   trackletsR.SetOwner();
  if (tracklets.GetEntries()==2) {
    AliTPCTracklet *t1=static_cast<AliTPCTracklet*>(tracklets[0]);
    AliTPCTracklet *t2=static_cast<AliTPCTracklet*>(tracklets[1]);
    if (t1->GetSector()>t2->GetSector()) {
      AliTPCTracklet* tmp=t1;
      t1=t2;
      t2=tmp;
    }
    AliExternalTrackParam *common1=0,*common2=0;
    if (AliTPCTracklet::PropagateToMeanX(*t1,*t2,common1,common2))
      ProcessTracklets(*common1,*common2,t1->GetSector(),t2->GetSector());
    delete common1;
    delete common2;
  }
  
}

void AliTPCcalibAlign::Analyze(){
  //
  // Analyze function 
  //
  EvalFitters();
}


void AliTPCcalibAlign::Terminate(){
  //
  // Terminate function
  // call base terminate + Eval of fitters
  //
  if (GetDebugLevel()>0) Info("AliTPCcalibAlign","Teminate");
  EvalFitters();
  AliTPCcalibBase::Terminate();
}




void AliTPCcalibAlign::ProcessTracklets(const AliExternalTrackParam &tp1,
					const AliExternalTrackParam &tp2,
					Int_t s1,Int_t s2) {

  //
  //
  //
  FillHisto(tp1,tp2,s1,s2);
  //
  // Process function to fill fitters
  //
  Double_t t1[5],t2[5];
  Double_t &x1=t1[0], &y1=t1[1], &z1=t1[2], &dydx1=t1[3], &dzdx1=t1[4];
  Double_t &x2=t2[0], &y2=t2[1], &z2=t2[2], &dydx2=t2[3], &dzdx2=t2[4];
  x1   =tp1.GetX();
  y1   =tp1.GetY();
  z1   =tp1.GetZ();
  Double_t snp1=tp1.GetSnp();
  dydx1=snp1/TMath::Sqrt(1.-snp1*snp1);
  Double_t tgl1=tp1.GetTgl();
  // dz/dx = 1/(cos(theta)*cos(phi))
  dzdx1=tgl1/TMath::Sqrt(1.-snp1*snp1);
  x2   =tp2.GetX();
  y2   =tp2.GetY();
  z2   =tp2.GetZ();
  Double_t snp2=tp2.GetSnp();
  dydx2=snp2/TMath::Sqrt(1.-snp2*snp2);
  Double_t tgl2=tp2.GetTgl();
  dzdx2=tgl2/TMath::Sqrt(1.-snp2*snp2);
  //
  //
  //
  if (fStreamLevel>1){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      static TVectorD vec1(5);
      static TVectorD vec2(5);
      vec1.SetElements(t1);
      vec2.SetElements(t2);
      AliExternalTrackParam *p1 = &((AliExternalTrackParam&)tp1);
      AliExternalTrackParam *p2 = &((AliExternalTrackParam&)tp2);
      (*cstream)<<"Tracklet"<<
	"tp1.="<<p1<<
	"tp2.="<<p2<<
	"v1.="<<&vec1<<
	"v2.="<<&vec2<<
	"s1="<<s1<<
	"s2="<<s2<<
	"\n";
    }
  }
  Process12(t1,t2,GetOrMakeFitter12(s1,s2));
  Process9(t1,t2,GetOrMakeFitter9(s1,s2));
  Process6(t1,t2,GetOrMakeFitter6(s1,s2));
  ++fPoints[72*s1+s2];
}

void AliTPCcalibAlign::Process12(const Double_t *t1,
				 const Double_t *t2,
				 TLinearFitter *fitter) {
  // x2    =  a00*x1 + a01*y1 + a02*z1 + a03
  // y2    =  a10*x1 + a11*y1 + a12*z1 + a13
  // z2    =  a20*x1 + a21*y1 + a22*z1 + a23
  // dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  //
  //       a00  a01 a02  a03     p[0]   p[1]  p[2]  p[9]
  //       a10  a11 a12  a13 ==> p[3]   p[4]  p[5]  p[10]
  //       a20  a21 a22  a23     p[6]   p[7]  p[8]  p[11] 
  const Double_t &x1=t1[0], &y1=t1[1], &z1=t1[2], &dydx1=t1[3], &dzdx1=t1[4];
  const Double_t /*&x2=t2[0],*/ &y2=t2[1], &z2=t2[2], &dydx2=t2[3], &dzdx2=t2[4];

  // TODO:
  Double_t sy    = 0.1;
  Double_t sz    = 0.1;
  Double_t sdydx = 0.001;
  Double_t sdzdx = 0.001;

  Double_t p[12];
  Double_t value;

  // x2  =  a00*x1 + a01*y1 + a02*z1 + a03
  // y2  =  a10*x1 + a11*y1 + a12*z1 + a13
  // y2' =  a10*x1 + a11*y1 + a12*z1 + a13 + (a01*y1 + a02*z1 + a03)*dydx2
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[3+0] = x1;          // a10
  p[3+1] = y1;          // a11
  p[3+2] = z1;          // a12
  p[9+1] = 1.;          // a13
  p[0+1] = y1*dydx2;    // a01
  p[0+2] = z1*dydx2;    // a02
  p[9+0] = dydx2;       // a03
  value  = y2;
  fitter->AddPoint(p,value,sy);

  // x2  =  a00*x1 + a01*y1 + a02*z1 + a03
  // z2  =  a20*x1 + a21*y1 + a22*z1 + a23
  // z2' =  a20*x1 + a21*y1 + a22*z1 + a23 + (a01*y1 + a02*z1 + a03)*dzdx2;
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[6+0] = x1;           // a20 
  p[6+1] = y1;           // a21
  p[6+2] = z1;           // a22
  p[9+2] = 1.;           // a23
  p[0+1] = y1*dzdx2;     // a01
  p[0+2] = z1*dzdx2;     // a02
  p[9+0] = dzdx2;        // a03
  value  = z2;
  fitter->AddPoint(p,value,sz);

  // dydx2 = (a10 + a11*dydx1 + a12*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a10 + a11*dydx1 + a12*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dydx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[3+0] = 1.;           // a10
  p[3+1] = dydx1;        // a11
  p[3+2] = dzdx1;        // a12
  p[0+0] = -dydx2;       // a00
  p[0+1] = -dydx1*dydx2; // a01
  p[0+2] = -dzdx1*dydx2; // a02
  value  = 0.;
  fitter->AddPoint(p,value,sdydx);

  // dzdx2 = (a20 + a21*dydx1 + a22*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a20 + a21*dydx1 + a22*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dzdx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[6+0] = 1;            // a20
  p[6+1] = dydx1;        // a21
  p[6+2] = dzdx1;        // a22
  p[0+0] = -dzdx2;       // a00
  p[0+1] = -dydx1*dzdx2; // a01
  p[0+2] = -dzdx1*dzdx2; // a02
  value  = 0.;
  fitter->AddPoint(p,value,sdzdx);
}

void AliTPCcalibAlign::Process9(Double_t *t1,
				Double_t *t2,
				TLinearFitter *fitter) {
  // x2    =  a00*x1 + a01*y1 + a02*z1 + a03
  // y2    =  a10*x1 + a11*y1 + a12*z1 + a13
  // z2    =  a20*x1 + a21*y1 + a22*z1 + a23
  // dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  //
  //       a00  a01 a02  a03     p[0]   p[1]  p[2]  p[9]
  //       a10  a11 a12  a13 ==> p[3]   p[4]  p[5]  p[10]
  //       a20  a21 a22  a23     p[6]   p[7]  p[8]  p[11] 
  Double_t &x1=t1[0], &y1=t1[1], &z1=t1[2], &dydx1=t1[3], &dzdx1=t1[4];
  Double_t /*&x2=t2[0],*/ &y2=t2[1], &z2=t2[2], &dydx2=t2[3], &dzdx2=t2[4];

  // TODO:
  Double_t sy    = 1.;
  Double_t sz    = 1.;
  Double_t sdydx = 1.;
  Double_t sdzdx = 1.;

  Double_t p[12];
  Double_t value;

  // x2  =  a00*x1 + a01*y1 + a02*z1 + a03
  // y2  =  a10*x1 + a11*y1 + a12*z1 + a13
  // y2' =  a10*x1 + a11*y1 + a12*z1 + a13 + (a01*y1 + a03)*dydx2
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[3+0] = x1;
  p[3+1] = y1;
  p[3+2] = z1;
  p[9+1] = 1.;
  p[0+1] = y1*dydx2;
  p[0+2] = z1*dydx2;
  p[9+0] = dydx2;
  value  = y2;
  fitter->AddPoint(p,value,sy);

  // x2  =  a00*x1 + a01*y1 + a02*z1 + a03
  // z2  =  a20*x1 + a21*y1 + a22*z1 + a23
  // z2' =  a20*x1 + a21*y1 + a22*z1 + a23 + (a01*y1 + a03)*dzdx2;
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[6+0] = x1;
  p[6+1] = y1;
  p[6+2] = z1;
  p[9+2] = 1.;
  p[0+1] = y1*dzdx2;
  p[0+2] = z1*dzdx2;
  p[9+0] = dzdx2;
  value  = z2;
  fitter->AddPoint(p,value,sz);

  // dydx2 = (a10 + a11*dydx1 + a12*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a10 + a11*dydx1 + a12*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dydx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[3+0] = 1.;
  p[3+1] = dydx1;
  p[3+2] = dzdx1;
  p[0+0] = -dydx2;
  p[0+1] = -dydx1*dydx2;
  p[0+2] = -dzdx1*dydx2;
  value  = 0.;
  fitter->AddPoint(p,value,sdydx);

  // dzdx2 = (a20 + a21*dydx1 + a22*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a20 + a21*dydx1 + a22*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dzdx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[6+0] = 1;
  p[6+1] = dydx1;
  p[6+2] = dzdx1;
  p[0+0] = -dzdx2;
  p[0+1] = -dydx1*dzdx2;
  p[0+2] = -dzdx1*dzdx2;
  value  = 0.;
  fitter->AddPoint(p,value,sdzdx);
}

void AliTPCcalibAlign::Process6(Double_t *t1,
				Double_t *t2,
				TLinearFitter *fitter) {
  // x2    =  1  *x1 +-a01*y1 + 0      +a03
  // y2    =  a01*x1 + 1  *y1 + 0      +a13
  // z2    =  a20*x1 + a21*y1 + 1  *z1 +a23
  // dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  //
  //       1   -a01 0    a03     x     -p[0]  x     p[3]
  //       a10  1   0    a13 ==> p[0]   x     x     p[4]
  //       a20  a21 1    a23     p[1]   p[2]  x     p[5] 
  Double_t &x1=t1[0], &y1=t1[1], &z1=t1[2], &dydx1=t1[3], &dzdx1=t1[4];
  Double_t /*&x2=t2[0],*/ &y2=t2[1], &z2=t2[2], &dydx2=t2[3], &dzdx2=t2[4];

  // TODO:
  Double_t sy    = 0.1;
  Double_t sz    = 0.1;
  Double_t sdydx = 0.001;
  Double_t sdzdx = 0.001;

  Double_t p[12];
  Double_t value;

  // x2  =  1  *x1 +-a01*y1 + 0      +a03
  // y2  =  a01*x1 + 1  *y1 + 0      +a13
  // y2' =  a10*x1 + a11*y1 + a12*z1 + a13 + (a01*y1 + a02*z1 + a03)*dydx2
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[3+0] = x1;
  p[3+1] = y1;
  p[3+2] = z1;
  p[9+1] = 1.;
  p[0+1] = y1*dydx2;
  p[0+2] = z1*dydx2;
  p[9+0] = dydx2;
  value  = y2;
  fitter->AddPoint(p,value,sy);

  // x2  =  1  *x1 +-a01*y1 + 0      + a03
  // z2  =  a20*x1 + a21*y1 + 1  *z1 + a23
  // z2' =  a20*x1 + a21*y1 + a22*z1 + a23 + (a01*y1 + a02*z1 +a03)*dzdx2;
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[6+0] = x1;
  p[6+1] = y1;
  p[6+2] = z1;
  p[9+2] = 1.;
  p[0+1] = y1*dzdx2;
  p[0+2] = z1*dzdx2;
  p[9+0] = dzdx2;
  value  = z2;
  fitter->AddPoint(p,value,sz);

  // dydx2 = (a10 + a11*dydx1 + a12*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a10 + a11*dydx1 + a12*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dydx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[3+0] = 1.;
  p[3+1] = dydx1;
  p[3+2] = dzdx1;
  p[0+0] = -dydx2;
  p[0+1] = -dydx1*dydx2;
  p[0+2] = -dzdx1*dydx2;
  value  = 0.;
  fitter->AddPoint(p,value,sdydx);

  // dzdx2 = (a20 + a21*dydx1 + a22*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a20 + a21*dydx1 + a22*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dzdx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[6+0] = 1;
  p[6+1] = dydx1;
  p[6+2] = dzdx1;
  p[0+0] = -dzdx2;
  p[0+1] = -dydx1*dzdx2;
  p[0+2] = -dzdx1*dzdx2;
  value  = 0.;
  fitter->AddPoint(p,value,sdzdx);
}




void AliTPCcalibAlign::EvalFitters() {
  //
  // Analyze function 
  // 
  // Perform the fitting using linear fitters
  //
  Int_t kMinPoints =50;
  TLinearFitter *f;
  TFile fff("alignDebug.root","recreate");
  for (Int_t s1=0;s1<72;++s1)
    for (Int_t s2=0;s2<72;++s2){
      if ((f=GetFitter12(s1,s2))&&fPoints[72*s1+s2]>kMinPoints) {
	//	cerr<<s1<<","<<s2<<": "<<fPoints[72*s1+s2]<<endl;
	if (f->Eval()!=0) {
	  cerr<<"Evaluation failed for "<<s1<<","<<s2<<endl;
	  f->Write(Form("f12_%d_%d",s1,s2));
	}else{
	  f->Write(Form("f12_%d_%d",s1,s2));
	}
      }
      if ((f=GetFitter9(s1,s2))&&fPoints[72*s1+s2]>kMinPoints) {
	//	cerr<<s1<<","<<s2<<": "<<fPoints[72*s1+s2]<<endl;
	if (f->Eval()!=0) {
	  cerr<<"Evaluation failed for "<<s1<<","<<s2<<endl;
	}else{
	  f->Write(Form("f9_%d_%d",s1,s2));
	}
      }
      if ((f=GetFitter6(s1,s2))&&fPoints[72*s1+s2]>kMinPoints) {
	//	cerr<<s1<<","<<s2<<": "<<fPoints[72*s1+s2]<<endl;
	if (f->Eval()!=0) {
	  cerr<<"Evaluation failed for "<<s1<<","<<s2<<endl;
	}else{
	  f->Write(Form("f6_%d_%d",s1,s2));
	}
      }
    }
  this->Write("align");
  /*
		    
  fitter->Eval();
  fitter->Eval();
  chi212 = align->GetChisquare()/(4.*entries);

  TMatrixD mat(13,13);
  TVectorD par(13);
  align->GetParameters(par);
  align->GetCovarianceMatrix(mat);

  //
  //
  for (Int_t i=0; i<12;i++){
    palign12(i)= par(i+1);
    for (Int_t j=0; j<12;j++){
      pcovar12(i,j)   = mat(i+1,j+1);
      pcovar12(i,j) *= chi212;
    }
  }
  //
  for (Int_t i=0; i<12;i++){
    psigma12(i)  = TMath::Sqrt(pcovar12(i,i));
    palignR12(i) = palign12(i)/TMath::Sqrt(pcovar12(i,i));
    for (Int_t j=0; j<12;j++){
      pcovarN12(i,j) = pcovar12(i,j)/TMath::Sqrt(pcovar12(i,i)*pcovar12(j,j));
    }
  }
  */
}

TLinearFitter* AliTPCcalibAlign::GetOrMakeFitter12(Int_t s1,Int_t s2) {
  //
  // get or make fitter - general linear transformation
  //
  TLinearFitter * fitter = GetFitter12(s1,s2);
  if (fitter) return fitter;
  fitter =new TLinearFitter(12,"x[0]++x[1]++x[2]++x[3]++x[4]++x[5]++x[6]++x[7]++x[8]++x[9]++x[10]++x[11]");
  fitter->StoreData(kFALSE);
  fFitterArray12.AddAt(fitter,GetIndex(s1,s2));
  return fitter;
}

TLinearFitter* AliTPCcalibAlign::GetOrMakeFitter9(Int_t s1,Int_t s2) {
  //
  //get or make fitter - general linear transformation - no scaling
  // 
  TLinearFitter * fitter = GetFitter9(s1,s2);
  if (fitter) return fitter;
  fitter =new TLinearFitter(9,"x[0]++x[1]++x[2]++x[3]++x[4]++x[5]++x[6]++x[7]++x[8]");
  fitter->StoreData(kFALSE);
  fFitterArray9.AddAt(fitter,GetIndex(s1,s2));
  return fitter;
}

TLinearFitter* AliTPCcalibAlign::GetOrMakeFitter6(Int_t s1,Int_t s2) {
  //
  // get or make fitter  - 6 paramater linear tranformation
  //                     - no scaling
  //                     - rotation x-y
  //                     - tilting x-z, y-z
  TLinearFitter * fitter = GetFitter6(s1,s2);
  if (fitter) return fitter;
  fitter=new TLinearFitter(6,"x[0]++x[1]++x[2]++x[3]++x[4]++x[5]");
  fitter->StoreData(kFALSE);
  fFitterArray6.AddAt(fitter,GetIndex(s1,s2));
  return fitter;
}





Bool_t AliTPCcalibAlign::GetTransformation12(Int_t s1,Int_t s2,TMatrixD &a) {
  //
  // GetTransformation matrix - 12 paramaters - generael linear transformation
  //
  if (!GetFitter12(s1,s2))
    return false;
  else {
    TVectorD p(12);
    GetFitter12(s1,s2)->GetParameters(p);
    a.ResizeTo(4,4);
    a[0][0]=p[0]; a[0][1]=p[1]; a[0][2]=p[2]; a[0][3]=p[9];
    a[1][0]=p[3]; a[1][1]=p[4]; a[1][2]=p[5]; a[1][3]=p[10];
    a[2][0]=p[6]; a[2][1]=p[7]; a[2][2]=p[8]; a[2][3]=p[11];
    a[3][0]=0.;   a[3][1]=0.;   a[3][2]=0.;   a[3][3]=1.;
    return true;
  } 
}

Bool_t AliTPCcalibAlign::GetTransformation9(Int_t s1,Int_t s2,TMatrixD &a) {
  //
  // GetTransformation matrix - 9 paramaters - general linear transformation
  //                            No scaling
  //
  if (!GetFitter9(s1,s2))
    return false;
  else {
    TVectorD p(9);
    GetFitter9(s1,s2)->GetParameters(p);
    a.ResizeTo(4,4);
    a[0][0]=p[0]; a[0][1]=p[1]; a[0][2]=p[2]; a[0][3]=p[9];
    a[1][0]=p[3]; a[1][1]=p[4]; a[1][2]=p[5]; a[1][3]=p[10];
    a[2][0]=p[6]; a[2][1]=p[7]; a[2][2]=p[8]; a[2][3]=p[11];
    a[3][0]=0.;   a[3][1]=0.;   a[3][2]=0.;   a[3][3]=1.;
    return true;
  } 
}

Bool_t AliTPCcalibAlign::GetTransformation6(Int_t s1,Int_t s2,TMatrixD &a) {
  //
  // GetTransformation matrix - 6  paramaters
  //                            3  translation
  //                            1  rotation -x-y  
  //                            2  tilting x-z y-z
  if (!GetFitter6(s1,s2))
    return false;
  else {
    TVectorD p(6);
    GetFitter6(s1,s2)->GetParameters(p);
    a.ResizeTo(4,4);
    a[0][0]=p[0];    a[0][1]=p[1]; a[0][2]=p[2]; a[0][3]=p[9];
    a[1][0]=p[3];    a[1][1]=p[4]; a[1][2]=p[5]; a[1][3]=p[10];
    a[2][0]=p[6];    a[2][1]=p[7]; a[2][2]=p[8]; a[2][3]=p[11];
    a[3][0]=0.;      a[3][1]=0.;   a[3][2]=0.;   a[3][3]=1.;
    return true;
  } 
}

void AliTPCcalibAlign::FillHisto(const AliExternalTrackParam &tp1,
					const AliExternalTrackParam &tp2,
					Int_t s1,Int_t s2) {
  //
  // Fill residual histograms
  //
  if (s2-s1==36) {//only inner-outer
    GetHisto(kPhi,s1,s2,kTRUE)->Fill(TMath::ASin(tp1.GetSnp())-TMath::ASin(tp2.GetSnp()));    
    GetHisto(kTheta,s1,s2,kTRUE)->Fill(TMath::ATan(tp1.GetTgl())-TMath::ATan(tp2.GetTgl()));
    GetHisto(kY,s1,s2,kTRUE)->Fill(tp1.GetY()-tp2.GetY());
    GetHisto(kZ,s1,s2,kTRUE)->Fill(tp1.GetZ()-tp2.GetZ());
  }  
}



TH1 * AliTPCcalibAlign::GetHisto(HistoType type, Int_t s1, Int_t s2, Bool_t force)
{
  //
  // return specified residual histogram - it is only QA
  // if force specified the histogram and given histogram is not existing 
  //  new histogram is created
  //
  if (GetIndex(s1,s2)>=72*72) return 0;
  TObjArray *histoArray=0;
  switch (type) {
  case kY:
    histoArray = &fDyHistArray; break;
  case kZ:
    histoArray = &fDzHistArray; break;
  case kPhi:
    histoArray = &fDphiHistArray; break;
  case kTheta:
    histoArray = &fDthetaHistArray; break;
  }
  TH1 * histo= (TH1*)histoArray->At(GetIndex(s1,s2));
  if (histo) return histo;
  if (force==kFALSE) return 0; 
  //
  stringstream name;
  stringstream title;
  switch (type) {    
  case kY:
    name<<"hist_y_"<<s1<<"_"<<s2;
    title<<"Y Missalignment for sectors "<<s1<<" and "<<s2;
    histo =new TH1D(name.str().c_str(),title.str().c_str(),512,-0.3,0.3); // +/- 3 mm
    break;
  case kZ:
    name<<"hist_z_"<<s1<<"_"<<s2;
    title<<"Z Missalignment for sectors "<<s1<<" and "<<s2;
    histo = new TH1D(name.str().c_str(),title.str().c_str(),512,-0.3,0.3); // +/- 3 mm
    break;
  case kPhi:
    name<<"hist_phi_"<<s1<<"_"<<s2;
    title<<"Phi Missalignment for sectors "<<s1<<" and "<<s2;
    histo =new TH1D(name.str().c_str(),title.str().c_str(),512,-0.01,0.01); // +/- 10 mrad
    break;
  case kTheta:
    name<<"hist_theta_"<<s1<<"_"<<s2;
    title<<"Theta Missalignment for sectors "<<s1<<" and "<<s2;
    histo =new TH1D(name.str().c_str(),title.str().c_str(),512,-0.01,0.01); // +/- 10 mrad
    break;
  }
  histo->SetDirectory(0);
  histoArray->AddAt(histo,GetIndex(s1,s2));
  return histo;
}
