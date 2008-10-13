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
//     Requierements - Warnings:
//     1. Before using this componenent the magnetic filed has to be set properly //
//     2. The systematic effects  - unlinearities has to be understood
//
//     If systematic and unlinearities are not under control
//     the alignment is just effective alignment. Not second order corrction
//     are calculated.
//    
//     The histograming of the edge effects and unlineratities integral part
//     of the component (currently only in debug stream)
//
//     3 general type of linear transformation investigated (see bellow)
//
//     By default only 6 parameter alignment to be used - other just for QA purposes

//     Different linear tranformation investigated
//     12 parameters - arbitrary linear transformation 
//                     a00  a01 a02  a03     p[0]   p[1]  p[2]  p[9]
//                     a10  a11 a12  a13 ==> p[3]   p[4]  p[5]  p[10]
//                     a20  a21 a22  a23     p[6]   p[7]  p[8]  p[11] 
//
//      9 parameters - scaling fixed to 1
//                     a00  a01  a02 a03     1      p[0]  p[1]   p[6]
//                     a10  a11  a12 a13 ==> p[2]   1     p[3]   p[7]
//                     a20  a21  a22 a23     p[4]   p[5]  1      p[8] 
//
//      6 parameters - x-y rotation x-z, y-z tiliting
//                     a00  a01  a02 a03     1     -p[0]  0     p[3]
//                     a10  a11  a12 a13 ==> p[0]   1     0     p[4]
//                     a20  a21  a22 a23     p[1]   p[2]  1     p[5] 
//
//
//      Debug stream supported
//      0. Align    - The main output of the Alignment component
//                  - Used for visualization of the misalignment between sectors
//                  - Results of the missalignment fit and the mean and sigmas of histograms
//                   stored there
//      1. Tracklet - StreamLevel >1
//                  - Dump all information about tracklet match from sector1 to sector 2
//                  - Default histogram residulas created in parallel
//                  - Check this streamer in case of suspicious content of these histograms
//      2. Track    - StreamLevel>5  
//                  - For debugging of the edge effects
//                  - All information  - extrapolation inside of one sectors
//                  - Created in order to distinguish between unlinearities inside of o
//                    sector and  missalignment 
   
//
//
/*
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("align.txt","Track",0,10200);
  chain->Lookup();
  TCut cutA("abs(tp1.fP[1]-tp2.fP[1])<0.3&&abs(tp1.fP[0]-tp2.fP[0])<0.15&&abs(tp1.fP[3]-tp2.fP[3])<0.01&&abs(tp1.fP[2]-tp2.fP[2])<0.01");
  TCut cutS("s1%36==s2%36");
  
  .x ~/UliStyle.C
  .x $ALICE_ROOT/macros/loadlibsREC.C

  gSystem->Load("$ROOTSYS/lib/libXrdClient.so");
  gSystem->Load("libProof");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");
  //
  // compare reference
  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");

  AliTPCcalibAlign * align = ( AliTPCcalibAlign *)array->FindObject("alignTPC");
  //
  //
  align->MakeTree("alignTree.root");
  TFile f("alignTree.root");
  TTree * tree = (TTree*)f.Get("Align");
  

*/

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
#include "TF1.h"
#include "TGraphErrors.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTracker.h"
#include "TClonesArray.h"


#include "TTreeStream.h"
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
  //
  // 
  //
  
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
      ProcessTracklets(*common1,*common2,seed, t1->GetSector(),t2->GetSector());
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
  Info("AliTPCcalibAlign","Terminate");
  EvalFitters();
  AliTPCcalibBase::Terminate();
}




void AliTPCcalibAlign::ProcessTracklets(const AliExternalTrackParam &tp1,
					const AliExternalTrackParam &tp2,
					const AliTPCseed * seed,
					Int_t s1,Int_t s2) {

  //
  //
  //
  //
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
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
	"tp1.="<<p1<<
	"tp2.="<<p2<<
	"v1.="<<&vec1<<
	"v2.="<<&vec2<<
	"s1="<<s1<<
	"s2="<<s2<<
	"\n";
    }
  }
  //
  // Aplly cut selection 
  /*
    TChain * chainalign = tool.MakeChain("align.txt","Tracklet",0,1000000)
    chainalign->Lookup();

    // Cuts to be justified with the debug streamer
    //
    TCut c1pt("abs((tp1.fP[4]+tp2.fP[4])*0.5)<3"); // pt cut  - OK
    TCut cdy("abs(tp1.fP[0]-tp2.fP[0])<2");
    TCut cdz("abs(tp1.fP[1]-tp2.fP[1])<2");
    TCut cdphi("abs(tp1.fP[2]-tp2.fP[2])<0.02");
    TCut cdt("abs(tp1.fP[3]-tp2.fP[3])<0.02");
    TCut cd1pt("abs(tp1.fP[4]-tp2.fP[4])<0.3");    // delta 1/pt cut  -OK   
    //
    //
    TCut acut =  c1pt+cdy+cdz+cdphi+cdt+cd1pt;

  */
  //   1. pt cut
  //   2. dy
  //   3. dz
  //   4. dphi
  //   5. dtheta
  //   6. d1pt
  if (GetDebugLevel()>50) printf("Process track\n");
  if (TMath::Abs(tp1.GetParameter()[0]-tp2.GetParameter()[0])>2)    return;
  if (TMath::Abs(tp1.GetParameter()[1]-tp2.GetParameter()[1])>2)    return;
  if (TMath::Abs(tp1.GetParameter()[2]-tp2.GetParameter()[2])>0.02) return;
  if (TMath::Abs(tp1.GetParameter()[3]-tp2.GetParameter()[3])>0.02) return;
  if (TMath::Abs(tp1.GetParameter()[4]-tp2.GetParameter()[4])>0.3)  return;
  if (TMath::Abs((tp1.GetParameter()[4]+tp2.GetParameter()[4])*0.5)>3)  return;
  if (TMath::Abs((tp1.GetParameter()[0]-tp2.GetParameter()[0]))<0.000000001)  return;
   if (GetDebugLevel()>50) printf("Filling track\n");
  //
  // fill resolution histograms - previous cut included
  FillHisto(tp1,tp2,s1,s2);  
  //
  Process12(t1,t2,GetOrMakeFitter12(s1,s2));
  Process9(t1,t2,GetOrMakeFitter9(s1,s2));
  Process6(t1,t2,GetOrMakeFitter6(s1,s2));
  ProcessDiff(tp1,tp2, seed,s1,s2);
  ++fPoints[GetIndex(s1,s2)];
}

void  AliTPCcalibAlign::ProcessDiff(const AliExternalTrackParam &t1,
				    const AliExternalTrackParam &t2,
				    const AliTPCseed *seed,
				    Int_t s1,Int_t s2)
{
  //
  // Process local residuals function
  // 
  TVectorD vecX(160);
  TVectorD vecY(160);
  TVectorD vecZ(160);
  TVectorD vecClY(160);
  TVectorD vecClZ(160);
  TClonesArray arrCl("AliTPCclusterMI",160);
  arrCl.ExpandCreateFast(160);
  Int_t count1=0, count2=0;
  for (Int_t i=0;i<160;++i) {
    AliTPCclusterMI *c=seed->GetClusterPointer(i);
    vecX[i]=0;
    vecY[i]=0;
    vecZ[i]=0;
    if (!c) continue;
    AliTPCclusterMI & cl = (AliTPCclusterMI&) (*arrCl[i]);
    if (c->GetDetector()!=s1 && c->GetDetector()!=s2) continue;
    vecClY[i] = c->GetY();
    vecClZ[i] = c->GetZ();
    cl=*c;
    const AliExternalTrackParam *par = (c->GetDetector()==s1)? &t1:&t2;
    if (c->GetDetector()==s1) ++count1;
    if (c->GetDetector()==s2) ++count2;
    Double_t gxyz[3],xyz[3];
    t1.GetXYZ(gxyz);
    Float_t bz = AliTracker::GetBz(gxyz);
    par->GetYAt(c->GetX(), bz, xyz[1]);
    par->GetZAt(c->GetX(), bz, xyz[2]);
    vecX[i] = c->GetX();
    vecY[i]= xyz[1];
    vecZ[i]= xyz[2];
  }
  //
  //
  if (fStreamLevel>5){
    //
    // huge output - cluster residuals to be investigated
    //
    TTreeSRedirector *cstream = GetDebugStreamer();
    AliTPCseed * t = (AliTPCseed*) seed;
    //AliExternalTrackParam *p0 = &((AliExternalTrackParam&)seed);
    AliExternalTrackParam *p1 = &((AliExternalTrackParam&)t1);
    AliExternalTrackParam *p2 = &((AliExternalTrackParam&)t2);
    /*
      
      Track->Draw("Cl[].fY-vtY.fElements:vtY.fElements-vtX.fElements*tan(pi/18.)>>his(100,-10,0)","Cl.fY!=0&&abs(Cl.fY-vtY.fElements)<1","prof");

    */

    if (cstream){
      (*cstream)<<"Track"<<
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
	"Cl.="<<&arrCl<<
	//"tp0.="<<p0<<
	"tp1.="<<p1<<
	"tp2.="<<p2<<
	"vtX.="<<&vecX<<
	"vtY.="<<&vecY<<
	"vtZ.="<<&vecZ<<
	"vcY.="<<&vecClY<<
	"vcZ.="<<&vecClZ<<
	"s1="<<s1<<
	"s2="<<s2<<
	"c1="<<count1<<
	"c2="<<count2<<
	"\n";
    }
  }
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
  //                     a00  a01 a02  a03     p[0]   p[1]  p[2]  p[9]
  //                     a10  a11 a12  a13 ==> p[3]   p[4]  p[5]  p[10]
  //                     a20  a21 a22  a23     p[6]   p[7]  p[8]  p[11] 



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
  //                     a00  a01  a02 a03     1      p[0]  p[1]   p[6]
  //                     a10  a11  a12 a13 ==> p[2]   1     p[3]   p[7]
  //                     a20  a21  a21 a23     p[4]   p[5]  1      p[8] 


  Double_t &x1=t1[0], &y1=t1[1], &z1=t1[2], &dydx1=t1[3], &dzdx1=t1[4];
  Double_t /*&x2=t2[0],*/ &y2=t2[1], &z2=t2[2], &dydx2=t2[3], &dzdx2=t2[4];

  // TODO:
  Double_t sy    = 0.1;
  Double_t sz    = 0.1;
  Double_t sdydx = 0.001;
  Double_t sdzdx = 0.001;
  //
  Double_t p[12];
  Double_t value;

  // x2  =  a00*x1 + a01*y1 + a02*z1 + a03
  // y2  =  a10*x1 + a11*y1 + a12*z1 + a13
  // y2' =  a10*x1 + a11*y1 + a12*z1 + a13 + (a01*y1 + a02*z1 + a03)*dydx2
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[2]   += x1;           // a10
  //p[]  +=1;             // a11
  p[3]   += z1;           // a12    
  p[7]   += 1;            // a13
  p[0]   += y1*dydx2;     // a01
  p[1]   += z1*dydx2;     // a02
  p[6]   += dydx2;        // a03
  value   = y2-y1;        //-a11
  fitter->AddPoint(p,value,sy);
  //
  // x2  =  a00*x1 + a01*y1 + a02*z1 + a03
  // z2  =  a20*x1 + a21*y1 + a22*z1 + a23
  // z2' =  a20*x1 + a21*y1 + a22*z1 + a23 + (a01*y1 + a02*z1 + a03)*dzdx2;
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[4]   += x1;           // a20 
  p[5]   += y1;           // a21
  //p[]  += z1;           // a22
  p[8]   += 1.;           // a23
  p[0]   += y1*dzdx2;     // a01
  p[1]   += z1*dzdx2;     // a02
  p[6]   += dzdx2;        // a03
  value  = z2-z1;         //-a22
  fitter->AddPoint(p,value,sz);

  // dydx2 = (a10 + a11*dydx1 + a12*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a10 + a11*dydx1 + a12*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dydx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[2]   += 1.;           // a10
  //p[]  += dydx1;      // a11
  p[3]   += dzdx1;        // a12
  //p[]  += -dydx2;       // a00
  p[0]   += -dydx1*dydx2; // a01
  p[1]   += -dzdx1*dydx2; // a02
  value  = -dydx1+dydx2;  // -a11 + a00
  fitter->AddPoint(p,value,sdydx);

  // dzdx2 = (a20 + a21*dydx1 + a22*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a20 + a21*dydx1 + a22*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dzdx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[4]   += 1;            // a20
  p[5]   += dydx1;        // a21
  //p[]  += dzdx1;        // a22
  //p[]  += -dzdx2;       // a00
  p[0]   += -dydx1*dzdx2; // a01
  p[1]   += -dzdx1*dzdx2; // a02
  value  = -dzdx1+dzdx2;  // -a22 + a00
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
  //                     a00  a01  a02 a03     1     -p[0]  0     p[3]
  //                     a10  a11  a12 a13 ==> p[0]   1     0     p[4]
  //                     a20  a21  a21 a23     p[1]   p[2]  1     p[5] 

  Double_t &x1=t1[0], &y1=t1[1], &z1=t1[2], &dydx1=t1[3], &dzdx1=t1[4];
  Double_t /*&x2=t2[0],*/ &y2=t2[1], &z2=t2[2], &dydx2=t2[3], &dzdx2=t2[4];

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
  p[0]   += x1;           // a10
  //p[]  +=1;             // a11
  //p[]  += z1;           // a12    
  p[4]   += 1;            // a13
  p[0]   += -y1*dydx2;    // a01
  //p[]  += z1*dydx2;     // a02
  p[3]   += dydx2;        // a03
  value   = y2-y1;        //-a11
  fitter->AddPoint(p,value,sy);
  //
  // x2  =  a00*x1 + a01*y1 + a02*z1 + a03
  // z2  =  a20*x1 + a21*y1 + a22*z1 + a23
  // z2' =  a20*x1 + a21*y1 + a22*z1 + a23 + (a01*y1 + a02*z1 + a03)*dzdx2;
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[1]   += x1;           // a20 
  p[2]   += y1;           // a21
  //p[]  += z1;           // a22
  p[5]   += 1.;           // a23
  p[0]   += -y1*dzdx2;    // a01
  //p[]   += z1*dzdx2;     // a02
  p[3]   += dzdx2;        // a03
  value  = z2-z1;         //-a22
  fitter->AddPoint(p,value,sz);

  // dydx2 = (a10 + a11*dydx1 + a12*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a10 + a11*dydx1 + a12*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dydx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[0]   += 1.;           // a10
  //p[]  += dydx1;      // a11
  //p[]   += dzdx1;        // a12
  //p[]  += -dydx2;       // a00
  p[0]   +=  dydx1*dydx2; // a01
  //p[]   += -dzdx1*dydx2; // a02
  value  = -dydx1+dydx2;  // -a11 + a00
  fitter->AddPoint(p,value,sdydx);

  // dzdx2 = (a20 + a21*dydx1 + a22*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a20 + a21*dydx1 + a22*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dzdx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[1]   += 1;            // a20
  p[2]   += dydx1;        // a21
  //p[]  += dzdx1;        // a22
  //p[]  += -dzdx2;       // a00
  p[0]   +=  dydx1*dzdx2; // a01
  //p[]  += -dzdx1*dzdx2; // a02
  value  = -dzdx1+dzdx2;  // -a22 + a00
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
      if ((f=GetFitter12(s1,s2))&&fPoints[GetIndex(s1,s2)]>kMinPoints) {
	//	cerr<<s1<<","<<s2<<": "<<fPoints[GetIndex(s1,s2)]<<endl;
	if (f->Eval()!=0) {
	  cerr<<"Evaluation failed for "<<s1<<","<<s2<<endl;
	  f->Write(Form("f12_%d_%d",s1,s2));
	}else{
	  f->Write(Form("f12_%d_%d",s1,s2));
	}
      }
      if ((f=GetFitter9(s1,s2))&&fPoints[GetIndex(s1,s2)]>kMinPoints) {
	//	cerr<<s1<<","<<s2<<": "<<fPoints[GetIndex(s1,s2)]<<endl;
	if (f->Eval()!=0) {
	  cerr<<"Evaluation failed for "<<s1<<","<<s2<<endl;
	}else{
	  f->Write(Form("f9_%d_%d",s1,s2));
	}
      }
      if ((f=GetFitter6(s1,s2))&&fPoints[GetIndex(s1,s2)]>kMinPoints) {
	//	cerr<<s1<<","<<s2<<": "<<fPoints[GetIndex(s1,s2)]<<endl;
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
  static Int_t counter12=0;
  static TF1 f12("f12","x[0]++x[1]++x[2]++x[3]++x[4]++x[5]++x[6]++x[7]++x[8]++x[9]++x[10]++x[11]");
  TLinearFitter * fitter = GetFitter12(s1,s2);
  if (fitter) return fitter;
  //  fitter =new TLinearFitter(12,"x[0]++x[1]++x[2]++x[3]++x[4]++x[5]++x[6]++x[7]++x[8]++x[9]++x[10]++x[11]");
  fitter =new TLinearFitter(&f12,"");
  fitter->StoreData(kFALSE);
  fFitterArray12.AddAt(fitter,GetIndex(s1,s2));	
  counter12++;
  if (GetDebugLevel()>0) cerr<<"Creating fitter12 "<<s1<<","<<s2<<"  :  "<<counter12<<endl;
  return fitter;
}

TLinearFitter* AliTPCcalibAlign::GetOrMakeFitter9(Int_t s1,Int_t s2) {
  //
  //get or make fitter - general linear transformation - no scaling
  // 
  static Int_t counter9=0;
  static TF1 f9("f9","x[0]++x[1]++x[2]++x[3]++x[4]++x[5]++x[6]++x[7]++x[8]");
  TLinearFitter * fitter = GetFitter9(s1,s2);
  if (fitter) return fitter;
  //  fitter =new TLinearFitter(9,"x[0]++x[1]++x[2]++x[3]++x[4]++x[5]++x[6]++x[7]++x[8]");
  fitter =new TLinearFitter(&f9,"");
  fitter->StoreData(kFALSE);
  fFitterArray9.AddAt(fitter,GetIndex(s1,s2));
  counter9++;
  if (GetDebugLevel()>0) cerr<<"Creating fitter12 "<<s1<<","<<s2<<"  :  "<<counter9<<endl;
  return fitter;
}

TLinearFitter* AliTPCcalibAlign::GetOrMakeFitter6(Int_t s1,Int_t s2) {
  //
  // get or make fitter  - 6 paramater linear tranformation
  //                     - no scaling
  //                     - rotation x-y
  //                     - tilting x-z, y-z
  static Int_t counter6=0;
  static TF1 f6("f6","x[0]++x[1]++x[2]++x[3]++x[4]++x[5]");
  TLinearFitter * fitter = GetFitter6(s1,s2);
  if (fitter) return fitter;
  //  fitter=new TLinearFitter(6,"x[0]++x[1]++x[2]++x[3]++x[4]++x[5]");
  fitter=new TLinearFitter(&f6,"");
  fitter->StoreData(kFALSE);
  fFitterArray6.AddAt(fitter,GetIndex(s1,s2));
  counter6++;
  if (GetDebugLevel()>0) cerr<<"Creating fitter6 "<<s1<<","<<s2<<"  :  "<<counter6<<endl;
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
    a[0][0]=1;    a[0][1]=p[0]; a[0][2]=p[1]; a[0][3]=p[6];
    a[1][0]=p[2]; a[1][1]=1;    a[1][2]=p[3]; a[1][3]=p[7];
    a[2][0]=p[4]; a[2][1]=p[5]; a[2][2]=1;    a[2][3]=p[8];
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
    a[0][0]=1;       a[0][1]=-p[0];a[0][2]=0;    a[0][3]=p[3];
    a[1][0]=p[0];    a[1][1]=1;    a[1][2]=0;    a[1][3]=p[4];
    a[2][0]=p[1];    a[2][1]=p[2]; a[2][2]=1;    a[2][3]=p[5];
    a[3][0]=0.;      a[3][1]=0.;   a[3][2]=0.;   a[3][3]=1.;
    return true;
  } 
}

void AliTPCcalibAlign::FillHisto(const AliExternalTrackParam &tp1,
					const AliExternalTrackParam &tp2,
					Int_t s1,Int_t s2) {
  //
  // Fill residual histograms
  // Innner-Outer
  // Left right - x-y
  // A-C side
  if (TMath::Abs(s2%36-s1%36)<2 || TMath::Abs(s2%18-s1%18)==0)  {  
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

TGraphErrors * AliTPCcalibAlign::MakeGraph(Int_t sec0, Int_t sec1, Int_t dsec, 
					   Int_t i0, Int_t i1, FitType type) 
{
  //
  //
  //
  TMatrixD mat;
  TObjArray *fitArray=0;
  Double_t xsec[1000];
  Double_t ysec[1000];
  Int_t npoints=0;
  for (Int_t isec = sec0; isec<=sec1; isec++){
    Int_t isec2 = (isec+dsec)%72;    
    switch (type) {
    case k6:
      GetTransformation6(isec,isec2,mat);break;
    case k9:
      GetTransformation9(isec,isec2,mat);break;
    case k12:
      GetTransformation12(isec,isec2,mat);break;
    }
    xsec[npoints]=isec;
    ysec[npoints]=mat(i0,i1);
    ++npoints;
  }
  TGraphErrors *gr = new TGraphErrors(npoints,xsec,ysec,0,0);
  Char_t name[1000];
  sprintf(name,"Mat[%d,%d]  Type=%d",i0,i1,type);
  gr->SetName(name);
  return gr;
}

void  AliTPCcalibAlign::MakeTree(const char *fname){
  //
  // make tree with alignment cosntant  -
  // For  QA visualization
  //
  /*
   TFile f("CalibObjects.root");
   TObjArray *array  = (TObjArray*)f.Get("TPCCalib");
   AliTPCcalibAlign   *alignTPC = (AliTPCcalibAlign   *)array->At(0);
   alignTPC->MakeTree("alignTree.root");
   TFile falign("alignTree.root");
   Align->Draw("dy")
   */
  const Int_t kMinPoints=50;
  TTreeSRedirector cstream(fname);
  for (Int_t s1=0;s1<72;++s1)
    for (Int_t s2=0;s2<72;++s2){
      if (fPoints[GetIndex(s1,s2)]<kMinPoints) continue;
      TMatrixD m6;
      TMatrixD m9;
      TMatrixD m12;
      GetTransformation6(s1,s2,m6);
      GetTransformation9(s1,s2,m9);
      GetTransformation12(s1,s2,m12);
      Double_t dy=0, dz=0, dphi=0,dtheta=0;
      Double_t sy=0, sz=0, sphi=0,stheta=0;
      Double_t ny=0, nz=0, nphi=0,ntheta=0;
      TH1 * his=0;
      his = GetHisto(kY,s1,s2);
      if (his) { dy = his->GetMean(); sy = his->GetRMS(); ny = his->GetEntries();}
      his = GetHisto(kZ,s1,s2);
      if (his) { dz = his->GetMean(); sz = his->GetRMS(); nz = his->GetEntries();}
      his = GetHisto(kPhi,s1,s2);
      if (his) { dphi = his->GetMean(); sphi = his->GetRMS(); nphi = his->GetEntries();}
      his = GetHisto(kTheta,s1,s2);
      if (his) { dtheta = his->GetMean(); stheta = his->GetRMS(); ntheta = his->GetEntries();}
      //
      cstream<<"Align"<<
	"s1="<<s1<<     // reference sector
	"s2="<<s2<<     // sector to align
	"m6.="<<&m6<<   // tranformation matrix
	"m9.="<<&m9<<   // 
	"m12.="<<&m12<<
	//               histograms mean RMS and entries
	"dy="<<dy<<  
	"sy="<<sy<<
	"ny="<<ny<<
	"dz="<<dz<<
	"sz="<<sz<<
	"nz="<<nz<<
	"dphi="<<dphi<<
	"sphi="<<sphi<<
	"nphi="<<nphi<<
	"dtheta="<<dtheta<<
	"stheta="<<stheta<<
	"ntheta="<<ntheta<<
	"\n";
    }

}


//_____________________________________________________________________
Long64_t AliTPCcalibAlign::Merge(TCollection* list) {
  //
  // merge function 
  //
  if (GetDebugLevel()>0) Info("AliTPCcalibAlign","Merge");
  if (!list)
    return 0;  
  if (list->IsEmpty())
    return 1;
  
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;  
  iter->Reset();
  Int_t count=0;
  while((obj = iter->Next()) != 0)
    {
      AliTPCcalibAlign* entry = dynamic_cast<AliTPCcalibAlign*>(obj);
      if (entry == 0) continue; 
      Add(entry);
      count++;
    } 
  return count;
}


void AliTPCcalibAlign::Add(AliTPCcalibAlign * align){
  //
  // Add entry
  //

  for (Int_t i=0; i<72;i++){
    for (Int_t j=0; j<72;j++){
      fPoints[GetIndex(i,j)]+=align->fPoints[GetIndex(i,j)];

      //
      // dy
      TH1* hdy0 = GetHisto(kY,i,j);
      TH1* hdy1 = align->GetHisto(kY,i,j);
      if (hdy1){
	if (hdy0) hdy0->Add(hdy1);
	else {
	  hdy0 = GetHisto(kY,i,j,kTRUE);
	  hdy0->Add(hdy1);
	}
      }      
      //
      // dz
      TH1* hdz0 = GetHisto(kZ,i,j);
      TH1* hdz1 = align->GetHisto(kZ,i,j);
      if (hdz1){
	if (hdz0) hdz0->Add(hdz1);
	else {
	  hdz0 = GetHisto(kZ,i,j,kTRUE);
	  hdz0->Add(hdz1);
	}
      }
      //
      // dphi
      TH1* hdphi0 = GetHisto(kPhi,i,j);
      TH1* hdphi1 = align->GetHisto(kPhi,i,j);
      if (hdphi1){
	if (hdphi0) hdphi0->Add(hdphi1);
	else {
	  hdphi0 = GetHisto(kPhi,i,j,kTRUE);
	  hdphi0->Add(hdphi1);
	}
      }      
      //
      // dtheta
      TH1* hdTheta0 = GetHisto(kTheta,i,j);
      TH1* hdTheta1 = align->GetHisto(kTheta,i,j);
      if (hdTheta1){
	if (hdTheta0) hdTheta0->Add(hdTheta1);
	else {
	  hdTheta0 = GetHisto(kTheta,i,j,kTRUE);
	  hdTheta0->Add(hdTheta1);
	}
      }           
    }
  }
  TLinearFitter *f0=0;
  TLinearFitter *f1=0;
  for (Int_t i=0; i<72;i++){
    for (Int_t j=0; j<72;j++){
      //
      // fitter12
      f0 =  GetFitter12(i,j);
      f1 =  GetFitter12(i,j);
      if (f1){
	if (f0) f0->Add(f1);
	else {
	  f0 = GetOrMakeFitter12(i,j);
	  f0->Add(f1);
	}
      }      
      //
      // fitter9
      f0 =  GetFitter9(i,j);
      f1 =  GetFitter9(i,j);
      if (f1){
	if (f0) f0->Add(f1);
	else {
	  f0 = GetOrMakeFitter9(i,j);
	  f0->Add(f1);
	}
      }      
      f0 =  GetFitter6(i,j);
      f1 =  GetFitter6(i,j);
      if (f1){
	if (f0) f0->Add(f1);
	else {
	  f0 = GetOrMakeFitter6(i,j);
	  f0->Add(f1);
	}
      }   
    }
  }
}





/*
  

gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
AliXRDPROOFtoolkit tool;
TChain * chainTr = tool.MakeChain("align.txt","Track",0,10200);
chainTr->Lookup();



TCut cutS("s1%36==s2%36");

TCut cutN("c1>32&&c2>60");
TCut cutC0("sqrt(tp2.fC[0])<1");

TCut cutP0("abs(tp1.fP[0]-tp2.fP[0])<0.4");
TCut cutP2("abs(tp1.fP[2]-tp2.fP[2])<0.01");
TCut cutP3("abs(tp1.fP[3]-tp2.fP[3])<0.01");
TCut cutP4("abs(tp1.fP[4]-tp2.fP[4])<0.5");
TCut cutP=cutP0+cutP2+cutP3+cutP4+cutC0;

TCut cutX("abs(tp2.fX-133.6)<2");

TCut cutA = cutP+cutN;


TCut cutY("abs(vcY.fElements-vtY.fElements)<0.3&&vcY.fElements!=0")
TCut cutZ("abs(vcZ.fElements-vtZ.fElements)<0.3&&vcZ.fElements!=0")

*/
