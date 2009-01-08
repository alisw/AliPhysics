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
  align->EvalFitters();
  align->MakeTree("alignTree.root");
  TFile falignTree("alignTree.root");
  TTree * treeAlign = (TTree*)falignTree.Get("Align");
   

*/

////
//// 

#include "TLinearFitter.h"
#include "AliTPCcalibAlign.h"
#include "AliExternalTrackParam.h"
#include "AliTPCTracklet.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TVectorD.h"
#include "TTreeStream.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTracker.h"
#include "TClonesArray.h"
#include "AliExternalComparison.h"


#include "TTreeStream.h"
#include <iostream>
#include <sstream>
using namespace std;

AliTPCcalibAlign* AliTPCcalibAlign::fgInstance = 0;
ClassImp(AliTPCcalibAlign)




AliTPCcalibAlign* AliTPCcalibAlign::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgInstance == 0){
    fgInstance = new AliTPCcalibAlign();
  }
  return fgInstance;
}




AliTPCcalibAlign::AliTPCcalibAlign()
  :  AliTPCcalibBase(),
     fDphiHistArray(72*72),
     fDthetaHistArray(72*72),
     fDyHistArray(72*72),
     fDzHistArray(72*72),
     //
     fDyPhiHistArray(72*72),      // array of residual histograms  y     -kYPhi
     fDzThetaHistArray(72*72),    // array of residual histograms  z-z   -kZTheta
     fDphiZHistArray(72*72),      // array of residual histograms  phi   -kPhiz
     fDthetaZHistArray(72*72),    // array of residual histograms  theta -kThetaz
     fDyZHistArray(72*72),        // array of residual histograms  y     -kYz
     fDzZHistArray(72*72),        // array of residual histograms  z     -kZz     
     fFitterArray12(72*72),
     fFitterArray9(72*72),
     fFitterArray6(72*72),
     fMatrixArray12(72*72),
     fMatrixArray9(72*72),
     fMatrixArray6(72*72),
     fCombinedMatrixArray6(72),
     fCompTracklet(0),             // tracklet comparison
     fNoField(kFALSE)
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
   fDyPhiHistArray(72*72),      // array of residual histograms  y     -kYPhi
   fDzThetaHistArray(72*72),    // array of residual histograms  z-z   -kZTheta
   fDphiZHistArray(72*72),      // array of residual histograms  phi   -kPhiz
   fDthetaZHistArray(72*72),    // array of residual histograms  theta -kThetaz
   fDyZHistArray(72*72),        // array of residual histograms  y     -kYz
   fDzZHistArray(72*72),        // array of residual histograms  z     -kZz     //
   fFitterArray12(72*72),
   fFitterArray9(72*72),
   fFitterArray6(72*72),
   fMatrixArray12(72*72),
   fMatrixArray9(72*72),
   fMatrixArray6(72*72),
   fCombinedMatrixArray6(72),
   fCompTracklet(0),             // tracklet comparison
   fNoField(kFALSE)
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


AliTPCcalibAlign::AliTPCcalibAlign(const AliTPCcalibAlign &align)
  :AliTPCcalibBase(align),  
   fDphiHistArray(align.fDphiHistArray),
   fDthetaHistArray(align.fDthetaHistArray),
   fDyHistArray(align.fDyHistArray),
   fDzHistArray(align.fDzHistArray),
   fDyPhiHistArray(align.fDyPhiHistArray),      // array of residual histograms  y     -kYPhi
   fDzThetaHistArray(align.fDzThetaHistArray),    // array of residual histograms  z-z   -kZTheta
   fDphiZHistArray(align.fDphiZHistArray),      // array of residual histograms  phi   -kPhiz
   fDthetaZHistArray(align.fDthetaZHistArray),    // array of residual histograms  theta -kThetaz
   fDyZHistArray(align.fDyZHistArray),        // array of residual histograms  y     -kYz
   fDzZHistArray(align.fDzZHistArray),        // array of residual histograms  z     -kZz     
   //
   fFitterArray12(align.fFitterArray12),
   fFitterArray9(align.fFitterArray9),
   fFitterArray6(align.fFitterArray6),
   
   fMatrixArray12(align.fMatrixArray12),
   fMatrixArray9(align.fMatrixArray9),
   fMatrixArray6(align.fMatrixArray6),
   fCombinedMatrixArray6(align.fCombinedMatrixArray6),
   fCompTracklet(align.fCompTracklet),             // tracklet comparison
   fNoField(align.fNoField)   
{
  //
  // copy constructor - copy also the content
  //
  TH1 * his = 0;
  TObjArray * arr0=0;
  const TObjArray *arr1=0;
  for (Int_t index =0; index<72*72; index++){
    for (Int_t iarray=0;iarray<10; iarray++){
      if (iarray==kY){
	arr0 = &fDyHistArray;
	arr1 = &align.fDyHistArray;
      }
      if (iarray==kZ){
	arr0 = &fDzHistArray;
	arr1 = &align.fDzHistArray;
      }
      if (iarray==kPhi){
	arr0 = &fDphiHistArray;
	arr1 = &align.fDphiHistArray;
      }
      if (iarray==kTheta){
	arr0 = &fDthetaHistArray;
	arr1 = &align.fDthetaHistArray;
      }
      if (iarray==kYz){
	arr0 = &fDyZHistArray;
	arr1 = &align.fDyZHistArray;
      }
      if (iarray==kZz){
	arr0 = &fDzZHistArray;
	arr1 = &align.fDzZHistArray;
      }
      if (iarray==kPhiZ){
	arr0 = &fDphiZHistArray;
	arr1 = &align.fDphiZHistArray;
      }
      if (iarray==kThetaZ){
	arr0 = &fDthetaZHistArray;
	arr1 = &align.fDthetaZHistArray;
      }

      if (iarray==kYPhi){
	arr0 = &fDyPhiHistArray;
	arr1 = &align.fDyPhiHistArray;
      }
      if (iarray==kZTheta){
	arr0 = &fDzThetaHistArray;
	arr1 = &align.fDzThetaHistArray;
      }

      if (arr1->At(index)) {
	his = (TH1*)arr1->At(index)->Clone();
	his->SetDirectory(0);
	arr0->AddAt(his,index);
      }    
    }
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
  Double_t t1[10],t2[10];
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
  Int_t accept =   AcceptTracklet(tp1,tp2);  
  //
  //
  //
  if (fStreamLevel>1 && seed){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      static TVectorD vec1(5);
      static TVectorD vec2(5);
      vec1.SetElements(t1);
      vec2.SetElements(t2);
      AliExternalTrackParam *p1 = &((AliExternalTrackParam&)tp1);
      AliExternalTrackParam *p2 = &((AliExternalTrackParam&)tp2);
      (*cstream)<<"Tracklet"<<
	"accept="<<accept<<
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
	"isOK="<<accept<<           //  flag - used for alignment
	"tp1.="<<p1<<
	"tp2.="<<p2<<
	"v1.="<<&vec1<<
	"v2.="<<&vec2<<
	"s1="<<s1<<
	"s2="<<s2<<
	"\n";
    }
  }
  if (accept>0) return;
  t1[0]-=134.;
  t2[0]-=134.;
  t1[5]=0; t2[5]=0;
  t1[6]=TMath::Sqrt(tp1.GetSigmaY2()+tp2.GetSigmaY2()); t2[6]=t1[6];
  t1[7]=TMath::Sqrt(tp1.GetSigmaZ2()+tp2.GetSigmaZ2()); t2[7]=t1[7];
  t1[8]=TMath::Sqrt(tp1.GetSigmaSnp2()+tp2.GetSigmaSnp2()); t2[8]=t1[8];
  t1[9]=TMath::Sqrt(tp1.GetSigmaTgl2()+tp2.GetSigmaTgl2()); t2[9]=t1[9];

  if (GetDebugLevel()>50) printf("Process track\n");
  if (GetDebugLevel()>50) printf("Filling track\n");
  //
  // fill resolution histograms - previous cut included
  if (seed) ProcessDiff(tp1,tp2, seed,s1,s2);
  FillHisto(tp1,tp2,s1,s2);  
  ProcessAlign(t1,t2,s1,s2);
}

void AliTPCcalibAlign::ProcessAlign(Double_t * t1,
				    Double_t * t2,
				    Int_t s1,Int_t s2){
  //
  // Do intersector alignment
  //
  Process12(t1,t2,GetOrMakeFitter12(s1,s2));
  Process9(t1,t2,GetOrMakeFitter9(s1,s2));
  Process6(t1,t2,GetOrMakeFitter6(s1,s2));
  ++fPoints[GetIndex(s1,s2)];
}

void AliTPCcalibAlign::ProcessTree(TTree * chainTracklet, AliExternalComparison *comp){
  //
  // Process the debug streamer tree
  // Possible to modify selection criteria
  // Used with entry list
  //
  TTreeSRedirector * cstream = new TTreeSRedirector("aligndump.root");

  AliTPCcalibAlign *align = this;
  //
  TVectorD * vec1 = 0;
  TVectorD * vec2 = 0;
  AliExternalTrackParam * tp1 = 0;
  AliExternalTrackParam * tp2 = 0;  
  Int_t      s1 = 0;
  Int_t      s2 = 0;				
  Int_t npoints =0;
  {
    Int_t entries=chainTracklet->GetEntries();
    for (Int_t i=0; i< entries; i++){
      chainTracklet->GetBranch("tp1.")->SetAddress(&tp1);
      chainTracklet->GetBranch("tp2.")->SetAddress(&tp2);
      chainTracklet->GetBranch("v1.")->SetAddress(&vec1);
      chainTracklet->GetBranch("v2.")->SetAddress(&vec2);
      chainTracklet->GetBranch("s1")->SetAddress(&s1);
      chainTracklet->GetBranch("s2")->SetAddress(&s2);      
      chainTracklet->GetEntry(i);
      if (!vec1) continue;
      if (!vec2) continue;
      if (!tp1) continue;
      if (!tp2) continue;
      if (!vec1->GetMatrixArray()) continue;
      if (!vec2->GetMatrixArray()) continue;
      // make a local copy
      AliExternalTrackParam par1(*tp1);
      AliExternalTrackParam par2(*tp2);
      TVectorD svec1(*vec1);
      TVectorD svec2(*vec2);
      //
      if (s1==s2) continue;
      if (i%100==0) printf("%d\t%d\t%d\t%d\t\n",i, npoints,s1,s2);
      AliExternalTrackParam  cpar1(par1);
      AliExternalTrackParam  cpar2(par2);      
      Constrain1Pt(cpar1,par2,fNoField);
      Constrain1Pt(cpar2,par1,fNoField);
      Bool_t acceptComp = kFALSE;
      if (comp) acceptComp=comp->AcceptPair(&par1,&par2);
      if (comp) acceptComp&=comp->AcceptPair(&cpar1,&cpar2);
      //
      Int_t reject =   align->AcceptTracklet(par1,par2);
      Int_t rejectC =align->AcceptTracklet(cpar1,cpar2); 

      if (1||fStreamLevel>0){
	(*cstream)<<"Tracklet"<<
	  "s1="<<s1<<
	  "s2="<<s2<<
	  "reject="<<reject<<
	  "rejectC="<<rejectC<<
	  "acceptComp="<<acceptComp<<
	  "tp1.="<<&par1<<
	  "tp2.="<<&par2<<	
	  "ctp1.="<<&cpar1<<
	  "ctp2.="<<&cpar2<<
	  "v1.="<<&svec1<<
	  "v2.="<<&svec2<<
	  "\n";
      }
      //
      if (fNoField){
	//
	//
      }
      if (acceptComp) comp->Process(&cpar1,&cpar2);
      //
      if (reject>0 || rejectC>0) continue;
      npoints++;
      align->ProcessTracklets(cpar1,cpar2,0,s1,s2);
      align->ProcessTracklets(cpar2,cpar1,0,s2,s1); 
    }
  }
  delete cstream;
}


Int_t AliTPCcalibAlign::AcceptTracklet(const AliExternalTrackParam &p1,
				       const AliExternalTrackParam &p2){

  //
  // Accept pair of tracklets?
  //
  /*
  // resolution cuts
  TCut cutS0("sqrt(tp2.fC[0]+tp1.fC[0])<0.2");
  TCut cutS1("sqrt(tp2.fC[2]+tp1.fC[2])<0.2");
  TCut cutS2("sqrt(tp2.fC[5]+tp1.fC[5])<0.01");
  TCut cutS3("sqrt(tp2.fC[9]+tp1.fC[9])<0.01");
  TCut cutS4("sqrt(tp2.fC[14]+tp1.fC[14])<0.25");
  TCut cutS=cutS0+cutS1+cutS2+cutS3+cutS4;
  //
  // parameters matching cuts
  TCut cutP0("abs(tp1.fP[0]-tp2.fP[0])<0.6");
  TCut cutP1("abs(tp1.fP[1]-tp2.fP[1])<0.6");
  TCut cutP2("abs(tp1.fP[2]-tp2.fP[2])<0.03");
  TCut cutP3("abs(tp1.fP[3]-tp2.fP[3])<0.03");
  TCut cutP4("abs(tp1.fP[4]-tp2.fP[4])<0.5");
  TCut cutPP4("abs(tp1.fP[4]-tp2.fP[4])/sqrt(tp2.fC[14]+tp1.fC[14])<3");
  TCut cutP=cutP0+cutP1+cutP2+cutP3+cutP4+cutPP4;
  */  
  //
  // resolution cuts
  Int_t reject=0;
  const Double_t *cp1 = p1.GetCovariance();
  const Double_t *cp2 = p2.GetCovariance();
  if (TMath::Sqrt(cp1[0]+cp2[0])>0.2)  reject|=1;;
  if (TMath::Sqrt(cp1[2]+cp2[2])>0.2)  reject|=2;
  if (TMath::Sqrt(cp1[5]+cp2[5])>0.01) reject|=4;
  if (TMath::Sqrt(cp1[9]+cp2[9])>0.01) reject|=8;
  if (TMath::Sqrt(cp1[14]+cp2[14])>0.2) reject|=16;

  //parameters difference
  const Double_t *tp1 = p1.GetParameter();
  const Double_t *tp2 = p2.GetParameter();
  if (TMath::Abs(tp1[0]-tp2[0])>0.6) reject|=32;
  if (TMath::Abs(tp1[1]-tp2[1])>0.6) reject|=64;
  if (TMath::Abs(tp1[2]-tp2[2])>0.03) reject|=128;
  if (TMath::Abs(tp1[3]-tp2[3])>0.03) reject|=526;
  if (TMath::Abs(tp1[4]-tp2[4])>0.4) reject|=1024;
  if (TMath::Abs(tp1[4]-tp2[4])/TMath::Sqrt(cp1[14]+cp2[14])>4) reject|=2048;
  
  //
  if (TMath::Abs(tp2[1])>235) reject|=2*4096;
  
  if (fNoField){
    
  }

  return reject;
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

  //
  Double_t sy    = t1[6];
  Double_t sz    = t1[7];
  Double_t sdydx = t1[8];
  Double_t sdzdx = t1[9];

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
  //
  Double_t sy    = t1[6];
  Double_t sz    = t1[7];
  Double_t sdydx = t1[8];
  Double_t sdzdx = t1[9];
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

  //
  Double_t sy    = t1[6];
  Double_t sz    = t1[7];
  Double_t sdydx = t1[8];
  Double_t sdzdx = t1[9];

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
  //p[0]   +=  dydx1*dydx2; // a01         FIXME- 0912 MI
  //p[]   += -dzdx1*dydx2; // a02
  value  = -dydx1+dydx2;  // -a11 + a00
  fitter->AddPoint(p,value,sdydx);

  // dzdx2 = (a20 + a21*dydx1 + a22*dzdx1)/( a00 + a01*dydx1 + a02*dzdx1)
  // (a20 + a21*dydx1 + a22*dzdx1) - (a00 + a01*dydx1 + a02*dzdx1)*dzdx2 = 0
  for (Int_t i=0; i<12;i++) p[i]=0.;
  p[1]   += 1;            // a20
  //  p[2]   += dydx1;        // a21   FIXME- 0912 MI
  //p[]  += dzdx1;        // a22
  //p[]  += -dzdx2;       // a00
  //p[0]   +=  dydx1*dzdx2; // a01     FIXME- 0912 MI
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
  TMatrixD mat(4,4);
  for (Int_t s1=0;s1<72;++s1)
    for (Int_t s2=0;s2<72;++s2){
      if (GetTransformation12(s1,s2,mat)){
	fMatrixArray12.AddAt(mat.Clone(), GetIndex(s1,s2));
      }
      if (GetTransformation9(s1,s2,mat)){
	fMatrixArray9.AddAt(mat.Clone(), GetIndex(s1,s2));
      }
      if (GetTransformation6(s1,s2,mat)){
	fMatrixArray6.AddAt(mat.Clone(), GetIndex(s1,s2));
      }
    }
  //this->Write("align");
  
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
  fitter->StoreData(kTRUE);
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
    //
    GetHisto(kPhiZ,s1,s2,kTRUE)->Fill(tp1.GetZ(),TMath::ASin(tp1.GetSnp())-TMath::ASin(tp2.GetSnp()));    
    GetHisto(kThetaZ,s1,s2,kTRUE)->Fill(tp1.GetZ(),TMath::ATan(tp1.GetTgl())-TMath::ATan(tp2.GetTgl()));
    GetHisto(kYz,s1,s2,kTRUE)->Fill(tp1.GetZ(),tp1.GetY()-tp2.GetY());
    GetHisto(kZz,s1,s2,kTRUE)->Fill(tp1.GetZ(),tp1.GetZ()-tp2.GetZ());
    //
    GetHisto(kYPhi,s1,s2,kTRUE)->Fill(tp1.GetSnp(),tp1.GetY()-tp2.GetY());
    GetHisto(kZTheta,s1,s2,kTRUE)->Fill(tp1.GetTgl(),tp1.GetZ()-tp2.GetZ());


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
  case kYPhi:
    histoArray = &fDyPhiHistArray; break;
  case kZTheta:
    histoArray = &fDzThetaHistArray; break;
  case kYz:
    histoArray = &fDyZHistArray; break;
  case kZz:
    histoArray = &fDzZHistArray; break;
  case kPhiZ:
    histoArray = &fDphiZHistArray; break;
  case kThetaZ:
    histoArray = &fDthetaZHistArray; break;
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
    //
    //
  case kYPhi:
    name<<"hist_yphi_"<<s1<<"_"<<s2;
    title<<"Y Missalignment for sectors Phi"<<s1<<" and "<<s2;
    histo =new TH2F(name.str().c_str(),title.str().c_str(),20,-1,1,128,-0.3,0.3); // +/- 3 mm
    break;
  case kZTheta:
    name<<"hist_ztheta_"<<s1<<"_"<<s2;
    title<<"Z Missalignment for sectors Theta"<<s1<<" and "<<s2;
    histo = new TH2F(name.str().c_str(),title.str().c_str(),128,20,-1,1,-0.3,0.3); // +/- 3 mm
    break;
    //
    //
    //
  case kYz:
    name<<"hist_yz_"<<s1<<"_"<<s2;
    title<<"Y Missalignment for sectors Z"<<s1<<" and "<<s2;
    histo =new TH2F(name.str().c_str(),title.str().c_str(),20,-250,250,128,-0.3,0.3); // +/- 3 mm
    break;
  case kZz:
    name<<"hist_zz_"<<s1<<"_"<<s2;
    title<<"Z Missalignment for sectors Z"<<s1<<" and "<<s2;
    histo = new TH2F(name.str().c_str(),title.str().c_str(),20,-250,250,128,-0.3,0.3); // +/- 3 mm
    break;
  case kPhiZ:
    name<<"hist_phiz_"<<s1<<"_"<<s2;
    title<<"Phi Missalignment for sectors Z"<<s1<<" and "<<s2;
    histo =new TH2F(name.str().c_str(),title.str().c_str(),20,-250,250,128,-0.01,0.01); // +/- 10 mrad
    break;
  case kThetaZ:
    name<<"hist_thetaz_"<<s1<<"_"<<s2;
    title<<"Theta Missalignment for sectors Z"<<s1<<" and "<<s2;
    histo =new TH2F(name.str().c_str(),title.str().c_str(),20,-250,250,128,-0.01,0.01); // +/- 10 mrad
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
  //TObjArray *fitArray=0;
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
      TMatrixD m6;
      TMatrixD m6FX;
      TMatrixD m9;
      TMatrixD m12;
      Double_t dy=0, dz=0, dphi=0,dtheta=0;
      Double_t sy=0, sz=0, sphi=0,stheta=0;
      Double_t ny=0, nz=0, nphi=0,ntheta=0;
      Double_t chi2v12=0, chi2v9=0, chi2v6=0;
      Int_t npoints=0;
      TLinearFitter * fitter = 0;      
      if (fPoints[GetIndex(s1,s2)]>kMinPoints){
	//
	//
	//
	fitter = GetFitter12(s1,s2);
	npoints = fitter->GetNpoints();
	chi2v12 = TMath::Sqrt(fitter->GetChisquare()/npoints);
	//
	fitter = GetFitter9(s1,s2);
	npoints = fitter->GetNpoints();
	chi2v9 = TMath::Sqrt(fitter->GetChisquare()/npoints);
	//
	fitter = GetFitter6(s1,s2);
	npoints = fitter->GetNpoints();
	chi2v6 = TMath::Sqrt(fitter->GetChisquare()/npoints);

	//
	GetTransformation6(s1,s2,m6);
	GetTransformation9(s1,s2,m9);
	GetTransformation12(s1,s2,m12);
	//
	fitter = GetFitter6(s1,s2);
	fitter->FixParameter(3,0);
	fitter->Eval();
	GetTransformation6(s1,s2,m6FX);
	//
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

      }

      // x2    =  a00*x1 + a01*y1 + a02*z1 + a03
      // y2    =  a10*x1 + a11*y1 + a12*z1 + a13
      // z2    =  a20*x1 + a21*y1 + a22*z1 + a23
      // dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
      // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
      //
      //                     a00  a01 a02  a03     p[0]   p[1]  p[2]  p[9]
      //                     a10  a11 a12  a13 ==> p[3]   p[4]  p[5]  p[10]
      //                     a20  a21 a22  a23     p[6]   p[7]  p[8]  p[11] 
      
      //
      // 
      // dy:-(134*m6.fElements[4]+m6.fElements[7])
      // 
      // dphi:-(m6.fElements[4])
      //
      // dz:134*m6.fElements[8]+m6.fElements[11]
      //
      // dtheta:m6.fElements[8]
      //
      cstream<<"Align"<<
	"s1="<<s1<<     // reference sector
	"s2="<<s2<<     // sector to align
	"m6FX.="<<&m6FX<<   // tranformation matrix
	"m6.="<<&m6<<   // tranformation matrix
	"m9.="<<&m9<<   // 
	"m12.="<<&m12<<
	"chi2v12="<<chi2v12<<
	"chi2v9="<<chi2v9<<
	"chi2v6="<<chi2v6<<
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
  //
  TString str1(GetName());
  while((obj = iter->Next()) != 0)
    {
      AliTPCcalibAlign* entry = dynamic_cast<AliTPCcalibAlign*>(obj);
      if (entry == 0) continue; 
      if (str1.CompareTo(entry->GetName())!=0) continue;
      Add(entry);
      count++;
    } 
  return count;
}


void AliTPCcalibAlign::Add(AliTPCcalibAlign * align){
  //
  // Add entry - used for merging of compoents
  //
  for (Int_t i=0; i<72;i++){
    for (Int_t j=0; j<72;j++){
      if (align->fPoints[GetIndex(i,j)]<10) continue;
      fPoints[GetIndex(i,j)]+=align->fPoints[GetIndex(i,j)];
      //
      //
      //
      for (Int_t itype=0; itype<10; itype++){
	TH1 * his0=0, *his1=0;
	his0 = GetHisto((HistoType)itype,i,j);
	his1 = align->GetHisto((HistoType)itype,i,j);
	if (his1){
	  if (his0) his0->Add(his1);
	  else {
	    his0 = GetHisto(kY,i,j,kTRUE);
	    his0->Add(his1);
	  }
	}   	
      }      
    }
  }
  TLinearFitter *f0=0;
  TLinearFitter *f1=0;
  for (Int_t i=0; i<72;i++){
    for (Int_t j=0; j<72;j++){     
      if (align->fPoints[GetIndex(i,j)]<20) continue;
      //
      //
      // fitter12
      f0 =  GetFitter12(i,j);
      f1 =  align->GetFitter12(i,j);
      if (f1 &&f1->Eval()){
	if (f0&&f0->GetNpoints()>10) f0->Add(f1);
	else {
	  fFitterArray12.AddAt(f1->Clone(),GetIndex(i,j));
	}
      }      
      //
      // fitter9
      f0 =  GetFitter9(i,j);
      f1 =  align->GetFitter9(i,j);
      if (f1&&f1->Eval()){
	if (f0&&f0->GetNpoints()>10) f0->Add(f1);
	else { 
	  fFitterArray9.AddAt(f1->Clone(),GetIndex(i,j));
	}
      }      
      f0 =  GetFitter6(i,j);
      f1 =  align->GetFitter6(i,j);
      if (f1 &&f1->Eval()){
	if (f0&&f0->GetNpoints()>10) f0->Add(f1);
	else {
	  fFitterArray6.AddAt(f1->Clone(),GetIndex(i,j));
	}
      }   
    }
  }
}

Double_t AliTPCcalibAlign::Correct(Int_t type, Int_t value, Int_t s1, Int_t s2, Double_t x1, Double_t y1, Double_t z1, Double_t dydx1,Double_t dzdx1){
  //
  // GetTransformed value
  //
  //
  // x2    =  a00*x1 + a01*y1 + a02*z1 + a03
  // y2    =  a10*x1 + a11*y1 + a12*z1 + a13
  // z2    =  a20*x1 + a21*y1 + a22*z1 + a23
  // dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  
  
  const TMatrixD * mat = GetTransformation(s1,s2,type);
  if (!mat) {
    if (value==0) return x1;
    if (value==1) return y1;
    if (value==2) return z1;
    if (value==3) return dydx1;
    if (value==4) return dzdx1;
    //
    if (value==5) return dydx1;
    if (value==6) return dzdx1;
    return 0;
  }
  Double_t valT=0;

  if (value==0){
    valT = (*mat)(0,0)*x1+(*mat)(0,1)*y1+(*mat)(0,2)*z1+(*mat)(0,3);
  }

  if (value==1){
    valT = (*mat)(1,0)*x1+(*mat)(1,1)*y1+(*mat)(1,2)*z1+(*mat)(1,3);
  }
  if (value==2){
    valT = (*mat)(2,0)*x1+(*mat)(2,1)*y1+(*mat)(2,2)*z1+(*mat)(2,3);
  }
  if (value==3){
    //    dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
    valT =  (*mat)(1,0)    +(*mat)(1,1)*dydx1 +(*mat)(1,2)*dzdx1;
    valT/= ((*mat)(0,0)    +(*mat)(0,1)*dydx1 +(*mat)(0,2)*dzdx1);
  }

  if (value==4){
    // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)    
    valT =  (*mat)(2,0)    +(*mat)(2,1)*dydx1 +(*mat)(2,2)*dzdx1;
    valT/= ((*mat)(0,0)    +(*mat)(0,1)*dydx1 +(*mat)(0,2)*dzdx1);
  }
  //
  if (value==5){
    // onlys shift in angle
    //    dydx2 =  (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
    valT =  (*mat)(1,0)    +(*mat)(1,1)*dydx1;
  }

  if (value==6){
    // only shift in angle
    // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)    
    valT =  (*mat)(2,0)    +(*mat)(2,1)*dydx1;
  }
  //
  return valT;
}


void  AliTPCcalibAlign::Constrain1Pt(AliExternalTrackParam &track1, const AliExternalTrackParam &track2, Bool_t noField){
  //
  // Update track parameters t1
  //
  TMatrixD vecXk(5,1);    // X vector
  TMatrixD covXk(5,5);    // X covariance 
  TMatrixD matHk(1,5);    // vector to mesurement
  TMatrixD measR(1,1);    // measurement error 
  //TMatrixD matQk(5,5);    // prediction noise vector
  TMatrixD vecZk(1,1);    // measurement
  //
  TMatrixD vecYk(1,1);    // Innovation or measurement residual
  TMatrixD matHkT(5,1);
  TMatrixD matSk(1,1);    // Innovation (or residual) covariance
  TMatrixD matKk(5,1);    // Optimal Kalman gain
  TMatrixD mat1(5,5);     // update covariance matrix
  TMatrixD covXk2(5,5);   // 
  TMatrixD covOut(5,5);
  //
  Double_t *param1=(Double_t*) track1.GetParameter();
  Double_t *covar1=(Double_t*) track1.GetCovariance();

  //
  // copy data to the matrix
  for (Int_t ipar=0; ipar<5; ipar++){
    vecXk(ipar,0) = param1[ipar];
    for (Int_t jpar=0; jpar<5; jpar++){
      covXk(ipar,jpar) = covar1[track1.GetIndex(ipar, jpar)];
    }
  }
  //
  //
  //
  vecZk(0,0) = track2.GetParameter()[4];   // 1/pt measurement from track 2
  measR(0,0) = track2.GetCovariance()[14];  // 1/pt measurement error
  if (noField) {
    measR(0,0)*=0.000000001;
    vecZk(0,0)=0.;
  }
  //
  matHk(0,0)=0; matHk(0,1)= 0; matHk(0,2)= 0;  
  matHk(0,3)= 0;    matHk(0,4)= 1;           // vector to measurement
  //
  //
  //
  vecYk = vecZk-matHk*vecXk;                 // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;      // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;              //  Optimal Kalman gain
  vecXk += matKk*vecYk;                      //  updated vector 
  mat1(0,0)=1; mat1(1,1)=1; mat1(2,2)=1; mat1(3,3)=1; mat1(4,4)=1;
  covXk2 = (mat1-(matKk*matHk));
  covOut =  covXk2*covXk; 
  //
  //
  //
  // copy from matrix to parameters
  if (0) {
    covOut.Print();
    vecXk.Print();
    covXk.Print();
    track1.Print();
    track2.Print();
  }

  for (Int_t ipar=0; ipar<5; ipar++){
    param1[ipar]= vecXk(ipar,0) ;
    for (Int_t jpar=0; jpar<5; jpar++){
      covar1[track1.GetIndex(ipar, jpar)]=covOut(ipar,jpar);
    }
  }

}

void AliTPCcalibAlign::GlobalAlign6(Int_t minPoints, Float_t sysError, Int_t niter){
  //
  //  Global Align -combine the partial alignment of pair of sectors
  //  minPoints - minimal number of points - don't use sector alignment wit smaller number
  //  sysError  - error added to the alignemnt error
  //
  AliTPCcalibAlign * align = this;
  TMatrixD * arrayAlign[72];
  TMatrixD * arrayAlignDiff[72];
  //
  for (Int_t i=0;i<72; i++) {
    TMatrixD * mat = new TMatrixD(4,4);
    mat->UnitMatrix();
    arrayAlign[i]=mat;
    arrayAlignDiff[i]=(TMatrixD*)(mat->Clone());
  }

  TTreeSRedirector *cstream = new TTreeSRedirector("galign6.root");
  for (Int_t iter=0; iter<niter;iter++){
    printf("Iter=\t%d\n",iter);
    for (Int_t is0=0;is0<72; is0++) {
      //
      //TMatrixD  *mati0 = arrayAlign[is0];
      TMatrixD matDiff(4,4);
      Double_t sumw=0;      
      for (Int_t is1=0;is1<72; is1++) {
	Bool_t invers=kFALSE;
	Int_t npoints=0;
	TMatrixD covar;
	TVectorD errors;
	const TMatrixD *mat = align->GetTransformation(is0,is1,0); 
	if (mat){
	  npoints = align->GetFitter6(is0,is1)->GetNpoints();
	  if (npoints>minPoints){
	    align->GetFitter6(is0,is1)->GetCovarianceMatrix(covar);
	    align->GetFitter6(is0,is1)->GetErrors(errors);
	  }
	}
	else{
	  invers=kTRUE;
	  mat = align->GetTransformation(is1,is0,0); 
	  if (mat) {
	    npoints = align->GetFitter6(is1,is0)->GetNpoints();
	    if (npoints>minPoints){
	      align->GetFitter6(is1,is0)->GetCovarianceMatrix(covar);
	      align->GetFitter6(is1,is0)->GetErrors(errors);
	    }
	  }
	}
	if (!mat) continue;
	if (npoints<minPoints) continue;
	//
	Double_t weight=1;
 	if (is1/36>is0/36) weight*=2./3.; //IROC-OROC
 	if (is1/36<is0/36) weight*=1./3.; //OROC-IROC
 	if (is1/36==is0/36) weight*=1/3.; //OROC-OROC
 	if (is1%36!=is0%36) weight*=1/2.; //Not up-down
	weight/=(errors[4]*errors[4]+sysError*sysError);
	//
	//
	TMatrixD matT = *mat;	
	if (invers)  matT.Invert();
	TMatrixD diffMat= (*(arrayAlign[is1]))*matT;
	diffMat-=(*arrayAlign[is0]);
	matDiff+=weight*diffMat;
	sumw+=weight;

	(*cstream)<<"LAlign"<<
	  "iter="<<iter<<
	  "s0="<<is0<<
	  "s1="<<is1<<
	  "npoints="<<npoints<<
	  "m60.="<<arrayAlign[is0]<<
	  "m61.="<<arrayAlign[is1]<<
	  "m01.="<<&matT<<
	  "diff.="<<&diffMat<<
	  "cov.="<<&covar<<
	  "err.="<<&errors<<
	  "\n";
      }
      if (sumw>0){
	matDiff*=1/sumw;
	matDiff(0,0)=0;
	matDiff(1,1)=0;
	matDiff(1,1)=0;
	matDiff(1,1)=0;	
	(*arrayAlignDiff[is0]) = matDiff;       
      }
    }
    for (Int_t is0=0;is0<72; is0++) {
      if (is0<36) (*arrayAlign[is0]) += 0.4*(*arrayAlignDiff[is0]);
      if (is0>=36) (*arrayAlign[is0]) += 0.2*(*arrayAlignDiff[is0]);
       //
      (*cstream)<<"GAlign"<<
	"iter="<<iter<<
	"s0="<<is0<<
	"m6.="<<arrayAlign[is0]<<
	"\n";
    }    
  }
  delete cstream;
  for (Int_t isec=0;isec<72;isec++){
    fCombinedMatrixArray6.AddAt(arrayAlign[isec],isec);
    delete arrayAlignDiff[isec];
  }
}


/*
  

gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
AliXRDPROOFtoolkit tool;
TChain * chainTr = tool.MakeChain("align.txt","Track",0,10200);
chainTr->Lookup();
TChain * chainTracklet = tool.MakeChain("align.txt","Tracklet",0,20);
chainTracklet->Lookup();


TCut cutS0("sqrt(tp2.fC[0]+tp2.fC[0])<0.6");
TCut cutS1("sqrt(tp2.fC[2]+tp2.fC[2])<0.6");
TCut cutS2("sqrt(tp2.fC[5]+tp2.fC[5])<0.04");
TCut cutS3("sqrt(tp2.fC[9]+tp2.fC[9])<0.02");
TCut cutS4("sqrt(tp2.fC[14]+tp2.fC[14])<0.5");
// resolution cuts
TCut cutS=cutS0+cutS1+cutS2+cutS3+cutS4;

TCut cutP0("abs(tp1.fP[0]-tp2.fP[0])<0.6");
TCut cutP1("abs(tp1.fP[1]-tp2.fP[1])<0.6");
TCut cutP2("abs(tp1.fP[2]-tp2.fP[2])<0.02");
TCut cutP3("abs(tp1.fP[3]-tp2.fP[3])<0.01");
TCut cutE("abs(tp2.fP[1])<235");
TCut cutP=cutP0+cutP1+cutP2+cutP3+cutE;


//
//
TCut cutA =  cutP+cutS;
chainTr->Draw(">>listEL",cutA,"entryList");
TEntryList *elist = (TEntryList*)gDirectory->Get("listEL");
chainTr->SetEntryList(elist);






TCut cutS("s1%36==s2%36");

TCut cutN("c1>32&&c2>60");
TCut cutC0("sqrt(tp2.fC[0])<1");

TCut cutX("abs(tp2.fX-133.6)<2");

TCut cutA = cutP+cutN;


TCut cutY("abs(vcY.fElements-vtY.fElements)<0.3&&vcY.fElements!=0")
TCut cutZ("abs(vcZ.fElements-vtZ.fElements)<0.3&&vcZ.fElements!=0")






TCut cutA =  cutP+cutS;
chainTracklet->Draw(">>listEL",cutA,"entryList");
TEntryList *elist = (TEntryList*)gDirectory->Get("listEL");
chainTracklet->SetEntryList(elist);

//
TVectorD * vec1 = 0;
TVectorD * vec2 = 0;
AliExternalTrackParam * tp1 = 0;
AliExternalTrackParam * tp2 = 0;

Int_t      s1 = 0;
Int_t      s2 = 0;
chainTracklet->GetBranch("v1.")->SetAddress(&vec1);
chainTracklet->GetBranch("v2.")->SetAddress(&vec2);
chainTracklet->GetBranch("s1")->SetAddress(&s1);
chainTracklet->GetBranch("s2")->SetAddress(&s2);


AliTPCcalibAlign align;
{
for (Int_t i=0; i< elist->GetN(); i++){
//for (Int_t i=0; i<100000; i++){
chainTracklet->GetBranch("tp1.")->SetAddress(&tp1);
chainTracklet->GetBranch("tp2.")->SetAddress(&tp2);
chainTracklet->GetBranch("v1.")->SetAddress(&vec1);
chainTracklet->GetBranch("v2.")->SetAddress(&vec2);
chainTracklet->GetBranch("s1")->SetAddress(&s1);
chainTracklet->GetBranch("s2")->SetAddress(&s2);

chainTracklet->GetEntry(i);
if (i%100==0) printf("%d\t%d\t%d\t%d\t\n",i,tentry, s1,s2);
//vec1.Print();
TLinearFitter * fitter = align.GetOrMakeFitter6(s1,s2);
if (fitter) align.Process6(vec1->GetMatrixArray(),vec2->GetMatrixArray(),fitter);
}
}

*/
