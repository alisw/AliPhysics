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
#include "AliTPCROC.h"
#include "AliTPCPointCorrection.h"
#include "AliTrackPointArray.h"

#include "AliExternalTrackParam.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"

#include "AliTPCTracklet.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THnSparse.h"
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
#include "AliLog.h"
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"


#include "TTreeStream.h"
#include "Riostream.h"
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
     fNoField(kFALSE),
     fXIO(0),
     fXmiddle(0),
     fXquadrant(0),
     fArraySectorIntParam(36), // array of sector alignment parameters
     fArraySectorIntCovar(36), // array of sector alignment covariances 
     //
     // Kalman filter for global alignment
     //
     fSectorParamA(0),     // Kalman parameter   for A side
     fSectorCovarA(0),     // Kalman covariance  for A side 
     fSectorParamC(0),     // Kalman parameter   for A side
     fSectorCovarC(0),     // Kalman covariance  for A side 
     fUseInnerOuter(kTRUE)// flag- use Inner Outer sector for left righ alignment
{
  //
  // Constructor
  //
  for (Int_t i=0;i<72*72;++i) {
    fPoints[i]=0;
  }
  AliTPCROC * roc = AliTPCROC::Instance();
  fXquadrant = roc->GetPadRowRadii(36,53);
  fXmiddle   = ( roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  fXIO       = ( roc->GetPadRowRadii(0,roc->GetNRows(0)-1)+roc->GetPadRowRadii(36,0))*0.5;
  fClusterDelta[0]=0;   // cluster residuals
  fClusterDelta[1]=0;   // cluster residuals
  fClusterDelta[2]=0;   // cluster residuals - vertex constrained
  fClusterDelta[3]=0;   // cluster residuals
  fClusterDelta[4]=0;   // cluster residuals - ITS constrained
  fClusterDelta[5]=0;   // cluster residuals
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
   fNoField(kFALSE),
   fXIO(0),
   fXmiddle(0),
   fXquadrant(0),
   fArraySectorIntParam(36), // array of sector alignment parameters
   fArraySectorIntCovar(36), // array of sector alignment covariances 
   //
   // Kalman filter for global alignment
   //
   fSectorParamA(0),     // Kalman parameter   for A side
   fSectorCovarA(0),     // Kalman covariance  for A side 
   fSectorParamC(0),     // Kalman parameter   for A side
   fSectorCovarC(0),     // Kalman covariance  for A side      
   fUseInnerOuter(kTRUE)// flag- use Inner Outer sector for left righ alignment

{
  //
  // Constructor
  //
  SetName(name);
  SetTitle(title);
  for (Int_t i=0;i<72*72;++i) {
    fPoints[i]=0;
  }
  AliTPCROC * roc = AliTPCROC::Instance();
  fXquadrant = roc->GetPadRowRadii(36,53);
  fXmiddle   = ( roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  fXIO       = ( roc->GetPadRowRadii(0,roc->GetNRows(0)-1)+roc->GetPadRowRadii(36,0))*0.5;
  fClusterDelta[0]=0;   // cluster residuals
  fClusterDelta[1]=0;   // cluster residuals
  fClusterDelta[2]=0;   // cluster residuals - vertex constrained
  fClusterDelta[3]=0;   // cluster residuals
  fClusterDelta[4]=0;   // cluster residuals - ITS constrained
  fClusterDelta[5]=0;   // cluster residuals



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
   fNoField(align.fNoField),
   fXIO(align.fXIO),   
   fXmiddle(align.fXmiddle),   
   fXquadrant(align.fXquadrant),   
   fArraySectorIntParam(align.fArraySectorIntParam), // array of sector alignment parameters
   fArraySectorIntCovar(align.fArraySectorIntCovar), // array of sector alignment covariances 
   fSectorParamA(0),     // Kalman parameter   for A side
   fSectorCovarA(0),     // Kalman covariance  for A side 
   fSectorParamC(0),     // Kalman parameter   for A side
   fSectorCovarC(0),      // Kalman covariance  for A side 
  fUseInnerOuter(kTRUE)// flag- use Inner Outer sector for left righ alignment
  
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
  //
  //
  //
  if (align.fSectorParamA){
    fSectorParamA = (TMatrixD*)align.fSectorParamA->Clone();
    fSectorParamA = (TMatrixD*)align.fSectorCovarA->Clone();
    fSectorParamC = (TMatrixD*)align.fSectorParamA->Clone();
    fSectorParamC = (TMatrixD*)align.fSectorCovarA->Clone();
  }
  fClusterDelta[0]=0;   // cluster residuals
  fClusterDelta[1]=0;   // cluster residuals
  fClusterDelta[2]=0;   // cluster residuals - vertex constrained
  fClusterDelta[3]=0;   // cluster residuals
  fClusterDelta[4]=0;   // cluster residuals - ITS constrained
  fClusterDelta[5]=0;   // cluster residuals

}


AliTPCcalibAlign::~AliTPCcalibAlign() {
  //
  // destructor
  //
  fDphiHistArray.SetOwner(kTRUE);    // array of residual histograms  phi      -kPhi
  fDthetaHistArray.SetOwner(kTRUE);  // array of residual histograms  theta    -kTheta
  fDyHistArray.SetOwner(kTRUE);      // array of residual histograms  y        -kY
  fDzHistArray.SetOwner(kTRUE);      // array of residual histograms  z        -kZ
  //
  fDyPhiHistArray.SetOwner(kTRUE);      // array of residual histograms  y     -kYPhi
  fDzThetaHistArray.SetOwner(kTRUE);    // array of residual histograms  z-z   -kZTheta
  //
  fDphiZHistArray.SetOwner(kTRUE);      // array of residual histograms  phi   -kPhiz
  fDthetaZHistArray.SetOwner(kTRUE);    // array of residual histograms  theta -kThetaz
  fDyZHistArray.SetOwner(kTRUE);        // array of residual histograms  y     -kYz
  fDzZHistArray.SetOwner(kTRUE);        // array of residual histograms  z     -kZz

  fDphiHistArray.Delete();    // array of residual histograms  phi      -kPhi
  fDthetaHistArray.Delete();  // array of residual histograms  theta    -kTheta
  fDyHistArray.Delete();      // array of residual histograms  y        -kY
  fDzHistArray.Delete();      // array of residual histograms  z        -kZ
  //
  fDyPhiHistArray.Delete();      // array of residual histograms  y     -kYPhi
  fDzThetaHistArray.Delete();    // array of residual histograms  z-z   -kZTheta
  //
  fDphiZHistArray.Delete();      // array of residual histograms  phi   -kPhiz
  fDthetaZHistArray.Delete();    // array of residual histograms  theta -kThetaz
  fDyZHistArray.Delete();        // array of residual histograms  y     -kYz
  fDzZHistArray.Delete();        // array of residual histograms  z     -kZz

  fFitterArray12.SetOwner(kTRUE);    // array of fitters
  fFitterArray9.SetOwner(kTRUE);     // array of fitters
  fFitterArray6.SetOwner(kTRUE);     // array of fitters
  //
  fMatrixArray12.SetOwner(kTRUE);    // array of transnformtation matrix
  fMatrixArray9.SetOwner(kTRUE);     // array of transnformtation matrix
  fMatrixArray6.SetOwner(kTRUE);     // array of transnformtation matrix 
  //
  fFitterArray12.Delete();    // array of fitters
  fFitterArray9.Delete();     // array of fitters
  fFitterArray6.Delete();     // array of fitters
  //
  fMatrixArray12.Delete();    // array of transnformtation matrix
  fMatrixArray9.Delete();     // array of transnformtation matrix
  fMatrixArray6.Delete();     // array of transnformtation matrix 

  if (fCompTracklet) delete fCompTracklet;

  fArraySectorIntParam.SetOwner(kTRUE); // array of sector alignment parameters
  fArraySectorIntCovar.SetOwner(kTRUE); // array of sector alignment covariances 
  fArraySectorIntParam.Delete(); // array of sector alignment parameters
  fArraySectorIntCovar.Delete(); // array of sector alignment covariances 
  for (Int_t i=0; i<6; i++){
    delete fClusterDelta[i];   // cluster residuals
  }
}

void AliTPCcalibAlign::Process(AliESDEvent *event) {
  //
  // Process pairs of cosmic tracks
  //
  if (!fClusterDelta[0])  MakeResidualHistos();
  //
  fCurrentEvent=event;
  ExportTrackPoints(event);  // export track points for external calibration 
  const Int_t kMaxTracks =6;
  const Int_t kminCl = 40;
  AliESDfriend *eESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!eESDfriend) return;
  Int_t ntracks=event->GetNumberOfTracks(); 
  Float_t dca0[2];
  Float_t dca1[2];
  //
  //
  // process seeds
  //
  for (Int_t i0=0;i0<ntracks;++i0) {
    AliESDtrack *track0 = event->GetTrack(i0);
    AliESDfriendTrack *friendTrack = 0;
    TObject *calibObject=0;
    AliTPCseed *seed0 = 0;
    //
    friendTrack = (AliESDfriendTrack *)eESDfriend->GetTrack(i0);;
    if (!friendTrack) continue;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed0=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }
    if (!seed0) continue;
    fCurrentTrack=track0;
    fCurrentSeed=seed0;
    fCurrentEvent=event;
    ProcessSeed(seed0);
  }
  //
  // process cosmic pairs
  //
  if (ntracks>kMaxTracks) return;  
  //
  //select pairs - for alignment  
  for (Int_t i0=0;i0<ntracks;++i0) {
    AliESDtrack *track0 = event->GetTrack(i0);
    //    if (track0->GetTPCNcls()<kminCl) continue;
    track0->GetImpactParameters(dca0[0],dca0[1]);
    //    if (TMath::Abs(dca0[0])>30) continue;
    //
    for (Int_t i1=0;i1<ntracks;++i1) {
      if (i0==i1) continue;
      AliESDtrack *track1 = event->GetTrack(i1);
      //      if (track1->GetTPCNcls()<kminCl) continue;
      track1->GetImpactParameters(dca1[0],dca1[1]);
      // fast cuts on dca and theta
      //      if (TMath::Abs(dca1[0]+dca0[0])>15) continue;
      //      if (TMath::Abs(dca1[1]-dca0[1])>15) continue;
      if (TMath::Abs(track0->GetParameter()[3]+track1->GetParameter()[3])>0.1) continue;
      //
      AliESDfriendTrack *friendTrack = 0;
      TObject *calibObject=0;
      AliTPCseed *seed0 = 0,*seed1=0;
      //
      friendTrack = (AliESDfriendTrack *)eESDfriend->GetTrack(i0);;
      if (!friendTrack) continue;
      for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
	if ((seed0=dynamic_cast<AliTPCseed*>(calibObject))) break;
      }
      friendTrack = (AliESDfriendTrack *)eESDfriend->GetTrack(i1);;
      if (!friendTrack) continue;
      for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
	if ((seed1=dynamic_cast<AliTPCseed*>(calibObject))) break;
      }
      if (!seed0) continue;
      //
      //
      if (!seed1) continue;
      Int_t nclsectors0[72], nclsectors1[72];
      for (Int_t isec=0;isec<72;isec++){
	nclsectors0[isec]=0;
	nclsectors1[isec]=0;
      }
      for (Int_t i=0;i<160;i++){
	AliTPCclusterMI *c0=seed0->GetClusterPointer(i);
	AliTPCclusterMI *c1=seed1->GetClusterPointer(i);
	if (c0)  nclsectors0[c0->GetDetector()]+=1;
	if (c1)  nclsectors1[c1->GetDetector()]+=1;
      }

      for (Int_t isec0=0; isec0<72;isec0++){
	if (nclsectors0[isec0]<kminCl) continue;
	for (Int_t isec1=0; isec1<72;isec1++){
	  if (nclsectors1[isec1]<kminCl) continue;
	  Int_t s0 = isec0;
	  Int_t s1 = isec1;
	  Double_t parLine0[10];
	  Double_t parLine1[10];
	  TMatrixD par0(4,1),cov0(4,4),par1(4,1),cov1(4,4);
	  Bool_t useInnerOuter = kFALSE;
	  if (s1%36!=s0%36) useInnerOuter = fUseInnerOuter;  // for left - right alignment both sectors refit can be used if specified
	  Int_t nl0 = RefitLinear(seed0,s0, parLine0, s0,par0,cov0,fXIO,useInnerOuter);
	  Int_t nl1 = RefitLinear(seed1,s1, parLine1, s0,par1,cov1,fXIO,useInnerOuter);
	  parLine0[0]=0;  // reference frame in IO boundary
	  parLine1[0]=0;
	  //      if (nl0<kminCl || nl1<kminCl) continue;
	  //
	  //
	  Bool_t isOK=kTRUE;
	  if (TMath::Min(nl0,nl1)<kminCl) isOK=kFALSE;
	  // apply selection criteria
	  //
	  Float_t dp0,dp1,dp3;
	  Float_t pp0,pp1,pp3;
	  dp0=par0(0,0)-par1(0,0); 
	  dp1=par0(1,0)-par1(1,0); 
	  dp3=par0(3,0)-par1(3,0); 
	  pp0=dp0/TMath::Sqrt(cov0(0,0)+cov1(0,0)+0.1*0.1);
	  pp1=dp1/TMath::Sqrt(cov0(1,1)+cov1(1,1)+0.0015*0.0015);
	  pp3=dp3/TMath::Sqrt(cov0(3,3)+cov1(3,3)+0.0015*0.0015);
	  //
	  if (TMath::Abs(dp0)>1.0)  isOK=kFALSE;
	  if (TMath::Abs(dp1)>0.02) isOK=kFALSE;
	  if (TMath::Abs(dp3)>0.02) isOK=kFALSE;
	  if (TMath::Abs(pp0)>6)  isOK=kFALSE;
	  if (TMath::Abs(pp1)>6) isOK=kFALSE;
	  if (TMath::Abs(pp3)>6) isOK=kFALSE;	  
	  //
	  if (isOK){
	    FillHisto(parLine0,parLine1,s0,s1);  
	    ProcessAlign(parLine0,parLine1,s0,s1);
	    UpdateKalman(s0,s1,par0, cov0, par1, cov1);
	  }
	  if (fStreamLevel>0){
	    TTreeSRedirector *cstream = GetDebugStreamer();
	    if (cstream){
	      (*cstream)<<"cosmic"<<
		"isOK="<<isOK<<
		"s0="<<s0<<
		"s1="<<s1<<
		"nl0="<<nl0<<
		"nl1="<<nl1<<
		"p0.="<<&par0<<
		"p1.="<<&par1<<
		"c0.="<<&cov0<<
		"c1.="<<&cov1<<
		"\n";
	    }
	  }
	}
      }
    }
  }
}

void  AliTPCcalibAlign::ExportTrackPoints(AliESDEvent *event){
  //
  // Export track points for alignment - calibration
  // export space points for pairs of tracks if possible
  //
  AliESDfriend *eESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!eESDfriend) return;
  Int_t ntracks=event->GetNumberOfTracks();
  Int_t kMaxTracks=4;   // maximal number of tracks for cosmic pairs
  Int_t kMinVertexTracks=5;   // maximal number of tracks for vertex mesurement

  //cuts
  const Int_t kminCl     = 60;
  const Int_t kminClSum  = 120;
  //const Double_t kDistY  = 5;
  // const Double_t kDistZ  = 40;
  const Double_t kDistTh = 0.05;
  const Double_t kDistThS = 0.002;
  const Double_t kDist1Pt = 0.1;
  const Double_t kMaxD0 =3;  // max distance to the primary vertex
  const Double_t kMaxD1 =5;  // max distance to the primary vertex
  const AliESDVertex *tpcVertex = 0;
  // get the primary vertex TPC
  if (ntracks>kMinVertexTracks) {
    tpcVertex = event->GetPrimaryVertexSPD();
    if (tpcVertex->GetNContributors()<kMinVertexTracks) tpcVertex=0;
  }
  //
  Float_t dca0[2];
  //  Float_t dca1[2];
  Int_t index0=0,index1=0;
  //
  for (Int_t i0=0;i0<ntracks;++i0) {
    AliESDtrack *track0 = event->GetTrack(i0);
    if (!track0) continue;
    if ((track0->GetStatus() & AliESDtrack::kTPCrefit)==0) continue;
    if (track0->GetOuterParam()==0) continue;
    if (track0->GetInnerParam()==0) continue;
    if (TMath::Abs(track0->GetInnerParam()->GetSigned1Pt()-track0->GetOuterParam()->GetSigned1Pt())>kDist1Pt) continue;
    if (TMath::Abs(track0->GetInnerParam()->GetSigned1Pt())>kDist1Pt) continue;
    if (TMath::Abs(track0->GetInnerParam()->GetTgl()-track0->GetOuterParam()->GetTgl())>kDistThS) continue;
    AliESDtrack *track1P = 0;
    if (track0->GetTPCNcls()<kminCl) continue;
    track0->GetImpactParameters(dca0[0],dca0[1]);
    index0=i0;
    index1=-1;
    //
    if (ntracks<kMaxTracks) for (Int_t i1=i0+1;i1<ntracks;++i1) {
      if (i0==i1) continue;
      AliESDtrack *track1 = event->GetTrack(i1);
      if (!track1) continue;
      if ((track1->GetStatus() & AliESDtrack::kTPCrefit)==0) continue;
      if (track1->GetOuterParam()==0) continue;
      if (track1->GetInnerParam()==0) continue;
      if (track1->GetTPCNcls()<kminCl) continue;
      if (TMath::Abs(track1->GetInnerParam()->GetSigned1Pt()-track1->GetOuterParam()->GetSigned1Pt())>kDist1Pt) continue;
      if (TMath::Abs(track1->GetInnerParam()->GetTgl()-track1->GetOuterParam()->GetTgl())>kDistThS) continue;
      if (TMath::Abs(track1->GetInnerParam()->GetSigned1Pt())>kDist1Pt) continue;
      //track1->GetImpactParameters(dca1[0],dca1[1]);
      //if (TMath::Abs(dca1[0]-dca0[0])>kDistY) continue;
      //if (TMath::Abs(dca1[1]-dca0[1])>kDistZ) continue;
      if (TMath::Abs(track0->GetTgl()+track1->GetTgl())>kDistTh) continue;
      if (TMath::Abs(track0->GetSigned1Pt()+track1->GetSigned1Pt())>kDist1Pt) continue;
      track1P = track1;
      index1=i1;
    }
    AliESDfriendTrack *friendTrack = 0;
    TObject *calibObject=0;
    AliTPCseed *seed0 = 0,*seed1=0;
    //
    friendTrack = (AliESDfriendTrack *)eESDfriend->GetTrack(index0);;
    if (!friendTrack) continue;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed0=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }
    if (index1>0){
      friendTrack = (AliESDfriendTrack *)eESDfriend->GetTrack(index1);;
      if (!friendTrack) continue;
      for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
	if ((seed1=dynamic_cast<AliTPCseed*>(calibObject))) break;
      }
    }
    //
    Int_t npoints=0, ncont=0;
    if (seed0) {npoints+=seed0->GetNumberOfClusters(); ncont++;}
    if (seed1) {npoints+=seed1->GetNumberOfClusters(); ncont++;}
    if (npoints<kminClSum) continue;    
    Int_t cpoint=0;
    AliTrackPointArray array(npoints);    
    if (tpcVertex){
      Double_t dxyz[3]={0,0,0};
      Double_t dc[6]={0,0,0};
      tpcVertex->GetXYZ(dxyz);
      tpcVertex->GetCovarianceMatrix(dc);
      Float_t xyz[3]={dxyz[0],dxyz[1],dxyz[2]};
      Float_t cov[6]={dc[0]+1,dc[1],dc[2]+1,dc[3],dc[4],dc[5]+100.}; 
      AliTrackPoint point(xyz,cov,73);  // add point to not existing volume
      Float_t dtpc[2],dcov[3];
      track0->GetImpactParametersTPC(dtpc,dcov);
      if (TMath::Abs(dtpc[0])>kMaxD0) continue;
      if (TMath::Abs(dtpc[1])>kMaxD1) continue;
    }

    if (seed0) for (Int_t icl = 0; icl<160; icl++){
      AliTPCclusterMI *cluster=seed0->GetClusterPointer(icl);
      if (!cluster) continue;
      Float_t xyz[3];
      Float_t cov[6];
      cluster->GetGlobalXYZ(xyz);
      cluster->GetGlobalCov(cov);
      AliTrackPoint point(xyz,cov,cluster->GetDetector());
      array.AddPoint(npoints, &point); 
      if (cpoint>=npoints) continue;  //shoul not happen
      array.AddPoint(cpoint, &point);
      cpoint++;
    }
    if (seed1) for (Int_t icl = 0; icl<160; icl++){
      AliTPCclusterMI *cluster=seed1->GetClusterPointer(icl);
      if (!cluster) continue;
      Float_t xyz[3];
      Float_t cov[6];
      cluster->GetGlobalXYZ(xyz);
      cluster->GetGlobalCov(cov);
      AliTrackPoint point(xyz,cov,cluster->GetDetector());
      array.AddPoint(npoints, &point);
      if (cpoint>=npoints) continue;  //shoul not happen
      array.AddPoint(cpoint, &point);
      cpoint++;
    }
    //
    //
    //
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      Bool_t isVertex=(tpcVertex)? kTRUE:kFALSE;
      Double_t tof0=track0->GetTOFsignal();
      Double_t tof1=(track1P) ?  track1P->GetTOFsignal(): 0;
      static AliExternalTrackParam dummy;
      AliExternalTrackParam *p0In  = &dummy;
      AliExternalTrackParam *p1In  = &dummy;
      AliExternalTrackParam *p0Out = &dummy;
      AliExternalTrackParam *p1Out = &dummy;
      AliESDVertex vdummy;
      AliESDVertex *pvertex= (tpcVertex)? (AliESDVertex *)tpcVertex: &vdummy;
      if (track0) {
	p0In= new AliExternalTrackParam(*track0);
	p0Out=new AliExternalTrackParam(*(track0->GetOuterParam()));
      }
      if (track1P) {
	p1In= new AliExternalTrackParam(*track1P);
	p1Out=new AliExternalTrackParam(*(track1P->GetOuterParam()));
      }

      (*cstream)<<"trackPoints"<<
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"triggerClass="<<&fTriggerClass<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
	"pvertex.="<<pvertex<<      // vertex
	//
	"isVertex="<<isVertex<<     // flag is with prim vertex
	"tof0="<<tof0<<             // tof signal 0
	"tof1="<<tof1<<             // tof signal 1
	"seed0.="<<seed0<<          // track info
	"ntracks="<<ntracks<<       // number of tracks
	"ncont="<<ncont<<           // number of contributors
	"p0In.="<<p0In<<              // track parameters0 
	"p1In.="<<p1In<<              // track parameters1
	"p0Out.="<<p0Out<<          // track parameters0
	"p1Out.="<<p1Out<<          // track parameters0
	"p.="<<&array<<
	"\n";
    }
  }  
}




void AliTPCcalibAlign::ProcessSeed(AliTPCseed *seed) {
  //
  // 
  //
  // make a kalman tracklets out of seed
  //
  TObjArray tracklets=
    AliTPCTracklet::CreateTracklets(seed,AliTPCTracklet::kKalman,
				    kFALSE,20,4);
  tracklets.SetOwner();
  Int_t ntracklets = tracklets.GetEntries();
  if (ntracklets<2) return;
  //
  //
  for (Int_t i1=0;i1<ntracklets;i1++)
    for (Int_t i2=0;i2<ntracklets;i2++){
      if (i1==i2) continue;
      AliTPCTracklet *t1=static_cast<AliTPCTracklet*>(tracklets[i1]);
      AliTPCTracklet *t2=static_cast<AliTPCTracklet*>(tracklets[i2]);
      AliExternalTrackParam *common1=0,*common2=0;
      if (AliTPCTracklet::PropagateToMeanX(*t1,*t2,common1,common2)){
	ProcessTracklets(*common1,*common2,seed, t1->GetSector(),t2->GetSector());	
	UpdateAlignSector(seed,t1->GetSector());
      }
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


void AliTPCcalibAlign::UpdatePointCorrection(AliTPCPointCorrection * correction){
  //
  // Update point correction with alignment coefficients
  //
  for (Int_t isec=0;isec<36;isec++){
    TMatrixD * matCorr = (TMatrixD*)(correction->fArraySectorIntParam.At(isec));
    TMatrixD * matAlign = (TMatrixD*)(fArraySectorIntParam.At(isec));
    TMatrixD * matAlignCovar = (TMatrixD*)(fArraySectorIntCovar.At(isec));
    if (!matAlign) continue;
    if (!matCorr) {
      correction->fArraySectorIntParam.AddAt(matAlign->Clone(),isec); 
      correction->fArraySectorIntCovar.AddAt(matAlignCovar->Clone(),isec); 
      continue;
    }
    (*matCorr)+=(*matAlign);
    correction->fArraySectorIntCovar.AddAt(matAlignCovar->Clone(),isec); 
  }
  //

}


void AliTPCcalibAlign::ProcessTracklets(const AliExternalTrackParam &tp1,
					const AliExternalTrackParam &tp2,
					const AliTPCseed * seed,
					Int_t s1,Int_t s2) {
  //
  // Process function to fill fitters
  //
  Double_t t1[10],t2[10];
  Double_t &x1=t1[0], &y1=t1[1], &z1=t1[3], &dydx1=t1[2], &dzdx1=t1[4];
  Double_t &x2=t2[0], &y2=t2[1], &z2=t2[3], &dydx2=t2[2], &dzdx2=t2[4];
  x1   =tp1.GetX();
  y1   =tp1.GetY();
  z1   =tp1.GetZ();
  Double_t snp1=tp1.GetSnp();
  dydx1=snp1/TMath::Sqrt((1.-snp1)*(1.+snp1));
  Double_t tgl1=tp1.GetTgl();
  // dz/dx = 1/(cos(theta)*cos(phi))
  dzdx1=tgl1/TMath::Sqrt((1.-snp1)*(1.+snp1));
  x2   =tp2.GetX();
  y2   =tp2.GetY();
  z2   =tp2.GetZ();
  Double_t snp2=tp2.GetSnp();
  dydx2=snp2/TMath::Sqrt((1.-snp2)*(1.+snp2));
  Double_t tgl2=tp2.GetTgl();
  dzdx2=tgl2/TMath::Sqrt((1.-snp2)*(1.+snp2));
  
  //
  // Kalman parameters
  //
  t1[0]-=fXIO;
  t2[0]-=fXIO;
  // errors
  t1[5]=0; t2[5]=0;
  t1[6]=TMath::Sqrt(tp1.GetSigmaY2());
  t1[7]=TMath::Sqrt(tp1.GetSigmaSnp2());
  t1[8]=TMath::Sqrt(tp1.GetSigmaZ2()); 
  t1[9]=TMath::Sqrt(tp1.GetSigmaTgl2());
  
  t2[6]=TMath::Sqrt(tp2.GetSigmaY2());
  t2[7]=TMath::Sqrt(tp2.GetSigmaSnp2()); 
  t2[8]=TMath::Sqrt(tp2.GetSigmaZ2());
  t2[9]=TMath::Sqrt(tp2.GetSigmaTgl2());
  //
  // linear parameters
  //
  Double_t parLine1[10];
  Double_t parLine2[10];
  TMatrixD par1(4,1),cov1(4,4),par2(4,1),cov2(4,4);
  Bool_t useInnerOuter = kFALSE;
  if (s1%36!=s2%36) useInnerOuter = fUseInnerOuter;  // for left - right alignment bot sectors refit can be used if specified
  Int_t nl1 = RefitLinear(seed,s1, parLine1, s1,par1,cov1,tp1.GetX(), useInnerOuter);
  Int_t nl2 = RefitLinear(seed,s2, parLine2, s1,par2,cov2,tp1.GetX(), useInnerOuter);
  parLine1[0]=tp1.GetX()-fXIO;   // parameters in  IROC-OROC boundary
  parLine2[0]=tp1.GetX()-fXIO;   // parameters in  IROC-OROC boundary
  //
  //
  //
  Int_t accept       =   AcceptTracklet(tp1,tp2);  
  Int_t acceptLinear =   AcceptTracklet(parLine1,parLine2);
 
  if (fStreamLevel>1 && seed){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      static TVectorD vec1(5);
      static TVectorD vec2(5);
      static TVectorD vecL1(9);
      static TVectorD vecL2(9);
      vec1.SetElements(t1);
      vec2.SetElements(t2);
      vecL1.SetElements(parLine1);
      vecL2.SetElements(parLine2);
      AliExternalTrackParam *p1 = &((AliExternalTrackParam&)tp1);
      AliExternalTrackParam *p2 = &((AliExternalTrackParam&)tp2);
      (*cstream)<<"Tracklet"<<
	"accept="<<accept<<
	"acceptLinear="<<acceptLinear<<  // accept linear tracklets
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"triggerClass="<<&fTriggerClass<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
	"isOK="<<accept<<           //  flag - used for alignment
	"tp1.="<<p1<<
	"tp2.="<<p2<<
	"v1.="<<&vec1<<
	"v2.="<<&vec2<<
	"s1="<<s1<<
	"s2="<<s2<<
	"nl1="<<nl1<<       // linear fit - n points
	"nl2="<<nl2<<       // linear fit - n points
	"vl1.="<<&vecL1<<   // linear fits
	"vl2.="<<&vecL2<<   // linear fits
	"\n";
    }
  }
  if (TMath::Abs(fMagF)<0.005){
    //
    // use Linear fit
    //
    if (nl1>10 && nl2>10 &&(acceptLinear==0)){
      if (seed) ProcessDiff(tp1,tp2, seed,s1,s2);
      if (TMath::Abs(parLine1[2])<0.8 &&TMath::Abs(parLine1[2])<0.8 ){ //angular cut
	FillHisto(parLine1,parLine2,s1,s2);  
	ProcessAlign(parLine1,parLine2,s1,s2);
	//UpdateKalman(s1,s2,par1, cov1, par2, cov2); - OBSOLETE to be removed - 50 % of time here
      }
    }
  }
  if (accept>0) return;
  //
  // fill resolution histograms - previous cut included
  if (TMath::Abs(fMagF)>0.005){
    //
    // use Kalman if mag field
    //
    if (seed) ProcessDiff(tp1,tp2, seed,s1,s2);
    FillHisto(t1,t2,s1,s2);  
    ProcessAlign(t1,t2,s1,s2);
  }
}

void AliTPCcalibAlign::ProcessAlign(Double_t * t1,
				    Double_t * t2,
				    Int_t s1,Int_t s2){
  //
  // Do intersector alignment
  //
  //Process12(t1,t2,GetOrMakeFitter12(s1,s2));
  //Process9(t1,t2,GetOrMakeFitter9(s1,s2));
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
				       const AliExternalTrackParam &p2) const
{

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


Int_t  AliTPCcalibAlign::AcceptTracklet(const Double_t *t1, const Double_t *t2) const
{
  //
  // accept tracklet  - 
  //  dist cut + 6 sigma cut 
  //
  Double_t dy     = t2[1]-t1[1];
  Double_t dphi   = t2[2]-t1[2];
  Double_t dz     = t2[3]-t1[3];
  Double_t dtheta = t2[4]-t1[4];
  //
  Double_t sy       = TMath::Sqrt(t1[6]*t1[6]+t2[6]*t2[6]+0.05*0.05);
  Double_t sdydx    = TMath::Sqrt(t1[7]*t1[7]+t2[7]*t2[7]+0.001*0.001);
  Double_t sz       = TMath::Sqrt(t1[8]*t1[8]+t2[8]*t2[8]+0.05*0.05);
  Double_t sdzdx    = TMath::Sqrt(t1[9]*t1[9]+t2[9]*t2[9]+0.001*0.001);
  //
  Int_t reject=0;
  if (TMath::Abs(dy)>1.)         reject|=2;
  if (TMath::Abs(dphi)>0.1)      reject|=4;
  if (TMath::Abs(dz)>1.)         reject|=8;
  if (TMath::Abs(dtheta)>0.1)    reject|=16;
  //
  if (TMath::Abs(dy/sy)>6)         reject|=32;
  if (TMath::Abs(dphi/sdydx)>6)    reject|=64;
  if (TMath::Abs(dz/sz)>6)         reject|=128;
  if (TMath::Abs(dtheta/sdzdx)>6)  reject|=256;
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
  if (fStreamLevel>5 && count1>10 && count2>10){
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
	"triggerClass="<<&fTriggerClass<<      //  trigger
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
				 TLinearFitter *fitter) const
{
  // x2    =  a00*x1 + a01*y1 + a02*z1 + a03
  // y2    =  a10*x1 + a11*y1 + a12*z1 + a13
  // z2    =  a20*x1 + a21*y1 + a22*z1 + a23
  // dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  //
  //                     a00  a01 a02  a03     p[0]   p[1]  p[2]  p[9]
  //                     a10  a11 a12  a13 ==> p[3]   p[4]  p[5]  p[10]
  //                     a20  a21 a22  a23     p[6]   p[7]  p[8]  p[11] 



  const Double_t &x1=t1[0], &y1=t1[1], &z1=t1[3], &dydx1=t1[2], &dzdx1=t1[4];
  const Double_t /*&x2=t2[0],*/ &y2=t2[1], &z2=t2[3], &dydx2=t2[2], &dzdx2=t2[4];

  //
  Double_t sy       = TMath::Sqrt(t1[6]*t1[6]+t2[6]*t2[6]);
  Double_t sdydx    = TMath::Sqrt(t1[7]*t1[7]+t2[7]*t2[7]);
  Double_t sz       = TMath::Sqrt(t1[8]*t1[8]+t2[8]*t2[8]);
  Double_t sdzdx    = TMath::Sqrt(t1[9]*t1[9]+t2[9]*t2[9]);

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

void AliTPCcalibAlign::Process9(const Double_t * const t1,
				const Double_t * const t2,
				TLinearFitter *fitter) const
{
  // x2    =  a00*x1 + a01*y1 + a02*z1 + a03
  // y2    =  a10*x1 + a11*y1 + a12*z1 + a13
  // z2    =  a20*x1 + a21*y1 + a22*z1 + a23
  // dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  //
  //                     a00  a01  a02 a03     1      p[0]  p[1]   p[6]
  //                     a10  a11  a12 a13 ==> p[2]   1     p[3]   p[7]
  //                     a20  a21  a21 a23     p[4]   p[5]  1      p[8] 


  const Double_t &x1=t1[0], &y1=t1[1], &z1=t1[3], &dydx1=t1[2], &dzdx1=t1[4];
  const Double_t /*&x2=t2[0],*/ &y2=t2[1], &z2=t2[3], &dydx2=t2[2], &dzdx2=t2[4];
  //
  Double_t sy       = TMath::Sqrt(t1[6]*t1[6]+t2[6]*t2[6]);
  Double_t sdydx    = TMath::Sqrt(t1[7]*t1[7]+t2[7]*t2[7]);
  Double_t sz       = TMath::Sqrt(t1[8]*t1[8]+t2[8]*t2[8]);
  Double_t sdzdx    = TMath::Sqrt(t1[9]*t1[9]+t2[9]*t2[9]);

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

void AliTPCcalibAlign::Process6(const Double_t *const t1,
				const Double_t *const t2,
				TLinearFitter *fitter) const
{
  // x2    =  1  *x1 +-a01*y1 + 0      +a03
  // y2    =  a01*x1 + 1  *y1 + 0      +a13
  // z2    =  a20*x1 + a21*y1 + 1  *z1 +a23
  // dydx2 = (a10    + a11*dydx1 + a12*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  // dzdx2 = (a20    + a21*dydx1 + a22*dzdx1)/(a00    + a01*dydx1 + a02*dzdx1)
  //
  //                     a00  a01  a02 a03     1     -p[0]  0     p[3]
  //                     a10  a11  a12 a13 ==> p[0]   1     0     p[4]
  //                     a20  a21  a21 a23     p[1]   p[2]  1     p[5] 

  const Double_t &x1=t1[0], &y1=t1[1], &z1=t1[3], &dydx1=t1[2], &dzdx1=t1[4];
  const Double_t /*&x2=t2[0],*/ &y2=t2[1], &z2=t2[3], &dydx2=t2[2], &dzdx2=t2[4];

  //
  Double_t sy       = TMath::Sqrt(t1[6]*t1[6]+t2[6]*t2[6]);
  Double_t sdydx    = TMath::Sqrt(t1[7]*t1[7]+t2[7]*t2[7]);
  Double_t sz       = TMath::Sqrt(t1[8]*t1[8]+t2[8]*t2[8]);
  Double_t sdzdx    = TMath::Sqrt(t1[9]*t1[9]+t2[9]*t2[9]);

 
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




void AliTPCcalibAlign::EvalFitters(Int_t minPoints) {
  //
  // Analyze function 
  // 
  // Perform the fitting using linear fitters
  //
  TLinearFitter *f;
  TFile fff("alignDebug.root","recreate");
  for (Int_t s1=0;s1<72;++s1)
    for (Int_t s2=0;s2<72;++s2){
      if ((f=GetFitter12(s1,s2))&&fPoints[GetIndex(s1,s2)]>minPoints) {
	//	cerr<<s1<<","<<s2<<": "<<fPoints[GetIndex(s1,s2)]<<endl;
	if (f->Eval()!=0) {
	  cerr<<"Evaluation failed for "<<s1<<","<<s2<<endl;
	  f->Write(Form("f12_%d_%d",s1,s2));
	}else{
	  f->Write(Form("f12_%d_%d",s1,s2));
	}
      }
      if ((f=GetFitter9(s1,s2))&&fPoints[GetIndex(s1,s2)]>minPoints) {
	//	cerr<<s1<<","<<s2<<": "<<fPoints[GetIndex(s1,s2)]<<endl;
	if (f->Eval()!=0) {
	  cerr<<"Evaluation failed for "<<s1<<","<<s2<<endl;
	}else{
	  f->Write(Form("f9_%d_%d",s1,s2));
	}
      }
      if ((f=GetFitter6(s1,s2))&&fPoints[GetIndex(s1,s2)]>minPoints) {
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
  //  if (GetDebugLevel()>0) cerr<<"Creating fitter12 "<<s1<<","<<s2<<"  :  "<<counter12<<endl;
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
  //  if (GetDebugLevel()>0) cerr<<"Creating fitter12 "<<s1<<","<<s2<<"  :  "<<counter9<<endl;
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
  //  if (GetDebugLevel()>0) cerr<<"Creating fitter6 "<<s1<<","<<s2<<"  :  "<<counter6<<endl;
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

void AliTPCcalibAlign::MakeResidualHistos(){
  //
  // Make cluster residual histograms
  //
  Double_t xminTrack[9], xmaxTrack[9];
  Int_t    binsTrack[9];
  TString  axisName[9],axisTitle[9];
  //
  // 0 - delta   of interest
  // 1 - global  phi in sector number  as float
  // 2 - local   x
  // 3 - local   ky
  // 4 - local   kz
  // 
  axisName[0]="delta";   axisTitle[0]="#Delta (cm)"; 
  binsTrack[0]=60;       xminTrack[0]=-0.6;        xmaxTrack[0]=0.6; 
  //
  axisName[1]="sector";   axisTitle[1]="Sector Number"; 
  binsTrack[1]=180;       xminTrack[1]=0;        xmaxTrack[1]=18; 
  //
  axisName[2]="localX";   axisTitle[2]="x (cm)"; 
  binsTrack[2]=53;       xminTrack[2]=85.;        xmaxTrack[2]=245.; 
  //
  axisName[3]="kY";      axisTitle[3]="dy/dx"; 
  binsTrack[3]=1;       xminTrack[3]=-0.16;        xmaxTrack[3]=0.16; 
  //
  axisName[4]="kZ";      axisTitle[4]="dz/dx"; 
  binsTrack[4]=22;       xminTrack[4]=-1.1;        xmaxTrack[4]=1.1; 
  //
  fClusterDelta[0] = new THnSparseF("#Delta_{Y} (cm)","#Delta_{Y} (cm)", 5, binsTrack,xminTrack, xmaxTrack);
  fClusterDelta[1] = new THnSparseF("#Delta_{Z} (cm)","#Delta_{Z} (cm)", 5, binsTrack,xminTrack, xmaxTrack);
  fClusterDelta[2] = new THnSparseF("#Delta_{Y} (cm) const","#Delta_{Y} (cm) const ", 5, binsTrack,xminTrack, xmaxTrack);
  fClusterDelta[3] = new THnSparseF("#Delta_{Z} (cm) const","#Delta_{Z} (cm) const", 5, binsTrack,xminTrack, xmaxTrack);
  fClusterDelta[4] = new THnSparseF("#Delta_{Y} (cm) ITS","#Delta_{Y} (cm) ITS", 5, binsTrack,xminTrack, xmaxTrack);
  fClusterDelta[5] = new THnSparseF("#Delta_{Z} (cm) ITS","#Delta_{Z} (cm) ITS", 5, binsTrack,xminTrack, xmaxTrack);
  //
  //
  //
  for (Int_t ivar=0;ivar<6;ivar++){
    for (Int_t ivar2=0;ivar2<5;ivar2++){
      fClusterDelta[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      fClusterDelta[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    }
  }

}

void AliTPCcalibAlign::FillHisto(const Double_t *t1,
				 const Double_t *t2,
				 Int_t s1,Int_t s2) {
  //
  // Fill residual histograms
  // Track2-Track1
  // Innner-Outer
  // Left right - x-y
  // A-C side
  if (1)  {  
    Double_t dy     = t2[1]-t1[1];
    Double_t dphi   = t2[2]-t1[2];
    Double_t dz     = t2[3]-t1[3];
    Double_t dtheta = t2[4]-t1[4];
    Double_t zmean = (t2[3]+t1[3])*0.5;
    //
    GetHisto(kPhi,s1,s2,kTRUE)->Fill(dphi);    
    GetHisto(kTheta,s1,s2,kTRUE)->Fill(dtheta);
    GetHisto(kY,s1,s2,kTRUE)->Fill(dy);
    GetHisto(kZ,s1,s2,kTRUE)->Fill(dz);
    //
    GetHisto(kPhiZ,s1,s2,kTRUE)->Fill(zmean,dphi);    
    GetHisto(kThetaZ,s1,s2,kTRUE)->Fill(zmean,dtheta);
    GetHisto(kYz,s1,s2,kTRUE)->Fill(zmean,dy);
    GetHisto(kZz,s1,s2,kTRUE)->Fill(zmean,dz);
    //
    GetHisto(kYPhi,s1,s2,kTRUE)->Fill(t2[2],dy);
    GetHisto(kZTheta,s1,s2,kTRUE)->Fill(t2[4],dz);
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
    histo =new TH1D(name.str().c_str(),title.str().c_str(),100,-0.5,0.5); // +/- 5 mm
    break;
  case kZ:
    name<<"hist_z_"<<s1<<"_"<<s2;
    title<<"Z Missalignment for sectors "<<s1<<" and "<<s2;
    histo = new TH1D(name.str().c_str(),title.str().c_str(),100,-0.3,0.3); // +/- 3 mm
    break;
  case kPhi:
    name<<"hist_phi_"<<s1<<"_"<<s2;
    title<<"Phi Missalignment for sectors "<<s1<<" and "<<s2;
    histo =new TH1D(name.str().c_str(),title.str().c_str(),100,-0.01,0.01); // +/- 10 mrad
    break;
  case kTheta:
    name<<"hist_theta_"<<s1<<"_"<<s2;
    title<<"Theta Missalignment for sectors "<<s1<<" and "<<s2;
    histo =new TH1D(name.str().c_str(),title.str().c_str(),100,-0.01,0.01); // +/- 10 mrad
    break;
    //
    //
  case kYPhi:
    name<<"hist_yphi_"<<s1<<"_"<<s2;
    title<<"Y Missalignment for sectors Phi"<<s1<<" and "<<s2;
    histo =new TH2F(name.str().c_str(),title.str().c_str(),20,-1,1,100,-0.5,0.5); // +/- 5 mm
    break;
  case kZTheta:
    name<<"hist_ztheta_"<<s1<<"_"<<s2;
    title<<"Z Missalignment for sectors Theta"<<s1<<" and "<<s2;
    histo = new TH2F(name.str().c_str(),title.str().c_str(),20,-1,1,100,-0.3,0.3); // +/- 3 mm
    break;
    //
    //
    //
  case kYz:
    name<<"hist_yz_"<<s1<<"_"<<s2;
    title<<"Y Missalignment for sectors Z"<<s1<<" and "<<s2;
    histo =new TH2F(name.str().c_str(),title.str().c_str(),20,-250,250,100,-0.5,0.5); // +/- 5 mm
    break;
  case kZz:
    name<<"hist_zz_"<<s1<<"_"<<s2;
    title<<"Z Missalignment for sectors Z"<<s1<<" and "<<s2;
    histo = new TH2F(name.str().c_str(),title.str().c_str(),20,-250,250,100,-0.3,0.3); // +/- 3 mm
    break;
  case kPhiZ:
    name<<"hist_phiz_"<<s1<<"_"<<s2;
    title<<"Phi Missalignment for sectors Z"<<s1<<" and "<<s2;
    histo =new TH2F(name.str().c_str(),title.str().c_str(),20,-250,250,100,-0.01,0.01); // +/- 10 mrad
    break;
  case kThetaZ:
    name<<"hist_thetaz_"<<s1<<"_"<<s2;
    title<<"Theta Missalignment for sectors Z"<<s1<<" and "<<s2;
    histo =new TH2F(name.str().c_str(),title.str().c_str(),20,-250,250,100,-0.01,0.01); // +/- 10 mrad
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

void  AliTPCcalibAlign::MakeTree(const char *fname, Int_t minPoints){
  //
  // make tree with alignment cosntant  -
  // For  QA visualization
  //
  /*
    TFile fcalib("CalibObjects.root");
    TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
    AliTPCcalibAlign * align = ( AliTPCcalibAlign *)array->FindObject("alignTPC");
    align->EvalFitters();
    align->MakeTree("alignTree.root");
    TFile falignTree("alignTree.root");
    TTree * treeAlign = (TTree*)falignTree.Get("Align");
   */
  TTreeSRedirector cstream(fname);
  for (Int_t s1=0;s1<72;++s1)
    for (Int_t s2=0;s2<72;++s2){
      TMatrixD m6;
      TMatrixD m6FX;
      TMatrixD m9;
      TMatrixD m12;
      TVectorD param6Diff;  // align parameters diff 
      TVectorD param6s1(6);    // align parameters sector1 
      TVectorD param6s2(6);    // align parameters sector2 

      //
      //
      TMatrixD * kpar = fSectorParamA;
      TMatrixD * kcov = fSectorCovarA;
      if (s1%36>=18){
	kpar = fSectorParamC;
	kcov = fSectorCovarC;
      }
      for (Int_t ipar=0;ipar<6;ipar++){
	Int_t isec1 = s1%18;
	Int_t isec2 = s2%18;
	if (s1>35) isec1+=18;
	if (s2>35) isec2+=18;	
	param6s1(ipar)=(*kpar)(6*isec1+ipar,0);
	param6s2(ipar)=(*kpar)(6*isec2+ipar,0);
      }


      Double_t dy=0, dz=0, dphi=0,dtheta=0;
      Double_t sy=0, sz=0, sphi=0,stheta=0;
      Double_t ny=0, nz=0, nphi=0,ntheta=0;
      Double_t chi2v12=0, chi2v9=0, chi2v6=0;
      //      Int_t npoints=0;
      // TLinearFitter * fitter = 0;      
      if (fPoints[GetIndex(s1,s2)]>minPoints){
	//
	//
	//
// 	fitter = GetFitter12(s1,s2);
// 	npoints = fitter->GetNpoints();
// 	chi2v12 = TMath::Sqrt(fitter->GetChisquare()/npoints);
	
// 	//
// 	fitter = GetFitter9(s1,s2);
// 	npoints = fitter->GetNpoints();
// 	chi2v9 = TMath::Sqrt(fitter->GetChisquare()/npoints);
// 	//
// 	fitter = GetFitter6(s1,s2);
// 	npoints = fitter->GetNpoints();
// 	chi2v6 = TMath::Sqrt(fitter->GetChisquare()/npoints);
// 	fitter->GetParameters(param6Diff);
// 	//
// 	GetTransformation6(s1,s2,m6);
// 	GetTransformation9(s1,s2,m9);
// 	GetTransformation12(s1,s2,m12);
// 	//
// 	fitter = GetFitter6(s1,s2);
// 	//fitter->FixParameter(3,0);
// 	//fitter->Eval();
// 	GetTransformation6(s1,s2,m6FX);
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
      AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
      if (!magF) AliError("Magneticd field - not initialized");
      Double_t bz = magF->SolenoidField()/10.; //field in T

      cstream<<"Align"<<
	"run="<<fRun<<  // run
	"bz="<<bz<<
	"s1="<<s1<<     // reference sector
	"s2="<<s2<<     // sector to align
	"m6FX.="<<&m6FX<<   // tranformation matrix
	"m6.="<<&m6<<   // tranformation matrix
	"m9.="<<&m9<<   // 
	"m12.="<<&m12<<
	"chi2v12="<<chi2v12<<
	"chi2v9="<<chi2v9<<
	"chi2v6="<<chi2v6<<
	//
	"p6.="<<&param6Diff<<
	"p6s1.="<<&param6s1<<
	"p6s2.="<<&param6s2<<
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
Long64_t AliTPCcalibAlign::Merge(TCollection* const list) {
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
      if (align->fPoints[GetIndex(i,j)]<1) continue;
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
	    his0 = GetHisto((HistoType)itype,i,j,kTRUE);
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
      if (align->fPoints[GetIndex(i,j)]<1) continue;
      //
      //
      // fitter12
      f0 =  GetFitter12(i,j);
      f1 =  align->GetFitter12(i,j);
      if (f1){
	if (f0) f0->Add(f1);
	else {
	  fFitterArray12.AddAt(f1->Clone(),GetIndex(i,j));
	}
      }      
      //
      // fitter9
      f0 =  GetFitter9(i,j);
      f1 =  align->GetFitter9(i,j);
      if (f1){
	if (f0) f0->Add(f1);
	else { 
	  fFitterArray9.AddAt(f1->Clone(),GetIndex(i,j));
	}
      }      
      f0 =  GetFitter6(i,j);
      f1 =  align->GetFitter6(i,j);
      if (f1){
	if (f0) f0->Add(f1);
	else {
	  fFitterArray6.AddAt(f1->Clone(),GetIndex(i,j));
	}
      }   
    }
  }
  //
  // Add Kalman filter
  //
  for (Int_t i=0;i<36;i++){
    TMatrixD *par0 = (TMatrixD*)fArraySectorIntParam.At(i);
    if (!par0){
      MakeSectorKalman();
      par0 = (TMatrixD*)fArraySectorIntParam.At(i);      
    }
    TMatrixD *par1 = (TMatrixD*)align->fArraySectorIntParam.At(i);
    if (!par1) continue;
    //
    TMatrixD *cov0 = (TMatrixD*)fArraySectorIntCovar.At(i);
    TMatrixD *cov1 = (TMatrixD*)align->fArraySectorIntCovar.At(i);
    UpdateSectorKalman(*par0,*cov0,*par1,*cov1);
  }
  if (!fSectorParamA){
    MakeKalman();
  }
  if (align->fSectorParamA){
    UpdateKalman(*fSectorParamA,*fSectorCovarA,*align->fSectorParamA,*align->fSectorCovarA);
    UpdateKalman(*fSectorParamC,*fSectorCovarC,*align->fSectorParamC,*align->fSectorCovarC);
  }
  if (!fClusterDelta[1]) MakeResidualHistos();

  for (Int_t i=0; i<6; i++){
    if (i==0 || i==3){
      delete fClusterDelta[i];   // memory problem do not fit into memory
      fClusterDelta[i]=0;        // 
      delete align->fClusterDelta[i];   // memory problem do not fit into memory
      align->fClusterDelta[i]=0;        // 
    }
    if (i==3) continue;  // skip constrained histo z
    if (i==0) continue;  // skip non constrained histo y
    if (align->fClusterDelta[i]) fClusterDelta[i]->Add(align->fClusterDelta[i]);
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
	weight/=(errors[4]*errors[4]+sysError*sysError);  // wieghting with error in Y
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


 Int_t  AliTPCcalibAlign::RefitLinear(const AliTPCseed * track, Int_t isec, Double_t *fitParam, Int_t refSector,  TMatrixD &tparam, TMatrixD&tcovar, Double_t xRef, Bool_t both){
  //
  // Refit tracklet linearly using clusters at  given sector isec
  // Clusters are rotated to the  reference frame of sector refSector
  // 
  // fit parameters and errors retruning in the fitParam
  //
  // seed      - acces to the original clusters
  // isec      - sector to be refited
  // fitParam  - 
  //           0  lx             
  //           1  ly
  //           2  dy/dz
  //           3  lz
  //           4  dz/dx
  //           5  sx 
  //           6  sy
  //           7  sdydx
  //           8  sz
  //           9  sdzdx
  // ref sector is the sector defining ref frame - rotation
  // return value - number of used clusters

  const Int_t      kMinClusterF=15;
  const  Int_t     kdrow1 =10;        // rows to skip at the end      
  const  Int_t     kdrow0 =3;        // rows to skip at beginning  
  const  Float_t   kedgeyIn=2.5;
  const  Float_t   kedgeyOut=4.0;
  const  Float_t   kMaxDist=5;       // max distance -in sigma 
  const  Float_t   kMaxCorrY=0.05;    // max correction
  //
  Double_t dalpha = 0;
  if ((refSector%18)!=(isec%18)){
    dalpha = -((refSector%18)-(isec%18))*TMath::TwoPi()/18.;
  }
  Double_t ca = TMath::Cos(dalpha);
  Double_t sa = TMath::Sin(dalpha);
  //
  //
  AliTPCPointCorrection * corr =  AliTPCPointCorrection::Instance();
  //
  // full track fit parameters
  // 
  static TLinearFitter fyf(2,"pol1");   // change to static - suggestion of calgrind - 30 % of time
  static TLinearFitter fzf(2,"pol1");   // relative to time of given class
  TVectorD pyf(2), peyf(2),pzf(2), pezf(2);
  TMatrixD covY(4,4),covZ(4,4);
  Double_t chi2FacY =1;
  Double_t chi2FacZ =1;
  Int_t nf=0;
  //
  //
  //
  Float_t erry=0.1;   // initial cluster error estimate
  Float_t errz=0.1;   // initial cluster error estimate
  for (Int_t iter=0; iter<2; iter++){
    fyf.ClearPoints();
    fzf.ClearPoints();
    for (Int_t irow=kdrow0;irow<159-kdrow1;irow++) {
      AliTPCclusterMI *c=track->GetClusterPointer(irow);
      if (!c) continue;      
      //
      if (c->GetDetector()%36!=(isec%36)) continue;
      if (!both && c->GetDetector()!=isec) continue;

      if (c->GetRow()<kdrow0) continue;
      //cluster position in reference frame 
      Double_t lxR   =   ca*c->GetX()-sa*c->GetY();  
      Double_t lyR   =  +sa*c->GetX()+ca*c->GetY();
      Double_t lzR   =  c->GetZ();

      Double_t dx = lxR -xRef;   // distance to reference X
      Double_t x[2]={dx, dx*dx};

      Double_t yfitR  =    pyf[0]+pyf[1]*dx;  // fit value Y in ref frame
      Double_t zfitR  =    pzf[0]+pzf[1]*dx;  // fit value Z in ref frame
      //
      Double_t yfit   =  -sa*lxR + ca*yfitR;  // fit value Y in local frame
      //
      if (iter==0 &&c->GetType()<0) continue;
      if (iter>0){	
	if (TMath::Abs(lyR-yfitR)>kMaxDist*erry) continue;
	if (TMath::Abs(lzR-zfitR)>kMaxDist*errz) continue;
	Double_t dedge =  c->GetX()*TMath::Tan(TMath::Pi()/18.)-TMath::Abs(yfit);
	if (isec<36 && dedge<kedgeyIn) continue;
	if (isec>35 && dedge<kedgeyOut) continue;
	Double_t corrtrY =  
	  corr->RPhiCOGCorrection(isec,c->GetRow(), c->GetPad(),
				  c->GetY(),yfit, c->GetZ(), pyf[1], c->GetMax(),2.5);
	Double_t corrclY =  
	  corr->RPhiCOGCorrection(isec,c->GetRow(), c->GetPad(),
				  c->GetY(),c->GetY(), c->GetZ(), pyf[1], c->GetMax(),2.5);
	if (TMath::Abs((corrtrY+corrclY)*0.5)>kMaxCorrY) continue;
	if (TMath::Abs(corrtrY)>kMaxCorrY) continue;
      }
      fyf.AddPoint(x,lyR,erry);
      fzf.AddPoint(x,lzR,errz);
    }
    nf = fyf.GetNpoints();
    if (nf<kMinClusterF) return 0;   // not enough points - skip 
    fyf.Eval(); 
    fyf.GetParameters(pyf); 
    fyf.GetErrors(peyf);
    fzf.Eval(); 
    fzf.GetParameters(pzf); 
    fzf.GetErrors(pezf);
    chi2FacY = TMath::Sqrt(fyf.GetChisquare()/(fyf.GetNpoints()-2.));
    chi2FacZ = TMath::Sqrt(fzf.GetChisquare()/(fzf.GetNpoints()-2.));
    peyf[0]*=chi2FacY;
    peyf[1]*=chi2FacY;
    pezf[0]*=chi2FacZ;
    pezf[1]*=chi2FacZ;
    erry*=chi2FacY;
    errz*=chi2FacZ;
    fyf.GetCovarianceMatrix(covY);
    fzf.GetCovarianceMatrix(covZ);
    for (Int_t i0=0;i0<2;i0++)
      for (Int_t i1=0;i1<2;i1++){
	covY(i0,i1)*=chi2FacY*chi2FacY;
	covZ(i0,i1)*=chi2FacZ*chi2FacZ;
      }
  }
  fitParam[0] = xRef;
  //
  fitParam[1] = pyf[0];
  fitParam[2] = pyf[1];
  fitParam[3] = pzf[0];
  fitParam[4] = pzf[1];
  //
  fitParam[5] = 0;
  fitParam[6] = peyf[0];
  fitParam[7] = peyf[1];
  fitParam[8] = pezf[0];
  fitParam[9] = pezf[1];
  //
  //
  tparam(0,0) = pyf[0];
  tparam(1,0) = pyf[1];
  tparam(2,0) = pzf[0];
  tparam(3,0) = pzf[1];
  //
  tcovar(0,0) = covY(0,0);
  tcovar(1,1) = covY(1,1);
  tcovar(1,0) = covY(1,0);
  tcovar(0,1) = covY(0,1); 
  tcovar(2,2) = covZ(0,0);
  tcovar(3,3) = covZ(1,1);
  tcovar(3,2) = covZ(1,0);
  tcovar(2,3) = covZ(0,1);
  return nf;
}

void  AliTPCcalibAlign::UpdateAlignSector(const AliTPCseed * track,Int_t isec){
  //
  // Update Kalman filter of Alignment 
  //       IROC - OROC quadrants
  //
  if (!fClusterDelta[0])  MakeResidualHistos();
  const Int_t kMinClusterF=40;
  const Int_t kMinClusterQ=10;
  //
  const  Int_t     kdrow1Fit =5;         // rows to skip from fit at the end      
  const  Int_t     kdrow0Fit =10;        // rows to skip from fit at beginning  
  const  Float_t   kedgey=3.0;
  const  Float_t   kMaxDist=1;
  const  Float_t   kMaxCorrY=0.05;
  const  Float_t   kPRFWidth = 0.6;   //cut 2 sigma of PRF
  isec = isec%36;     // use the hardware numbering
  //
  //
  AliTPCPointCorrection * corr =  AliTPCPointCorrection::Instance();
  //
  // full track fit parameters
  // 
  static TLinearFitter fyf(2,"pol1");   // make it static - too much time for comiling of formula
  static TLinearFitter fzf(2,"pol1");   // calgrind recomendation
  TVectorD pyf(2), peyf(2),pzf(2), pezf(2);
  TVectorD pyfc(2),pzfc(2);
  Int_t nf=0;
  //
  // make full fit as reference
  //
  for (Int_t iter=0; iter<2; iter++){
    fyf.ClearPoints();
    fzf.ClearPoints();
    for (Int_t irow=kdrow0Fit;irow<159-kdrow1Fit;irow++) {
      AliTPCclusterMI *c=track->GetClusterPointer(irow);
      if (!c) continue;
      if ((c->GetDetector()%36)!=isec) continue;
      if (c->GetRow()<kdrow0Fit) continue;
      Double_t dx = c->GetX()-fXmiddle;
      Double_t x[2]={dx, dx*dx};
      if (iter==0 &&c->GetType()<0) continue;
      if (iter==1){	
	Double_t yfit  =  pyf[0]+pyf[1]*dx;
	Double_t zfit  =  pzf[0]+pzf[1]*dx;
	Double_t dedge =  c->GetX()*TMath::Tan(TMath::Pi()/18.)-TMath::Abs(yfit);
	if (TMath::Abs(c->GetY()-yfit)>kMaxDist) continue;
	if (TMath::Abs(c->GetZ()-zfit)>kMaxDist) continue;
	if (dedge<kedgey) continue;
	Double_t corrtrY =  
	  corr->RPhiCOGCorrection(c->GetDetector(),c->GetRow(), c->GetPad(),
				  c->GetY(),yfit, c->GetZ(), pyf[1], c->GetMax(),2.5);
	if (TMath::Abs(corrtrY)>kMaxCorrY) continue;
      }
      fyf.AddPoint(x,c->GetY(),0.1);
      fzf.AddPoint(x,c->GetZ(),0.1);
    }
    nf = fyf.GetNpoints();
    if (nf<kMinClusterF) return;   // not enough points - skip 
    fyf.Eval(); 
    fyf.GetParameters(pyf); 
    fyf.GetErrors(peyf);
    fzf.Eval(); 
    fzf.GetParameters(pzf); 
    fzf.GetErrors(pezf);
  }
  //
  //
  //
  TVectorD vecX(2*nf+kdrow0Fit+kdrow1Fit+5);          // x         vector
  TVectorD vecY(2*nf+kdrow0Fit+kdrow1Fit+5);          // residuals vector
  TVectorD vecZ(2*nf+kdrow0Fit+kdrow1Fit+5);                              // residuals vector
  TVectorD vPosG(3);                  //vertex position
  TVectorD vPosL(3);                 // vertex position in the TPC local system
  TVectorF vImpact(2);               //track impact parameter
  Double_t tofSignal= fCurrentTrack->GetTOFsignal();      // tof signal
  TVectorF tpcPosG(3);                                    // global position of track at the middle of fXmiddle
  Double_t lphi = TMath::ATan2(pyf[0],fXmiddle);          // expected local phi angle - if vertex at 0
  Double_t gphi = 2.*TMath::Pi()*(isec%18+0.5)/18.+lphi;  // expected global phi if vertex at 0
  Double_t ky   = pyf[0]/fXmiddle;
  Double_t kyE  =0, kzE=0;    // ky and kz expected
  Double_t alpha =2.*TMath::Pi()*(isec%18+0.5)/18.;
  Double_t scos=TMath::Cos(alpha);
  Double_t ssin=TMath::Sin(alpha);
  const AliESDVertex* vertex = fCurrentEvent->GetPrimaryVertexTracks();
  vertex->GetXYZ(vPosG.GetMatrixArray());
  fCurrentTrack->GetImpactParameters(vImpact[0],vImpact[1]);  // track impact parameters
  //
  tpcPosG[0]= scos*fXmiddle-ssin*pyf[0];
  tpcPosG[1]=+ssin*fXmiddle+scos*pyf[0];
  tpcPosG[2]=pzf[0];
  vPosL[0]= scos*vPosG[0]+ssin*vPosG[1];
  vPosL[1]=-ssin*vPosG[0]+scos*vPosG[1];
  vPosL[2]= vPosG[2];
  kyE = (pyf[0]-vPosL[1])/(fXmiddle-vPosL[0]);
  kzE = (pzf[0]-vPosL[2])/(fXmiddle-vPosL[0]);
  //
  // get constrained parameters
  //
  Double_t xvertex=vPosL[0]-fXmiddle;
  fyf.AddPoint(&xvertex,vPosL[1], 0.1+TMath::Abs(vImpact[0]));
  fzf.AddPoint(&xvertex,vPosL[2], 0.1+TMath::Abs(vImpact[1]));
  fyf.Eval(); 
  fyf.GetParameters(pyfc); 
  fzf.Eval(); 
  fzf.GetParameters(pzfc); 
  //
  //
  // Make Fitters and params for 5 fitters
  // 1-4 OROC quadrants 
  //   0 IROC
  //
  static TLinearFitter *fittersY[5]={0,0,0,0,0};   // calgrind recomendation - fater to clear points
  static TLinearFitter *fittersZ[5]={0,0,0,0,0};   // than create the fitter
  if (fittersY[0]==0){
    for (Int_t i=0;i<5;i++) {
      fittersY[i] = new TLinearFitter(2,"pol1");
      fittersZ[i] = new TLinearFitter(2,"pol1");
    }
  }
  //
  Int_t npoints[5];
  TVectorD paramsY[5];
  TVectorD errorsY[5];
  TMatrixD covY[5];
  Double_t chi2CY[5];
  TVectorD paramsZ[5];
  TVectorD errorsZ[5];
  TMatrixD covZ[5];
  Double_t chi2CZ[5];
  for (Int_t i=0;i<5;i++) {
    npoints[i]=0;
    paramsY[i].ResizeTo(2);
    errorsY[i].ResizeTo(2);
    covY[i].ResizeTo(2,2);
    paramsZ[i].ResizeTo(2);
    errorsZ[i].ResizeTo(2);
    covZ[i].ResizeTo(2,2);
    fittersY[i]->ClearPoints();
    fittersZ[i]->ClearPoints();
  }
  //
  // Update fitters
  //
  Int_t countRes=0;
  for (Int_t irow=0;irow<159;irow++) {
    AliTPCclusterMI *c=track->GetClusterPointer(irow);
    if (!c) continue;
    if ((c->GetDetector()%36)!=isec) continue;
    Double_t dx = c->GetX()-fXmiddle;
    Double_t x[2]={dx, dx*dx};
    Double_t yfit  =  pyf[0]+pyf[1]*dx;
    Double_t zfit  =  pzf[0]+pzf[1]*dx;
    Double_t yfitC  =  pyfc[0]+pyfc[1]*dx;
    Double_t zfitC  =  pzfc[0]+pzfc[1]*dx;
    Double_t dedge =  c->GetX()*TMath::Tan(TMath::Pi()/18.)-TMath::Abs(yfit);
    if (TMath::Abs(c->GetY()-yfit)>kMaxDist) continue;
    if (TMath::Abs(c->GetZ()-zfit)>kMaxDist) continue;
    if (dedge<kedgey) continue;
    Double_t corrtrY =  
      corr->RPhiCOGCorrection(c->GetDetector(),c->GetRow(), c->GetPad(),
			      c->GetY(),yfit, c->GetZ(), pyf[1], c->GetMax(),2.5);
    if (TMath::Abs(corrtrY)>kMaxCorrY) continue;  
    //
    vecX[countRes]=c->GetX();
    vecY[countRes]=c->GetY()-yfit;
    vecZ[countRes]=c->GetZ()-zfit;
    countRes++;
    //
    // Fill THnSparse cluster residuals
    // use only primary candidates with ITS signal
    if (nf>100&&fCurrentTrack->IsOn(0x4)&&TMath::Abs(vImpact[0])<1&&TMath::Abs(vImpact[1])<1){    
      Double_t resVector[5];
      resVector[1]= 9.*gphi/TMath::Pi();
      resVector[2]= c->GetX();
      resVector[3]= c->GetY()/c->GetX();
      resVector[4]= c->GetZ()/c->GetX();
      //
      resVector[0]= c->GetY()-yfit;
      //fClusterDelta[0]->Fill(resVector);
      resVector[0]= c->GetZ()-zfit;
      fClusterDelta[1]->Fill(resVector);
      //
      resVector[0]= c->GetY()-yfitC;
      fClusterDelta[2]->Fill(resVector);
      resVector[0]= c->GetZ()-zfitC;
      //fClusterDelta[3]->Fill(resVector);

    }
    if (c->GetRow()<kdrow0Fit) continue;      
    if (c->GetRow()>159-kdrow1Fit) continue;      
    //

    if (c->GetDetector()>35){      
      if (c->GetX()<fXquadrant){
	if (yfit<-kPRFWidth) fittersY[1]->AddPoint(x,c->GetY(),0.1);
	if (yfit<-kPRFWidth) fittersZ[1]->AddPoint(x,c->GetZ(),0.1);
	if (yfit>kPRFWidth) fittersY[2]->AddPoint(x,c->GetY(),0.1);
	if (yfit>kPRFWidth) fittersZ[2]->AddPoint(x,c->GetZ(),0.1);
      }
      if (c->GetX()>fXquadrant){
	if (yfit<-kPRFWidth) fittersY[3]->AddPoint(x,c->GetY(),0.1);
	if (yfit<-kPRFWidth) fittersZ[3]->AddPoint(x,c->GetZ(),0.1);
	if (yfit>kPRFWidth) fittersY[4]->AddPoint(x,c->GetY(),0.1);
	if (yfit>kPRFWidth) fittersZ[4]->AddPoint(x,c->GetZ(),0.1);
      }
    }
    if (c->GetDetector()<36){
      fittersY[0]->AddPoint(x,c->GetY(),0.1);
      fittersZ[0]->AddPoint(x,c->GetZ(),0.1);
    }
  }
  //
  //
  //
  for (Int_t i=0;i<5;i++) {
    npoints[i] = fittersY[i]->GetNpoints();
    if (npoints[i]>=kMinClusterQ){
      // Y fit 
      fittersY[i]->Eval();
      Double_t chi2FacY = TMath::Sqrt(fittersY[i]->GetChisquare()/(fittersY[i]->GetNpoints()-2));
      chi2CY[i]=chi2FacY;
      fittersY[i]->GetParameters(paramsY[i]);
      fittersY[i]->GetErrors(errorsY[i]);
      fittersY[i]->GetCovarianceMatrix(covY[i]);
      // renormalize errors
      errorsY[i][0]*=chi2FacY;
      errorsY[i][1]*=chi2FacY;
      covY[i](0,0)*=chi2FacY*chi2FacY;
      covY[i](0,1)*=chi2FacY*chi2FacY;
      covY[i](1,0)*=chi2FacY*chi2FacY;
      covY[i](1,1)*=chi2FacY*chi2FacY;
      // Z fit
      fittersZ[i]->Eval();
      Double_t chi2FacZ = TMath::Sqrt(fittersZ[i]->GetChisquare()/(fittersZ[i]->GetNpoints()-2));
      chi2CZ[i]=chi2FacZ;
      fittersZ[i]->GetParameters(paramsZ[i]);
      fittersZ[i]->GetErrors(errorsZ[i]);
      fittersZ[i]->GetCovarianceMatrix(covZ[i]);
      // renormalize errors
      errorsZ[i][0]*=chi2FacZ;
      errorsZ[i][1]*=chi2FacZ;
      covZ[i](0,0)*=chi2FacZ*chi2FacZ;
      covZ[i](0,1)*=chi2FacZ*chi2FacZ;
      covZ[i](1,0)*=chi2FacZ*chi2FacZ;
      covZ[i](1,1)*=chi2FacZ*chi2FacZ;      
    }
  }
  //
  // void UpdateSectorKalman
  //
  for (Int_t i0=0;i0<5;i0++){
    for (Int_t i1=i0+1;i1<5;i1++){
      if(npoints[i0]<kMinClusterQ) continue;
      if(npoints[i1]<kMinClusterQ) continue;
      TMatrixD v0(4,1),v1(4,1);        // measurement
      TMatrixD cov0(4,4),cov1(4,4);    // covariance
      //
      v0(0,0)= paramsY[i0][0];         v1(0,0)= paramsY[i1][0]; 
      v0(1,0)= paramsY[i0][1];         v1(1,0)= paramsY[i1][1]; 
      v0(2,0)= paramsZ[i0][0];         v1(2,0)= paramsZ[i1][0]; 
      v0(3,0)= paramsZ[i0][1];         v1(3,0)= paramsZ[i1][1]; 
      //covariance i0
      cov0(0,0) = covY[i0](0,0);
      cov0(1,1) = covY[i0](1,1);
      cov0(1,0) = covY[i0](1,0);
      cov0(0,1) = covY[i0](0,1); 
      cov0(2,2) = covZ[i0](0,0);
      cov0(3,3) = covZ[i0](1,1);
      cov0(3,2) = covZ[i0](1,0);
      cov0(2,3) = covZ[i0](0,1);
      //covariance i1
      cov1(0,0) = covY[i1](0,0);
      cov1(1,1) = covY[i1](1,1);
      cov1(1,0) = covY[i1](1,0);
      cov1(0,1) = covY[i1](0,1); 
      cov1(2,2) = covZ[i1](0,0);
      cov1(3,3) = covZ[i1](1,1);
      cov1(3,2) = covZ[i1](1,0);
      cov1(2,3) = covZ[i1](0,1);
      //
      // And now update
      //
      if (TMath::Abs(pyf[1])<0.8){ //angular cut	
	UpdateSectorKalman(isec, i0,i1, &v0,&cov0,&v1,&cov1);
      }
    }
  }

  //
  // Dump debug information
  //
  if (fStreamLevel>0){
    // get vertex position
     //
    TTreeSRedirector *cstream = GetDebugStreamer();  

    
    if (cstream){      
      for (Int_t i0=0;i0<5;i0++){
	for (Int_t i1=i0;i1<5;i1++){
	  if (i0==i1) continue;
	  if(npoints[i0]<kMinClusterQ) continue;
	  if(npoints[i1]<kMinClusterQ) continue;
	  (*cstream)<<"sectorAlign"<<
	    "run="<<fRun<<              //  run number
	    "event="<<fEvent<<          //  event number
	    "time="<<fTime<<            //  time stamp of event
	    "trigger="<<fTrigger<<      //  trigger
	    "triggerClass="<<&fTriggerClass<<      //  trigger
	    "mag="<<fMagF<<             //  magnetic field	  
	    "isec="<<isec<<             //  current sector 
	    //
	    "xref="<<fXmiddle<<         // reference X
	    "vPosG.="<<&vPosG<<        // vertex position in global system
	    "vPosL.="<<&vPosL<<        // vertex position in local  system
	    "vImpact.="<<&vImpact<<   // track impact parameter
	    "tofSignal="<<tofSignal<<   // tof signal
	    "tpcPosG.="<<&tpcPosG<<   // global position of track at the middle of fXmiddle
	    "lphi="<<lphi<<             // expected local phi angle - if vertex at 0
	    "gphi="<<gphi<<             // expected global phi if vertex at 0
	    "ky="<<ky<<
	    "kyE="<<kyE<<               // expect ky - assiming pirmary track
	    "kzE="<<kzE<<               // expected kz - assuming primary tracks
	    "salpha="<<alpha<<          // sector alpha
	    "scos="<<scos<<              // tracking cosinus
	    "ssin="<<ssin<<              // tracking sinus
	    //
	    // full fit
	    //
	    "nf="<<nf<<                 //  total number of points
	    "pyf.="<<&pyf<<             //  full OROC fit y
	    "pzf.="<<&pzf<<             //  full OROC fit z
	    "vX.="<<&vecX<<              //  x cluster
	    "vY.="<<&vecY<<              //  y residual cluster
	    "vZ.="<<&vecZ<<              //  z residual cluster
	    // quadrant and IROC fit
	    "i0="<<i0<<                 // quadrant number
	    "i1="<<i1<<                 
	    "n0="<<npoints[i0]<<        // number of points
	    "n1="<<npoints[i1]<<
	    //
	    "py0.="<<&paramsY[i0]<<       // parameters
	    "py1.="<<&paramsY[i1]<< 
	    "ey0.="<<&errorsY[i0]<<       // errors
	    "ey1.="<<&errorsY[i1]<<
	    "chiy0="<<chi2CY[i0]<<       // chi2s
	    "chiy1="<<chi2CY[i1]<<
	    //
	    "pz0.="<<&paramsZ[i0]<<       // parameters
	    "pz1.="<<&paramsZ[i1]<< 
	    "ez0.="<<&errorsZ[i0]<<       // errors
	    "ez1.="<<&errorsZ[i1]<<
	    "chiz0="<<chi2CZ[i0]<<       // chi2s
	    "chiz1="<<chi2CZ[i1]<<
	    "\n";
	}    
      }
    }
  }
}

void AliTPCcalibAlign::UpdateSectorKalman(Int_t sector, Int_t quadrant0, Int_t quadrant1,  TMatrixD *const p0, TMatrixD *const c0, TMatrixD *const p1, TMatrixD *const c1 ){
  //
  //
  // tracks are refitted at sector middle
  //
  if (fArraySectorIntParam.At(0)==NULL) MakeSectorKalman();
  //
  //
  static TMatrixD matHk(4,30);    // vector to mesurement
  static TMatrixD measR(4,4);     // measurement error 
  //  static TMatrixD matQk(2,2);     // prediction noise vector
  static TMatrixD vecZk(4,1);     // measurement
  //
  static TMatrixD vecYk(4,1);     // Innovation or measurement residual
  static TMatrixD matHkT(30,4);   // helper matrix Hk transpose
  static TMatrixD matSk(4,4);     // Innovation (or residual) covariance
  static TMatrixD matKk(30,4);    // Optimal Kalman gain
  static TMatrixD mat1(30,30);    // update covariance matrix
  static TMatrixD covXk2(30,30);  // helper matrix
  //
  TMatrixD *vOrig = (TMatrixD*)(fArraySectorIntParam.At(sector));
  TMatrixD *cOrig = (TMatrixD*)(fArraySectorIntCovar.At(sector));
  //
  TMatrixD vecXk(*vOrig);             // X vector
  TMatrixD covXk(*cOrig);             // X covariance
  //
  //Unit matrix
  //
  for (Int_t i=0;i<30;i++)
    for (Int_t j=0;j<30;j++){
      mat1(i,j)=0;
      if (i==j) mat1(i,j)=1;
    }
  //
  //
  // matHk - vector to measurement
  //
  for (Int_t i=0;i<4;i++)
    for (Int_t j=0;j<30;j++){
      matHk(i,j)=0;
    }
  //
  // Measurement  
  // 0  - y
  // 1  - ky
  // 2  - z
  // 3  - kz
  
  matHk(0,6*quadrant1+4)  =  1.;            // delta y
  matHk(1,6*quadrant1+0)  =  1.;            // delta ky
  matHk(2,6*quadrant1+5)  =  1.;            // delta z
  matHk(3,6*quadrant1+1)  =  1.;            // delta kz
  // bug fix 24.02  - aware of sign in dx
  matHk(0,6*quadrant1+3)  =  -(*p0)(1,0);    // delta x to delta y  - through ky
  matHk(2,6*quadrant1+3)  =  -(*p0)(3,0);    // delta x to delta z  - thorugh kz
  matHk(2,6*quadrant1+2)  =  ((*p0)(0,0));  // y       to delta z  - through psiz
  //
  matHk(0,6*quadrant0+4)  =  -1.;           // delta y
  matHk(1,6*quadrant0+0)  =  -1.;           // delta ky
  matHk(2,6*quadrant0+5)  =  -1.;           // delta z
  matHk(3,6*quadrant0+1)  =  -1.;           // delta kz
  // bug fix 24.02 be aware of sign in dx
  matHk(0,6*quadrant0+3)  =  ((*p0)(1,0)); // delta x to delta y  - through ky
  matHk(2,6*quadrant0+3)  =  ((*p0)(3,0)); // delta x to delta z  - thorugh kz
  matHk(2,6*quadrant0+2)  =  -((*p0)(0,0)); // y       to delta z  - through psiz

  //
  //
  
  vecZk =(*p1)-(*p0);                 // measurement
  measR =(*c1)+(*c0);                 // error of measurement
  vecYk = vecZk-matHk*vecXk;          // Innovation or measurement residual
  //
  //
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  vecXk += matKk*vecYk;                    //  updated vector 
  covXk2= (mat1-(matKk*matHk));
  covXk =  covXk2*covXk;    
  //
  //
  (*cOrig)=covXk;
  (*vOrig)=vecXk;
}

void AliTPCcalibAlign::MakeSectorKalman(){
  //
  // Make a initial Kalman paramaters for IROC - Quadrants alignment
  //
  TMatrixD param(5*6,1);
  TMatrixD covar(5*6,5*6);
  //
  // Set inital parameters
  //
  for (Int_t ip=0;ip<5*6;ip++) param(ip,0)=0;  // mean alignment to 0
  //
  for (Int_t iq=0;iq<5;iq++){
    // Initial uncertinty
    covar(iq*6+0,iq*6+0) = 0.002*0.002;        // 2 mrad  
    covar(iq*6+1,iq*6+1) = 0.002*0.002;        // 2 mrad  rotation
    covar(iq*6+2,iq*6+2) = 0.002*0.002;        // 2 mrad 
    //
    covar(iq*6+3,iq*6+3) = 0.02*0.02;        // 0.2 mm  
    covar(iq*6+4,iq*6+4) = 0.02*0.02;        // 0.2 mm  translation
    covar(iq*6+5,iq*6+5) = 0.02*0.02;        // 0.2 mm 
  }

  for (Int_t isec=0;isec<36;isec++){
    fArraySectorIntParam.AddAt(param.Clone(),isec);
    fArraySectorIntCovar.AddAt(covar.Clone(),isec);
  }
}

void AliTPCcalibAlign::UpdateSectorKalman(TMatrixD &par0, TMatrixD &cov0, TMatrixD &par1, TMatrixD &cov1){
  //
  // Update kalman vector para0 with vector par1
  // Used for merging
  //
  static TMatrixD matHk(30,30);    // vector to mesurement
  static TMatrixD measR(30,30);     // measurement error 
  static TMatrixD vecZk(30,1);     // measurement
  //
  static TMatrixD vecYk(30,1);     // Innovation or measurement residual
  static TMatrixD matHkT(30,30);   // helper matrix Hk transpose
  static TMatrixD matSk(30,30);     // Innovation (or residual) covariance
  static TMatrixD matKk(30,30);    // Optimal Kalman gain
  static TMatrixD mat1(30,30);    // update covariance matrix
  static TMatrixD covXk2(30,30);  // helper matrix
  //
  TMatrixD vecXk(par0);             // X vector
  TMatrixD covXk(cov0);             // X covariance

  //
  //Unit matrix
  //
  for (Int_t i=0;i<30;i++)
    for (Int_t j=0;j<30;j++){
      mat1(i,j)=0;
      if (i==j) mat1(i,j)=1;
    }
  matHk = mat1;                       // unit matrix 
  //
  vecZk = par1;                       // measurement
  measR = cov1;                        // error of measurement
  vecYk = vecZk-matHk*vecXk;          // Innovation or measurement residual
  //
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  //matKk.Print();
  vecXk += matKk*vecYk;                    //  updated vector 
  covXk2= (mat1-(matKk*matHk));
  covXk =  covXk2*covXk;   
  CheckCovariance(covXk);
  CheckCovariance(cov1);
  //
  par0  = vecXk;                      // update measurement param
  cov0  = covXk;                      // update measurement covar
}

Double_t AliTPCcalibAlign::GetCorrectionSector(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t /*lz*/){
  //
  // Get position correction for given sector
  //

  TMatrixD * param = (TMatrixD*)fArraySectorIntParam.At(sector%36);
  if (!param) return 0;
  Int_t quadrant=0;
  if(lx>fXIO){
    if (lx<fXquadrant) {
      if (ly<0) quadrant=1;
      if (ly>0) quadrant=2;
    }
    if (lx>fXquadrant) {
      if (ly<0) quadrant=3;
      if (ly>0) quadrant=4;
    }
  }
  Double_t a10 = (*param)(quadrant*6+0,0);
  Double_t a20 = (*param)(quadrant*6+1,0);
  Double_t a21 = (*param)(quadrant*6+2,0);
  Double_t dx  = (*param)(quadrant*6+3,0);
  Double_t dy  = (*param)(quadrant*6+4,0);
  Double_t dz  = (*param)(quadrant*6+5,0);
  Double_t deltaX = lx-fXIO;
  if (coord==0) return dx;
  if (coord==1) return dy+deltaX*a10;
  if (coord==2) return dz+deltaX*a20+ly*a21;
  return 0;
}

Double_t AliTPCcalibAlign::SGetCorrectionSector(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz){
  //
  //
  //
  if (!Instance()) return 0;
  return Instance()->GetCorrectionSector(coord,sector,lx,ly,lz);
}

void AliTPCcalibAlign::MakeKalman(){
  //
  // Make a initial Kalman paramaters for sector Alignemnt
  //
  fSectorParamA = new TMatrixD(6*36+6,1);
  fSectorParamC = new TMatrixD(6*36+6,1);
  fSectorCovarA = new TMatrixD(6*36+6,6*36+6);
  fSectorCovarC = new TMatrixD(6*36+6,6*36+6);
  //
  // set starting parameters at 0
  //
  for (Int_t isec=0;isec<37;isec++)
    for (Int_t ipar=0;ipar<6;ipar++){
      (*fSectorParamA)(isec*6+ipar,0) =0;
      (*fSectorParamC)(isec*6+ipar,0) =0;
  }
  //
  // set starting covariance
  //
  for (Int_t isec=0;isec<36;isec++)
    for (Int_t ipar=0;ipar<6;ipar++){
      if (ipar<3){
	(*fSectorCovarA)(isec*6+ipar,isec*6+ipar) =0.002*0.002;   // 2 mrad
	(*fSectorCovarC)(isec*6+ipar,isec*6+ipar) =0.002*0.002;
      }
      if (ipar>=3){
	(*fSectorCovarA)(isec*6+ipar,isec*6+ipar) =0.02*0.02;     // 0.2 mm
	(*fSectorCovarC)(isec*6+ipar,isec*6+ipar) =0.02*0.02;
      }
  }
  (*fSectorCovarA)(36*6+0,36*6+0) =0.04;   // common shift y  up-up
  (*fSectorCovarA)(36*6+1,36*6+1) =0.04;   // common shift y  down-down
  (*fSectorCovarA)(36*6+2,36*6+2) =0.04;   // common shift y  up-down
  (*fSectorCovarA)(36*6+3,36*6+3) =0.004;   // common shift phi  up-up
  (*fSectorCovarA)(36*6+4,36*6+4) =0.004;   // common shift phi down-down
  (*fSectorCovarA)(36*6+5,36*6+5) =0.004;   // common shift phi up-down
  //
  (*fSectorCovarC)(36*6+0,36*6+0) =0.04;   // common shift y  up-up
  (*fSectorCovarC)(36*6+1,36*6+1) =0.04;   // common shift y  down-down
  (*fSectorCovarC)(36*6+2,36*6+2) =0.04;   // common shift y  up-down
  (*fSectorCovarC)(36*6+3,36*6+3) =0.004;   // common shift phi  up-up
  (*fSectorCovarC)(36*6+4,36*6+4) =0.004;   // common shift phi down-down
  (*fSectorCovarC)(36*6+5,36*6+5) =0.004;   // common shift phi up-down
}

void AliTPCcalibAlign::UpdateKalman(Int_t sector0, Int_t sector1,  TMatrixD &p0, TMatrixD &c0, TMatrixD &p1, TMatrixD &c1){
  //
  // Update Kalman parameters
  // Note numbering from 0..36  0..17 IROC 18..35 OROC
  // 
  //
  if (fSectorParamA==NULL) MakeKalman();
  if (CheckCovariance(c0)>0) return;
  if (CheckCovariance(c1)>0) return;
  const Int_t nelem = 6*36+6;
  //
  //
  static TMatrixD matHk(4,nelem);    // vector to mesurement
  static TMatrixD measR(4,4);     // measurement error 
  static TMatrixD vecZk(4,1);     // measurement
  //
  static TMatrixD vecYk(4,1);     // Innovation or measurement residual
  static TMatrixD matHkT(nelem,4);   // helper matrix Hk transpose
  static TMatrixD matSk(4,4);     // Innovation (or residual) covariance
  static TMatrixD matKk(nelem,4);    // Optimal Kalman gain
  static TMatrixD mat1(nelem,nelem);    // update covariance matrix
  static TMatrixD covXk2(nelem,nelem);  // helper matrix
  //
  TMatrixD *vOrig = 0;
  TMatrixD *cOrig = 0;
  vOrig = (sector0%36>=18) ? fSectorParamA:fSectorParamC;
  cOrig = (sector0%36>=18) ? fSectorCovarA:fSectorCovarC;
  //
  Int_t sec0= sector0%18;
  Int_t sec1= sector1%18;
  if (sector0>35) sec0+=18;
  if (sector1>35) sec1+=18;
  //
  TMatrixD vecXk(*vOrig);             // X vector
  TMatrixD covXk(*cOrig);             // X covariance
  //
  //Unit matrix
  //
  for (Int_t i=0;i<nelem;i++)
    for (Int_t j=0;j<nelem;j++){
      mat1(i,j)=0;
      if (i==j) mat1(i,j)=1;
    }
  //
  //
  // matHk - vector to measurement
  //
  for (Int_t i=0;i<4;i++)
    for (Int_t j=0;j<nelem;j++){
      matHk(i,j)=0;
    }
  //
  // Measurement  
  // 0  - y
  // 1  - ky
  // 2  - z
  // 3  - kz
  
  matHk(0,6*sec1+4)  =  1.;            // delta y
  matHk(1,6*sec1+0)  =  1.;            // delta ky
  matHk(2,6*sec1+5)  =  1.;            // delta z
  matHk(3,6*sec1+1)  =  1.;            // delta kz
  matHk(0,6*sec1+3)  =  p0(1,0);    // delta x to delta y  - through ky
  matHk(2,6*sec1+3)  =  p0(3,0);    // delta x to delta z  - thorugh kz
  matHk(2,6*sec1+2)  =  p0(0,0);    // y       to delta z  - through psiz
  //
  matHk(0,6*sec0+4)  =  -1.;           // delta y
  matHk(1,6*sec0+0)  =  -1.;           // delta ky
  matHk(2,6*sec0+5)  =  -1.;           // delta z
  matHk(3,6*sec0+1)  =  -1.;           // delta kz
  matHk(0,6*sec0+3)  =  -p0(1,0); // delta x to delta y  - through ky
  matHk(2,6*sec0+3)  =  -p0(3,0); // delta x to delta z  - thorugh kz
  matHk(2,6*sec0+2)  =  -p0(0,0); // y       to delta z  - through psiz

  Int_t dsec = (sector1%18)-(sector0%18);
  if (dsec<-2) dsec+=18; 
  if (TMath::Abs(dsec)==1){
    //
    // Left right systematic fit part
    //
    Double_t dir = 0;
    if (dsec>0) dir= 1.;
    if (dsec<0) dir=-1.;
    if (sector0>35&&sector1>35){
      matHk(0,36*6+0)=dir; 
      matHk(1,36*6+3+0)=dir; 
    }
    if (sector0<36&&sector1<36){
      matHk(0,36*6+1)=dir; 
      matHk(1,36*6+3+1)=dir; 
    }
    if (sector0<36&&sector1>35){
      matHk(0,36*6+2)=dir; 
      matHk(1,36*6+3+2)=dir; 
    }
    if (sector0>35&&sector1<36){
      matHk(0,36*6+2)=-dir; 
      matHk(1,36*6+3+2)=-dir; 
    }
  }
  //
  //  
  vecZk =(p1)-(p0);                 // measurement
  measR =(c1)+(c0);                 // error of measurement
  vecYk = vecZk-matHk*vecXk;          // Innovation or measurement residual
  //
  //
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  vecXk += matKk*vecYk;                    //  updated vector 
  covXk2= (mat1-(matKk*matHk));
  covXk =  covXk2*covXk;    

  if (CheckCovariance(covXk)>0) return;

  //
  //
  (*cOrig)=covXk;
  (*vOrig)=vecXk;
}


void AliTPCcalibAlign::UpdateKalman(TMatrixD &par0, TMatrixD &cov0, TMatrixD &par1, TMatrixD &cov1){
  //
  // Update kalman vector para0 with vector par1
  // Used for merging
  //
  Int_t nelem = 6*36+6;
  static TMatrixD matHk(nelem,nelem);    // vector to mesurement
  static TMatrixD measR(nelem,nelem);     // measurement error 
  static TMatrixD vecZk(nelem,1);     // measurement
  //
  static TMatrixD vecYk(nelem,1);     // Innovation or measurement residual
  static TMatrixD matHkT(nelem,nelem);   // helper matrix Hk transpose
  static TMatrixD matSk(nelem,nelem);     // Innovation (or residual) covariance
  static TMatrixD matKk(nelem,nelem);    // Optimal Kalman gain
  static TMatrixD mat1(nelem,nelem);    // update covariance matrix
  static TMatrixD covXk2(nelem,nelem);  // helper matrix
  //
  TMatrixD vecXk(par0);             // X vector
  TMatrixD covXk(cov0);             // X covariance

  //
  //Unit matrix
  //
  for (Int_t i=0;i<nelem;i++)
    for (Int_t j=0;j<nelem;j++){
      mat1(i,j)=0;
      if (i==j) mat1(i,j)=1;
    }
  matHk = mat1;                       // unit matrix 
  //
  vecZk = par1;                       // measurement
  measR = cov1;                        // error of measurement
  vecYk = vecZk-matHk*vecXk;          // Innovation or measurement residual
  //
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  //matKk.Print();
  vecXk += matKk*vecYk;                    //  updated vector 
  covXk2= (mat1-(matKk*matHk));
  covXk =  covXk2*covXk;
  //
  CheckCovariance(cov0);
  CheckCovariance(cov1);
  CheckCovariance(covXk);
  //
  par0  = vecXk;                      // update measurement param
  cov0  = covXk;                      // update measurement covar
}




Int_t  AliTPCcalibAlign::CheckCovariance(TMatrixD &covar){
  //
  // check the consistency of covariance matrix
  //
  Int_t ncols = covar.GetNcols();
  Int_t nrows= covar.GetNrows();
  const Float_t kEpsilon = 0.0001;
  Int_t nerrors =0;
  //
  //
  //
  if (nrows!=ncols) {
    printf("Error 0 - wrong matrix\n");
    nerrors++;
  }
  //
  // 1. Check that the non diagonal elements
  //
  for (Int_t i0=0;i0<nrows;i0++)
    for (Int_t i1=i0+1;i1<ncols;i1++){      
      Double_t r0 = covar(i0,i1)/TMath::Sqrt(covar(i0,i0)*covar(i1,i1));
      Double_t r1 = covar(i1,i0)/TMath::Sqrt(covar(i0,i0)*covar(i1,i1));
      if (TMath::Abs(r0-r1)>kEpsilon){
	printf("Error 1 - non symetric matrix %d\t%d\t%f",i0,i1,r1-r0);	    
	nerrors++;
      }
      if (TMath::Abs(r0)>=1){
	printf("Error 2 - Wrong correlation %d\t%d\t%f\n",i0,i1,r0);
	nerrors++;	    
      }     
      if (TMath::Abs(r1)>=1){
	printf("Error 3 - Wrong correlation %d\t%d\t%f\n",i0,i1,r1);
	nerrors++;
      }     
    }
  return nerrors;
}



void AliTPCcalibAlign::MakeReportDy(TFile *output){
  //
  // Draw histogram of dY
  //
  Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};

  AliTPCcalibAlign *align = this;
  TVectorD vecOOP(36);
  TVectorD vecOOM(36);
  TVectorD vecOIP(36);
  TVectorD vecOIM(36);
  TVectorD vecOIS(36);
  TVectorD vecSec(36);
  TCanvas * cOROCdy = new TCanvas("OROC dy","OROC dy",900,600);
  cOROCdy->Divide(6,6);
  TCanvas * cIROCdy = new TCanvas("IROC dy","IROC dy",900,600);
  cIROCdy->Divide(6,6);
  TCanvas * cDy = new TCanvas("Dy","Dy",600,700);
  cDy->Divide(1,2);
  for (Int_t isec=0;isec<36;isec++){
    Bool_t isDraw=kFALSE;
    vecSec(isec)=isec;
    cOROCdy->cd(isec+1);
    Int_t secPlus = (isec%18==17)? isec-17:isec+1;
    Int_t secMinus= (isec%18==0) ? isec+17:isec-1;
    printf("%d\t%d\t%d\n",isec,secPlus,secMinus);
    TH1 * hisOOP= align->GetHisto(AliTPCcalibAlign::kY,isec+36,secPlus+36);
    TH1 * hisOOM= align->GetHisto(AliTPCcalibAlign::kY,isec+36,secMinus+36);    
    TH1 * hisOIS= align->GetHisto(AliTPCcalibAlign::kY,isec+36,isec);    

    if (hisOIS) {
      hisOIS = (TH1F*)(hisOIS->Clone());
      hisOIS->SetDirectory(0);
      hisOIS->Scale(1./(hisOIS->GetMaximum()+1));
      hisOIS->SetLineColor(kmicolors[0]);
      hisOIS->Draw();
      isDraw = kTRUE;
      vecOIS(isec)=10*hisOIS->GetMean();
    }
    if (hisOOP) {
      hisOOP = (TH1F*)(hisOOP->Clone());
      hisOOP->SetDirectory(0);
      hisOOP->Scale(1./(hisOOP->GetMaximum()+1));
      hisOOP->SetLineColor(kmicolors[1]);      
      if (isDraw) hisOOP->Draw("same");
      if (!isDraw) {hisOOP->Draw(""); isDraw=kTRUE;}
      vecOOP(isec)=10*hisOOP->GetMean();
    }
    if (hisOOM) {
      hisOOM = (TH1F*)(hisOOM->Clone());
      hisOOM->SetDirectory(0);
      hisOOM->Scale(1/(hisOOM->GetMaximum()+1));
      hisOOM->SetLineColor(kmicolors[3]);
      if (isDraw) hisOOM->Draw("same");
      if (!isDraw) {hisOOM->Draw(""); isDraw=kTRUE;}
      vecOOM(isec)=10*hisOOM->GetMean();
    }
  }
  //
  //
  for (Int_t isec=0;isec<36;isec++){
    Bool_t isDraw=kFALSE;
    cIROCdy->cd(isec+1);
    Int_t secPlus = (isec%18==17)? isec-17:isec+1;
    Int_t secMinus= (isec%18==0) ? isec+17:isec-1;
    printf("%d\t%d\t%d\n",isec,secPlus,secMinus);
    TH1 * hisOIP= align->GetHisto(AliTPCcalibAlign::kY,isec+36,secPlus);
    TH1 * hisOIM= align->GetHisto(AliTPCcalibAlign::kY,isec+36,secMinus);    
    TH1 * hisOIS= align->GetHisto(AliTPCcalibAlign::kY,isec+36,isec);    
    if (hisOIS) {
      hisOIS = (TH1F*)(hisOIS->Clone());
      hisOIS->SetDirectory(0);
      hisOIS->Scale(1./(hisOIS->GetMaximum()+1));
      hisOIS->SetLineColor(kmicolors[0]);
      hisOIS->Draw();
      isDraw = kTRUE;
      vecOIS(isec)=10*hisOIS->GetMean();
    }
    if (hisOIP) {
      hisOIP = (TH1F*)(hisOIP->Clone());
      hisOIP->SetDirectory(0);
      hisOIP->Scale(1./(hisOIP->GetMaximum()+1));
      hisOIP->SetLineColor(kmicolors[1]);
      if (isDraw) hisOIP->Draw("same");
      if (!isDraw) {hisOIP->Draw(""); isDraw=kTRUE;}
      hisOIP->Draw("same");
      vecOIP(isec)=10*hisOIP->GetMean();
    }
    if (hisOIM) {
      hisOIM = (TH1F*)(hisOIM->Clone());
      hisOIM->SetDirectory(0);
      hisOIM->Scale(1/(hisOIM->GetMaximum()+1));
      hisOIM->SetLineColor(kmicolors[3]);
      if (isDraw) hisOIM->Draw("same");
      if (!isDraw) {hisOIM->Draw(""); isDraw=kTRUE;}
      vecOIM(isec)=10*hisOIM->GetMean();
    }
  }
  TGraph* grOIM = new TGraph(36,vecSec.GetMatrixArray(),vecOIM.GetMatrixArray());
  TGraph* grOIP = new TGraph(36,vecSec.GetMatrixArray(),vecOIP.GetMatrixArray());
  TGraph* grOIS = new TGraph(36,vecSec.GetMatrixArray(),vecOIS.GetMatrixArray());  
  TGraph* grOOM = new TGraph(36,vecSec.GetMatrixArray(),vecOOM.GetMatrixArray());
  TGraph* grOOP = new TGraph(36,vecSec.GetMatrixArray(),vecOOP.GetMatrixArray());
  //
  grOIS->SetMarkerStyle(kmimarkers[0]);
  grOIP->SetMarkerStyle(kmimarkers[1]);
  grOIM->SetMarkerStyle(kmimarkers[3]);
  grOOP->SetMarkerStyle(kmimarkers[1]);
  grOOM->SetMarkerStyle(kmimarkers[3]);
  grOIS->SetMarkerColor(kmicolors[0]);
  grOIP->SetMarkerColor(kmicolors[1]);
  grOIM->SetMarkerColor(kmicolors[3]);
  grOOP->SetMarkerColor(kmicolors[1]);
  grOOM->SetMarkerColor(kmicolors[3]);
  grOIS->SetLineColor(kmicolors[0]);
  grOIP->SetLineColor(kmicolors[1]);
  grOIM->SetLineColor(kmicolors[3]);
  grOOP->SetLineColor(kmicolors[1]);
  grOOM->SetLineColor(kmicolors[3]);
  grOIS->SetMaximum(1.5);
  grOIS->SetMinimum(-1.5);
  grOIS->GetXaxis()->SetTitle("Sector number");
  grOIS->GetYaxis()->SetTitle("#Delta_{y} (mm)");

  cDy->cd(1);
  grOIS->Draw("apl");
  grOIM->Draw("pl");
  grOIP->Draw("pl");
  cDy->cd(2);
  grOIS->Draw("apl");
  grOOM->Draw("pl");
  grOOP->Draw("pl");
  cOROCdy->SaveAs("picAlign/OROCOROCdy.eps");
  cOROCdy->SaveAs("picAlign/OROCOROCdy.gif");
  cOROCdy->SaveAs("picAlign/OROCOROCdy.root");
  //
  cIROCdy->SaveAs("picAlign/OROCIROCdy.eps");
  cIROCdy->SaveAs("picAlign/OROCIROCdy.gif");
  cIROCdy->SaveAs("picAlign/OROCIROCdy.root");
  //
  cDy->SaveAs("picAlign/Sectordy.eps");
  cDy->SaveAs("picAlign/Sectordy.gif");
  cDy->SaveAs("picAlign/Sectordy.root");
  if (output){
    output->cd();
    cOROCdy->Write("OROCOROCDy");
    cIROCdy->Write("OROCIROCDy");
    cDy->Write("SectorDy");
  }
}

void AliTPCcalibAlign::MakeReportDyPhi(TFile */*output*/){
  //
  //
  //
  Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};

  AliTPCcalibAlign *align = this;
  TCanvas * cOROCdyPhi = new TCanvas("OROC dyphi","OROC dyphi",900,600);
  cOROCdyPhi->Divide(6,6);
  for (Int_t isec=0;isec<36;isec++){
    cOROCdyPhi->cd(isec+1);
    Int_t secPlus = (isec%18==17)? isec-17:isec+1;
    Int_t secMinus= (isec%18==0) ? isec+17:isec-1;
    //printf("%d\t%d\t%d\n",isec,secPlus,secMinus);
    TH2F *htemp=0;
    TProfile * profdyphiOOP=0,*profdyphiOOM=0,*profdyphiOOS=0;
    htemp = (TH2F*) (align->GetHisto(AliTPCcalibAlign::kYPhi,isec+36,secPlus+36));
    if (htemp) profdyphiOOP= htemp->ProfileX();
    htemp = (TH2F*)(align->GetHisto(AliTPCcalibAlign::kYPhi,isec+36,secMinus+36));
    if (htemp) profdyphiOOM= htemp->ProfileX();
    htemp = (TH2F*)(align->GetHisto(AliTPCcalibAlign::kYPhi,isec+36,isec));
    if (htemp) profdyphiOOS= htemp->ProfileX();
    
    if (profdyphiOOS){
      profdyphiOOS->SetLineColor(kmicolors[0]);
      profdyphiOOS->SetMarkerStyle(kmimarkers[0]);
      profdyphiOOS->SetMarkerSize(0.2);
      profdyphiOOS->SetMaximum(0.5);
      profdyphiOOS->SetMinimum(-0.5);
      profdyphiOOS->SetXTitle("tan(#phi)");
      profdyphiOOS->SetYTitle("#DeltaY (cm)");
    }
    if (profdyphiOOP){
      profdyphiOOP->SetLineColor(kmicolors[1]);
      profdyphiOOP->SetMarkerStyle(kmimarkers[1]);
      profdyphiOOP->SetMarkerSize(0.2);
    }
    if (profdyphiOOM){
      profdyphiOOM->SetLineColor(kmicolors[3]);
      profdyphiOOM->SetMarkerStyle(kmimarkers[3]);
      profdyphiOOM->SetMarkerSize(0.2);
    }
    if (profdyphiOOS){
      profdyphiOOS->Draw();
    }else{
      if (profdyphiOOM) profdyphiOOM->Draw("");
      if (profdyphiOOP) profdyphiOOP->Draw("");
    }
    if (profdyphiOOM) profdyphiOOM->Draw("same");
    if (profdyphiOOP) profdyphiOOP->Draw("same");
    
  }
}

