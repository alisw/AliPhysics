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

/*
  laser track clasification;
  TCut cutT("cutT","abs(Tr.fP[3])<0.06");
  TCut cutPt("cutPt","abs(Tr.fP[4])<0.1");
  TCut cutN("cutN","fTPCncls>100");

  TCut cutFi("cutZB","");
  TCut cutFi("cutFi","abs((180*atan2(x1,x0)/pi-20)%40)<5");
  
*/



#include "TLinearFitter.h"
#include "AliTPCcalibLaser.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
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
#include "AliTPCLaserTrack.h"

using namespace std;

ClassImp(AliTPCcalibLaser)

AliTPCcalibLaser::AliTPCcalibLaser():
  AliTPCcalibBase()
{
  //
  // Constructor
  //
}

AliTPCcalibLaser::AliTPCcalibLaser(const Text_t *name, const Text_t *title):
  AliTPCcalibBase()  
{
  SetName(name);
  SetTitle(title);
  //
  // Constructor
  //
  
}

AliTPCcalibLaser::~AliTPCcalibLaser() {
  //
  // destructor
  //
}

void AliTPCcalibLaser::Process(AliESDtrack *track, Int_t run) {
  //
  // 
  //
  // 1. Propagate track to the mirror radius
  Float_t kRadius = 271.6;
  if (!track->GetOuterParam()) return;
  AliExternalTrackParam param(*(track->GetOuterParam()));
  AliTracker::PropagateTrackTo(&param,270.,0.0005,3,kTRUE);
  AliTracker::PropagateTrackTo(&param,kRadius,0.0005,0.1,kTRUE);
  AliTPCLaserTrack ltr;
  AliTPCLaserTrack *ltrp=0x0;
  Int_t id = AliTPCLaserTrack::IdentifyTrack(&param);
  if (id!=-1) ltrp=(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
  else ltrp=&ltr;
  //
  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    Double_t xyz[3];
    Double_t pxyz[3];
    param.GetXYZ(xyz);
    param.GetPxPyPz(pxyz);
    Int_t side = (param.GetZ()>0) ? 0:1;
    Int_t beam = 0;
    if (side==0) beam = TMath::Nint((180*atan2(xyz[1],xyz[0])/TMath::Pi()+20)/60.);
    if (side==1) beam = TMath::Nint((180*atan2(xyz[1],xyz[0])/TMath::Pi()-20)/60.);
    //Int_t id(180*atan2(x1,x0)/pi+20)/60.;
    Int_t bundle=TMath::Nint(param.GetZ()/80.);
    if (cstream){
      (*cstream)<<"Track"<<
	"run="<<run<<
	"id="<<id<<
	"fSide="<<side<<    // side A-C
	"fBeam="<<beam<<   // phi id
	"fBundle="<<bundle<< // laser Z
	//
        "LTr.="<<ltrp<<
	"Esd.="<<track<<
	"Tr.="<<&param<<
	"x0="<<xyz[0]<<
	"x1="<<xyz[1]<<
	"x2="<<xyz[2]<<
	"px0="<<pxyz[0]<<
	"px1="<<pxyz[1]<<
	"px2="<<pxyz[2]<<
	"\n";
    }
  }
}

void AliTPCcalibLaser::Analyze(){
  //
  //
  //
}




