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
  TCut cutFi("cutFi","abs((180*atan2(x1,x0)/pi-20)%60)<5");
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

void AliTPCcalibLaser::Process(AliESDtrack *track) {
  //
  // 
  //
  // 1. Propagate track to the mirror radius
  Float_t kRadius = 271.6;
  if (!track->GetOuterParam()) return;
  AliExternalTrackParam param(*(track->GetOuterParam()));
  AliTracker::PropagateTrackTo(&param,270.,0.0005,3,kTRUE);
  AliTracker::PropagateTrackTo(&param,kRadius,0.0005,0.1,kTRUE);
  //
  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    Double_t xyz[3];
    Double_t pxyz[3];
    param.GetXYZ(xyz);
    param.GetPxPyPz(pxyz);
    if (cstream){
      (*cstream)<<"Track"<<
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


void AliTPCcalibLaser::Terminate(){
  //
  //
  //
}




