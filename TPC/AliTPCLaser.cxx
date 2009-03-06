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

//
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Laser for the TPChamber version 2 -- detailed TPC and slow simulation    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TLorentzVector.h>
#include "AliMC.h"


#include "AliTPCLaser.h"


ClassImp(AliTPCLaser)
 
//_____________________________________________________________________________
AliTPCLaser::AliTPCLaser(const char *name, const char *title) :
  AliTPCv2(name, title),
  fNelPerCollision(10),
  fLaserPID(13), // muons
  fCollisionsPerCm(20)
{

}
//______________________________________________________________
void AliTPCLaser::StepManager()
{
  // laser tracks are particles with PID fLaserPID (default PID=13) 
  // stopped in the the TPC inner containment vessel (14) 

  if (gMC->TrackPid() != fLaserPID) {
    // in this way we can prevent delta-electrons
    gMC->StopTrack();
    return;
  }
  
  Int_t copy;
  Int_t vol[2];
  vol[0] = gMC->CurrentVolID(copy);
  
  if (gMC->TrackPid() == fLaserPID
      && vol[0] == 14) {// 14 = TIIN (inner containment vessel)
    gMC->StopTrack();
    return;
  }
  
  TLorentzVector p;
  Float_t hits[5]={0,0,0,0,0};
  gMC->TrackPosition(p);
  hits[0]=p[0];
  hits[1]=p[1];
  hits[2]=p[2];
  hits[3]=fNelPerCollision;
  hits[4]=gMC->TrackTime();

  Int_t index[3];  
  vol[0]=fTPCParam->Transform0to1(hits,index);
  AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits);
  
  Double_t rnd = gMC->GetRandom()->Rndm();  
  gMC->SetMaxStep(-TMath::Log(rnd)/fCollisionsPerCm);
}
