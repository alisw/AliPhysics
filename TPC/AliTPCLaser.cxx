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

#include <stdlib.h>

#include <TLorentzVector.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliTPCDigitsArray.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCTrackHitsV2.h"
#include "AliTPCLaser.h"

ClassImp(AliTPCLaser)
 
//_____________________________________________________________________________
AliTPCLaser::AliTPCLaser(const char *name, const char *title) :
  AliTPCv2(name, title) 
{
  // only use the AliTPCv2 constructor
}
//______________________________________________________________
void AliTPCLaser::StepManager()
{
  // laser tracks are muons (PID=13) 
  // stopped in the the inner containment vessel (PID=14) 

   TVirtualMC* mc = TVirtualMC::GetMC();
   Int_t copy, vol;
   vol = mc->CurrentVolID(copy);
   // Debug
   // printf("Vol name %s\n",mc->CurrentVolName());
   if (mc->TrackPid() == 13 // muons
       && vol == 14) {// 14 = TIIN (inner containment vessel)
     mc->StopTrack();
     return;
   }
   AliTPCv2::StepManager();   
}
