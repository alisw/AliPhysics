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
// TOF tender: 
//             load updated calibration if needed
//             reapply TOF pid on the fly                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliESDpid.h>
#include <AliTender.h>

#include <AliTOFcalib.h>
#include <AliTOFT0maker.h>

#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliT0CalibSeasonTimeShift.h>

#include "AliTOFTenderSupply.h"

AliTOFTenderSupply::AliTOFTenderSupply() :
  AliTenderSupply(),
  fESDpid(0x0),
  fIsMC(kFALSE),
  fApplyT0(kFALSE),
  fTimeZeroType(3),
  fCorrectExpTimes(kTRUE),
  fLHC10dPatch(kFALSE),
  fTOFCalib(0x0),
  fTOFT0maker(0x0),
  fTOFres(100.)


{
  //
  // default ctor
  //

  fT0shift[0] = 0;
  fT0shift[1] = 0;
  fT0shift[2] = 0;
  fT0shift[3] = 0;
}

//_____________________________________________________
AliTOFTenderSupply::AliTOFTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fESDpid(0x0),
  fIsMC(kFALSE),
  fApplyT0(kFALSE),
  fTimeZeroType(3),
  fCorrectExpTimes(kTRUE),
  fLHC10dPatch(kFALSE),
  fTOFCalib(0x0),
  fTOFT0maker(0x0),
  fTOFres(100.) 
 
{
  //
  // named ctor
  //

  fT0shift[0] = 0;
  fT0shift[1] = 0;
  fT0shift[2] = 0;
  fT0shift[3] = 0;
}

//_____________________________________________________
void AliTOFTenderSupply::Init()
{
  //
  // Initialise TOF tender
  //


  //
  // Setup PID object
  //

  // Check if another detector already created the esd pid object
  // if not we create it and set it to the ESD input handler
  fESDpid=fTender->GetESDhandler()->GetESDpid();
  if (!fESDpid) {
    fESDpid=new AliESDpid;
    fTender->GetESDhandler()->SetESDpid(fESDpid);
  }

  //Set proper resolution in case of MC
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  if (mgr->GetMCtruthEventHandler()) fIsMC=kTRUE;


  //
  // Create TOF calibration classes
  //
  if (!fTOFCalib){
      fTOFCalib=new AliTOFcalib();
      fTOFCalib->SetCorrectTExp(fCorrectExpTimes); // apply a fine tuning on the expected times at low momenta
      if(fIsMC) fTOFCalib->SetCalibrateTOFsignal(kFALSE); // no new calibration
//      fTOFCalib->Init();
  }
  if (!fTOFT0maker) {
      fTOFT0maker = new AliTOFT0maker(fESDpid,fTOFCalib);
      fTOFT0maker->SetTimeResolution(fTOFres); // set TOF resolution for the PID
      printf("tof time res = %f\n",fTOFres);
  }
}

//_____________________________________________________
void AliTOFTenderSupply::ProcessEvent()
{
  //
  // Reapply pid information
  //

  //no corrections for MC
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  if (mgr->GetMCtruthEventHandler()) return;

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
    
  //recalculate TOF signal
  if (fTender->RunChanged()){
    fTOFCalib->Init(fTender->GetRun());
    
    if(event->GetT0TOF()){ // read T0 detector correction from OCDB
      // OCDB instance
      AliCDBManager* ocdbMan = AliCDBManager::Instance();
      ocdbMan->SetRun(fTender->GetRun());    
      AliCDBEntry *entry = ocdbMan->Get("T0/Calib/TimeAdjust/");
      if(entry) {
	AliT0CalibSeasonTimeShift *clb = (AliT0CalibSeasonTimeShift*) entry->GetObject();
	Float_t *t0means= clb->GetT0Means();
	//      Float_t *t0sigmas = clb->GetT0Sigmas();
	fT0shift[0] = t0means[0];
	fT0shift[1] = t0means[1];
	fT0shift[2] = t0means[2];
	fT0shift[3] = t0means[3];
      }   
    }
  }
  fTOFCalib->CalibrateESD(event);

  if (fLHC10dPatch) RecomputeTExp(event);

  if(fIsMC){
    Float_t t0true = fTOFT0maker->TuneForMC(event);
    
    if(event->GetT0TOF()){// add the t0 smearing also to the T0 detector information
      event->SetT0TOF(0,event->GetT0TOF(0) + t0true);
      event->SetT0TOF(1,event->GetT0TOF(1) + t0true);
      event->SetT0TOF(2,event->GetT0TOF(2) + t0true);  
    }
  }
  
  if(event->GetT0TOF()){ // write the T0 detector corrected times
    if(event->GetT0TOF(0) == 0) event->SetT0TOF(0, 9999999.);
    if(event->GetT0TOF(1) == 0) event->SetT0TOF(1, 99999.);
    if(event->GetT0TOF(2) == 0) event->SetT0TOF(2, 99999.);

    event->SetT0TOF(0,event->GetT0TOF(0) - fT0shift[0]);
    event->SetT0TOF(1,event->GetT0TOF(1) - fT0shift[1]);
    event->SetT0TOF(2,event->GetT0TOF(2) - fT0shift[2]);      

    if(event->GetT0TOF(0) > 9000000) event->SetT0TOF(0, 9999999.);
    if(event->GetT0TOF(1) > 90000) event->SetT0TOF(1, 99999.);
    if(event->GetT0TOF(2) > 90000) event->SetT0TOF(2, 99999.);
  }
  
  //Calculate event time zero
  fTOFT0maker->ComputeT0TOF(event);
  fTOFT0maker->WriteInESD(event);

  // subtract the T0-TOF information to the TOF times
  if(fApplyT0){
    fTOFT0maker->ApplyT0TOF(event);
    event->SetT0(0.0);
  }

  //
  // recalculate PID probabilities
  //

  fESDpid->SetTOFResponse(event, (AliESDpid::EStartTimeType_t)fTimeZeroType);
  
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    fESDpid->MakeTOFPID(event->GetTrack(itrack),0);
  }
  
  
}


//_____________________________________________________
void AliTOFTenderSupply::RecomputeTExp(AliESDEvent *event) const
{
  /*
   * calibrate TExp
   */

  
  /* loop over tracks */
  AliESDtrack *track = NULL;
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
    /* get track and calibrate */
    track = event->GetTrack(itrk);
    RecomputeTExp(track);
  }
  
}

//_____________________________________________________
void AliTOFTenderSupply::RecomputeTExp(AliESDtrack *track) const
{
  /*** 
       THIS METHOD IS BASED ON THEORETICAL EXPECTED TIME COMPUTED
       USING AVERAGE MOMENTUM BETWEEN INNER/OUTER TRACK PARAMS 
       IT IS A ROUGH APPROXIMATION APPLIED TO FIX LHC10d-pass2 DATA
       WHERE A WRONG GEOMETRY (FULL TRD) WAS INSERTED
  ***/

  Double_t texp[AliPID::kSPECIES];
  if (!track || !(track->GetStatus() & AliESDtrack::kTOFout)) return;


  /* get track params */
  Float_t l = track->GetIntegratedLength();
  Float_t p = track->P();
  if (track->GetInnerParam() && track->GetOuterParam()) {
    Float_t pin = track->GetInnerParam()->P();
    Float_t pout = track->GetOuterParam()->P();
    p = 0.5 * (pin + pout);
  }
  /* loop over particle types and compute expected time */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    texp[ipart] = GetExpTimeTh(AliPID::ParticleMass(ipart), p, l) - 37.; 
  // 37 is a final semiempirical offset to further adjust (calibrations were
  // done with "standard" integratedTimes)
  /* set integrated times */
  track->SetIntegratedTimes(texp);

}
