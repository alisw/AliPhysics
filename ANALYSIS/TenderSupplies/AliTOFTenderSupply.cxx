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
// TOF tender: reapply TOF pid on the fly                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliESDpid.h>
#include <AliTender.h>

#include "AliTOFcalib.h"
#include "AliTOFT0maker.h"

#include "AliTOFTenderSupply.h"


AliTOFTenderSupply::AliTOFTenderSupply() :
  AliTenderSupply(),
  fESDpid(0x0),
  fTOFCalib(0x0),
  fTOFT0maker(0x0),
  fTOFres(100.),
  fIsMC(kFALSE)
{
  //
  // default ctor
  //
}

//_____________________________________________________
AliTOFTenderSupply::AliTOFTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fESDpid(0x0),
  fTOFCalib(0x0),
  fTOFT0maker(0x0),
  fTOFres(100.),
  fIsMC(kFALSE)
{
  //
  // named ctor
  //
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
      fTOFCalib->SetCorrectTExp(kTRUE); // apply a fine tuning on the expected times at low momenta
      if(fIsMC) fTOFCalib->SetCalibrateTOFsignal(kFALSE); // no new calibration
      fTOFCalib->Init();
  }
  if (!fTOFT0maker) {
      fTOFT0maker = new AliTOFT0maker(fESDpid,fTOFCalib);
      fTOFT0maker->SetTimeResolution(fTOFres); // set TOF resolution for the PID
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
    fTOFCalib->Init();
  }
  fTOFCalib->CalibrateESD(event);

  if(fIsMC) fTOFT0maker->TuneForMC(event);

  //Calculate event time zero
  fTOFT0maker->ComputeT0TOF(event);
  fTOFT0maker->ApplyT0TOF(event);

  event->SetT0(0.0);

  //
  // recalculate PID probabilities
  //

  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    fESDpid->MakeTOFPID(event->GetTrack(itrack),0);
  }

}


