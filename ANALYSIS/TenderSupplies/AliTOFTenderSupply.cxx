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

#include "AliTOFcalibESD.h"
#include "AliTOFT0makerANA.h"

#include "AliTOFTenderSupply.h"


AliTOFTenderSupply::AliTOFTenderSupply() :
  AliTenderSupply(),
  fESDpid(0x0),
  fTOFesdCalib(0x0),
  fTOFT0maker(0x0),
  fTOFres(130.)
{
  //
  // default ctor
  //
}

//_____________________________________________________
AliTOFTenderSupply::AliTOFTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fESDpid(0x0),
  fTOFesdCalib(0x0),
  fTOFT0maker(0x0),
  fTOFres(130.)
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
  if (mgr->GetMCtruthEventHandler())   fESDpid->GetTOFResponse().SetTimeResolution(80.);
  
  
  //
  // Create TOF calibration classes
  //
  if (!fTOFesdCalib) fTOFesdCalib=new AliTOFcalibESD;
  if (!fTOFT0maker) {
    fTOFT0maker = new AliTOFT0makerANA(fESDpid);
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
    fTOFesdCalib->Init(fTender->GetRun());
  }
  fTOFesdCalib->CalibrateESD(event);
  
  //Calculate event time zero
  Double_t* calcolot0;
  calcolot0=fTOFT0maker->RemakePID(event); // calculate T0-TOF(T0-FILL) and
  Double_t t0best=calcolot0[0];  // T0-Event = (T0-TOF .OR. T0-FILL) <- This is what you asked me
  event->SetT0(t0best);
  
  //
  // recalculate PID probabilities
  //
  
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    fESDpid->MakeTOFPID(event->GetTrack(itrack),0);
  }
  
}
