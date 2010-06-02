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

////////////////////////////////////////////////////////////////////
//                                                                //
// recalculate the TOF signal from the TOF raw signal             //
// using the updates in the OCDB                                  //
//                                                                //
////////////////////////////////////////////////////////////////////

#include "AliTOFChannelOnlineStatusArray.h"
#include "TObjArray.h"
#include "AliTOFDeltaBCOffset.h"
#include "AliTOFCTPLatency.h"
#include "AliTOFRunParams.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTOFChannelOffline.h"
#include "AliTOFGeometry.h"

#include "AliTOFcalibESD.h"

ClassImp(AliTOFcalibESD)

//______________________________________________________________-

AliTOFcalibESD::AliTOFcalibESD() :
  TObject(),
  fInitFlag(kFALSE),
  fChannelStatusArray(NULL),
  fParOfflineArray(NULL),
  fDeltaBCOffsetObj(NULL),
  fCTPLatencyObj(NULL),
  fRunParamsObj(NULL),
  fTimeZero(0),
  fTOFResolution(130)
{
  /*
   * default constructor
   */
}

//______________________________________________________________-

AliTOFcalibESD::~AliTOFcalibESD()
{
  /*
   * default destructor
   */
}

//______________________________________________________________-

Bool_t
AliTOFcalibESD::Init(Int_t run)
{
  /*
   * init
   */

  /* get cdb instance */
  AliCDBManager *cdb = AliCDBManager::Instance();
  AliCDBEntry *entry = NULL;

  /* get channel status array */
  entry = cdb->Get("TOF/Calib/Status", run);
  if (!entry || !entry->GetObject()) return kFALSE;
  fChannelStatusArray = (AliTOFChannelOnlineStatusArray *)entry->GetObject();
  /* get par offline array */
  entry = cdb->Get("TOF/Calib/ParOffline", run);
  if (!entry || !entry->GetObject()) return kFALSE;
  fParOfflineArray = (TObjArray *)entry->GetObject();
  /* get deltaBC offset obj */
  entry = cdb->Get("TOF/Calib/DeltaBCOffset", run);
  if (!entry || !entry->GetObject()) return kFALSE;
  fDeltaBCOffsetObj = (AliTOFDeltaBCOffset *)entry->GetObject();
  /* get CTP latency obj */
  entry = cdb->Get("TOF/Calib/CTPLatency", run);
  if (!entry || !entry->GetObject()) return kFALSE;
  fCTPLatencyObj = (AliTOFCTPLatency *)entry->GetObject();
  /* get run params obj */
  entry = cdb->Get("TOF/Calib/RunParams", run);
  if (!entry || !entry->GetObject()) return kFALSE;
  fRunParamsObj = (AliTOFRunParams *)entry->GetObject();

  /* all done */
  fInitFlag = kTRUE;
  return kTRUE;
}

//______________________________________________________________-

void
AliTOFcalibESD::CalibrateESD(AliESDEvent *event)
{
  /*
   * calibrate ESD
   */
  
  /* get global calibration params */
  AliTOFChannelOffline *parOffline = NULL;
//   Int_t deltaBCOffset = fDeltaBCOffsetObj->GetDeltaBCOffset();  //see below
  Float_t ctpLatency = fCTPLatencyObj->GetCTPLatency();
  Float_t tdcLatencyWindow;
  Float_t timezero = fRunParamsObj->EvalT0(event->GetTimeStamp());
  fTimeZero=timezero;
  fTOFResolution=fRunParamsObj->EvalTOFResolution(event->GetTimeStamp());
  
  /* loop over tracks */
  AliESDtrack *track = NULL;
  Int_t index, l0l1, deltaBC;
  Double_t time, tot;
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {

    /* get track */
    track = event->GetTrack(itrk);
    if (!track || !(track->GetStatus() & AliESDtrack::kTOFout)) continue;
    
    /* get info */
    index = track->GetTOFCalChannel();
    time = track->GetTOFsignalRaw();
    tot = track->GetTOFsignalToT();
    l0l1 = track->GetTOFL0L1();
    deltaBC = track->GetTOFDeltaBC();
    
    /* get channel dependent calibration params */
    parOffline = (AliTOFChannelOffline *)fParOfflineArray->At(index);
    tdcLatencyWindow = fChannelStatusArray->GetLatencyWindow(index) * 1.e3;
    
    /* deltaBC correction (inhibited for the time being) */
    //  time -= (deltaBC - deltaBCOffset) * AliTOFGeometry::BunchCrossingBinWidth();
    /* L0-L1 latency correction */
    time += l0l1 * AliTOFGeometry::BunchCrossingBinWidth();
    /* CTP latency correction */
    time += ctpLatency;
    /* TDC latency window correction */
    time -= tdcLatencyWindow;
    /* time-zero correction */
    time -= timezero;
    /* time calibration correction */
    if (tot < AliTOFGeometry::SlewTOTMin()) 
      tot = AliTOFGeometry::SlewTOTMin();
    if (tot > AliTOFGeometry::SlewTOTMax()) 
      tot = AliTOFGeometry::SlewTOTMax();
    for (Int_t islew = 0; islew < 6; islew++)
      time -= parOffline->GetSlewPar(islew) * TMath::Power(tot, islew) * 1.e3;
    
    /* set new TOF signal */
    track->SetTOFsignal(time);

  }

}
