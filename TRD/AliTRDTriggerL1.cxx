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

///////////////////////////////////////////////////////
//                                                   //
//                                                   //
//  TRD main trigger class for L1                    //
//                                                   //
//                                                   //
///////////////////////////////////////////////////////

#include <TMath.h>

#include "AliRun.h"
#include "AliTriggerInput.h"

#include "AliTRDTriggerL1.h"
#include "AliTRDtrigParam.h"
#include "AliTRDtrigger.h"
#include "AliTRDgtuTrack.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDTriggerL1)

//_____________________________________________________________________________
AliTRDTriggerL1::AliTRDTriggerL1():AliTriggerDetector()
{

  SetName("TRD");

}

//_____________________________________________________________________________
void AliTRDTriggerL1::CreateInputs()
{

  fInputs.AddLast(new AliTriggerInput( "TRD_HadrLPt_L1",          "Single hadron low pt ",        0x01 ));
  fInputs.AddLast(new AliTriggerInput( "TRD_HadrHPt_L1",          "Single hadron high pt",        0x02 ));
  fInputs.AddLast(new AliTriggerInput( "TRD_Unlike_EPair_L1",     "Unlike electron pair",         0x04 ));
  fInputs.AddLast(new AliTriggerInput( "TRD_Unlike_EPair_HPt_L1", "Unlike electron pair high pt", 0x08 ));
  fInputs.AddLast(new AliTriggerInput( "TRD_Like_EPair_L1",       "Like electron pair",           0x10 ));
  fInputs.AddLast(new AliTriggerInput( "TRD_Like_EPair_HPt_L1",   "Like electron pair high pt",   0x20 ));
  fInputs.AddLast(new AliTriggerInput( "TRD_Electron_L1",         "Single electron",              0x40 ));
  fInputs.AddLast(new AliTriggerInput( "TRD_Electron_HPt_L1",     "Single electron high pt",      0x80 ));

}

//_____________________________________________________________________________
void AliTRDTriggerL1::Trigger()
{

  AliRunLoader* runLoader = gAlice->GetRunLoader();

  AliLoader *loader=runLoader->GetLoader("TRDLoader");

  Int_t nEvents = runLoader->GetNumberOfEvents();

  // Trigger (tracklets, LTU)

  AliTRDgeometry *geo = AliTRDgeometry::GetGeometry(runLoader);

  AliTRDtrigger trdTrigger("Trigger","Trigger class"); 

  AliTRDtrigParam *trigp = new AliTRDtrigParam("TRDtrigParam","TRD Trigger parameters");

  gAlice = runLoader->GetAliRun();
  Double_t x[3] = { 0.0, 0.0, 0.0 };
  Double_t b[3];
  gAlice->Field(x,b);  // b[] is in kilo Gauss
  Float_t field = b[2] * 0.1; // Tesla
  Info("Trigger","Trigger set for magnetic field = %f Tesla \n",field);

  trigp->SetField(field);
  trigp->Init();

  trdTrigger.SetParameter(trigp);
  trdTrigger.SetRunLoader(runLoader);
  trdTrigger.Init();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {

    trdTrigger.Open(runLoader->GetFileName(), iEvent);
    trdTrigger.ReadDigits();
    trdTrigger.MakeTracklets(kTRUE);
    trdTrigger.WriteTracklets(-1);

  }

  // Trigger (tracks, GTU)

  Float_t highPt = trigp->GetHighPt();

  Float_t pid, pt;
  Int_t   det, sec;
  Bool_t  isElectron;

  const Int_t maxEle = 1000;

  Int_t   electronPlus,          electronMinus;
  Int_t   sectorElePlus[maxEle], sectorEleMinus[maxEle];
  Float_t ptElePlus[maxEle],     ptEleMinus[maxEle];
  Int_t   hadronLowPt, hadronHighPt;

  hadronLowPt  = 0;
  hadronHighPt = 0;

  electronPlus  = 0;
  electronMinus = 0;

  AliTRDgtuTrack *gtuTrack;
  Int_t nTracks = trdTrigger.GetNumberOfTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

    gtuTrack = trdTrigger.GetTrack(iTrack);

    pid        = gtuTrack->GetPID();
    isElectron = gtuTrack->IsElectron();
    pt         = gtuTrack->GetPt();
    det        = gtuTrack->GetDetector();

    sec        = geo->GetSector(det);

    if (isElectron) {

      if (pt < 0.0) {
	sectorEleMinus[electronMinus] = sec;
	ptEleMinus[electronMinus]     = pt;
	electronMinus++;
      } else {
	sectorElePlus[electronPlus] = sec;
	ptElePlus[electronPlus]     = pt;
	electronPlus++;
      }

    } else {

      if (TMath::Abs(pt) < highPt) {
	hadronLowPt++;
      } else {
	hadronHighPt++;
      }

    }

  }

  loader->UnloadTracks();

  // hadrons

  if (hadronLowPt)  SetInput("TRD_Hadr_LPt_L1");
  if (hadronHighPt) SetInput("TRD_Hadr_HPt_L1");

  // electron-positron pairs (open angle > 80 deg)

  Int_t  secPlus, secMinus, secDiff;
  Bool_t electronUnlikePair    = kFALSE;
  Bool_t electronUnlikePairHPt = kFALSE;

  if (electronMinus > 0 && electronPlus > 0) {
    for (Int_t iPlus = 0; iPlus < electronPlus; iPlus++) {
      secPlus = sectorElePlus[iPlus];
      for (Int_t iMinus = 0; iMinus < electronMinus; iMinus++) {
	secMinus = sectorEleMinus[iMinus];
	secDiff = TMath::Abs(secPlus-secMinus);
	if (secDiff >  9) secDiff = 18 - secDiff;
	if (secDiff >= 5) {
	  electronUnlikePair = kTRUE;
	  if (TMath::Abs(ptElePlus[iPlus]) > highPt && TMath::Abs(ptEleMinus[iMinus]) > highPt) {
	    electronUnlikePairHPt = kTRUE;
	  }
	}
      }
    }
  }

  if (electronUnlikePair)    SetInput("TRD_Unlike_EPair_L1");
  if (electronUnlikePairHPt) SetInput("TRD_Unlike_EPair_HPt_L1");

  // like electron/positron pairs

  Bool_t ele1, ele1HPt;
  Bool_t ele2, ele2HPt;

  // positive

  ele1    = kFALSE;
  ele2    = kFALSE;
  ele1HPt = kFALSE;
  ele2HPt = kFALSE;
  if (electronPlus > 1) {
    for (Int_t iPlus = 0; iPlus < electronPlus; iPlus++) {
      if (!ele1) {
	ele1 = kTRUE;
      } else if (!ele2) {
	ele2 = kTRUE;
      }
      if (TMath::Abs(ptElePlus[iPlus]) > highPt) {
	if (!ele1HPt) {
	  ele1HPt = kTRUE;
	} else if (!ele2HPt) {
	  ele2HPt = kTRUE;
	}
      }
    }
  }

  if (ele1    && ele2   ) SetInput("TRD_Like_EPair_L1");
  if (ele1HPt && ele2HPt) SetInput("TRD_Like_EPair_HPt_L1");
  
  // negative

  ele1    = kFALSE;
  ele2    = kFALSE;
  ele1HPt = kFALSE;
  ele2HPt = kFALSE;
  if (electronMinus > 1) {
    for (Int_t iMinus = 0; iMinus < electronMinus; iMinus++) {
      if (!ele1) {
	ele1 = kTRUE;
      } else if (!ele2) {
	ele2 = kTRUE;
      }
      if (TMath::Abs(ptEleMinus[iMinus]) > highPt) {
	if (!ele1HPt) {
	  ele1HPt = kTRUE;
	} else if (!ele2HPt) {
	  ele2HPt = kTRUE;
	}
      }
    }
  }

  if (ele1    && ele2   ) SetInput("TRD_Like_EPair_L1");
  if (ele1HPt && ele2HPt) SetInput("TRD_Like_EPair_HPt_L1");
  
  // single electron/positron

  if (electronPlus > 0 || electronMinus > 0) {
    SetInput("TRD_Electron_L1");
    for (Int_t iPlus = 0; iPlus < electronPlus; iPlus++) {
      if (TMath::Abs(ptElePlus[iPlus]) > highPt) SetInput("TRD_Electron_HPt_L1");
      break;
    }
    for (Int_t iMinus = 0; iMinus < electronMinus; iMinus++) {
      if (TMath::Abs(ptEleMinus[iMinus]) > highPt) SetInput("TRD_Electron_HPt_L1");
      break;
    }
  }

}

