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

/* $Id$ */
//-------------------------------------------
// Class : AliVZEROTriggerMask
//
// Fill up the trigger mask word.
//

#include <Riostream.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TF1.h>
#include <TMath.h>

#include <AliGeomManager.h>
#include "AliLog.h"
#include "AliVZEROTriggerMask.h"
#include "AliVZEROConst.h"
#include "AliVZEROCalibData.h"
#include "AliESDVZERO.h"
#include "AliVZEROReconstructor.h"

//______________________________________________________________________
ClassImp(AliVZEROTriggerMask)

//______________________________________________________________________

AliVZEROTriggerMask::AliVZEROTriggerMask()
  :TObject(),
   fAdcThresHold(0.0),
   fTimeWindowBBALow(-9.5),
   fTimeWindowBBAUp(22.5),
   fTimeWindowBGALow(-2.5),
   fTimeWindowBGAUp(5.0),
   fTimeWindowFakeALow(-17.5),
   fTimeWindowFakeAUp(-9.5),
   fTimeWindowBBCLow(-2.5),
   fTimeWindowBBCUp(22.5),
   fTimeWindowBGCLow(-2.5),
   fTimeWindowBGCUp(2.5),
   fTimeWindowFakeCLow(-22.5),
   fTimeWindowFakeCUp(-8.5),
   fV0ADist(0),
   fV0CDist(0)
{
  // Default constructor
  //
  Float_t zV0A = TMath::Abs(GetZPosition("VZERO/V0A"));
  Float_t zV0C = TMath::Abs(GetZPosition("VZERO/V0C"));

  // distance in time units from nominal vertex to V0
  fV0ADist = zV0A/TMath::Ccgs()*1e9;
  fV0CDist = zV0C/TMath::Ccgs()*1e9;
}

//________________________________________________________________________________
Double_t AliVZEROTriggerMask::GetZPosition(const char* symname){
// Get the global z coordinate of the given V0 alignable volume
//
  Double_t *tr;
  TGeoPNEntry *pne = gGeoManager->GetAlignableEntry(symname);
  if (!pne) {
    AliFatalClass(Form("TGeoPNEntry with symbolic name %s does not exist!",symname));
    return 0;
  }

  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if(pnode){
          TGeoHMatrix* hm = pnode->GetMatrix();
           tr = hm->GetTranslation();
  }else{
          const char* path = pne->GetTitle();
          if(!gGeoManager->cd(path)){
                  AliFatalClass(Form("Volume path %s not valid!",path));
                  return 0;
          }
         tr = gGeoManager->GetCurrentMatrix()->GetTranslation();
  }
  return tr[2];

}

//________________________________________________________________________________


void AliVZEROTriggerMask::FillMasks(AliESDVZERO *esdV0,
				    AliVZEROCalibData *cal,
				    TF1 *slewing)
{
  // Fill up the trigger mask word
  // using the TDC data (already corrected for
  // slewing and misalignment between channels)

  esdV0->SetBit(AliESDVZERO::kTriggerBitsFilled,kTRUE);
  esdV0->SetBit(AliESDVZERO::kDecisionFilled,kTRUE);

  UInt_t aBBtriggerV0A = 0; // bit mask for Beam-Beam trigger in V0A
  UInt_t aBGtriggerV0A = 0; // bit mask for Beam-Gas trigger in V0A
  UInt_t aBBtriggerV0C = 0; // bit mask for Beam-Beam trigger in V0C
  UInt_t aBGtriggerV0C = 0; // bit mask for Beam-Gas trigger in V0C

  const Float_t p1 = 2.50; // photostatistics term in the time resolution
  const Float_t p2 = 3.00; // slewing related term in the time resolution

  // loop over vzero channels
  Float_t timeAW = 0,timeCW = 0;
  Float_t weightA = 0,weightC = 0;
  Int_t ntimeA = 0, ntimeC = 0;
  for (Int_t i = 0; i < 64; ++i) {
    Float_t adc = esdV0->GetAdc(i);
    if (adc > fAdcThresHold) {
      Float_t tdc = esdV0->GetTime(i);
      if (tdc > (AliVZEROReconstructor::kInvalidTime + 1e-6)) {
	Float_t nphe = adc*kChargePerADC/(cal->GetGain(i)*TMath::Qe());
	Float_t timeErr = TMath::Sqrt(kIntTimeRes*kIntTimeRes+
				      p1*p1/nphe+
				      p2*p2*(slewing->GetParameter(0)*slewing->GetParameter(1))*(slewing->GetParameter(0)*slewing->GetParameter(1))*
				      TMath::Power(adc/cal->GetDiscriThr(i),2.*(slewing->GetParameter(1)-1.))/
				      (cal->GetDiscriThr(i)*cal->GetDiscriThr(i)));

	if (i < 32) { // in V0C
	  ntimeC++;
	  timeCW += tdc/(timeErr*timeErr);
	  weightC += 1.0/(timeErr*timeErr);

	  if (tdc > (fV0CDist + fTimeWindowBBCLow) &&
	      tdc < (fV0CDist + fTimeWindowBBCUp))
	    aBBtriggerV0C |= (1 << i);
	  if (tdc > (-fV0CDist + fTimeWindowBGCLow) &&
	      tdc < (-fV0CDist + fTimeWindowBGCUp))
	    aBGtriggerV0C |= (1 << i); 
	}
	else { // in V0A
	  ntimeA++;
	  timeAW += tdc/(timeErr*timeErr);
	  weightA += 1.0/(timeErr*timeErr);

	  Int_t shift = i - 32;
	  if (tdc > (fV0ADist + fTimeWindowBBALow) &&
	      tdc < (fV0ADist + fTimeWindowBBAUp)) 
	    aBBtriggerV0A |= (1 << shift);
	  if (tdc > (-fV0ADist + fTimeWindowBGALow) &&
	      tdc < (-fV0ADist + fTimeWindowBGAUp))
	    aBGtriggerV0A |= (1 << shift); 
	}
      }
    }
  } // end of loop over channels

  esdV0->SetBBtriggerV0A(aBBtriggerV0A);
  esdV0->SetBGtriggerV0A(aBGtriggerV0A);
  esdV0->SetBBtriggerV0C(aBBtriggerV0C);
  esdV0->SetBGtriggerV0C(aBGtriggerV0C);

  if (weightA > 0) timeAW = timeAW/weightA;
  else timeAW = AliVZEROReconstructor::kInvalidTime;

  if (weightC > 0) timeCW = timeCW/weightC;
  else timeCW = AliVZEROReconstructor::kInvalidTime;

  esdV0->SetV0ATime(timeAW);
  esdV0->SetV0CTime(timeCW);
  esdV0->SetV0ATimeError((weightA > 0) ? (1./TMath::Sqrt(weightA)) : 0);
  esdV0->SetV0CTimeError((weightC > 0) ? (1./TMath::Sqrt(weightC)) : 0);

  esdV0->SetV0ADecision(AliESDVZERO::kV0Empty);
  esdV0->SetV0CDecision(AliESDVZERO::kV0Empty);

  if (timeAW > (fV0ADist + fTimeWindowBBALow) &&
      timeAW < (fV0ADist + fTimeWindowBBAUp)) 
    esdV0->SetV0ADecision(AliESDVZERO::kV0BB);
  else if (timeAW > (-fV0ADist + fTimeWindowBGALow) &&
	   timeAW < (-fV0ADist + fTimeWindowBGAUp))
    esdV0->SetV0ADecision(AliESDVZERO::kV0BG);
  else if (timeAW > (fV0ADist + fTimeWindowFakeALow) &&
	   timeAW < (fV0ADist + fTimeWindowFakeAUp))
    esdV0->SetV0ADecision(AliESDVZERO::kV0Fake);

  if (timeCW > (fV0CDist + fTimeWindowBBCLow) &&
      timeCW < (fV0CDist + fTimeWindowBBCUp)) 
    esdV0->SetV0CDecision(AliESDVZERO::kV0BB);
  else if (timeCW > (-fV0CDist + fTimeWindowBGCLow) &&
	   timeCW < (-fV0CDist + fTimeWindowBGCUp))
    esdV0->SetV0CDecision(AliESDVZERO::kV0BG);
  else if (timeCW > (fV0CDist + fTimeWindowFakeCLow) &&
	   timeCW < (fV0CDist + fTimeWindowFakeCUp))
    esdV0->SetV0CDecision(AliESDVZERO::kV0Fake);

}
