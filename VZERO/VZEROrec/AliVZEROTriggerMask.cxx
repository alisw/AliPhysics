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
   fV0ADist(0),
   fV0CDist(0),
   fRecoParam(NULL)
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

  const Float_t p1 = 2.50; // photostatistics term in the time resolution
  const Float_t p2 = 3.00; // slewing related term in the time resolution

  // loop over vzero channels
  Float_t timeAW = 0,timeCW = 0;
  Float_t weightA = 0,weightC = 0;
  Int_t ntimeA = 0, ntimeC = 0;
  Double_t timesA[32], timesC[32];
  Double_t wA[32],wC[32];
  Int_t indA[32], indC[32];
  for (Int_t i = 0; i < 64; ++i) {
    Float_t adc = esdV0->GetAdc(i);
    if (adc > GetRecoParam()->GetAdcThresHold()) {
      Float_t tdc = esdV0->GetTime(i);
      if (tdc > (AliVZEROReconstructor::kInvalidTime + 1e-6)) {
	Float_t nphe = adc*kChargePerADC/(cal->GetGain(i)*TMath::Qe());
	Float_t timeErr = TMath::Sqrt(kIntTimeRes*kIntTimeRes+
				      p1*p1/nphe+
				      p2*p2*(slewing->GetParameter(0)*slewing->GetParameter(1))*(slewing->GetParameter(0)*slewing->GetParameter(1))*
				      TMath::Power(adc/cal->GetCalibDiscriThr(i,kTRUE),2.*(slewing->GetParameter(1)-1.))/
				      (cal->GetCalibDiscriThr(i,kTRUE)*cal->GetCalibDiscriThr(i,kTRUE)));

	if (i < 32) { // in V0C
	  timesC[ntimeC] = tdc;
	  wC[ntimeC] = 1.0/(timeErr*timeErr);
	  indC[ntimeC] = i;
	  ntimeC++;
	  timeCW += tdc/(timeErr*timeErr);
	  weightC += 1.0/(timeErr*timeErr);
	}
	else { // in V0A
	  timesA[ntimeA] = tdc;
	  wA[ntimeA] = 1.0/(timeErr*timeErr);
	  indA[ntimeA] = i - 32;
	  ntimeA++;
	  timeAW += tdc/(timeErr*timeErr);
	  weightA += 1.0/(timeErr*timeErr);
	}
      }
    }
  } // end of loop over channels

  if (weightA > 0) timeAW = timeAW/weightA;
  else timeAW = AliVZEROReconstructor::kInvalidTime;

  if (weightC > 0) timeCW = timeCW/weightC;
  else timeCW = AliVZEROReconstructor::kInvalidTime;

  esdV0->SetBit(AliESDVZERO::kRobustMeanTime,kTRUE);

  Double_t medianTimeA = AliVZEROReconstructor::kInvalidTime;
  if (ntimeA > 0) medianTimeA = TMath::Median(ntimeA,timesA,wA);
  Double_t medianTimeC = AliVZEROReconstructor::kInvalidTime;
  if (ntimeC > 0) medianTimeC = TMath::Median(ntimeC,timesC,wC);

  Float_t robTimeAW = 0,robTimeCW = 0;
  Float_t robWeightA = 0,robWeightC = 0;
  Int_t nrobTimeA = 0, nrobTimeC = 0;
  Int_t robIndA[32], robIndC[32];
  for(Int_t i = 0; i < ntimeA; ++i) {
    AliDebug(1,Form("ChannelsAResiduals %f %f %d",timesA[i]-medianTimeA,1./TMath::Sqrt(wA[i]),ntimeA));
    if (TMath::Abs(timesA[i]-medianTimeA) < GetRecoParam()->GetMaxResid()/TMath::Sqrt(wA[i])) {
      robIndA[nrobTimeA] = indA[i];
      nrobTimeA++;
      robTimeAW += timesA[i]*wA[i];
      robWeightA += wA[i];
    }
  }
  for(Int_t i = 0; i < ntimeC; ++i) {
    AliDebug(1,Form("ChannelsCResiduals %f %f %d",timesC[i]-medianTimeC,1./TMath::Sqrt(wC[i]),ntimeC));
    if (TMath::Abs(timesC[i]-medianTimeC) < GetRecoParam()->GetMaxResid()/TMath::Sqrt(wC[i])) {
      robIndC[nrobTimeC] = indC[i];
      nrobTimeC++;
      robTimeCW += timesC[i]*wC[i];
      robWeightC += wC[i];
    }
  }

  if (robWeightA > 0) robTimeAW = robTimeAW/robWeightA;
  else robTimeAW = AliVZEROReconstructor::kInvalidTime;

  if (robWeightC > 0) robTimeCW = robTimeCW/robWeightC;
  else robTimeCW = AliVZEROReconstructor::kInvalidTime;

  AliDebug(1,Form("V0timesA %f %f %f %f %d",timeAW,(weightA > 0) ? (1./TMath::Sqrt(weightA)) : 0,
		  medianTimeA,robTimeAW,ntimeA));
  AliDebug(1,Form("V0timesC %f %f %f %f %d",timeCW,(weightC > 0) ? (1./TMath::Sqrt(weightC)) : 0,
		  medianTimeC,robTimeCW,ntimeC));

  esdV0->SetV0ATime(robTimeAW);
  esdV0->SetV0CTime(robTimeCW);
  esdV0->SetV0ATimeError((robWeightA > 0) ? (1./TMath::Sqrt(robWeightA)) : 0);
  esdV0->SetV0CTimeError((robWeightC > 0) ? (1./TMath::Sqrt(robWeightC)) : 0);

  esdV0->SetV0ADecision(AliESDVZERO::kV0Empty);
  esdV0->SetV0CDecision(AliESDVZERO::kV0Empty);

  if (robTimeAW > (fV0ADist + GetRecoParam()->GetTimeWindowBBALow()) &&
      robTimeAW < (fV0ADist + GetRecoParam()->GetTimeWindowBBAUp())) 
    esdV0->SetV0ADecision(AliESDVZERO::kV0BB);
  else if (robTimeAW > (-fV0ADist + GetRecoParam()->GetTimeWindowBGALow()) &&
	   robTimeAW < (-fV0ADist + GetRecoParam()->GetTimeWindowBGAUp()))
    esdV0->SetV0ADecision(AliESDVZERO::kV0BG);
  else if (robTimeAW > (AliVZEROReconstructor::kInvalidTime + 1e-6))
    esdV0->SetV0ADecision(AliESDVZERO::kV0Fake);

  if (robTimeCW > (fV0CDist + GetRecoParam()->GetTimeWindowBBCLow()) &&
      robTimeCW < (fV0CDist + GetRecoParam()->GetTimeWindowBBCUp())) 
    esdV0->SetV0CDecision(AliESDVZERO::kV0BB);
  else if (robTimeCW > (-fV0CDist + GetRecoParam()->GetTimeWindowBGCLow()) &&
	   robTimeCW < (-fV0CDist + GetRecoParam()->GetTimeWindowBGCUp()))
    esdV0->SetV0CDecision(AliESDVZERO::kV0BG);
  else if (robTimeCW > (AliVZEROReconstructor::kInvalidTime + 1e-6))
    esdV0->SetV0CDecision(AliESDVZERO::kV0Fake);

  UInt_t aBBtriggerV0A = 0; // bit mask for Beam-Beam trigger in V0A
  UInt_t aBGtriggerV0A = 0; // bit mask for Beam-Gas trigger in V0A
  UInt_t aBBtriggerV0C = 0; // bit mask for Beam-Beam trigger in V0C
  UInt_t aBGtriggerV0C = 0; // bit mask for Beam-Gas trigger in V0C

  for(Int_t i = 0; i < nrobTimeA; ++i) {
    if (esdV0->GetV0ADecision() == AliESDVZERO::kV0BB)
      aBBtriggerV0A |= (1 << (robIndA[i]));
    else if (esdV0->GetV0ADecision() == AliESDVZERO::kV0BG)
      aBGtriggerV0A |= (1 << (robIndA[i]));
  }

  for(Int_t i = 0; i < nrobTimeC; ++i) {
    if (esdV0->GetV0CDecision() == AliESDVZERO::kV0BB)
      aBBtriggerV0C |= (1 << (robIndC[i]));
    else if (esdV0->GetV0CDecision() == AliESDVZERO::kV0BG)
      aBGtriggerV0C |= (1 << (robIndC[i]));
  }

  esdV0->SetBBtriggerV0A(aBBtriggerV0A);
  esdV0->SetBGtriggerV0A(aBGtriggerV0A);
  esdV0->SetBBtriggerV0C(aBBtriggerV0C);
  esdV0->SetBGtriggerV0C(aBGtriggerV0C);
}
