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
// Class : AliADDecision
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
#include "AliADDecision.h"
#include "AliADConst.h"
#include "AliADCalibData.h"
#include "AliESDAD.h"
#include "AliADReconstructor.h"

//______________________________________________________________________
ClassImp(AliADDecision)

//______________________________________________________________________

AliADDecision::AliADDecision()
  :TObject(),
   fADADist(0),
   fADCDist(0),
   fRecoParam(NULL),
   fEarlyHitCutShape(NULL)
{
  // Default constructor
  //AD has two layers, filling average 
  Float_t zADA = (TMath::Abs(GetZPosition("AD/ADA1")) + TMath::Abs(GetZPosition("AD/ADA2")))/2; 
  Float_t zADC = (TMath::Abs(GetZPosition("AD/ADC1")) + TMath::Abs(GetZPosition("AD/ADC2")))/2;

  // distance in time units from nominal vertex to AD
  fADADist = zADA/TMath::Ccgs()*1e9; 
  fADCDist = zADC/TMath::Ccgs()*1e9;
  
  //fADADist = 56.6958892608299081;
  //fADCDist = 65.1917667655268360;
  
  fEarlyHitCutShape = new TF1("fEarlyHitCutShape", " [0]+(x>[2])*[1]*(x-[2])**2");
  
}

//______________________________________________________________________
AliADDecision::~AliADDecision()
{
  // d-tor
  delete fEarlyHitCutShape;
}


//________________________________________________________________________________
Double_t AliADDecision::GetZPosition(const char* symname)
{
// Get the global z coordinate of the given AD alignable volume
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
void AliADDecision::FillDecisions(AliESDAD *esdAD)
{
  // Fill up offline trigger decisions
  // using the TDC data (already corrected for
  // slewing and misalignment between channels)
  
  esdAD->SetBit(AliESDAD::kDecisionFilled,kTRUE);
  esdAD->SetBit(AliESDAD::kTriggerBitsFilled,kTRUE);

  // loop over AD channels
  Double_t timeADA =0., timeADC = 0.;
  Double_t weightADA =0., weightADC = 0.;
  UInt_t   ntimeADA=0, ntimeADC=0;
  UInt_t   itimeADA[8], itimeADC[8];
  
  //Compute average time with basic method(charge weighted average)
  Double_t timeBasicADA=0, timeBasicADC=0;
  Double_t weightBasicADA =0., weightBasicADC = 0.;
  UInt_t   ntimeBasicADA=0, ntimeBasicADC=0;
  UInt_t   itimeBasicADA[8], itimeBasicADC[8];
  
  for (Int_t i = 0; i < 16; ++i) {
    Float_t adc = esdAD->GetAdc(i);
    if (adc > GetRecoParam()->GetAdcThresHold()) {
      Float_t time = esdAD->GetTime(i);
	if(time > (AliADReconstructor::kInvalidTime+1.e-6)){
		Float_t timeErr = 1;
		if (adc>1) timeErr = 1/adc;

		if (i<8) {
			itimeBasicADC[ntimeBasicADC] = i;
	    		ntimeBasicADC++;
	    		timeBasicADC += time/(timeErr*timeErr);
	    		weightBasicADC += 1./(timeErr*timeErr);
	  		}
		else{
			itimeBasicADA[ntimeBasicADA] = i-8;
	    		ntimeBasicADA++;
	    		timeBasicADA += time/(timeErr*timeErr);
	    		weightBasicADA += 1./(timeErr*timeErr);
	  		}
	
      		}
    	}
  } // end of loop over channels
  
  
  if(weightBasicADA > 1) timeBasicADA /= weightBasicADA; 
  else timeBasicADA = -1024.;
  if(weightBasicADC > 1) timeBasicADC /= weightBasicADC;
  else timeBasicADC = -1024.;
  
  //Robust time: Pad coincidence and early hit removal
  Double_t timeRobustADA=0, timeRobustADC=0;
  Double_t weightRobustADA =0., weightRobustADC = 0.;
  UInt_t   ntimeRobustADA=0, ntimeRobustADC=0;
  UInt_t   itimeRobustADA[8], itimeRobustADC[8];
  
  fEarlyHitCutShape->SetParameter(0,GetRecoParam()->GetMaxResid());
  fEarlyHitCutShape->SetParameter(1,GetRecoParam()->GetResidRise());
  
  //C-side
  fEarlyHitCutShape->SetParameter(2,2*fADCDist - 10);
  for (Int_t i = 0; i < 4; ++i) {
    Float_t adc1 = esdAD->GetAdc(i);
    Float_t adc2 = esdAD->GetAdc(i+4);
    if (adc1 > GetRecoParam()->GetAdcThresHold() && adc2 > GetRecoParam()->GetAdcThresHold()) {
      Float_t time1 = esdAD->GetTime(i);
      Float_t time2 = esdAD->GetTime(i+4);
	if(time1 > (AliADReconstructor::kInvalidTime+1.e-6) && time2 > (AliADReconstructor::kInvalidTime+1.e-6)){
		Float_t timeErr1 = 1/adc1;
		Float_t timeErr2 = 1/adc2;
		Float_t timeDiff = TMath::Abs(time1-time2);
		Float_t timeSum = time1+time2;
		Float_t earlyHitCut = 1000;
		if(TMath::Abs(timeSum - 2*fADCDist) < 20) earlyHitCut = fEarlyHitCutShape->Eval(timeSum);
		if(timeDiff < earlyHitCut){
			itimeRobustADC[ntimeRobustADC] = i;
			ntimeRobustADC++;
			timeRobustADC += time1/(timeErr1*timeErr1);
			weightRobustADC += 1./(timeErr1*timeErr1);
			
			itimeRobustADC[ntimeRobustADC] = i+4;
			ntimeRobustADC++;
			timeRobustADC += time2/(timeErr2*timeErr2);
			weightRobustADC += 1./(timeErr2*timeErr2);
			}
		}
	}
		
  }
  
  //A-side
  fEarlyHitCutShape->SetParameter(2,2*fADADist - 10);
  for (Int_t i = 8; i < 12; ++i) {
    Float_t adc1 = esdAD->GetAdc(i);
    Float_t adc2 = esdAD->GetAdc(i+4);
    if (adc1 > GetRecoParam()->GetAdcThresHold() && adc2 > GetRecoParam()->GetAdcThresHold()) {
      Float_t time1 = esdAD->GetTime(i);
      Float_t time2 = esdAD->GetTime(i+4);
	if(time1 > (AliADReconstructor::kInvalidTime+1.e-6) && time2 > (AliADReconstructor::kInvalidTime+1.e-6)){
		Float_t timeErr1 = 1/adc1;
		Float_t timeErr2 = 1/adc2;
		Float_t timeDiff = TMath::Abs(time1-time2);
		Float_t timeSum = time1+time2;
		Float_t earlyHitCut = 1000;
		if(TMath::Abs(timeSum - 2*fADADist) < 20) earlyHitCut = fEarlyHitCutShape->Eval(timeSum);
		if(timeDiff < earlyHitCut){
			itimeRobustADA[ntimeRobustADA] = i-8;
			ntimeRobustADA++;
			timeRobustADA += time1/(timeErr1*timeErr1);
			weightRobustADA += 1./(timeErr1*timeErr1);
			
			itimeRobustADA[ntimeRobustADA] = i-4;
			ntimeRobustADA++;
			timeRobustADA += time2/(timeErr2*timeErr2);
			weightRobustADA += 1./(timeErr2*timeErr2);
			}
		}
	}
		
  }
  
  if(weightRobustADA > 1) timeRobustADA /= weightRobustADA; 
  else timeRobustADA = -1024.;
  if(weightRobustADC > 1) timeRobustADC /= weightRobustADC;
  else timeRobustADC = -1024.;
  
  /*/Use basic time for decisions
  timeADA = timeBasicADA; timeADC = timeBasicADC;
  weightADA = weightBasicADA; weightADC = weightBasicADC;
  ntimeADA = ntimeBasicADA; ntimeADC = ntimeBasicADC;
  for(Int_t i=0; i<8; ++i){itimeADA[i] = itimeBasicADA[i]; itimeADC[i] = itimeBasicADC[i];}/*/
  
  //Use Robust time for decisions
  esdAD->SetBit(AliESDAD::kRobustMeanTime,kTRUE);
  timeADA = timeRobustADA; timeADC = timeRobustADC;
  weightADA = weightRobustADA; weightADC = weightRobustADC;
  ntimeADA = ntimeRobustADA; ntimeADC = ntimeRobustADC;
  for(Int_t i=0; i<8; ++i){itimeADA[i] = itimeRobustADA[i]; itimeADC[i] = itimeRobustADC[i];}

  esdAD->SetADATime(timeADA);
  esdAD->SetADCTime(timeADC);
  esdAD->SetADATimeError((weightADA > 1) ? (1./TMath::Sqrt(weightADA)) : 666);
  esdAD->SetADCTimeError((weightADC > 1) ? (1./TMath::Sqrt(weightADC)) : 666);
  
  esdAD->SetADADecision(AliESDAD::kADEmpty);
  esdAD->SetADCDecision(AliESDAD::kADEmpty);

  if (timeADA > (fADADist + GetRecoParam()->GetTimeWindowBBALow()) &&
      timeADA < (fADADist + GetRecoParam()->GetTimeWindowBBAUp())) 
    esdAD->SetADADecision(AliESDAD::kADBB);
  else if (timeADA > (-fADADist + GetRecoParam()->GetTimeWindowBGALow()) &&
	   timeADA < (-fADADist + GetRecoParam()->GetTimeWindowBGAUp()))
    esdAD->SetADADecision(AliESDAD::kADBG);
  else if (timeADA > (AliADReconstructor::kInvalidTime + 1e-6))
    esdAD->SetADADecision(AliESDAD::kADFake);

  if (timeADC > (fADCDist + GetRecoParam()->GetTimeWindowBBCLow()) &&
      timeADC < (fADCDist + GetRecoParam()->GetTimeWindowBBCUp())) 
    esdAD->SetADCDecision(AliESDAD::kADBB);
  else if (timeADC > (-fADCDist + GetRecoParam()->GetTimeWindowBGCLow()) &&
	   timeADC < (-fADCDist + GetRecoParam()->GetTimeWindowBGCUp()))
    esdAD->SetADCDecision(AliESDAD::kADBG);
  else if (timeADC > (AliADReconstructor::kInvalidTime + 1e-6))
    esdAD->SetADCDecision(AliESDAD::kADFake);

  UInt_t aBBtriggerADA = 0; // bit mask for Beam-Beam trigger in ADA
  UInt_t aBGtriggerADA = 0; // bit mask for Beam-Gas trigger in ADA
  UInt_t aBBtriggerADC = 0; // bit mask for Beam-Beam trigger in ADC
  UInt_t aBGtriggerADC = 0; // bit mask for Beam-Gas trigger in ADC

  for(Int_t i = 0; i < ntimeADA; ++i) {
    if (esdAD->GetADADecision() == AliESDAD::kADBB)
      aBBtriggerADA |= (1 << (itimeADA[i]));
    else if (esdAD->GetADADecision() == AliESDAD::kADBG)
      aBGtriggerADA |= (1 << (itimeADA[i]));
  }

  for(Int_t i = 0; i < ntimeADC; ++i) {
    if (esdAD->GetADCDecision() == AliESDAD::kADBB)
      aBBtriggerADC |= (1 << (itimeADC[i]));
    else if (esdAD->GetADCDecision() == AliESDAD::kADBG)
      aBGtriggerADC |= (1 << (itimeADC[i]));
  }

  esdAD->SetBBtriggerADA(aBBtriggerADA);
  esdAD->SetBGtriggerADA(aBGtriggerADA);
  esdAD->SetBBtriggerADC(aBBtriggerADC);
  esdAD->SetBGtriggerADC(aBGtriggerADC);

 
}
