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
   fRecoParam(NULL)
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
	if(time > 1.e-6){
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
  
  //Use basic time for decisions
  timeADA = timeBasicADA; timeADC = timeBasicADC;
  weightADA = weightBasicADA; weightADC = weightBasicADC;
  ntimeADA = ntimeBasicADA; ntimeADC = ntimeBasicADC;
  for(Int_t i=0; i<8; ++i){itimeADA[i] = itimeBasicADA[i]; itimeADC[i] = itimeBasicADC[i];} 

  esdAD->SetADATime(timeADA);
  esdAD->SetADCTime(timeADC);
  esdAD->SetADATimeError((weightADA > 1) ? (1./TMath::Sqrt(weightADA)) : 666);
  esdAD->SetADCTimeError((weightADC > 1) ? (1./TMath::Sqrt(weightADC)) : 666);
  
  esdAD->SetADADecision(AliESDAD::kADEmpty);
  esdAD->SetADCDecision(AliESDAD::kADEmpty);

  if (timeADA > (fADADist + GetRecoParam()->GetTimeWindowBBALow()) &&
      timeADA < (fADADist + GetRecoParam()->GetTimeWindowBBAUp())) 
    esdAD->SetADADecision(AliESDAD::kADBB);
  else if (timeADA > (fADADist + GetRecoParam()->GetTimeWindowBGALow()) &&
	   timeADA < (fADADist + GetRecoParam()->GetTimeWindowBGAUp()))
    esdAD->SetADADecision(AliESDAD::kADBG);
  else if (timeADA > (AliADReconstructor::kInvalidTime + 1e-6))
    esdAD->SetADADecision(AliESDAD::kADFake);

  if (timeADC > (fADCDist + GetRecoParam()->GetTimeWindowBBCLow()) &&
      timeADC < (fADCDist + GetRecoParam()->GetTimeWindowBBCUp())) 
    esdAD->SetADCDecision(AliESDAD::kADBB);
  else if (timeADC > (fADCDist + GetRecoParam()->GetTimeWindowBGCLow()) &&
	   timeADC < (fADCDist + GetRecoParam()->GetTimeWindowBGCUp()))
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
