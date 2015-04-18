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
  /*/AD has two layers, filling average 
  Float_t zADA = (TMath::Abs(GetZPosition("AD/ADA1")) + TMath::Abs(GetZPosition("AD/ADA2")))/2; 
  Float_t zADC = (TMath::Abs(GetZPosition("AD/ADC1")) + TMath::Abs(GetZPosition("AD/ADC2")))/2;

  // distance in time units from nominal vertex to V0
  fADADist = zADA/TMath::Ccgs()*1e9;
  fADCDist = zADC/TMath::Ccgs()*1e9;/*/
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
  // Fill up the trigger mask word
  // using the TDC data (already corrected for
  // slewing and misalignment between channels)

  // loop over AD channels
  Double_t timeADA =0., timeADC = 0.;
  Double_t weightADA =0., weightADC = 0.;
  UInt_t   itimeADA=0, itimeADC=0;
  
  for (Int_t i = 0; i < 16; ++i) {
    Float_t adc = esdAD->GetAdc(i);
    if (adc > GetRecoParam()->GetAdcThresHold()) {
      Float_t time = esdAD->GetTime(i);
	if(time > 1.e-6){
		Float_t timeErr = 1;
		if (adc>1) timeErr = 1/adc;

		if (i<8) {
	    		itimeADC++;
	    		timeADC += time/(timeErr*timeErr);
	    		weightADC += 1./(timeErr*timeErr);
	  		}
		else{
	    		itimeADA++;
	    		timeADA += time/(timeErr*timeErr);
	    		weightADA += 1./(timeErr*timeErr);
	  		}
	
      		}
    	}
  } // end of loop over channels
  
  if(weightADA > 1) timeADA /= weightADA; 
  else timeADA = -1024.;
  if(weightADC > 1) timeADC /= weightADC;
  else timeADC = -1024.;

  esdAD->SetADATime(timeADA);
  esdAD->SetADCTime(timeADC);
  esdAD->SetADATimeError((weightADA > 1) ? (1./TMath::Sqrt(weightADA)) : 666);
  esdAD->SetADCTimeError((weightADC > 1) ? (1./TMath::Sqrt(weightADC)) : 666);
 
}
