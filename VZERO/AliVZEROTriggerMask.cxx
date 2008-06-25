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
#include <Riostream.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <AliGeomManager.h>
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliVZEROTriggerMask.h"
#include "AliVZEROdigit.h"

//______________________________________________________________________
ClassImp(AliVZEROTriggerMask)

//______________________________________________________________________

AliVZEROTriggerMask::AliVZEROTriggerMask()
  :TObject(),
   fAdcThresHold(0.0),
   fTimeWindowWidthBBA(50.0),
   fTimeWindowWidthBGA(20.0),
   fTimeWindowWidthBBC(50.0),
   fTimeWindowWidthBGC(20.0),
   fBBtriggerV0A(0),
   fBGtriggerV0A(0),
   fBBtriggerV0C(0),
   fBGtriggerV0C(0)
   
{
   SetAdcThreshold();
}

//________________________________________________________________________________
Double_t AliVZEROTriggerMask::GetZPosition(const char* symname){
// Get the global z coordinate of the given V0 alignable volume
//
  Double_t *tr;
  TGeoPNEntry *pne = gGeoManager->GetAlignableEntry(symname);
  if (!pne) return 0;


  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if(pnode){
          TGeoHMatrix* hm = pnode->GetMatrix();
           tr = hm->GetTranslation();
  }else{
          const char* path = pne->GetTitle();
          if(!gGeoManager->cd(path)){
                  AliErrorClass(Form("Volume path %s not valid!",path));
                  return 0;
          }
         tr = gGeoManager->GetCurrentMatrix()->GetTranslation();
  }
  return tr[2];

}

//________________________________________________________________________________


void AliVZEROTriggerMask::FillMasks(TTree* vzeroDigitsTree,
				TClonesArray* vzeroDigits)
{
  const Double_t LightSpeed = 2.9979245800; // cm/100 ps
  Float_t Z_V0A = TMath::Abs(GetZPosition("VZERO/V0A"));
  Float_t Z_V0C = TMath::Abs(GetZPosition("VZERO/V0C"));

  // distance in time units from nominal vertex to V0
  Float_t v0a_dist = Z_V0A/LightSpeed; // 100 of picoseconds
  Float_t v0c_dist = Z_V0C/LightSpeed; // 100 of picoseconds
  Float_t bunch_separation = 1000.0; // 100 of picoseconds

  // mask
  UInt_t one=1;
 
  // loop over vzero entries
  Int_t nEntries = (Int_t)vzeroDigitsTree->GetEntries();
  for (Int_t e=0; e<nEntries; e++) {
    vzeroDigitsTree->GetEvent(e);

    Int_t nDigits = vzeroDigits->GetEntriesFast();
    
    for (Int_t d=0; d<nDigits; d++) {
      //      vzeroDigitsTree->GetEvent(d);
      AliVZEROdigit* digit = (AliVZEROdigit*)vzeroDigits->At(d);
      
      Int_t   PMNumber   = digit->PMNumber();
      Float_t adc        = digit->ADC();
      Float_t tdc        = digit->Time(); // in 100 of picoseconds

      if (adc>fAdcThresHold) {
	if (PMNumber<32) { // in V0C
	  if (tdc>(v0c_dist-fTimeWindowWidthBBC/2.0) &&
	      tdc<(v0c_dist+fTimeWindowWidthBBC/2.0))
	    fBBtriggerV0C+=(one<<PMNumber);
	  if (tdc>(bunch_separation-v0c_dist-fTimeWindowWidthBGC/2.0) &&
	      tdc<(bunch_separation-v0c_dist+fTimeWindowWidthBGC/2.0))
	   fBGtriggerV0C+=(one<<PMNumber); 
	}
	if (PMNumber>31) { // in V0A
	  Int_t shift = PMNumber-32;
	  if (tdc>(v0a_dist-fTimeWindowWidthBBA/2.0) &&
	      tdc<(v0a_dist+fTimeWindowWidthBBA/2.0)) 
	    fBBtriggerV0A+=(one<<shift);
	  if (tdc>(bunch_separation-v0a_dist-fTimeWindowWidthBGA/2.0) &&
	      tdc<(bunch_separation-v0a_dist+fTimeWindowWidthBGA/2.0))
	    fBGtriggerV0A+=(one<<shift); 
	}
      }
    } // end of loop over digits
  } // end of loop over events in digits tree
}
