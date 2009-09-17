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
//____________________________________________________________________
//                                                                          
// This class implements the FMD offline trigger as requested for 
// ALICE first physics.
// 
// 
// 
//
#include "AliFMDOfflineTrigger.h"	
#include <iostream>
#include "AliESDFMD.h"

//____________________________________________________________________
ClassImp(AliFMDOfflineTrigger)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDOfflineTrigger::AliFMDOfflineTrigger() : 
  fLowCut(0.2),
  fHitCut(0.5)
{
  // CTOR 

}

//____________________________________________________________________
AliFMDOfflineTrigger::AliFMDOfflineTrigger(const AliFMDOfflineTrigger& o)
  : TObject(o),
    fLowCut(o.fLowCut),
    fHitCut(o.fHitCut)
{
  
  // Copy Ctor 
}

//____________________________________________________________________
AliFMDOfflineTrigger&
AliFMDOfflineTrigger::operator=(const AliFMDOfflineTrigger& o)
{
  // Assignment operator 
  return (*this);
}
//_____________________________________________________________________
Bool_t AliFMDOfflineTrigger::ASideHasHit(AliESDFMD* fmd) {

  Float_t totalMult = 0;
  for(UShort_t det=1;det<=2;det++) {
    Int_t nRings = (det == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  if(mult == AliESDFMD::kInvalidMult) continue;
	  
	  if(mult > fLowCut)
	    totalMult = totalMult + mult;
	  else
	    {
	      if( totalMult > fHitCut) {
		return kTRUE;
	      }
	      else totalMult = 0 ;
	    }
	}
      }
    }
  }
  return kFALSE;
  
}
//_____________________________________________________________________
Bool_t AliFMDOfflineTrigger::CSideHasHit(AliESDFMD* fmd) {
  
  Float_t totalMult = 0;
  UShort_t det = 3;
  Int_t nRings = 2;
  for (UShort_t ir = 0; ir < nRings; ir++) {
    Char_t   ring = (ir == 0 ? 'I' : 'O');
    UShort_t nsec = (ir == 0 ? 20  : 40);
    UShort_t nstr = (ir == 0 ? 512 : 256);
    for(UShort_t sec =0; sec < nsec;  sec++)  {
      for(UShort_t strip = 0; strip < nstr; strip++) {
	Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	if(mult == AliESDFMD::kInvalidMult) continue;
	
	if(mult > fLowCut)
	  totalMult = totalMult + mult;
	else
	  {
	    if( totalMult > fHitCut) {
	      return kTRUE;
	    }
	    else totalMult = 0 ;
	  }
      }
    }
  }
 
  return kFALSE;
}
//____________________________________________________________________
//
// EOF
//
