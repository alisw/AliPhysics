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

/* $Id: AliTRDCalDCS.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS parameters                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCS.h"
#include "AliTRDCalDCSFEE.h"

ClassImp(AliTRDCalDCS)
  
//_____________________________________________________________________________
AliTRDCalDCS::AliTRDCalDCS()
  :TNamed()
  ,fGNumberOfTimeBins(-1)
  ,fGConfigTag(-1)
  ,fGSingleHitThres(-1)
  ,fGThreePadClustThres(-1)
  ,fGSelNoZS(-1)
  ,fGTCFilterWeight(-1)
  ,fGTCFilterShortDecPar(-1)
  ,fGTCFilterLongDecPar(-1)
  ,fGFastStatNoise(-1)
  ,fGConfigVersion(0)
  ,fGConfigName(0)
  ,fGFilterType(0)
  ,fGReadoutParam(0)
  ,fGTestPattern(0)
  ,fGTrackletMode(0)
  ,fGTrackletDef(0)
  ,fGTriggerSetup(0)
  ,fGAddOptions(0)
  ,fFEEArr(new TObjArray(540))
  ,fPTRArr(new TObjArray(6))
  ,fGTUArr(new TObjArray(19))
{
  //
  // AliTRDCalDCS default constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCS::AliTRDCalDCS(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
  ,fGNumberOfTimeBins(-1)
  ,fGConfigTag(-1)
  ,fGSingleHitThres(-1)
  ,fGThreePadClustThres(-1)
  ,fGSelNoZS(-1)
  ,fGTCFilterWeight(-1)
  ,fGTCFilterShortDecPar(-1)
  ,fGTCFilterLongDecPar(-1)
  ,fGFastStatNoise(-1)
  ,fGConfigVersion(0)
  ,fGConfigName(0)
  ,fGFilterType(0)
  ,fGReadoutParam(0)
  ,fGTestPattern(0)
  ,fGTrackletMode(0)
  ,fGTrackletDef(0)
  ,fGTriggerSetup(0)
  ,fGAddOptions(0)
  ,fFEEArr(new TObjArray(540))
  ,fPTRArr(new TObjArray(6))
  ,fGTUArr(new TObjArray(19))
{
  //
  // AliTRDCalDCS constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCS::AliTRDCalDCS(const AliTRDCalDCS &cd)
  :TNamed(cd)
  ,fGNumberOfTimeBins(-1)
  ,fGConfigTag(-1)
  ,fGSingleHitThres(-1)
  ,fGThreePadClustThres(-1)
  ,fGSelNoZS(-1)
  ,fGTCFilterWeight(-1)
  ,fGTCFilterShortDecPar(-1)
  ,fGTCFilterLongDecPar(-1)
  ,fGFastStatNoise(-1)
  ,fGConfigVersion(0)
  ,fGConfigName(0)
  ,fGFilterType(0)
  ,fGReadoutParam(0)
  ,fGTestPattern(0)
  ,fGTrackletMode(0)
  ,fGTrackletDef(0)
  ,fGTriggerSetup(0)
  ,fGAddOptions(0)
  ,fFEEArr(0)
  ,fPTRArr(0)
  ,fGTUArr(0)
{
  //
  // AliTRDCalDCS copy constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCS &AliTRDCalDCS::operator=(const AliTRDCalDCS &cd)
{
  //
  // Assignment operator
  //
  if (&cd == this) return *this;

  new (this) AliTRDCalDCS(cd);
  return *this;
}

//_____________________________________________________________________________
void AliTRDCalDCS::EvaluateGlobalParameters()
{
  for(Int_t i=0; i<540; i++) {
    AliTRDCalDCSFEE *iDCSFEEObj;
    iDCSFEEObj = GetCalDCSFEEObj(i);
    if(iDCSFEEObj != NULL) {
      if(iDCSFEEObj->GetStatusBit() == 0) {
	// first, set the parameters of the first good ROC as global
	fGNumberOfTimeBins    = iDCSFEEObj->GetNumberOfTimeBins();
	fGConfigTag           = iDCSFEEObj->GetConfigTag();
	fGSingleHitThres      = iDCSFEEObj->GetSingleHitThres();
	fGThreePadClustThres  = iDCSFEEObj->GetThreePadClustThres();
	fGSelNoZS             = iDCSFEEObj->GetSelectiveNoZS();
	fGTCFilterWeight      = iDCSFEEObj->GetTCFilterWeight();
	fGTCFilterShortDecPar = iDCSFEEObj->GetTCFilterShortDecPar();
	fGTCFilterLongDecPar  = iDCSFEEObj->GetTCFilterLongDecPar();
	fGFastStatNoise       = iDCSFEEObj->GetFastStatNoise();
	fGConfigVersion       = iDCSFEEObj->GetConfigVersion();
	fGConfigName          = iDCSFEEObj->GetConfigName();
	fGFilterType          = iDCSFEEObj->GetFilterType();
	fGReadoutParam        = iDCSFEEObj->GetReadoutParam();
	fGTestPattern         = iDCSFEEObj->GetTestPattern();
	fGTrackletMode        = iDCSFEEObj->GetTrackletMode();
	fGTrackletDef         = iDCSFEEObj->GetTrackletDef();
	fGTriggerSetup        = iDCSFEEObj->GetTriggerSetup();
	fGAddOptions          = iDCSFEEObj->GetAddOptions();
	break;
      }
    }
  }

  for(Int_t i=0; i<540; i++) {
    AliTRDCalDCSFEE *iDCSFEEObj;
    iDCSFEEObj = GetCalDCSFEEObj(i);
    if(iDCSFEEObj != NULL) {
      if(iDCSFEEObj->GetStatusBit() == 0) {
	// second, if any of the other good chambers differ, set the global value to -1/""
	if(fGNumberOfTimeBins    != iDCSFEEObj->GetNumberOfTimeBins())    fGNumberOfTimeBins    = -2;
	if(fGConfigTag           != iDCSFEEObj->GetConfigTag())           fGConfigTag           = -2;
	if(fGSingleHitThres      != iDCSFEEObj->GetSingleHitThres())      fGSingleHitThres      = -2;
	if(fGThreePadClustThres  != iDCSFEEObj->GetThreePadClustThres())  fGThreePadClustThres  = -2;
	if(fGSelNoZS             != iDCSFEEObj->GetSelectiveNoZS())       fGSelNoZS             = -2;
	if(fGTCFilterWeight      != iDCSFEEObj->GetTCFilterWeight())      fGTCFilterWeight      = -2;
	if(fGTCFilterShortDecPar != iDCSFEEObj->GetTCFilterShortDecPar()) fGTCFilterShortDecPar = -2;
	if(fGTCFilterLongDecPar  != iDCSFEEObj->GetTCFilterLongDecPar())  fGTCFilterLongDecPar  = -2;
	if(fGFastStatNoise       != iDCSFEEObj->GetFastStatNoise())       fGFastStatNoise       = -2;
	if(fGConfigVersion       != iDCSFEEObj->GetConfigVersion())       fGConfigVersion       = "mixed";
	if(fGConfigName          != iDCSFEEObj->GetConfigName())          fGConfigName          = "mixed";
	if(fGFilterType          != iDCSFEEObj->GetFilterType())          fGFilterType          = "mixed";
	if(fGReadoutParam        != iDCSFEEObj->GetReadoutParam())        fGReadoutParam        = "mixed";
	if(fGTestPattern         != iDCSFEEObj->GetTestPattern())         fGTestPattern         = "mixed";
	if(fGTrackletMode        != iDCSFEEObj->GetTrackletMode())        fGTrackletMode        = "mixed";
	if(fGTrackletDef         != iDCSFEEObj->GetTrackletDef())         fGTrackletDef         = "mixed";
	if(fGTriggerSetup        != iDCSFEEObj->GetTriggerSetup())        fGTriggerSetup        = "mixed";
	if(fGAddOptions          != iDCSFEEObj->GetAddOptions())          fGAddOptions          = "mixed";
      }
    }
  }
}

