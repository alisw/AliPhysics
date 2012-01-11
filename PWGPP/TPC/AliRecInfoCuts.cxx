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

//------------------------------------------------------------------------------
// Implementation of the AliRecInfoCuts class. It keeps selection cuts for 
// reconstructed tracks. 
//
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

#include "AliRecInfoCuts.h"

ClassImp(AliRecInfoCuts)

//_____________________________________________________________________________
AliRecInfoCuts::AliRecInfoCuts(const Char_t* name,const Char_t *title) : AliESDtrackCuts(name, title)
, fMinTPCsignalN(0)
, fMaxAbsTanTheta(0)
, fMinNClustersTRD(0)
, fTPCITSMatchingRadius(0)
, fTPCTRDMatchingRadius(0)
, fMinNTrackletsTRD(0)
{
  // init data members with defaults
  InitME();
}

//_____________________________________________________________________________
void AliRecInfoCuts::InitME()
{
  // set default values 
  SetMinTPCsignalN();
  SetMaxAbsTanTheta();
  SetMinNClustersTRD();
  SetMinNTrackletsTRD();
  SetTPCITSMatchingRadius();
  SetTPCTRDMatchingRadius();
}


//_____________________________________________________________________________
/*
Long64_t AliRecInfoCuts::Merge(TCollection* list) const 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
  AliRecInfoCuts* entry = dynamic_cast<AliRecInfoCuts*>(obj);
  if (entry == 0) 
   continue;

  count++;
  }

return count;
}
*/

