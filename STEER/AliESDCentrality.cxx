/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//*****************************************************
//   Class AliCentralitySelectionTask
//   author: Alberica Toia
//*****************************************************
/// A container for the centrality stored in the ESD.
 
#include "AliESDCentrality.h"

ClassImp(AliESDCentrality)

AliESDCentrality::AliESDCentrality() : TNamed("ESDCentrality", "Centrality"),
fCentrality(0)
{
  /// constructor
}

AliESDCentrality::AliESDCentrality(const AliESDCentrality& cnt) : 
  TNamed(cnt), fCentrality(cnt.fCentrality)
{
  /// Copy constructor
}

AliESDCentrality& AliESDCentrality::operator=(const AliESDCentrality& c)
{
  /// Assignment operator
  if (this!=&c) {
    TNamed::operator=(c);
    fCentrality = c.fCentrality;
  }

  return *this;
}

AliESDCentrality::~AliESDCentrality()
{
  /// destructor
}

Float_t AliESDCentrality::GetCentralityPercentile()
{
  return fCentrality;
}

Int_t AliESDCentrality::GetCentralityClass10()
{
  return (Int_t) fCentrality / 10.0;
}

Int_t AliESDCentrality::GetCentralityClass5()
{
  return (Int_t) fCentrality / 5.0;
}

Bool_t AliESDCentrality::IsEventInCentralityClass(Float_t a, Float_t b)
{
  if (fCentrality >=a && fCentrality < b) return kTRUE;
  else return kFALSE;
}

