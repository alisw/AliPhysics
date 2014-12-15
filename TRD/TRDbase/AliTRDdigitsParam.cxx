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

/* $Id: AliTRDdigitsParam.cxx 34070 2009-08-04 15:34:53Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class containing parameters for digits                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"

#include "AliTRDdigitsParam.h"

ClassImp(AliTRDdigitsParam)

//_____________________________________________________________________________
AliTRDdigitsParam::AliTRDdigitsParam()
  :TObject()
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 540; i++) {
    fNTimeBins[i]       = 0;
    fPretriggerPhase[i] = 0;
    fADCbaseline[i]     = 0;
  }

}

//_____________________________________________________________________________
AliTRDdigitsParam::~AliTRDdigitsParam() 
{
  //
  // Destructor
  //

}

//_____________________________________________________________________________
AliTRDdigitsParam::AliTRDdigitsParam(const AliTRDdigitsParam &p)
  :TObject(p)
{
  //
  // Copy constructor
  //

  for (Int_t i = 0; i < 540; i++) {
    fNTimeBins[i]       = p.fNTimeBins[i];
    fPretriggerPhase[i] = p.fPretriggerPhase[i];
    fADCbaseline[i]     = p.fADCbaseline[i];
  }

}

//_____________________________________________________________________________
AliTRDdigitsParam &AliTRDdigitsParam::operator=(const AliTRDdigitsParam &p)
{
  //
  // Assignment operator
  //

  if (this != &p) {
    ((AliTRDdigitsParam &) p).Copy(*this);
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDdigitsParam::Copy(TObject &p) const
{
  //
  // Copy function
  //
  
  AliTRDdigitsParam *target = dynamic_cast<AliTRDdigitsParam*> (&p);
  if (!target) {
    return;
  }  

  for (Int_t i = 0; i < 540; i++) {
    target->fNTimeBins[i]       = fNTimeBins[i];
    target->fPretriggerPhase[i] = fPretriggerPhase[i];
    target->fADCbaseline[i]     = fADCbaseline[i];
  }

}
