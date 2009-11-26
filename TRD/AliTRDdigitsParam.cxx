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
  ,fNTimeBins(0)
  ,fADCbaseline(0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 540; i++) {
    fPretriggerPhase[i] = 0;
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
  ,fNTimeBins(p.fNTimeBins)
  ,fADCbaseline(p.fADCbaseline)
{
  //
  // Copy constructor
  //

  for (Int_t i = 0; i < 540; i++) {
    fPretriggerPhase[i] = p.fPretriggerPhase[i];
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

  target->fNTimeBins   = fNTimeBins;
  target->fADCbaseline = fADCbaseline;

  for (Int_t i = 0; i < 540; i++) {
    target->fPretriggerPhase[i] = fPretriggerPhase[i];
  }

}
