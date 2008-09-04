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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class containing constant common parameters                               //
//                                                                           //
// Request an instance with AliTRDCommonParam::Instance()                    //
// Then request the needed values                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>

#include "AliTracker.h"
#include "AliRun.h"

#include "AliTRDCommonParam.h"

ClassImp(AliTRDCommonParam)

AliTRDCommonParam *AliTRDCommonParam::fgInstance = 0;
Bool_t AliTRDCommonParam::fgTerminated = kFALSE;

//_ singleton implementation __________________________________________________
AliTRDCommonParam* AliTRDCommonParam::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  // 
  
  if (fgTerminated != kFALSE) {
    return 0;
  }

  if (fgInstance == 0) {
    fgInstance = new AliTRDCommonParam();
  }  

  return fgInstance;

}

//_____________________________________________________________________________
void AliTRDCommonParam::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag,
  // instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != 0) {
    delete fgInstance;
    fgInstance = 0;
  }

}

//_____________________________________________________________________________
AliTRDCommonParam::AliTRDCommonParam()
  :TObject()
  ,fExBOn(kFALSE)
  ,fSamplingFrequency(0.0)
{
  //
  // Default constructor
  //
  
  Init();

}

//_____________________________________________________________________________
AliTRDCommonParam::AliTRDCommonParam(TRootIoCtor *)
  :TObject()
  ,fExBOn(0)
  ,fSamplingFrequency(0.0)
{
  //
  // IO constructor
  //

}

//_____________________________________________________________________________
void AliTRDCommonParam::Init()
{
  //
  // Initialization
  //
  
  // E x B effects
  fExBOn             = kTRUE;

  // Sampling Frequency in MHz
  fSamplingFrequency = 10.0;

}

//_____________________________________________________________________________
AliTRDCommonParam::~AliTRDCommonParam() 
{
  //
  // Destructor
  //
  
}

//_____________________________________________________________________________
AliTRDCommonParam::AliTRDCommonParam(const AliTRDCommonParam &p)
  :TObject(p)
  ,fExBOn(p.fExBOn)
  ,fSamplingFrequency(p.fSamplingFrequency)
{
  //
  // Copy constructor
  //

}

//_____________________________________________________________________________
AliTRDCommonParam &AliTRDCommonParam::operator=(const AliTRDCommonParam &p)
{
  //
  // Assignment operator
  //

  if (this != &p) {
    ((AliTRDCommonParam &) p).Copy(*this);
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDCommonParam::Copy(TObject &p) const
{
  //
  // Copy function
  //
  
  AliTRDCommonParam *target = dynamic_cast<AliTRDCommonParam*> (&p);
  if (!target) {
    return;
  }  

  target->fExBOn             = fExBOn;
  target->fSamplingFrequency = fSamplingFrequency;

}
