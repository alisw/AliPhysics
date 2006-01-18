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
//  TRD parameter class                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// WARNING: This class is obsolete. As soon as all function calls are replaced, the class will be removed.

#include <iostream>

#include "AliRun.h"

#include "AliTRDparameter.h"

ClassImp(AliTRDparameter)

//_____________________________________________________________________________
AliTRDparameter::AliTRDparameter():TNamed()
{
  //
  // AliTRDparameter default constructor
  //

  fDriftVelocity      = 0.0;

}

//_____________________________________________________________________________
AliTRDparameter::AliTRDparameter(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDparameter constructor
  //

  fDriftVelocity      = 0.0;
  
  Init();

}


//_____________________________________________________________________________
AliTRDparameter::AliTRDparameter(const AliTRDparameter &p):TNamed(p)
{
  //
  // AliTRDparameter copy constructor
  //

  ((AliTRDparameter &) p).Copy(*this);

}

///_____________________________________________________________________________
AliTRDparameter::~AliTRDparameter()
{
  //
  // AliTRDparameter destructor
  //
}

//_____________________________________________________________________________
AliTRDparameter &AliTRDparameter::operator=(const AliTRDparameter &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDparameter &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDparameter::Copy(TObject &p) const
{
  //
  // Copy function
  //

  /*((AliTRDparameter &) p).fDiffusionT         = fDiffusionT;
  ((AliTRDparameter &) p).fDiffusionL         = fDiffusionL;*/
  //((AliTRDparameter &) p).fOmegaTau           = fOmegaTau;
  //((AliTRDparameter &) p).fLorentzFactor      = fLorentzFactor;
  ((AliTRDparameter &) p).fDriftVelocity      = fDriftVelocity;


  /*((AliTRDparameter &) p).fTimeStructOn       = fTimeStructOn;
  if (((AliTRDparameter &) p).fTimeStruct1) 
    delete [] ((AliTRDparameter &) p).fTimeStruct1;
  ((AliTRDparameter &) p).fTimeStruct1 = new Float_t[38*11];
  for (Int_t i = 0; i < 38*11; i++) {
    ((AliTRDparameter &) p).fTimeStruct1[i] = fTimeStruct1[i];
  }
  if (((AliTRDparameter &) p).fTimeStruct2) 
    delete [] ((AliTRDparameter &) p).fTimeStruct2;
  ((AliTRDparameter &) p).fTimeStruct2 = new Float_t[38*11];
  for (Int_t i = 0; i < 38*11; i++) {
    ((AliTRDparameter &) p).fTimeStruct2[i] = fTimeStruct2[i];
  }*/

}

//_____________________________________________________________________________
void AliTRDparameter::Init()
{
  //
  // Initializes the parameter
  //



  //
  // ----------------------------------------------------------------------------
  // The digitization parameters
  // ----------------------------------------------------------------------------
  //

  // The drift velocity (in drift region) and the time structure of the
  // drift cells. Default is 1.5 cm/mus
  fDriftVelocity = 1.5;

  // Additional time bins before and after the drift region.
  // Default is to only sample the drift region
  fTimeMax = 22;
  SetExpandTimeBin(0,0);

  //
  // ----------------------------------------------------------------------------
  // The clusterization parameter
  // ----------------------------------------------------------------------------
  //

  ReInit();

}

//_____________________________________________________________________________
void AliTRDparameter::ReInit()
{
  //
  // Reinitializes the parameter class after a change
  //
}


//_____________________________________________________________________________
void AliTRDparameter::PrintDriftVelocity()
{
  //
  // Prints the used drift velocity
  //

  printf("<AliTRDparameter::PrintDriftVelocity> Driftvelocity = %.3f\n"
        ,fDriftVelocity);

}
