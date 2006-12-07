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
#include "AliTRDpadPlane.h"


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
  // Deletes the instance of this class and sets the terminated flag, instances cannot be requested anymore
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
  ,fField(0)
  ,fExBOn(kFALSE)
  ,fPadPlaneArray(0)
{
  //
  // Default constructor
  //
  
  Init();

}

//_____________________________________________________________________________
void AliTRDCommonParam::Init()
{
  //
  // Initialization
  //
  
  // E x B effects
  fExBOn          = kTRUE;

  // The magnetic field strength in Tesla
  fField = AliTracker::GetBz() * 0.1; 

  if (TMath::Abs(fField) < 1e-5) {
    Info("Init", "MC B field ... ");
    Double_t x[3] = { 0.0, 0.0, 0.0 };
    Double_t b[3]; 	 
    gAlice->Field(x,b);  // b[] is in kilo Gauss 	 
    fField = b[2] * 0.1; // Tesla
  }

  
  // ----------------------------------------------------------------------------
  // The pad planes
  // ----------------------------------------------------------------------------
  
  fPadPlaneArray = new TObjArray(kNplan * kNcham);
  
  for (Int_t iplan = 0; iplan < kNplan; iplan++) {
    for (Int_t icham = 0; icham < kNcham; icham++) {
      Int_t ipp = iplan + icham * kNplan;
      fPadPlaneArray->AddAt(new AliTRDpadPlane(iplan,icham),ipp);
    }
  }

}

//_____________________________________________________________________________
AliTRDCommonParam::~AliTRDCommonParam() 
{
  //
  // Destructor
  //
  
  if (fPadPlaneArray) {
    fPadPlaneArray->Delete();
    delete fPadPlaneArray;
    fPadPlaneArray = 0;
  }

}

//_____________________________________________________________________________
AliTRDCommonParam::AliTRDCommonParam(const AliTRDCommonParam &p)
  :TObject(p)
  ,fField(p.fField)
  ,fExBOn(p.fExBOn)
  ,fPadPlaneArray(0)
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

  target->fExBOn = fExBOn;
  target->fField = fField;

}

//_____________________________________________________________________________
AliTRDpadPlane *AliTRDCommonParam::GetPadPlane(Int_t p, Int_t c) const
{
  //
  // Returns the pad plane for a given plane <p> and chamber <c> number
  //

  Int_t ipp = p + c * kNplan;
  return ((AliTRDpadPlane *) fPadPlaneArray->At(ipp));

}

//_____________________________________________________________________________
Int_t AliTRDCommonParam::GetRowMax(Int_t p, Int_t c, Int_t /*s*/) const
{
  //
  // Returns the number of rows on the pad plane
  //

  return GetPadPlane(p,c)->GetNrows();

}

//_____________________________________________________________________________
Int_t AliTRDCommonParam::GetColMax(Int_t p) const
{
  //
  // Returns the number of rows on the pad plane
  //

  return GetPadPlane(p,0)->GetNcols();

}

//_____________________________________________________________________________
Double_t AliTRDCommonParam::GetRow0(Int_t p, Int_t c, Int_t /*s*/) const
{
  //
  // Returns the position of the border of the first pad in a row
  //

  return GetPadPlane(p,c)->GetRow0();

}

//_____________________________________________________________________________
Double_t AliTRDCommonParam::GetCol0(Int_t p) const
{
  //
  // Returns the position of the border of the first pad in a column
  //

  return GetPadPlane(p,0)->GetCol0();

}
