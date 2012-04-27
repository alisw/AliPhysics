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

// *
// *
// *
// * this class defines the TOF object to be stored
// * in OCDB on a run-by-run basis in order to have the status
// * of TOF trigger inputs. it stores 32 bit masks for each crate
// * 
// *
// *
// *

#include "AliTOFTriggerMask.h"

ClassImp(AliTOFTriggerMask)

//_________________________________________________________

AliTOFTriggerMask::AliTOFTriggerMask() :
  TObject(),
  fTriggerMask()
{
  /*
   * default constructor
   */

  for (Int_t iddl = 0; iddl < 72; iddl++) fTriggerMask[iddl] = 0;
}

//_________________________________________________________

AliTOFTriggerMask::~AliTOFTriggerMask()
{
  /*
   * default destructor
   */

}

//_________________________________________________________

AliTOFTriggerMask::AliTOFTriggerMask(const AliTOFTriggerMask &source) :
  TObject(source),
  fTriggerMask()
{
  /*
   * copy constructor
   */

  for (Int_t iddl = 0; iddl < 72; iddl++) fTriggerMask[iddl] = source.fTriggerMask[iddl];
}

//_________________________________________________________

AliTOFTriggerMask &
AliTOFTriggerMask::operator=(const AliTOFTriggerMask &source)
{
  /*
   * operator=
   */

  if (this == &source) return *this;
  TObject::operator=(source);
  
  for (Int_t iddl = 0; iddl < 72; iddl++) fTriggerMask[iddl] = source.fTriggerMask[iddl];

  return *this;
}

//_________________________________________________________

void
AliTOFTriggerMask::SetTriggerMaskArray(UInt_t *array)
{
  /*
   * set trigger mask array
   */

  for (Int_t iddl = 0; iddl < 72; iddl++) fTriggerMask[iddl] = array[iddl];
}
//_________________________________________________________

Int_t AliTOFTriggerMask::GetNumberMaxiPadOn() {
  Int_t n=0;
  for(Int_t j=0;j<72;j++) 
    for(Int_t i=23;i>0;i--) 
      n += (fTriggerMask[j]%Int_t(TMath::Power(2.,i+1.)))/Int_t(TMath::Power(2.,i+0.));
  return n;
};
