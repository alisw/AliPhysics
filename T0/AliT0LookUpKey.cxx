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
/*************************************************************************
 * class for logical adress of chanells:
 * can be call by nubmer or by name
 *
 *  Alla.Maevskaya@cern.ch
 ************************************************************************/ 

#include "AliT0LookUpKey.h"

ClassImp(AliT0LookUpKey)

 AliT0LookUpKey::AliT0LookUpKey():
   TObject(),
   fKey(0),
   fName("")
 {
 
 }

//--------------------------------------------------------------------------------------
 AliT0LookUpKey::AliT0LookUpKey(Int_t key):
   TObject(),
   fKey(key),
   fName("")
 {
 
 }

//--------------------------------------------------------------------------------------
 AliT0LookUpKey::AliT0LookUpKey(TString name):
   TObject(),
   fKey(),
   fName(name)
 {
 
 }
//________________________________________________________________
AliT0LookUpKey::AliT0LookUpKey(const AliT0LookUpKey& o) :		
   TObject(),
   fKey(),
   fName("")

{
// copy constructor
 ((AliT0LookUpKey&) o).Copy(*this);
}
//--------------------------------------------------------------------------------------

Bool_t AliT0LookUpKey:: IsEqual(const TObject* obj) const
{
  Int_t k;
//  printf("IsEqual\n");
    k=((AliT0LookUpKey*)obj)->GetKey();
    if (k==fKey) return kTRUE; else return kFALSE;
}
//--------------------------------------------------------------------------------------
void AliT0LookUpKey::Print(Option_t *) const
{
  printf(" AliT0LookUpKey key %i name %s\n",fKey,fName.Data());

}
