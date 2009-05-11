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

// $Id$

#include "AliMUONErrorCounter.h"

#include <Riostream.h>

//-----------------------------------------------------------------------------
/// \class AliMUONErrorCounter
/// 
/// add 
/// 
///
/// \author Alberto Baldisseri, JL Charvet (05/05/2009)
//-----------------------------------------------------------------------------

//______________________________________________________________________________
AliMUONErrorCounter::AliMUONErrorCounter(Int_t bp, Int_t manu, Int_t ev) 
  : fBusPatch(bp), 
    fManuId(manu), 
    fEvents(ev) 
{
}

//______________________________________________________________________________
Int_t AliMUONErrorCounter::Compare(const TObject* obj) const
{
  /// Compare function
  Int_t patch1, patch2, manu1, manu2;
  patch1 = fBusPatch;
  manu1 = fManuId;
  patch2 = ((AliMUONErrorCounter*)obj)->BusPatch();
  manu2 = ((AliMUONErrorCounter*)obj)->ManuId();

  if (patch1 == patch2)
  {
    if (manu1 == manu2)
    {
      return 0;
    }
    else
      return (manu1 >= manu2) ? 1 : -1;
  }
  else
    return (patch1 >= patch2) ? 1 : -1;
}

//______________________________________________________________________________
void AliMUONErrorCounter::Print(const Option_t* option) const
{
  TNamed::Print(option);
  cout<<"bp "<<fBusPatch<<" events "<<fEvents<<endl;
}

//______________________________________________________________________________
void AliMUONErrorCounter::Print_uncal(const Option_t* option) const
{
  TNamed::Print(option);
  cout<<"bp =  "<<fBusPatch<< "  manu = " << fManuId << " uncal = "<< fEvents <<endl;
}
