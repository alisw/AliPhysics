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

// ---------------------------
// Class AliMUONDigitMapA1
// ---------------------------
// Implements cluster Map as a 1-dim array
// Author: Christian Finck

#include "AliMUONDigitMapA1.h"
#include "AliMUONDigit.h"

#include "AliLog.h"

#include <TObjArray.h>
#include <TMath.h>

/// \cond CLASSIMP
ClassImp(AliMUONDigitMapA1)
/// \endcond

AliMUONDigitMapA1::AliMUONDigitMapA1() 
  : fIdDE(0),
    fNpx(0),
    fNpy(0),
    fDigits(0),
    fMaxIndex(0),
    fHitMap(0)
{
/// Default constructor

}
//________________________________________________________________________________
AliMUONDigitMapA1::AliMUONDigitMapA1(Int_t idDE, Int_t npx, Int_t npy )
  : fIdDE(idDE),
    fNpx(npx),
    fNpy(npy),
    fDigits(0),
    fMaxIndex(fNpx * fNpy + fNpx),
    fHitMap(new Int_t[fMaxIndex])
{
/// Standard constructor
   Clear();
}

//_________________________________
AliMUONDigitMapA1::~AliMUONDigitMapA1()
{
/// Destructor
    if (fHitMap) delete[] fHitMap;
}
//______________________________________
void AliMUONDigitMapA1::Clear(const char *)
{
/// Clear hitmap
    memset(fHitMap,0,sizeof(int)*fMaxIndex);
}

//_________________________________________________________
Int_t AliMUONDigitMapA1::CheckedIndex(Int_t ix, Int_t iy) const
{
/// Return checked indices ix, iy
  Int_t index =  iy * fNpx + ix;
    if ( index < 0 || index >= fMaxIndex ) {
      AliWarning(Form("index outside array idDE %d ix %d iy %d MaxIndex %d index %d Npx %d Npy %d",
		      fIdDE, ix,iy, fMaxIndex, index, fNpx, fNpy));
	return  fMaxIndex-1;
    } else {
	return index;
    }
}
//_____________________________
void  AliMUONDigitMapA1::FillHits(TObjArray* digits)
{
/// Fill hits from digits list  
    fDigits = digits;
    Int_t ndigits = fDigits->GetEntriesFast();
    if (!ndigits) return;
    AliMUONDigit *dig;
    for (Int_t ndig = 0; ndig < ndigits; ndig++) {
	dig = (AliMUONDigit*)fDigits->UncheckedAt(ndig);
	SetHit(dig->PadX(),dig->PadY(),ndig);
    }
}
//___________________________________________________________
void  AliMUONDigitMapA1::SetHit(Int_t ix, Int_t iy, Int_t idigit)
{
/// Assign digit to hit cell ix,iy

    fHitMap[CheckedIndex(ix, iy)] = idigit + 1;
}
//_______________________________________________
void AliMUONDigitMapA1::DeleteHit(Int_t ix, Int_t iy)
{
/// Delete hit at cell ix,iy

    fHitMap[CheckedIndex(ix, iy)]=0;
}
//_____________________________________________
void AliMUONDigitMapA1::FlagHit(Int_t ix, Int_t iy)
{
/// Flag hit as used
    fHitMap[CheckedIndex(ix, iy)]=
	-TMath::Abs(fHitMap[CheckedIndex(ix, iy)]);
}
//________________________________________________________
Int_t AliMUONDigitMapA1::GetHitIndex(Int_t ix, Int_t iy) const
{
/// Get absolute value of contents of hit cell ix,iy
    return TMath::Abs(fHitMap[CheckedIndex(ix, iy)])-1;
}
//_______________________________________________________
TObject* AliMUONDigitMapA1::GetHit(Int_t ix, Int_t iy) const
{
/// Get pointer to object at hit cell ix, iy
/// Force crash if index does not exist ! (Manu)
    Int_t index = GetHitIndex(ix,iy);
    return (index <0) ? 0 : fDigits->UncheckedAt(GetHitIndex(ix,iy));
}
//_________________________________________________
FlagType AliMUONDigitMapA1::TestHit(Int_t ix, Int_t iy) const
{
/// Check if hit cell is empty, used or unused

    Int_t index = CheckedIndex(ix, iy);
    if (index<0 || index >= fMaxIndex) return kEmpty;

    Int_t inf = fHitMap[index];
    if (inf < 0) {
	return kUsed;
    } else if (inf == 0) {
	return kEmpty;
    } else {
	return kUnused;
    }
}
