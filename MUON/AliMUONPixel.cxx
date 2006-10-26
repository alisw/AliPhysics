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

// -------------------------------------
// Class AliMUONPixel
// -------------------------------------
// Basic object of the cluster / rec. point finder based 
// on Expectation-Minimization approach (AZ cluster finder)
// Author: Alexander Zinchenko, JINR Dubna

#include "AliMUONPixel.h"

/// \cond CLASSIMP
ClassImp(AliMUONPixel) // Class implementation in ROOT context
/// \endcond

//_____________________________________________________________________________
AliMUONPixel::AliMUONPixel()
  : TObject(),
    fCharge(0),
    fFlag(0)
{
/// Default constructor
  fXY[0] = fXY[1] = fSize[0] = fSize[1] = 0;
} 

//_____________________________________________________________________________
AliMUONPixel::AliMUONPixel(Double_t xc, Double_t yc, Double_t wx, Double_t wy, Double_t charge)
  : TObject(),
    fCharge(charge),
    fFlag(0)
{
/// Constructor
  fXY[0] = xc; fXY[1] = yc; fSize[0] = wx; fSize[1] = wy;
}

//_____________________________________________________________________________
AliMUONPixel::~AliMUONPixel()
{
/// Destructor
}

//__________________________________________________________________________
Int_t AliMUONPixel::Compare(const TObject* pixel) const
{
/// "Compare" function to sort with decreasing pixel charge.
/// Returns -1 (0, +1) if charge of current pixel
/// is greater than (equal to, less than) charge of pixel
  if (fCharge > ((AliMUONPixel*)pixel)->Charge()) return(-1);
  else if (fCharge == ((AliMUONPixel*)pixel)->Charge()) return( 0);
  else return(+1);
}

//__________________________________________________________________________
void AliMUONPixel::Print(const char* /*opt*/) const
{
/// Print function
  printf("%9.4f %9.4f %9.4f %9.4f %9.4f %1d\n", fXY[0], fXY[1], fSize[0], fSize[1], fCharge, fFlag);
}
