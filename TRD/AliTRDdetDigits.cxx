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

/*
$Log$
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Container class for the TRD digits of one detector segment (chamber).    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDdetDigits.h"

ClassImp(AliTRDdetDigits)

//_____________________________________________________________________________
AliTRDdetDigits::AliTRDdetDigits():AliDigits()
{
  //
  // Default constructor
  //

  fNrowTRD = 0;
  fNcolTRD = 0;

}

//_____________________________________________________________________________
AliTRDdetDigits::Allocate(Int_t nrow, Int_t ncol, Int_t ntime)
{
  //
  // Allocates an empty buffer for the digits data of the size
  // <nrow> * <ncol> * <ntime>
  //

  fNrowTRD = nrow;
  fNcolTRD = ncol;
  
  // The two-dimensional row/column structure of the TRD gets mapped into 
  // an one-dimensional row array which can be used inside AliDigits.
  // The TRD timebins correspond to the columns in AliDigits.
  AliDigits::Allocate(fNrowTRD*fNcolTRD,ntime);

}

//_____________________________________________________________________________
AliTRDdetDigits::SetDigits(Int_t row, Int_t col, Int_t time, Short_t value)
{
  //
  // Sets the value of one given digit
  //

  Int_t index = row * fNrowTRD + col;
  AliDigits::SetDigitFast(value,index,time);

}
