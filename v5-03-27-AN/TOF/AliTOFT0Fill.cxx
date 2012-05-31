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
// * this class defines the T0Fill object to be stored
// * in OCDB in order to apply T0Fill correction during 
// * reconstruction. 
// *
// *
// *

#include "AliTOFT0Fill.h"

ClassImp(AliTOFT0Fill)

//_________________________________________________________

AliTOFT0Fill::AliTOFT0Fill() :
  TObject(),
  fT0Fill(0.)
{
  /*
   * default constructor
   */
}

//_________________________________________________________

AliTOFT0Fill::~AliTOFT0Fill()
{
  /*
   * default destructor
   */
}

//_________________________________________________________

AliTOFT0Fill::AliTOFT0Fill(const AliTOFT0Fill &source) :
  TObject(source),
  fT0Fill(source.fT0Fill)
{
  /*
   * copy constructor
   */
}

//_________________________________________________________

AliTOFT0Fill &
AliTOFT0Fill::operator=(const AliTOFT0Fill &source)
{
  /*
   * operator=
   */

  if (this == &source) return *this;
  TObject::operator=(source);
  fT0Fill = source.fT0Fill;
  return *this;
}

