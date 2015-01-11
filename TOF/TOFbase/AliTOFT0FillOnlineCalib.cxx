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
// * this class defines the T0Fill online calibration object to be stored
// * in OCDB in order to obtain a corrected T0Fill from online algorithm
// *
// *
// *

#include "AliTOFT0FillOnlineCalib.h"

ClassImp(AliTOFT0FillOnlineCalib)

//_________________________________________________________

AliTOFT0FillOnlineCalib::AliTOFT0FillOnlineCalib() :
  TObject(),
  fOffset(0.),
  fCoefficient(0.)
{
  /*
   * default constructor
   */
}

//_________________________________________________________

AliTOFT0FillOnlineCalib::~AliTOFT0FillOnlineCalib()
{
  /*
   * default destructor
   */
}

//_________________________________________________________

AliTOFT0FillOnlineCalib::AliTOFT0FillOnlineCalib(const AliTOFT0FillOnlineCalib &source) :
  TObject(source),
  fOffset(source.fOffset),
  fCoefficient(source.fCoefficient)
{
  /*
   * copy constructor
   */
}

//_________________________________________________________

AliTOFT0FillOnlineCalib &
AliTOFT0FillOnlineCalib::operator=(const AliTOFT0FillOnlineCalib &source)
{
  /*
   * operator=
   */

  if (this == &source) return *this;
  TObject::operator=(source);
  fOffset = source.fOffset;
  fCoefficient = source.fCoefficient;
  return *this;
}

