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
***************************************************************************/

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/


//////////////////////////////////////////////////////////////////////
//                                                                  //
//                                                                  //
//        This class provides a definition for TDC errors.          //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliTOFTDCError.h"

ClassImp(AliTOFTDCError)

AliTOFTDCError::AliTOFTDCError() :
  TObject(),
  fErrorFlags(0),
  fTDCID(0)
{
  /* default constructor */
}

//_________________________________________________________________

AliTOFTDCError::AliTOFTDCError(const AliTOFTDCError &source) :
  TObject(),
  fErrorFlags(source.fErrorFlags),
  fTDCID(source.fTDCID)
{
  /* copy contructor */
}

//_________________________________________________________________

AliTOFTDCError &
AliTOFTDCError::operator = (const AliTOFTDCError &source)
{
  /* operator = */
  fErrorFlags = source.fErrorFlags;
  fTDCID = source.fTDCID;
  return *this;
}

//_________________________________________________________________

AliTOFTDCError::~AliTOFTDCError()
{
  /* default destructor */
}
