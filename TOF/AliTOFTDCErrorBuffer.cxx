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
//        This class provides a buffer for TDC errors.              //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliTOFTDCErrorBuffer.h"
#include "AliLog.h"

ClassImp(AliTOFTDCErrorBuffer)

AliTOFTDCErrorBuffer::AliTOFTDCErrorBuffer() :
  TObject(),
  fBuffer("AliTOFTDCError")
{
  /* default constructor */
  fBuffer.SetOwner(kTRUE);
}

//_________________________________________________________________

AliTOFTDCErrorBuffer::~AliTOFTDCErrorBuffer()
{
  /* destructor */
}

//_________________________________________________________________

void
AliTOFTDCErrorBuffer::Add(const AliTOFTDCError &Error)
{
  /* add function */

  new (fBuffer[GetEntries()]) AliTOFTDCError(Error);
}
