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
//        This class provides a buffer for TDC hits.                //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliTOFTDCHitBuffer.h"
#include "AliLog.h"

ClassImp(AliTOFTDCHitBuffer)

AliTOFTDCHitBuffer::AliTOFTDCHitBuffer() :
  TObject(),
  fBuffer("AliTOFTDCHit")
{
  /* default constructor */
  fBuffer.SetOwner(kTRUE);
}

//_________________________________________________________________

AliTOFTDCHitBuffer::~AliTOFTDCHitBuffer()
{
  /* destructr */
}

//_________________________________________________________________

void
AliTOFTDCHitBuffer::Add(const AliTOFTDCHit &Hit)
{
  /* add function */

  new (fBuffer[GetEntries()]) AliTOFTDCHit(Hit);
}
