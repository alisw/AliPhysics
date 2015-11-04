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

// Header for EPOS v3.111 generated event.
//
// Author: Natalia Zhigareva <Natalia.Zhigareva@cern.ch>

#include "AliGenEpos3EventHeader.h" 

ClassImp(AliGenEpos3EventHeader)

AliGenEpos3EventHeader::AliGenEpos3EventHeader():
  fIversn(0),
  fLaproj(0),
  fMaproj(0),
  fLatarg(0),
  fMatarg(0),
  fEngy(0.),
  fNfull(0),
  fNfreeze(0),
  fBim(0)

{
  // default constructor
}

AliGenEpos3EventHeader::AliGenEpos3EventHeader(const char*name):
  AliGenEventHeader(name),
  fIversn(0),
  fLaproj(0),
  fMaproj(0),
  fLatarg(0),
  fMatarg(0),
  fEngy(0.),
  fNfull(0),
  fNfreeze(0),
  fBim(0)
{
    //constructor
}
// end
