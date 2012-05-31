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


#include "AliACORDEdigit.h"

ClassImp(AliACORDEdigit)

//_____________________________________________________________________________
AliACORDEdigit::AliACORDEdigit()
  : AliDigit(),
    fModule(0),
    fTime(0)
{
  //
  // Default constructor
  //
}
//_____________________________________________________________________________
AliACORDEdigit::AliACORDEdigit(Int_t module, Float_t pulse_time)
  : AliDigit(),
    fModule(module),
    fTime(pulse_time)
{
	//
}
//_____________________________________________________________________________
AliACORDEdigit::AliACORDEdigit(Int_t* tracks, Int_t module, Float_t pulse_time)
  : AliDigit(tracks),
    fModule(module),
    fTime(pulse_time)
{
	//
}

//_____________________________________________________________________________
AliACORDEdigit::AliACORDEdigit(Int_t* modules, Float_t pulse_time)
  : AliDigit(modules),
    fModule(0),
    fTime(pulse_time)
{
	//MRC's part
}
//_____________________________________________________________________________
AliACORDEdigit::~AliACORDEdigit()
{
  //
  //
  //
}

//_____________________________________________________________________________
void AliACORDEdigit::Print(const Option_t*) const
{
	Dump();
}

