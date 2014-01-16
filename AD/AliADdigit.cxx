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

#include "AliADdigit.h"

ClassImp(AliADdigit)

//__________________________________________________________________________
AliADdigit::AliADdigit()
   :AliDigit(),
    fModule(0),
    fCell(0)

{
  // Standard default
  // constructor 
}

//__________________________________________________________________________
AliADdigit::AliADdigit(Int_t module, Float_t cellPad)
	: AliDigit(),
	fModule(module),
	fCell(cellPad)
{
}

//__________________________________________________________________________
AliADdigit::AliADdigit(Int_t* tracks, Int_t  module, Float_t cellPad)
	:AliDigit(tracks),
	fModule(module),
	fCell(cellPad)
{
}
//__________________________________________________________________________
AliADdigit::AliADdigit(Int_t* module, Float_t cellPad)
	: AliDigit(module),
	  fModule(0),
	  fCell(cellPad)
{
}
//__________________________________________________________________________
AliADdigit::~AliADdigit()
{
  //
  //
  //
}


//__________________________________________________________________________
void AliADdigit::Print(const Option_t*) const
{
    // Dumps digit object
    Dump();
}
