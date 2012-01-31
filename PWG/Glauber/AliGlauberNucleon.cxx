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

////////////////////////////////////////////////////////////////////////////////
//
//  AliGlauberNucleon implementation
//  support class for Glauber MC
//
//  origin: PHOBOS experiment
//  alification: Mikolaj Krzewicki, Nikhef, mikolaj.krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include "AliGlauberNucleon.h"

ClassImp(AliGlauberNucleon)

//______________________________________________________________________________
AliGlauberNucleon::AliGlauberNucleon():
  fX(0),
  fY(0),
  fZ(0),
  fInNucleusA(0),
  fNColl(0)
{
  //def ctor
}


