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
Revision 1.7  2000/10/02 21:28:14  fca
Removal of useless dependecies via forward declarations

Revision 1.6  2000/07/11 18:24:59  fca
Coding convention corrections + few minor bug fixes

Revision 1.5  2000/06/09 19:55:18  morsch
Introduce new class AliMagFDM - Galina Chabratova


Revision 1.4  2000/03/28 12:40:24  fca
Introduce factor for magnetic field


Revision 1.3  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/


#include "AliMagF.h"

#include <stdlib.h>
#include <stdio.h>


ClassImp(AliMagF)

//________________________________________
AliMagF::AliMagF(const char *name, const char *title, const Int_t integ, 
		 const Float_t factor, const Float_t fmax)
  : TNamed(name,title)
{
  //
  // Standard constructor
  //
    if(integ<0 || integ > 2) {
	Warning("SetField",
		"Invalid magnetic field flag: %5d; Helix tracking chosen instead\n"
		,integ);
	fInteg = 2;
    } else {
	fInteg = integ;
    }
    fType = kUndef;
    fFactor = factor;
    fMax = fmax;
}

//________________________________________
void AliMagF::Field(Float_t*, Float_t *b)
{
  //
  // Method to return the field in one point -- dummy in this case
  //
  printf("Undefined MagF Field called, returning 0\n");
  b[0]=b[1]=b[2]=0;
}
      


  




