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
Revision 1.9  2001/05/16 14:57:22  alibrary
New files for folders and Stack

Revision 1.8  2000/12/18 10:44:01  morsch
Possibility to set field map by passing pointer to objet of type AliMagF via
SetField().
Example:
gAlice->SetField(new AliMagFCM("Map2", "$(ALICE_ROOT)/data/field01.dat",2,1.,10.));

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

ClassImp(AliMagF)

//_______________________________________________________________________
AliMagF::AliMagF():
  fMap(0),
  fType(0),
  fInteg(0),
  fFactor(0),
  fMax(0),
  fDebug(0)
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________
AliMagF::AliMagF(const char *name, const char *title, const Int_t integ, 
                 const Float_t factor, const Float_t fmax):
  TNamed(name,title),
  fMap(0),
  fType(0),
  fInteg(0),
  fFactor(factor),
  fMax(fmax),
  fDebug(0)
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
    //
    fDebug = 0;
}

//_______________________________________________________________________
void AliMagF::Field(Float_t*, Float_t *b)
{
  //
  // Method to return the field in one point -- dummy in this case
  //
  Warning("Field","Undefined MagF Field called, returning 0\n");
  b[0]=b[1]=b[2]=0;
}
