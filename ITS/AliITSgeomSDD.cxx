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
Revision 1.2.4.2  2000/03/04 23:55:35  nilsen
Fixed up comments/documentation.

Revision 1.2.4.1  2000/01/12 19:03:32  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

Revision 1.2  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

#include "AliITSgeomSDD.h"

ClassImp(AliITSgeomSDD)
AliITSgeomSDD::AliITSgeomSDD(){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////

    Float_t fDx = 3.5;  // cm. (Geant 3.12 units) Orthonormal to y and z
    Float_t fDy = 0.014;  // cm. (Geant 3.12 units) Radialy from the Beam Pipe
    Float_t fDz = 3.763;  // cm. (Geant 3.12 units) Allong the Beam Pipe

    fShapeSDD = new TBRIK("ActiveSDD","Active volume of SDD","SDD SI DET",
			    fDx,fDy,fDz);
}
//________________________________________________________________________
AliITSgeomSDD::AliITSgeomSDD(AliITSgeomSDD &source){
  // Copy constructor
   if(this==&source) return;
   this->fShapeSDD = new TBRIK(*(source.fShapeSDD));
   this->fDx = source.fDx;
   this->fDy = source.fDy;
   this->fDz = source.fDz;
}
//________________________________________________________________________
AliITSgeomSDD& AliITSgeomSDD::operator=(AliITSgeomSDD &source){
  // = operator
   if(this==&source) return *this;
   this->fShapeSDD = new TBRIK(*(source.fShapeSDD));
   this->fDx = source.fDx;
   this->fDy = source.fDy;
   this->fDz = source.fDz;
   return *this;
}
