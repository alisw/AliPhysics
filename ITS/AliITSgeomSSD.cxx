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
Revision 1.2.4.2  2000/03/04 23:55:59  nilsen
Fixed up the comments/documentation

Revision 1.2.4.1  2000/01/12 19:03:32  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

Revision 1.2  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

#include "AliITSgeomSSD.h"

ClassImp(AliITSgeomSSD)
AliITSgeomSSD::AliITSgeomSSD(){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////

    Float_t dx = 3.65;  // cm. (Geant 3.12 units) Orthonormal to y and z
    Float_t dy = 0.015; // cm. (Geant 3.12 units) Radialy from the Beam Pipe
    Float_t dz = 2.0;   // cm. (Geant 3.12 units) Allong the Beam Pipe

    fShapeSSD = new TBRIK("ActiveSSD","Active volume of SSD","SSD SI DET",
		  	    dx,dy,dz);
}
