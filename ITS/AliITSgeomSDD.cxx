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
*/

#include "AliITSgeomSDD.h"

ClassImp(AliITSgeomSDD)
AliITSgeomSDD::AliITSgeomSDD(){
    //
    // default constructor
    //
    fShapeSDD = new TBRIK("ActiveSDD","Active volume of SDD","SDD SI CHIP",
			    3.0E-2/2.,7.25/2.,7.53/2.);
}
