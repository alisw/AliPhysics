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

#include "AliITSgeomSPD.h"

ClassImp(AliITSgeomSPD)
AliITSgeomSPD::AliITSgeomSPD(){
    //
    // default constructor
    //
    fShapeSPD = new TBRIK("ActiveSPD","Active volume of SPD","SPD SI CHIP",
			  2.5E-2/2.0,1.38/2.0,8.2/2.0);
}
