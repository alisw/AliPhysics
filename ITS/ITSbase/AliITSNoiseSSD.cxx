/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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


#include "AliITSNoiseSSD.h"

//////////////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
//                                                  //
//////////////////////////////////////////////////////

ClassImp(AliITSNoiseSSD)

//______________________________________________________________________
AliITSNoiseSSD::AliITSNoiseSSD():
fMod(0),
fNoisP(0),
fNoisN(0) {
    // Default Constructor
}

//______________________________________________________________________
AliITSNoiseSSD::AliITSNoiseSSD(const AliITSNoiseSSD &source): TObject(source),
fMod(source.fMod),
fNoisP(source.fNoisP),
fNoisN(source.fNoisN) {
    // copy Constructor
}
//______________________________________________________________________
AliITSNoiseSSD::~AliITSNoiseSSD(){
    // destructor

}

//______________________________________________________________________
AliITSNoiseSSD& AliITSNoiseSSD::operator=(const AliITSNoiseSSD &source) {
 // ass. op.
    if (this == &source)
      return *this;
    fMod = source.fMod;
    fNoisP =  source.fMod;
    fNoisN =  source.fMod;
    return *this;
}
