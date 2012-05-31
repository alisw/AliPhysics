/**************************************************************************
 * Copyright(c) 2007-9, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

#include "AliITSPedestalSSD.h"

//////////////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 12/12/2007
//                                                  //
//////////////////////////////////////////////////////

ClassImp(AliITSPedestalSSD)

//______________________________________________________________________
AliITSPedestalSSD::AliITSPedestalSSD():
fMod(0),
fPedP(0),
fPedN(0) {
    // Default Constructor
}

//______________________________________________________________________
AliITSPedestalSSD::AliITSPedestalSSD(const AliITSPedestalSSD &source): TObject(source),
fMod(source.fMod),
fPedP(source.fPedP),
fPedN(source.fPedN) {
    // copy Constructor
}
//______________________________________________________________________
AliITSPedestalSSD::~AliITSPedestalSSD(){
    // destructor

}

//______________________________________________________________________
AliITSPedestalSSD& AliITSPedestalSSD::operator=(const AliITSPedestalSSD &source) {
 // ass. op.
    if (this == &source)
      return *this;
    fMod = source.fMod;
    fPedP =  source.fMod;
    fPedN =  source.fMod;
    return *this;
}
