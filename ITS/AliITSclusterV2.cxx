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

//-------------------------------------------------------------------------
//              Implementation of the ITS cluster class
//
//         Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "AliITSclusterV2.h"

ClassImp(AliITSclusterV2)
//_______________________________________________________
AliITSclusterV2::AliITSclusterV2() : AliCluster(),
fIndex(0),
fQ(0),
fLayer(0),
fNz(0),
fNy(0),
fChargeRatio(0),
fType(0),
fDeltaProb(0) {
  //default constructor
}

//_______________________________________________________
AliITSclusterV2::AliITSclusterV2(Int_t *lab,Float_t *hit, Int_t *info) : AliCluster(lab,hit),
fIndex(lab[3]),
fQ(hit[4]),
fLayer(info[2]),
fNz(info[1]),
fNy(info[0]),
fChargeRatio(0),
fType(0),
fDeltaProb(0){
  //standard constructor
}
