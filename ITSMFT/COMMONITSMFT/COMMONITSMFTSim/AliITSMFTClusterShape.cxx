/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

#include <iostream>
#include <TBits.h>
#include <TRandom.h>

#include "AliITSMFTClusterShape.h"

ClassImp(AliITSMFTClusterShape)

//______________________________________________________________________
AliITSMFTClusterShape::AliITSMFTClusterShape() :
fNrows(0),
fNcols(0),
fNFPix(0),
fShape(0)
{}


//______________________________________________________________________
AliITSMFTClusterShape::AliITSMFTClusterShape(UInt_t Nrows, UInt_t Ncols, UInt_t NFPix) {
  fNrows = Nrows;
  fNcols = Ncols;
  fNFPix = NFPix;
  fShape = new UInt_t[fNFPix];
}


//______________________________________________________________________
AliITSMFTClusterShape::~AliITSMFTClusterShape() {}


//______________________________________________________________________
Long64_t AliITSMFTClusterShape::GetShapeID() {
  // DJBX33X
  Long64_t id = 5381;
  for (UInt_t i = 0; i < fNFPix; ++i) {
    id = ((id << 5) + id) ^ fShape[i];
  }
  return id;
}
