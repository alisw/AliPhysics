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
#include "AliITSMFTSimuClusterShaper.h"

ClassImp(AliITSMFTSimuClusterShaper)

//______________________________________________________________________
AliITSMFTSimuClusterShaper::AliITSMFTSimuClusterShaper() :
fNrows(0),
fNcols(0),
fNpixOn(0)
{  
}

//______________________________________________________________________
AliITSMFTSimuClusterShaper::AliITSMFTSimuClusterShaper(const Int_t &cs) :
fNrows(0),
fNcols(0),
fNpixOn(cs)
{
  while(fNrows*fNcols < fNpixOn) {
    fNrows += 1;
    fNcols += 1;
  }
}

//______________________________________________________________________
AliITSMFTSimuClusterShaper::~AliITSMFTSimuClusterShaper()
{
    // destructor
}

//______________________________________________________________________
void AliITSMFTSimuClusterShaper::FillClusterRandomly(Int_t *clusterConf)
{
  Int_t matrixSize = fNrows*fNcols;

  // generate UNIQUE random numbers
  Int_t i = 0;
  TBits *bits = new TBits(fNpixOn);
  while (i < fNpixOn) {
    Int_t j = gRandom->Integer(matrixSize); // [0, matrixSize-1]
    if (bits->TestBitNumber(j)) continue; 
    bits->SetBitNumber(j);
    i++;
  }

  Int_t bit = 0;
  for (i = 0; i < fNpixOn; ++i) {
    Int_t j = bits->FirstSetBit(bit);
    clusterConf[i] = j;
    bit = j+1;
  }
}


//______________________________________________________________________ 
void AliITSMFTSimuClusterShaper::AddNoisePixel() 
{
  
}
