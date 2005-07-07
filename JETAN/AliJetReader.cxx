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
 
//---------------------------------------------------------------------
// Jet reader base class
// manages the reading of input for jet algorithms
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <TClonesArray.h>

#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliESD.h"
#include "AliHeader.h"

ClassImp(AliJetReader)

////////////////////////////////////////////////////////////////////////

AliJetReader::AliJetReader()
{
  // Constructor
  fChain = 0; 
  fChainMC = 0;
  fESD = 0;
  fMomentumArray = new TClonesArray("TLorentzVector",2000);
  fArrayMC = 0;
  fAliHeader = 0;
  fSignalFlag = TArrayI();
}

////////////////////////////////////////////////////////////////////////

AliJetReader::~AliJetReader()
{
  // Destructor
  delete fChain;
  delete fChainMC;
  delete fESD;
  delete fAliHeader;
  if (fMomentumArray) {
      fMomentumArray->Delete();
      delete fMomentumArray;
  }
  delete fArrayMC;
}


////////////////////////////////////////////////////////////////////////

void AliJetReader::ClearArray()

{
  if (fMomentumArray)  fMomentumArray->Clear();
}
