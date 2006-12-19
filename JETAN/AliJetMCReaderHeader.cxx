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
 
// Jet MC Reader Header
// Header for the MC reader in the jet analysis
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)

 
#include "AliJetMCReaderHeader.h"

ClassImp(AliJetMCReaderHeader)

//____________________________________________________________________________

AliJetMCReaderHeader::AliJetMCReaderHeader():
    AliJetESDReaderHeader()
{
  // Constructor
    fDCA = 3.;
}

//____________________________________________________________________________
AliJetMCReaderHeader::~AliJetMCReaderHeader()
{
  // destructor
}
