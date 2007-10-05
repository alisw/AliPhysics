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
//  FastJet header class
// Stores parameters of particle algoritm
// Author: Rafael.Diaz.Valdes@cern.ch
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TMath.h>
#include "AliFastJetHeader.h"

ClassImp(AliFastJetHeader)

////////////////////////////////////////////////////////////////////////

AliFastJetHeader::AliFastJetHeader():
    AliJetHeader("AliFastJetHeader"),
    fLegoNbinEta(60),
    fLegoNbinPhi(210),
    fLegoEtaMin(-0.9),
    fLegoEtaMax(0.9),
    fLegoPhiMin(0.),
    fLegoPhiMax(2. * TMath::Pi()),
    fRadius(1.0),
    fMinJetEt(10.0),
    fSoftBackg(kTRUE),
    fPrecBg(0.035) 
{
  // Constructor
}

////////////////////////////////////////////////////////////////////////

void AliFastJetHeader::PrintParameters() const
{
  // prints out parameters of jet algorithm

  cout << " FastJet algorithm " << endl;

}
