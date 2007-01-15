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
 
/* $Id$ */

//---------------------------------------------------------------------
//  Jet Gen header class
// Stores parameters of particle algoritm
// Author: Rafael.Diaz.Valdes@cern.ch
//---------------------------------------------------------------------

#include <Riostream.h> 
#include <TMath.h>

#include "AliUA1JetHeaderV1.h"
ClassImp(AliUA1JetHeaderV1)

////////////////////////////////////////////////////////////////////////

AliUA1JetHeaderV1::AliUA1JetHeaderV1():

    AliJetHeader("AliUA1JetHeaderV1"),
    fConeRadius(0.3),
    fEtSeed(3.0),
    fMinJetEt(10.),
    fMinMove(0.05),
    fMaxMove(0.15),
    fBackgMode(1),
    fPrecBg(0.035),
    fBackgStat(0.0), 
    fBackgCutRatio(1.0),
    fNAcceptJets(3),
    fLegoNbinEta(36),
    fLegoNbinPhi(124),
    fLegoEtaMin(-0.9),
    fLegoEtaMax(0.9),
    fLegoPhiMin(0.),
    fLegoPhiMax(2.*TMath::Pi())
{
  // Constructor
}

////////////////////////////////////////////////////////////////////////

void AliUA1JetHeaderV1::PrintParameters() const
{
  // prints out parameters of jet algorithm

  cout << " UA version1 jet algorithm " << endl;
  cout << " * Jet parameters: " << endl;
  cout << "   Cone size: " << fConeRadius<< endl;
  cout << "   Minimum energy for a seed: " << fEtSeed << endl;
  cout << "   Minumum energy for a jet: " << fMinJetEt << endl;
  cout << "   Minimum allowed move: " << fMinMove << endl;
  cout << "   Maximum allowed move: " << fMaxMove << endl;
  cout << " * Lego parameters: " << endl;
  cout << "   Number of bins in eta: " << fLegoNbinEta<< endl;
  cout << "   Number of bins in phi: " << fLegoNbinPhi<< endl;
  cout << "   Minimum azimuthal angle: " << fLegoPhiMin<< endl;
  cout << "   Maximum azimuthal angle: " << fLegoPhiMax<< endl;
  cout << "   Minimum rapidity angle: " << fLegoEtaMin<< endl;
  cout << "   Maximum rapidity angle: " << fLegoEtaMax<< endl;
}
