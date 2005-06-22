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
// UA1 Jet header class 
// Stores parameters of UA1 jet algo
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <Riostream.h> 
#include <TMath.h>
#include "AliUA1JetHeader.h"
ClassImp(AliUA1JetHeader)
 
 
////////////////////////////////////////////////////////////////////////

AliUA1JetHeader::AliUA1JetHeader():
  AliJetHeader("AliUA1JetHeader")
{
  //
  // Constructor
  //
  fConeRadius =  0.3;
  fEtSeed     =  3.0;
  fMinJetEt   = 10.0;
  fMinCellEt  =  0.0;
  fMode       =  1;
  fMinMove    =  0.05;
  fMaxMove    =  0.15;
  fPrecBg     =  0.035;
  fNbinEta    =  36;
  fNbinPhi    = 124;
  fPhiMin     =   0.;
  fPhiMax     = 2. * TMath::Pi();
  fEtaMin     = -0.9;
  fEtaMax     =  0.9;
}
 

////////////////////////////////////////////////////////////////////////
 
void AliUA1JetHeader::PrintParameters() const

{
  //
  // prints out parameters of jet algorithm
  //

  cout << " UA1 jet algorithm " << endl;
  cout << " * Jet parameters: " << endl;
  cout << "   Cone size: " << fConeRadius<< endl;
  cout << "   Minimum energy for a seed: " << fEtSeed << endl;
  cout << "   Minumum energy for a jet: " << fMinJetEt << endl;
  cout << "   Minumum energy for a cell: " << fMinCellEt << endl;
  cout << " * Background substraction parameters: " << endl;
  cout << "   Substraction mode: " << fMode << endl;
  cout << "   Minimum allowed move: " << fMinMove << endl;
  cout << "   Maximum allowed move: " << fMaxMove << endl;
  cout << "   Precision for background: " << fPrecBg << endl;
  cout << " * Lego parameters: " << endl;
  cout << "   Number of bins in eta: " << fNbinEta<< endl;
  cout << "   Number of bins in phi: " << fNbinPhi<< endl;
  cout << "   Minimum azimuthal angle: " << fPhiMin<< endl;
  cout << "   Maximum azimuthal angle: " << fPhiMax<< endl;
  cout << "   Minimum rapidity angle: " << fEtaMin<< endl;
  cout << "   Maximum rapidity angle: " << fEtaMax<< endl;
}
