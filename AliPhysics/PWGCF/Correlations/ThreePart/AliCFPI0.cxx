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
/* $Id: $ */
 
//_________________________________________________________________________
// Class to collect neutral pion kinematics and information
//
//-- Author: Paul Batzing

#include "AliCFPI0.h"
ClassImp(AliCFPI0)
//
AliCFPI0::AliCFPI0():AliVParticle()
  , fPx(0.)
  , fPy(0.)
  , fPz(0.)
  , fPt(0.)
  , fP(0.)
  , fPhi(0.)
  , fTheta(0.)
  , fEta(0.)
  , fE(0.)
  , fIsPhos(0)
  , fIsEmcal(0)
{
  
}

AliCFPI0::AliCFPI0(TLorentzVector V):AliVParticle()

  , fPx(0.)
  , fPy(0.)
  , fPz(0.)
  , fPt(0.)
  , fPhi(0.)
  , fTheta(0.)
  , fEta(0.)
  , fE(0.)
  , fIsPhos(0)
  , fIsEmcal(0)
{
  fPx = V.Px();
  fPy = V.Py();
  fPz = V.Pz();
  fPt = V.Pt();
  fP = V.P();
  fE = V.E();
  fTheta = V.Theta();
  fPhi = V.Phi();
  fEta = V.Eta();
}


