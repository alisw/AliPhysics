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
//////////////////////////////////////////////////////
//  Response class for set:ITS                      //
//  Specific subdetector implementation for         //
//  Silicon pixels                                  //
//  An alternative version "SPDdubna"               //
//  is also available                               //
//////////////////////////////////////////////////////

#include "AliITSresponseSPD.h"

const Double_t AliITSresponseSPD::fgkDiffCoeffDefault = 0.;
const Double_t AliITSresponseSPD::fgkThreshDefault = 2000.;
const Double_t AliITSresponseSPD::fgkSigmaDefault = 280.;

ClassImp(AliITSresponseSPD)	
//______________________________________________________________________
AliITSresponseSPD::AliITSresponseSPD():
AliITSresponse(),
fBaseline(0.0),
fNoise(0.0),
fThresh(fgkThreshDefault),
fSigma(fgkSigmaDefault),
fCouplCol(0.0),
fCouplRow(0.0),
fDeadPixels(0.01){
  // constructor

   SetThresholds(fgkThreshDefault,fgkSigmaDefault);
   //SetDiffCoeff(fgkDiffCoeffDefault,0.);
   SetNoiseParam(0.,0.);
   SetDataType("simulated");
   SetFractionDead();
}
//_________________________________________________________________________
Bool_t AliITSresponseSPD::IsPixelDead(Int_t mod,Int_t ix,Int_t iz) const {
  // Returns kTRUE if pixel is dead
  // Inputs:
  //    Int_t mod      module number
  //    Int_t ix       x pixel number
  //    Int_t iz       z pixel number
  // Outputs:
  //    none.
  // Return:
  //    kFALSE if pixel is alive, or kTRUE if pixel is dead.
  Bool_t  dead = kFALSE;
  Int_t   seed;
  static TRandom ran; // don't use gRandom. This must not be a true randome
  // sequence. These sequence must be random one and then fully repetable.

  seed = mod*256*256+iz*256+ix;
  ran.SetSeed(seed);
  if(ran.Rndm(0)<fDeadPixels) dead = kTRUE;
  return dead;
}
