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

const Float_t AliITSresponseSPD::fgkDiffCoeffDefault = 0.;
const Float_t AliITSresponseSPD::fgkThreshDefault = 2000.;
const Float_t AliITSresponseSPD::fgkSigmaDefault = 280.;

ClassImp(AliITSresponseSPD)	
//______________________________________________________________________
AliITSresponseSPD::AliITSresponseSPD(){
  // constructor

   SetThresholds(fgkThreshDefault,fgkSigmaDefault);
   SetDiffCoeff(fgkDiffCoeffDefault,0.);
   SetNoiseParam(0.,0.);
   SetDataType();
   SetFractionDead();
}

