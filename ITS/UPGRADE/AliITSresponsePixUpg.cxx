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
///////////////////////////////////////////////////////////////////////////
//  Base Response for pixels                      
//  It is used to set static data members           
//  connected to parameters equal for all           
//  the SPD modules                                 
//
//  Modified by D. Elia, G.E. Bruno
//  March-April 2006
//  September   2007: Coupling params taken out
//                    left in AliITSCalibrationSPD only
//
///////////////////////////////////////////////////////////////////////////

#include "AliITSresponsePixUpg.h"

const Float_t AliITSresponsePixUpg::fgkDiffCoeffDefault = 0.; //change this
const TString AliITSresponsePixUpg::fgkCouplingOptDefault = "old";
const Float_t AliITSresponsePixUpg::fgkEccentricityDiffDefault = 0.85;

ClassImp(AliITSresponsePixUpg)	
//______________________________________________________________________
AliITSresponsePixUpg::AliITSresponsePixUpg():
  AliITSresponse(),
fCouplOpt(0),
fEccDiff(0){

  // constructor
  SetCouplingOption(fgkCouplingOptDefault);
  SetDiffCoeff(fgkDiffCoeffDefault,0.);
  SetSigmaDiffusionAsymmetry(fgkEccentricityDiffDefault);

}
