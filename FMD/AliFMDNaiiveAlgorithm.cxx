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

//____________________________________________________________________
// 
// Base class for FMD naiive algorithms. 
//
// Derived classes will implement various ways of reconstructing the
// charge particle multiplicity in the FMD.  
// 
#include "AliFMDNaiiveAlgorithm.h"	// ALIFMDNAIIVEALGORITHM_H
#include "AliFMDDigit.h"		// ALIFMDDIGIT_H

//____________________________________________________________________
ClassImp(AliFMDNaiiveAlgorithm);

//____________________________________________________________________
AliFMDNaiiveAlgorithm::AliFMDNaiiveAlgorithm()
  : AliFMDReconstructionAlgorithm("Naiive", "Naiive")
{}


//____________________________________________________________________
void
AliFMDNaiiveAlgorithm::Reset() 
{
  // Reset internal data 
  fTotal = 0;
  fEdep.Clear(0);
  fMultiplicity.Clear(0);
}

//____________________________________________________________________
void
AliFMDNaiiveAlgorithm::ProcessDigit(AliFMDDigit* digit, 
				    Float_t eta, 
				    Float_t phi, 
				    UShort_t count)
{
  // Process one digit. 
  // 
  // Parameters: 
  //    
  //   digit		Digit to process 
  //   eta              Pseudo-rapidity of digit
  //   phi              Azimuthal angle of digit
  //   count            ADC (corrected for the pedestal)
  //
  // This calculates the energy deposited and the number of MIPs that
  // this energy deposition corresponds to 
  // 
  //   EnergyDeposited = cos(theta) * gain * count 
  //   Multiplicity    = EnergyDeposited / EnergyDepositedPerMIP
  // 
  // where gain is a conversion factor from number of counts to an
  // energy:
  //          Pre_Amp_MIP_Range           1
  //   gain = ----------------- * ---------------------
  //          ADC_channel_size    EnergyDepositedPerMip
  // 
  // and theta is the particles incident angle on the strip, given by 
  //
  //   theta = 2 * atan  * exp(-eta)
  //
  // The cos(theta) factor corrects for the fact that the particle may
  // traverse the strip at an angle, and therefor have a longer flight
  // length, leading to a larger energy deposition. 
  // 
  if (!digit) return;
  Double_t theta = 2 * TMath::Tan(TMath::Exp(- eta));
  Double_t edep  = TMath::Cos(theta) * fGain * count;
  Double_t mult  = edep / fEdepMip;
  fEdep(digit->Detector() - 1, digit->Ring(), 
	digit->Setctor(), digit->Strip()) = edep;
  fMultiplicity(digit->Detector() - 1, digit->Ring(), 
		digit->Setctor(), digit->Strip()) = mult;
}


//____________________________________________________________________
// 
// EOF
//
