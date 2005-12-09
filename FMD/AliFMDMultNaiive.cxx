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
// Reconstruct charged particle multiplicity in the FMD 
// 
// [See also the AliFMDReconstructor class]
// 
// This class reconstructs the multiplicity based on the assumption
// that all particles are minimum ionizing.   
// Hence, the name `naiive'
// 
#include "AliFMD.h"			// ALIFMD_H
#include "AliFMDMultNaiive.h"		// ALIFMDMULTNAIIVE_H
#include "AliFMDParameters.h"           // ALIFMDPARAMETERS_H
#include "AliFMDMultStrip.h"		// ALIFMDMULTNAIIVE_H
#include "AliFMDDigit.h"		// ALIFMDDIGIT_H
#include <TClonesArray.h>               // ROOT_TClonesArray
#include <TTree.h>               	// ROOT_TTree

//____________________________________________________________________
ClassImp(AliFMDMultNaiive)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDMultNaiive::AliFMDMultNaiive()
  : AliFMDMultAlgorithm("Naiive", "Naiive")
{
  // Default CTOR
  fMult = new TClonesArray("AliFMDMultStrip", 1000);
}

//____________________________________________________________________
void
AliFMDMultNaiive::PreRun(AliFMD* fmd) 
{
  // Initialise before a run 
  AliFMDMultAlgorithm::PreRun(fmd);
  AliFMDParameters* pars = AliFMDParameters::Instance();
  fEdepMip = pars->GetEdepMip();
  fGain = (Float_t(pars->GetVA1MipRange()) / pars->GetAltroChannelSize() 
	   * fEdepMip);
}

//____________________________________________________________________
void
AliFMDMultNaiive::PreEvent(TTree* treeR, Float_t ipZ) 
{
  // Reset internal data 
  AliFMDMultAlgorithm::PreEvent(treeR, ipZ);  
  fTreeR->Branch("FMDNaiive", &fMult);
}

//____________________________________________________________________
void
AliFMDMultNaiive::ProcessDigit(AliFMDDigit* digit, 
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
  //   theta = 2 * atan(exp(-eta))
  //
  // The cos(theta) factor corrects for the fact that the particle may
  // traverse the strip at an angle, and therefor have a longer flight
  // length, leading to a larger energy deposition. 
  // 
  if (!digit) return;
  Double_t edep  = Adc2Energy(digit, eta, count);
  Double_t mult  = Energy2Multiplicity(digit, edep);
  
  AliFMDMultStrip* m = 
    new ((*fMult)[fNMult]) AliFMDMultStrip(digit->Detector(), 
					   digit->Ring(), 
					   digit->Sector(),
					   digit->Strip(),
					   eta, phi, 
					   edep, mult,
				   AliFMDMult::kNaiive);
  (void)m;
  fNMult++;
}
//____________________________________________________________________
Float_t
AliFMDMultNaiive::Adc2Energy(AliFMDDigit* /* digit */, 
			     Float_t      eta, 
			     UShort_t     count) 
{
  // Converts number of ADC counts to energy deposited. 
  // Note, that this member function can be overloaded by derived
  // classes to do strip-specific look-ups in databases or the like,
  // to find the proper gain for a strip. 
  // 
  // In this simple version, we calculate the energy deposited as 
  // 
  //    EnergyDeposited = cos(theta) * gain * count
  // 
  // where 
  // 
  //           Pre_amp_MIP_Range
  //    gain = ----------------- * Energy_deposited_per_MIP
  //           ADC_channel_size    
  // 
  // is constant and the same for all strips. 
  Double_t theta = 2 * TMath::ATan(TMath::Exp(-eta));
  Double_t edep  = TMath::Abs(TMath::Cos(theta)) * fGain * count;
  return edep;
}

//____________________________________________________________________
Float_t
AliFMDMultNaiive::Energy2Multiplicity(AliFMDDigit* /* digit */, 
				      Float_t      edep)
{
  // Converts an energy signal to number of particles. 
  // Note, that this member function can be overloaded by derived
  // classes to do strip-specific look-ups in databases or the like,
  // to find the proper gain for a strip. 
  // 
  // In this simple version, we calculate the multiplicity as 
  // 
  //   multiplicity = Energy_deposited / Energy_deposited_per_MIP
  // 
  // where 
  //
  //   Energy_deposited_per_MIP = 1.664 * SI_density * SI_thickness 
  // 
  // is constant and the same for all strips 
  return edep / fEdepMip;
}



//____________________________________________________________________
// 
// EOF
//
