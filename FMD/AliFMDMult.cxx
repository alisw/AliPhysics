/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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
//  Forward Multiplicity Detector have to be reconstructed number of
//  particles in fixed pseudorapidity interval from fNumOfMinRing
//  to fNumOfMaxRing and phi interval from fNumOfMinSector to
//  fNumOfMaxSector
//
#include "AliFMDMult.h"		// ALIFMDMULT_H
#include <TString.h>            // ROOT_TString 
#include <Riostream.h>		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDMult);

//____________________________________________________________________
AliFMDMult::AliFMDMult(Float_t  particles, UShort_t method)
  : fParticles(particles),
    fMethod(method)
{
  // Constructor
  switch (fMethod) {
  case kPoission: 
  case kIterative: 
  case kNaiive:    
    break;    
  default:        
    Warning("AliFMDMult", "unknown method: %d", method);
    break;
  }
}


//____________________________________________________________________
void
AliFMDMult::Print(Option_t* option) const
{
  // Print information 
  //
  // Options:
  //
  //    V      Be verbose
  // 
  TString opt(option);
  if (!opt.Contains("v", TString::kIgnoreCase)) return;
  cout << "  Method:        " << flush;
  switch (fMethod) {
  case kPoission:  cout << "Poission"  << endl; break;
  case kIterative: cout << "Iterative" << endl; break;
  case kNaiive:    cout << "Naive"     << endl; break; 
  default:         cout << "Unknown"   << endl; break;
  }
}

    
//____________________________________________________________________
//
// EOF
//
