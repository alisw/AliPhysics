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
// Base class for reconstructed charged particle multiplicty in the
// FMD.  The class contains the field fMethod which is a flag set by
// the AliFMDMultAlgorithm that created the object. The flag tells us
// which algorithm was used to create the data stored in the object. 
//
#include "AliFMDMult.h"		// ALIFMDMULT_H
#include <TString.h>            // ROOT_TString 
#include <Riostream.h>		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDMult)

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
