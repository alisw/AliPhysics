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
/* $Id: AliFMDDigitizer.cxx 22496 2007-11-26 13:50:44Z cholm $ */
/** @file    AliFMDDigitizer.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:38:26 2006
    @brief   FMD Digitizers implementation
    @ingroup FMD_sim
*/
//////////////////////////////////////////////////////////////////////////////
//
//  This class contains the procedures simulation ADC  signal for the
//  Forward Multiplicity detector  : SDigits->Digits
// 
//  Digits consists of
//   - Detector #
//   - Ring ID                                             
//   - Sector #     
//   - Strip #
//   - ADC count in this channel                                  
//
//  Digits consists of
//   - Detector #
//   - Ring ID                                             
//   - Sector #     
//   - Strip #
//   - Total energy deposited in the strip
//   - ADC count in this channel                                  
//
// As the Digits and SDigits have so much in common, the classes
// AliFMDDigitizer and AliFMDSDigitizer are implemented via a base
// class AliFMDBaseDigitizer.
//
//                 +---------------------+
//                 | AliFMDBaseDigitizer |
//                 +---------------------+
//                           ^
//                           |
//                +----------+---------+
//                |                    |
//      +-----------------+     +------------------+
//      | AliFMDDigitizer |	| AliFMDSDigitizer |
//      +-----------------+	+------------------+
//                |
//     +-------------------+
//     | AliFMDSSDigitizer |
//     +-------------------+
//
// These classes has several paramters: 
//
//     fPedestal
//     fPedestalWidth
//         (Only AliFMDDigitizer)
//         Mean and width of the pedestal.  The pedestal is simulated
//         by a Guassian, but derived classes my override MakePedestal
//         to simulate it differently (or pick it up from a database).
//
//     fVA1MipRange
//         The dymamic MIP range of the VA1_ALICE pre-amplifier chip 
//
//     fAltroChannelSize
//         The largest number plus one that can be stored in one
//         channel in one time step in the ALTRO ADC chip. 
//
//     fSampleRate
//         How many times the ALTRO ADC chip samples the VA1_ALICE
//         pre-amplifier signal.   The VA1_ALICE chip is read-out at
//         10MHz, while it's possible to drive the ALTRO chip at
//         25MHz.  That means, that the ALTRO chip can have time to
//         sample each VA1_ALICE signal up to 2 times.  Although it's
//         not certain this feature will be used in the production,
//         we'd like have the option, and so it should be reflected in
//         the code.
//
//
// The shaping function of the VA1_ALICE is generally given by 
//
//      f(x) = A(1 - exp(-Bx))
//
// where A is the total charge collected in the pre-amp., and B is a
// paramter that depends on the shaping time of the VA1_ALICE circut.
// 
// When simulating the shaping function of the VA1_ALICe
// pre-amp. chip, we have to take into account, that the shaping
// function depends on the previous value of read from the pre-amp. 
//
// That results in the following algorithm:
//
//    last = 0;
//    FOR charge IN pre-amp. charge train DO 
//      IF last < charge THEN 
//        f(t) = (charge - last) * (1 - exp(-B * t)) + last
//      ELSE
//        f(t) = (last - charge) * exp(-B * t) + charge)
//      ENDIF
//      FOR i IN # samples DO 
//        adc_i = f(i / (# samples))
//      DONE
//      last = charge
//   DONE
//
// Here, 
//
//   pre-amp. charge train 
//       is a series of 128 charges read from the VA1_ALICE chip
//
//   # samples
//       is the number of times the ALTRO ADC samples each of the 128
//       charges from the pre-amp. 
//
// Where Q is the total charge collected by the VA1_ALICE
// pre-amplifier.   Q is then given by 
//
//           E S 
//      Q =  - -
//           e R
//
// where E is the total energy deposited in a silicon strip, R is the
// dynamic range of the VA1_ALICE pre-amp (fVA1MipRange), e is the
// energy deposited by a single MIP, and S ALTRO channel size in each
// time step (fAltroChannelSize).  
//
// The energy deposited per MIP is given by 
//
//      e = M * rho * w 
//
// where M is the universal number 1.664, rho is the density of
// silicon, and w is the depth of the silicon sensor. 
//
// The final ADC count is given by 
//
//      C' = C + P
//
// where P is the (randomized) pedestal (see MakePedestal)
//
// This class uses the class template AliFMDMap<Type> to make an
// internal cache of the energy deposted of the hits.  The class
// template is instantasized as 
//
//  typedef AliFMDMap<std::pair<Float_t, UShort_t> > AliFMDEdepMap;
//
// The first member of the values is the summed energy deposition in a
// given strip, while the second member of the values is the number of
// hits in a given strip.  Using the second member, it's possible to
// do some checks on just how many times a strip got hit, and what
// kind of error we get in our reconstructed hits.  Note, that this
// information is currently not written to the digits tree.  I think a
// QA (Quality Assurance) digit tree is better suited for that task.
// However, the information is there to be used in the future. 
//
//
// Latest changes by Christian Holm Christensen
//
//////////////////////////////////////////////////////////////////////////////

//      /1
//      |           A(-1 + B + exp(-B))
//      | f(x) dx = ------------------- = 1
//      |                    B
//      / 0
//
// and B is the a parameter defined by the shaping time (fShapingTime).  
//
// Solving the above equation, for A gives
//
//                 B
//      A = ----------------
//          -1 + B + exp(-B)
//
// So, if we define the function g: [0,1] -> [0:1] by 
//
//               / v
//               |              Bu + exp(-Bu) - Bv - exp(-Bv) 
//      g(u,v) = | f(x) dx = -A -----------------------------
//               |                            B
//               / u
//
// we can evaluate the ALTRO sample of the VA1_ALICE pre-amp between
// any two times (u, v), by 
//       
//
//                                B	    Bu + exp(-Bu) - Bv - exp(-Bv)
//      C = Q g(u,v) = - Q ---------------- -----------------------------
//		           -1 + B + exp(-B)              B	            
//
//               Bu + exp(-Bu) - Bv - exp(-Bv) 
//        = -  Q -----------------------------
//                    -1 + B + exp(-B)
//

#include <TTree.h>		// ROOT_TTree
#include "AliFMDDebug.h" 	// Better debug macros
#include "AliFMDSSDigitizer.h"	// ALIFMDSSDIGITIZER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDSDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDParameters.h"   // ALIFMDPARAMETERS_H
#include <AliRunDigitizer.h>	// ALIRUNDIGITIZER_H
#include <AliRun.h>		// ALIRUN_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRunLoader.h>	// ALIRUNLOADER_H
    
//====================================================================
ClassImp(AliFMDSSDigitizer)

//____________________________________________________________________
void
AliFMDSSDigitizer::SumContributions(AliFMD* fmd) 
{
  AliFMDDebug(1, ("Runnin our version of SumContributions"));

  // Sum energy deposited contributions from each hit in a cache
  // (fEdep).  
  if (!fRunLoader) 
    Fatal("SumContributions", "no run loader");
  
  // Clear array of deposited energies 
  fEdep.Reset();
  
  // Get the FMD loader 
  AliLoader* inFMD = fRunLoader->GetLoader("FMDLoader");
  // And load the hits 
  inFMD->LoadSDigits("READ");
  
  // Get the tree of hits 
  TTree* sdigitsTree = inFMD->TreeS();
  if (!sdigitsTree)  {
    // Try again 
    // inFMD->LoadSDigits("READ");
    // sdigitsTree = inFMD->TreeH();
    AliError("No sdigit tree from manager");
  }
  
  // Get the FMD branch 
  TBranch* sdigitsBranch = sdigitsTree->GetBranch("FMD");
  if (sdigitsBranch) fmd->SetSDigitsAddressBranch(sdigitsBranch);
  else            AliFatal("Branch FMD hit not found");
  
  // Get a list of hits from the FMD manager 
  TClonesArray *fmdSDigits = fmd->SDigits();
  
  // Get number of entries in the tree 
  Int_t nevents  = Int_t(sdigitsTree->GetEntries());
  
  AliFMDParameters* param = AliFMDParameters::Instance();
  Int_t read = 0;
  // Loop over the events in the 
  for (Int_t event = 0; event < nevents; event++)  {
    // Read in entry number `event' 
    read += sdigitsBranch->GetEntry(event);
    
    // Get the number of sdigits 
    Int_t nsdigits = fmdSDigits->GetEntries ();
    AliFMDDebug(1, ("Got %5d SDigits", nsdigits));
    for (Int_t sdigit = 0; sdigit < nsdigits; sdigit++) {
      // Get the sdigit number `sdigit'
      AliFMDSDigit* fmdSDigit = 
	static_cast<AliFMDSDigit*>(fmdSDigits->UncheckedAt(sdigit));
      
      // Extract parameters 
      UShort_t detector = fmdSDigit->Detector();
      Char_t   ring     = fmdSDigit->Ring();
      UShort_t sector   = fmdSDigit->Sector();
      UShort_t strip    = fmdSDigit->Strip();
      Float_t  edep     = fmdSDigit->Edep();
      // UShort_t minstrip = param->GetMinStrip(detector, ring, sector, strip);
      // UShort_t maxstrip = param->GetMaxStrip(detector, ring, sector, strip);
      // Check if strip is `dead' 
      AliFMDDebug(10, ("SDigit in FMD%d%c[%2d,%3d]=%f",
		      detector, ring, sector, strip, edep));
      if (param->IsDead(detector, ring, sector, strip)) { 
	AliFMDDebug(5, ("FMD%d%c[%2d,%3d] is marked as dead", 
			 detector, ring, sector, strip));
	continue;
      }
      // Check if strip is out-side read-out range 
      // if (strip < minstrip || strip > maxstrip) {
      //   AliFMDDebug(5, ("FMD%d%c[%2d,%3d] is outside range [%3d,%3d]", 
      //		    detector,ring,sector,strip,minstrip,maxstrip));
      //   continue;
      // }
	
      // Give warning in case of double sdigit 
      if (fEdep(detector, ring, sector, strip).fEdep != 0)
	AliFMDDebug(5, ("Double sdigit in %d%c(%d,%d)", 
			 detector, ring, sector, strip));
      
      // Sum energy deposition
      fEdep(detector, ring, sector, strip).fEdep  += edep;
      fEdep(detector, ring, sector, strip).fN     += 1;
      // Add this to the energy deposited for this strip
    }  // sdigit loop
  } // event loop
  AliFMDDebug(3, ("Size of cache: %d bytes, read %d bytes", 
		   sizeof(fEdep), read));
}

//____________________________________________________________________
//
// EOF
// 




