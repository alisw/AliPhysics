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
/** @file    AliFMDDigitizer.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:38:26 2006
    @brief   FMD Digitizers implementation
    @ingroup FMD_sim
*/
//////////////////////////////////////////////////////////////////////////////
//
//  This class contains the procedures simulation ADC  signal for the
//  Forward Multiplicity detector  : Hits->Digits
// 
//  Digits consists of
//   - Detector #
//   - Ring ID                                             
//   - Sector #     
//   - Strip #
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
// These parameters are fetched from OCDB via the mananger AliFMDParameters.
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
#include <TRandom.h>		// ROOT_TRandom
#include "AliFMDDebug.h"        // Better debug macros
#include "AliFMDDigitizer.h"	// ALIFMDDIGITIZER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDParameters.h"   // ALIFMDPARAMETERS_H
#include <AliRunDigitizer.h>	// ALIRUNDIGITIZER_H
#include <AliRun.h>		// ALIRUN_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRunLoader.h>	// ALIRUNLOADER_H
    
//====================================================================
ClassImp(AliFMDDigitizer)

//____________________________________________________________________
void 
AliFMDDigitizer::OutputTree(AliLoader* outFMD, AliFMD* fmd)
{
  // Load digits from the tree 
  outFMD->LoadDigits("update");

  // Get the tree of digits 
  TTree* digitTree = outFMD->TreeD();
  if (!digitTree) {
    outFMD->MakeTree("D");
    digitTree = outFMD->TreeD();
  }
  digitTree->Reset();
  
  // Get the digits 
  TClonesArray* digits =  fmd->Digits();
  if (!digits) { 
    AliError("Failed to get digits");
    return;
  }
  AliFMDDebug(1, ("Got a total of %5d digits", digits->GetEntries()));

  // Make a branch in the tree 
  fmd->MakeBranchInTree(digitTree, fmd->GetName(), &(digits), 4000, 0);
  // TBranch* digitBranch = digitTree->GetBranch(fmd->GetName());
  // Fill the tree 
  Int_t write = 0;
  write = digitTree->Fill();
  AliFMDDebug(1, ("Wrote %d bytes to digit tree", write));
  
  // Write the digits to disk 
  outFMD->WriteDigits("OVERWRITE");
  outFMD->UnloadHits();
  outFMD->UnloadDigits();

  // Reset the digits in the AliFMD object 
  fmd->ResetDigits();
}

//____________________________________________________________________
UShort_t
AliFMDDigitizer::MakePedestal(UShort_t  detector, 
			      Char_t    ring, 
			      UShort_t  sector, 
			      UShort_t  strip) const 
{
  // Make a pedestal 
  AliFMDParameters* param =AliFMDParameters::Instance();
  Float_t           mean  =param->GetPedestal(detector,ring,sector,strip);
  Float_t           width =param->GetPedestalWidth(detector,ring,sector,strip);
  return UShort_t(TMath::Max(gRandom->Gaus(mean, width), 0.));
}

//____________________________________________________________________
void
AliFMDDigitizer::AddDigit(AliFMD*  fmd,
			  UShort_t detector, 
			  Char_t   ring,
			  UShort_t sector, 
			  UShort_t strip, 
			  Float_t  /* edep */, 
			  UShort_t count1, 
			  Short_t  count2, 
			  Short_t  count3,
			  Short_t  count4) const
{
  // Add a digit
  fmd->AddDigitByFields(detector, ring, sector, strip, 
			count1, count2, count3, count4);
}

//____________________________________________________________________
void
AliFMDDigitizer::CheckDigit(AliFMDDigit*    digit,
			    UShort_t        nhits,
			    const TArrayI&  counts) 
{
  // Check that digit is consistent
  AliFMDParameters* param = AliFMDParameters::Instance();
  UShort_t          det   = digit->Detector();
  Char_t            ring  = digit->Ring();
  UShort_t          sec   = digit->Sector();
  UShort_t          str   = digit->Strip();
  Float_t           mean  = param->GetPedestal(det,ring,sec,str);
  Float_t           width = param->GetPedestalWidth(det,ring,sec,str);
  UShort_t          range = param->GetVA1MipRange();
  UShort_t          size  = param->GetAltroChannelSize();
  Int_t             integral = counts[0];
  if (counts[1] >= 0) integral += counts[1];
  if (counts[2] >= 0) integral += counts[2];
  if (counts[3] >= 0) integral += counts[3];
  integral -= Int_t(mean + 2 * width);
  if (integral < 0) integral = 0;
  
  Float_t convF = Float_t(range) / size;
  Float_t mips  = integral * convF;
  if (mips > Float_t(nhits) + .5 || mips < Float_t(nhits) - .5) 
    Warning("CheckDigit", "Digit -> %4.2f MIPS != %d +/- .5 hits", 
	    mips, nhits);
}

//____________________________________________________________________
//
// EOF
// 




