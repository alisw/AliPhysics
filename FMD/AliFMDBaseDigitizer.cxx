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
/** @file    AliFMDBaseDigitizer.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:38:26 2006
    @brief   FMD Digitizers implementation
    @ingroup FMD_sim
*/
//////////////////////////////////////////////////////////////////////////////
//
//  This class contains the procedures simulation ADC  signal for the
//  Forward Multiplicity detector  : Hits->Digits and Hits->SDigits
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

#include <TMath.h>
#include <TTree.h>		// ROOT_TTree
//#include <TRandom.h>		// ROOT_TRandom
#include <AliLog.h>		// ALILOG_H
#include "AliFMDBaseDigitizer.h" // ALIFMDDIGITIZER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"	        // ALIFMDRING_H
#include "AliFMDHit.h"		// ALIFMDHIT_H
// #include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDParameters.h"   // ALIFMDPARAMETERS_H
// #include <AliRunDigitizer.h>	// ALIRUNDIGITIZER_H
//#include <AliRun.h>		// ALIRUN_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRunLoader.h>	// ALIRUNLOADER_H
    
//====================================================================
ClassImp(AliFMDBaseDigitizer)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDBaseDigitizer::AliFMDBaseDigitizer()  
  : fRunLoader(0),
    fEdep(AliFMDMap::kMaxDetectors, 
	  AliFMDMap::kMaxRings, 
	  AliFMDMap::kMaxSectors, 
	  AliFMDMap::kMaxStrips),
    fShapingTime(0)
{
  // Default ctor - don't use it
}

//____________________________________________________________________
AliFMDBaseDigitizer::AliFMDBaseDigitizer(AliRunDigitizer* manager) 
  : AliDigitizer(manager, "AliFMDBaseDigitizer", "FMD Digitizer base class"), 
    fRunLoader(0),
    fEdep(AliFMDMap::kMaxDetectors, 
	  AliFMDMap::kMaxRings, 
	  AliFMDMap::kMaxSectors, 
	  AliFMDMap::kMaxStrips), 
    fShapingTime(0)
{
  // Normal CTOR
  AliDebug(1," processed");
  SetShapingTime();
}

//____________________________________________________________________
AliFMDBaseDigitizer::AliFMDBaseDigitizer(const Char_t* name, 
					 const Char_t* title) 
  : AliDigitizer(name, title),
    fRunLoader(0),
    fEdep(AliFMDMap::kMaxDetectors, 
	  AliFMDMap::kMaxRings, 
	  AliFMDMap::kMaxSectors, 
	  AliFMDMap::kMaxStrips)
{
  // Normal CTOR
  AliDebug(1," processed");
  SetShapingTime();
}

//____________________________________________________________________
AliFMDBaseDigitizer::~AliFMDBaseDigitizer()
{
  // Destructor
}

//____________________________________________________________________
Bool_t 
AliFMDBaseDigitizer::Init()
{
  // Initialization
  AliFMDParameters::Instance()->Init();
  return kTRUE;
}
 

//____________________________________________________________________
UShort_t
AliFMDBaseDigitizer::MakePedestal(UShort_t, 
				  Char_t, 
				  UShort_t, 
				  UShort_t) const 
{ 
  // Make a pedestal
  return 0; 
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::SumContributions(AliFMD* fmd) 
{
  // Sum energy deposited contributions from each hit in a cache
  // (fEdep).  
  if (!fRunLoader) 
    Fatal("SumContributions", "no run loader");
  
  // Clear array of deposited energies 
  fEdep.Reset();
  
  // Get the FMD loader 
  AliLoader* inFMD = fRunLoader->GetLoader("FMDLoader");
  // And load the hits 
  inFMD->LoadHits("READ");
  
  // Get the tree of hits 
  TTree* hitsTree = inFMD->TreeH();
  if (!hitsTree)  {
    // Try again 
    inFMD->LoadHits("READ");
    hitsTree = inFMD->TreeH();
  }
  
  // Get the FMD branch 
  TBranch* hitsBranch = hitsTree->GetBranch("FMD");
  if (hitsBranch) fmd->SetHitsAddressBranch(hitsBranch);
  else            AliFatal("Branch FMD hit not found");
  
  // Get a list of hits from the FMD manager 
  TClonesArray *fmdHits = fmd->Hits();
  
  // Get number of entries in the tree 
  Int_t ntracks  = Int_t(hitsTree->GetEntries());
  
  AliFMDParameters* param = AliFMDParameters::Instance();
  Int_t read = 0;
  // Loop over the tracks in the 
  for (Int_t track = 0; track < ntracks; track++)  {
    // Read in entry number `track' 
    read += hitsBranch->GetEntry(track);
    
    // Get the number of hits 
    Int_t nhits = fmdHits->GetEntries ();
    for (Int_t hit = 0; hit < nhits; hit++) {
      // Get the hit number `hit'
      AliFMDHit* fmdHit = 
	static_cast<AliFMDHit*>(fmdHits->UncheckedAt(hit));
      
      // Extract parameters 
      UShort_t detector = fmdHit->Detector();
      Char_t   ring     = fmdHit->Ring();
      UShort_t sector   = fmdHit->Sector();
      UShort_t strip    = fmdHit->Strip();
      Float_t  edep     = fmdHit->Edep();
      // UShort_t minstrip = param->GetMinStrip(detector, ring, sector, strip);
      // UShort_t maxstrip = param->GetMaxStrip(detector, ring, sector, strip);
      // Check if strip is `dead' 
      if (param->IsDead(detector, ring, sector, strip)) { 
	AliDebug(5, Form("FMD%d%c[%2d,%3d] is marked as dead", 
			 detector, ring, sector, strip));
	continue;
      }
      // Check if strip is out-side read-out range 
      // if (strip < minstrip || strip > maxstrip) {
      //   AliDebug(5, Form("FMD%d%c[%2d,%3d] is outside range [%3d,%3d]", 
      //		    detector,ring,sector,strip,minstrip,maxstrip));
      //   continue;
      // }
	
      // Give warning in case of double hit 
      if (fEdep(detector, ring, sector, strip).fEdep != 0)
	AliDebug(5, Form("Double hit in %d%c(%d,%d)", 
			 detector, ring, sector, strip));
      
      // Sum energy deposition
      fEdep(detector, ring, sector, strip).fEdep  += edep;
      fEdep(detector, ring, sector, strip).fN     += 1;
      // Add this to the energy deposited for this strip
    }  // hit loop
  } // track loop
  AliDebug(1, Form("Size of cache: %d bytes, read %d bytes", 
		   sizeof(fEdep), read));
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::DigitizeHits(AliFMD* fmd) const
{
  // For the stored energy contributions in the cache (fEdep), convert
  // the energy signal to ADC counts, and store the created digit in
  // the digits array (AliFMD::fDigits)
  //
  AliFMDGeometry* geometry = AliFMDGeometry::Instance();
  
  TArrayI counts(3);
  for (UShort_t detector=1; detector <= 3; detector++) {
    // Get pointer to subdetector 
    AliFMDDetector* det = geometry->GetDetector(detector);
    if (!det) continue;
    for (UShort_t ringi = 0; ringi <= 1; ringi++) {
      Char_t ring = ringi == 0 ? 'I' : 'O';
      // Get pointer to Ring
      AliFMDRing* r = det->GetRing(ring);
      if (!r) continue;
      
      // Get number of sectors 
      UShort_t nSectors = UShort_t(360. / r->GetTheta());
      // Loop over the number of sectors 
      for (UShort_t sector = 0; sector < nSectors; sector++) {
	// Get number of strips 
	UShort_t nStrips = r->GetNStrips();
	// Loop over the stips 
	Float_t last = 0;
	for (UShort_t strip = 0; strip < nStrips; strip++) {
	  // Reset the counter array to the invalid value -1 
	  counts.Reset(-1);
	  // Reset the last `ADC' value when we've get to the end of a
	  // VA1_ALICE channel. 
	  if (strip % 128 == 0) last = 0;
	  
	  Float_t edep = fEdep(detector, ring, sector, strip).fEdep;
	  ConvertToCount(edep, last, detector, ring, sector, strip, counts);
	  last = edep;
	  AddDigit(fmd, detector, ring, sector, strip, edep, 
		   UShort_t(counts[0]), Short_t(counts[1]), 
		   Short_t(counts[2]));
#if 0
	  // This checks if the digit created will give the `right'
	  // number of particles when reconstructed, using a naiive
	  // approach.  It's here only as a quality check - nothing
	  // else. 
	  CheckDigit(digit, fEdep(detector, ring, sector, strip).fN,
		     counts);
#endif
	} // Strip
      } // Sector 
    } // Ring 
  } // Detector 
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::ConvertToCount(Float_t   edep, 
				    Float_t   last,
				    UShort_t  detector, 
				    Char_t    ring, 
				    UShort_t  sector, 
				    UShort_t  strip,
				    TArrayI&  counts) const
{
  // Convert the total energy deposited to a (set of) ADC count(s). 
  // 
  // This is done by 
  // 
  //               Energy_Deposited      ALTRO_Channel_Size
  //    ADC = -------------------------- ------------------- + pedestal
  //          Energy_Deposition_Of_1_MIP VA1_ALICE_MIP_Range
  //
  //               Energy_Deposited             fAltroChannelSize
  //        = --------------------------------- ----------------- + pedestal 
  //          1.664 * Si_Thickness * Si_Density   fVA1MipRange	 
  //          
  // 
  //        = Energy_Deposited * ConversionFactor + pedestal
  // 
  // However, this is modified by the response function of the
  // VA1_ALICE pre-amp. chip in case we are doing oversampling of the
  // VA1_ALICE output. 
  // 
  // In that case, we get N=fSampleRate values of the ADC, and the
  // `EnergyDeposited' is a function of which sample where are
  // calculating the ADC for 
  // 
  //     ADC_i = f(EnergyDeposited, i/N, Last) * ConversionFactor + pedestal 
  // 
  // where Last is the Energy deposited in the previous strip. 
  // 
  // Here, f is the shaping function of the VA1_ALICE.   This is given
  // by 
  //                       
  //                    |   (E - l) * (1 - exp(-B * t) + l   if E > l
  //       f(E, t, l) = <
  //                    |   (l - E) * exp(-B * t) + E        otherwise
  //                       
  // 
  //                  = E + (l - E) * ext(-B * t)
  // 
  AliFMDParameters* param = AliFMDParameters::Instance();
  Float_t  convF          = 1./param->GetPulseGain(detector,ring,sector,strip);
  Int_t    ped            = MakePedestal(detector,ring,sector,strip);
  Int_t    maxAdc         = param->GetAltroChannelSize()-1;
  if (maxAdc < 0) {
    AliWarning(Form("Maximum ADC is %d < 0, forcing it to 1023", maxAdc));
    maxAdc = 1023;
  }
  UShort_t rate           = param->GetSampleRate(detector,ring,sector,strip);
  if (rate < 1 || rate > 3) rate = 1;
  
  // In case we don't oversample, just return the end value. 
  if (rate == 1) {
    Float_t    a = edep * convF + ped;
    if (a < 0) a = 0;
    counts[0]    = UShort_t(TMath::Min(a, Float_t(maxAdc)));
    AliDebug(2, Form("FMD%d%c[%2d,%3d]: converting ELoss %f to "
		     "ADC %4d (%f,%d)",
		     detector,ring,sector,strip,edep,counts[0],convF,ped));
    return;
  }
  
  // Create a pedestal 
  Float_t b = fShapingTime;
  for (Ssiz_t i = 0; i < rate;  i++) {
    Float_t t  = Float_t(i) / rate;
    Float_t s  = edep + (last - edep) * TMath::Exp(-b * t);
    Float_t a  = Int_t(s * convF + ped);
    if (a < 0) a = 0;
    counts[i]  = UShort_t(TMath::Min(a, Float_t(maxAdc)));
  }
}



//____________________________________________________________________
//
// EOF
// 




