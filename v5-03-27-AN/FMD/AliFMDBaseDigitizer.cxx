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
// #include <AliLog.h>		// ALILOG_H
#include "AliFMDDebug.h" // Better debug macros
#include "AliFMDBaseDigitizer.h" // ALIFMDDIGITIZER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"	        // ALIFMDRING_H
#include "AliFMDHit.h"		// ALIFMDHIT_H
// #include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDParameters.h"   // ALIFMDPARAMETERS_H
// #include <AliDigitizationInput.h>	// ALIRUNDIGITIZER_H
//#include <AliRun.h>		// ALIRUN_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRun.h>		// ALILOADER_H
#include <AliRunLoader.h>	// ALIRUNLOADER_H
#include <TRandom.h>
    
//====================================================================
ClassImp(AliFMDBaseDigitizer)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDBaseDigitizer::AliFMDBaseDigitizer()  
  : fFMD(0),
    fRunLoader(0),
    fEdep(AliFMDMap::kMaxDetectors, 
	  AliFMDMap::kMaxRings, 
	  AliFMDMap::kMaxSectors, 
	  AliFMDMap::kMaxStrips),
    fShapingTime(6),
    fStoreTrackRefs(kTRUE), 
    fIgnoredLabels(0)
{
  AliFMDDebug(1, ("Constructed"));
  // Default ctor - don't use it
}

//____________________________________________________________________
AliFMDBaseDigitizer::AliFMDBaseDigitizer(AliDigitizationInput* digInput) 
  : AliDigitizer(digInput, "AliFMDBaseDigitizer", "FMD Digitizer base class"), 
    fFMD(0),
    fRunLoader(0),
    fEdep(0),        // nDet==0 means 51200 slots
    fShapingTime(6),
    fStoreTrackRefs(kTRUE), 
    fIgnoredLabels(0)
{
  // Normal CTOR
  AliFMDDebug(1, ("Constructed"));
  SetShapingTime();
}

//____________________________________________________________________
AliFMDBaseDigitizer::AliFMDBaseDigitizer(const Char_t* name, 
					 const Char_t* title) 
  : AliDigitizer(name, title),
    fFMD(0),
    fRunLoader(0),
    fEdep(0),        // nDet==0 means 51200 slots
    fShapingTime(6),
    fStoreTrackRefs(kTRUE), 
    fIgnoredLabels(0)
{
  // Normal CTOR
  AliFMDDebug(1, (" Constructed"));
  SetShapingTime();
}

//____________________________________________________________________
AliFMDBaseDigitizer::~AliFMDBaseDigitizer()
{
  // Destructor
}

//____________________________________________________________________
AliFMDBaseDigitizer&
AliFMDBaseDigitizer::operator=(const AliFMDBaseDigitizer& o) 
{ 
  // 
  // Assignment operator
  // 
  // Return:
  //    Reference to this object 
  //
  if (&o == this) return *this; 
  AliDigitizer::operator=(o);
  fRunLoader      = o.fRunLoader;
  fEdep           = o.fEdep;
  fShapingTime    = o.fShapingTime;
  fStoreTrackRefs = o.fStoreTrackRefs;
  fIgnoredLabels  = o.fIgnoredLabels;
  return *this; 
}

//____________________________________________________________________
Bool_t 
AliFMDBaseDigitizer::Init()
{
  // Initialization.   Get a pointer to the parameter manager, and
  // initialize it.  
  AliFMDParameters::Instance()->Init();
  if (AliLog::GetDebugLevel("FMD","") >= 15) 
    AliFMDParameters::Instance()->Print("");
  return kTRUE;
}

//____________________________________________________________________
UShort_t
AliFMDBaseDigitizer::MakePedestal(UShort_t detector, 
				  Char_t   ring, 
				  UShort_t sector, 
				  UShort_t strip) const 
{ 
  // Make a pedestal.  The pedestal value is drawn from a Gaussian
  // distribution.  The mean of the distribution is the measured
  // pedestal, and the width is the measured noise. 
  AliFMDParameters* param =AliFMDParameters::Instance();
  Float_t           mean  =param->GetPedestal(detector,ring,sector,strip);
  Float_t           width =param->GetPedestalWidth(detector,ring,sector,strip);
  return UShort_t(TMath::Max(gRandom->Gaus(mean, width), 0.));
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::AddContribution(UShort_t detector, 
				     Char_t   ring, 
				     UShort_t sector, 
				     UShort_t strip, 
				     Float_t  edep, 
				     Bool_t   isPrimary,
				     Int_t    nTrack,
				     Int_t*   tracknos)
{
  // Add edep contribution from (detector,ring,sector,strip) to cache
  AliFMDParameters* param = AliFMDParameters::Instance();
  AliFMDDebug(10, ("Adding contribution %7.5f for FMD%d%c[%2d,%3d] "
		  " from %d tracks (%s)", 
		  edep,
		  detector, 
		  ring,
		  sector, 
		  strip, 
		  nTrack, 
		  (isPrimary ? "primary" : "secondary")));
  // Check if strip is `dead' 
  if (param->IsDead(detector, ring, sector, strip)) { 
    AliFMDDebug(5, ("FMD%d%c[%2d,%3d] is marked as dead", 
		    detector, ring, sector, strip));
    return;
  }
  // Check if strip is out-side read-out range 
  // if (strip < minstrip || strip > maxstrip) {
  //   AliFMDDebug(5, ("FMD%d%c[%2d,%3d] is outside range [%3d,%3d]", 
  //		    detector,ring,sector,strip,minstrip,maxstrip));
  //   continue;
  // }
  
  AliFMDEdepHitPair& entry = fEdep(detector, ring, sector, strip);

  // Give warning in case of double sdigit 
  if (entry.fEdep != 0)
    AliFMDDebug(5, ("Double digit in FMD%d%c[%2d,%3d]", 
		    detector, ring, sector, strip));
      
  // Sum energy deposition
  Int_t oldN  =  entry.fN;
  entry.fEdep += edep;
  entry.fN    += nTrack;
  if (isPrimary) entry.fNPrim += nTrack;
  if (fStoreTrackRefs) { 
    if (entry.fLabels.fN < entry.fN) {
      AliFMDDebug(15, ("== New label array size %d, was %d, added %d", 
		       entry.fN, entry.fLabels.fN, nTrack));
      entry.fLabels.Set(entry.fN);
    }
    for (Int_t i = 0; i < nTrack; i++) {
      AliFMDDebug(15, ("=> Setting track label # %d", oldN+i));
      entry.fLabels[oldN + i] = tracknos[i];
      AliFMDDebug(15, ("<= Setting track label # %d", oldN+i));
    }
  }
  AliFMDDebug(15,("Adding contribution %f to FMD%d%c[%2d,%3d] (%f) track %d", 
		  edep, detector, ring, sector, strip,
		  entry.fEdep, (nTrack > 0 ? tracknos[0] : -1)));
  
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::DigitizeHits() const
{
  // For the stored energy contributions in the cache (fEdep), convert
  // the energy signal to ADC counts, and store the created digit in
  // the digits array (AliFMD::fDigits)
  //
  AliFMDDebug(5, ("Will now digitize all the summed signals"));
  fIgnoredLabels = 0;
  AliFMDGeometry* geometry = AliFMDGeometry::Instance();
  
  TArrayI counts(4);
  for (UShort_t detector=1; detector <= 3; detector++) {
    AliFMDDebug(10, ("Processing hits in FMD%d", detector));
    // Get pointer to subdetector 
    AliFMDDetector* det = geometry->GetDetector(detector);
    if (!det) continue;
    for (UShort_t ringi = 0; ringi <= 1; ringi++) {
      Char_t ring = ringi == 0 ? 'I' : 'O';
      AliFMDDebug(10, (" Processing hits in FMD%d%c", detector,ring));
      // Get pointer to Ring
      AliFMDRing* r = det->GetRing(ring);
      if (!r) continue;
      
      // Get number of sectors 
      UShort_t nSectors = UShort_t(360. / r->GetTheta());
      // Loop over the number of sectors 
      for (UShort_t sector = 0; sector < nSectors; sector++) {
	AliFMDDebug(10, ("  Processing hits in FMD%d%c[%2d]", 
			detector,ring,sector));
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
	  
	  const AliFMDEdepHitPair& entry  = fEdep(detector,ring,sector,strip);
	  Float_t                  edep   = entry.fEdep;
	  UShort_t                 ntot   = entry.fN;
	  UShort_t                 nprim  = entry.fNPrim;
	  const TArrayI&           labels = entry.fLabels;
	  if (edep > 0)
	    AliFMDDebug(15, ("Edep = %f for FMD%d%c[%2d,%3d]", 
			     edep, detector, ring, sector, strip));
	  ConvertToCount(edep, last, detector, ring, sector, strip, counts);
	  last = edep;
	  

	  // The following line was introduced - wrongly - by Peter
	  // Hristov.  It _will_ break the digitisation and the
	  // following reconstruction.  The behviour of the
	  // digitisation models exactly the front-end as it should
	  // (no matter what memory concuption it may entail).  The
	  // check should be on zero suppression, since that's what
	  // models the front-end - if zero suppression is turned on
	  // in the front-end, then we can suppress empty digits -
	  // otherwise we shoud never do that.  Note, that the line
	  // affects _both_ normal digitisation and digitisation for
	  // summable digits, since the condition is on the energy
	  // deposition and not on the actual number of counts.  If
	  // this line should go anywhere, it should be in the
	  // possible overloaded AliFMDSDigitizer::AddDigit - not
	  // here. 
	  // 
	  //   if (edep<=0) continue;
	  AddDigit(detector, ring, sector, strip, edep, 
		   UShort_t(counts[0]), Short_t(counts[1]), 
		   Short_t(counts[2]), Short_t(counts[3]), 
		   ntot, nprim, labels);
	  AliFMDDebug(15, ("   Adding digit in FMD%d%c[%2d,%3d]=%d", 
			   detector,ring,sector,strip,counts[0]));
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
  if (fIgnoredLabels > 0) 
    AliWarning(Form("%d track labels could not be associated with digits "
		    "due to limited storage facilities in AliDigit", 
		    fIgnoredLabels));
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
  Float_t  convF          = (param->GetDACPerMIP() / param->GetEdepMip() *
			     param->GetPulseGain(detector,ring,sector,strip));
  Int_t    ped            = MakePedestal(detector,ring,sector,strip);
  Int_t    maxAdc         = param->GetAltroChannelSize()-1;
  if (maxAdc < 0) {
    AliWarning(Form("Maximum ADC is %d < 0, forcing it to 1023", maxAdc));
    maxAdc = 1023;
  }
  UShort_t rate           = param->GetSampleRate(detector,ring,sector,strip);
  AliFMDDebug(15, ("Sample rate for FMD%d%c[%2d,%3d] = %d", 
		   detector, ring, sector, strip, rate));
  if (rate < 1 || rate > 4) {
    AliWarning(Form("Invalid sample rate for for FMD%d%c[%2d,%3d] = %d", 
		    detector, ring, sector, strip, rate));
    rate = 1;
  }

  // In case we don't oversample, just return the end value. 
  if (rate == 1) {
    Float_t    a = edep * convF + ped;
    if (a < 0) a = 0;
    counts[0]    = UShort_t(TMath::Min(a, Float_t(maxAdc)));
    AliFMDDebug(15, ("FMD%d%c[%2d,%3d]: converting ELoss %f to "
		     "ADC %4d (%f,%d)",
		     detector,ring,sector,strip,edep,counts[0],convF,ped));
    return;
  }

  
  // Create a pedestal 
  Float_t b = fShapingTime;
  for (Ssiz_t i = 0; i < rate;  i++) {
    Float_t t  = Float_t(i) / rate + 1./rate;
    Float_t s  = edep + (last - edep) * TMath::Exp(-b * t);
    Float_t a  = Int_t(s * convF + ped);
    if (a < 0) a = 0;
    counts[i]  = UShort_t(TMath::Min(a, Float_t(maxAdc)));
  }
  AliFMDDebug(15, ("Converted edep = %f to ADC (%x,%x,%x,%x) "
		   "[gain: %f=(%f/%f*%f), pedestal: %d, rate: %d]", 
		   edep, counts[0], counts[1], counts[2], counts[3], 
		   convF, param->GetDACPerMIP(),param->GetEdepMip(),
		   param->GetPulseGain(detector,ring,sector,strip), 
		   ped, rate));
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::AddDigit(UShort_t        detector, 
			      Char_t          ring,
			      UShort_t        sector, 
			      UShort_t        strip, 
			      Float_t         /* edep */, 
			      UShort_t        count1, 
			      Short_t         count2, 
			      Short_t         count3,
			      Short_t         count4,
			      UShort_t        ntot, 
			      UShort_t        /* nprim */,
			      const TArrayI&  refs) const
{
  // Add a digit or summable digit
  fFMD->AddDigitByFields(detector, ring, sector, strip, 
			 count1, count2, count3, count4, 
			 ntot, fStoreTrackRefs ? refs.fArray : 0);
  if (fStoreTrackRefs && ntot > 3) fIgnoredLabels += ntot - 3;
}

//____________________________________________________________________
TTree*
AliFMDBaseDigitizer::MakeOutputTree(AliLoader* loader)
{
  // Create output tree using loader.   If the passed loader differs
  // from the currently set loader in the FMD object, reset the FMD
  // loader to be the passed loader.   This is for the cases wher the
  // output is different from the output. 
  AliFMDDebug(5, ("Making digits tree"));
  loader->LoadDigits("UPDATE"); // "RECREATE");
  TTree* out = loader->TreeD();
  if (!out) loader->MakeTree("D");
  out = loader->TreeD(); 
  if (out) { 
    out->Reset();
    if (loader != fFMD->GetLoader()) 
      fFMD->SetLoader(loader);
    fFMD->MakeBranch("D");
  }
  return out;
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::StoreDigits(const AliLoader* loader)
{
  // Write the digits to disk 
  AliFMDDebug(5, ("Storing %d digits",   fFMD->Digits()->GetEntries()));
  loader->WriteDigits("OVERWRITE");
  loader->UnloadDigits();
  // Reset the digits in the AliFMD object 
  fFMD->ResetDigits();
}

//____________________________________________________________________
//
// EOF
// 




