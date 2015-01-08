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
/* $Id: AliFMDHitDigitizer.cxx 28055 2008-08-18 00:33:20Z cholm $ */
/** @file    AliFMDHitDigitizer.cxx
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
// AliFMDHitDigitizer and AliFMDSDigitizer are implemented via a base
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
//      | AliFMDHitDigitizer |	| AliFMDSDigitizer |
//      +-----------------+	+------------------+
//
// These classes has several paramters: 
//
//     fPedestal
//     fPedestalWidth
//         (Only AliFMDHitDigitizer)
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
#include "AliFMDDebug.h"        // Better debug macros
#include "AliFMDHitDigitizer.h"	// ALIFMDDIGITIZER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDParameters.h"   // ALIFMDPARAMETERS_H
#include <AliRun.h>		// ALIRUN_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRunLoader.h>	// ALIRUNLOADER_H
#include <AliFMDHit.h>
#include <AliStack.h>
#include <TFile.h>
#include <TParticle.h>

//====================================================================
ClassImp(AliFMDHitDigitizer)    
#if 0
;
#endif

//____________________________________________________________________
AliFMDHitDigitizer::AliFMDHitDigitizer(AliFMD* fmd, Output_t  output)
  : AliFMDBaseDigitizer("FMD", (output == kDigits ? 
				"FMD Hit->Digit digitizer" :
				"FMD Hit->SDigit digitizer")),
    fOutput(output), 
    fHoldTime(2e-6),
    fStack(0)
{
  fFMD = fmd;
}

//____________________________________________________________________
AliFMDHitDigitizer& 
AliFMDHitDigitizer::operator=(const AliFMDHitDigitizer& o) 
{
  /** 
   * Assignment operator
   *
   * @param o Object to assign from 
   * @return Reference to this 
   */
  if (&o == this) return *this; 
  AliFMDBaseDigitizer::operator=(o);
  fHoldTime    = o.fHoldTime;
  fOutput      = o.fOutput;
  fStack       = o.fStack;
  return *this;
}

//____________________________________________________________________
void
AliFMDHitDigitizer::Digitize(Option_t* /*option*/)
{
  // Run this digitizer 
  // Get an inititialize parameter manager
  AliFMDParameters::Instance()->Init();
  if (AliLog::GetDebugLevel("FMD","") >= 10) 
    AliFMDParameters::Instance()->Print("ALL");

  // Get loader, and ask it to read in the hits 
  AliLoader* loader = fFMD->GetLoader();
  if (!loader) { 
    AliError("Failed to get loader from detector object");
    return;
  }
  loader->LoadHits("READ");
  
  // Get the run loader 
  AliRunLoader* runLoader = loader->GetRunLoader();
  if (!runLoader) {
    AliError("Failed to get run loader from loader");
    return;
  }
  
  // Now loop over events
  Int_t nEvents = runLoader->GetNumberOfEvents();
  for (Int_t event = 0; event < nEvents; event++) { 
    // Get the current event folder. 
    TFolder* folder = loader->GetEventFolder();
    if (!folder) { 
      AliError("Failed to get event folder from loader");
      return;
    }

    // Get the run-loader of this event. 
    const char* loaderName = AliRunLoader::GetRunLoaderName();
    AliRunLoader* thisLoader = 
      static_cast<AliRunLoader*>(folder->FindObject(loaderName));
    if (!thisLoader) { 
      AliError(Form("Failed to get loader '%s' from event folder", loaderName));
      return;
    }
    
    // Read in the event
    AliFMDDebug(1, ("Now digitizing (Hits->%s) event # %d", 
		    (fOutput == kDigits ? "digits" : "sdigits"), event));
    thisLoader->GetEvent(event);
    
    // Load kinematics to get primary information for SDigits
    fStack = 0;
    if (fOutput == kSDigits) {
      if (thisLoader->LoadKinematics("READ")) {
	AliError("Failed to get kinematics from event loader");
	return;
      }
      AliFMDDebug(5, ("Loading stack of kinematics"));
      fStack = thisLoader->Stack();
    }

    // Check that we have the hits 
    if (!loader->TreeH() && loader->LoadHits()) {
      AliError("Failed to load hits");
      return;
    }
    TTree*   hitsTree = loader->TreeH();
    TBranch* hitsBranch = hitsTree->GetBranch(fFMD->GetName());
    if (!hitsBranch) { 
      AliError("Failed to get hits branch in tree");
      return;
    }
    // Check that we can make the output digits - This must come
    // before AliFMD::SetBranchAddress
    TTree* outTree = MakeOutputTree(loader);
    if (!outTree) { 
      AliError("Failed to get output tree");
      return;
    }
    AliFMDDebug(5, ("Output tree name for %s is '%s'", 
		    (fOutput == kDigits ? "digits" : "sdigits"),
		    outTree->GetName()));
    if (AliLog::GetDebugLevel("FMD","") >= 5) {
      TFile* file = outTree->GetCurrentFile();
      if (!file) {
	AliWarning("Output tree has no file!");
      }
      else { 
	AliFMDDebug(5, ("Output tree file %s content:", file->GetName()));
	file->ls();
      }
    }

    // Set-up the branch addresses 
    fFMD->SetTreeAddress();
    
    // Now sum all contributions in cache 
    SumContributions(hitsBranch);
    loader->UnloadHits();

    // And now digitize the hits 
    DigitizeHits();
    
    // Write digits to tree
    Int_t write = outTree->Fill();
    AliFMDDebug(5, ("Wrote %d bytes to digit tree", write));

    // Store the digits
    StoreDigits(loader);

  }  
}

//____________________________________________________________________
TTree*
AliFMDHitDigitizer::MakeOutputTree(AliLoader* loader)
{
  /** 
   * Make the output tree using the passed loader 
   *
   * @param loader 
   * @return The generated tree. 
   */
  if (fOutput == kDigits) 
    return AliFMDBaseDigitizer::MakeOutputTree(loader);
  
  AliFMDDebug(5, ("Making sdigits tree"));
  loader->LoadSDigits("UPDATE"); // RECREATE");
  TTree* out = loader->TreeS();
  if (!out) loader->MakeTree("S");
  out = loader->TreeS(); 
  if (out) { 
    out->Reset();
    fFMD->MakeBranch("S");
  }
  return out;
}


//____________________________________________________________________
void
AliFMDHitDigitizer::SumContributions(TBranch* hitsBranch) 
{
  // Sum energy deposited contributions from each hit in a cache
  // (fEdep).  
  
  // Clear array of deposited energies 
  fEdep.Reset();
  
  // Get a list of hits from the FMD manager 
  AliFMDDebug(5, ("Get array of FMD hits"));
  TClonesArray *fmdHits = fFMD->Hits();
  

  // Get number of entries in the tree 
  AliFMDDebug(5, ("Get # of tracks"));
  Int_t ntracks  = Int_t(hitsBranch->GetEntries());
  AliFMDDebug(5, ("We got %d tracks", ntracks));

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

      // Ignore hits that arrive too late
      if (fmdHit->Time() > fHoldTime) continue;
      

      // Check if this is a primary particle
      Bool_t isPrimary = kTRUE;
      Int_t  trackno   = -1;
      if (fStack) {
	trackno = fmdHit->Track();
	AliFMDDebug(10, ("Will get track # %d/%d from entry # %d", 
			trackno, fStack->GetNtrack(), track));
	if (fStack->GetNtrack() < trackno) {
	  AliError(Form("Track number %d/%d out of bounds", 
			trackno, fStack->GetNtrack()));
	  continue;
	}
#if 1
	isPrimary = fStack->IsPhysicalPrimary(trackno);
#else // This is our hand-crafted code.  We use the ALICE definition
	TParticle* part    = fStack->Particle(trackno);
	isPrimary          = part->IsPrimary();
	if (!isPrimary) { 
	  // Extended testing of mother status - this is for Pythia6.
	  Int_t      mother1   = part->GetFirstMother();
	  TParticle* mother    = fStack->Particle(mother1);
	  if (!mother || mother->GetStatusCode() > 1)
	    isPrimary = kTRUE;
	  AliFMDDebug(15,
		      ("Track %d secondary, mother: %d - %s - status %d: %s", 
		       trackno, mother1, 
		       (mother ? "found"                 : "not found"), 
		       (mother ? mother->GetStatusCode() : -1),
		       (isPrimary ? "primary" : "secondary")));
	}
#endif
      }
    
      // Extract parameters 
      AliFMDDebug(15,("Adding contribution %7.5f for FMD%d%c[%2d,%3d] "
		      " for trackno %6d (%s)", 
		      fmdHit->Edep(),
		      fmdHit->Detector(), 
		      fmdHit->Ring(),
		      fmdHit->Sector(), 
		      fmdHit->Strip(), 
		      trackno, 
		      (isPrimary ? "primary" : "secondary")));
      AddContribution(fmdHit->Detector(),
		      fmdHit->Ring(),
		      fmdHit->Sector(),
		      fmdHit->Strip(),
		      fmdHit->Edep(), 
		      isPrimary, 
		      1, 
		      &trackno);
    }  // hit loop
  } // track loop
  AliFMDDebug(5, ("Size of cache: %d bytes, read %d bytes", 
		  int(sizeof(fEdep)), read));
}


//____________________________________________________________________
UShort_t
AliFMDHitDigitizer::MakePedestal(UShort_t  detector, 
				 Char_t    ring, 
				 UShort_t  sector, 
				 UShort_t  strip) const 
{
  // Make a pedestal 
  if (fOutput == kSDigits) return 0;
  return AliFMDBaseDigitizer::MakePedestal(detector, ring, sector, strip);
}



//____________________________________________________________________
void
AliFMDHitDigitizer::AddDigit(UShort_t        detector, 
			     Char_t          ring,
			     UShort_t        sector, 
			     UShort_t        strip, 
			     Float_t         edep, 
			     UShort_t        count1, 
			     Short_t         count2, 
			     Short_t         count3,
			     Short_t         count4, 
			     UShort_t        ntotal,
			     UShort_t        nprim, 
			     const TArrayI&  refs) const
{
  // Add a digit or summable digit
  if (fOutput == kDigits) { 
    AliFMDDebug(15,("Adding digit for FMD%d%c[%2d,%3d] = (%x,%x,%x,%x)",
		    detector, ring, sector, strip, 
		    count1, count2, count3, count4));
    AliFMDBaseDigitizer::AddDigit(detector, ring, sector, strip, 0,
				  count1, count2, count3, count4, 
				  ntotal, nprim, refs);
    return;
  }
  if (edep <= 0) { 
    AliFMDDebug(15, ("Digit edep = %f <= 0 for FMD%d%c[%2d,%3d]", 
		    edep, detector, ring, sector, strip));
    return;
  }
  if (count1 == 0 && count2 <= 0 && count3 <= 0 && count4 <= 0) {
    AliFMDDebug(15, ("Digit counts = (%x,%x,%x,%x) <= 0 for FMD%d%c[%2d,%3d]", 
		    count1, count2, count3, count4, 
		    detector, ring, sector, strip));
    return;
  }
  AliFMDDebug(15, ("Adding sdigit for FMD%d%c[%2d,%3d] = "
		   "(%x,%x,%x,%x) [%d/%d] %d",
		   detector, ring, sector, strip, 
		   count1, count2, count3, count4, nprim, ntotal, refs.fN));
  fFMD->AddSDigitByFields(detector, ring, sector, strip, edep,
			  count1, count2, count3, count4, 
			  ntotal, nprim, fStoreTrackRefs ? refs.fArray : 0);
  if (fStoreTrackRefs && nprim > 3) fIgnoredLabels += nprim - 3;
}

//____________________________________________________________________
void
AliFMDHitDigitizer::CheckDigit(AliFMDDigit*    digit,
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
void
AliFMDHitDigitizer::StoreDigits(const AliLoader* loader)
{
  /** 
   * Store the data using the loader 
   *
   * @param loader The loader 
   */
  if (fOutput == kDigits) { 
    AliFMDBaseDigitizer::StoreDigits(loader);
    return;
  }
  AliFMDDebug(5, ("Storing %d sdigits",   fFMD->SDigits()->GetEntries()));
  // Write the digits to disk 
  loader->WriteSDigits("OVERWRITE");
  loader->UnloadSDigits();
  // Reset the digits in the AliFMD object 
  fFMD->ResetSDigits();
}


//____________________________________________________________________
//
// EOF
// 




