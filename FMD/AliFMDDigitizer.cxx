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

#include <TTree.h>		// ROOT_TTree
#include <TRandom.h>		// ROOT_TRandom
#include <AliLog.h>		// ALILOG_H
#include "AliFMDDigitizer.h"	// ALIFMDDIGITIZER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"	        // ALIFMDRING_H
#include "AliFMDHit.h"		// ALIFMDHIT_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDParameters.h"   // ALIFMDPARAMETERS_H
#include <AliRunDigitizer.h>	// ALIRUNDIGITIZER_H
#include <AliRun.h>		// ALIRUN_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRunLoader.h>	// ALIRUNLOADER_H
    
//====================================================================
ClassImp(AliFMDBaseDigitizer)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDBaseDigitizer::AliFMDBaseDigitizer()  
  : fRunLoader(0)
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
	  AliFMDMap::kMaxStrips)
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

      // Check if strip is `dead' 
      if (param->IsDead(detector, ring, sector, strip)) continue;
      
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
  Float_t  convF          = 1/param->GetPulseGain(detector,ring,sector,strip);
  UShort_t ped            = MakePedestal(detector,ring,sector,strip);
  UInt_t   maxAdc         = param->GetAltroChannelSize();
  UShort_t rate           = param->GetSampleRate(AliFMDParameters::kBaseDDL);
  UShort_t size           = param->GetAltroChannelSize();
  
  // In case we don't oversample, just return the end value. 
  if (rate == 1) {
    counts[0] = UShort_t(TMath::Min(edep * convF + ped, Float_t(size)));
    return;
  }
  
  // Create a pedestal 
  Float_t b = fShapingTime;
  for (Ssiz_t i = 0; i < rate;  i++) {
    Float_t t = Float_t(i) / rate;
    Float_t s = edep + (last - edep) * TMath::Exp(-b * t);
    counts[i] = UShort_t(TMath::Min(s * convF + ped, Float_t(maxAdc)));
  }
}


//====================================================================
ClassImp(AliFMDDigitizer)

//____________________________________________________________________
AliFMDDigitizer::AliFMDDigitizer()  
  : AliFMDBaseDigitizer()
{
  // Default ctor - don't use it
}

//____________________________________________________________________
AliFMDDigitizer::AliFMDDigitizer(AliRunDigitizer* manager) 
  : AliFMDBaseDigitizer(manager)
{
  // Normal CTOR
  AliDebug(1," processed");
}

//____________________________________________________________________
void
AliFMDDigitizer::Exec(Option_t*) 
{
  // Get the output manager 
  TString outFolder(fManager->GetOutputFolderName());
  AliRunLoader* out = 
    AliRunLoader::GetRunLoader(outFolder.Data());
  // Get the FMD output manager 
  AliLoader* outFMD = out->GetLoader("FMDLoader");

  // Get the input loader 
  TString inFolder(fManager->GetInputFolderName(0));
  fRunLoader = 
    AliRunLoader::GetRunLoader(inFolder.Data());
  if (!fRunLoader) {
    AliError("Can not find Run Loader for input stream 0");
    return;
  }
  // Get the AliRun object 
  if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();

  // Get the AliFMD object 
  AliFMD* fmd = 
    static_cast<AliFMD*>(fRunLoader->GetAliRun()->GetDetector("FMD"));
  if (!fmd) {
    AliError("Can not get FMD from gAlice");
    return;
  }

  Int_t nFiles= fManager->GetNinputs();
  for (Int_t inputFile = 0; inputFile < nFiles; inputFile++) {
    AliDebug(1,Form(" Digitizing event number %d",
		    fManager->GetOutputEventNr()));
    // Get the current loader 
    fRunLoader = 
      AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    if (!fRunLoader) Fatal("Exec", "no run loader");
    // Cache contriutions 
    SumContributions(fmd);
  }
  // Digitize the event 
  DigitizeHits(fmd);

  // Load digits from the tree 
  outFMD->LoadDigits("update");
  // Get the tree of digits 
  TTree* digitTree = outFMD->TreeD();
  if (!digitTree) {
    outFMD->MakeTree("D");
    digitTree = outFMD->TreeD();
  }
  digitTree->Reset();
  // Make a branch in the tree 
  TClonesArray* digits = fmd->Digits();
  fmd->MakeBranchInTree(digitTree, fmd->GetName(), &(digits), 4000, 0);
  // TBranch* digitBranch = digitTree->GetBranch(fmd->GetName());
  // Fill the tree 
  Int_t write = 0;
  write = digitTree->Fill();
  AliDebug(1, Form("Wrote %d bytes to digit tree", write));
  
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
			  Short_t  count3) const
{
  // Add a digit
  fmd->AddDigitByFields(detector, ring, sector, strip, count1, count2, count3);
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
  integral -= Int_t(mean + 2 * width);
  if (integral < 0) integral = 0;
  
  Float_t convF = Float_t(range) / size;
  Float_t mips  = integral * convF;
  if (mips > Float_t(nhits) + .5 || mips < Float_t(nhits) - .5) 
    Warning("CheckDigit", "Digit -> %4.2f MIPS != %d +/- .5 hits", 
	    mips, nhits);
}

//====================================================================
ClassImp(AliFMDSDigitizer)

//____________________________________________________________________
AliFMDSDigitizer::AliFMDSDigitizer()  
{
  // Default ctor - don't use it
}

//____________________________________________________________________
AliFMDSDigitizer::AliFMDSDigitizer(const Char_t* headerFile, 
				   const Char_t* /* sdigfile */)
  : AliFMDBaseDigitizer("FMDSDigitizer", "FMD SDigitizer")
{
  // Normal CTOR
  AliDebug(1," processed");

  fRunLoader = AliRunLoader::GetRunLoader(); // Open(headerFile);
  if (!fRunLoader) 
    Fatal("AliFMDSDigitizer", "cannot open session, header file '%s'",
	  headerFile);
  AliLoader* loader = fRunLoader->GetLoader("FMDLoader");
  if (!loader) 
    Fatal("AliFMDSDigitizer", "cannot find FMD loader in specified event");

  // Add task to tasks folder 
  loader->PostSDigitizer(this);

}

//____________________________________________________________________
AliFMDSDigitizer::~AliFMDSDigitizer() 
{
  // Destructor
  AliLoader* loader = fRunLoader->GetLoader("FMDLoader");
  loader->CleanSDigitizer();
}

//____________________________________________________________________
void
AliFMDSDigitizer::Exec(Option_t*) 
{
  // Get the output manager 
  if (!fRunLoader) {
    Error("Exec", "Run loader is not set");
    return;
  }
  if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
  if (!fRunLoader->TreeE())     fRunLoader->LoadHeader();
  
  AliLoader* fmdLoader = fRunLoader->GetLoader("FMDLoader");
  if (!fmdLoader) Fatal("Exec", "no FMD loader");
  
  // Get the AliFMD object 
  AliFMD* fmd = 
    static_cast<AliFMD*>(fRunLoader->GetAliRun()->GetDetector("FMD"));
  if (!fmd) {
    AliError("Can not get FMD from gAlice");
    return;
  }

  Int_t nEvents = Int_t(fRunLoader->TreeE()->GetEntries());
  for (Int_t event = 0; event < nEvents; event++) {
    AliDebug(1,Form(" Digitizing event number %d", event));
    // Get the current loader 
    fRunLoader->GetEvent(event);

    if (!fmdLoader->TreeS()) fmdLoader->MakeTree("S");
    // Make a branch
    fmd->MakeBranch("S");
    
    // Cache contriutions 
    SumContributions(fmd);

    // Digitize the event 
    DigitizeHits(fmd);

    fmdLoader->TreeS()->Reset();
    fmdLoader->TreeS()->Fill();
    fmdLoader->WriteSDigits("OVERWRITE");
  }
}

//____________________________________________________________________
void
AliFMDSDigitizer::AddDigit(AliFMD*  fmd,
			   UShort_t detector, 
			   Char_t   ring,
			   UShort_t sector, 
			   UShort_t strip, 
			   Float_t  edep, 
			   UShort_t count1, 
			   Short_t  count2, 
			   Short_t  count3) const
{
  // Add a summable digit
  fmd->AddSDigitByFields(detector, ring, sector, strip, edep, 
			 count1, count2, count3); 
}



//____________________________________________________________________
//
// EOF
// 




