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
// The shaping function of the VA1_ALICE is given by 
//
//      f(x) = A(1 - exp(-Bx))
//
// Where A is a normalization constant, tuned so that the integral 
//
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

#ifndef ROOT_TTree
# include <TTree.h>
#endif
#ifndef ROOT_TRandom
# include <TRandom.h>
#endif
#ifndef ALILOG_H
# include "AliLog.h"
#endif
#ifndef ALIFMDDIGITIZER_H
# include "AliFMDDigitizer.h"
#endif
#ifndef ALIFMD_H
# include "AliFMD.h"
#endif
#ifndef ALIFMDHIT_H
# include "AliFMDHit.h"
#endif
#ifndef ALIFMDDIGIT_H
# include "AliFMDDigit.h"
#endif
#ifndef ALIFMDDIGIT_H
# include "AliFMDSDigit.h"
#endif
#ifndef ALIRUNDIGITIZER_H
# include "AliRunDigitizer.h"
#endif
#ifndef ALIRUN_H
# include "AliRun.h"
#endif
#ifndef ALILOADER_H
# include "AliLoader.h"
#endif
#ifndef ALIRUNLOADER_H
# include "AliRunLoader.h"
#endif
    
//____________________________________________________________________
ClassImp(AliFMDEdepMap);

//====================================================================
ClassImp(AliFMDBaseDigitizer);

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
    fEdep(kMaxDetectors, kMaxRings, kMaxSectors, kMaxStrips)
{
  // Normal CTOR
  AliDebug(1," processed");
  SetVA1MipRange();
  SetAltroChannelSize();
  SetSampleRate();
  SetShapingTime();
}

//____________________________________________________________________
AliFMDBaseDigitizer::AliFMDBaseDigitizer(const Char_t* name, 
					 const Char_t* title) 
  : AliDigitizer(name, title),
    fRunLoader(0),
    fEdep(kMaxDetectors, kMaxRings, kMaxSectors, kMaxStrips)
{
  // Normal CTOR
  AliDebug(1," processed");
  SetVA1MipRange();
  SetAltroChannelSize();
  SetSampleRate();
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
  return kTRUE;
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
  fEdep.Clear();
  
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
  
  // Loop over the tracks in the 
  for (Int_t track = 0; track < ntracks; track++)  {
    // Read in entry number `track' 
    hitsBranch->GetEntry(track);
    
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
      if (fEdep(detector, ring, sector, strip).first != 0)
	AliDebug(1, Form("Double hit in %d%c(%d,%d)", 
			 detector, ring, sector, strip));
      
      fEdep(detector, ring, sector, strip).first  += edep;
      fEdep(detector, ring, sector, strip).second += 1;
      // Add this to the energy deposited for this strip
    }  // hit loop
  } // track loop
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::DigitizeHits(AliFMD* fmd) const
{
  // For the stored energy contributions in the cache (fEdep), convert
  // the energy signal to ADC counts, and store the created digit in
  // the digits array (AliFMD::fDigits)
  //
  TArrayI counts(3);
  for (UShort_t detector=1; detector <= 3; detector++) {
    // Get pointer to subdetector 
    AliFMDSubDetector* det = 0;
    switch (detector) {
    case 1: det = fmd->GetFMD1(); break;
    case 2: det = fmd->GetFMD2(); break;
    case 3: det = fmd->GetFMD3(); break;
    }
    if (!det) continue;
    for (UShort_t ringi = 0; ringi <= 1; ringi++) {
      // Get pointer to Ring
      AliFMDRing* r = 0;
      switch (ringi) {
      case 0: if (det->GetInner()) r = det->GetInner(); break;
      case 1: if (det->GetOuter()) r = det->GetOuter(); break;
      }
      if (!r) continue;
      
      // Get number of sectors 
      UShort_t nSectors = UShort_t(360. / r->GetTheta());
      // Loop over the number of sectors 
      for (UShort_t sector = 0; sector < nSectors; sector++) {
	// Get number of strips 
	UShort_t nStrips = r->GetNStrips();
	// Loop over the stips 
	for (UShort_t strip = 0; strip < nStrips; strip++) {
	  counts.Reset(-1);
	  Float_t edep = fEdep(detector, r->GetId(), sector, strip).first;
	  ConvertToCount(edep, r->GetSiThickness(), fmd->GetSiDensity(), 
			 counts);
	  AddDigit(fmd, detector, r->GetId(), sector, strip, 
		   edep, UShort_t(counts[0]), 
		   Short_t(counts[1]), Short_t(counts[2]));
#if 0
	  // This checks if the digit created will give the `right'
	  // number of particles when reconstructed, using a naiive
	  // approach.  It's here only as a quality check - nothing
	  // else. 
	  CheckDigit(fEdep(detector, r->GetId(), sector, strip).first,
		     fEdep(detector, r->GetId(), sector, strip).second,
		     counts);
#endif
	} // Strip
      } // Sector 
    } // Ring 
  } // Detector 
}

//____________________________________________________________________
Float_t
AliFMDBaseDigitizer::ShapeIntegral(Float_t u, Float_t v) const
{
  // Calculates the integral 
  // 
  //      / v
  //      |               Bu + exp(-Bu) - Bv - exp(-Bv) 
  //      | f(x) dx = - A -----------------------------
  //      |                             B
  //      / u
  // 
  // of the shaping function of the VA1_ALICE between times u and v
  //
  //      f(x) = A(1 - exp(-Bx))
  //
  // where A is a normalization constant, tuned so that the integral 
  //
  //      /1
  //      |                              B	    
  //      | f(x) dx = 1   =>  A = ----------------
  //      |           	          -1 + B + exp(-B)
  //      / 0
  //
  // and B is the a parameter defined by the shaping time (fShapingTime).  
  // 
  // That is, the function return the value 
  // 
  //        Bu + exp(-Bu) - Bv - exp(-Bv) 
  //      - -----------------------------
  //               -1 + B + exp(-B)
  // 
  // u,v should lie in the interval [0,1], and u < v
  if (u == 0 && v == 1) return 1;
  Float_t B = fShapingTime;
  
  // Calculate the integral 
  Float_t res = - ((B * u + TMath::Exp(-B * u) - B * v - TMath::Exp(-B * v)) /
		 (-1 + B + TMath::Exp(-B)));
  return res;
}

//____________________________________________________________________
void
AliFMDBaseDigitizer::ConvertToCount(Float_t   edep, 
				    Float_t   siThickness, 
				    Float_t   siDensity, 
				    TArrayI&  counts) const
{
  // Put noise and make ADC signal
  // This is calculated as the product 
  // 
  //   DeltaEmip * SiThickness * SiDensity / Number 
  //
  // Where 
  //  
  //   DeltaEmip     is the energy loss of a MIP 
  //   SiThickness   is the thickness of the silicon 
  //   SiDensity     is the Silicon density 
  //   Number        is # of e^- per MIP
  //
  // Note: Need to check this is correct. 
  // 
  const Float_t mipI = 1.664 * siThickness * siDensity;
  // const Float_t mipI = 1.664 * 0.04 * 2.33 / 22400; // = 6.923e-6;
  
  // Create a pedestal 
  UShort_t ped = MakePedestal();
  
  Float_t convf = 1 / mipI * Float_t(fAltroChannelSize) / fVA1MipRange;
  Int_t n = fSampleRate;
  for (Ssiz_t i = 0; i < n;  i++) {
    Float_t w = ShapeIntegral(Float_t(i)/n, Float_t(i+1)/n);
    counts[i] = UShort_t(TMath::Min(w * edep * convf + ped, 
				    Float_t(fAltroChannelSize))); 
  }
}


//====================================================================
ClassImp(AliFMDDigitizer);

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
  SetPedestal();
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
    AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
  if (!fRunLoader) {
    AliError("Can not find Run Loader for input stream 0");
    return;
  }
  // Get the AliRun object 
  if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();

  // Get the AliFMD object 
  AliFMD* fmd = static_cast<AliFMD*>(fRunLoader->GetAliRun()->GetDetector("FMD"));
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
  digitTree->Fill();
  
  // Write the digits to disk 
  outFMD->WriteDigits("OVERWRITE");
  outFMD->UnloadHits();
  outFMD->UnloadDigits();

  // Reset the digits in the AliFMD object 
  fmd->ResetDigits();
}


//____________________________________________________________________
UShort_t
AliFMDDigitizer::MakePedestal() const 
{
  return UShort_t(TMath::Max(gRandom->Gaus(fPedestal, fPedestalWidth), 0.));
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
  fmd->AddDigit(detector, ring, sector, strip, count1, count2, count3);
}

//____________________________________________________________________
void
AliFMDDigitizer::CheckDigit(Float_t         /* edep */, 
			    UShort_t        nhits,
			    const TArrayI&  counts) 
{
  Int_t integral = counts[0];
  if (counts[1] >= 0) integral += counts[1];
  if (counts[2] >= 0) integral += counts[2];
  integral -= Int_t(fPedestal + 2 * fPedestalWidth);
  if (integral < 0) integral = 0;
  
  Float_t convf = Float_t(fVA1MipRange) / fAltroChannelSize;
  Float_t mips  = integral * convf;
  if (mips > Float_t(nhits) + .5 || mips < Float_t(nhits) - .5) 
    Warning("CheckDigit", "Digit -> %4.2f MIPS != %d +/- .5 hits", 
	    mips, nhits);
}

//====================================================================
ClassImp(AliFMDSDigitizer);

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
  fmd->AddSDigit(detector, ring, sector, strip, edep, count1, count2, count3);
}



//____________________________________________________________________
//
// EOF
// 




