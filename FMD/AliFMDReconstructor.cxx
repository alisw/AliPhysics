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
// This is a class that constructs ReconstParticles (reconstructed
// particles) out of Digits
//
//
// This class reads either digits from a TClonesArray or raw data from
// a DDL file (or similar), and stores the read ADC counts in an
// internal cache (fAdcs). 
//
// From the cached values it then calculates the number of particles
// that hit a region of the FMDs, as specified by the user. 
//
// The reconstruction can be done in two ways: Either via counting the
// number of empty strips (Poisson method), or by converting the ADC
// signal to an energy deposition, and then dividing by the typical
// energy loss of a particle.
// 
// Currently, this class only reads the digits from a TClonesArray,
// and the Poission method for reconstruction. 
//
// 
//-- Authors: Evgeny Karpechev(INR) and Alla Maevsksia
//  Latest changes by Christian Holm Christensen <cholm@nbi.dk>
//
//
//____________________________________________________________________

#include "AliFMD.h"
#include "AliFMDDigit.h"
#include "AliFMDParticles.h"
#include "AliFMDReconstructor.h"
#include "AliAltroBuffer.h"
#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliFMDRawStream.h"
#include "AliRawReader.h"

//____________________________________________________________________
ClassImp(AliFMDReconstructor);

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor() 
  : AliReconstructor(),
    fAdcs(kMaxDetectors, kMaxRings, kMaxSectors, kMaxStrips),
    fDeltaEta(0), 
    fDeltaPhi(0), 
    fThreshold(0),
    fPedestal(0), 
    fPedestalWidth(0)
{
  SetDeltaEta();
  SetDeltaPhi();
  SetThreshold();
  SetPedestal();

  fParticles = new TClonesArray("AliFMDParticles", 1000);
  fFMDLoader = 0;
  fRunLoader = 0;
  fFMD       = 0;
}
  
//____________________________________________________________________
void 
AliFMDReconstructor::SetPedestal(Float_t mean, Float_t width) 
{
  // Set the pedestal, and pedestal width 
  fPedestal      = mean;
  fPedestalWidth = width;
}

//____________________________________________________________________
void 
AliFMDReconstructor::Reconstruct(AliRunLoader* runLoader, 
				 AliRawReader* rawReader) const
{ 
  // Collects all digits in the same active volume into number of
  // particles
  //
  // Reconstruct number of particles in given group of pads for given
  // FMDvolume determined by numberOfVolume,
  // numberOfMinSector, numberOfMaxSector, numberOfMinRing,
  // numberOgMaxRing 
  //
  // The reconstruction method is choosen based on the number of empty
  // strips. 
  fParticles->Clear();
  if (!runLoader) {
    Error("Exec","Run Loader loader is NULL - Session not opened");
    return;
  }
  fRunLoader = runLoader;
  fFMDLoader = runLoader->GetLoader("FMDLoader");
  if (!fFMDLoader) 
    Fatal("AliFMDReconstructor","Can not find FMD (loader) "
	  "in specified event");

  // Get the AliRun object
  if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
  
  // Get the AliFMD object
  fFMD = static_cast<AliFMD*>(fRunLoader->GetAliRun()->GetDetector("FMD"));
  if (!fFMD) {
    AliError("Can not get FMD from gAlice");
    return;
  }
  fFMDLoader->LoadRecPoints("RECREATE");

  if (!fRunLoader->TreeE())     fRunLoader->LoadHeader();

  TClonesArray* digits = fFMD->Digits();
  if (rawReader) {
    Int_t event = 0;
    while (rawReader->NextEvent()) {
      ProcessEvent(event, rawReader, digits);
      event++;
    }
  }
  else {
    Int_t nEvents= Int_t(fRunLoader->TreeE()->GetEntries()); 
    for(Int_t event = 0; event < nEvents; event++) 
      ProcessEvent(event, 0, digits);
  }


  fFMDLoader->UnloadRecPoints();
  fFMDLoader = 0;
  fRunLoader = 0;
  fFMD       = 0;
}

//____________________________________________________________________
void 
AliFMDReconstructor::Reconstruct(AliRunLoader* runLoader) const
{ 
  // Collects all digits in the same active volume into number of
  // particles
  //
  // Reconstruct number of particles in given group of pads for given
  // FMDvolume determined by numberOfVolume,
  // numberOfMinSector, numberOfMaxSector, numberOfMinRing,
  // numberOgMaxRing 
  //
  // The reconstruction method is choosen based on the number of empty
  // strips. 
  Reconstruct(runLoader, 0);
}


//____________________________________________________________________
void 
AliFMDReconstructor::ProcessEvent(Int_t event, 
				  AliRawReader* reader, 
				  TClonesArray* digits) const
{
  fRunLoader->GetEvent(event) ;
  //event z-vertex for correction eta-rad dependence      
  AliHeader *header            = fRunLoader->GetHeader();
  if (!header) 
    Warning("ProcessEvent", "no AliHeader found!");
  AliGenEventHeader* genHeader = (header ? header->GenEventHeader() : 0);

  // Get the Z--coordinate from the event header 
  TArrayF o(3); 
  if (genHeader) genHeader->PrimaryVertex(o);
  Float_t zVertex = o.At(2);

  // If the recontruction tree isn't loaded, load it
  if(fFMDLoader->TreeR()==0) fFMDLoader->MakeTree("R");
  
  //Make branches to hold the reconstructed particles 
  const Int_t kBufferSize = 16000;
  fFMDLoader->TreeR()->Branch("FMD", &fParticles, kBufferSize);

  // Load or recreate the digits 
  if (fFMDLoader->LoadDigits((reader ? "UPDATE" : "READ"))) {
    if (!reader) {
      Error("Exec","Error occured while loading digits. Exiting.");
      return;
    }
    
  }
  // Get the digits tree 
  TTree* digitTree = fFMDLoader->TreeD();
  if (!digitTree) { 
    if (!reader) {
      Error("Exec","Can not get Tree with Digits. "
	    "Nothing to reconstruct - Exiting");
      return;
    }
    fFMDLoader->MakeTree("D");
    digitTree = fFMDLoader->TreeD();
    
  }
  // Get the FMD branch holding the digits. 
  TBranch *digitBranch = digitTree->GetBranch("FMD");
  if (!digitBranch) {
    if (!reader) {
      Error("Exec", "No digit branch for the FMD found");
      return;
    }
    fFMD->MakeBranchInTree(digitTree, fFMD->GetName(), &(digits), 4000, 0);
  }
  if (!reader) digitBranch->SetAddress(&digits);

  fEmptyStrips = 0;
  fTotalStrips = 0;
  Bool_t ok = kFALSE;
  if      (reader) ok = ReadAdcs(reader);
  else if (digits) ok = ReadAdcs(digits);
  if (!ok) return;
  
  ReconstructFromCache(zVertex);

  if (reader) {
    digitTree->Fill();
    fFMDLoader->WriteDigits("OVERWRITE");
  }
  fFMDLoader->UnloadDigits();
  fFMDLoader->TreeR()->Reset();
  fFMDLoader->TreeR()->Fill(); 
  fFMDLoader->WriteRecPoints("OVERWRITE");
}

//____________________________________________________________________
Bool_t
AliFMDReconstructor::ReadAdcs(TClonesArray* digits) const
{
  AliDebug(10, "Reading ADCs from Digits array");
  // read Digits, and reconstruct the particles
  if (!fFMDLoader->TreeD()->GetEvent(0)) return kFALSE;

  // Reads the digits from the array, and fills up the cache (fAdcs) 
  fAdcs.Clear();
  Int_t nDigits = digits->GetEntries();
  for (Int_t digit = 0; digit < nDigits; digit++) {
    AliFMDDigit* fmdDigit = 
      static_cast<AliFMDDigit*>(digits->UncheckedAt(digit));    
    
    ProcessDigit(fmdDigit);
  } //digit loop
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDReconstructor::ReadAdcs(AliRawReader* reader) const
{
  AliDebug(10, "Reading ADCs from RawReader");
  // Reads the digits from a RAW data 
  fAdcs.Clear();
  // reader->Reset();

  if (!reader->ReadHeader()) {
    Error("ReadAdcs", "Couldn't read header");
    return kFALSE;
  }

  // Use AliAltroRawStream to read the ALTRO format.  No need to
  // reinvent the wheel :-) 
  AliFMDRawStream input(reader);
  // Select FMD DDL's 
  reader->Select(AliFMD::kBaseDDL >> 8);
  
  Int_t    oldDDL      = -1;
  Int_t    count       = 0;
  UShort_t detector    = 1; // Must be one here
  UShort_t oldDetector = 0;
  // Loop over data in file 
  Bool_t   next       = kTRUE;

  // local Cache 
  TArrayI counts(10);
  counts.Reset(-1);
  Int_t offset = 0;
  
  while (next) {
    next = input.Next();


    count++; 
    Int_t ddl = reader->GetDDLID();
    if (ddl != oldDDL 
	|| input.IsNewStrip()
	|| !next) {
      // Make a new digit, if we have some data! 
      if (counts[0] >= 0) {
	// Got a new strip. 
	AliDebug(10, Form("Add a new strip: FMD%d%c[%2d,%3d] "
			  "(current: FMD%d%c[%2d,%3d])", 
			  oldDetector, input.PrevRing(), 
			  input.PrevSector() , input.PrevStrip(),
			  detector , input.Ring(), input.Sector(), 
			  input.Strip()));
	fFMD->AddDigit(oldDetector, 
		       input.PrevRing(), 
		       input.PrevSector(), 
		       input.PrevStrip(), 
		       counts[0], counts[1], counts[2]);
	AliFMDDigit* digit = 
	  static_cast<AliFMDDigit*>(fFMD->Digits()->
				    UncheckedAt(fFMD->GetNdigits()-1));
	ProcessDigit(digit);
      }
	
      if (!next) { 
	AliDebug(10, Form("Read %d channels for FMD%d", 
			  count + 1, detector));
	break;
      }
    
    
      // If we got a new DDL, it means we have a new detector. 
      if (ddl != oldDDL) {
	if (detector != 0) 
	  AliDebug(10, Form("Read %d channels for FMD%d", 
			    count + 1, detector));
	// Reset counts, and update the DDL cache 
	count       = 0;
	oldDDL      = ddl;
	// Check that we're processing a FMD detector 
	Int_t detId = reader->GetDetectorID();
	if (detId != (AliFMD::kBaseDDL >> 8)) {
	  Error("ReadAdcs", "Detector ID %d != %d",
		detId, (AliFMD::kBaseDDL >> 8));
	  break;
	}
	// Figure out what detector we're deling with 
	oldDetector = detector;
	switch (ddl) {
	case 0: detector = 1; break;
	case 1: detector = 2; break;
	case 2: detector = 3; break;
	default:
	  Error("ReadAdcs", "Unknown DDL 0x%x for FMD", ddl);
	  return kFALSE;
	}
	AliDebug(10, Form("Reading ADCs for 0x%x  - That is FMD%d",
			  reader->GetEquipmentId(), detector));
      }
      counts.Reset(-1);
      offset = 0;
    }
    
    counts[offset] = input.Count();
    offset++;
    
    AliDebug(10, Form("ADC of FMD%d%c[%2d,%3d] += %d",
		      detector, input.Ring(), input.Sector(), 
		      input.Strip(), input.Count()));
    oldDetector = detector;
  }
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDReconstructor::ProcessDigit(AliFMDDigit* digit) const
{
  // Process a digit.  Derived classes can overload this member
  // function to do stuff to the digit.  However, it should write the
  // ADC count to the internal cache 
  //
  //    fAdcs(detector - 1, ring, sector, strip) = counts;
  //
  // In this implementation, we count the number of strips below
  // threshold.   This we do to later choose what kind of
  // reconstruction algorithm we'd like to use. 
  // 
  UShort_t detector = digit->Detector();
  Char_t   ring     = digit->Ring();
  UShort_t sector   = digit->Sector();
  UShort_t strip    = digit->Strip();    
  
  UShort_t counts  = SubtractPedestal(digit);
  
  fAdcs(detector - 1, ring, sector, strip) = counts;
  if (counts < fThreshold) fEmptyStrips++;
  fTotalStrips++;
}

//____________________________________________________________________
UShort_t
AliFMDReconstructor::SubtractPedestal(AliFMDDigit* digit) const
{
  // Member function to subtract the pedestal from a digit
  // This implementation does nothing, but a derived class could over
  // load this to subtract a pedestal that was given in a database or
  // something like that. 

  Int_t counts = 
    TMath::Max(Int_t(digit->Count1() - fPedestal - 3 * fPedestalWidth), 0);
  if (digit->Count2() >= 0) 
    counts += 
      TMath::Max(Int_t(digit->Count2() - fPedestal - 3 * fPedestalWidth), 0);
  if (digit->Count3() >= 0) 
    counts += 
      TMath::Max(Int_t(digit->Count3() - fPedestal - 3 * fPedestalWidth), 0);
      
  if (counts < 0) counts = 0;
  return  UShort_t(counts);
}

//____________________________________________________________________
void
AliFMDReconstructor::ReconstructFromCache(Float_t zVertex) const
{
  // Based on the information in the cache, do the reconstruction. 
  Int_t nRecon = 0;
  // Loop over the detectors 
  for (Int_t i = 1; i <= 3; i++) {
    AliFMDSubDetector* sub = 0;
    switch (i) {
    case 1: sub = fFMD->GetFMD1(); break;
    case 2: sub = fFMD->GetFMD2(); break;
    case 3: sub = fFMD->GetFMD3(); break;
    }
    if (!sub) continue;
	
    // Loop over the rings in the detector
    for (Int_t j = 0; j < 2; j++) {
      Float_t     rZ = 0;
      AliFMDRing* r  = 0;
      switch (j) {
      case 0: r  = sub->GetInner(); rZ = sub->GetInnerZ(); break;
      case 1: r  = sub->GetOuter(); rZ = sub->GetOuterZ(); break;
      }
      if (!r) continue;
      
      // Calculate low/high theta and eta 
      // FIXME: Is this right? 
      Float_t realZ    = zVertex + rZ;
      Float_t thetaOut = TMath::ATan2(r->GetHighR(), realZ);
      Float_t thetaIn  = TMath::ATan2(r->GetLowR(), realZ);
      Float_t etaOut   = - TMath::Log(TMath::Tan(thetaOut / 2));
      Float_t etaIn    = - TMath::Log(TMath::Tan(thetaIn / 2));
      if (TMath::Abs(etaOut) > TMath::Abs(etaIn)) {
	Float_t tmp = etaIn;
	etaIn       = etaOut;
	etaOut      = tmp;
      }

      //-------------------------------------------------------------
      //
      // Here starts poisson method 
      //
      // Calculate eta step per strip, number of eta steps, number of
      // phi steps, and check the sign of the eta increment 
      Float_t stripEta = (Float_t(r->GetNStrips()) / (etaIn - etaOut));
      Int_t   nEta     = Int_t(TMath::Abs(etaIn - etaOut) / fDeltaEta); 
      Int_t   nPhi     = Int_t(360. / fDeltaPhi);
      Float_t sign     = TMath::Sign(Float_t(1.), etaIn);

      AliDebug(10, Form("FMD%d%c Eta range: %f, %f %d Phi steps",
			sub->GetId(), r->GetId(), etaOut, etaIn, nPhi));

      // Loop over relevant phi values 
      for (Int_t p = 0; p < nPhi; p++) {
	Float_t  minPhi    = p * fDeltaPhi;
	Float_t  maxPhi    = minPhi + fDeltaPhi;
	UShort_t minSector = UShort_t(minPhi / 360) * r->GetNSectors();
	UShort_t maxSector = UShort_t(maxPhi / 360) * r->GetNSectors();
	
	AliDebug(10, Form(" Now in phi range %f, %f (sectors %d,%d)",
			  minPhi, maxPhi, minSector, maxSector));
	// Loop over relevant eta values 
	for (Int_t e = nEta; e >= 0; --e) {
	  Float_t  maxEta   = etaIn  - sign * e * fDeltaEta;
	  Float_t  minEta   = maxEta - sign * fDeltaEta;
	  if (sign > 0)  minEta = TMath::Max(minEta, etaOut);
	  else           minEta = TMath::Min(minEta, etaOut);
	  Float_t  theta1   = 2 * TMath::ATan(TMath::Exp(-minEta));
	  Float_t  theta2   = 2 * TMath::ATan(TMath::Exp(-maxEta));
	  Float_t  minR     = TMath::Abs(realZ * TMath::Tan(theta2));
	  Float_t  maxR     = TMath::Abs(realZ * TMath::Tan(theta1));
	  UShort_t minStrip = UShort_t((etaIn - maxEta) * stripEta + 0.5);
	  UShort_t maxStrip = UShort_t((etaIn - minEta) * stripEta + 0.5);

	  AliDebug(10, Form("  Now in eta range %f, %f (strips %d, %d)\n"
			    "    [radii %f, %f, thetas %f, %f, sign %d]", 
			    minEta, maxEta, minStrip, maxStrip, 
			    minR, maxR, theta1, theta2, sign));

	  // Count number of empty strips
	  Int_t   emptyStrips = 0;
	  for (Int_t sector = minSector; sector < maxSector; sector++) 
	    for (Int_t strip = minStrip; strip < maxStrip; strip++) 
	      if (fAdcs(sub->GetId() - 1, r->GetId(), sector, strip) 
		  < fThreshold) emptyStrips++;
	  
	  // The total number of strips 
	  Float_t nTotal = (maxSector - minSector) * (maxStrip - minStrip);
	  
	  // Log ratio of empty to total number of strips 
	  AliDebug(10, Form("Lambda= %d / %d = %f", 
			    emptyStrips, nTotal, 
			    Float_t(emptyStrips) / nTotal));
	  
	  Double_t lambda = (emptyStrips > 0 ? 
			     - TMath::Log(Double_t(emptyStrips) / nTotal) :
			     1);

	  // The reconstructed number of particles is then given by 
	  Int_t reconstructed = Int_t(lambda * nTotal + 0.5);
	    
	  // Add a AliFMDParticles to the reconstruction tree. 
	  new((*fParticles)[nRecon])   
	    AliFMDParticles(sub->GetId(), r->GetId(),
			    minSector, maxSector, minStrip, maxStrip,
			    minEta, maxEta, minPhi, maxPhi,
			    reconstructed, AliFMDParticles::kPoission);
	  nRecon++;
	} // phi 
      } // eta
    } // ring 
  } // detector 
}

 
//____________________________________________________________________
void 
AliFMDReconstructor::FillESD(AliRunLoader* /*fRunLoader*/, 
			     AliESD* /*esd*/) const
{
// nothing to be done

}

