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
//
// This is a class that constructs AliFMDMult (reconstructed
// multiplicity) from of Digits
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
//      +---------------------+       +---------------------+
//      | AliFMDReconstructor |<>-----| AliFMDMultAlgorithm |
//      +---------------------+       +---------------------+
//                                               ^
//                                               |
//                                   +-----------+---------+
//                                   |                     |
//                         +-------------------+   +------------------+
//                         | AliFMDMultPoisson |   | AliFMDMultNaiive |
//                         +-------------------+   +------------------+
//
// AliFMDReconstructor acts as a manager class.  It contains a list of
// AliFMDMultAlgorithm objects.  The call graph looks something like 
//
//
//       +----------------------+            +----------------------+
//       | :AliFMDReconstructor |            | :AliFMDMultAlgorithm |
//       +----------------------+            +----------------------+
//                  |                                  |
//    Reconstruct  +-+                                 |
//    ------------>| |                         PreRun +-+
//                 | |------------------------------->| |   
//                 | |                                +-+
//                 | |-----+ (for each event)          |
//                 | |     | *ProcessEvent             |
//                 |+-+    |                           |
//                 || |<---+                 PreEvent +-+
//                 || |------------------------------>| |      
//                 || |                               +-+
//                 || |-----+                          |
//                 || |     | ProcessDigits            |
//                 ||+-+    |                          |
//                 ||| |<---+                          |
//                 ||| |         *ProcessDigit(digit) +-+
//                 ||| |----------------------------->| |
//                 ||| |                              +-+
//                 ||+-+                               |
//                 || |                     PostEvent +-+
//                 || |------------------------------>| |
//                 || |                               +-+
//                 |+-+                                |
//                 | |                        PostRun +-+
//                 | |------------------------------->| |
//                 | |                                +-+
//                 +-+                                 |
//                  |                                  |
//
//
// 
//-- Authors: Evgeny Karpechev(INR) and Alla Maevsksia
//  Latest changes by Christian Holm Christensen <cholm@nbi.dk>
//
//
//____________________________________________________________________

#include <AliLog.h>                        // ALILOG_H
#include <AliRun.h>                        // ALIRUN_H
#include <AliRunLoader.h>                  // ALIRUNLOADER_H
#include <AliLoader.h>                     // ALILOADER_H
#include <AliHeader.h>                     // ALIHEADER_H
#include <AliRawReader.h>                  // ALIRAWREADER_H
#include <AliGenEventHeader.h>             // ALIGENEVENTHEADER_H
#include "AliFMD.h"                        // ALIFMD_H
#include "AliFMDDigit.h"                   // ALIFMDDIGIT_H
#include "AliFMDReconstructor.h"           // ALIFMDRECONSTRUCTOR_H
#include "AliFMDRawStream.h"               // ALIFMDRAWSTREAM_H
#include "AliFMDRawReader.h"               // ALIFMDRAWREADER_H
#include "AliFMDMultAlgorithm.h" 	   // ALIFMDMULTALGORITHM_H
#include "AliFMDMultPoisson.h"		   // ALIFMDMULTPOISSON_H
#include "AliFMDMultNaiive.h"		   // ALIFMDMULTNAIIVE_H

//____________________________________________________________________
ClassImp(AliFMDReconstructor);

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor() 
  : AliReconstructor(),
    fPedestal(0), 
    fPedestalWidth(0),
    fPedestalFactor(0)
{
  // Make a new FMD reconstructor object - default CTOR.
  SetPedestal();

  fFMDLoader = 0;
  fRunLoader = 0;
  fFMD       = 0;
  fAlgorithms.Add(new AliFMDMultNaiive);
  fAlgorithms.Add(new AliFMDMultPoisson);
}
  

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor(const AliFMDReconstructor& other) 
  : AliReconstructor(),
    fPedestal(0), 
    fPedestalWidth(0),
    fPedestalFactor(0)
{
  // Copy constructor 
  SetPedestal(other.fPedestal, other.fPedestalWidth, other.fPedestalFactor);

  fFMDLoader = other.fFMDLoader;
  fRunLoader = other.fRunLoader;
  fFMD       = other.fFMD;
  
  fAlgorithms.Delete();
  TIter next(&(other.fAlgorithms));
  AliFMDMultAlgorithm* algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    fAlgorithms.Add(algorithm);
  fAlgorithms.SetOwner(kFALSE);
}
  

//____________________________________________________________________
AliFMDReconstructor&
AliFMDReconstructor::operator=(const AliFMDReconstructor& other) 
{
  // Assignment operator
  SetPedestal(other.fPedestal, other.fPedestalWidth, other.fPedestalFactor);

  fFMDLoader = other.fFMDLoader;
  fRunLoader = other.fRunLoader;
  fFMD       = other.fFMD;

  fAlgorithms.Delete();
  TIter next(&(other.fAlgorithms));
  AliFMDMultAlgorithm* algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    fAlgorithms.Add(algorithm);
  fAlgorithms.SetOwner(kFALSE);

  return *this;
}

//____________________________________________________________________
AliFMDReconstructor::~AliFMDReconstructor() 
{
  // Destructor 
  fAlgorithms.Delete();
}
  
//____________________________________________________________________
void 
AliFMDReconstructor::SetPedestal(Float_t mean, Float_t width, Float_t factor) 
{
  // Set the pedestal, and pedestal width 
  fPedestal       = mean;
  fPedestalWidth  = width;
  fPedestalFactor = factor;
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

  TIter next(&fAlgorithms);
  AliFMDMultAlgorithm* algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    algorithm->PreRun(fFMD);

  if (rawReader) {
    Int_t event = 0;
    while (rawReader->NextEvent()) {
      ProcessEvent(event, rawReader);
      event++;
    }
  }
  else {
    Int_t nEvents= Int_t(fRunLoader->TreeE()->GetEntries()); 
    for(Int_t event = 0; event < nEvents; event++) 
      ProcessEvent(event, 0);
  }

  next.Reset();
  algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    algorithm->PostRun();

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
				  AliRawReader* reader) const
{
  // Process one event read from either a clones array or from a a raw
  // data reader. 
  fRunLoader->GetEvent(event) ;
  //event z-vertex for correction eta-rad dependence      
  AliHeader *header            = fRunLoader->GetHeader();
  if (!header) Warning("ProcessEvent", "no AliHeader found!");
  AliGenEventHeader* genHeader = (header ? header->GenEventHeader() : 0);

  // Get the Z--coordinate from the event header 
  TArrayF o(3); 
  if (genHeader) genHeader->PrimaryVertex(o);
  else           Warning("ProcessEvent", "No GenEventHeader Found");
  fCurrentVertex = o.At(2);

  // If the recontruction tree isn't loaded, load it
  if(fFMDLoader->TreeR()==0) fFMDLoader->MakeTree("R");
  
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
  TClonesArray* digits = fFMD->Digits();
  if (!digitBranch) {
    if (!reader) {
      Error("Exec", "No digit branch for the FMD found");
      return;
    }
    fFMD->MakeBranchInTree(digitTree, fFMD->GetName(), &(digits), 4000, 0);
  }
  if (!reader) digitBranch->SetAddress(&digits);

  if  (reader) {
    AliFMDRawReader rawRead(fFMD, reader);
    AliDebug(10, Form("Making raw reader with sample rate: %d",
		      fFMD->GetSampleRate()));
    rawRead.SetSampleRate(fFMD->GetSampleRate());
    rawRead.Exec();
  }
  else {
    // Read the ADC values from a clones array. 
    AliDebug(10, "Reading ADCs from Digits array");
    // read Digits, and reconstruct the particles
    if (!fFMDLoader->TreeD()->GetEvent(0)) return;
  }
  
  TIter next(&fAlgorithms);
  AliFMDMultAlgorithm* algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    algorithm->PreEvent(fFMDLoader->TreeR(), fCurrentVertex);

  ProcessDigits(digits);

  next.Reset();
  algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    algorithm->PostEvent();
  
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
UShort_t
AliFMDReconstructor::SubtractPedestal(AliFMDDigit* digit) const
{
  // Member function to subtract the pedestal from a digit
  // This implementation does nothing, but a derived class could over
  // load this to subtract a pedestal that was given in a database or
  // something like that. 

  Int_t counts = 0;
  Float_t ped = fPedestal + fPedestalFactor * fPedestalWidth;
  if (digit->Count3() > 0)      counts = digit->Count3();
  else if (digit->Count2() > 0) counts = digit->Count2();
  else                          counts = digit->Count1();
  counts = TMath::Max(Int_t(counts - ped), 0);
  return  UShort_t(counts);
}

//____________________________________________________________________
void
AliFMDReconstructor::ProcessDigits(TClonesArray* digits) const
{
  Int_t nDigits = digits->GetEntries();
  for (Int_t i = 0; i < nDigits; i++) {
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    AliFMDSubDetector* subDetector = 0;
    switch (digit->Detector()) {
    case 1: subDetector = fFMD->GetFMD1(); break;
    case 2: subDetector = fFMD->GetFMD2(); break;
    case 3: subDetector = fFMD->GetFMD3(); break;
    }
    if (!subDetector) { 
      Warning("ProcessDigits", "Unknown detector: FMD%d" , digit->Detector());
      continue;
    }
    
    AliFMDRing* ring  = 0;
    Float_t     ringZ = 0;
    switch(digit->Ring()) {
    case 'i':
    case 'I':  
      ring  = subDetector->GetInner(); 
      ringZ = subDetector->GetInnerZ();
      break;
    case 'o':
    case 'O':  
      ring  = subDetector->GetOuter(); 
      ringZ = subDetector->GetOuterZ();
      break;
    }
    if (!ring) {
      Warning("ProcessDigits", "Unknown ring: FMD%d%c", digit->Detector(), 
	      digit->Ring());
      break;
    }
    
    Float_t  realZ    = fCurrentVertex + ringZ;
    Float_t  stripR   = ((ring->GetHighR() - ring->GetLowR()) 
			 / ring->GetNStrips() * (digit->Strip() + .5) 
			 + ring->GetLowR());
    Float_t  theta    = TMath::ATan2(stripR, realZ);
    Float_t  phi      = (2 * TMath::Pi() / ring->GetNSectors() 
			 * (digit->Sector() + .5));
    Float_t  eta      = -TMath::Log(TMath::Tan(theta / 2));
    UShort_t counts   = SubtractPedestal(digit);
    
    TIter next(&fAlgorithms);
    AliFMDMultAlgorithm* algorithm = 0;
    while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
      algorithm->ProcessDigit(digit, eta, phi, counts);
  }
}
      
 
//____________________________________________________________________
void 
AliFMDReconstructor::FillESD(AliRunLoader* /*fRunLoader*/, 
			     AliESD* /*esd*/) const
{
  // nothing to be done

}

//____________________________________________________________________
//
// EOF
//
