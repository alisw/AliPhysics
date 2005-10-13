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
#include "AliFMD.h"         		   // ALIFMD_H
#include "AliFMDGeometry.h"                // ALIFMDGEOMETRY_H
#include "AliFMDParameters.h"              // ALIFMDPARAMETERS_H
#include "AliFMDDetector.h"                // ALIFMDDETECTOR_H
#include "AliFMDRing.h"                    // ALIFMDRING_H
#include "AliFMDDigit.h"                   // ALIFMDDIGIT_H
#include "AliFMDReconstructor.h"           // ALIFMDRECONSTRUCTOR_H
#include "AliFMDRawStream.h"               // ALIFMDRAWSTREAM_H
#include "AliFMDRawReader.h"               // ALIFMDRAWREADER_H
#include "AliFMDMultAlgorithm.h" 	   // ALIFMDMULTALGORITHM_H
#include "AliFMDMultPoisson.h"		   // ALIFMDMULTPOISSON_H
#include "AliFMDMultNaiive.h"		   // ALIFMDMULTNAIIVE_H
#include "AliESD.h"			   // ALIESD_H

//____________________________________________________________________
ClassImp(AliFMDReconstructor)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor() 
  : AliReconstructor(),
    fPedestal(0), 
    fPedestalWidth(0),
    fPedestalFactor(0)
{
  // Make a new FMD reconstructor object - default CTOR.
  AliFMDParameters* pars = AliFMDParameters::Instance();
  fPedestal       = pars->GetPedestal();
  fPedestalWidth  = pars->GetPedestalWidth();
  fPedestalFactor = pars->GetPedestalFactor();
  
  fAlgorithms.Add(new AliFMDMultNaiive);
  fAlgorithms.Add(new AliFMDMultPoisson);

  TIter next(&fAlgorithms);
  AliFMDMultAlgorithm* algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    algorithm->PreRun(0);
}
  

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor(const AliFMDReconstructor& other) 
  : AliReconstructor(),
    fPedestal(0), 
    fPedestalWidth(0),
    fPedestalFactor(0)
{
  // Copy constructor 
  fPedestal       = other.fPedestal;
  fPedestalWidth  = other.fPedestalWidth;
  fPedestalFactor = other.fPedestalFactor;

  
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
  fPedestal       = other.fPedestal;
  fPedestalWidth  = other.fPedestalWidth;
  fPedestalFactor = other.fPedestalFactor;

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
AliFMDReconstructor::Init(AliRunLoader* runLoader) 
{
  // Initialize the reconstructor 
  AliDebug(1, Form("Init called with runloader 0x%x", runLoader));
  fCurrentVertex = 0;
  if (!runLoader) { 
    Warning("Init", "No run loader");
    return;
  }
  AliHeader* header = runLoader->GetHeader();
  if (!header) {
    Warning("Init", "No header");
    return;
  }
  AliGenEventHeader* eventHeader = header->GenEventHeader();
  if (eventHeader) {
    TArrayF vtx;
    eventHeader->PrimaryVertex(vtx);
    fCurrentVertex = vtx[2];
    AliDebug(1, Form("Primary vertex Z coordinate for event # %d/%d is %f", 
		     header->GetRun(), header->GetEvent(), fCurrentVertex));
    Warning("Init", "no generator event header");
  }
  else {
    Warning("Init", "No generator event header - "
	    "perhaps we get the vertex from ESD?");
  }
  // Get the ESD tree 
  SetESD(new AliESD);
}

//____________________________________________________________________
void 
AliFMDReconstructor::ConvertDigits(AliRawReader* reader, 
				   TTree* digitsTree) const
{
  // Convert Raw digits to AliFMDDigit's in a tree 
  AliFMDRawReader rawRead(reader, digitsTree);
  // rawRead.SetSampleRate(fFMD->GetSampleRate());
  rawRead.Exec();
}

//____________________________________________________________________
void 
AliFMDReconstructor::Reconstruct(TTree* digitsTree, 
				 TTree* clusterTree) const 
{
  // Reconstruct event from digits in tree 
  // Get the FMD branch holding the digits. 
  // FIXME: The vertex may not be known yet, so we may have to move
  // some of this to FillESD. 
  AliDebug(1, "Reconstructing from digits in a tree");

  if (fESD) {
    const AliESDVertex* vertex = fESD->GetVertex();
    // if (vertex) {
    //   AliDebug(1, Form("Got vertex from ESD: %f", vertex->GetZv()));
    //   fCurrentVertex = vertex->GetZv();
    // }
  }
  
  TBranch *digitBranch = digitsTree->GetBranch("FMD");
  if (!digitBranch) {
    Error("Exec", "No digit branch for the FMD found");
    return;
  }
  TClonesArray* digits = new TClonesArray("AliFMDDigit");
  digitBranch->SetAddress(&digits);

  TIter next(&fAlgorithms);
  AliFMDMultAlgorithm* algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    algorithm->PreEvent(clusterTree, fCurrentVertex);
  digitBranch->GetEntry(0);
  
  ProcessDigits(digits);

  next.Reset();
  algorithm = 0;
  while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) 
    algorithm->PostEvent();
  clusterTree->Fill();
  digits->Delete();
  delete digits;
}
 
//____________________________________________________________________
void
AliFMDReconstructor::ProcessDigits(TClonesArray* digits) const
{
  // For each digit, find the pseudo rapdity, azimuthal angle, and
  // number of corrected ADC counts, and pass it on to the algorithms
  // used. 
  Int_t nDigits = digits->GetEntries();
  AliDebug(1, Form("Got %d digits", nDigits));
  for (Int_t i = 0; i < nDigits; i++) {
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    AliFMDGeometry* fmd = AliFMDGeometry::Instance();
    AliFMDDetector* subDetector = fmd->GetDetector(digit->Detector());
    if (!subDetector) { 
      Warning("ProcessDigits", "Unknown detector: FMD%d" , digit->Detector());
      continue;
    }
    
    AliFMDRing* ring  = subDetector->GetRing(digit->Ring());
    Float_t     ringZ = subDetector->GetRingZ(digit->Ring());
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
AliFMDReconstructor::FillESD(TTree*  /* digitsTree */, 
			     TTree*  /* clusterTree */,
			     AliESD* /* esd*/) const
{
  // nothing to be done
  // FIXME: The vertex may not be known when Reconstruct is executed,
  // so we may have to move some of that member function here. 
#if 0
  TClonesArray* multStrips  = 0;
  TClonesArray* multRegions = 0;
  TTree*        treeR  = fmdLoader->TreeR();
  TBranch*      branchRegions = treeR->GetBranch("FMDPoisson");
  TBranch*      branchStrips  = treeR->GetBranch("FMDNaiive");
  branchRegions->SetAddress(&multRegions);
  branchStrips->SetAddress(&multStrips);
  
  Int_t total = 0;
  Int_t nEntries  = clusterTree->GetEntries();
  for (Int_t entry = 0; entry < nEntries; entry++) {
    AliDebug(5, Form("Entry # %d in cluster tree", entry));
    treeR->GetEntry(entry);
    
    
    Int_t nMults = multRegions->GetLast();
    for (Int_t i = 0; i <= nMults; i++) {
      AliFMDMultRegion* multR =
	static_cast<AliFMDMultRegion*>(multRegions->UncheckedAt(i));
      Int_t nParticles=multR->Particles();
      if (i>=0 && i<=13)   hEtaPoissonI1->AddBinContent(i+1,nParticles);
      if (i>=14 && i<=27 ) hEtaPoissonI2->AddBinContent(i-13,nParticles);
      if (i>=28 && i<=33 );
      if (i>=34 && i<=47 ) hEtaPoissonI3->AddBinContent(48-i,nParticles);
      if (i>=48 && i<=53)  hEtaPoissonO3->AddBinContent(54-i,nParticles);
    }
  }
#endif   
}

//____________________________________________________________________
//
// EOF
//
