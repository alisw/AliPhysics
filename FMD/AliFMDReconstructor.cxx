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
#include "AliFMDRecPoint.h"	   	   // ALIFMDMULTNAIIVE_H
#include "AliESD.h"			   // ALIESD_H
#include <AliESDFMD.h>			   // ALIESDFMD_H
#include <TFile.h>

//____________________________________________________________________
ClassImp(AliFMDReconstructor)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor() 
  : AliReconstructor(),
    fMult(0), 
    fESDObj(0)
{
  // Make a new FMD reconstructor object - default CTOR.  
}
  

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor(const AliFMDReconstructor& other) 
  : AliReconstructor(), 
    fMult(other.fMult),
    fESDObj(other.fESDObj)
{
  // Copy constructor 
}
  

//____________________________________________________________________
AliFMDReconstructor&
AliFMDReconstructor::operator=(const AliFMDReconstructor& other) 
{
  // Assignment operator
  fMult   = other.fMult;
  fESDObj = other.fESDObj;
  return *this;
}

//____________________________________________________________________
AliFMDReconstructor::~AliFMDReconstructor() 
{
  // Destructor 
  if (fMult)   fMult->Delete();
  if (fMult)   delete fMult;
  if (fESDObj) delete fESDObj;
}

//____________________________________________________________________
void 
AliFMDReconstructor::Init(AliRunLoader* runLoader) 
{
  // Initialize the reconstructor 
  AliDebug(1, Form("Init called with runloader 0x%x", runLoader));
  // Initialize the geometry 
  AliFMDGeometry* fmd = AliFMDGeometry::Instance();
  fmd->Init();
  fmd->InitTransformations();
  
  // Current vertex position
  fCurrentVertex = 0;
  // Create array of reconstructed strip multiplicities 
  fMult = new TClonesArray("AliFMDRecPoint", 51200);
  // Create ESD output object 
  fESDObj = new AliESDFMD;
  
  // Check that we have a run loader
  if (!runLoader) { 
    Warning("Init", "No run loader");
    return;
  }

  // Check if we can get the header tree 
  AliHeader* header = runLoader->GetHeader();
  if (!header) {
    Warning("Init", "No header");
    return;
  }

  // Check if we can get a simulation header 
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
}

//____________________________________________________________________
void 
AliFMDReconstructor::ConvertDigits(AliRawReader* reader, 
				   TTree* digitsTree) const
{
  // Convert Raw digits to AliFMDDigit's in a tree 
  AliDebug(1, "Reading raw data into digits tree");
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
#if 1
  if (fESD) {
    const AliESDVertex* vertex = fESD->GetVertex();
    if (vertex) {
      AliDebug(1, Form("Got vertex from ESD: %f", vertex->GetZv()));
      fCurrentVertex = vertex->GetZv();
    }
  }
#endif  
  TBranch *digitBranch = digitsTree->GetBranch("FMD");
  if (!digitBranch) {
    Error("Exec", "No digit branch for the FMD found");
    return;
  }
  TClonesArray* digits = new TClonesArray("AliFMDDigit");
  digitBranch->SetAddress(&digits);

  // TIter next(&fAlgorithms);
  // AliFMDMultAlgorithm* algorithm = 0;
  // while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) {
  //   AliDebug(10, Form("PreEvent called for algorithm %s", 
  //                     algorithm->GetName()));    
  //   algorithm->PreEvent(clusterTree, fCurrentVertex);
  // }
  if (fMult)   fMult->Clear();
  if (fESDObj) fESDObj->Clear();
  
  fNMult = 0;
  fTreeR = clusterTree;
  fTreeR->Branch("FMD", &fMult);
  
  AliDebug(5, "Getting entry 0 from digit branch");
  digitBranch->GetEntry(0);
  
  AliDebug(5, "Processing digits");
  ProcessDigits(digits);

  // next.Reset();
  // algorithm = 0;
  // while ((algorithm = static_cast<AliFMDMultAlgorithm*>(next()))) {
  //   AliDebug(10, Form("PostEvent called for algorithm %s", 
  //                     algorithm->GetName()));
  // algorithm->PostEvent();
  // }
  Int_t written = clusterTree->Fill();
  AliDebug(10, Form("Filled %d bytes into cluster tree", written));
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
    AliFMDParameters* param  = AliFMDParameters::Instance();
    // Check that the strip is not marked as dead 
    if (param->IsDead(digit->Detector(), digit->Ring(), 
		      digit->Sector(), digit->Strip())) continue;

    // digit->Print();
    // Get eta and phi 
    Float_t eta, phi;
    PhysicalCoordinates(digit, eta, phi);
    
    // Substract pedestal. 
    UShort_t counts   = SubtractPedestal(digit);
    
    // Gain match digits. 
    Double_t edep     = Adc2Energy(digit, eta, counts);
    
    // Make rough multiplicity 
    Double_t mult     = Energy2Multiplicity(digit, edep);
    
    AliDebug(10, Form("FMD%d%c[%2d,%3d]: "
		      "ADC: %d, Counts: %d, Energy: %f, Mult: %f", 
		      digit->Detector(), digit->Ring(), digit->Sector(),
		      digit->Strip(), digit->Counts(), counts, edep, mult));
    
    // Create a `RecPoint' on the output branch. 
    AliFMDRecPoint* m = 
      new ((*fMult)[fNMult]) AliFMDRecPoint(digit->Detector(), 
					    digit->Ring(), 
					    digit->Sector(),
					    digit->Strip(),
					    eta, phi, 
					    edep, mult);
    (void)m; // Suppress warnings about unused variables. 
    fNMult++;

    fESDObj->SetMultiplicity(digit->Detector(), digit->Ring(), 
			     digit->Sector(),  digit->Strip(), mult);
    fESDObj->SetEta(digit->Detector(), digit->Ring(), 
		    digit->Sector(),  digit->Strip(), eta);
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

  Int_t             counts = 0;
  AliFMDParameters* param  = AliFMDParameters::Instance();
  Float_t           pedM   = param->GetPedestal(digit->Detector(), 
						digit->Ring(), 
						digit->Sector(), 
						digit->Strip());
  AliDebug(10, Form("Subtracting pedestal %f from signal %d", 
		   pedM, digit->Counts()));
  if (digit->Count3() > 0)      counts = digit->Count3();
  else if (digit->Count2() > 0) counts = digit->Count2();
  else                          counts = digit->Count1();
  counts = TMath::Max(Int_t(counts - pedM), 0);
  if (counts > 0) AliDebug(10, "Got a hit strip");
  
  return  UShort_t(counts);
}

//____________________________________________________________________
Float_t
AliFMDReconstructor::Adc2Energy(AliFMDDigit* digit, 
				Float_t      /* eta */, 
				UShort_t     count) const
{
  // Converts number of ADC counts to energy deposited. 
  // Note, that this member function can be overloaded by derived
  // classes to do strip-specific look-ups in databases or the like,
  // to find the proper gain for a strip. 
  // 
  // In this simple version, we calculate the energy deposited as 
  // 
  //    EnergyDeposited = cos(theta) * gain * count
  // 
  // where 
  // 
  //           Pre_amp_MIP_Range
  //    gain = ----------------- * Energy_deposited_per_MIP
  //           ADC_channel_size    
  // 
  // is constant and the same for all strips.

  // Double_t theta = 2 * TMath::ATan(TMath::Exp(-eta));
  // Double_t edep  = TMath::Abs(TMath::Cos(theta)) * fGain * count;
  AliFMDParameters* param = AliFMDParameters::Instance();
  Float_t           gain  = param->GetPulseGain(digit->Detector(), 
						digit->Ring(), 
						digit->Sector(), 
						digit->Strip());
  Double_t          edep  = count * gain;
  AliDebug(10, Form("Converting counts %d to energy via factor %f", 
		    count, gain));
  return edep;
}

//____________________________________________________________________
Float_t
AliFMDReconstructor::Energy2Multiplicity(AliFMDDigit* /* digit */, 
					 Float_t      edep) const
{
  // Converts an energy signal to number of particles. 
  // Note, that this member function can be overloaded by derived
  // classes to do strip-specific look-ups in databases or the like,
  // to find the proper gain for a strip. 
  // 
  // In this simple version, we calculate the multiplicity as 
  // 
  //   multiplicity = Energy_deposited / Energy_deposited_per_MIP
  // 
  // where 
  //
  //   Energy_deposited_per_MIP = 1.664 * SI_density * SI_thickness 
  // 
  // is constant and the same for all strips 
  AliFMDParameters* param   = AliFMDParameters::Instance();
  Double_t          edepMIP = param->GetEdepMip();
  Float_t           mult    = edep / edepMIP;
  if (edep > 0) 
    AliDebug(10, Form("Translating energy %f to multiplicity via "
		     "divider %f->%f", edep, edepMIP, mult));
  return mult;
}

//____________________________________________________________________
void
AliFMDReconstructor::PhysicalCoordinates(AliFMDDigit* digit, 
					 Float_t& eta, 
					 Float_t& phi) const
{
  // Get the eta and phi of a digit 
  // 
  // Get geometry. 
  AliFMDGeometry* fmd = AliFMDGeometry::Instance();
#if 0
  AliFMDDetector* subDetector = fmd->GetDetector(digit->Detector());
  if (!subDetector) { 
    Warning("ProcessDigits", "Unknown detector: FMD%d" , digit->Detector());
    return;
  }
    
  // Get the ring - we should really use geometry->Detector2XYZ(...)
  // here. 
  AliFMDRing* ring  = subDetector->GetRing(digit->Ring());
  Float_t     ringZ = subDetector->GetRingZ(digit->Ring());
  if (!ring) {
    Warning("ProcessDigits", "Unknown ring: FMD%d%c", digit->Detector(), 
	    digit->Ring());
    return;
  }
  Float_t  realZ    = fCurrentVertex + ringZ;
  Float_t  stripR   = ((ring->GetHighR() - ring->GetLowR()) 
		       / ring->GetNStrips() * (digit->Strip() + .5) 
		       + ring->GetLowR());
  Float_t  theta    = TMath::ATan2(stripR, realZ);
#endif    
  Double_t x, y, z, r, theta;
  fmd->Detector2XYZ(digit->Detector(), digit->Ring(), digit->Sector(), 
		    digit->Strip(), x, y, z);
  // Correct for vertex offset. 
  z     += fCurrentVertex;
  phi   =  TMath::ATan2(y, x);
  r     =  TMath::Sqrt(y * y + x * x);
  theta =  TMath::ATan2(r, z);
  eta   = -TMath::Log(TMath::Tan(theta / 2));
}

      

//____________________________________________________________________
void 
AliFMDReconstructor::FillESD(TTree*  /* digitsTree */, 
			     TTree*  /* clusterTree */,
			     AliESD* esd) const
{
  // nothing to be done
  // FIXME: The vertex may not be known when Reconstruct is executed,
  // so we may have to move some of that member function here. 
  AliDebug(1, Form("Calling FillESD with two trees and one ESD"));
  // fESDObj->Print();

  if (esd) { 
    AliDebug(1, Form("Writing FMD data to ESD tree"));
    esd->SetFMDData(fESDObj);
    // Let's check the data in the ESD
#if 0
    AliESDFMD* fromEsd = esd->GetFMDData();
    if (!fromEsd) {
      AliWarning("No FMD object in ESD!");
      return;
    }
    for (UShort_t det = 1; det <= fESDObj->MaxDetectors(); det++) {
      for (UShort_t ir = 0; ir < fESDObj->MaxRings(); ir++) {
	Char_t ring = (ir == 0 ? 'I' : 'O');
	for (UShort_t sec = 0; sec < fESDObj->MaxSectors(); sec++) {
	  for (UShort_t str = 0; str < fESDObj->MaxStrips(); str++) {
	    if (fESDObj->Multiplicity(det, ring, sec, str) != 
		fromEsd->Multiplicity(det, ring, sec, str))
	      AliWarning(Form("Mult for FMD%d%c[%2d,%3d]",det,ring,sec,str));
	    if (fESDObj->Eta(det, ring, sec, str) != 
		fromEsd->Eta(det, ring, sec, str))
	      AliWarning(Form("Eta for FMD%d%c[%2d,%3d]", det,ring,sec,str));
	    if (fESDObj->Multiplicity(det, ring, sec, str) > 0 && 
		fESDObj->Multiplicity(det, ring, sec, str) 
		!= AliESDFMD::kInvalidMult) 
	      AliInfo(Form("Mult in FMD%d%c[%2d,%3d] is %f", det,ring,sec,str,
			   fESDObj->Multiplicity(det, ring, sec, str)));
	  }
	}
      }
    }
#endif
  }

#if 0  
  static Int_t evNo = -1;
  evNo++;
  if (esd) evNo = esd->GetEventNumber();
  TString fname(Form("FMD.ESD.%03d.root", evNo));
  TFile* file = TFile::Open(fname.Data(), "RECREATE");
  if (!file) {
    AliError(Form("Failed to open file %s", fname.Data()));
    return;
  }
  fESDObj->Write();
  file->Close();
#endif
    
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
void 
AliFMDReconstructor::Reconstruct(AliRawReader*,TTree*) const 
{
  // Cannot be used.  See member function with same name but with 2
  // TTree arguments.   Make sure you do local reconstrucion 
  AliDebug(1, Form("Calling FillESD with loader and tree"));
  AliError("MayNotUse");
}
//____________________________________________________________________
void 
AliFMDReconstructor::Reconstruct(AliRunLoader*) const 
{
  // Cannot be used.  See member function with same name but with 2
  // TTree arguments.   Make sure you do local reconstrucion 
  AliDebug(1, Form("Calling FillESD with loader"));
  AliError("MayNotUse");
}
//____________________________________________________________________
void 
AliFMDReconstructor::Reconstruct(AliRunLoader*, AliRawReader*) const 
{
  // Cannot be used.  See member function with same name but with 2
  // TTree arguments.   Make sure you do local reconstrucion 
  AliDebug(1, Form("Calling FillESD with loader and raw reader"));
  AliError("MayNotUse");
}
//____________________________________________________________________
void 
AliFMDReconstructor::FillESD(AliRawReader*,TTree*,AliESD*) const 
{
  // Cannot be used.  See member function with same name but with 2
  // TTree arguments.   Make sure you do local reconstrucion 
  AliDebug(1, Form("Calling FillESD with raw reader, tree, and ESD"));
  AliError("MayNotUse");
}
//____________________________________________________________________
void 
AliFMDReconstructor::FillESD(AliRunLoader*,AliESD*) const
{
  // Cannot be used.  See member function with same name but with 2
  // TTree arguments.   Make sure you do local reconstrucion 
  AliDebug(1, Form("Calling FillESD with loader and ESD"));
  AliError("MayNotUse");
}
//____________________________________________________________________
void 
AliFMDReconstructor::FillESD(AliRunLoader*,AliRawReader*,AliESD*) const 
{
  // Cannot be used.  See member function with same name but with 2
  // TTree arguments.   Make sure you do local reconstrucion 
  AliDebug(1, Form("Calling FillESD with loader, raw reader, and ESD"));
  AliError("MayNotUse");
}

//____________________________________________________________________
//
// EOF
//
