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
/** @file    AliFMDReconstructor.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:47:09 2006
    @brief   FMD reconstruction 
*/
//____________________________________________________________________
//
// This is a class that constructs AliFMDRecPoint objects from of Digits
// This class reads either digits from a TClonesArray or raw data from 
// a DDL file (or similar), and stores the read ADC counts in an
// internal cache (fAdcs).   The rec-points are made via the naiive
// method. 
//
//-- Authors: Evgeny Karpechev(INR) and Alla Maevsksia
//  Latest changes by Christian Holm Christensen <cholm@nbi.dk>
//
//
//____________________________________________________________________

#include <AliLog.h>                        // ALILOG_H
// #include <AliRun.h>                        // ALIRUN_H
#include <AliRunLoader.h>                  // ALIRUNLOADER_H
#include <AliHeader.h>                     // ALIHEADER_H
#include <AliGenEventHeader.h>             // ALIGENEVENTHEADER_H
#include "AliFMDGeometry.h"                // ALIFMDGEOMETRY_H
#include "AliFMDParameters.h"              // ALIFMDPARAMETERS_H
#include "AliFMDDigit.h"                   // ALIFMDDIGIT_H
#include "AliFMDReconstructor.h"           // ALIFMDRECONSTRUCTOR_H
#include "AliFMDRawReader.h"               // ALIFMDRAWREADER_H
#include "AliFMDRecPoint.h"	   	   // ALIFMDMULTNAIIVE_H
#include "AliESD.h"			   // ALIESD_H
#include <AliESDFMD.h>			   // ALIESDFMD_H
class AliRawReader;

//____________________________________________________________________
ClassImp(AliFMDReconstructor)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor() 
  : AliReconstructor(),
    fMult(0x0), 
    fNMult(0),
    fTreeR(0x0),
    fCurrentVertex(0),
    fESDObj(0x0),
    fESD(0x0)
{
  // Make a new FMD reconstructor object - default CTOR.  
  SetNoiseFactor();
  SetAngleCorrect();
}
  

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor(const AliFMDReconstructor& other) 
  : AliReconstructor(), 
    fMult(other.fMult),
    fNMult(other.fNMult),
    fTreeR(other.fTreeR),
    fCurrentVertex(other.fCurrentVertex),
    fESDObj(other.fESDObj),
    fESD(other.fESD), 
    fNoiseFactor(other.fNoiseFactor),
    fAngleCorrect(other.fAngleCorrect)
{
  // Copy constructor 
}
  

//____________________________________________________________________
AliFMDReconstructor&
AliFMDReconstructor::operator=(const AliFMDReconstructor& other) 
{
  // Assignment operator
  fMult          = other.fMult;
  fNMult         = other.fNMult;
  fTreeR         = other.fTreeR;
  fCurrentVertex = other.fCurrentVertex;
  fESDObj        = other.fESDObj;
  fESD           = other.fESD;
  fNoiseFactor   = other.fNoiseFactor;
  fAngleCorrect  = other.fAngleCorrect;
   
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
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();

  // Initialize the parameters
  AliFMDParameters* param = AliFMDParameters::Instance();
  param->Init();
  
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

  if (fMult)   fMult->Clear();
  if (fESDObj) fESDObj->Clear();
  
  fNMult = 0;
  fTreeR = clusterTree;
  fTreeR->Branch("FMD", &fMult);
  
  AliDebug(5, "Getting entry 0 from digit branch");
  digitBranch->GetEntry(0);
  
  AliDebug(5, "Processing digits");
  ProcessDigits(digits);

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
  fESDObj->SetNoiseFactor(fNoiseFactor);
  fESDObj->SetAngleCorrected(fAngleCorrect);
  for (Int_t i = 0; i < nDigits; i++) {
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    AliFMDParameters* param  = AliFMDParameters::Instance();
    // Check that the strip is not marked as dead 
    if (param->IsDead(digit->Detector(), digit->Ring(), 
		      digit->Sector(), digit->Strip())) {
      AliDebug(10, Form("FMD%d%c[%2d,%3d] is dead", digit->Detector(), 
			digit->Ring(), digit->Sector(), digit->Strip()));
      continue;
    }

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
  Float_t           ped    = param->GetPedestal(digit->Detector(), 
						digit->Ring(), 
						digit->Sector(), 
						digit->Strip());
  Float_t           noise  = param->GetPedestalWidth(digit->Detector(), 
						     digit->Ring(), 
						     digit->Sector(), 
						     digit->Strip());
  AliDebug(15, Form("Subtracting pedestal %f from signal %d", 
		   ped, digit->Counts()));
  if (digit->Count3() > 0)      counts = digit->Count3();
  else if (digit->Count2() > 0) counts = digit->Count2();
  else                          counts = digit->Count1();
  counts = TMath::Max(Int_t(counts - ped), 0);
  if (counts < noise * fNoiseFactor) counts = 0;
  if (counts > 0) AliDebug(15, "Got a hit strip");
  
  return  UShort_t(counts);
}

//____________________________________________________________________
Float_t
AliFMDReconstructor::Adc2Energy(AliFMDDigit* digit, 
				Float_t      eta, 
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

  if (count <= 0) return 0;
  AliFMDParameters* param = AliFMDParameters::Instance();
  Float_t           gain  = param->GetPulseGain(digit->Detector(), 
						digit->Ring(), 
						digit->Sector(), 
						digit->Strip());
  AliDebug(15, Form("Converting counts %d to energy via factor %f", 
		    count, gain));

  Double_t edep  = count * gain;
  if (fAngleCorrect) {
    Double_t theta =  2 * TMath::ATan(TMath::Exp(-eta));
    edep           *= TMath::Abs(TMath::Cos(theta));
  }
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
    AliDebug(15, Form("Translating energy %f to multiplicity via "
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
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  Double_t x, y, z, r, theta;
  geom->Detector2XYZ(digit->Detector(), digit->Ring(), digit->Sector(), 
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
  }
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
