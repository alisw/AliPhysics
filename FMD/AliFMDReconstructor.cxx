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

// #include <AliLog.h>                        // ALILOG_H
// #include <AliRun.h>                        // ALIRUN_H
#include "AliFMDDebug.h"
#include "AliFMDGeometry.h"                // ALIFMDGEOMETRY_H
#include "AliFMDParameters.h"              // ALIFMDPARAMETERS_H
#include "AliFMDAltroMapping.h"            // ALIFMDALTROMAPPING_H
#include "AliFMDDigit.h"                   // ALIFMDDIGIT_H
#include "AliFMDSDigit.h"                  // ALIFMDDIGIT_H
#include "AliFMDReconstructor.h"           // ALIFMDRECONSTRUCTOR_H
#include "AliFMDRawReader.h"               // ALIFMDRAWREADER_H
#include "AliFMDRecPoint.h"	   	   // ALIFMDMULTNAIIVE_H
#include "AliESDEvent.h"		   // ALIESDEVENT_H
#include "AliESDVertex.h"		   // ALIESDVERTEX_H
#include "AliESDTZERO.h"		   // ALIESDVERTEX_H
#include <AliESDFMD.h>			   // ALIESDFMD_H
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <climits>
// Import revertexer into a private namespace (to prevent conflicts) 
namespace { 
# include "AliFMDESDRevertexer.h"
}


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
    fNoiseFactor(0),
    fAngleCorrect(kTRUE),
    fVertexType(kNoVertex),
    fESD(0x0),
    fDiagnostics(kTRUE),
    fDiagStep1(0), 
    fDiagStep2(0),
    fDiagStep3(0),
    fDiagStep4(0),
    fDiagAll(0)
{
  // Make a new FMD reconstructor object - default CTOR.  
  SetNoiseFactor();
  SetAngleCorrect();
  if (AliDebugLevel() > 0) fDiagnostics = kTRUE;
  for(Int_t det = 1; det<=3; det++) {
    fZS[det-1]       = kFALSE;
    fZSFactor[det-1] = 0;
  }
}
  

//____________________________________________________________________
AliFMDReconstructor::AliFMDReconstructor(const AliFMDReconstructor& other) 
  : AliReconstructor(), 
    fMult(other.fMult),
    fNMult(other.fNMult),
    fTreeR(other.fTreeR),
    fCurrentVertex(other.fCurrentVertex),
    fESDObj(other.fESDObj),
    fNoiseFactor(other.fNoiseFactor),
    fAngleCorrect(other.fAngleCorrect),
    fVertexType(other.fVertexType),
    fESD(other.fESD),
    fDiagnostics(other.fDiagnostics),
    fDiagStep1(other.fDiagStep1), 
    fDiagStep2(other.fDiagStep2),
    fDiagStep3(other.fDiagStep3),
    fDiagStep4(other.fDiagStep4),
    fDiagAll(other.fDiagAll) 
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
  fNoiseFactor   = other.fNoiseFactor;
  fAngleCorrect  = other.fAngleCorrect;
  fVertexType    = other.fVertexType;
  fESD           = other.fESD;
  fDiagnostics   = other.fDiagnostics;
  fDiagStep1     = other.fDiagStep1;
  fDiagStep2     = other.fDiagStep2;
  fDiagStep3     = other.fDiagStep3;
  fDiagStep4     = other.fDiagStep4;
  fDiagAll       = other.fDiagAll;
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
AliFMDReconstructor::Init() 
{
  // Initialize the reconstructor 

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
  
  // Check if we need diagnostics histograms 
  if (!fDiagnostics) return;
  AliInfo("Making diagnostics histograms");
  fDiagStep1   = new TH2I("diagStep1", "Read ADC vs. Noise surpressed ADC",
			1024, -.5, 1023.5, 1024, -.5, 1023.5);
  fDiagStep1->SetDirectory(0);
  fDiagStep1->GetXaxis()->SetTitle("ADC (read)");
  fDiagStep1->GetYaxis()->SetTitle(Form("ADC (noise surpressed %4.f)", 
					fNoiseFactor));
  fDiagStep2  = new TH2F("diagStep2",  "ADC vs Edep deduced",
			1024, -.5, 1023.5, 100, 0, 2);
  fDiagStep2->SetDirectory(0);
  fDiagStep2->GetXaxis()->SetTitle("ADC (noise surpressed)");
  fDiagStep2->GetYaxis()->SetTitle("#Delta E [GeV]");
  fDiagStep3  = new TH2F("diagStep3",  "Edep vs Edep path corrected",
			100, 0., 2., 100, 0., 2.);
  fDiagStep3->SetDirectory(0);
  fDiagStep3->GetXaxis()->SetTitle("#Delta E [GeV]");
  fDiagStep3->GetYaxis()->SetTitle("#Delta E/#Delta x #times #delta x [GeV]");
  fDiagStep4  = new TH2F("diagStep4",  "Edep vs Multiplicity deduced", 
			100, 0., 2., 100, -.1, 19.9);
  fDiagStep4->SetDirectory(0);
  fDiagStep4->GetXaxis()->SetTitle("#Delta E/#Delta x #times #delta x [GeV]");
  fDiagStep4->GetYaxis()->SetTitle("Multiplicity");
  fDiagAll    = new TH2F("diagAll",    "Read ADC vs Multiplicity deduced", 
			 1024, -.5, 1023.5, 100, -.1, 19.9);
  fDiagAll->SetDirectory(0);
  fDiagAll->GetXaxis()->SetTitle("ADC (read)");
  fDiagAll->GetYaxis()->SetTitle("Multiplicity");
}

//____________________________________________________________________
void 
AliFMDReconstructor::ConvertDigits(AliRawReader* reader, 
				   TTree* digitsTree) const
{
  // Convert Raw digits to AliFMDDigit's in a tree 
  AliFMDDebug(1, ("Reading raw data into digits tree"));
  AliFMDRawReader rawRead(reader, digitsTree);
  // rawRead.SetSampleRate(fFMD->GetSampleRate());
  rawRead.Exec();
  AliFMDAltroMapping* map = AliFMDParameters::Instance()->GetAltroMap();
  for (size_t i = 1; i <= 3; i++) { 
    fZS[i]       = rawRead.IsZeroSuppressed(map->Detector2DDL(i));
    fZSFactor[i] = rawRead.NoiseFactor(map->Detector2DDL(i));
  }
}

//____________________________________________________________________
void 
AliFMDReconstructor::GetVertex(AliESDEvent* esd) const
{
  // Return the vertex to use. 
  // This is obtained from the ESD object. 
  // If not found, a warning is issued.
  fVertexType    = kNoVertex;
  fCurrentVertex = 0;
  if (!esd) return;
  
  const AliESDVertex* vertex = esd->GetPrimaryVertex();
  if (!vertex)        vertex = esd->GetPrimaryVertexSPD();
  if (!vertex)        vertex = esd->GetPrimaryVertexTPC();
  if (!vertex)        vertex = esd->GetVertex();

  if (vertex) {
    AliFMDDebug(2, ("Got %s (%s) from ESD: %f", 
		    vertex->GetName(), vertex->GetTitle(), vertex->GetZv()));
    fCurrentVertex = vertex->GetZv();
    fVertexType    = kESDVertex;
    return;
  }
  else if (esd->GetESDTZERO()) { 
    AliFMDDebug(2, ("Got primary vertex from T0: %f", esd->GetT0zVertex()));
    fCurrentVertex = esd->GetT0zVertex();
    fVertexType    = kESDVertex;
    return;
  }
  AliWarning("Didn't get any vertex from ESD or generator");
}
  

//____________________________________________________________________
void 
AliFMDReconstructor::Reconstruct(AliRawReader* reader, TTree*) const
{
  // Reconstruct directly from raw data (no intermediate output on
  // digit tree or rec point tree).  
  // 
  // Parameters: 
  //   reader	Raw event reader 
  //   ctree    Not used. 
  AliFMDRawReader rawReader(reader, 0);

  UShort_t det, sec, str, fac;
  Short_t  adc, oldDet = -1;
  Bool_t   zs;
  Char_t   rng;
    
  while (rawReader.NextSignal(det, rng, sec, str, adc, zs, fac)) { 
    if (det != oldDet) { 
      fZS[det-1]       = zs;
      fZSFactor[det-1] = fac;
      oldDet           = det;
    }
    ProcessSignal(det, rng, sec, str, adc);
  }
}

//____________________________________________________________________
void 
AliFMDReconstructor::Digitize(AliRawReader* reader, TClonesArray* sdigits) const
{
  // Reconstruct directly from raw data (no intermediate output on
  // digit tree or rec point tree).  
  // 
  // Parameters: 
  //   reader	Raw event reader 
  //   ctree    Not used. 
  AliFMDRawReader rawReader(reader, 0);

  UShort_t det, sec, str, sam, rat, fac;
  Short_t  adc, oldDet = -1;
  Bool_t   zs;
  Char_t   rng;
    
  while (rawReader.NextSample(det, rng, sec, str, sam, rat, adc, zs, fac)) { 
    if (!rawReader.SelectSample(sam, rat)) continue;
    if (det != oldDet) { 
      fZS[det-1]       = zs;
      fZSFactor[det-1] = fac;
      oldDet           = det;
    }
    DigitizeSignal(sdigits, det, rng, sec, str, sam, adc);
  }
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
  // 
  // Parameters: 
  //   digitsTree	Pointer to a tree containing digits 
  //   clusterTree	Pointer to output tree 
  // 
  AliFMDDebug(2, ("Reconstructing from digits in a tree"));
  GetVertex(fESD);
  
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
  
  AliFMDDebug(5, ("Getting entry 0 from digit branch"));
  digitBranch->GetEntry(0);
  
  AliFMDDebug(1, ("Processing digits"));
  ProcessDigits(digits);

  Int_t written = clusterTree->Fill();
  AliFMDDebug(10, ("Filled %d bytes into cluster tree", written));
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
  // 
  // Parameters: 
  //    digits	Array of digits
  // 
  Int_t nDigits = digits->GetEntries();
  AliFMDDebug(1, ("Got %d digits", nDigits));
  fESDObj->SetNoiseFactor(fNoiseFactor);
  fESDObj->SetAngleCorrected(fAngleCorrect);
  for (Int_t i = 0; i < nDigits; i++) {
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    if (!digit) continue;
    ProcessDigit(digit);
  }
}

//____________________________________________________________________
void
AliFMDReconstructor::ProcessDigit(AliFMDDigit* digit) const
{
  UShort_t det = digit->Detector();
  Char_t   rng = digit->Ring();
  UShort_t sec = digit->Sector();
  UShort_t str = digit->Strip();
  Short_t  adc = digit->Counts();
  
  ProcessSignal(det, rng, sec, str, adc);
}

//____________________________________________________________________
void
AliFMDReconstructor::ProcessSignal(UShort_t det, 
				   Char_t   rng, 
				   UShort_t sec, 
				   UShort_t str, 
				   Short_t  adc) const
{
  // Process the signal from a single strip 
  // 
  // Parameters: 
  //    det	Detector ID
  //    rng	Ring ID
  //    sec	Sector ID
  //    rng	Strip ID
  //    adc     ADC counts
  // 
  AliFMDParameters* param  = AliFMDParameters::Instance();
  // Check that the strip is not marked as dead 
  if (param->IsDead(det, rng, sec, str)) {
    AliFMDDebug(10, ("FMD%d%c[%2d,%3d] is dead", det, rng, sec, str));
    return;
  }
  
  // digit->Print();
  // Get eta and phi 
  Float_t eta, phi;
  PhysicalCoordinates(det, rng, sec, str, eta, phi);
    
  // Substract pedestal. 
  UShort_t counts   = SubtractPedestal(det, rng, sec, str, adc);
  if(counts == USHRT_MAX) return;
  
    // Gain match digits. 
  Double_t edep     = Adc2Energy(det, rng, sec, str, eta, counts);
  // Get rid of nonsense energy
  if(edep < 0)  return;
  
  // Make rough multiplicity 
  Double_t mult     = Energy2Multiplicity(det, rng, sec, str, edep);
  // Get rid of nonsense mult
  if (mult < 0)  return; 
  AliFMDDebug(5, ("FMD%d%c[%2d,%3d]: "
		    "ADC: %d, Counts: %d, Energy: %f, Mult: %f",
		  det, rng, sec, str, adc, counts, edep, mult));
  
  // Create a `RecPoint' on the output branch. 
  if (fMult) {
    AliFMDRecPoint* m = 
      new ((*fMult)[fNMult]) AliFMDRecPoint(det, rng, sec, str, 
					    eta, phi, edep, mult);
    (void)m; // Suppress warnings about unused variables. 
    fNMult++;
  }
  
  fESDObj->SetMultiplicity(det, rng, sec, str, mult);
  fESDObj->SetEta(det, rng, sec, str, eta);
  
  if (fDiagAll) fDiagAll->Fill(adc, mult);  

}

//____________________________________________________________________
void
AliFMDReconstructor::DigitizeSignal(TClonesArray* sdigits, 
				    UShort_t det, 
				    Char_t   rng, 
				    UShort_t sec, 
				    UShort_t str, 
				    UShort_t /* sam */,
				    Short_t  adc) const
{
  // Process the signal from a single strip 
  // 
  // Parameters: 
  //    det	Detector ID
  //    rng	Ring ID
  //    sec	Sector ID
  //    rng	Strip ID
  //    adc     ADC counts
  // 
  AliFMDParameters* param  = AliFMDParameters::Instance();
  // Check that the strip is not marked as dead 
  if (param->IsDead(det, rng, sec, str)) {
    AliFMDDebug(10, ("FMD%d%c[%2d,%3d] is dead", det, rng, sec, str));
    return;
  }
  
  // Substract pedestal. 
  UShort_t counts   = SubtractPedestal(det, rng, sec, str, adc);
  if(counts == USHRT_MAX || counts == 0) return;
  
    // Gain match digits. 
  Double_t edep     = Adc2Energy(det, rng, sec, str, counts);
  // Get rid of nonsense energy
  if(edep < 0)  return;

  Int_t n = sdigits->GetEntriesFast();
  // AliFMDSDigit* sdigit = 
  new ((*sdigits)[n]) 
    AliFMDSDigit(det, rng, sec, str, edep, counts, counts, counts, counts);
  // sdigit->SetCount(sam, counts);
}

//____________________________________________________________________
UShort_t
AliFMDReconstructor::SubtractPedestal(UShort_t det, 
				      Char_t   rng, 
				      UShort_t sec, 
				      UShort_t str, 
				      UShort_t adc, 
				      Float_t  noiseFactor,
				      Bool_t   zsEnabled, 
				      UShort_t zsNoiseFactor) const
{
  AliFMDParameters* param  = AliFMDParameters::Instance();
  Float_t           ped    = (zsEnabled ? 0 : 
				param->GetPedestal(det, rng, sec, str));
  Float_t           noise  = param->GetPedestalWidth(det, rng, sec, str);
  if(ped < 0 || noise < 0) { 
    AliWarningClass(Form("Invalid pedestal (%f) or noise (%f) "
			 "for FMD%d%c[%02d,%03d]", 
		    ped, noise, det, rng, sec, str));
    return USHRT_MAX;
  }
  AliDebugClass(15, Form("Subtracting pedestal for FMD%d%c[%2d,%3d]=%4d "
			 "(%s w/factor %d, noise factor %f, "
			 "pedestal %8.2f+/-%8.2f)",
			 det, rng, sec, str, adc, 
			 (zsEnabled ? "zs'ed" : "straight"), 
			 zsNoiseFactor, noiseFactor, ped, noise));

  Int_t counts = adc + Int_t(zsEnabled ? zsNoiseFactor * noise : - ped);
  counts =  TMath::Max(Int_t(counts), 0);
  // Calculate the noise factor for suppressing remenants of the noise
  // peak.  If we have done on-line zero suppression, we only check
  // for noise signals that are larger than the suppressed noise.  If
  // the noise factor used on line is larger than the factor used
  // here, we do not do this check at all.  
  // 
  // For example:
  //    Online factor  |  Read factor |  Result 
  //    ---------------+--------------+-------------------------------
  //           2       |      3       | Check if signal > 1 * noise
  //           3       |      3       | Check if signal > 0
  //           3       |      2       | Check if signal > 0
  //
  // In this way, we make sure that we do not suppress away too much
  // data, and that the read-factor is the most stringent cut. 
  Float_t nf = TMath::Max(0.F, noiseFactor - (zsEnabled ? zsNoiseFactor : 0));
  if (counts < noise * nf) counts = 0;
  if (counts > 0) AliDebugClass(15, "Got a hit strip");

  return counts;
}


//____________________________________________________________________
UShort_t
AliFMDReconstructor::SubtractPedestal(UShort_t det, 
				      Char_t   rng, 
				      UShort_t sec, 
				      UShort_t str, 
				      Short_t  adc) const
{
  // Member function to subtract the pedestal from a digit
  //
  // Parameters: 
  //    det	Detector ID
  //    rng	Ring ID
  //    sec	Sector ID
  //    rng	Strip ID
  //    adc     # of ADC counts
  // Return:
  //    Pedestal subtracted signal or USHRT_MAX in case of problems 
  //
  UShort_t counts = SubtractPedestal(det, rng, sec, str, adc, 
				     fNoiseFactor, fZS[det-1], 
				     fZSFactor[det-1]);
  if (fDiagStep1) fDiagStep1->Fill(adc, counts);
  
  return counts;
}

//____________________________________________________________________
Float_t
AliFMDReconstructor::Adc2Energy(UShort_t det, 
				Char_t   rng, 
				UShort_t sec, 
				UShort_t str, 
				UShort_t count) const
{
  // Converts number of ADC counts to energy deposited. 
  // Note, that this member function can be overloaded by derived
  // classes to do strip-specific look-ups in databases or the like,
  // to find the proper gain for a strip. 
  // 
  // In the first simple version, we calculate the energy deposited as 
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
  //
  // For the production we use the conversion measured in the NBI lab.
  // The total conversion is then:
  // 
  //    gain = ADC / DAC
  // 
  //                  EdepMip * count
  //      => energy = ----------------
  //                  gain * DACPerADC
  // 
  // Parameters: 
  //    det	Detector ID
  //    rng	Ring ID
  //    sec	Sector ID
  //    rng	Strip ID
  //    counts  Number of ADC counts over pedestal
  // Return 
  //    The energy deposited in a single strip, or -1 in case of problems
  //
  if (count <= 0) return 0;
  AliFMDParameters* param = AliFMDParameters::Instance();
  Float_t           gain  = param->GetPulseGain(det, rng, sec, str);
  // 'Tagging' bad gains as bad energy
  if (gain < 0) { 
    AliWarning(Form("Invalid gain (%f) for FMD%d%c[%02d,%03d]", 
		    gain, det, rng, sec, str));
    return -1;
  }
  AliFMDDebug(5, ("Converting counts %d to energy (factor=%f, DAC2MIP=%f)", 
		  count, gain,param->GetDACPerMIP()));

  Double_t edep  = ((count * param->GetEdepMip()) 
		    / (gain * param->GetDACPerMIP()));
  return edep;
}

//____________________________________________________________________
Float_t
AliFMDReconstructor::Adc2Energy(UShort_t det, 
				Char_t   rng, 
				UShort_t sec, 
				UShort_t str, 
				Float_t  eta, 
				UShort_t count) const
{
  // Converts number of ADC counts to energy deposited. 
  // Note, that this member function can be overloaded by derived
  // classes to do strip-specific look-ups in databases or the like,
  // to find the proper gain for a strip. 
  // 
  // In the first simple version, we calculate the energy deposited as 
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
  //
  // For the production we use the conversion measured in the NBI lab.
  // The total conversion is then:
  // 
  //    gain = ADC / DAC
  // 
  //                  EdepMip * count
  //      => energy = ----------------
  //                  gain * DACPerADC
  // 
  // Parameters: 
  //    det	Detector ID
  //    rng	Ring ID
  //    sec	Sector ID
  //    rng	Strip ID
  //    eta     Psuedo-rapidity
  //    counts  Number of ADC counts over pedestal
  // Return 
  //    The energy deposited in a single strip, or -1 in case of problems
  //
  Double_t edep = Adc2Energy(det, rng, sec, str, count);
  
  if (fDiagStep2) fDiagStep2->Fill(count, edep);  
  if (fAngleCorrect) {
    Double_t theta = 2 * TMath::ATan(TMath::Exp(-eta));
    Double_t corr  = TMath::Abs(TMath::Cos(theta));
    Double_t cedep = corr * edep;
    AliFMDDebug(10, ("correcting for path %f * %f = %f (eta=%f, theta=%f)",
		      edep, corr, cedep, eta, theta));
    if (fDiagStep3) fDiagStep3->Fill(edep, cedep);  
    edep = cedep;
  }
  return edep;
}

//____________________________________________________________________
Float_t
AliFMDReconstructor::Energy2Multiplicity(UShort_t /*det*/, 
					 Char_t   /*rng*/, 
					 UShort_t /*sec*/, 
					 UShort_t /*str*/, 
					 Float_t  edep) const
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
  //
  // Parameters: 
  //    det	Detector ID
  //    rng	Ring ID
  //    sec	Sector ID
  //    rng	Strip ID
  //    edep    Energy deposited in a single strip
  // Return 
  //    The "bare" multiplicity corresponding to the energy deposited
  AliFMDParameters* param   = AliFMDParameters::Instance();
  Double_t          edepMIP = param->GetEdepMip();
  Float_t           mult    = edep / edepMIP;
#if 0
  if (edep > 0) 
    AliFMDDebug(15, ("Translating energy %f to multiplicity via "
		     "divider %f->%f", edep, edepMIP, mult));
#endif
  if (fDiagStep4) fDiagStep4->Fill(edep, mult);  
  return mult;
}

//____________________________________________________________________
void
AliFMDReconstructor::PhysicalCoordinates(UShort_t det, 
					 Char_t   rng, 
					 UShort_t sec, 
					 UShort_t str, 
					 Float_t& eta, 
					 Float_t& phi) const
{
  // Get the eta and phi of a digit 
  // 
  // Parameters: 
  //    det	Detector ID
  //    rng	Ring ID
  //    sec	Sector ID
  //    rng	Strip ID
  //    eta	On return, contains the psuedo-rapidity of the strip
  //    phi     On return, contains the azimuthal angle of the strip
  // 
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  Double_t x, y, z, r, theta;
  geom->Detector2XYZ(det, rng, sec, str, x, y, z);
  // Correct for vertex offset. 
  z     -= fCurrentVertex;
  phi   =  TMath::ATan2(y, x);
  r     =  TMath::Sqrt(y * y + x * x);
  theta =  TMath::ATan2(r, z);
  eta   = -TMath::Log(TMath::Tan(theta / 2));
}

      

//____________________________________________________________________
void 
AliFMDReconstructor::FillESD(TTree*  /* digitsTree */, 
			     TTree*  /* clusterTree */,
			     AliESDEvent* esd) const
{
  // nothing to be done
  // FIXME: The vertex may not be known when Reconstruct is executed,
  // so we may have to move some of that member function here. 
  AliFMDDebug(2, ("Calling FillESD with two trees and one ESD"));
  // fESDObj->Print();

  Double_t oldVz = fCurrentVertex;
  GetVertex(esd);
  if (fVertexType != kNoVertex) { 
    AliFMDDebug(2, ("Revertexing the ESD data to vz=%f (was %f)",
		    fCurrentVertex, oldVz));
    AliFMDESDRevertexer revertexer;
    revertexer.Revertex(fESDObj, fCurrentVertex);
  }

  if (esd) { 
    AliFMDDebug(2, ("Writing FMD data to ESD tree"));
    esd->SetFMDData(fESDObj);
  }

  if (!fDiagnostics || !esd) return;
  static bool first = true;
  // This is most likely NOT the event number you'd like to use. It
  // has nothing to do with the 'real' event number. 
  // - That's OK.  We just use it for the name of the directory -
  // nothing else.  Christian
  Int_t evno = esd->GetEventNumberInFile(); 
  AliFMDDebug(1, ("Writing diagnostics histograms to FMD.Diag.root/%03d",evno));
  TFile f("FMD.Diag.root", (first ? "RECREATE" : "UPDATE"));
  first = false;
  f.cd(); 
  TDirectory* d = f.mkdir(Form("%03d", evno), 
			  Form("Diagnostics histograms for event # %d", evno));
  d->cd();
  if (fDiagStep1) fDiagStep1->Write();
  if (fDiagStep2) fDiagStep2->Write();
  if (fDiagStep3) fDiagStep3->Write();
  if (fDiagStep4) fDiagStep4->Write();
  if (fDiagAll)   fDiagAll->Write();
  d->Write();
  f.Write();
  f.Close();

  if (fDiagStep1) fDiagStep1->Reset();
  if (fDiagStep2) fDiagStep2->Reset();
  if (fDiagStep3) fDiagStep3->Reset();
  if (fDiagStep4) fDiagStep4->Reset();
  if (fDiagAll)   fDiagAll->Reset();
}

//____________________________________________________________________
void
AliFMDReconstructor::FillESD(AliRawReader*, TTree* clusterTree, 
			     AliESDEvent* esd) const
{
  TTree* dummy = 0;
  FillESD(dummy, clusterTree, esd);
}

//____________________________________________________________________
//
// EOF
//
