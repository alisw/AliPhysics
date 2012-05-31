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
/** @file    AliFMD.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 17:59:18 2006
    @brief   Implementation of AliFMD base class 
*/
//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. This class
// is the driver for especially simulation. 
//
// The Forward Multiplicity Detector consists of 3 sub-detectors FMD1,
// FMD2, and FMD3, each of which has 1 or 2 rings of silicon sensors. 
//                                                       
// This is the base class for all FMD manager classes. 
//                    
// The actual code is done by various separate classes.   Below is
// diagram showing the relationship between the various FMD classes
// that handles the simulation
//
//
//       +----------+   +----------+   
//       | AliFMDv1 |	| AliFMDv0 |   
//       +----------+   +----------+   
//            |              |                    +-----------------+
//       +----+--------------+                 +--| AliFMDDigitizer |
//       |                                     |  +-----------------+
//       |           +---------------------+   |
//       |        +--| AliFMDBaseDigitizer |<--+
//       V     1  |  +---------------------+   |
//  +--------+<>--+                            |  +------------------+
//  | AliFMD |                                 +--| AliFMDSDigitizer |    
//  +--------+<>--+                               +------------------+       
//	       1  |  +---------------------+
//	          +--| AliFMDReconstructor |
//	  	     +---------------------+
//
// *  AliFMD 
//    This defines the interface for the various parts of AliROOT that
//    uses the FMD, like AliFMDSimulator, AliFMDDigitizer, 
//    AliFMDReconstructor, and so on. 
//
// *  AliFMDv0
//    This is a concrete implementation of the AliFMD interface. 
//    It is the responsibility of this class to create the FMD
//    geometry.
//
// *  AliFMDv1 
//    This is a concrete implementation of the AliFMD interface. 
//    It is the responsibility of this class to create the FMD
//    geometry, process hits in the FMD, and serve hits and digits to
//    the various clients. 
//  
// *  AliFMDSimulator
//    This is the base class for the FMD simulation tasks.   The
//    simulator tasks are responsible to implment the geoemtry, and
//    process hits. 
//                                                                          
// *  AliFMDReconstructor
//    This is a concrete implementation of the AliReconstructor that
//    reconstructs pseudo-inclusive-multiplicities from digits (raw or
//    from simulation)
//
// Calibration and geometry parameters are managed by separate
// singleton managers.  These are AliFMDGeometry and
// AliFMDParameters.  Please refer to these classes for more
// information on these.
//

// These files are not in the same directory, so there's no reason to
// ask the preprocessor to search in the current directory for these
// files by including them with `#include "..."' 
#include <TBrowser.h>		// ROOT_TBrowser
#include <TClonesArray.h>	// ROOT_TClonesArray
#include <TGeoGlobalMagField.h> // ROOT_TGeoGlobalMagField
#include <TGeoManager.h>        // ROOT_TGeoManager
#include <TRotMatrix.h>		// ROOT_TRotMatrix
#include <TTree.h>		// ROOT_TTree
#include <TVector2.h>           // ROOT_TVector2 
#include <TVirtualMC.h>	        // ROOT_TVirtualMC
#include <cmath>                // __CMATH__

#include <AliDigitizationInput.h>	// ALIRUNDIGITIZER_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRun.h>		// ALIRUN_H
#include <AliMC.h>		// ALIMC_H
#include <AliMagF.h>		// ALIMAGF_H
// #include <AliLog.h>		// ALILOG_H
#include "AliFMDDebug.h" // Better debug macros
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDSDigit.h"	// ALIFMDSDIGIT_H
#include "AliFMDHit.h"		// ALIFMDHIT_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMDDigitizer.h"	// ALIFMDDIGITIZER_H
#include "AliFMDHitDigitizer.h"	// ALIFMDSDIGITIZER_H
// #define USE_SSDIGITIZER 
//#ifdef USE_SSDIGITIZER
//# include "AliFMDSSDigitizer.h"	// ALIFMDSDIGITIZER_H
//#endif
// #include "AliFMDGeometryBuilder.h"
#include "AliFMDRawWriter.h"	// ALIFMDRAWWRITER_H
#include "AliFMDRawReader.h"	// ALIFMDRAWREADER_H
#include "AliTrackReference.h" 
#include "AliFMDStripIndex.h"
#include "AliFMDParameters.h"
#include "AliFMDReconstructor.h"

//____________________________________________________________________
ClassImp(AliFMD)
#if 0
  ; // This is to keep Emacs from indenting the next line 
#endif 

//____________________________________________________________________
AliFMD::AliFMD()
  : AliDetector(),
    fSDigits(0), 
    fNsdigits(0),
    fDetailed(kTRUE),
    fUseOld(kFALSE),
    fUseAssembly(kTRUE),
    fBad(0) 
{
  //
  // Default constructor for class AliFMD
  //
  AliFMDDebug(10, ("\tDefault CTOR"));
  fHits        = 0;
  fDigits      = 0;
  fIshunt      = 0;
  // fBad         = new TClonesArray("AliFMDHit");
}

//____________________________________________________________________
AliFMD::AliFMD(const char *name, const char *title)
  : AliDetector (name, title),
    fSDigits(0),
    fNsdigits(0),
    fDetailed(kTRUE),
    fUseOld(kFALSE),
    fUseAssembly(kFALSE),
    fBad(0)
{
  //
  // Standard constructor for Forward Multiplicity Detector
  //
  AliFMDDebug(10, ("\tStandard CTOR"));
  // fBad         = new TClonesArray("AliFMDHit");
  
  // Initialise Hit array
  // HitsArray();
  // gAlice->GetMCApp()->AddHitList(fHits);

  // (S)Digits for the detectors disk
  // DigitsArray();
  // SDigitsArray();
  
  // CHC: What is this?
  fIshunt = 0;
  //PH  SetMarkerColor(kRed);
  //PH  SetLineColor(kYellow);
}

//____________________________________________________________________
AliFMD::~AliFMD ()
{
  // Destructor for base class AliFMD
  if (fHits) {
    fHits->Delete();
    delete fHits;
    fHits = 0;
  }
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits = 0;
  }
  if (fSDigits) {
    fSDigits->Delete();
    delete fSDigits;
    fSDigits = 0;
  }
  if (fBad) {
    fBad->Delete();
    delete fBad;
    fBad = 0;
  }
}


//====================================================================
//
// GEometry ANd Traking
//
//____________________________________________________________________
void 
AliFMD::CreateGeometry()
{
  //
  // Create the geometry of Forward Multiplicity Detector.  The actual
  // construction of the geometry is delegated to the class
  // AliFMDGeometryBuilder, invoked by the singleton manager
  // AliFMDGeometry. 
  //
  AliFMDGeometry*  fmd = AliFMDGeometry::Instance();
  fmd->SetDetailed(fDetailed);
  fmd->UseAssembly(fUseAssembly);
  fmd->Build();
}    

//____________________________________________________________________
void AliFMD::CreateMaterials() 
{
  // Define the materials and tracking mediums needed by the FMD
  // simulation.   These mediums are made by sending the messages
  // AliMaterial, AliMixture, and AliMedium to the passed AliModule
  // object module.   The defined mediums are 
  // 
  //	FMD Si$		Silicon (active medium in sensors)
  //	FMD C$		Carbon fibre (support cone for FMD3 and vacuum pipe)
  //	FMD Al$		Aluminium (honeycomb support plates)
  //	FMD PCB$	Printed Circuit Board (FEE board with VA1_3)
  //	FMD Chip$	Electronics chips (currently not used)
  //	FMD Air$	Air (Air in the FMD)
  //	FMD Plastic$	Plastic (Support legs for the hybrid cards)
  //
  // The geometry builder should really be the one that creates the
  // materials, but the architecture of AliROOT makes that design
  // akward.  What should happen, was that the AliFMDGeometryBuilder
  // made the mediums, and that this class retrives pointers from the
  // TGeoManager, and registers the mediums here.  Alas, it's not
  // really that easy. 
  //
  AliFMDDebug(10, ("\tCreating materials"));
  // Get pointer to geometry singleton object. 
  AliFMDGeometry* geometry = AliFMDGeometry::Instance();
  geometry->Init();
#if 0
  if (gGeoManager && gGeoManager->GetMedium("FMD Si$")) {
    // We need to figure out the some stuff about the geometry
    fmd->ExtractGeomInfo();
    return;
  }
#endif  
  Int_t    id;
  Double_t a                = 0;
  Double_t z                = 0;
  Double_t density          = 0;
  Double_t radiationLength  = 0;
  Double_t absorbtionLength = 999;
  Int_t    fieldType        = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();     // Field type 
  Double_t maxField         = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();     // Field max.
  Double_t maxBending       = 0;     // Max Angle
  Double_t maxStepSize      = 0.001; // Max step size 
  Double_t maxEnergyLoss    = 1;     // Max Delta E
  Double_t precision        = 0.001; // Precision
  Double_t minStepSize      = 0.001; // Minimum step size 
 
  // Silicon 
  a                = 28.0855;
  z                = 14.;
  density          = geometry->GetSiDensity();
  radiationLength  = 9.36;
  maxBending       = 1;
  maxStepSize      = .001;
  precision        = .001;
  minStepSize      = .001;
  id               = kSiId;
  AliMaterial(id, "Si$", a, z, density, radiationLength, absorbtionLength);
  AliMedium(kSiId, "Si$", id,1,fieldType,maxField,maxBending,
	    maxStepSize,maxEnergyLoss,precision,minStepSize);
  

  // Carbon 
  a                = 12.011;
  z                = 6.;
  density          = 2.265;
  radiationLength  = 18.8;
  maxBending       = 10;
  maxStepSize      = .01;
  precision        = .003;
  minStepSize      = .003;
  id               = kCarbonId;
  AliMaterial(id, "Carbon$", a, z, density, radiationLength, absorbtionLength);
  AliMedium(kCarbonId, "Carbon$", id,0,fieldType,maxField,maxBending,
		    maxStepSize,maxEnergyLoss,precision,minStepSize);

  // Aluminum
  a                = 26.981539;
  z                = 13.;
  density          = 2.7;
  radiationLength  = 8.9;
  id               = kAlId;
  AliMaterial(id, "Aluminum$",a,z, density, radiationLength, absorbtionLength);
  AliMedium(kAlId, "Aluminum$", id, 0, fieldType, maxField, maxBending,
	    maxStepSize, maxEnergyLoss, precision, minStepSize);
  
  
  // Copper 
  a                = 63.546;
  z                = 29;
  density          =  8.96;
  radiationLength  =  1.43;
  id               = kCopperId;
  AliMaterial(id, "Copper$", 
		      a, z, density, radiationLength, absorbtionLength);
  AliMedium(kCopperId, "Copper$", id, 0, fieldType, maxField, maxBending,
	    maxStepSize, maxEnergyLoss, precision, minStepSize);
  

  // Silicon chip 
  {
    Float_t as[] = { 12.0107,      14.0067,      15.9994,
		      1.00794,     28.0855,     107.8682 };
    Float_t zs[] = {  6.,           7.,           8.,
		      1.,          14.,          47. };
    Float_t ws[] = {  0.039730642,  0.001396798,  0.01169634,
		      0.004367771,  0.844665,     0.09814344903 };
    density          = 2.36436;
    maxBending       = 10;
    maxStepSize      = .01;
    precision        = .003;
    minStepSize      = .003;
    id               = kSiChipId;
    AliMixture(id, "Si Chip$", as, zs, density, 6, ws);
    AliMedium(kSiChipId, "Si Chip$",  id, 0, fieldType, maxField, maxBending, 
	      maxStepSize, maxEnergyLoss, precision, minStepSize);
  }
  
  // Kaption
  {
    Float_t as[] = { 1.00794,  12.0107,  14.010,   15.9994};
    Float_t zs[] = { 1.,        6.,       7.,       8.};
    Float_t ws[] = { 0.026362,  0.69113,  0.07327,  0.209235};
    density          = 1.42;
    maxBending       = 1;
    maxStepSize      = .001;
    precision        = .001;
    minStepSize      = .001;
    id               = kKaptonId;
    AliMixture(id, "Kaption$", as, zs, density, 4, ws);
    AliMedium(kKaptonId, "Kaption$", id,0,fieldType,maxField,maxBending,
	      maxStepSize,maxEnergyLoss,precision,minStepSize);
  }

  // Air
  {
    Float_t as[] = { 12.0107, 14.0067,   15.9994,  39.948 };
    Float_t zs[] = {  6.,      7.,       8.,       18. };
    Float_t ws[] = { 0.000124, 0.755267, 0.231781, 0.012827 }; 
    density      = .00120479;
    maxBending   = 1;
    maxStepSize  = .001;
    precision    = .001;
    minStepSize  = .001;
    id           = kAirId;
    AliMixture(id, "Air$", as, zs, density, 4, ws);
    AliMedium(kAirId, "Air$", id,0,fieldType,maxField,maxBending,
	      maxStepSize,maxEnergyLoss,precision,minStepSize);
  }
  
  // PCB
  {
    Float_t zs[] = { 14.,         20.,         13.,         12.,
		      5.,         22.,         11.,         19.,
		     26.,          9.,          8.,          6.,
		      7.,          1.};
    Float_t as[] = { 28.0855,     40.078,      26.981538,   24.305, 
		     10.811,      47.867,      22.98977,    39.0983,
		     55.845,      18.9984,     15.9994,     12.0107,
		     14.0067,      1.00794};
    Float_t ws[] = {  0.15144894,  0.08147477,  0.04128158,  0.00904554, 
		      0.01397570,  0.00287685,  0.00445114,  0.00498089,
		      0.00209828,  0.00420000,  0.36043788,  0.27529426,
		      0.01415852,  0.03427566};
    density      = 1.8;
    maxBending   = 1;
    maxStepSize  = .001;
    precision    = .001;
    minStepSize  = .001;
    id           = kPcbId;
    AliMixture(id, "PCB$", as, zs, density, 14, ws);
    AliMedium(kPcbId, "PCB$", id,0,fieldType,maxField,maxBending,
	      maxStepSize,maxEnergyLoss,precision,minStepSize);
  }
  
  // Stainless steel
  {
    Float_t as[] = { 55.847, 51.9961, 58.6934, 28.0855 };
    Float_t zs[] = { 26.,    24.,     28.,     14.     };
    Float_t ws[] = { .715,   .18,     .1,      .005    };
    density      = 7.88;
    id           = kSteelId;
    AliMixture(id, "Steel$", as, zs, density, 4, ws);
    AliMedium(kSteelId, "Steel$", id, 0, fieldType, maxField, maxBending, 
	      maxStepSize, maxEnergyLoss, precision, minStepSize);
  }
  // Plastic 
  {
    Float_t as[] = { 1.01, 12.01 };
    Float_t zs[] = { 1.,   6.    };
    Float_t ws[] = { 1.,   1.    };
    density      = 1.03;
    maxBending   = 10;
    maxStepSize  = .01;
    precision    = .003;
    minStepSize  = .003;
    id           = kPlasticId;
    AliMixture(id, "Plastic$", as, zs, density, -2, ws);
    AliMedium(kPlasticId, "Plastic$", id,0,fieldType,maxField,maxBending,
	      maxStepSize,maxEnergyLoss,precision,minStepSize);
  }

}

#if 0
//____________________________________________________________________
void  
AliFMD::SetTrackingParameters(Int_t imed, 
			      Float_t gamma,                 
			      Float_t electron, 
			      Float_t neutral_hadron, 
			      Float_t charged_hadron, 
			      Float_t muon,
			      Float_t electron_bremstrahlung, 
			      Float_t muon__bremstrahlung, 
			      Float_t electron_delta,
			      Float_t muon_delta,
			      Float_t muon_pair,
			      Int_t   annihilation, 
			      Int_t   bremstrahlung, 
			      Int_t   compton_scattering, 
			      Int_t   decay,
			      Int_t   delta_ray, 
			      Int_t   hadronic, 
			      Int_t   energy_loss, 
			      Int_t   multiple_scattering, 
			      Int_t   pair_production, 
			      Int_t   photon_production, 
			      Int_t   rayleigh_scattering)
{
  // Disabled by request of FCA, kept for reference only
  if (!gMC) return;
  TArrayI& idtmed = *(GetIdtmed());
  Int_t    iimed  = idtmed[imed];
  // gMC->Gstpar(iimed, "CUTGAM",	gamma);
  // gMC->Gstpar(iimed, "CUTELE",	electron);
  // gMC->Gstpar(iimed, "CUTNEU",	neutral_hadron);
  // gMC->Gstpar(iimed, "CUTHAD",	charged_hadron);
  // gMC->Gstpar(iimed, "CUTMUO",	muon);
  // gMC->Gstpar(iimed, "BCUTE",	electron_bremstrahlung);
  // gMC->Gstpar(iimed, "BCUTM",	muon__bremstrahlung);
  // gMC->Gstpar(iimed, "DCUTE",	electron_delta);
  // gMC->Gstpar(iimed, "DCUTM",	muon_delta);
  // gMC->Gstpar(iimed, "PPCUTM",	muon_pair);
  // gMC->Gstpar(iimed, "ANNI",	Float_t(annihilation));
  // gMC->Gstpar(iimed, "BREM",	Float_t(bremstrahlung));
  // gMC->Gstpar(iimed, "COMP",	Float_t(compton_scattering));
  // gMC->Gstpar(iimed, "DCAY",	Float_t(decay));
  // gMC->Gstpar(iimed, "DRAY",	Float_t(delta_ray));
  // gMC->Gstpar(iimed, "HADR",	Float_t(hadronic));
  // gMC->Gstpar(iimed, "LOSS",	Float_t(energy_loss));
  // gMC->Gstpar(iimed, "MULS",	Float_t(multiple_scattering));
  // gMC->Gstpar(iimed, "PAIR",	Float_t(pair_production));
  // gMC->Gstpar(iimed, "PHOT",	Float_t(photon_production));
  // gMC->Gstpar(iimed, "RAYL",	Float_t(rayleigh_scattering));
}
#endif

//____________________________________________________________________
void  
AliFMD::Init()
{
  // Initialize the detector 
  // 
  AliFMDDebug(1, ("Initialising FMD detector object"));
  TVirtualMC*      mc     = TVirtualMC::GetMC();
  AliFMDGeometry*  fmd    = AliFMDGeometry::Instance();
  TArrayI          actGeo = fmd->ActiveIds();
  bool             valid  = true;
  if (actGeo.fN <= 0) valid = false;
  else { 
    for (int i = 0; i < actGeo.fN; i++) {
      if (actGeo[i] < 0) { 
	valid = false;
	break;
      }
    }
  }
  if (!valid) { 
    AliFMDDebug(1, ("Extracting geometry info from loaded geometry"));
    fmd->ExtractGeomInfo();
    actGeo = fmd->ActiveIds();
  }
  TArrayI          actVmc(actGeo.fN);
  for (Int_t i = 0; i < actGeo.fN; i++) {
    if (actGeo[i] < 0) { 
      AliError(Form("Invalid id: %d", actGeo[i]));
      continue;
    }
    TGeoVolume *sens = gGeoManager->GetVolume(actGeo[i]);
    if (!sens) {
      AliError(Form("No TGeo volume for sensitive volume ID=%d",actGeo[i]));
      continue;
    }   
    actVmc[i] = mc->VolId(sens->GetName());
    AliFMDDebug(1, ("Active vol id # %d: %d changed to %d", 
		    i, actGeo[i], actVmc[i]));
  }
  fmd->SetActive(actVmc.fArray, actVmc.fN);
  // fmd->InitTransformations();
}

//____________________________________________________________________
void
AliFMD::FinishEvent()
{
  // Called at the end of the an event in simulations.  If the debug
  // level is high enough, then the `bad' hits are printed.
  // 
  if (AliLog::GetDebugLevel("FMD", "AliFMD") < 10) return;
  if (fBad && fBad->GetEntries() > 0) {
    AliWarning(Form("got %d 'bad' hits", fBad->GetEntries()));
    TIter next(fBad);
    AliFMDHit* hit;
    while ((hit = static_cast<AliFMDHit*>(next()))) hit->Print("D");
    fBad->Clear();
  }
}



//====================================================================
//
// Hit and Digit managment 
//
//____________________________________________________________________
void 
AliFMD::MakeBranch(Option_t * option)
{
  // Create Tree branches for the FMD.
  //
  // Options:
  //
  //    H          Make a branch of TClonesArray of AliFMDHit's
  //    D          Make a branch of TClonesArray of AliFMDDigit's
  //    S          Make a branch of TClonesArray of AliFMDSDigit's
  // 
  const Int_t kBufferSize = 16000;
  TString branchname(GetName());
  TString opt(option);
  
  if (opt.Contains("H", TString::kIgnoreCase)) {
    HitsArray();
    AliDetector::MakeBranch(option); 
  }
  if (opt.Contains("D", TString::kIgnoreCase)) { 
    DigitsArray();
    MakeBranchInTree(fLoader->TreeD(), branchname.Data(),
		     &fDigits, kBufferSize, 0);
  }
  if (opt.Contains("S", TString::kIgnoreCase)) { 
    SDigitsArray();
    MakeBranchInTree(fLoader->TreeS(), branchname.Data(),
		     &fSDigits, kBufferSize, 0);
  }
}

//____________________________________________________________________
void 
AliFMD::SetTreeAddress()
{
  // Set branch address for the Hits, Digits, and SDigits Tree.
  if (fLoader->TreeH()) HitsArray();
  AliDetector::SetTreeAddress();

  TTree *treeD = fLoader->TreeD();
  if (treeD) {
    DigitsArray();
    TBranch* branch = treeD->GetBranch ("FMD");
    if (branch) branch->SetAddress(&fDigits);
  }

  TTree *treeS = fLoader->TreeS();
  if (treeS) {
    SDigitsArray();
    TBranch* branch = treeS->GetBranch ("FMD");
    if (branch) branch->SetAddress(&fSDigits);
  }
}

//____________________________________________________________________
void 
AliFMD::SetHitsAddressBranch(TBranch *b)
{
  // Set the TClonesArray to read hits into. 
  b->SetAddress(&fHits);
}
//____________________________________________________________________
void 
AliFMD::SetSDigitsAddressBranch(TBranch *b)
{
  // Set the TClonesArray to read hits into. 
  b->SetAddress(&fSDigits);
}

//____________________________________________________________________
void 
AliFMD::AddHit(Int_t track, Int_t *vol, Float_t *hits) 
{
  // Add a hit to the hits tree 
  // 
  // The information of the two arrays are decoded as 
  // 
  // Parameters
  //    track	 	     Track #
  //    ivol[0]  [UShort_t ] Detector # 
  //    ivol[1]	 [Char_t   ] Ring ID 
  //    ivol[2]	 [UShort_t ] Sector #
  //    ivol[3]	 [UShort_t ] Strip # 
  //    hits[0]	 [Float_t  ] Track's X-coordinate at hit 
  //    hits[1]	 [Float_t  ] Track's Y-coordinate at hit
  //    hits[3]  [Float_t  ] Track's Z-coordinate at hit
  //    hits[4]  [Float_t  ] X-component of track's momentum 	       	 
  //    hits[5]	 [Float_t  ] Y-component of track's momentum	       	 
  //    hits[6]	 [Float_t  ] Z-component of track's momentum	       	
  //    hits[7]	 [Float_t  ] Energy deposited by track		       	
  //    hits[8]	 [Int_t    ] Track's particle Id # 
  //    hits[9]	 [Float_t  ] Time when the track hit
  // 
  // 
  AddHitByFields(track, 
		 UShort_t(vol[0]),  // Detector # 
		 Char_t(vol[1]),    // Ring ID
		 UShort_t(vol[2]),  // Sector # 
		 UShort_t(vol[3]),  // Strip # 
		 hits[0],           // X
		 hits[1],           // Y
		 hits[2],           // Z
		 hits[3],           // Px
		 hits[4],           // Py
		 hits[5],           // Pz
		 hits[6],           // Energy loss 
		 Int_t(hits[7]),    // PDG 
		 hits[8]);          // Time
}

//____________________________________________________________________
AliFMDHit*
AliFMD::AddHitByFields(Int_t    track, 
		       UShort_t detector, 
		       Char_t   ring, 
		       UShort_t sector, 
		       UShort_t strip, 
		       Float_t  x, 
		       Float_t  y, 
		       Float_t  z,
		       Float_t  px, 
		       Float_t  py, 
		       Float_t  pz,
		       Float_t  edep,
		       Int_t    pdg,
		       Float_t  t, 
		       Float_t  l, 
		       Bool_t   stop)
{
  // Add a hit to the list
  //
  // Parameters:
  // 
  //    track	  Track #
  //    detector  Detector # (1, 2, or 3)                      
  //    ring	  Ring ID ('I' or 'O')
  //    sector	  Sector # (For inner/outer rings: 0-19/0-39)
  //    strip	  Strip # (For inner/outer rings: 0-511/0-255)
  //    x	  Track's X-coordinate at hit
  //    y	  Track's Y-coordinate at hit
  //    z	  Track's Z-coordinate at hit
  //    px	  X-component of track's momentum 
  //    py	  Y-component of track's momentum
  //    pz	  Z-component of track's momentum
  //    edep	  Energy deposited by track
  //    pdg	  Track's particle Id #
  //    t	  Time when the track hit 
  //    l         Track length through the material. 
  //    stop      Whether track was stopped or disappeared
  // 
  TClonesArray& a = *(HitsArray());
  // Search through the list of already registered hits, and see if we
  // find a hit with the same parameters.  If we do, then don't create
  // a new hit, but rather update the energy deposited in the hit.
  // This is done, so that a FLUKA based simulation will get the
  // number of hits right, not just the enerrgy deposition. 
  AliFMDHit* hit = 0;
  for (Int_t i = 0; i < fNhits; i++) {
    if (!a.At(i)) continue;
    hit = static_cast<AliFMDHit*>(a.At(i));
    if (hit->Detector() == detector 
	&& hit->Ring() == ring
	&& hit->Sector() == sector 
	&& hit->Strip() == strip
	&& hit->Track() == track) {
      AliFMDDebug(1, ("already had a hit in FMD%d%c[%2d,%3d] for track # %d,"
		       " adding energy (%f) to that hit (%f) -> %f", 
		       detector, ring, sector, strip, track, edep, hit->Edep(),
		       hit->Edep() + edep));
      hit->SetEdep(hit->Edep() + edep);
      return hit;
    }
  }
  // If hit wasn't already registered, do so know. 
  hit = new (a[fNhits]) AliFMDHit(fIshunt, track, detector, ring, sector, 
				  strip, x, y, z, px, py, pz, edep, pdg, t, 
				  l, stop);
  // gMC->AddTrackReference(track, 12);
  fNhits++;
  
  //Reference track

  AliMC *mcApplication = (AliMC*)gAlice->GetMCApp();
  
  AliTrackReference* trackRef = 
    AddTrackReference(mcApplication->GetCurrentTrackNumber(), 
		      AliTrackReference::kFMD); 
  UInt_t stripId = AliFMDStripIndex::Pack(detector,ring,sector,strip);
  trackRef->SetUserId(stripId);
  
  
  
  return hit;
}

//____________________________________________________________________
void 
AliFMD::AddDigit(Int_t* digits, Int_t*)
{
  // Add a digit to the Digit tree 
  // 
  // Paramters 
  //
  //    digits[0]  [UShort_t] Detector #
  //    digits[1]  [Char_t]   Ring ID
  //    digits[2]  [UShort_t] Sector #
  //    digits[3]  [UShort_t] Strip #
  //    digits[4]  [UShort_t] ADC Count 
  //    digits[5]  [Short_t]  ADC Count, -1 if not used
  //    digits[6]  [Short_t]  ADC Count, -1 if not used 
  // 
  AddDigitByFields(UShort_t(digits[0]),  // Detector #
		   Char_t(digits[1]),    // Ring ID
		   UShort_t(digits[2]),  // Sector #
		   UShort_t(digits[3]),  // Strip #
		   UShort_t(digits[4]),  // ADC Count1 
		   Short_t(digits[5]), 	 // ADC Count2 
		   Short_t(digits[6]),   // ADC Count3 
		   Short_t(digits[7])); 
}

//____________________________________________________________________
void 
AliFMD::AddDigitByFields(UShort_t       detector, 
			 Char_t         ring, 
			 UShort_t       sector, 
			 UShort_t       strip, 
			 UShort_t       count1, 
			 Short_t        count2,
			 Short_t        count3, 
			 Short_t        count4,
			 UShort_t	nrefs,
			 Int_t*		refs)
{
  // add a real digit - as coming from data
  // 
  // Parameters 
  //
  //    detector  Detector # (1, 2, or 3)                      
  //    ring	  Ring ID ('I' or 'O')
  //    sector	  Sector # (For inner/outer rings: 0-19/0-39)
  //    strip	  Strip # (For inner/outer rings: 0-511/0-255)
  //    count1    ADC count (a 10-bit word)
  //    count2    ADC count (a 10-bit word), or -1 if not used
  //    count3    ADC count (a 10-bit word), or -1 if not used
  TClonesArray& a = *(DigitsArray());
  
  AliFMDDebug(15, ("Adding digit # %5d/%5d for FMD%d%c[%2d,%3d]"
		   "=(%d,%d,%d,%d) with %d tracks",
		   fNdigits-1, a.GetEntriesFast(),
		   detector, ring, sector, strip, 
		   count1, count2, count3, count4, nrefs));
  new (a[fNdigits++]) 
    AliFMDDigit(detector, ring, sector, strip, 
		count1, count2, count3, count4, nrefs, refs);
  
}

//____________________________________________________________________
void 
AliFMD::AddSDigit(Int_t* digits)
{
  // Add a digit to the SDigit tree 
  // 
  // Paramters 
  //
  //    digits[0]  [UShort_t] Detector #
  //    digits[1]  [Char_t]   Ring ID
  //    digits[2]  [UShort_t] Sector #
  //    digits[3]  [UShort_t] Strip #
  //    digits[4]  [Float_t]  Total energy deposited 
  //    digits[5]  [UShort_t] ADC Count 
  //    digits[6]  [Short_t]  ADC Count, -1 if not used
  //    digits[7]  [Short_t]  ADC Count, -1 if not used 
  // 
  AddSDigitByFields(UShort_t(digits[0]),   // Detector #
		    Char_t(digits[1]),     // Ring ID
		    UShort_t(digits[2]),   // Sector #
		    UShort_t(digits[3]),   // Strip #
		    Float_t(digits[4]),    // Edep
		    UShort_t(digits[5]),   // ADC Count1 
		    Short_t(digits[6]),    // ADC Count2 
		    Short_t(digits[7]),    // ADC Count3 
		    Short_t(digits[8]),    // ADC Count4
		    UShort_t(digits[9]),   // N particles
		    UShort_t(digits[10])); // N primaries
}

//____________________________________________________________________
void 
AliFMD::AddSDigitByFields(UShort_t       detector, 
			  Char_t         ring, 
			  UShort_t       sector, 
			  UShort_t       strip, 
			  Float_t        edep,
			  UShort_t       count1, 
			  Short_t        count2,
			  Short_t        count3, 
			  Short_t        count4, 
			  UShort_t       ntot, 
			  UShort_t       nprim,
			  Int_t*	 refs)
{
  // add a summable digit
  // 
  // Parameters 
  //
  //    detector  Detector # (1, 2, or 3)                      
  //    ring	  Ring ID ('I' or 'O')
  //    sector	  Sector # (For inner/outer rings: 0-19/0-39)
  //    strip	  Strip # (For inner/outer rings: 0-511/0-255)
  //    edep      Total energy deposited
  //    count1    ADC count (a 10-bit word)
  //    count2    ADC count (a 10-bit word), or -1 if not used
  //    count3    ADC count (a 10-bit word), or -1 if not used
  //
  TClonesArray& a = *(SDigitsArray());
  // AliFMDDebug(0, ("Adding sdigit # %d", fNsdigits));
  
  AliFMDDebug(15, ("Adding sdigit # %5d/%5d for FMD%d%c[%2d,%3d]"
		   "=(%d,%d,%d,%d) with %d tracks %d primaries (%p)",
		   fNsdigits-1, a.GetEntriesFast(),
		   detector, ring, sector, strip, 
		   count1, count2, count3, count4, ntot, nprim, refs));
  new (a[fNsdigits++]) 
    AliFMDSDigit(detector, ring, sector, strip, edep, 
		 count1, count2, count3, count4, ntot, nprim, refs);
}

//____________________________________________________________________
void 
AliFMD::ResetSDigits()
{
  // Reset number of digits and the digits array for this detector. 
  //
  fNsdigits   = 0;
  if (fSDigits) fSDigits->Clear();
}


//____________________________________________________________________
TClonesArray*
AliFMD::HitsArray() 
{
  // Initialize hit array if not already, and return pointer to it. 
  if (!fHits) { 
    fHits = new TClonesArray("AliFMDHit", 1000);
    fNhits = 0;
    if (gAlice && gAlice->GetMCApp() && gAlice->GetMCApp()->GetHitLists()) 
      gAlice->GetMCApp()->AddHitList(fHits);
  }
  return fHits;
}

//____________________________________________________________________
TClonesArray*
AliFMD::DigitsArray() 
{
  // Initialize digit array if not already, and return pointer to it. 
  if (!fDigits) { 
    fDigits = new TClonesArray("AliFMDDigit", 1000);
    fNdigits = 0;
  }
  return fDigits;
}

//____________________________________________________________________
TClonesArray*
AliFMD::SDigitsArray() 
{
  // Initialize digit array if not already, and return pointer to it. 
  if (!fSDigits) { 
    fSDigits = new TClonesArray("AliFMDSDigit", 1000);
    fNsdigits = 0;
  }
  return fSDigits;
}

//====================================================================
//
// Digitization 
//
//____________________________________________________________________
void 
AliFMD::Hits2Digits() 
{
  // Create AliFMDDigit's from AliFMDHit's.  This is done by making a
  // AliFMDDigitizer, and executing that code.
  // 
  AliFMDHitDigitizer digitizer(this, AliFMDHitDigitizer::kDigits);
  digitizer.Init();
  digitizer.Digitize("");
}

//____________________________________________________________________
void 
AliFMD::Hits2SDigits() 
{
  // Create AliFMDSDigit's from AliFMDHit's.  This is done by creating
  // an AliFMDSDigitizer object, and executing it. 
  // 
  AliFMDHitDigitizer digitizer(this, AliFMDHitDigitizer::kSDigits);
  digitizer.Init();
  digitizer.Digitize("");
}

  
//____________________________________________________________________
AliDigitizer* 
AliFMD::CreateDigitizer(AliDigitizationInput* digInput) const
{
  // Create a digitizer object 
  
  /* This is what we probably _should_ do */
  AliFMDBaseDigitizer* digitizer = 0;
  
#ifdef USE_SSDIGITIZER
  digitizer = new AliFMDSSDigitizer(digInput);
#else 
  /* This is what we actually do, and will work */
#if 0
  AliInfo("SDigit->Digit conversion not really supported, "
	  "doing Hit->Digit conversion instead");
#endif
  digitizer = new AliFMDDigitizer(digInput);
#endif
  return digitizer;
}

//====================================================================
//
// Raw data simulation 
//
//__________________________________________________________________
void 
AliFMD::Digits2Raw() 
{
  // Turn digits into raw data. 
  // 
  // This uses the class AliFMDRawWriter to do the job.   Please refer
  // to that class for more information. 
  AliFMDRawWriter writer(this);
  writer.Exec();
}

//====================================================================
//
// Raw data reading 
//
//__________________________________________________________________
Bool_t
AliFMD::Raw2SDigits(AliRawReader* reader) 
{
  // Turn digits into raw data. 
  // 
  // This uses the class AliFMDRawWriter to do the job.   Please refer
  // to that class for more information. 
  AliFMDParameters::Instance()->Init();
  MakeTree("S");
  MakeBranch("S");
  
  TClonesArray*       sdigits = SDigits();
  AliFMDReconstructor rec;
  
  // The two boolean arguments
  //   Make sdigits instead of digits 
  //   Subtract the pedestal off the signal
  rec.Digitize(reader, sdigits);
  // 
  // Bool_t ret = fmdReader.ReadAdcs(sdigits, kTRUE, kTRUE);
  // sdigits->ls();
  UShort_t ns = sdigits->GetEntriesFast();
  if (AliLog::GetDebugLevel("FMD", 0) > 5) {
    for (UShort_t i = 0; i < ns; i++) 
      sdigits->At(i)->Print("pl");
  } 
  AliFMDDebug(1, ("Got a total of %d SDigits", ns));

  fLoader->TreeS()->Fill();
  ResetSDigits();
  fLoader->WriteSDigits("OVERWRITE");

  return kTRUE;
}


//====================================================================
//
// Utility 
//
//__________________________________________________________________
void 
AliFMD::Browse(TBrowser* b) 
{
  // Browse this object. 
  //
  AliFMDDebug(30, ("\tBrowsing the FMD"));
  AliDetector::Browse(b);
  b->Add(AliFMDGeometry::Instance());
}

//____________________________________________________________________	
void
AliFMD::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  // 
  // This code was made by Raffaele Grosso <rgrosso@mail.cern.ch>.  I
  // (cholm) will probably want to change it.   For one, I think it
  // should be the job of the geometry manager to deal with this. 
  AliInfo("Add FMD alignable volumes");
  AliFMDGeometry::Instance()->SetAlignableVolumes();
#if 0  
  for(size_t f = 1; f <= 3; f++){ // Detector 1,2,3
    for(size_t tb =  0; tb <2 ; tb++){ // Top/Bottom 
      char     stb = tb == 0 ? 'T' : 'B';
      unsigned min = tb == 0 ? 0   : 5;

      TString halfVol(Form("/ALIC_1/F%dM%c_%d", f, stb, f));
      TString halfSym(halfVol);
      if(!gGeoManager->SetAlignableEntry(halfSym.Data(),halfVol.Data()))
	AliFatal(Form("Alignable entry %s not created. "
		      "Volume path %s not valid", 
		      halfSym.Data(),halfVol.Data()));
      for(size_t io = 0; io < 2; io++){ // inner, outer 
	if (f==1 && io==1) continue; // Only one ring in FMD1 
	if(tb == 1 && io==1) min=10;
	char     sio = (io == 0 ? 'I' : 'O');
	unsigned nio = (io == 0 ? 3   : 9);
	unsigned max = (io == 0 ? 5   : 10) + min;
	
	for(size_t i = min; i < max; i++) { // Modules
	  TString modVol(Form("%s/F%c%cV_7%d/F%cSE_%d", halfVol.Data(), 
			      sio, stb, nio, sio, i));
	  TString modSym(modVol);
	  if(!gGeoManager->SetAlignableEntry(modSym.Data(),modVol.Data()))
	    AliFatal(Form("Alignable entry %s not created. "
			  "Volume path %s not valid", 
			  modSym.Data(), modVol.Data()));
	}
      }
    }
  }
#endif
}
//___________________________________________________________________
//
// EOF
//
