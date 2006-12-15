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
#include <cmath>                // __CMATH__
#include <TClonesArray.h>	// ROOT_TClonesArray
#include <TGeometry.h>		// ROOT_TGeomtry
#include <TNode.h>		// ROOT_TNode
#include <TXTRU.h>		// ROOT_TXTRU
#include <TRotMatrix.h>		// ROOT_TRotMatrix
#include <TTUBE.h>		// ROOT_TTUBE
#include <TTree.h>		// ROOT_TTree
#include <TBrowser.h>		// ROOT_TBrowser
// #include <TVirtualMC.h>	// ROOT_TVirtualMC
#include <TVector2.h>           // ROOT_TVector2 
#include <TGeoManager.h>        // ROOT_TGeoManager

#include <AliRunDigitizer.h>	// ALIRUNDIGITIZER_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRun.h>		// ALIRUN_H
#include <AliMC.h>		// ALIMC_H
#include <AliMagF.h>		// ALIMAGF_H
#include <AliLog.h>		// ALILOG_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDSDigit.h"	// ALIFMDSDIGIT_H
#include "AliFMDHit.h"		// ALIFMDHIT_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMDDigitizer.h"	// ALIFMDDIGITIZER_H
#include "AliFMDSDigitizer.h"	// ALIFMDSDIGITIZER_H
// #include "AliFMDGeometryBuilder.h"
#include "AliFMDRawWriter.h"	// ALIFMDRAWWRITER_H
#include "AliFMDPoints.h"       // ALIFMDPOINTS_H

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
  AliDebug(10, "\tDefault CTOR");
  fHits        = 0;
  fDigits      = 0;
  fIshunt      = 0;
  fBad         = new TClonesArray("AliFMDHit");
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
  AliDebug(10, "\tStandard CTOR");
  fBad         = new TClonesArray("AliFMDHit");
  
  // Initialise Hit array
  HitsArray();
  gAlice->GetMCApp()->AddHitList(fHits);

  // (S)Digits for the detectors disk
  DigitsArray();
  SDigitsArray();
  
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
  AliDebug(10, "\tCreating materials");
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
  Int_t    fieldType        = gAlice->Field()->Integ();     // Field type 
  Double_t maxField         = gAlice->Field()->Max();     // Field max.
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

//____________________________________________________________________
void  
AliFMD::Init()
{
  // Initialize the detector 
  // 
  AliDebug(1, "Initialising FMD detector object");
  // AliFMDGeometry*  fmd = AliFMDGeometry::Instance();
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
    AliWarning((Form("EndEvent", "got %d 'bad' hits", fBad->GetEntries())));
    TIter next(fBad);
    AliFMDHit* hit;
    while ((hit = static_cast<AliFMDHit*>(next()))) hit->Print("D");
    fBad->Clear();
  }
}


//====================================================================
//
// Graphics and event display
//
//____________________________________________________________________
void 
AliFMD::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display. With the new
  // geometry modeller, TGeoManager, this seems rather redundant. 
  AliDebug(10, "\tCreating a simplified geometry");

  AliFMDGeometry* fmd = AliFMDGeometry::Instance();
  
  static TXTRU*     innerShape = 0;
  static TXTRU*     outerShape = 0;
  static TObjArray* innerRot   = 0;
  static TObjArray* outerRot   = 0;

  if (!innerShape || !outerShape) {
    // Make the shapes for the modules 
    for (Int_t i = 0; i < 2; i++) {
      AliFMDRing* r = 0;
      switch (i) {
      case 0: r = fmd->GetRing('I'); break;
      case 1: r = fmd->GetRing('O'); break;
      }
      if (!r) {
	AliError(Form("no ring found for i=%d", i));
	return;
      }
      Double_t    siThick  = r->GetSiThickness();
      const Int_t knv      = r->GetNVerticies();
      Double_t    theta    = r->GetTheta();
      Int_t       nmod     = r->GetNModules();
      
      TXTRU* shape = new TXTRU(r->GetName(), r->GetTitle(), "void", knv, 2);
      for (Int_t j = 0; j < knv; j++) {
	TVector2* vv = r->GetVertex(knv - 1 - j);
	shape->DefineVertex(j, vv->X(), vv->Y());
      }
      shape->DefineSection(0, -siThick / 2, 1, 0, 0);
      shape->DefineSection(1, +siThick / 2, 1, 0, 0);
      shape->SetLineColor(kYellow); //PH kYellow is the default line color in FMD
      
      TObjArray* rots = new TObjArray(nmod);
      for (Int_t j = 0; j < nmod; j++) {
	Double_t th = (j + .5) * theta * 2;
	TString name(Form("FMD_ring_%c_rot_%02d", r->GetId(), j));
	TString title(Form("FMD Ring %c Rotation # %d", r->GetId(), j));
	TRotMatrix* rot = new TRotMatrix(name.Data(), title.Data(),
					 90, th, 90, fmod(90+th,360), 0, 0);
	rots->AddAt(rot, j);
      }
      
      switch (r->GetId()) {
      case 'i':
      case 'I': innerShape = shape; innerRot = rots; break;
      case 'o':
      case 'O': outerShape = shape; outerRot = rots; break;
      }
    }
  }
  
  TNode* top = gAlice->GetGeometry()->GetNode("alice");
  
  for (Int_t i = 1; i <= 3; i++) {
    AliFMDDetector* det = fmd->GetDetector(i);
    if (!det) {
      Warning("BuildGeometry", "FMD%d seems to be disabled", i);
      continue;
    }
    Double_t w  = 0;
    Double_t rh = det->GetRing('I')->GetHighR();
    Char_t   id = 'I';
    if (det->GetRing('O')) {
      w  = TMath::Abs(det->GetRingZ('O') - det->GetRingZ('I'));
      id = (TMath::Abs(det->GetRingZ('O')) 
	    > TMath::Abs(det->GetRingZ('I')) ? 'O' : 'I');
      rh = det->GetRing('O')->GetHighR();
    }
    w += (det->GetRing(id)->GetModuleSpacing() +
	  det->GetRing(id)->GetSiThickness());
    TShape* shape = new TTUBE(det->GetName(), det->GetTitle(), "void",
			      det->GetRing('I')->GetLowR(), rh, w / 2);
    Double_t z = (det->GetRingZ('I') - w / 2);
    if (z > 0) z += det->GetRing(id)->GetModuleSpacing();
    top->cd();
    TNode* node = new TNode(det->GetName(), det->GetTitle(), shape, 
			    0, 0, z, 0);
    fNodes->Add(node);
    
    for (Int_t j = 0; j < 2; j++) {
      AliFMDRing* r      = 0;
      TShape*     rshape = 0;
      TObjArray*  rots   = 0;
      switch (j) {
      case 0: 
	r = det->GetRing('I'); rshape = innerShape; rots = innerRot; break;
      case 1: 
	r = det->GetRing('O'); rshape = outerShape; rots = outerRot; break;
      }
      if (!r) continue;
      
      Double_t    siThick  = r->GetSiThickness();
      Int_t       nmod     = r->GetNModules();
      Double_t    modspace = r->GetModuleSpacing();
      Double_t    rz       = - (z - det->GetRingZ(r->GetId()));
      
      for (Int_t k = 0; k < nmod; k++) {
	node->cd();
	Double_t    offz    = (k % 2 == 1 ? modspace : 0);
	TRotMatrix* rot     = static_cast<TRotMatrix*>(rots->At(k));
	TString name(Form("%s%c_module_%02d", det->GetName(), r->GetId(),k));
	TString title(Form("%s%c Module %d", det->GetName(), r->GetId(),k));
	TNode* mnod = new TNode(name.Data(), title.Data(), rshape, 
				0, 0, rz - siThick / 2 
				+ TMath::Sign(offz,z), rot);
	mnod->SetLineColor(kYellow); //PH kYellow is the default line color in FMD
	fNodes->Add(mnod);
      } // for (Int_t k = 0 ; ...)
    } // for (Int_t j = 0 ; ...)
  } // for (Int_t i = 1 ; ...)
}

//____________________________________________________________________
void 
AliFMD::LoadPoints(Int_t /* track */) 
{
  // Store x, y, z of all hits in memory for display. 
  // 
  // Normally, the hits are drawn using TPolyMarker3D - however, that
  // is not very useful for the FMD.  Therefor, this member function
  // is overloaded to make TMarker3D, via the class AliFMDPoints.
  // AliFMDPoints is a local class. 
  //
  if (!fHits) {
    AliError(Form("fHits == 0. Name is %s",GetName()));
    return;
  }
  Int_t nHits = fHits->GetEntriesFast();
  if (nHits == 0) {
    return;
  }
  Int_t tracks = gAlice->GetMCApp()->GetNtrack();
  if (fPoints == 0) fPoints = new TObjArray(2 * tracks);

  // Get geometry 
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();

  // Now make markers for each hit  
  // AliInfo(Form("Drawing %d hits (have %d points) for track %d", 
  //              nHits, fPoints->GetEntriesFast(), track));
  for (Int_t ihit = 0; ihit < nHits; ihit++) {
    AliFMDHit* hit = static_cast<AliFMDHit*>(fHits->At(ihit));
    if (!hit) continue;
    Double_t edep    = hit->Edep();
    Double_t m       = hit->M();
    Double_t poverm  = (m == 0 ? 0 : hit->P());
    Double_t absQ    = TMath::Abs(hit->Q());
    Bool_t   bad     = kFALSE;
    // This `if' is to debug abnormal energy depositions.  We trigger on
    // p/m approx larger than or equal to a MIP, and a large edep - more 
    // than 1 keV - a MIP is 100 eV. 
    if (edep > absQ * absQ && poverm > 1) bad = kTRUE;

    AliFMDPoints* p1 = new AliFMDPoints(hit, kRed); //PH kRed is the default marker color in FMD
    // AliPoints* p1 = new AliPoints();
    // p1->SetMarkerColor(GetMarkerColor());
    // p1->SetMarkerSize(GetMarkerSize());
    // p1->SetPoint(0, hit->X(), hit->Y(), hit->Z());
    p1->SetDetector(this);
    p1->SetParticle(hit->GetTrack());
    fPoints->AddAt(p1, hit->GetTrack());
    if (bad) {
      p1->SetMarkerColor(4);
      // p1->SetMarkerSize(2 * GetMarkerSize());
    }
    
    Double_t x, y, z;
    geom->Detector2XYZ(hit->Detector(), hit->Ring(), hit->Sector(), 
		       hit->Strip(), x, y, z);
    AliFMDPoints* p = new AliFMDPoints(hit, 3);
    // AliPoints* p = new AliPoints();
    // p->SetMarkerColor(3);
    // p->SetMarkerSize(GetMarkerSize());
    // p->SetPoint(0, x, y, z);
    p->SetDetector(this);
    p->SetParticle(hit->GetTrack());
    p->SetXYZ(x, y, z);
    p->SetMarkerColor(3);
    fPoints->AddAt(p, tracks+hit->GetTrack());
    if (bad) {
      p->SetMarkerColor(5);
      // p->SetMarkerSize(2 * GetMarkerSize());
    }
    // AliInfo(Form("Adding point at %d", tracks+hit->GetTrack()));
  }
}

//____________________________________________________________________
void 
AliFMD::DrawDetector()
{
  // Draw a shaded view of the Forward multiplicity detector.  This
  // isn't really useful anymore. 
  AliDebug(10, "\tDraw detector");
}

//____________________________________________________________________
Int_t 
AliFMD::DistancetoPrimitive(Int_t, Int_t)
{
  // Calculate the distance from the mouse to the FMD on the screen
  // Dummy routine.
  //
  return 9999;
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
      AliDebug(1, Form("already had a hit in FMD%d%c[%2d,%3d] for track # %d,"
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
  fNhits++;
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
		   Short_t(digits[6]));  // ADC Count3 
}

//____________________________________________________________________
void 
AliFMD::AddDigitByFields(UShort_t detector, 
			 Char_t   ring, 
			 UShort_t sector, 
			 UShort_t strip, 
			 UShort_t count1, 
			 Short_t  count2,
			 Short_t  count3)
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
  
  new (a[fNdigits++]) 
    AliFMDDigit(detector, ring, sector, strip, count1, count2, count3);
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
  AddSDigitByFields(UShort_t(digits[0]),  // Detector #
		    Char_t(digits[1]),    // Ring ID
		    UShort_t(digits[2]),  // Sector #
		    UShort_t(digits[3]),  // Strip #
		    Float_t(digits[4]),   // Edep
		    UShort_t(digits[5]),  // ADC Count1 
		    Short_t(digits[6]),   // ADC Count2 
		    Short_t(digits[7]));  // ADC Count3 
}

//____________________________________________________________________
void 
AliFMD::AddSDigitByFields(UShort_t detector, 
			  Char_t   ring, 
			  UShort_t sector, 
			  UShort_t strip, 
			  Float_t  edep,
			  UShort_t count1, 
			  Short_t  count2,
			  Short_t  count3)
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
  
  new (a[fNsdigits++]) 
    AliFMDSDigit(detector, ring, sector, strip, edep, count1, count2, count3);
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
  Warning("Hits2Digits", "Try not to use this method.\n"
	  "Instead, use AliSimulator");
  AliRunDigitizer* manager = new AliRunDigitizer(1, 1);
  manager->SetInputStream(0, "galice.root");
  manager->SetOutputFile("H2Dfile");
  
  /* AliDigitizer* dig =*/ CreateDigitizer(manager);
  manager->Exec("");
  delete manager;
}

//____________________________________________________________________
void 
AliFMD::Hits2SDigits() 
{
  // Create AliFMDSDigit's from AliFMDHit's.  This is done by creating
  // an AliFMDSDigitizer object, and executing it. 
  // 
  AliFMDSDigitizer* digitizer = new AliFMDSDigitizer("galice.root");
  digitizer->Exec("");
  delete digitizer;
}

  
//____________________________________________________________________
AliDigitizer* 
AliFMD::CreateDigitizer(AliRunDigitizer* manager) const
{
  // Create a digitizer object 
  AliFMDDigitizer* digitizer = new AliFMDDigitizer(manager);
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
// Utility 
//
//__________________________________________________________________
void 
AliFMD::Browse(TBrowser* b) 
{
  // Browse this object. 
  //
  AliDebug(30, "\tBrowsing the FMD");
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
