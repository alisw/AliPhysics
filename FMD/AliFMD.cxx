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
// Forward Multiplicity Detector based on Silicon wafers. This class
// contains the base procedures for the Forward Multiplicity detector
// Detector consists of 5 Si volumes covered pseudorapidity interval
// from 1.7 to 5.1.
//                                                       
// This is the base class for all FMD manager classes. 
//                    
// The actual code is done by various separate classes.   Below is
// diagram showing the relationship between the various FMD classes
// that handles the geometry 
//
//
//       +----------+   +----------+   
//       | AliFMDv1 |	| AliFMDv1 |   
//       +----------+   +----------+   
//            |              |
//       +----+--------------+
//       |
//       |           +------------+ 1  +---------------+
//       |        +- | AliFMDRing |<>--| AliFMDPolygon | 
//       V     2  |  +------------+    +---------------+   
//  +--------+<>--+        |
//  | AliFMD |             ^                       
//  +--------+<>--+        V 1..2                     
//	       3  | +-------------------+ 
//	          +-| AliFMDSubDetector | 
//	  	    +-------------------+
//                           ^              
//                           |
//             +-------------+-------------+
//             |             |             |	      
//        +---------+   +---------+   +---------+
//        | AliFMD1 |   | AliFMD2 |   | AliFMD3 |
//        +---------+   +---------+   +---------+
//      
//
// *  AliFMD 
//    This defines the interface for the various parts of AliROOT that
//    uses the FMD, like AliFMDDigitizer, AliFMDReconstructor, and so
//    on. 
//
// *  AliFMDv1 
//    This is a concrete implementation of the AliFMD interface. 
//    It is the responsibility of this class to create the FMD
//    geometry, process hits in the FMD, and serve hits and digits to
//    the various clients. 
//  
//    It uses the objects of class AliFMDSubDetector to do the various
//    stuff for FMD1, 2, and 3 
//
// *  AliFMDRing 
//    This class contains all stuff needed to do with a ring.  It's
//    used by the AliFMDSubDetector objects to instantise inner and
//    outer rings.  The AliFMDRing objects are shared by the
//    AliFMDSubDetector objects, and owned by the AliFMDv1 object. 
//
// *  AliFMDPolygon 
//    The code I lifted from TGeoPolygon to help with the geometry of
//    the modules, as well as to decide wether a hit is actually with
//    in the real module shape.  The point is, that the shape of the
//    various ring modules are really polygons (much like the lid of a
//    coffin), but it's segmented at constant radius.  That is very
//    hard to implement using GEANT 3.21 shapes, so instead the
//    modules are implemented as TUBS (tube sections), and in the step
//    procedure we do the test whether the track was inside the real
//    shape of the module.  
//
// *  AliFMD1, AliFMD2, and AliFMD3 
//    These are specialisation of AliFMDSubDetector, that contains the
//    particularities of each of the sub-detector system.  It is
//    envisioned that the classes should also define the support
//    volumes and material for each of the detectors.                          
//                                                                          
// The responsible person for this module is Alla Maevskaia
// <Alla.Maevskaia@cern.ch>.
//
// Many modifications by Christian Holm Christensen <cholm@nbi.dk>
//

#include "TClonesArray.h"	// ROOT_TClonesArray
#include "TGeometry.h"		// ROOT_TGeomtry
#include "TNode.h"		// ROOT_TNode
#include "TTUBE.h"		// ROOT_TTUBE
#include "TTree.h"		// ROOT_TTree
#include "TVirtualMC.h"		// ROOT_TVirtualMC
#include "TBrowser.h"		// ROOT_TBrowser
#include "TMath.h"		// ROOT_TMath

#include "AliRunDigitizer.h"	// ALIRUNDIGITIZER_H
#include "AliLoader.h"		// ALILOADER_H
#include "AliRun.h"		// ALIRUN_H
#include "AliMC.h"		// ALIMC_H
#include "AliLog.h"		// ALILOG_H
#include "AliMagF.h"		// ALIMAGF_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDDigit.h"	// ALIFMDDIGIG_H
#include "AliFMDHit.h"		// ALIFMDHIT_H
#include "AliFMDDigitizer.h"	// ALIFMDDIGITIZER_H
#include "AliFMD1.h"		// ALIFMD1_H
#include "AliFMD2.h"		// ALIFMD2_H
#include "AliFMD3.h"		// ALIFMD3_H
#include "AliFMDRawWriter.h"	// ALIFMDRAWWRITER_H

//____________________________________________________________________
ClassImp(AliFMD);

//____________________________________________________________________
AliFMD::AliFMD()
  : fInner(0), 
    fOuter(0),
    fFMD1(0),
    fFMD2(0), 
    fFMD3(0), 
    fSDigits(0), 
    fNsdigits(0),
    fSiDensity(0),
    fPrintboardRotationId(0),
    fIdentityRotationId(0),
    fShortLegId(0),
    fLongLegId(0),
    fLegLength(0),
    fLegRadius(0),
    fModuleSpacing(0)
{
  //
  // Default constructor for class AliFMD
  //
  AliDebug(0, "Default CTOR");
  fHits     = 0;
  fDigits   = 0;
  fIshunt   = 0;
}

//____________________________________________________________________
AliFMD::AliFMD(const char *name, const char *title, bool detailed)
  : AliDetector (name, title),
    fInner(0), 
    fOuter(0),
    fFMD1(0),
    fFMD2(0), 
    fFMD3(0),
    fSDigits(0),
    fNsdigits(0),
    fSiDensity(0),
    fPrintboardRotationId(0),
    fIdentityRotationId(0),
    fShortLegId(0),
    fLongLegId(0),
    fLegLength(0),
    fLegRadius(0),
    fModuleSpacing(0)
{
  //
  // Standard constructor for Forward Multiplicity Detector
  //
  AliDebug(0, "Standard CTOR");

  // Initialise Hit array
  HitsArray();
  gAlice->GetMCApp()->AddHitList(fHits);

  // (S)Digits for the detectors disk
  DigitsArray();
  SDigitsArray();
  
  // CHC: What is this?
  fIshunt = 0;
  SetMarkerColor(kRed);
  SetLineColor(kYellow);
  SetSiDensity();

  // Create sub-volume managers 
  fInner = new AliFMDRing('I', detailed);
  fOuter = new AliFMDRing('O', detailed);
  fFMD1  = new AliFMD1();
  fFMD2  = new AliFMD2();
  fFMD3  = new AliFMD3();

  // Specify parameters of sub-volume managers 
  fFMD1->SetInner(fInner);
  fFMD1->SetOuter(0);

  fFMD2->SetInner(fInner);
  fFMD2->SetOuter(fOuter);
  
  fFMD3->SetInner(fInner);
  fFMD3->SetOuter(fOuter);

  SetLegLength();
  SetLegRadius();
  SetLegOffset();
  SetModuleSpacing();
  
  fInner->SetLowR(4.3);
  fInner->SetHighR(17.2);
  fInner->SetWaferRadius(13.4/2);
  fInner->SetTheta(36/2);
  fInner->SetNStrips(512);
  fInner->SetSiThickness(.03);
  fInner->SetPrintboardThickness(.11);
  fInner->SetBondingWidth(.5);

  fOuter->SetLowR(15.6);
  fOuter->SetHighR(28.0);
  fOuter->SetWaferRadius(13.4/2);
  fOuter->SetTheta(18/2);
  fOuter->SetNStrips( 256);
  fOuter->SetSiThickness(.03);
  fOuter->SetPrintboardThickness(.1);
  fOuter->SetBondingWidth(.5);
  
  
  fFMD1->SetHoneycombThickness(1);
  fFMD1->SetInnerZ(340.0);
  
  fFMD2->SetHoneycombThickness(1);
  fFMD2->SetInnerZ(83.4);
  fFMD2->SetOuterZ(75.2);

  fFMD3->SetHoneycombThickness(1);
  fFMD3->SetInnerZ(-62.8);
  fFMD3->SetOuterZ(-75.2);
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
  // construction of the geometry is delegated to the class AliFMDRing
  // and AliFMDSubDetector and the relevant derived classes. 
  //
  // The flow of this member function is:
  //
  //   FOR rings fInner and fOuter DO  
  //     AliFMDRing::Init();
  //   END FOR
  // 
  //   Set up hybrud card support (leg) volume shapes  
  // 
  //   FOR rings fInner and fOuter DO  
  //     AliFMDRing::SetupGeometry();
  //   END FOR
  // 
  //   FOR subdetectors fFMD1, fFMD2, and fFMD3 DO 
  //     AliFMDSubDetector::SetupGeomtry();
  //   END FOR
  // 
  //   FOR subdetectors fFMD1, fFMD2, and fFMD3 DO 
  //     AliFMDSubDetector::Geomtry();
  //   END FOR
  //   

  // DebugGuard guard("AliFMD::CreateGeometry");
  AliDebug(10, "Creating geometry");

  fInner->Init();
  fOuter->Init();

  TString name;
  Double_t par[3];

  par[0]      =  fLegRadius - .1;
  par[1]      =  fLegRadius;
  par[2]      =  fLegLength / 2;
  name        =  "FSL";
  fShortLegId =  gMC->Gsvolu(name.Data(),"TUBE",(*fIdtmed)[kPlasticId],par,3);
  
  par[2]      += fModuleSpacing / 2;
  name        = "FLL";
  fLongLegId  =  gMC->Gsvolu(name.Data(),"TUBE",(*fIdtmed)[kPlasticId],par,3);

  fInner->SetupGeometry((*fIdtmed)[kAirId], 
			(*fIdtmed)[kSiId], 
			(*fIdtmed)[kPcbId], 
			fPrintboardRotationId, 
			fIdentityRotationId);
  fOuter->SetupGeometry((*fIdtmed)[kAirId], 
			(*fIdtmed)[kSiId], 
			(*fIdtmed)[kPcbId], 
			fPrintboardRotationId, 
			fIdentityRotationId);

  fFMD1->SetupGeometry((*fIdtmed)[kAirId], (*fIdtmed)[kKaptionId]);
  fFMD2->SetupGeometry((*fIdtmed)[kAirId], (*fIdtmed)[kKaptionId]);
  fFMD3->SetupGeometry((*fIdtmed)[kAirId], (*fIdtmed)[kKaptionId]);
  
  fFMD1->Geometry("ALIC", fPrintboardRotationId, fIdentityRotationId);
  fFMD2->Geometry("ALIC", fPrintboardRotationId, fIdentityRotationId);
  fFMD3->Geometry("ALIC", fPrintboardRotationId, fIdentityRotationId);    
}    

//____________________________________________________________________
void AliFMD::CreateMaterials() 
{
  // Register various materials and tracking mediums with the
  // backend.   
  // 
  // Currently defined materials and mediums are 
  // 
  //    FMD Air		Normal air 
  //    FMD Si          Active silicon of sensors 
  //    FMD Carbon      Normal carbon used in support, etc. 
  //    FMD Kapton      Carbon used in Honeycomb
  //    FMD PCB         Printed circuit board material 
  //    FMD Plastic     Material for support legs 
  // 
  // Also defined are two rotation matricies. 
  //
  // DebugGuard guard("AliFMD::CreateMaterials");
  AliDebug(10, "Creating materials");
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
  density          = fSiDensity;
  radiationLength  = 9.36;
  maxBending       = 1;
  maxStepSize      = .001;
  precision        = .001;
  minStepSize      = .001;
  id               = kSiId;
  AliMaterial(id, "FMD Si$", a, z, density, radiationLength, absorbtionLength);
  AliMedium(kSiId, "FMD Si$",id,1,fieldType,maxField,maxBending,
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
  AliMaterial(id, "FMD Carbon$", a, z, density, radiationLength, 
	      absorbtionLength);
  AliMedium(kCarbonId, "FMD Carbon$",id,0,fieldType,maxField,maxBending,
	    maxStepSize,maxEnergyLoss,precision,minStepSize);

  // Silicon chip 
  {
    Float_t as[] = { 12.0107,      14.0067,      15.9994,
		     1.00794,      28.0855,     107.8682 };
    Float_t zs[] = {  6.,           7.,           8.,
		      1.,          14.,          47. };
    Float_t ws[] = {  0.039730642,  0.001396798,  0.01169634,
		      0.004367771,  0.844665,     0.09814344903 };
    density = 2.36436;
    maxBending       = 10;
    maxStepSize      = .01;
    precision        = .003;
    minStepSize      = .003;
    id = kSiChipId;
    AliMixture(id, "FMD Si Chip$", as, zs, density, 6, ws);
    AliMedium(kSiChipId, "FMD Si Chip$", id, 0, fieldType, maxField, 
	      maxBending, maxStepSize, maxEnergyLoss, precision, minStepSize);
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
    id               = kKaptionId;
    AliMixture(id, "FMD Kaption$", as, zs, density, 4, ws);
    AliMedium(kKaptionId, "FMD Kaption$",id,0,fieldType,maxField,maxBending,
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
    AliMixture(id, "FMD Air$", as, zs, density, 4, ws);
    AliMedium(kAirId, "FMD Air$", id,0,fieldType,maxField,maxBending,
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
    AliMixture(id, "FMD PCB$", as, zs, density, 14, ws);
    AliMedium(kPcbId, "FMD PCB$", id,1,fieldType,maxField,maxBending,
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
    AliMixture(id, "FMD Plastic$", as, zs, density, -2, ws);
    AliMedium(kPlasticId, "FMD Plastic$", id,0,fieldType,maxField,maxBending,
		maxStepSize,maxEnergyLoss,precision,minStepSize);
  }
  AliMatrix(fPrintboardRotationId, 90, 90, 0, 90, 90, 0);
  AliMatrix(fIdentityRotationId, 90, 0, 90, 90, 0, 0);
}

//____________________________________________________________________
void  
AliFMD::Init()
{
  //
  // Initialis the FMD after it has been built
  Int_t i;
  //
  if (fDebug) {
    cout << "\n" << ClassName() << ": " << flush;
    for (i = 0; i < 35; i++) cout << "*";
    cout << " FMD_INIT ";
    for (i = 0; i < 35; i++) cout << "*";
    cout << "\n" << ClassName() << ": " << flush;
    //
    // Here the FMD initialisation code (if any!)
    for (i = 0; i < 80; i++) cout << "*";
    cout << endl;
  }
  //
  //
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
  // Build simple ROOT TNode geometry for event display
  //
  // Build a simplified geometry of the FMD used for event display  
  // 
  // The actual building of the TNodes is done by
  // AliFMDSubDetector::SimpleGeometry. 
  AliDebug(10, "Creating a simplified geometry");

  TNode* top = gAlice->GetGeometry()->GetNode("alice");
  
  fFMD1->SimpleGeometry(fNodes, top, GetLineColor(), 0);
  fFMD2->SimpleGeometry(fNodes, top, GetLineColor(), 0);
  fFMD3->SimpleGeometry(fNodes, top, GetLineColor(), 0);
}

//____________________________________________________________________
void 
AliFMD::DrawDetector()
{
  //
  // Draw a shaded view of the Forward multiplicity detector
  //
  // DebugGuard guard("AliFMD::DrawDetector");
  AliDebug(10, "Draw detector");
  
  //Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);

  //Set volumes visible
  fFMD1->Gsatt();
  fFMD2->Gsatt();
  fFMD3->Gsatt();
  fInner->Gsatt();
  fOuter->Gsatt();

  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 12, .055, .055);
  gMC->Gdhead(1111, "Forward Multiplicity Detector");
  gMC->Gdman(16, 10, "MAN");
  gMC->Gdopt("hide", "off");
}

//____________________________________________________________________
const Int_t 
AliFMD::DistanceToPrimitive(Int_t, Int_t)
{
  //
  // Calculate the distance from the mouse to the FMD on the screen
  // Dummy routine
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
  AddHit(track, 
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
void 
AliFMD::AddHit(Int_t    track, 
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
	       Float_t  t)
{
  //
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
  // 
  TClonesArray& a = *(HitsArray());
  // Search through the list of already registered hits, and see if we
  // find a hit with the same parameters.  If we do, then don't create
  // a new hit, but rather update the energy deposited in the hit.
  // This is done, so that a FLUKA based simulation will get the
  // number of hits right, not just the enerrgy deposition. 
  for (Int_t i = 0; i < fNhits; i++) {
    if (!a.At(i)) continue;
    AliFMDHit* hit = static_cast<AliFMDHit*>(a.At(i));
    if (hit->Detector() == detector 
	&& hit->Ring() == ring
	&& hit->Sector() == sector 
	&& hit->Strip() == strip
	&& hit->Track() == track) {
      Warning("AddHit", "already had a hit in FMD%d%c[%2d,%3d] for track # %d,"
	      " adding energy (%f) to that hit (%f) -> %f", 
	      detector, ring, sector, strip, track, edep, hit->Edep(),
	      hit->Edep() + edep);
      hit->SetEdep(hit->Edep() + edep);
      return;
    }
  }
  // If hit wasn't already registered, do so know. 
  new (a[fNhits]) AliFMDHit(fIshunt, track, detector, ring, sector, strip, 
			    x, y, z, px, py, pz, edep, pdg, t);
  fNhits++;
}

//____________________________________________________________________
void 
AliFMD::AddDigit(Int_t* digits)
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
  AddDigit(UShort_t(digits[0]),  // Detector #
	   Char_t(digits[1]),    // Ring ID
	   UShort_t(digits[2]),  // Sector #
	   UShort_t(digits[3]),  // Strip #
	   UShort_t(digits[4]),  // ADC Count1 
	   Short_t(digits[5]), 	 // ADC Count2 
	   Short_t(digits[6]));  // ADC Count3 
}

//____________________________________________________________________
void 
AliFMD::AddDigit(UShort_t detector, 
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
  AddSDigit(UShort_t(digits[0]),  // Detector #
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
AliFMD::AddSDigit(UShort_t detector, 
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
  //
  // Reset number of digits and the digits array for this detector
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
  AliRunDigitizer* manager = new AliRunDigitizer(1, 1);
  manager->SetInputStream(0, "galice.root");
  manager->SetOutputFile("H2Dfile");
  
  /* AliDigitizer* dig =*/ CreateDigitizer(manager);
  manager->Exec("");
}

//____________________________________________________________________
void 
AliFMD::Hits2SDigits() 
{
  // Create AliFMDSDigit's from AliFMDHit's.  This is done by creating
  // an AliFMDSDigitizer object, and executing it. 
  // 
  AliDigitizer* sdig = new AliFMDSDigitizer("galice.root");
  sdig->Exec("");
}

  
//____________________________________________________________________
AliDigitizer* 
AliFMD::CreateDigitizer(AliRunDigitizer* manager) const
{
  // Create a digitizer object 
  return new AliFMDDigitizer(manager);
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
  
#if 0
  // Digits are read from the Digit branch, and processed to make
  // three DDL files, one for each of the sub-detectors FMD1, FMD2,
  // and FMD3. 
  //
  // The raw data files consists of a header, followed by ALTRO
  // formatted blocks.  
  // 
  //          +-------------+
  //          | Header      |
  //          +-------------+
  //          | ALTRO Block |
  //          | ...         |
  //          +-------------+
  //          DDL file 
  // 
  // An ALTRO formatted block, in the FMD context, consists of a
  // number of counts followed by a trailer. 
  // 
  //          +------------------+
  //          | Count            |
  //          | ...              |
  //          | possible fillers |
  //          +------------------+
  //          | Trailer          |
  //          +------------------+
  //          ALTRO block 
  // 
  // The counts are listed backwards, that is, starting with the
  // latest count, and ending in the first. 
  // 
  // Each count consist of 1 or more ADC samples of the VA1_ALICE
  // pre-amp. signal.  Just how many samples are used depends on
  // whether the ALTRO over samples the pre-amp.  Each sample is a
  // 10-bit word, and the samples are grouped into 40-bit blocks 
  //
  //          +------------------------------------+
  //          |  S(n)   | S(n-1) | S(n-2) | S(n-3) |
  //          |  ...    | ...    | ...    | ...    |
  //          |  S(2)   | S(1)   | AA     | AA     |
  //          +------------------------------------+
  //          Counts + possible filler 
  //
  // The trailer of the number of words of signales, the starting
  // strip number, the sector number, and the ring ID; each 10-bit
  // words,  packed into 40-bits. 
  // 
  //          +------------------------------------+
  //          | # words | start  | sector | ring   |
  //          +------------------------------------+
  //          Trailer
  // 
  // Note, that this method assumes that the digits are ordered. 
  //
  AliFMD* fmd = static_cast<AliFMD*>(gAlice->GetDetector(GetName()));
  fLoader->LoadDigits();
  TTree* digitTree = fLoader->TreeD();
  if (!digitTree) {
    Error("Digits2Raw", "no digit tree");
    return;
  }
  
  TClonesArray* digits = new TClonesArray("AliFMDDigit", 1000);
  fmd->SetTreeAddress();
  TBranch* digitBranch = digitTree->GetBranch(GetName());
  if (!digitBranch) {
    Error("Digits2Raw", "no branch for %s", GetName());
    return;
  }
  digitBranch->SetAddress(&digits);
  
  Int_t nEvents = Int_t(digitTree->GetEntries());
  for (Int_t event = 0; event < nEvents; event++) {
    fmd->ResetDigits();
    digitTree->GetEvent(event);
    
    Int_t nDigits = digits->GetEntries();
    if (nDigits < 1) continue;


    UShort_t prevDetector = 0;
    Char_t   prevRing     = '\0';
    UShort_t prevSector   = 0;
    // UShort_t prevStrip    = 0;

    // The first seen strip number for a channel 
    UShort_t startStrip   = 0;
    
    // Which channel number in the ALTRO channel we're at 
    UShort_t offset       = 0;

    // How many times the ALTRO Samples one VA1_ALICE channel 
    Int_t sampleRate = 1;

    // A buffer to hold 1 ALTRO channel - Normally, one ALTRO channel
    // holds 128 VA1_ALICE channels, sampled at a rate of `sampleRate' 
    TArrayI channel(128 * sampleRate);
    
    // The Altro buffer 
    AliAltroBuffer* altro = 0;
    
    // Loop over the digits in the event.  Note, that we assume the
    // the digits are in order in the branch.   If they were not, we'd
    // have to cache all channels before we could write the data to
    // the ALTRO buffer, or we'd have to set up a map of the digits. 
    for (Int_t i = 0; i < nDigits; i++) {
      // Get the digit
      AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));

      UShort_t det    = digit->Detector();
      Char_t   ring   = digit->Ring();
      UShort_t sector = digit->Sector();
      UShort_t strip  = digit->Strip();
      if (det != prevDetector) {
	AliDebug(10, Form("FMD: New DDL, was %d, now %d",
			  kBaseDDL + prevDetector - 1,
			  kBaseDDL + det - 1));
	// If an altro exists, delete the object, flushing the data to
	// disk, and closing the file. 
	if (altro) { 
	  // When the first argument is false, we write the real
	  // header. 
	  AliDebug(10, Form("New altro: Write channel at %d Strip: %d "
			    "Sector: %d  Ring: %d", 
			    i, startStrip, prevSector, prevRing));
	  // TPC to FMD translations 
	  // 
	  //    TPC                FMD
	  //    ----------+-----------
	  //    pad       |      strip
	  //    row       |     sector
	  //    sector    |       ring
	  // 
	  altro->WriteChannel(Int_t(startStrip), 
			      Int_t(prevSector), 
			      Int_t((prevRing == 'I' ? 0 : 1)), 
			      channel.fN, channel.fArray, 0);
	  altro->Flush();
	  altro->WriteDataHeader(kFALSE, kFALSE);
	  delete altro;
	  altro = 0;
	}

	prevDetector = det;
	// Need to open a new DDL! 
	Int_t ddlId = kBaseDDL + det - 1;
	TString filename(Form("%s_%d.ddl", GetName(),  ddlId));

	AliDebug(10, Form("New altro buffer with DDL file %s", 
			  filename.Data()));
	AliDebug(10, Form("New altro at %d", i));
	// Create a new altro buffer - a `1' as the second argument
	// means `write mode' 
	altro = new AliAltroBuffer(filename.Data(), 1);
	
	// Write a dummy (first argument is true) header to the DDL
	// file - later on, when we close the file, we write the real
	// header
	altro->WriteDataHeader(kTRUE, kFALSE);

	// Figure out the sample rate 
	if (digit->Count2() > 0) sampleRate = 2;
	if (digit->Count3() > 0) sampleRate = 3;

	channel.Set(128 * sampleRate);
	offset     = 0;
	prevRing   = ring;
	prevSector = sector;
	startStrip = strip;
      }
      else if (offset == 128                        
	       || digit->Ring() != prevRing 
	       || digit->Sector() != prevSector) {
	// Force a new Altro channel
	AliDebug(10, Form("Flushing channel to disk because %s",
			  (offset == 128 ? "channel is full" :
			   (ring != prevRing ? "new ring up" :
			    "new sector up"))));
	AliDebug(10, Form("New Channel: Write channel at %d Strip: %d "
			  "Sector: %d  Ring: %d", 
			  i, startStrip, prevSector, prevRing));
	altro->WriteChannel(Int_t(startStrip), 
			    Int_t(prevSector), 
			    Int_t((prevRing == 'I' ? 0 : 1)), 
			    channel.fN, channel.fArray, 0);
	// Reset and update channel variables 
	channel.Reset(0);
	offset     = 0; 
	startStrip = strip;
	prevRing   = ring;
	prevSector = sector;
      }

      // Store the counts of the ADC in the channel buffer 
      channel[offset * sampleRate] = digit->Count1();
      if (sampleRate > 1) 
	channel[offset * sampleRate + 1] = digit->Count2();
      if (sampleRate > 2) 
	channel[offset * sampleRate + 2] = digit->Count3();
      offset++;
    }
    // Finally, we need to close the final ALTRO buffer if it wasn't
    // already 
    if (altro) {
      altro->Flush();
      altro->WriteDataHeader(kFALSE, kFALSE);
      delete altro;
    }
  }
  fLoader->UnloadDigits();
#endif
}

//==================================================================
//
// Various setter functions for the common paramters 
//

//__________________________________________________________________
void 
AliFMD::SetLegLength(Double_t length) 
{
  // Set lenght of plastic legs that hold the hybrid (print board and
  // silicon sensor) onto the honeycomp support
  //
  // DebugGuard guard("AliFMD::SetLegLength");
  AliDebug(10, "AliFMD::SetLegLength");
  fLegLength = length;
  fInner->SetLegLength(fLegLength);
  fOuter->SetLegLength(fLegLength);
}

//__________________________________________________________________
void 
AliFMD::SetLegOffset(Double_t offset) 
{
  // Set offset from edge of hybrid to plastic legs that hold the
  // hybrid (print board and silicon sensor) onto the honeycomp
  // support 
  //
  // DebugGuard guard("AliFMD::SetLegOffset");
  AliDebug(10, "AliFMD::SetLegOffset");
  fInner->SetLegOffset(offset);
  fOuter->SetLegOffset(offset);
}

//__________________________________________________________________
void 
AliFMD::SetLegRadius(Double_t radius) 
{
  // Set the diameter of the plastic legs that hold the hybrid (print
  // board and silicon sensor) onto the honeycomp support
  //
  // DebugGuard guard("AliFMD::SetLegRadius");
  AliDebug(10, "AliFMD::SetLegRadius");
  fLegRadius = radius;
  fInner->SetLegRadius(fLegRadius);
  fOuter->SetLegRadius(fLegRadius);
}

//__________________________________________________________________
void 
AliFMD::SetModuleSpacing(Double_t spacing) 
{
  // Set the distance between the front and back sensor modules
  // (module staggering). 
  //
  // DebugGuard guard("AliFMD::SetModuleSpacing");
  AliDebug(10, "AliFMD::SetModuleSpacing");  
  fModuleSpacing = spacing;
  fInner->SetModuleSpacing(fModuleSpacing);
  fOuter->SetModuleSpacing(fModuleSpacing);
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
  AliDebug(10, "AliFMD::Browse");
  AliDetector::Browse(b);
  if (fInner) b->Add(fInner, "Inner Ring");
  if (fOuter) b->Add(fOuter, "Outer Ring");
  if (fFMD1)  b->Add(fFMD1,  "FMD1 SubDetector");
  if (fFMD2)  b->Add(fFMD2,  "FMD2 SubDetector");
  if (fFMD3)  b->Add(fFMD3,  "FMD3 SubDetector");
}


//___________________________________________________________________
//
// EOF
//
