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
// Detector consists of 3 sub-detectors FMD1, FMD2, and FMD3, each of
// which has 1 or 2 rings of silicon sensors. 
//                                                       
// This is the base class for all FMD manager classes. 
//                    
// The actual code is done by various separate classes.   Below is
// diagram showing the relationship between the various FMD classes
// that handles the simulation
//
//      +--------+ 1     +-----------------+ 
//      | AliFMD |<>-----| AliFMDSimulator |
//      +--------+	 +-----------------+
//                               ^              
//                               |
//                 +-------------+-------------+
//                 |                           |	      
//        +--------------------+   +-------------------+
//        | AliFMDGeoSimulator |   | AliFMDG3Simulator | 
//        +--------------------+   +-------------------+
//                                           ^
//                                           |
//                                 +----------------------+
//				   | AliFMDG3OldSimulator |
//				   +----------------------+
//      
// *  AliFMD 
//    This defines the interface for the various parts of AliROOT that
//    uses the FMD, like AliFMDSimulator, AliFMDDigitizer, 
//    AliFMDReconstructor, and so on. 
//
// *  AliFMDSimulator
//    This is the base class for the FMD simulation tasks.   The
//    simulator tasks are responsible to implment the geoemtry, and
//    process hits. 
//                                                                          
// *  AliFMDGeoSimulator
//    This is a concrete implementation of the AliFMDSimulator that
//    uses the TGeo classes directly only.  This defines the active
//    volume as an ONLY XTRU shape with a divided MANY TUBS shape
//    inside to implement the particular shape of the silicon
//    sensors. 
//
// *  AliFMDG3Simulator
//    This is a concrete implementation of the AliFMDSimulator that
//    uses the TVirtualMC interface with GEANT 3.21-like messages.
//    This implements the active volume as a divided TUBS shape.  Hits
//    in the corners should be cut away at run time (but currently
//    isn't). 
//
// *  AliFMDG3OldSimulator
//    This is a concrete implementation of AliFMDSimulator.   It
//    approximates the of the rings as segmented disks. 
// 
#include "AliFMDSimulator.h"	// ALIFMDSIMULATOR_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMD1.h"		// ALIFMD1_H
#include "AliFMD2.h"		// ALIFMD2_H
#include "AliFMD3.h"		// ALIFMD3_H
#include "AliFMD.h"		// ALIFMD_H
#include <AliRun.h>		// ALIRUN_H
#include <AliMC.h>		// ALIMC_H
#include <AliMagF.h>		// ALIMAGF_H
#include <AliLog.h>		// ALILOG_H
#include <TGeoVolume.h>		// ROOT_TGeoVolume
#include <TGeoTube.h>		// ROOT_TGeoTube
#include <TGeoPcon.h>		// ROOT_TGeoPcon
#include <TGeoMaterial.h>	// ROOT_TGeoMaterial
#include <TGeoMedium.h>		// ROOT_TGeoMedium
#include <TGeoXtru.h>		// ROOT_TGeoXtru
#include <TGeoPolygon.h>	// ROOT_TGeoPolygon
#include <TGeoTube.h>		// ROOT_TGeoTube
#include <TGeoManager.h>	// ROOT_TGeoManager
#include <TTree.h>		// ROOT_TTree
#include <TParticle.h>		// ROOT_TParticle
#include <TLorentzVector.h>	// ROOT_TLorentzVector
#include <TVector2.h>		// ROOT_TVector2
#include <TVector3.h>		// ROOT_TVector3
#include <TVirtualMC.h>		// ROOT_TVirtualMC
#include <TArrayD.h>		// ROOT_TArrayD

//====================================================================
ClassImp(AliFMDSimulator)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
const Char_t* AliFMDSimulator::fgkActiveName	= "F%cAC";
const Char_t* AliFMDSimulator::fgkSectorName	= "F%cSE";
const Char_t* AliFMDSimulator::fgkStripName	= "F%cST";
const Char_t* AliFMDSimulator::fgkModuleName	= "F%cMO";
const Char_t* AliFMDSimulator::fgkPCBName	= "F%cP%c";
const Char_t* AliFMDSimulator::fgkLongLegName	= "F%cLL";
const Char_t* AliFMDSimulator::fgkShortLegName	= "F%cSL";
const Char_t* AliFMDSimulator::fgkFrontVName	= "F%cFV";
const Char_t* AliFMDSimulator::fgkBackVName	= "F%cBV";
const Char_t* AliFMDSimulator::fgkRingName	= "FMD%c";
const Char_t* AliFMDSimulator::fgkTopHCName	= "F%d%cI";
const Char_t* AliFMDSimulator::fgkBotHCName	= "F%d%cJ";
const Char_t* AliFMDSimulator::fgkTopIHCName	= "F%d%cK";
const Char_t* AliFMDSimulator::fgkBotIHCName	= "F%d%cL";
const Char_t* AliFMDSimulator::fgkNoseName      = "F3SN";
const Char_t* AliFMDSimulator::fgkBackName      = "F3SB";
const Char_t* AliFMDSimulator::fgkBeamName      = "F3SL";
const Char_t* AliFMDSimulator::fgkFlangeName    = "F3SF";

//____________________________________________________________________
AliFMDSimulator::AliFMDSimulator() 
  : fFMD(0), 
    fDetailed(kFALSE),
    fInnerId(-1),
    fOuterId(-1), 
    fActiveId(4), 
    fUseDivided(kFALSE),
    fUseAssembly(kTRUE)
{
  // Default constructor
}

//____________________________________________________________________
AliFMDSimulator::AliFMDSimulator(AliFMD* fmd, Bool_t detailed) 
  : TTask("FMDSimulator", "Forward Multiplicity Detector Simulator"), 
    fFMD(fmd), 
    fDetailed(detailed),
    fInnerId(-1),
    fOuterId(-1),
    fActiveId(4),
    fUseDivided(kFALSE),
    fUseAssembly(kTRUE)
{
  // Normal constructor
  // 
  // Parameters: 
  // 
  //      fmd		Pointer to AliFMD object 
  //      detailed      Whether to make a detailed simulation or not 
  // 
}


//____________________________________________________________________
void
AliFMDSimulator::DefineMaterials() 
{
  // Define the materials and tracking mediums needed by the FMD
  // simulation.   These mediums are made by sending the messages
  // AliMaterial, AliMixture, and AliMedium to the passed AliModule
  // object module.   The defined mediums are 
  // 
  //	FMD Si$		Silicon (active medium in sensors)
  //	FMD C$		Carbon fibre (support cone for FMD3 and vacuum pipe)
  //	FMD Al$		Aluminium (honeycomb support plates)
  //	FMD PCB$	Printed Circuit Board (FEE board with VA1_ALICE)
  //	FMD Chip$	Electronics chips (currently not used)
  //	FMD Air$	Air (Air in the FMD)
  //	FMD Plastic$	Plastic (Support legs for the hybrid cards)
  //
  // Pointers to TGeoMedium objects are retrived from the TGeoManager
  // singleton.  These pointers are later used when setting up the
  // geometry 
  AliDebug(10, "\tCreating materials");
  // Get pointer to geometry singleton object. 
  AliFMDGeometry* geometry = AliFMDGeometry::Instance();
  geometry->Init();
  
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
  fFMD->AliMaterial(id, "Si$", 
		      a, z, density, radiationLength, absorbtionLength);
  fFMD->AliMedium(kSiId, "Si$",
		    id,1,fieldType,maxField,maxBending,
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
  fFMD->AliMaterial(id, "Carbon$", 
		      a, z, density, radiationLength, absorbtionLength);
  fFMD->AliMedium(kCarbonId, "Carbon$",
		    id,0,fieldType,maxField,maxBending,
		    maxStepSize,maxEnergyLoss,precision,minStepSize);

  // Aluminum
  a                = 26.981539;
  z                = 13.;
  density          = 2.7;
  radiationLength  = 8.9;
  id               = kAlId;
  fFMD->AliMaterial(id, "Aluminum$", 
		      a, z, density, radiationLength, absorbtionLength);
  fFMD->AliMedium(kAlId, "Aluminum$", 
		    id, 0, fieldType, maxField, maxBending,
		    maxStepSize, maxEnergyLoss, precision, minStepSize);
  
  
  // Copper 
  a                = 63.546;
  z                = 29;
  density          =  8.96;
  radiationLength  =  1.43;
  id               = kCopperId;
  fFMD->AliMaterial(id, "Copper$", 
		      a, z, density, radiationLength, absorbtionLength);
  fFMD->AliMedium(kCopperId, "Copper$", 
		    id, 0, fieldType, maxField, maxBending,
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
    fFMD->AliMixture(id, "Si Chip$", as, zs, density, 6, ws);
    fFMD->AliMedium(kSiChipId, "Si Chip$", 
		      id, 0, fieldType, maxField, maxBending, 
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
    fFMD->AliMixture(id, "Kaption$", as, zs, density, 4, ws);
    fFMD->AliMedium(kKaptonId, "Kaption$",
		      id,0,fieldType,maxField,maxBending,
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
    fFMD->AliMixture(id, "Air$", as, zs, density, 4, ws);
    fFMD->AliMedium(kAirId, "Air$", 
		      id,0,fieldType,maxField,maxBending,
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
    fFMD->AliMixture(id, "PCB$", as, zs, density, 14, ws);
    fFMD->AliMedium(kPcbId, "PCB$", 
		      id,0,fieldType,maxField,maxBending,
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
    fFMD->AliMixture(id, "Plastic$", as, zs, density, -2, ws);
    fFMD->AliMedium(kPlasticId, "Plastic$", 
		      id,0,fieldType,maxField,maxBending,
		      maxStepSize,maxEnergyLoss,precision,minStepSize);
  }
}

//____________________________________________________________________
Bool_t
AliFMDSimulator::IsActive(Int_t volId) const
{
  for (Int_t i = 0; i < fActiveId.fN; i++) 
    if  (volId == fActiveId[i]) return kTRUE;
  return kFALSE;
}

//____________________________________________________________________
Bool_t
AliFMDSimulator::VMC2FMD(TLorentzVector& v, UShort_t& detector,
                         Char_t& ring, UShort_t& sector, UShort_t& strip)
{
  TVirtualMC* mc = TVirtualMC::GetMC();

  // Get track position
  mc->TrackPosition(v);
  Int_t moduleno; mc->CurrentVolOffID(fModuleOff, moduleno);
  Int_t iring;    mc->CurrentVolOffID(fRingOff, iring);   ring = Char_t(iring);
  Int_t det;      mc->CurrentVolOffID(fDetectorOff, det); detector = det;
  

  // Get the ring geometry
  AliFMDGeometry*  fmd = AliFMDGeometry::Instance();
  //Int_t     nsec = fmd->GetDetector(detector)->GetRing(ring)->GetNSectors();
  Int_t     nstr = fmd->GetDetector(detector)->GetRing(ring)->GetNStrips();
  Double_t  lowr = fmd->GetDetector(detector)->GetRing(ring)->GetLowR();
  Double_t  highr= fmd->GetDetector(detector)->GetRing(ring)->GetHighR();
  Double_t  theta= fmd->GetDetector(detector)->GetRing(ring)->GetTheta();

  // Figure out the strip number
  Double_t r     = TMath::Sqrt(v.X() * v.X() + v.Y() * v.Y());
  Double_t pitch = (highr - lowr) / nstr;
  Int_t    str   = Int_t((r - lowr) / pitch);
  if (str < 0 || str >= nstr) return kFALSE;
  strip          = str;

  // Figure out the sector number
  Double_t phi    = TMath::ATan2(v.Y(), v.X()) * 180. / TMath::Pi();
  if (phi < 0) phi = 360. + phi;
  Double_t t      = phi - 2 * moduleno * theta;
  sector          = 2 * moduleno;
  if (t < 0 || t > 2 * theta) return kFALSE;
  else if (t > theta)         sector += 1;

  AliDebug(40, Form("<1> Inside an active FMD volume FMD%d%c[%2d,%3d] %s",
		    detector, ring, sector, strip, mc->CurrentVolPath()));
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDSimulator::VMC2FMD(Int_t copy, TLorentzVector& v,
                         UShort_t& detector, Char_t& ring,
                         UShort_t& sector, UShort_t& strip)
{
  TVirtualMC* mc = TVirtualMC::GetMC();

  strip = copy - 1;
  Int_t sectordiv; mc->CurrentVolOffID(fSectorOff, sectordiv);
  if (fModuleOff >= 0) {
    Int_t module;    mc->CurrentVolOffID(fModuleOff, module);
    sector = 2 * module + sectordiv;
  }
  else 
    sector = sectordiv;
  Int_t iring;     mc->CurrentVolOffID(fRingOff, iring); ring = Char_t(iring);
  Int_t det;       mc->CurrentVolOffID(fDetectorOff, det); detector = det;

  AliFMDGeometry* fmd = AliFMDGeometry::Instance();
  Double_t  rz  = fmd->GetDetector(detector)->GetRingZ(ring);
  Int_t     n   = fmd->GetDetector(detector)->GetRing(ring)->GetNSectors();
#if 0
  if (rz < 0) {
    Int_t s = ((n - sector + n / 2) % n) + 1;
    AliDebug(1, Form("Recalculating sector to %d (=%d-%d+%d/2%%%d+1 z=%f)",
		     s, n, sector, n, n, rz));
    sector = s;
  }
#endif
  if (sector < 1 || sector > n) {
    Warning("Step", "sector # %d out of range (0-%d)", sector-1, n-1);
    return kFALSE;
  }
  sector--;
  // Get track position
  mc->TrackPosition(v);
  AliDebug(40, Form("<2> Inside an active FMD volume FMD%d%c[%2d,%3d] %s",
		    detector, ring, sector, strip, mc->CurrentVolPath()));

  return kTRUE;
}

//____________________________________________________________________
void
AliFMDSimulator::Exec(Option_t* /* option */) 
{
  // Member function that is executed each time a hit is made in the
  // FMD.  None-charged particles are ignored.   Dead tracks  are
  // ignored. 
  //
  // The procedure is as follows: 
  // 
  //   - IF NOT track is alive THEN RETURN ENDIF
  //   - IF NOT particle is charged THEN RETURN ENDIF
  //   - IF NOT volume name is "STRI" or "STRO" THEN RETURN ENDIF 
  //   - Get strip number (volume copy # minus 1)
  //   - Get phi division number (mother volume copy #)
  //   - Get module number (grand-mother volume copy #)
  //   - section # = 2 * module # + phi division # - 1
  //   - Get ring Id from volume name 
  //   - Get detector # from grand-grand-grand-mother volume name 
  //   - Get pointer to sub-detector object. 
  //   - Get track position 
  //   - IF track is entering volume AND track is inside real shape THEN
  //   -   Reset energy deposited 
  //   -   Get track momentum 
  //   -   Get particle ID # 
  ///  - ENDIF
  //   - IF track is inside volume AND inside real shape THEN 
  ///  -   Update energy deposited 
  //   - ENDIF 
  //   - IF track is inside real shape AND (track is leaving volume,
  //         or it died, or it is stopped  THEN
  //   -   Create a hit 
  //   - ENDIF
  //     
  TVirtualMC* mc = TVirtualMC::GetMC();
  if (!mc->IsTrackAlive()) return;
  Double_t absQ = TMath::Abs(mc->TrackCharge());
  if (absQ <= 0) return;
  
  Int_t copy;
  Int_t vol = mc->CurrentVolID(copy);
  if (!IsActive(vol)) {
    AliDebug(50, Form("Not an FMD volume %d '%s' (%d or %d)", 
		      vol, mc->CurrentVolName(), fInnerId, fOuterId));
    return;
  }
  TLorentzVector v;
  UShort_t       detector;
  Char_t         ring;
  UShort_t       sector;
  UShort_t       strip;
  
  if (fUseDivided) {
    if (!VMC2FMD(copy, v, detector, ring, sector, strip)) return;
  } else {
    if (!VMC2FMD(v, detector, ring, sector, strip)) return;
  }
  TLorentzVector p;
  mc->TrackMomentum(p);
  Int_t    trackno = gAlice->GetMCApp()->GetCurrentTrackNumber();
  Int_t    pdg     = mc->TrackPid();
  Double_t mass    = mc->TrackMass();
  Double_t edep    = mc->Edep() * 1000; // keV
  Double_t poverm  = (mass == 0 ? 0 : p.P() / mass);

  // This `if' is to debug abnormal energy depositions.  We trigger on
  // p/m approx larger than or equal to a MIP, and a large edep - more 
  // than 1 keV - a MIP is 100 eV. 
  if (mc->Edep() * 1000 > absQ * absQ && poverm > 1) {
    TArrayI procs;
    mc->StepProcesses(procs);
    TString processes;
    for (Int_t ip = 0; ip < procs.fN; ip++) {
      if (ip != 0) processes.Append(",");
      processes.Append(TMCProcessName[procs.fArray[ip]]);
    }
    TDatabasePDG* pdgDB        = TDatabasePDG::Instance();
    TParticlePDG* particleType = pdgDB->GetParticle(pdg);
    TString pname(particleType ? particleType->GetName() : "???");
    TString what;
    if (mc->IsTrackEntering())    what.Append("entering ");
    if (mc->IsTrackExiting())     what.Append("exiting ");
    if (mc->IsTrackInside())      what.Append("inside ");
    if (mc->IsTrackDisappeared()) what.Append("disappeared ");
    if (mc->IsTrackStop())        what.Append("stopped ");
    if (mc->IsNewTrack())         what.Append("new ");
    if (mc->IsTrackAlive())       what.Append("alive ");
    if (mc->IsTrackOut())         what.Append("out ");
    
    Int_t mother = gAlice->GetMCApp()->GetPrimary(trackno);
    Warning("Step", "Track # %5d deposits a lot of energy\n" 
	    "  Volume:    %s\n" 
	    "  Momentum:  (%8.4f,%8.4f,%8.4f)\n"
	    "  PDG:       %d (%s)\n" 
	    "  Edep:      %-16.8f keV (mother %d)\n"
	    "  p/m:       %-16.8f\n"
	    "  Processes: %s\n"
	    "  What:      %s\n",
	    trackno, mc->CurrentVolPath(), p.X(), p.Y(), p.Z(),
	    pdg, pname.Data(), edep, mother, poverm, processes.Data(), 
	    what.Data());
  }
  
  // Check that the track is actually within the active area 
  Bool_t entering = mc->IsTrackEntering();
  Bool_t inside   = mc->IsTrackInside();
  Bool_t out      = (mc->IsTrackExiting()|| mc->IsTrackDisappeared()||
		     mc->IsTrackStop());
  // Reset the energy deposition for this track, and update some of
  // our parameters.
  if (entering) {
    AliDebug(15, Form("Track # %8d entering active FMD volume %s: "
		      "Edep=%f", 
		      gAlice->GetMCApp()->GetCurrentTrackNumber(),
		      mc->CurrentVolPath(), 1000 * mc->Edep()));
    fCurrentP      = p;
    fCurrentV      = v;    
    fCurrentDeltaE = edep;
    fCurrentPdg    = mc->IdFromPDG(pdg);
  }
  // If the track is inside, then update the energy deposition
  if (inside && fCurrentDeltaE >= 0) {
    fCurrentDeltaE += edep;
    AliDebug(15, Form("Track # %8d inside active FMD volume %s: Edep=%f, "
		      "Accumulated Edep=%f", 
		      trackno, mc->CurrentVolPath(), edep, 
		      fCurrentDeltaE));
  }
  // The track exits the volume, or it disappeared in the volume, or
  // the track is stopped because it no longer fulfills the cuts
  // defined, then we create a hit. 
  if (out && fCurrentDeltaE >= 0) {
    fCurrentDeltaE += edep;
    fFMD->AddHitByFields(trackno, detector, ring, sector, strip,
			 fCurrentV.X(),  fCurrentV.Y(), fCurrentV.Z(),
			 fCurrentP.X(),  fCurrentP.Y(), fCurrentP.Z(), 
			 fCurrentDeltaE, fCurrentPdg,   fCurrentV.T());
    fCurrentDeltaE = -1;
  }
}



//____________________________________________________________________
//
// EOF
//
