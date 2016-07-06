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
/** @file    AliFMDv1.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:48:51 2006
    @brief   Concrete implementation of FMD detector driver - detailed
    version 
    @ingroup FMD_sim
*/
//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. This class
// contains the base procedures for the Forward Multiplicity detector
// Detector consists of 3 sub-detectors FMD1, FMD2, and FMD3, each of
// which has 1 or 2 rings of silicon sensors. 
// This class contains the detailed version of the FMD - that is, hits
// are produced during simulation. 
//                                                                           
// See also the class AliFMD for a more detailed explanation of the
// various componets. 
//
#include <TVirtualMC.h>		// ROOT_TVirtualMC
#include <AliRun.h>		// ALIRUN_H
#include <AliMC.h>		// ALIMC_H
// #include <AliLog.h>		// ALILOG_H
#include "AliFMDDebug.h" // Better debug macros
#include "AliFMDv1.h"		// ALIFMDV1_H
// #include "AliFMDGeometryBuilder.h"
#include "AliFMDGeometry.h"
#include "AliFMDDetector.h"
#include "AliFMDRing.h"
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include "AliFMDHit.h"

//____________________________________________________________________
ClassImp(AliFMDv1)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif


//____________________________________________________________________
Bool_t
AliFMDv1::VMC2FMD(TLorentzVector& v, UShort_t& detector,
		  Char_t& ring, UShort_t& sector, UShort_t& strip) const
{
  // Convert VMC coordinates to detector coordinates 
  TVirtualMC*      mc  = TVirtualMC::GetMC();
  AliFMDGeometry*  fmd = AliFMDGeometry::Instance();

  // Get track position
  mc->TrackPosition(v);
  Int_t moduleno; mc->CurrentVolOffID(fmd->GetModuleOff(), moduleno);
  Int_t iring;    mc->CurrentVolOffID(fmd->GetRingOff(), iring);   
  ring = Char_t(iring);
  Int_t det;      mc->CurrentVolOffID(fmd->GetDetectorOff(), det); 
  detector = det;
  

  // Get the ring geometry
  //Int_t     nsec = fmd->GetDetector(detector)->GetRing(ring)->GetNSectors();
  Int_t     nstr  = fmd->GetDetector(detector)->GetRing(ring)->GetNStrips();
  Double_t  lowr  = fmd->GetDetector(detector)->GetRing(ring)->GetMinR();
  Double_t  theta = fmd->GetDetector(detector)->GetRing(ring)->GetTheta();
  Double_t  pitch = fmd->GetDetector(detector)->GetRing(ring)->GetPitch();

  // Figure out the strip number
  Double_t r     = TMath::Sqrt(v.X() * v.X() + v.Y() * v.Y());
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

  AliFMDDebug(40, ("<1> Inside an active FMD volume FMD%d%c[%2d,%3d] %s",
		    detector, ring, sector, strip, mc->CurrentVolPath()));
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDv1::VMC2FMD(Int_t copy, TLorentzVector& v,
		  UShort_t& detector, Char_t& ring,
		  UShort_t& sector, UShort_t& strip) const
{
  // Convert VMC coordinates to detector coordinates 
  TVirtualMC*     mc  = TVirtualMC::GetMC();
  AliFMDGeometry* fmd = AliFMDGeometry::Instance();

  strip = copy - 1;
  Int_t sectordiv; mc->CurrentVolOffID(fmd->GetSectorOff(), sectordiv);
  if (fmd->GetModuleOff() >= 0) {
    Int_t module;    mc->CurrentVolOffID(fmd->GetModuleOff(), module);
    sector = 2 * module + sectordiv;
  }
  else 
    sector = sectordiv;
  AliFMDDebug(30, ("Getting ring volume with offset %d -> %s", 
		    fmd->GetRingOff(), 
		    mc->CurrentVolOffName(fmd->GetRingOff())));
  Int_t iring;     mc->CurrentVolOffID(fmd->GetRingOff(), iring); 
  ring = Char_t(iring);
  Int_t det;       mc->CurrentVolOffID(fmd->GetDetectorOff(), det); 
  detector = det;

  //Double_t  rz  = fmd->GetDetector(detector)->GetRingZ(ring);
  AliFMDDetector* gdet  = fmd->GetDetector(detector);
  AliFMDRing*     gring = gdet->GetRing(ring);
  if (!gring) {
    AliFatal(Form("Ring %c not found (volume was %s at offset %d in path %s)", 
		  ring, 
		  mc->CurrentVolOffName(fmd->GetRingOff()),
		  fmd->GetRingOff(), 
		  mc->CurrentVolPath()));
  }
  Int_t n = gring->GetNSectors();
#if 0
  if (rz < 0) {
    Int_t s = ((n - sector + n / 2) % n) + 1;
    AliFMDDebug(1, ("Recalculating sector to %d (=%d-%d+%d/2%%%d+1 z=%f)",
		     s, n, sector, n, n, rz));
    sector = s;
  }
#endif
  if (sector < 1 || sector > n) {
    AliWarning(Form("sector # %d out of range (0-%d)", sector-1, n-1));
    return kFALSE;
  }
  sector--;
  // Get track position
  mc->TrackPosition(v);
  AliFMDDebug(40, ("<2> Inside an active FMD volume FMD%d%c[%2d,%3d] %s",
		    detector, ring, sector, strip, mc->CurrentVolPath()));

  return kTRUE;
}


//____________________________________________________________________
Bool_t
AliFMDv1::CheckHit(Int_t trackno, Int_t pdg, Float_t absQ, 
		   const TLorentzVector& p, Float_t edep) const
{
  // Check that a hit is good 
  if (AliLog::GetDebugLevel("FMD", "AliFMD") < 5) return kFALSE;
  TVirtualMC* mc   = TVirtualMC::GetMC();
  Double_t mass    = mc->TrackMass();
  Double_t poverm  = (mass == 0 ? 0 : p.P() / mass);
    
  // This `if' is to debug abnormal energy depositions.  We trigger on
  // p/m approx larger than or equal to a MIP, and a large edep - more 
  // than 1 keV - a MIP is 100 eV. 
  if (!(edep > absQ * absQ && poverm > 1)) return kFALSE;
  
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
  AliFMDDebug(15, ("Track # %5d deposits a lot of energy\n" 
		    "  Volume:    %s\n" 
		    "  Momentum:  (%7.4f,%7.4f,%7.4f)\n"
		    "  PDG:       %d (%s)\n" 
		    "  Edep:      %-14.7f keV (mother %d)\n"
		    "  p/m:       %-7.4f/%-7.4f = %-14.7f\n"
		    "  Processes: %s\n"
		    "  What:      %s\n",
		    trackno, mc->CurrentVolPath(), p.X(), p.Y(), p.Z(),
		    pdg, pname.Data(), edep, mother, p.P(), mass, 
		    poverm, processes.Data(), what.Data()));
  return kTRUE;
}


//____________________________________________________________________
void 
AliFMDv1::StepManager()
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
  TVirtualMC* mc = fMC;
  if (!mc->IsTrackAlive()) return;
  Double_t absQ = TMath::Abs(mc->TrackCharge());
  if (absQ <= 0) return;
  
  Int_t copy;
  Int_t vol = mc->CurrentVolID(copy);
  AliFMDGeometry*  fmd = AliFMDGeometry::Instance();
  if (!fmd->IsActive(vol)) {
    AliFMDDebug(50, ("Not an FMD volume %d '%s'",vol,mc->CurrentVolName()));
    return;
  }
  TLorentzVector v;
  UShort_t       detector;
  Char_t         ring;
  UShort_t       sector;
  UShort_t       strip;
  
  if (fmd->IsDetailed()) {
    if (!VMC2FMD(copy, v, detector, ring, sector, strip)) return;
  } else {
    if (!VMC2FMD(v, detector, ring, sector, strip)) return;
  }
  TLorentzVector p;
  mc->TrackMomentum(p);
  Int_t    trackno = gAlice->GetMCApp()->GetCurrentTrackNumber();
  Int_t    pdg     = mc->TrackPid();
  Double_t edep    = mc->Edep() * 1000; // keV
  Bool_t   isBad   = CheckHit(trackno, pdg, absQ, p, edep);
  
  // Check that the track is actually within the active area 
  Bool_t entering = mc->IsTrackEntering();
  Bool_t inside   = mc->IsTrackInside();
  Bool_t out      = (mc->IsTrackExiting()|| mc->IsTrackDisappeared()||
		     mc->IsTrackStop());
  // Reset the energy deposition for this track, and update some of
  // our parameters.
  if (entering) {
    AliFMDDebug(15, ("Track # %8d entering active FMD volume %s: "
		      "Edep=%f (%f,%f,%f)", trackno, mc->CurrentVolPath(),
		      edep, v.X(), v.Y(), v.Z()));
    fCurrentP      = p;
    fCurrentV      = v;    
    fCurrentDeltaE = edep;
    fCurrentPdg    = pdg; // mc->IdFromPDG(pdg);
  }
  // If the track is inside, then update the energy deposition
  if (inside && fCurrentDeltaE >= 0) {
    fCurrentDeltaE += edep;
    AliFMDDebug(15, ("Track # %8d inside active FMD volume %s: Edep=%f, "
		      "Accumulated Edep=%f  (%f,%f,%f)", trackno, 
		      mc->CurrentVolPath(), edep, fCurrentDeltaE, 
		      v.X(), v.Y(), v.Z()));
  }
  // The track exits the volume, or it disappeared in the volume, or
  // the track is stopped because it no longer fulfills the cuts
  // defined, then we create a hit. 
  if (out) {
    if (fCurrentDeltaE >= 0) {
      fCurrentDeltaE += edep;
      AliFMDDebug(15, ("Track # %8d exiting active FMD volume %s: Edep=%g, "
			"Accumulated Edep=%g (%f,%f,%f)", trackno, 
			mc->CurrentVolPath(), edep, fCurrentDeltaE, 
			v.X(), v.Y(), v.Z()));
      TVector3 cur(v.Vect());
      cur -= fCurrentV.Vect();
      Double_t len = cur.Mag();
      AliFMDHit* h = 
	AddHitByFields(trackno, detector, ring, sector, strip,
		       fCurrentV.X(),  fCurrentV.Y(), fCurrentV.Z(),
		       fCurrentP.X(),  fCurrentP.Y(), fCurrentP.Z(), 
		       fCurrentDeltaE, fCurrentPdg,   fCurrentV.T(), 
		       len, mc->IsTrackDisappeared()||mc->IsTrackStop());
      // Add a copy 
      if (isBad && fBad) { 
	new ((*fBad)[fBad->GetEntries()]) AliFMDHit(*h);
      }
      // Check the geometry that we can get back the coordinates. 
#ifdef CHECK_TRANS
      Double_t x, y, z;
      fmd->Detector2XYZ(detector, ring, sector, strip, x, y ,z);
      AliFMDDebug(1, ("Hit at (%f,%f,%f), geometry says (%f,%f,%f)", 
		       fCurrentV.X(), fCurrentV.Y(), fCurrentV.Z(), x, y, z));
#endif
    }
    fCurrentDeltaE = -1;
  }
}
//___________________________________________________________________
//
// EOF
//
