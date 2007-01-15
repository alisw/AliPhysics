/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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
/** @file    AliFMDHit.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:41:58 2006
    @brief   Hit in the FMD
    @ingroup FMD_sim
*/
//____________________________________________________________________
//
//  Hits in the FMD 
//  Contains information on:
//	Position of hit
//	Momentum of track
//	PID of track
//	Energy loss of track
//	Track #
//	Track path length
//	Track stopping status. 
//  Latest changes by Christian Holm Christensen
//

#include "Riostream.h"		// ROOT_Riostream
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TString.h>

#include "AliFMDHit.h"		// ALIFMDHIT_H
// #include "AliLog.h"		// ALILOG_H

//____________________________________________________________________
ClassImp(AliFMDHit)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif


//____________________________________________________________________
AliFMDHit::AliFMDHit()
  : fDetector(0), 
    fRing(0), 
    fSector(0), 
    fStrip('\0'), 
    fPx(0),
    fPy(0),
    fPz(0),
    fPdg(0),
    fEdep(0), 
    fTime(0), 
    fLength(0), 
    fStop(0)
{
  // Default CTOR
  fX = fY = fZ = 0;
}
  

//____________________________________________________________________
AliFMDHit::AliFMDHit(Int_t    shunt, 
		     Int_t    track, 
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
  : AliHit(shunt, track),
    fDetector(detector), 
    fRing(ring), 
    fSector(sector), 
    fStrip(strip), 
    fPx(px),
    fPy(py),
    fPz(pz),
    fPdg(pdg),
    fEdep(edep), 
    fTime(t), 
    fLength(l), 
    fStop(stop)
{
  // Normal FMD hit ctor
  // 
  // Parameters:
  // 
  //    shunt     ???
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
  fX = x;
  fY = y;
  fZ = z;
}

//____________________________________________________________________
const char*
AliFMDHit::GetName() const 
{ 
  // Get the name 
  static TString n;
  n = Form("FMD%d%c[%2d,%3d]", fDetector,fRing,fSector,fStrip);
  return n.Data();
}

//____________________________________________________________________
const char*
AliFMDHit::GetTitle() const 
{ 
  // Get the title 
  static TString t;
  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  TParticlePDG* pdg   = pdgDB->GetParticle(fPdg);
  t = Form("%s (%d): %f MeV / %f cm", (pdg ? pdg->GetName() : "?"), 
	   fTrack, fEdep, fLength);
  return t.Data();
}

//____________________________________________________________________
Float_t
AliFMDHit::P() const 
{
  // Get the momentum of the particle of the particle that made this
  // hit. 
  return TMath::Sqrt(fPx * fPx + fPy * fPy + fPz * fPz);
}

//____________________________________________________________________
Float_t
AliFMDHit::M() const 
{
  // Get the mass of the particle that made this hit. 
  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  TParticlePDG* pdg   = pdgDB->GetParticle(fPdg);
  return (pdg ? pdg->Mass() : -1);
}

//____________________________________________________________________
Float_t
AliFMDHit::Q() const
{
  // Get the charge of the particle that made this hit. 
  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  TParticlePDG* pdg   = pdgDB->GetParticle(fPdg);
  return (pdg ? pdg->Charge() : 0);
}


//____________________________________________________________________
void
AliFMDHit::Print(Option_t* option) const 
{
  // Print Hit to standard out 
  cout << "AliFMDHit: FMD" 
       << fDetector << fRing << "[" 
       << setw(3) << fSector << ","
       << setw(3) << fStrip << "] = " 
       << fEdep << endl;
  TString opt(option);
  if (opt.Contains("D", TString::kIgnoreCase)) {
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    TParticlePDG* pdg   = pdgDB->GetParticle(fPdg);
    cout << "\tPDG:\t" << fPdg << " " << (pdg ? pdg->GetName() : "?") << "\n"
	 << "\tP:\t(" << fPx << "," << fPy << "," << fPz << ") "<<P() << "\n" 
	 << "\tX:\t" << fX << "," << fY << "," << fZ << "\n" 
	 << "\tTrack #:\t" << fTrack << "\tLength:\t" 
	 << fLength << "cm\t" << (IsStop() ? "stopped" : "") << std::endl;
  }
}

//____________________________________________________________________
//
// EOF
//
