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
// This class contains the detailed version of the FMD - that is, hits
// are produced during simulation. 
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
// See also the class AliFMD for a more detailed explanation of the
// various componets. 
#include "TVirtualMC.h"		// ROOT_TVirtualMC
#include "AliFMDv1.h"		// ALIFMDV1_H
#include "AliRun.h"		// ALIRUN_H
#include "AliMC.h"		// ALIMC_H
#include "AliLog.h"		// ALILOG_H

//____________________________________________________________________
ClassImp(AliFMDv1);


//____________________________________________________________________
void 
AliFMDv1::StepManager()
{
  //
  // Called for every step in the Forward Multiplicity Detector
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
  //
  // DebugGuard guard("AliFMDv1::StepManager");
  AliDebug(10, "AliFMDv1::StepManager");
  // return;

  // If the track is gone, return
  if (!gMC->IsTrackAlive()) return;
  
  // Only process charged particles 
  if(TMath::Abs(gMC->TrackCharge()) <= 0) return; 

  // TString vol(gMC->CurrentVolName());
  // std::cout << "Is inside " << vol << " ... " << std::endl;
  // Only do stuff is the track is in one of the strips. 
  // TString vol(gMC->CurrentVolName());
  // if (!vol.Contains("STR")) return;
  Int_t copy;
  Int_t volumeId = gMC->CurrentVolID(copy);
  // The ring ID is encoded in the volume name 
  Char_t ring = '\0';
  if (volumeId == fInner->GetStripId())      ring = 'I';
  else if (volumeId == fOuter->GetStripId()) ring = 'O'; 
  else                                       return;

  // Get the strip number.  Note, that GEANT numbers divisions from 1,
  // so we subtract one 
  Int_t strip = copy - 1;

  // Get the phi division of the module 
  Int_t phiDiv;                         // * The phi division number (1 or 2)
  gMC->CurrentVolOffID(1, phiDiv);      //   in the module  

  // Active volume number - not used. 
  // Int_t active;                         
  // gMC->CurrentVolOffID(2, active);      

  // Get the module number in the ring. 
  Int_t module;                    
  gMC->CurrentVolOffID(3, module); 
  
  // Ring copy number - the same as the detector number - not used
  // Int_t ringCopy;                       // * Ring copy number
  // gMC->CurrentVolOffID(4, ringCopy);    //   Same as detector number 
  
  // Get the detector number from the path name 
  Int_t detector = Int_t((gMC->CurrentVolOffName(5)[3]) - 48);

  // The sector number, calculated from module and phi division # 
  Int_t  sector =  2 * module + phiDiv - 1;

  
  // Get a pointer to the sub detector structure 
  AliFMDSubDetector* det = 0;
  switch (detector) {
  case 1: det = fFMD1; break;
  case 2: det = fFMD2; break;
  case 3: det = fFMD3; break;
  }
  if (!det) return;

  // Get the current track position 
  TLorentzVector v;
  gMC->TrackPosition(v);
  // Check that the track is actually within the active area 
  Bool_t isWithin = det->CheckHit(ring, module, v.X(), v.Y());
  Bool_t entering = gMC->IsTrackEntering() && isWithin;
  Bool_t inside   = gMC->IsTrackInside()   && isWithin;
  Bool_t out      = (gMC->IsTrackExiting() 
		     || gMC->IsTrackDisappeared() 
		     || gMC->IsTrackStop() 
		     || !isWithin);
// Reset the energy deposition for this track, and update some of
  // our parameters.
  if (entering) {
    fCurrentDeltaE = 0;

    // Get production vertex and momentum of the track 
    fCurrentV = v;
    gMC->TrackMomentum(fCurrentP);
    fCurrentPdg = gMC->IdFromPDG(gMC->TrackPid());

    // if (fAnalyser) 
    //   fAnalyser->Update(detector, ring, isWithin, v.X(), v.Y());
  }
  
  // If the track is inside, then update the energy deposition
  if (inside && fCurrentDeltaE >= 0) 
    fCurrentDeltaE += 1000 * gMC->Edep();

  // The track exits the volume, or it disappeared in the volume, or
  // the track is stopped because it no longer fulfills the cuts
  // defined, then we create a hit. 
  if (out && fCurrentDeltaE >= 0) {
    fCurrentDeltaE += 1000 * gMC->Edep();

    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),
	   detector, ring,  sector, strip,
	   fCurrentV.X(), fCurrentV.Y(), fCurrentV.Z(),
	   fCurrentP.X(), fCurrentP.Y(), fCurrentP.Z(), 
	   fCurrentDeltaE, fCurrentPdg, fCurrentV.T());
    fCurrentDeltaE = -1;
  }
}
//___________________________________________________________________
//
// EOF
//
