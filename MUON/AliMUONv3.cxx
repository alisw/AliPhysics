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
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// AliMUONv3 Addapted for AliMUONv1 
// This was the last revision of AliMUONv1
// $Log$
// Revision 1.40  2003/01/28 13:21:06  morsch
// Improved response simulation for station 1.
// (M. Mac Cormick, I. Hrivnacova, D. Guez)
//
// Gines Martinez (Subatech) jan 2003

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 2    //
/////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TLorentzVector.h> 
#include <TNode.h> 
#include <TRandom.h> 
#include <TTUBE.h>

#include "AliMUONv3.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliCallf77.h"
#include "AliConst.h" 
#include "AliMUONChamber.h"
#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliMUONConstants.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONFactory.h"

ClassImp(AliMUONv3)
 
//___________________________________________
AliMUONv3::AliMUONv3() : AliMUONv1()
{
// Constructor
    fChambers = 0;
    fStations = 0;
}
 
//___________________________________________
AliMUONv3::AliMUONv3(const char *name, const char *title)
       : AliMUONv1(name,title)
{
// Constructor
    // By default include all stations
    fStations = new Int_t[5];
    for (Int_t i=0; i<5; i++) fStations[i] = 1;

    AliMUONFactory factory;
    factory.Build(this, title);
}
//___________________________________________
void AliMUONv3::StepManager()
{

  // Volume id
  Int_t   copy, id;
  Int_t   idvol;
  Int_t   iChamber=0;
  // Particule id, pos and mom vectors, 
  // theta, phi angles with respect the normal of the chamber, 
  // spatial step, delta_energy and time of flight
  Int_t          ipart;
  TLorentzVector pos, mom;
  Float_t        theta, phi, tof;
  Float_t        destep, step;

  TClonesArray &lhits = *fHits;

  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 

  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  idvol=-1;
  id=gMC->CurrentVolID(copy);
  for (Int_t i = 1; i <= AliMUONConstants::NCh(); i++) {
    if(id==((AliMUONChamber*)(*fChambers)[i-1])->GetGid()) {
      iChamber = i;
      idvol  = i-1;
    }
  }
  if (idvol == -1) return;

  // Get current particle id (ipart), track position (pos)  and momentum (mom)
  gMC->TrackPosition(pos);
  gMC->TrackMomentum(mom);
  ipart    = gMC->TrackPid();
  theta    = mom.Theta()*kRaddeg;     // theta of track
  phi      = mom.Phi()  *kRaddeg;     // phi of the track
  tof      = gMC->TrackTime();        // Time of flight
  //
  // momentum loss and steplength in last step
  destep = gMC->Edep();
  step   = gMC->TrackStep();
  //    new hit       
  new(lhits[fNhits++]) 
    AliMUONHit(fIshunt, gAlice->CurrentTrack(), iChamber, ipart, pos.X(), pos.Y(), pos.Z(), tof, mom.P(), theta, phi, step, destep);
}


