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

#include "AliRICHv0.h"
#include <TVirtualMC.h>
#include <TPDGCode.h>
#include "AliRICHConst.h" 
#include <AliRun.h>
#include <TLorentzVector.h>

ClassImp(AliRICHv0)
//__________________________________________________________________________________________________
AliRICHv0::AliRICHv0(const char *name, const char *title)
          :AliRICH(name,title)
{
  if(GetDebug())Info("named ctor","Start.");
  if(GetDebug())Info("named ctor","Stop.");
}//name ctor
//__________________________________________________________________________________________________
void AliRICHv0::StepManager()
{//
//  if(!gMC->IsNewTrack()) return;
 
  char *sParticle;
  switch(gMC->TrackPid()){
    case kProton:
      sParticle="p";break;
    case kNeutron:
      sParticle="n";break;
    case kGamma:
      sParticle="gamma";break;
    case 50000050:
      sParticle="photon";break;
    default:
      sParticle="not known";break;
  }

  Info("StepManager","Event=%i hunt=%i TID=%i PID=%s Mass=%f Charge=%i",
                      gMC->CurrentEvent(),
                                fIshunt,
                                        gAlice->GetCurrentTrackNumber(),
                                               sParticle,
                                                      gMC->TrackMass(),
                                                              gMC->TrackCharge());
  Info("StepManager","Flags:Alive(%i) Disap(%i) Enter(%i) Exit(%i) Inside(%i) Out(%i) Stop(%i) New(%i)",
                            gMC->IsTrackAlive(),
                                      gMC->IsTrackDisappeared(),
                                                gMC->IsTrackEntering(),
                                                          gMC->IsTrackExiting(),
                                                                    gMC->IsTrackInside(),
                                                                          gMC->IsTrackOut(),
                                                                                  gMC->IsTrackStop(),
                                                                                       gMC->IsNewTrack());
  Info("StepManager","Volume=%s of volume=%s",
                      gMC->CurrentVolName(),gMC->CurrentVolOffName(1));

//  Info("StepManager","TrackPID %i Particle %i",
//                     gMC->TrackPid(),gAlice->Particles()[gAlice->CurrentTrack()]
  TLorentzVector x4;
  gMC->TrackPosition(x4);
  Info("StepManager","x=%f y=%f z=%f r=%f theta=%f phi=%f\n",
                      x4.X(),x4.Y(),x4.Z(),x4.Rho(),x4.Theta()*r2d,x4.Phi()*r2d);  
}//AliRICHv0::StepManager()
//__________________________________________________________________________________________________
