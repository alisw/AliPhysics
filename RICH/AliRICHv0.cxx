//**************************************************************************
//  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//                                                                         *
//  Author: The ALICE Off-line Project.                                    *
//  Contributors are mentioned in the code where appropriate.              *
//                                                                         *
//  Permission to use, copy, modify and distribute this software and its   *
//  documentation strictly for non-commercial purposes is hereby granted   *
//  without fee, provided that the above copyright notice appears in all   *
//  copies and that both the copyright notice and this permission notice   *
//  appear in the supporting documentation. The authors make no claims     *
//  about the suitability of this software for any purpose. It is          *
//  provided "as is" without express or implied warranty.                  *
//**************************************************************************

#include "AliRICHv0.h"
#include "AliRICHConst.h"
#include "AliRICHChamber.h" 
#include <AliRun.h>
#include <AliMC.h>
#include <TVirtualMC.h>
#include <TPDGCode.h>
#include <TLorentzVector.h>
#include <TMath.h>

ClassImp(AliRICHv0)

void AliRICHv0::StepManager()
{
//This StepManager is a provision for different test-learn activities on the current MC layer  
  char *sParticle;
  switch(gMC->TrackPid()){
    case kProton:
      sParticle="proton";break;
    case kNeutron:
      sParticle="neutron";break;
    case kGamma:
      sParticle="gamma";break;
    case kCerenkov:
      sParticle="photon";break;
    case kPi0:
      sParticle="Pi0";break;  
    case kElectron:
      sParticle="electron";break;  
    default:
      sParticle="not known";break;
  }

  Info("","event=%i hunt=%i tid=%i pid=%i(%s) m=%f q=%3.1f",
                            gMC->CurrentEvent(),
                            fIshunt,
                            gAlice->GetMCApp()->GetCurrentTrackNumber(),
                            gMC->TrackPid(),
                            sParticle,
                            gMC->TrackMass(),
                            gMC->TrackCharge());
  Info("","Flags:alive(%i) disap(%i) enter(%i) exit(%i) inside(%i) out(%i) stop(%i) new(%i)",
                            gMC->IsTrackAlive(),
                            gMC->IsTrackDisappeared(),
                            gMC->IsTrackEntering(),
                            gMC->IsTrackExiting(),
                            gMC->IsTrackInside(),
                            gMC->IsTrackOut(),
                            gMC->IsTrackStop(),
                            gMC->IsNewTrack());
  Int_t copy0,copy1,copy2,copy3;
  Int_t vid0=gMC->CurrentVolID(copy0);
  Int_t vid1=gMC->CurrentVolOffID(1,copy1);
  Int_t vid2=gMC->CurrentVolOffID(2,copy2);
  Int_t vid3=gMC->CurrentVolOffID(3,copy3);
  Info("","vid0=%i(%s)c%i vid1=%i(%s)c%i vid2=%i(%s)c%i vid3=%i(%s)c%i   %s-%s-%s-%s",
                      vid0,gMC->VolName(vid0),copy0, 
                      vid1,gMC->VolName(vid1),copy1, 
                      vid2,gMC->VolName(vid2),copy2, 
                      vid3,gMC->VolName(vid3),copy3, 
                      gMC->CurrentVolName(),
                      gMC->CurrentVolOffName(1),
                      gMC->CurrentVolOffName(2),
                      gMC->CurrentVolOffName(3));
  
  Float_t a,z,den,rad,abs; a=z=den=rad=abs=kBad;
  Int_t mid=gMC->CurrentMaterial(a,z,den,rad,abs);
  Info("","mid=%i a=%7.2f z=%7.2f den=%7.2f rad=%7.2f abs=%7.2f",mid,a,z,den,rad,abs);
  
  TLorentzVector x4;
  gMC->TrackPosition(x4);
  Float_t glo[3],loc[3];
  glo[0]=x4.X();glo[1]=x4.Y();glo[2]=x4.Z();  
  gMC->Gmtod(glo,loc,1);
  Info("","glo(%+8.3f,%+8.3f,%+8.3f) r=%8.3f theta=%8.3f phi=%8.3f",
                      glo[0],glo[1],glo[2],x4.Rho(),x4.Theta()*TMath::RadToDeg(),x4.Phi()*TMath::RadToDeg());  
  Info("","loc(%+8.3f,%+8.3f,%8.3f) by gMC->Gmtod()",         loc[0],loc[1],loc[2]);  
  if(gMC->VolId("CSI ")==gMC->CurrentVolID(copy0)){
    Int_t iChamber;
    gMC->CurrentVolOffID(2,iChamber);
    TVector3 x3=C(iChamber)->G2L(x4);
    Info("","loc(%+8.3f,%+8.3f,%8.3f) by G2L",         x3.X(),x3.Y(),x3.Z());  
    x3=C(iChamber)->Glob2Loc(x4);
    Info("","loc(%+8.3f,%+8.3f,%8.3f) by Global2Local",         x3.X(),x3.Y(),x3.Z());  
  }
  Info("","end of current step\n");
}//StepManager()
//__________________________________________________________________________________________________
