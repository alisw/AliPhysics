// **************************************************************************
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: The ALICE Off-line Project.                                    *
// * Contributors are mentioned in the code where appropriate.              *
// *                                                                        *
// * Permission to use, copy, modify and distribute this software and its   *
// * documentation strictly for non-commercial purposes is hereby granted   *
// * without fee, provided that the above copyright notice appears in all   *
// * copies and that both the copyright notice and this permission notice   *
// * appear in the supporting documentation. The authors make no claims     *
// * about the suitability of this software for any purpose. It is          *
// * provided "as is" without express or implied warranty.                  *
// **************************************************************************


#include "AliRICHv1.h"
#include "AliRICHParam.h"
#include "AliRICHChamber.h"
#include <TParticle.h> 
#include <TRandom.h> 
#include <TVirtualMC.h>
#include <TPDGCode.h>

#include <AliConst.h>
#include <AliPDG.h>
#include <AliRun.h>
#include <AliMC.h>

ClassImp(AliRICHv1)    
//__________________________________________________________________________________________________
void AliRICHv1::StepManager()
{
//Full Step Manager

  Int_t          copy;
  static Int_t   iCurrentChamber;
        
//Treat photons    
  static TLorentzVector cerX4;
  if((gMC->TrackPid()==kCerenkov||gMC->TrackPid()==kFeedback)&&gMC->CurrentVolID(copy)==gMC->VolId("CSI ")){//photon in CSI
    if(gMC->Edep()>0.){//CF+CSI+DE
      if(IsLostByFresnel()){ gMC->StopTrack(); return;}        
      gMC->TrackPosition(cerX4); gMC->CurrentVolOffID(2,iCurrentChamber);
	
      AddHit(iCurrentChamber,gAlice->GetMCApp()->GetCurrentTrackNumber(),cerX4.Vect(),cerX4.Vect());//HIT for PHOTON in conditions CF+CSI+DE
      GenerateFeedbacks(iCurrentChamber);
    }//CF+CSI+DE
  }//CF in CSI
  
//Treat charged particles  
  static Float_t eloss;
  static TLorentzVector mipInX4,mipOutX4;
  if(gMC->TrackCharge() && gMC->CurrentVolID(copy)==gMC->VolId("GAP ")){//MIP in GAP
    gMC->CurrentVolOffID(3,iCurrentChamber);
    if(gMC->IsTrackEntering()||gMC->IsNewTrack()) {//MIP in GAP entering or newly created
      eloss=0;                                                           
      gMC->TrackPosition(mipInX4);
    }else if(gMC->IsTrackExiting()||gMC->IsTrackStop()||gMC->IsTrackDisappeared()){//MIP in GAP exiting or disappeared
      eloss+=gMC->Edep();//take into account last step dEdX
      gMC->TrackPosition(mipOutX4);  
      AddHit(iCurrentChamber,gAlice->GetMCApp()->GetCurrentTrackNumber(),mipInX4.Vect(),mipOutX4.Vect(),eloss);//HIT for MIP: MIP in GAP Exiting
      GenerateFeedbacks(iCurrentChamber,eloss);//MIP+GAP+Exit
    }else//MIP in GAP going inside
      eloss   += gMC->Edep();
  }//MIP in GAP
}//StepManager()
//__________________________________________________________________________________________________
Bool_t AliRICHv1::IsLostByFresnel()
{
  TLorentzVector p4;
  Double_t mom[3],localMom[3];
  gMC->TrackMomentum(p4); mom[0]=p4(0);   mom[1]=p4(1);   mom[2]=p4(2);   mom[3]=p4(3);
  gMC->Gmtod(mom,localMom,2);
  Double_t localTc    = localMom[0]*localMom[0]+localMom[2]*localMom[2];
  Double_t localTheta = TMath::ATan2(TMath::Sqrt(localTc),localMom[1]);
  Double_t cotheta = TMath::Abs(TMath::Cos(localTheta));
  if(gMC->GetRandom()->Rndm() < Fresnel(p4.E()*1e9,cotheta,1)){
    if(GetDebug()) Info("IsLostByFresnel","");
    return kTRUE;
  }else
    return kFALSE;
}//IsLostByFresnel()
//__________________________________________________________________________________________________
