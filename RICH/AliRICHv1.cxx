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
//______________________________________________________________________________
void AliRICHv1::StepManager()
{
//Full Step Manager

  Int_t          copy;
  static Int_t   iCurrentChamber;
  static TLorentzVector x4,p4,mipInX4,mipOutX4;
  Float_t        pos[3],mom[4],localPos[3],localMom[4];
  Float_t        coscerenkov;
       
  TParticle *current = (TParticle*)(*gAlice->GetMCApp()->Particles())[gAlice->GetMCApp()->GetCurrentTrackNumber()];
 
  Float_t cherenkovLoss=0;
    
  if(gMC->TrackPid()==kCerenkov){//C
    Float_t ckovEnergy = current->Energy();
    if(ckovEnergy > 5.6e-09 && ckovEnergy < 7.8e-09 ){//C+E
        if(gMC->IsTrackEntering()){                                     //is track entering?
		    
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy)){             //is it in csi?      
			gMC->TrackMomentum(p4);	mom[0]=p4(0);	mom[1]=p4(1);	mom[2]=p4(2);	mom[3]=p4(3);
			gMC->Gmtod(mom,localMom,2);
			Double_t localTc    = localMom[0]*localMom[0]+localMom[2]*localMom[2];
			Double_t localTheta = TMath::ATan2(TMath::Sqrt(localTc),localMom[1]);
			Double_t cotheta = TMath::Abs(TMath::Cos(localTheta));
                        if(gMC->GetRandom()->Rndm() < Fresnel(p4.E()*1e9,cotheta,1))  gMC->StopTrack();
				
		      }//C+E+produced in Freon
		  } //track entering?
    }//C+E
  }//C
    
  if((gMC->TrackPid()==kCerenkov||gMC->TrackPid()==kFeedback)&&gMC->CurrentVolID(copy)==gMC->VolId("CSI ")){//photon in CSI     
      if(gMC->Edep()>0.){//CF+CSI+DE
        gMC->TrackPosition(x4);   pos[0]=x4(0);   pos[1]=x4(1);   pos[2]=x4(2);
        gMC->TrackMomentum(p4);   mom[0]=p4(0);   mom[1]=p4(1);   mom[2]=p4(2);   mom[3]=p4(3);
        if(IsFresnelLoss()){ gMC->StopTrack(); return;}        
	gMC->CurrentVolOffID(2,copy);iCurrentChamber=copy;
	
        gMC->Gmtod(pos,localPos,1);     gMC->Gmtod(mom,localMom,2);

	cherenkovLoss  += gMC->Edep();
		    
        AliRICHhit *mipHit =  (AliRICHhit*) (fHits->UncheckedAt(0));
        if(mipHit){
          mom[0] = current->Px();   mom[1] = current->Py();   mom[2] = current->Pz();
          Float_t mipPx = mipHit->MomX();   Float_t mipPy = mipHit->MomY();   Float_t mipPz = mipHit->MomZ();
			
          Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
          Float_t rt = TMath::Sqrt(r);
          Float_t mipR = mipPx*mipPx + mipPy*mipPy + mipPz*mipPz;	
          Float_t mipRt = TMath::Sqrt(mipR);
          if((rt*mipRt) > 0)
            coscerenkov = (mom[0]*mipPx + mom[1]*mipPy + mom[2]*mipPz)/(rt*mipRt);
          else
            coscerenkov = 0;
        }
        AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),x4.Vect());//HIT for PHOTON in conditions CF+CSI+DE
	GenerateFeedbacks(iCurrentChamber,cherenkovLoss);//CF+CSI+DE
      }//CF+CSI+DE
  }//CF in CSI
  
//Treat charged particles  
  static Float_t eloss;
  if(gMC->TrackCharge() && gMC->CurrentVolID(copy)==gMC->VolId("GAP ")){//MIP in GAP
    gMC->CurrentVolOffID(3,copy); iCurrentChamber=copy;
    if(gMC->IsTrackEntering()||gMC->IsNewTrack()) {//MIP in GAP Entering
      eloss=0;                                                           
      gMC->TrackPosition(mipInX4);
    }else if(gMC->IsTrackExiting()||gMC->IsTrackStop()||gMC->IsTrackDisappeared()){//MIP in GAP Exiting
      eloss+=gMC->Edep();//take into account last step dEdX
      gMC->TrackPosition(mipOutX4);  
      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),mipInX4.Vect(),mipOutX4.Vect(),eloss);//HIT for MIP for conditions: MIP in GAP Exiting
      GenerateFeedbacks(iCurrentChamber,eloss);//MIP+GAP+Exit
    }else//MIP in GAP going inside
      eloss   += gMC->Edep();
  }//MIP in GAP
}//void AliRICHv1::StepManager()

Bool_t AliRICHv1::IsFresnelLoss()
{
  return kFALSE;
}
