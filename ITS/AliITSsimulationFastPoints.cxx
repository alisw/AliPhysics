/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                          *
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

/*
$Log$
*/
#include <TRandom.h>

#include "AliITS.h"
#include "AliRun.h"
#include "AliITSsimulationFastPoints.h"


ClassImp(AliITSsimulationFastPoints)

AliITSsimulationFastPoints::AliITSsimulationFastPoints()
{
  //constructor
  fSigmaRPhi[0] = fSigmaRPhi[1] = 12e-4;
  fSigmaRPhi[2] = fSigmaRPhi[3] = 38e-4;
  fSigmaRPhi[4] = fSigmaRPhi[5] = 20e-4;
  fSigmaZ[0] = fSigmaZ[1] = 70e-4;
  fSigmaZ[2] = fSigmaZ[3] = 28e-4;
  fSigmaZ[4] = fSigmaZ[5] = 830e-4;
  fSigmaDe[0] = fSigmaDe[1] = 0.72e-6;
  fSigmaDe[2] = fSigmaDe[3] = 0.90e-6;
  fSigmaDe[4] = fSigmaDe[5] =  5e-6;
  fThrDe[0] = fThrDe[1] = 7.2e-6;
  fThrDe[2] = fThrDe[3] = 2.70e-6;
  fThrDe[4] = fThrDe[5] = 10e-6;
}

//-------------------------------------------------------------
void AliITSsimulationFastPoints::CreateFastRecPoints(AliITSmodule *mod, Int_t module, TRandom *random){
  // Fast points simulator 

   AliITS *aliITS  = (AliITS*)gAlice->GetModule("ITS");
   AliITSgeom *gm = aliITS->GetITSgeom();

   const Float_t kdEdXtoQ = 2.778e+8; 

   Int_t ihit,flag,numofhits;
   Float_t xg,yg,zg,xl,yl,zl;
   Float_t px,py,pz;
   //Double_t p, theta, pt, ps;
   Float_t locals[3];
   Float_t globals[3];
   //Float_t xg1,yg1,zg1;
   Double_t sigmarphi=0., sigmaz=0., sigmade=0., thrde=0.;
   Float_t deltaXl,deltaZl,deltaDe;

   Int_t hitlay, hitlad, hitdet, hitstatus, hittrack;
   Float_t hitx, hity, hitz, hitpx, hitpy, hitpz, hitdestep;
   

   Int_t   hitstatus1, hittrack1;
   Float_t hitx1, hity1, hitz1;
   Float_t hitdestep1;

   Float_t xMg,yMg,zMg;
   //Float_t dx,dy,dz,ds;


   numofhits = mod->GetNhits();
   flag = 1;
   for(ihit=0;ihit<numofhits;ihit++){
     AliITShit *hit=mod->GetHit(ihit);
     hit->GetPositionG(hitx1,hity1,hitz1);
     hitstatus1 = hit->GetTrackStatus();
     hitdestep1 = hit->GetIonization();

     hittrack1 = hit->GetTrack();

     mod->MedianHit(module,hitx1,hity1,hitz1,hitstatus1,xMg,yMg,zMg,flag);
     if (flag!=1) {
       hitx      = xMg;
       hity      = yMg;
       hitz      = zMg;
       hit->GetDetectorID(hitlay,hitlad,hitdet);
       hit->GetMomentumG(hitpx,hitpy,hitpz);            
       hitdestep = hit->GetIonization();
       hitstatus = hitstatus1;
       hittrack  = hit->GetTrack();

       if (hitdestep > 0) {
	   xg = hitx;
	   yg = hity;
	   zg = hitz;
           // Transform to the module local frame
	   globals[0] = hitx;
	   globals[1] = hity;
	   globals[2] = hitz;
	   gm->GtoL(hitlay,hitlad,hitdet,globals,locals);
	   xl = locals[0];
	   yl = locals[1];
	   zl = locals[2];
	   px = hitpx;
	   py = hitpy;
	   pz = hitpz;
	   /*
           // Calculate transverse momentum and pseudorapidity
           // to allow pt and eta dependence in sigma values
           // of the spatial resolution
	   p  = TMath::Sqrt((px*px)+(py*py)+(pz*pz));
	   theta = TMath::ACos(pz/p);
	   pt = p * TMath::Sin(theta);
	   ps = -TMath::Log(TMath::Tan(theta/2));
	   */

           // Retrieve sigma values for position and energy, and energy
           // threshold 
	   
	   sigmarphi = SigmaRPhi(hitlay);
	   sigmaz = SigmaZ(hitlay);
	   sigmade = SigmaDe(hitlay);
	   thrde = ThrDe(hitlay);
	   // Randomize position and deposited energy
           Int_t k=3*(Int_t)((hitlay-1)/2);

	   deltaXl = (float)(random[k].Gaus(0,sigmarphi));
	   deltaZl = (float)(random[k+1].Gaus(0,sigmaz));
	   deltaDe = (float)(random[k+2].Gaus(0,sigmade));
           // Apply energy threshold and trasform back to global reference 
           // system
	   if ( (hitdestep+deltaDe) > thrde ){
	       locals[0] = xl + deltaXl;
	       locals[1] = yl;
	       locals[2] = zl + deltaZl;
	       AliITSRecPoint rp;
	       rp.fTracks[0]=hittrack;
	       rp.fTracks[1]=0;
	       rp.fTracks[2]=0;
	       rp.SetX(locals[0]);
	       rp.SetZ(locals[2]);
	       rp.SetdEdX(hitdestep+deltaDe);
	       rp.SetQ(kdEdXtoQ*(hitdestep+deltaDe));  // number of e
	       rp.SetSigmaX2(sigmarphi*sigmarphi);
	       rp.SetSigmaZ2(sigmaz*sigmaz);
	       rp.SetProbability(1.0);
	       aliITS->AddRecPoint(rp);
	       /*
	       gm->LtoG(hitlay,hitlad,hitdet,locals,globals);
	       xg1 = globals[0];
	       yg1 = globals[1];
	       zg1 = globals[2];
	       dx = TMath::Abs(xg1-hitx);
	       dy = TMath::Abs(yg1-hity);
	       dz = TMath::Abs(zg1-hitz);
	       ds = TMath::Abs(deltaDe);
	       */
	   } // if ( (hitdestep+deltaDe)
	   else flag=1;
       } // if (hitdestep > 0)
       else flag=1;
     } // if (flag!=1)
   }   

}
//_______________________________________________________________________
void AliITSsimulationFastPoints::SetSigmaRPhi(Double_t  srphi[6])
{
  // set sigmas in rphi

    Int_t i;
    for (i=0; i<6; i++) {
	fSigmaRPhi[i]=srphi[i];
    }
}
//_______________________________________________________________________
void AliITSsimulationFastPoints::SetSigmaZ(Double_t  sz[6])
{
  // set sigmas in z

    Int_t i;
    for (i=0; i<6; i++) {
	fSigmaZ[i]=sz[i];
    }
}
//_______________________________________________________________________
void AliITSsimulationFastPoints::SetSigmaDe(Double_t  sde[6])
{
  // set sigmas in energy

    Int_t i;
    for (i=0; i<6; i++) {
	fSigmaDe[i]=sde[i];
    }
}
//_______________________________________________________________________
void AliITSsimulationFastPoints::SetThrDe(Double_t  thrde[6])
{
  // set energy thersholds

    Int_t i;
    for (i=0; i<6; i++) {
	fThrDe[i]=thrde[i];
    }
}

