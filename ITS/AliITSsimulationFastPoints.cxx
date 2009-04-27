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

/* $Id$ */
//////////////////////////////////////////////////////////
// implements fast simulation                           //
//                                                      //
//                                                      //
//////////////////////////////////////////////////////////


#include <TRandom.h>

#include "AliITS.h"
#include "AliITShit.h"
#include "AliITSRecPoint.h"
#include "AliITSmodule.h"
#include "AliITSgeom.h"
#include "AliRun.h"
#include "AliITSsimulationFastPoints.h"


ClassImp(AliITSsimulationFastPoints)

AliITSsimulationFastPoints::AliITSsimulationFastPoints()
{
  //constructor
  fSigmaRPhi[0] = fSigmaRPhi[1] = 12e-4;
  fSigmaRPhi[2] = fSigmaRPhi[3] = 38e-4;
  fSigmaRPhi[4] = fSigmaRPhi[5] = 20e-4;
  fSigmaZ[0] = fSigmaZ[1] = 120e-4;        // resolution for 425 micron pixels
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
void AliITSsimulationFastPoints::CreateFastRecPoints(Int_t module, TClonesArray* recp){
    // Fast points simulator
    AliITS *aliITS  = (AliITS*)gAlice->GetModule("ITS");

    CreateFastRecPoints((AliITSmodule *)(aliITS->GetModule(module)),
			module,gRandom,recp);
}
//-------------------------------------------------------------
void AliITSsimulationFastPoints::CreateFastRecPoints(AliITSmodule *mod,
						     Int_t module,
						     TRandom *random,
						     TClonesArray* recp) {
  // Fast points simulator 

  TClonesArray &pt=*recp;
  AliITS *aliITS  = (AliITS*)gAlice->GetModule("ITS");
  AliITSgeom *gm = aliITS->GetITSgeom();
  const Float_t kdEdXtoQ = 1.0e+6;  // GeV->KeV

  Int_t lay,lad,det;
  gm->GetModuleId(module,lay,lad,det);
  Int_t ind=(lad-1)*gm->GetNdetectors(lay)+(det-1);
  Int_t lyr=(lay-1);


  Int_t ihit,flag,numofhits;
  Float_t locals[3];
  Float_t globals[3];
  Double_t sigmarphi=0., sigmaz=0., sigmade=0., thrde=0.;
  Float_t deltaXl,deltaZl,deltaDe;

  Int_t hitlay, hitlad, hitdet, hitstatus;
  Float_t hitpx, hitpy, hitpz, hitdestep;

  Int_t   hitstatus1, hittrack1;
  Float_t hitx1, hity1, hitz1;
  Float_t hitdestep1;
  Float_t xMg,yMg,zMg;
  Int_t irecp=0;
  numofhits = mod->GetNhits();
  //printf("numofhits %d \n",numofhits);
  for(ihit=0;ihit<numofhits;ihit++){
    AliITShit *hit=mod->GetHit(ihit);
    hit->GetPositionG(hitx1,hity1,hitz1);
    hitstatus1 = hit->GetTrackStatus();
    hitdestep1 = hit->GetIonization();
    hittrack1 = hit->GetTrack();
    
    mod->MedianHit(module,hitx1,hity1,hitz1,hitstatus1,xMg,yMg,zMg,flag);
    if (flag!=1) {
      hitdestep = hit->GetIonization();
      
      if (hitdestep > 0) {
	hit->GetDetectorID(hitlay,hitlad,hitdet);
	hit->GetMomentumG(hitpx,hitpy,hitpz);            
	hitstatus = hitstatus1;
		// Transform to the module local frame
	globals[0] = xMg; 
	globals[1] = yMg;
	globals[2] = zMg;
	gm->GtoL(hitlay,hitlad,hitdet,globals,locals);
	// Retrieve sigma values for position and energy, and energy
	// threshold
	sigmarphi = SigmaRPhi(hitlay);
	sigmaz = SigmaZ(hitlay);
	sigmade = SigmaDe(hitlay);
	thrde = ThrDe(hitlay);
	deltaXl = random->Gaus(0,sigmarphi);
	deltaZl = random->Gaus(0,sigmaz);
	deltaDe = random->Gaus(0,sigmade);
	
	// Apply energy threshold and trasform back to global reference
	// system
	
	if ( (hitdestep+deltaDe) > thrde ){
	  locals[0] += deltaXl;
	  locals[2] += deltaZl;
	  Int_t lab[4] = {hit->GetTrack(),-3,-3,ind};
	  Float_t q=kdEdXtoQ*(hitdestep+deltaDe);
	  if(hitlay<3) q=1.; // SPD binary readout
	  Float_t hitv[5] = {locals[0],locals[2],sigmarphi*sigmarphi,sigmaz*sigmaz,q};
	  Int_t info[3] = {0,0,lyr};
	  AliITSRecPoint rp(lab,hitv,info,kTRUE);

	  new (pt[irecp]) AliITSRecPoint(rp);
	  irecp++;
	} // end if ( (hitdestep+deltaDe)
      } // end if (hitdestep > 0)
    } // end if (flag!=1)
  } // end for ihit
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

