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
Revision 1.9.8.1  2002/07/24 09:27:50  alibrary
Updating on VirtualMC

Revision 1.11  2002/07/11 10:24:21  barbera
Fixes to make tracking V2 working with the HEAD and with fast points. Waiting for a fix in slow reconstruction of SPD, fast points are temporarily made the default for tracking V2.

Revision 1.10  2002/06/10 17:30:24  nilsen
A new CreateFastRecPoints has been made and the old one made compatible.

Revision 1.9  2001/10/01 19:36:03  nilsen
fixed a compilation warning about unused variable.

Revision 1.8  2001/07/27 08:06:49  hristov
Use global gRandom generator (M.Ivanov)

Revision 1.7  2001/05/11 09:15:21  barbera
Corrected to make fast point creation working with PPR geometry

Revision 1.6  2000/10/29 18:30:14  barbera
Z resolution of pixel changed according with the default lenght of 425 microns

Revision 1.1.2.6  2000/10/29 18:29:51  barbera
Z resolution of pixel changed according with the default lenght of 425 microns

Revision 1.1.2.5  2000/10/02 16:03:20  barbera
Forward declarations added

Revision 1.4  2000/09/22 12:43:59  nilsen
Default track number set to -3 and not 0.

*/
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
void AliITSsimulationFastPoints::CreateFastRecPoints(Int_t module){
    // Fast points simulator
    AliITS *aliITS  = (AliITS*)gAlice->GetModule("ITS");

    CreateFastRecPoints((AliITSmodule *)(aliITS->GetModule(module)),
			module,gRandom);
}
//-------------------------------------------------------------
void AliITSsimulationFastPoints::CreateFastRecPoints(AliITSmodule *mod,
						     Int_t module,
						     TRandom *random){
    // Fast points simulator 
    AliITS *aliITS  = (AliITS*)gAlice->GetModule("ITS");
    AliITSgeom *gm = aliITS->GetITSgeom();

    const Float_t kdEdXtoQ = 2.778e+8; 

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
		    AliITSRecPoint rp;
		    rp.fTracks[0]=hit->GetTrack();
		    //		    rp.fTracks[0]=mod->GetHitTrackIndex(ihit);
		    rp.fTracks[1]=-3;
		    rp.fTracks[2]=-3;
		    rp.SetX(locals[0]);
		    rp.SetZ(locals[2]);
		    rp.SetdEdX(hitdestep+deltaDe);
		    rp.SetQ(kdEdXtoQ*(hitdestep+deltaDe));  // number of e
		    rp.SetSigmaX2(sigmarphi*sigmarphi);
		    rp.SetSigmaZ2(sigmaz*sigmaz);
		    aliITS->AddRecPoint(rp);
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

