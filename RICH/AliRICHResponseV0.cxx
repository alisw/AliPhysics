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

/*
  $Log$
  Revision 1.1  2000/06/12 15:29:37  jbarbosa
  Cleaned up version.

*/

#include "AliRICHResponseV0.h"
#include "AliRICHSegmentation.h"
#include "AliRun.h"
#include "AliMC.h"

#include <TMath.h>
#include <TRandom.h>
#include <TParticle.h>
//___________________________________________
ClassImp(AliRICHResponseV0)

Float_t AliRICHResponseV0::IntPH(Float_t eloss)
{
    // Get number of electrons and return charge
    
    Int_t nel;
    nel= Int_t(eloss/fEIonisation);
    
    Float_t charge=0;
    if (nel == 0) nel=1;
    for (Int_t i=1;i<=nel;i++) {
	charge -= fChargeSlope*TMath::Log(gRandom->Rndm());    
    }
    return charge;
}

Float_t AliRICHResponseV0::IntPH()
{

//  Get number of electrons and return charge, for a single photon

    Float_t charge = -fChargeSlope*TMath::Log(gRandom->Rndm());
    return charge;
}



// -------------------------------------------
Float_t AliRICHResponseV0::IntXY(AliRICHSegmentation * segmentation)
{
    
    const Float_t kInversePitch = 1/fPitch;
    Float_t response;
//
//  Integration limits defined by segmentation model
//  
    
    Float_t xi1, xi2, yi1, yi2;
    segmentation->IntegrationLimits(xi1,xi2,yi1,yi2);

    xi1=xi1*kInversePitch;
    xi2=xi2*kInversePitch;
    yi1=yi1*kInversePitch;
    yi2=yi2*kInversePitch;

    //printf("Integration Limits: %f-%f, %f-%f\n",xi1,xi2,yi1,yi2);
    
    //printf("KInversePitch:%f\n",kInversePitch);

    //
// The Mathieson function 
    Double_t ux1=fSqrtKx3*TMath::TanH(fKx2*xi1);
    Double_t ux2=fSqrtKx3*TMath::TanH(fKx2*xi2);
    
    Double_t uy1=fSqrtKy3*TMath::TanH(fKy2*yi1);
    Double_t uy2=fSqrtKy3*TMath::TanH(fKy2*yi2);

    //printf("Integration Data: %f-%f, %f-%f\n",ux1,ux2,uy1,uy2);
    
    //printf("%f %f %f %f\n",fSqrtKx3,fKx2,fKy4,fKx4);
    
    response=4.*fKx4*(TMath::ATan(ux2)-TMath::ATan(ux1))*fKy4*(TMath::ATan(uy2)-TMath::ATan(uy1));

    //printf("Response:%f\n",response);

    return response;       
    
}

Int_t AliRICHResponseV0::FeedBackPhotons(Float_t *source, Float_t qtot)
{
  //
  // Generate FeedBack photons
  //
  Int_t j, ipart, nt;
    
  Int_t sNfeed=0;
  
  
  // Local variables 
  Float_t cthf, ranf[2], phif, enfp = 0, sthf;
  Int_t i, ifeed;
  Float_t e1[3], e2[3], e3[3];
  Float_t vmod, uswop;
  Float_t fp, random;
  Float_t dir[3], phi;
  Int_t nfp;
  Float_t pol[3], mom[3];
  TLorentzVector position;
  //
  // Determine number of feedback photons

  //  Get weight of current particle
  TParticle *current = (TParticle*) 
    (*gAlice->Particles())[gAlice->CurrentTrack()];
    
  ifeed = Int_t(current->GetWeight()/100+0.5);
  ipart = gMC->TrackPid();
  fp = fAlphaFeedback * qtot;
  nfp = gRandom->Poisson(fp);
  
  // This call to fill the time of flight
  gMC->TrackPosition(position);
  //
  // Generate photons
  for (i = 0; i <nfp; i++) {
	
    // Direction
    gMC->Rndm(ranf, 2);
    cthf = ranf[0] * 2 - 1.;
    if (cthf < 0)  continue;
    sthf = TMath::Sqrt((1 - cthf) * (1 + cthf));
    phif = ranf[1] * 2 * TMath::Pi();
    //
    gMC->Rndm(&random, 1);
    if (random <= .57) {
      enfp = 7.5e-9;
    } else if (random <= .7) {
      enfp = 6.4e-9;
    } else {
      enfp = 7.9e-9;
    }

    dir[0] = sthf * TMath::Sin(phif);
    dir[1] = cthf;
    dir[2] = sthf * TMath::Cos(phif);
    gMC->Gdtom(dir, mom, 2);
    mom[0]*=enfp;
    mom[1]*=enfp;
    mom[2]*=enfp;
    
    // Polarisation
    e1[0] = 0;
    e1[1] = -dir[2];
    e1[2] = dir[1];
    
    e2[0] = -dir[1];
    e2[1] = dir[0];
    e2[2] = 0;
    
    e3[0] = dir[1];
    e3[1] = 0;
    e3[2] = -dir[0];
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e1[j];
      e1[j]=e3[j];
      e3[j]=uswop;
    }
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e2[j];
      e2[j]=e3[j];
      e3[j]=uswop;
    }
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    vmod=TMath::Sqrt(1/vmod);
    for(j=0;j<3;j++) e1[j]*=vmod;
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    vmod=TMath::Sqrt(1/vmod);
    for(j=0;j<3;j++) e2[j]*=vmod;
    
    gMC->Rndm(ranf, 1);
    phi = ranf[0] * 2 * TMath::Pi();
    for(j=0;j<3;j++) pol[j]=e1[j]*TMath::Sin(phi)+e2[j]*TMath::Cos(phi);
    gMC->Gdtom(pol, pol, 2);
    
    // Put photon on the stack and label it as feedback (51, 52) 
    ++sNfeed;

    gAlice->SetTrack(Int_t(1), gAlice->CurrentTrack(), Int_t(50000051),
		     mom,source,pol,position[3],
		     "Feedback", nt, 1.);
  }
  return(sNfeed);
}




