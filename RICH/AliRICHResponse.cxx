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


#include <TMath.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TVirtualMC.h>

#include "AliRICHResponse.h"
#include "AliRun.h"
#include "AliSegmentation.h"
#include "AliMC.h"

ClassImp(AliRICHResponse)
//__________________________________________________________________________________________________
AliRICHResponse::AliRICHResponse()
{
   SetSigmaIntegration(5.);
   SetChargeSlope(27.);
   SetChargeSpread(0.18, 0.18);
   SetMaxAdc(4096);
   SetAlphaFeedback(0.036);
   SetEIonisation(26.e-9);
   SetSqrtKx3(0.77459667);
   SetKx2(0.962);
   SetKx4(0.379);
   SetSqrtKy3(0.77459667);
   SetKy2(0.962);
   SetKy4(0.379);
   SetPitch(0.25);
   SetWireSag(1);		      // 1->On, 0->Off
   SetVoltage(2150);		      // Should only be 2000, 2050, 2100 or 2150
}//AliRICHResponse::ctor()
//__________________________________________________________________________________________________
Float_t AliRICHResponse::IntPH(Float_t eloss, Float_t yhit)
{
    // Get number of electrons and return charge
    
    Int_t nel;
    nel= Int_t(eloss/fEIonisation);
    
    Float_t charge=0;
    Double_t gain_var=1;

    if (nel == 0) nel=1;

    if (fWireSag)
      {
	if (fVoltage==2150)
	  {
	    gain_var = 9e-6*TMath::Power(yhit,4) + 2e-7*TMath::Power(yhit,3) - 0.0316*TMath::Power(yhit,2) - 3e-4*yhit + 25.367;
	  }
	if (fVoltage==2100)
	    gain_var = 8e-6*TMath::Power(yhit,4) + 2e-7*TMath::Power(yhit,3) - 0.0283*TMath::Power(yhit,2) - 2e-4*yhit + 23.015;
	if (fVoltage==2050)
	    gain_var = 7e-6*TMath::Power(yhit,4) + 1e-7*TMath::Power(yhit,3) - 0.0254*TMath::Power(yhit,2) - 2e-4*yhit + 20.888;
	if (fVoltage==2000)
	    gain_var = 6e-6*TMath::Power(yhit,4) + 8e-8*TMath::Power(yhit,3) - 0.0227*TMath::Power(yhit,2) - 1e-4*yhit + 18.961;
		
	gain_var = gain_var/100;
	//printf("Yhit:%f, Gain variation:%f\n",yhit,gain_var);

	Float_t gain = (fChargeSlope + fChargeSlope*gain_var)*.9; 
	//printf(" Yhit:%f, Gain variation:%f\n",yhit, gain);

	for (Int_t i=1;i<=nel;i++) {
	  charge -= gain*TMath::Log(gRandom->Rndm());    
	}
      }
    else
      {
	for (Int_t i=1;i<=nel;i++) {
	  charge -= fChargeSlope*TMath::Log(gRandom->Rndm());    
	}
      }

    return charge;
}//InitPH()
//__________________________________________________________________________________________________
Float_t AliRICHResponse::IntPH(Float_t yhit)
{//  Get number of electrons and return charge, for a single photon

  Float_t charge=0;
  Double_t gain_var=1;

   if (fWireSag)
      {
	if (fVoltage==2150)
	  {
	    gain_var = 9e-6*TMath::Power(yhit,4) + 2e-7*TMath::Power(yhit,3) - 0.0316*TMath::Power(yhit,2) - 3e-4*yhit + 25.367;
	    //gain_var = 9e-5*TMath::Power(yhit,4) + 2e-6*TMath::Power(yhit,3) - 0.316*TMath::Power(yhit,2) - 3e-3*yhit + 253.67;
	  }
	if (fVoltage==2100)
	    gain_var = 8e-6*TMath::Power(yhit,4) + 2e-7*TMath::Power(yhit,3) - 0.0283*TMath::Power(yhit,2) - 2e-4*yhit + 23.015;
	if (fVoltage==2050)
	    gain_var = 7e-6*TMath::Power(yhit,4) + 1e-7*TMath::Power(yhit,3) - 0.0254*TMath::Power(yhit,2) - 2e-4*yhit + 20.888;
	if (fVoltage==2000)
	    gain_var = 6e-6*TMath::Power(yhit,4) + 8e-8*TMath::Power(yhit,3) - 0.0227*TMath::Power(yhit,2) - 1e-4*yhit + 18.961;

	gain_var = gain_var/100;
	//printf(" Yhit:%f, Gain variation:%f\n",yhit, gain_var);
	
	Float_t gain = (fChargeSlope + fChargeSlope*gain_var)*.9; 
	
	charge -= gain*TMath::Log(gRandom->Rndm());
	//printf(" Yhit:%f, Gain variation:%f\n",yhit, gain);
      }
   else
     {
       charge -= fChargeSlope*TMath::Log(gRandom->Rndm());
     }
    return charge;
}//IntPH()
//__________________________________________________________________________________________________
Float_t AliRICHResponse::IntXY(AliSegmentation * segmentation)
{
    
    const Float_t kInversePitch = 1/fPitch;
    Float_t response;

    //  Integration limits defined by segmentation model
    
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
    
}//IntXY
//__________________________________________________________________________________________________
Int_t AliRICHResponse::FeedBackPhotons(Float_t *source, Float_t qtot)
{  // Generate FeedBack photons
  Int_t j, ipart, nt;
    
  Int_t sNfeed=0;
  
  
  // Local variables 
  Double_t ranf[2];
  Float_t cthf, phif, enfp = 0, sthf;
  Int_t i, ifeed;
  Float_t e1[3], e2[3], e3[3];
  Float_t vmod, uswop;
  Float_t fp;
  Double_t random;
  Float_t dir[3], phi;
  Int_t nfp;
  Float_t pol[3], mom[4];
  TLorentzVector position;
  //
  // Determine number of feedback photons

  //  Get weight of current particle
  TParticle *current = (TParticle*) 
    (*gAlice->GetMCApp()->Particles())[gAlice->GetMCApp()->GetCurrentTrackNumber()];
    
  ifeed = Int_t(current->GetWeight()/100+0.5);
  ipart = gMC->TrackPid();
  fp = fAlphaFeedback * qtot;
  nfp = gRandom->Poisson(fp);
  
  // This call to fill the time of flight
  gMC->TrackPosition(position);
  //printf("Track position: %f %f %f %15.12f\n", position[0],position[1],position[2],position[3]);
  //
  // Generate photons
  for (i = 0; i <nfp; i++) {
	
    // Direction
    gMC->GetRandom()->RndmArray(2,ranf);
    cthf = ranf[0] * 2 - 1.;
    if (cthf < 0)  continue;
    sthf = TMath::Sqrt((1 - cthf) * (1 + cthf));
    phif = ranf[1] * 2 * TMath::Pi();
    //
    //gMC->Rndm(&random, 1);
    gMC->GetRandom()->RndmArray(1, &random);
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
    mom[3] = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
    //printf("Dir %f %f %f\n",dir[0],dir[1],dir[2]);
    //printf("Momentum %15.12f %15.12f %15.12f\n",mom[0],mom[1],mom[2]);
    //printf("Energy %e\n", mom[3]);
    
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
    
    //gMC->Rndm(ranf, 1);
    gMC->GetRandom()->RndmArray(1,ranf);
    phi = ranf[0] * 2 * TMath::Pi();
    for(j=0;j<3;j++) pol[j]=e1[j]*TMath::Sin(phi)+e2[j]*TMath::Cos(phi);
    gMC->Gdtom(pol, pol, 2);
    
    // Put photon on the stack and label it as feedback (51, 52) 
    ++sNfeed;

    gAlice->GetMCApp()->PushTrack(Int_t(1), gAlice->GetMCApp()->GetCurrentTrackNumber(), Int_t(50000051),
		     mom[0],mom[1],mom[2],mom[3],source[0],source[1],source[2],position[3],pol[0],pol[1],pol[2],
		     kPFeedBackPhoton, nt, 1.);
    
  }
  return(sNfeed);
}//FeedBackPhotons()
