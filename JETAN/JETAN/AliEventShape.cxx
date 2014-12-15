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
   
/* $Id$ */
//---------------------------------------------------------------------
// Event shape utility class
// Circularity, Thrust, ... 
// Authors: Antonio Ortiz Velasquez <Antonio.Ortiz.Velasquez@cern.ch>
//          
//---------------------------------------------------------------------


#include "AliEventShape.h"

#include "AliStack.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TMatrixDSymEigen.h>
#include <TParticle.h>
#include <TParticlePDG.h>




//____________________________________________________________________
ClassImp(AliEventShape)

//___________________________________________________________________
TArrayD * AliEventShape::GetThrustParamMC(AliMCEvent* mcEvent, Int_t  nstudymin, Double_t ptcutoff, Double_t etacutoff, Bool_t chom)
{
  /*
    This function returns an array of values of thrust. To get these values you have to do:
    TArrayD* eventshapes = 0;
    eventshapes = AliShapeEvent::GetThrustParamMC(mcEvent, 3, 1, 1, kFALSE);
    Double_t thrust=eventshapes->GetAt(0);
    Double_t thrustmin=eventshapes->GetAt(1);
    Double_t recoil=eventshapes->GetAt(2);  
    The calculus uses  primary  particles. The input parameters:
    1. nstudymin, is the minumum number of particles which you want that participate in the calculus.(default:nstudymin =3)
    2. ptcutoff, is the cut in pt applied to participants to calculate the variables.(default: ptcutoff=1)
    3. etacutoff, is the cut in acceptance applied to participants to calculate the variables.(default: etacutoff=1)
    4. if chom=kTRUE, then the calculus includes neutral particles (rejecting photons and neutrinos).
       if chom=kFALSE, then the calculus includes only charged particles (rejecting photons and neutrinos).
       Returned values: thrust->0: 2-jet event, thrust->0.5: isotropic event
       Recoil is a term which is sensitive to radiation outside from acceptance, 1>=recoil>=0,
       thrustmin, is a measure of the radiation which is perpendicular to the plane formed by beam axis and thrust axis, 2/TMath::Pi()>thrustmin>0. In the limit of 2 back-to-back jets thrusmin->0, while in the case of a uniformly distributed event thrustmin->2/TMath::Pi();
  */
  
  AliStack* stack = 0;

  stack = mcEvent->Stack();
  Double_t * ptT = 0;
  Double_t * pxT = 0;
  Double_t * pyT = 0;
  Double_t ptsuma = 0;
  Double_t pxsuma = 0;
  Double_t pysuma = 0;

  TArrayD* evsh = new TArrayD(3);  
  Int_t nPrim  = stack->GetNprimary();
  Int_t  nmctracks = 0;  
  for (Int_t iMCTracks = 0; iMCTracks < nPrim; iMCTracks++) {    
      TParticle* trackmc = stack->Particle(iMCTracks);
      if (!trackmc) continue;
      Double_t etamc =trackmc ->Eta();
      Double_t ptmc=trackmc->Pt();
      Int_t pdgCode = TMath::Abs(trackmc->GetPdgCode());
      if (TMath::Abs(etamc) > etacutoff) continue; //only particles in |eta|<=etacutoff
      if(ptmc < ptcutoff) continue;  // PT cut
      Bool_t isprimary = stack->IsPhysicalPrimary(iMCTracks);    // Check if particle is charged, and primary
      if(isprimary == 0) continue;  // only primary particles
      TParticlePDG* pdgPart =trackmc ->GetPDG();
      if(chom == 1){//include neutral particles
	  // skip photons and neutrinos
	  if (pdgCode == 22 || pdgCode == 12 || pdgCode == 14 || pdgCode == 16) continue;
	  nmctracks++;  
      }
      else{ //only charged particles
	  if (pdgPart->Charge() == 0)continue;
	  nmctracks++;     
      }
  }
  // Minimum number of particles used in the analysis
  if(nmctracks < nstudymin){
      evsh->AddAt(-2,0);
      evsh->AddAt(-2,1);
      evsh->AddAt(-2,2);
      return evsh;
  }
  
  Int_t j=0;
  pxT = new Double_t[nmctracks];
  pyT = new Double_t[nmctracks];
  ptT = new Double_t[nmctracks];
  for (Int_t i = 0; i < nmctracks; i++)
    {
      pxT[i] = 0;
      pyT[i] = 0;
      ptT[i] = 0;
    }
  for (Int_t iMCTracks = 0; iMCTracks < nPrim; ++iMCTracks) {    
      TParticle* trackmc = stack->Particle(iMCTracks);
      if (!trackmc) continue;
      Double_t etamc = trackmc ->Eta();
      Double_t pxmc  = trackmc->Px();
      Double_t pymc  = trackmc->Py();
      Double_t ptmc  = trackmc->Pt();
      Int_t pdgCode  = TMath::Abs(trackmc->GetPdgCode());
      if (TMath::Abs(etamc) > etacutoff) continue;
      if(ptmc < ptcutoff) continue;
      Bool_t isprimary = stack->IsPhysicalPrimary(iMCTracks); 
      if(isprimary==0) continue;
      TParticlePDG* pdgPart =trackmc ->GetPDG();

      if(chom==1){
	  if (pdgCode == 22 || pdgCode == 12 || pdgCode == 14 || pdgCode == 16)continue;
      } else {
	  if (pdgPart->Charge() == 0) continue;
      }
      
      ptT[j] = ptmc;
      pxT[j] = pxmc;
      pyT[j] = pymc;
      ptsuma += ptmc;
      pxsuma+=pxmc;
      pysuma+=pymc;
      j++;    
  }

  Double_t numerador = 0;
  Double_t numerador2 = 0;
  Double_t phimax = -1;  
  Double_t pFull = -1;
  Double_t pMax = 0;
  Double_t phi = 0;
  Double_t thrust = 80;
  Double_t thrustminor = 80;
  Double_t nx = 0;
  Double_t ny = 0;
  Double_t phiparam = 0;
  //Getting thrust
  for(Int_t i = 0; i < 360; ++i){
      numerador = 0;
      phiparam  = 0;
      nx = 0;
      ny = 0;
      phiparam=((TMath::Pi()) * i) / 180; // parametrization of the angle
      nx = TMath::Cos(phiparam);            // x component of an unitary vector n
      ny = TMath::Sin(phiparam);            // y component of an unitary vector n
      for(Int_t i1 = 0; i1 < nmctracks; ++i1){
	  numerador += TMath::Abs(nx * pxT[i1] + ny * pyT[i1]);//product between momentum proyection in XY plane and the unitari vector.
      }
      pFull=numerador / ptsuma;
      if(pFull > pMax)//maximization of pFull
      {
	  pMax = pFull;
	  phi = phiparam;
      }
  }

  phimax=(phi * 180) / TMath::Pi();//angular parameter of the unitary vector which maximiza thrust
  //if n vector and beam axis form a plane, then we can calculate a second unitary vector perpendicular to that plane
  Double_t nx1 = TMath::Cos(phi);
  Double_t ny1 = TMath::Sin(phi);
  for(Int_t i2 =0; i2 < nmctracks; ++i2){
      numerador2 += TMath::Abs(pxT[i2] * ny1 - nx1 * pyT[i2]);//cross product: P_{i} X n, P_{i}=(px_{i},py_{i})
  }
  thrust = 1 - pMax;//this is the value of thrust
  thrustminor = numerador2 / ptsuma;//this is the value of thrust minor
  Double_t recoil = TMath::Abs(TMath::Sqrt(pxsuma * pxsuma + pysuma * pysuma)) / (ptsuma);//factor sentsitive to radiation outside from acceptance 

  evsh->AddAt(thrust, 0);
  evsh->AddAt(thrustminor, 1);
  evsh->AddAt(recoil, 2);


  delete [] ptT;
  delete [] pxT;
  delete [] pyT;

  return evsh;  
}


Double_t AliEventShape::GetCircularityMC(AliMCEvent* mcEvent, Int_t  nstudymin, Double_t ptcutoff, Double_t etacutoff, Bool_t chom)
{
  /*
    This function returns the circularity value of the event 

    The calculus uses  primary  particles. The input parameters:
    1. nstudymin, is the minumum number of particles which you want that participate in the calculus.(default:nstudymin =3)
    2. ptcutoff, is the cut in pt applied to participants to calculate the variables.(default: ptcutoff=1)
    3. etacutoff, is the cut in acceptance applied to participants to calculate the variables.(default: etacutoff=1)
    4. if chom=kTRUE, then the calculus includes neutral particles (rejecting photons and neutrinos).
       if chom=kFALSE, then the calculus includes only charged particles (rejecting photons and neutrinos).
       1>=circularity>=0
  */


  AliStack* stack = 0;

  stack = mcEvent->Stack();

  TMatrixDSym s(2);

  Double_t s00 =  0;
  Double_t s01 =  0;
  Double_t s10 =  0;
  Double_t s11 =  0;
  Double_t ptot = 0;
  Double_t circularity = -2;
  Int_t  nmctracks = 0;  
  Int_t nPrim  = stack->GetNprimary();

  for (Int_t iMCTracks = 0; iMCTracks < nPrim; ++iMCTracks) {
    TParticle* trackmc = stack->Particle(iMCTracks);
    if (!trackmc) continue;
    Double_t etamc = trackmc ->Eta();
    Double_t ptmc  = trackmc->Pt();
    Double_t pxmc  = trackmc->Px();
    Double_t pymc  = trackmc->Py();
    Int_t pdgCode = TMath::Abs(trackmc->GetPdgCode());
    if (TMath::Abs(etamc) > etacutoff) continue;
    if (ptmc < ptcutoff) continue;
    Bool_t isprimary = stack->IsPhysicalPrimary(iMCTracks);
    if (isprimary == 0) continue;
    TParticlePDG* pdgPart = trackmc ->GetPDG();
    if(chom == kTRUE){
      // skip photons and neutrinos
      if (pdgCode == 22 || pdgCode == 12 || pdgCode == 14 || pdgCode == 16) continue; 
    }
    else{
      if (pdgPart->Charge() == 0)continue;
    }

    ptot = ptot + (ptmc * ptmc);
    s00 = s00 + (pxmc * pxmc);
    s01 = s01 + (pxmc * pymc);
    s10 = s10 + (pymc * pxmc);
    s11 = s11 + (pymc * pymc); 
    nmctracks++;
  } //track loop 


  if (nmctracks < nstudymin) {
      Printf("Too few particles, stopping");
      return -2;
  }


  
  if(ptot != 0){
    s(0,0) = s00 / ptot;
    s(0,1) = s01 / ptot;
    s(1,0) = s10 / ptot;
    s(1,1) = s11 / ptot;
    const TMatrixDSymEigen eigen(s);
    const TVectorD eigenVal=eigen.GetEigenValues();
    circularity = 2 * (1 - eigenVal(0));
  }
  return circularity;
}
  
  































