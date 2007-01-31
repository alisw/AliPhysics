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

// Base class for generators using external MC generators.
// For example AliGenPythia using Pythia.
// Provides basic functionality: setting of kinematic cuts on 
// decay products and particle selection.
// andreas.morsch@cern.ch

#include <TClonesArray.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TParticle.h>

#include "AliGenMC.h"
#include "AliGeometry.h"

ClassImp(AliGenMC)

AliGenMC::AliGenMC()
    :AliGenerator(),
     fParticles(0),
     fParentSelect(8),
     fChildSelect(8),
     fCutOnChild(0),
     fChildPtMin(0.),
     fChildPtMax(1.e10),
     fChildPMin(0.),
     fChildPMax(1.e10),
     fChildPhiMin(0.),
     fChildPhiMax(2. * TMath::Pi()),
     fChildThetaMin(0.),
     fChildThetaMax(TMath::Pi()),
     fChildYMin(-12.),
     fChildYMax(12.),
     fXingAngleX(0.),
     fXingAngleY(0.),
     fForceDecay(kAll),
     fMaxLifeTime(1.e-15),
     fAProjectile(1),
     fZProjectile(1),
     fATarget(1),
     fZTarget(1),
     fProjectile("P"),
     fTarget("P"),    
     fDyBoost(0.),
     fGeometryAcceptance(0),
     fPdgCodeParticleforAcceptanceCut(0),
     fNumberOfAcceptedParticles(2)
{
// Default Constructor
}

AliGenMC::AliGenMC(Int_t npart)
    :AliGenerator(npart),
     fParticles(0),
     fParentSelect(8),
     fChildSelect(8),
     fCutOnChild(0),
     fChildPtMin(0.),
     fChildPtMax(1.e10),
     fChildPMin(0.),
     fChildPMax(1.e10),
     fChildPhiMin(0.),
     fChildPhiMax(2. * TMath::Pi()),
     fChildThetaMin(0.),
     fChildThetaMax(TMath::Pi()),
     fChildYMin(-12.),
     fChildYMax(12.),
     fXingAngleX(0.),
     fXingAngleY(0.),
     fForceDecay(kAll),
     fMaxLifeTime(1.e-15),
     fAProjectile(1),
     fZProjectile(1),
     fATarget(1),
     fZTarget(1),
     fProjectile("P"),
     fTarget("P"),    
     fDyBoost(0.),
     fGeometryAcceptance(0),
     fPdgCodeParticleforAcceptanceCut(0),
     fNumberOfAcceptedParticles(2)
{
//  Constructor
// 
    for (Int_t i=0; i<8; i++) fParentSelect[i]=fChildSelect[i]=0;
}

AliGenMC::~AliGenMC()
{
// Destructor
}

void AliGenMC::Init()
{
//
//  Initialization
    switch (fForceDecay) {
    case kSemiElectronic:
    case kDiElectron:
    case kBJpsiDiElectron:
    case kBPsiPrimeDiElectron:
	fChildSelect[0] = kElectron;	
	break;
    case kHardMuons:	
    case kSemiMuonic:
    case kDiMuon:
    case kBJpsiDiMuon:
    case kBPsiPrimeDiMuon:
    case kPiToMu:
    case kKaToMu:
    case kWToMuon:
    case kWToCharmToMuon:
    case kZDiMuon:
	fChildSelect[0]=kMuonMinus;
	break;
    case kWToCharm:
	break;
    case kHadronicD:
	fChildSelect[0]=kPiPlus;
	fChildSelect[1]=kKPlus;
	break;
    case kPhiKK:
	fChildSelect[0]=kKPlus;
	break;
    case kBJpsi:
	fChildSelect[0]=443;
	break;
    case kOmega:	
    case kAll:
    case kNoDecay:
    case kNoDecayHeavy:
	break;
    }

    if (fZTarget != 0 && fAProjectile != 0) 
    {
	fDyBoost    = - 0.5 * TMath::Log(Double_t(fZProjectile) * Double_t(fATarget) / 
					 (Double_t(fZTarget)    * Double_t(fAProjectile)));
    }
}


Bool_t AliGenMC::ParentSelected(Int_t ip) const
{
// True if particle is in list of parent particles to be selected
    for (Int_t i=0; i<8; i++)
    {
	if (fParentSelect.At(i) == ip) return kTRUE;
    }
    return kFALSE;
}

Bool_t AliGenMC::ChildSelected(Int_t ip) const
{
// True if particle is in list of decay products to be selected
    for (Int_t i=0; i<5; i++)
    {
	if (fChildSelect.At(i) == ip) return kTRUE;
    }
    return kFALSE;
}

Bool_t AliGenMC::KinematicSelection(TParticle *particle, Int_t flag) const
{
// Perform kinematic selection
    Float_t pz    = particle->Pz();
    Float_t  e    = particle->Energy();
    Float_t pt    = particle->Pt();
    Float_t p     = particle->P();
    Float_t theta = particle->Theta();
    Float_t mass  = particle->GetCalcMass();
    Float_t mt2   = pt * pt + mass * mass;
    Float_t phi   = particle->Phi();
    
    Double_t y, y0;

    if (TMath::Abs(pz) <  e) {
	y = 0.5*TMath::Log((e+pz)/(e-pz));
    } else {
	y = 1.e10;
    }
    
    if (mt2) {
	y0 = 0.5*TMath::Log((e+TMath::Abs(pz))*(e+TMath::Abs(pz))/mt2);
    } else {
	if (TMath::Abs(y) < 1.e10) {
	    y0 = y;
	} else {
	    y0 = 1.e10;
	}
    }
      
    y = (pz < 0) ? -y0 : y0;
    
    if (flag == 0) {
//
// Primary particle cuts
//
//  transverse momentum cut    
	if (pt > fPtMax || pt < fPtMin) {
//	    printf("\n failed pt cut %f %f %f \n",pt,fPtMin,fPtMax);
	    return kFALSE;
	}
//
// momentum cut
	if (p > fPMax || p < fPMin) {
//	    printf("\n failed p cut %f %f %f \n",p,fPMin,fPMax);
	    return kFALSE;
	}
//
// theta cut
	if (theta > fThetaMax || theta < fThetaMin) {
//	    printf("\n failed theta cut %f %f %f \n",theta,fThetaMin,fThetaMax);
	    return kFALSE;
	}
//
// rapidity cut
	if (y > fYMax || y < fYMin) {
//	    printf("\n failed y cut %f %f %f \n",y,fYMin,fYMax);
	    return kFALSE;
	}
//
// phi cut
	if (phi > fPhiMax || phi < fPhiMin) {
//	    printf("\n failed phi cut %f %f %f \n",phi,fPhiMin,fPhiMax);
	    return kFALSE;
	}
    } else {
//
// Decay product cuts
//
//  transverse momentum cut    
	if (pt > fChildPtMax || pt < fChildPtMin) {
//	    printf("\n failed pt cut %f %f %f \n",pt,fChildPtMin,fChildPtMax);
	    return kFALSE;
	}
//
// momentum cut
	if (p > fChildPMax || p < fChildPMin) {
//	    printf("\n failed p cut %f %f %f \n",p,fChildPMin,fChildPMax);
	    return kFALSE;
	}
//
// theta cut
	if (theta > fChildThetaMax || theta < fChildThetaMin) {
//	    printf("\n failed theta cut %f %f %f \n",theta,fChildThetaMin,fChildThetaMax);
	    return kFALSE;
	}
//
// rapidity cut
	if (y > fChildYMax || y < fChildYMin) {
//	    printf("\n failed y cut %f %f %f \n",y,fChildYMin,fChildYMax);
	    return kFALSE;
	}
//
// phi cut
	if (phi > fChildPhiMax || phi < fChildPhiMin) {
//	    printf("\n failed phi cut %f %f %f \n",phi,fChildPhiMin,fChildPhiMax);
	    return kFALSE;
	}
    }
    
    return kTRUE;
}

Bool_t AliGenMC::CheckAcceptanceGeometry(Int_t np, TClonesArray* particles)
{
// Check the geometrical acceptance for particle.

  Bool_t check ;  
  Int_t numberOfPdgCodeParticleforAcceptanceCut = 0;
  Int_t numberOfAcceptedPdgCodeParticleforAcceptanceCut = 0;
  TParticle * particle;
  Int_t i;
  for (i = 0; i < np; i++) {
    particle =  (TParticle *) particles->At(i);
    if( TMath::Abs( particle->GetPdgCode() ) == TMath::Abs( fPdgCodeParticleforAcceptanceCut ) ) {
      numberOfPdgCodeParticleforAcceptanceCut++;
      if (fGeometryAcceptance->Impact(particle)) numberOfAcceptedPdgCodeParticleforAcceptanceCut++;
    }   
  }
  if ( numberOfAcceptedPdgCodeParticleforAcceptanceCut > (fNumberOfAcceptedParticles-1) )
    check = kTRUE;
  else
    check = kFALSE;

  return check;
}

Int_t AliGenMC::CheckPDGCode(Int_t pdgcode) const
{
//
//  If the particle is in a diffractive state, then take action accordingly
  switch (pdgcode) {
  case 91:
    return 92;
  case 110:
    //rho_diff0 -- difficult to translate, return rho0
    return 113;
  case 210:
    //pi_diffr+ -- change to pi+
    return 211;
  case 220:
    //omega_di0 -- change to omega0
    return 223;
  case 330:
    //phi_diff0 -- return phi0
    return 333;
  case 440:
    //J/psi_di0 -- return J/psi
    return 443;
  case 2110:
    //n_diffr -- return neutron
    return 2112;
  case 2210:
    //p_diffr+ -- return proton
    return 2212;
  }
  //non diffractive state -- return code unchanged
  return pdgcode;
}

void AliGenMC::Boost()
{
//
// Boost cms into LHC lab frame
//

    Double_t beta  = TMath::TanH(fDyBoost);
    Double_t gamma = 1./TMath::Sqrt(1.-beta*beta);
    Double_t gb    = gamma * beta;

    //    printf("\n Boosting particles to lab frame %f %f %f", fDyBoost, beta, gamma);
    
    Int_t i;
    Int_t np = fParticles->GetEntriesFast();
    for (i = 0; i < np; i++) 
    {
	TParticle* iparticle = (TParticle*) fParticles->At(i);

	Double_t e   = iparticle->Energy();
	Double_t px  = iparticle->Px();
	Double_t py  = iparticle->Py();
	Double_t pz  = iparticle->Pz();

	Double_t eb  = gamma * e -      gb * pz;
	Double_t pzb =   -gb * e +   gamma * pz;

	iparticle->SetMomentum(px, py, pzb, eb);
    }
}
