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
Revision 1.13  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.5.4.2  2002/07/24 08:56:28  alibrary
Updating EVGEN on TVirtulaMC

Revision 1.12  2002/07/19 11:42:33  morsch
Use CalcMass()

Revision 1.11  2002/06/06 15:26:24  morsch
Correct child-selection for kPhiKK

Revision 1.10  2002/06/05 14:05:46  morsch
Decayer option kPhiKK for forced phi->K+K- decay added.

Revision 1.9  2002/05/30 14:58:29  morsch
Add pointer to AliGeometry to handle geometrical acceptance. (G. MArtinez)

Revision 1.8  2002/04/26 10:42:35  morsch
Case kNoDecayHeavy added. (N. Carrer)

Revision 1.7  2002/04/17 10:32:32  morsch
Coding Rule violations corrected.

Revision 1.6  2002/03/26 14:19:36  morsch
Saver calculation of rapdity.

Revision 1.5  2002/03/12 17:02:20  morsch
Change in calculation of rapidity, include case in which numerically e == pz.

Revision 1.4  2001/11/27 13:13:07  morsch
Maximum lifetime for long-lived particles to be put on the stack is parameter.
It can be set via SetMaximumLifetime(..).

Revision 1.3  2001/10/16 08:48:56  morsch
Common vertex related code moved to base class AliGenerator.

Revision 1.2  2001/10/15 08:15:51  morsch
Event vertex and vertex truncation setting moved into AliMC.

Revision 1.1  2001/07/13 10:56:00  morsch
AliGenMC base class for AliGenParam and AliGenPythia commonalities.

*/

// Base class for generators using external MC generators.
// For example AliGenPythia using Pythia.
// Provides basic functionality: setting of kinematic cuts on 
// decay products and particle selection.
// andreas.morsch@cern.ch

#include <TMath.h>
#include <TPDGCode.h>
#include <TParticle.h>

#include "AliGenMC.h"

 ClassImp(AliGenMC)

AliGenMC::AliGenMC()
                 :AliGenerator()
{
// Default Constructor
    SetCutOnChild();
    SetChildMomentumRange();
    SetChildPtRange();
    SetChildPhiRange();
    SetChildThetaRange(); 
    SetChildYRange(); 
    SetMaximumLifetime();
    SetGeometryAcceptance();
    SetPdgCodeParticleforAcceptanceCut();
    SetNumberOfAcceptedParticles();
}

AliGenMC::AliGenMC(Int_t npart)
                 :AliGenerator(npart)
{
//  Constructor
    SetCutOnChild();
    SetChildMomentumRange();
    SetChildPtRange();
    SetChildPhiRange();
    SetChildThetaRange();
    SetChildYRange(); 
// 
    fParentSelect.Set(8);
    fChildSelect.Set(8);
    for (Int_t i=0; i<8; i++) fParentSelect[i]=fChildSelect[i]=0;
    SetMaximumLifetime();
    SetGeometryAcceptance();
    SetPdgCodeParticleforAcceptanceCut();
    SetNumberOfAcceptedParticles();
}

AliGenMC::AliGenMC(const AliGenMC & mc)
{
// copy constructor
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
    case kSemiMuonic:
    case kDiMuon:
    case kBJpsiDiMuon:
    case kBPsiPrimeDiMuon:
    case kPiToMu:
    case kKaToMu:
	fChildSelect[0]=kMuonMinus;
	break;
    case kHadronicD:
	fChildSelect[0]=kPiPlus;
	fChildSelect[1]=kKPlus;
	break;
    case kPhiKK:
	fChildSelect[0]=kKPlus;
    case kOmega:	
    case kAll:
    case kNoDecay:
    case kNoDecayHeavy:
	break;
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
    Float_t px    = particle->Px();
    Float_t py    = particle->Py();
    Float_t pz    = particle->Pz();
    Float_t  e    = particle->Energy();
    Float_t pt    = particle->Pt();
    Float_t p     = particle->P();
    Float_t theta = particle->Theta();
    Float_t mass  = particle->GetCalcMass();
    Float_t mt2   = pt * pt + mass * mass;
    
    Float_t phi   = Float_t(TMath::ATan2(Double_t(py),Double_t(px)));
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
  Bool_t Check ;  // All fPdgCodeParticleforAcceptanceCut particles are in in the fGeometryAcceptance acceptance
  Int_t NumberOfPdgCodeParticleforAcceptanceCut=0;
  Int_t NumberOfAcceptedPdgCodeParticleforAcceptanceCut=0;
  TParticle * particle;
  Int_t i;
  for (i=0; i<np; i++) {
    particle =  (TParticle *) particles->At(i);
    if( TMath::Abs( particle->GetPdgCode() ) == TMath::Abs( fPdgCodeParticleforAcceptanceCut ) ) {
      NumberOfPdgCodeParticleforAcceptanceCut++;
      if (fGeometryAcceptance->Impact(particle)) NumberOfAcceptedPdgCodeParticleforAcceptanceCut++;
    }   
  }
  if ( NumberOfAcceptedPdgCodeParticleforAcceptanceCut > (fNumberOfAcceptedParticles-1) )
    Check = kTRUE;
  else
    Check = kFALSE;

  return Check;
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
	  
AliGenMC& AliGenMC::operator=(const  AliGenMC& rhs)
{
// Assignment operator
    return *this;
}

