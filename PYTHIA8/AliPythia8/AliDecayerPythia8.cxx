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

// Implementation of AliDecayer using Pythia8
// Author: andreas.morsch@cern.ch
#include <TMath.h>
#include <TRandom.h>
#include <TPDGCode.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include "AliTPythia8.h"
#include "AliDecayerPythia8.h"
#include "ParticleData.h"

ClassImp(AliDecayerPythia8)

Bool_t AliDecayerPythia8::fgInit = kFALSE;

AliDecayerPythia8::AliDecayerPythia8():
  TVirtualMCDecayer(),
  fPythia8(AliTPythia8::Instance()),
  fDebug(0),
  fDecay(kAll),
  fHeavyFlavour(kTRUE)
{
    // Constructor
   fPythia8->Pythia8()->readString("SoftQCD:elastic = on");
   fPythia8->Pythia8()->init();
}

//___________________________________________________________________________
void AliDecayerPythia8::Decay(Int_t pdg, TLorentzVector* p)
{
   // Decay a single particle
   ClearEvent();
   AppendParticle(pdg, p);
   Int_t idPart = fPythia8->Pythia8()->event[0].id();
   fPythia8->Pythia8()->particleData.mayDecay(idPart,kTRUE);
   fPythia8->Pythia8()->moreDecays();
   if (fDebug > 0) fPythia8->EventListing();
}

//___________________________________________________________________________
Int_t AliDecayerPythia8::ImportParticles(TClonesArray *particles)
{
  //import the decay products into particles array
  const Float_t kconvT=0.001/2.999792458e8; // mm/c to seconds conversion
  const Float_t kconvL=1./10; // mm to cm conversion
  int np = fPythia8->ImportParticles(particles, "All");
  // pythia assigns decay time in mm/c, convert to seconds
  for (int ip=1;ip<np;ip++) {
    TParticle* prod = (TParticle*)particles->At(ip);
    if (!prod) continue;
    prod->SetProductionVertex(prod->Vx()*kconvL,prod->Vy()*kconvL,prod->Vz()*kconvL,kconvT*prod->T());
  }
  return np;
}


void AliDecayerPythia8::Init()
{
// Initialisation
//
    if (!fgInit) {
	fgInit = kTRUE;
	// fPythia->SetDecayTable();
    }

// Switch on heavy flavor decays
    
    Int_t j;
    Int_t heavy[14] = {411, 421, 431, 4122, 4132, 4232, 4332, 511, 521, 531, 5122, 5132, 5232, 5332};
//    fPythia->ResetDecayTable();
    for (j=0; j < 14; j++) {
	if (fDecay == kNoDecayHeavy) {
	    fPythia8->ReadString(Form("%d:onMode = off", heavy[j]));
	} else {
	    fPythia8->ReadString(Form("%d:onMode = on", heavy[j]));
	}
    }
 
    fPythia8->ReadString("111:onMode = on");
    
//...Switch off decay of K0S, Lambda, Sigma+-, Xi0-, Omega-.
    fPythia8->ReadString("310:onMode  = off");
    fPythia8->ReadString("3122:onMode = off");
    fPythia8->ReadString("3112:onMode = off");
    fPythia8->ReadString("3222:onMode = off");
    fPythia8->ReadString("3312:onMode = off");
    fPythia8->ReadString("3322:onMode = off");
    fPythia8->ReadString("3334:onMode = off");
// .. Force decay channels
    ForceDecay();
}

void AliDecayerPythia8::ForceDecay()
{
// 
// Force a particle decay mode
// Switch heavy flavour production off if requested
    if (!fHeavyFlavour) SwitchOffHeavyFlavour();
//
    Decay_t decay = fDecay;
    fPythia8->ReadString("HadronLevel:Decay = on");
    
    if (decay == kNoDecayHeavy) return;
//
// select mode    
    switch (decay) 
    {
    case kHardMuons:
//	B0 -> mu X 
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 13 443 100443");
//	B+/- -> mu X 
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 13 443 100443");
//	Bs -> mu X 
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 13 443 100443");
//	Lambda_b -> mu X 
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 13 443 100443");
//      Sigma_b- -> mu X 
	fPythia8->ReadString("5132:onMode = off");
	fPythia8->ReadString("5132:onIfAny = 13 443 100443");
//      Sigma_b0 -> mu X
	fPythia8->ReadString("5232:onMode = off");
	fPythia8->ReadString("5232:onIfAny = 13 443 100443");
//      Omega_b  -> mu X
	fPythia8->ReadString("5332:onMode = off");
	fPythia8->ReadString("5332:onIfAny = 13 443 100443");
//      Psi' -> mu+ mu-
	fPythia8->ReadString("100443:onMode = off");
	fPythia8->ReadString("100443:onIfAny = 443");
//      Psi  -> mu+ mu-
	fPythia8->ReadString("443:onMode = off");
	fPythia8->ReadString("443:onIfAll = 13 13");
//      D+/- -> mu X
	fPythia8->ReadString("411:onMode = off");
	fPythia8->ReadString("411:onIfAll = 13");
//      D0   -> mu X
	fPythia8->ReadString("421:onMode = off");
	fPythia8->ReadString("421:onIfAll = 13");
//      D_s  -> mu X
	fPythia8->ReadString("431:onMode = off");
	fPythia8->ReadString("431:onIfAll = 13");
//      Lambda_c -> mu X
	fPythia8->ReadString("4122:onMode = off");
	fPythia8->ReadString("4122:onIfAll = 13");
//      Sigma_c  -> mu X
	fPythia8->ReadString("4132:onMode = off");
	fPythia8->ReadString("4132:onIfAll = 13");
//      Sigma_c+ -> mu X
	fPythia8->ReadString("4232:onMode = off");
	fPythia8->ReadString("4232:onIfAll = 13");
//      Omega_c  -> mu X
	fPythia8->ReadString("4332:onMode = off");
	fPythia8->ReadString("4332:onIfAll = 13");

	break;
   case kChiToJpsiGammaToMuonMuon:
// Chi_1c  -> J/Psi  Gamma
	fPythia8->ReadString("20443:onMode = off");
	fPythia8->ReadString("20443:onIfAll = 443 22");
// Chi_2c  -> J/Psi  Gamma
	fPythia8->ReadString("445:onMode = off");
	fPythia8->ReadString("445:onIfAll = 443 22");
// J/Psi -> mu+ mu-
	fPythia8->ReadString("443:onMode = off");
	fPythia8->ReadString("443:onIfAll = 13 13");
	break;
    case kChiToJpsiGammaToElectronElectron:
// Chi_1c  -> J/Psi  Gamma
	fPythia8->ReadString("20443:onMode = off");
	fPythia8->ReadString("20443:onIfAll = 443 22");
// Chi_2c  -> J/Psi  Gamma
	fPythia8->ReadString("445:onMode = off");
	fPythia8->ReadString("445:onIfAll = 443 22");
// J/Psi -> e+ e-
	fPythia8->ReadString("443:onMode = off");
	fPythia8->ReadString("443:onIfAll = 11 11");
	break;

    case kBSemiMuonic:
//      B0 -> mu X 
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 13");
//	B+/- -> mu X 
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 13");
//	B_s -> mu X 
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 13");
//	Lambda_b -> mu X 
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 13");
//	Sigma_b -> mu X 
	fPythia8->ReadString("5132:onMode = off");
	fPythia8->ReadString("5132:onIfAny = 13");
//	Sigma_b0 -> mu X 
	fPythia8->ReadString("5232:onMode = off");
	fPythia8->ReadString("5232:onIfAny = 13");
//	Omega_b  -> mu X 
	fPythia8->ReadString("5332:onMode = off");
	fPythia8->ReadString("5332:onIfAny = 13");
	break;
    case kDSemiMuonic:
//      D+- -> mu X
	fPythia8->ReadString("411:onMode = off");
	fPythia8->ReadString("411:onIfAll = 13");
//      D0  -> mu X
	fPythia8->ReadString("421:onMode = off");
	fPythia8->ReadString("421:onIfAll = 13");
//      D_s  -> mu X
	fPythia8->ReadString("431:onMode = off");
	fPythia8->ReadString("431:onIfAll = 13");
//      Lambda_c -> mu X
	fPythia8->ReadString("4122:onMode = off");
	fPythia8->ReadString("4122:onIfAll = 13");
//      Sigma_c  -> mu X
	fPythia8->ReadString("4132:onMode = off");
	fPythia8->ReadString("4132:onIfAll = 13");
//      Sigma  -> mu X
	fPythia8->ReadString("4232:onMode = off");
	fPythia8->ReadString("4232:onIfAll = 13");
//      Omega_c  -> mu X
	fPythia8->ReadString("4332:onMode = off");
	fPythia8->ReadString("4332:onIfAll = 13");
	break;
    case kSemiMuonic:
//      D+- -> mu X
	fPythia8->ReadString("411:onMode = off");
	fPythia8->ReadString("411:onIfAll = 13");
//      D0  -> mu X
	fPythia8->ReadString("421:onMode = off");
	fPythia8->ReadString("421:onIfAll = 13");
//      D_s  -> mu X
	fPythia8->ReadString("431:onMode = off");
	fPythia8->ReadString("431:onIfAll = 13");
//      Lambda_c -> mu X
	fPythia8->ReadString("4122:onMode = off");
	fPythia8->ReadString("4122:onIfAll = 13");
//      Sigma_c  -> mu X
	fPythia8->ReadString("4132:onMode = off");
	fPythia8->ReadString("4132:onIfAll = 13");
//      Sigma  -> mu X
	fPythia8->ReadString("4232:onMode = off");
	fPythia8->ReadString("4232:onIfAll = 13");
//      Omega_c  -> mu X
	fPythia8->ReadString("4332:onMode = off");
	fPythia8->ReadString("4332:onIfAll = 13");
//      B0       -> mu X
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 13");
//      B+/-     -> mu X
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 13");
//      B_s      -> mu X
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 13");
//      Lambda_c -> mu X
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 13");
//      Sigma_c  -> mu X
	fPythia8->ReadString("5132:onMode = off");
	fPythia8->ReadString("5132:onIfAny = 13");
//      Sigma_c  -> mu X
	fPythia8->ReadString("5232:onMode = off");
	fPythia8->ReadString("5232:onIfAny = 13");
//      Omega_c  -> mu X
	fPythia8->ReadString("5332:onMode = off");
	fPythia8->ReadString("5332:onIfAny = 13");

	break;
    case kJpsiDiMuon:
//      J/Psi-> mu+ mu-
	fPythia8->ReadString("443:onMode = off");
	fPythia8->ReadString("443:onIfAll = 13 13");
	break;
    case kDiMuon:
//      Rho -> mu+ mu-
	fPythia8->ReadString("113:onMode = off");
	fPythia8->ReadString("113:onIfAll = 13 13");
//      Eta-> mu+ mu-
	fPythia8->ReadString("221:onMode = off");
	fPythia8->ReadString("221:onIfAll = 13 13");
//      omega-> mu+ mu-
	fPythia8->ReadString("223:onMode = off");
	fPythia8->ReadString("223:onIfAll = 13 13");
//      phi-> mu+ mu-
	fPythia8->ReadString("333:onMode = off");
	fPythia8->ReadString("333:onIfAll = 13 13");
//      J/Psi-> mu+ mu-
	fPythia8->ReadString("443:onMode = off");
	fPythia8->ReadString("443:onIfAll = 13 13");
//      Psi'-> mu+ mu-
	fPythia8->ReadString("100443:onMode = off");
	fPythia8->ReadString("100443:onIfAll = 13 13");
//      Ups-> mu+ mu-
	fPythia8->ReadString("553:onMode = off");
	fPythia8->ReadString("553:onIfAll = 13 13");
//      Ups'-> mu+ mu-
	fPythia8->ReadString("100553:onMode = off");
	fPythia8->ReadString("100553:onIfAll = 13 13"); 
//      Ups''-> mu+ mu-
	fPythia8->ReadString("200553:onMode = off");
	fPythia8->ReadString("200553:onIfAll = 13 13");
	break;
    case kBSemiElectronic:
//      B0 - > e+ e-
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 11");
//      B+- -> e+ e-
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 11");
//      B_s -> e+ e-
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 11");
//      Lambda_b -> e+ e-
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 11");
//      Sigma_b -> e+ e-
	fPythia8->ReadString("5132:onMode = off");
	fPythia8->ReadString("5132:onIfAny = 11");
//      Sigma_b -> e+ e-
	fPythia8->ReadString("5232:onMode = off");
	fPythia8->ReadString("5232:onIfAny = 11");
//      Omega_b ->e+ e-
	fPythia8->ReadString("5332:onMode = off");
	fPythia8->ReadString("5332:onIfAny = 11");
	break;
    case kSemiElectronic:
//      D+/- -> e X
	fPythia8->ReadString("411:onMode = off");
	fPythia8->ReadString("411:onIfAll = 11");
//      D0 -> e X
	fPythia8->ReadString("421:onMode = off");
	fPythia8->ReadString("421:onIfAll = 11");
//      D_s ->e X
	fPythia8->ReadString("431:onMode = off");
	fPythia8->ReadString("431:onIfAll = 11");
//      Lambda_c -> e X
	fPythia8->ReadString("4122:onMode = off");
	fPythia8->ReadString("4122:onIfAll = 11");
//      Sigma_c -> e X
	fPythia8->ReadString("4132:onMode = off");
	fPythia8->ReadString("4132:onIfAll = 11");
//      Sigma_c -> e X
	fPythia8->ReadString("4232:onMode = off");
	fPythia8->ReadString("4232:onIfAll = 11");
//      Omega_c -> e X
	fPythia8->ReadString("4332:onMode = off");
	fPythia8->ReadString("4332:onIfAll = 11");
//      B0 -> e X
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 11");
//      B+/- -> e X
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 11");
//      B_s -> e X
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 11");
//      Lambda_b -> e X
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 11");
//      Sigma_b -> e X
	fPythia8->ReadString("5132:onMode = off");
	fPythia8->ReadString("5132:onIfAny = 11");
//      Sigma_b -> e X
	fPythia8->ReadString("5232:onMode = off");
	fPythia8->ReadString("5232:onIfAny = 11");
//      Omega_b -> e X
	fPythia8->ReadString("5332:onMode = off");
	fPythia8->ReadString("5332:onIfAny = 11");
	break;
    case kDiElectron:
//      Rho -> e+e-
	fPythia8->ReadString("113:onMode = off");
	fPythia8->ReadString("113:onIfAll = 11 11");
//      Eta -> e+e-
	fPythia8->ReadString("221:onMode = off");
	fPythia8->ReadString("221:onIfAll = 11 11");
//      omega -> e+e-
	fPythia8->ReadString("223:onMode = off");
	fPythia8->ReadString("223:onIfAll = 11 11");
//      phi -> e+e-
	fPythia8->ReadString("333:onMode = off");
	fPythia8->ReadString("333:onIfAll = 11 11");
//      J/Psi -> e+e-
	fPythia8->ReadString("443:onMode = off");
	fPythia8->ReadString("443:onIfAll = 11 11");
//      Psi' -> e+e-
	fPythia8->ReadString("100443:onMode = off");
	fPythia8->ReadString("100443:onIfAll = 11 11");
//      Ups -> e+e-
	fPythia8->ReadString("553:onMode = off");
	fPythia8->ReadString("553:onIfAll = 11 11");
//      Ups' -> e+e-
	fPythia8->ReadString("100553:onMode = off");
	fPythia8->ReadString("100553:onIfAll = 11 11");
//      Ups'' -> e+e-
	fPythia8->ReadString("200553:onMode = off");
	fPythia8->ReadString("200553:onIfAll = 11 11");
	break;
    case kBJpsiDiMuon:
//      B0   -> J/Psi (Psi') X   
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 443 100443");
//      B+/-   -> J/Psi (Psi') X   
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 443 100443");
//      B_s   -> J/Psi (Psi') X   
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 443 100443");
//      Lambda_b -> J/Psi (Psi') X   
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 443 100443");
//
//      J/Psi -> mu+ mu-
	fPythia8->ReadString("443:onMode = off");
	fPythia8->ReadString("443:onIfAll = 13 13");
//      Psi' -> mu+ mu-
	fPythia8->ReadString("100443:onMode = off");
	fPythia8->ReadString("100443:onIfAll = 13 13");
	break;
    case kBPsiPrimeDiMuon:
//      B0   -> Psi' X   
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 100443");
//      B+/-   -> Psi' X   
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 100443");
//      B_s   -> Psi'  X   
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 100443");
//      Lambda_b -> Psi' X   
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 100443");
//
//      Psi' -> mu+ mu-
	fPythia8->ReadString("100443:onMode = off");
	fPythia8->ReadString("100443:onIfAll = 13 13");
	break;
    case kBJpsiDiElectron:
//      B0   -> Psi X   
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 443");
//      B+/-   -> Psi X   
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 443");
//      B_s   -> Psi  X   
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 443");
//      Lambda_b -> Psi X   
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 443");
//
//      Psi -> mu+ mu-
	fPythia8->ReadString("443:onMode = off");
	fPythia8->ReadString("443:onIfAll = 11 11");

	break;
    case kBJpsi:
//      B0   -> Psi X   
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 443");
//      B+/-   -> Psi X   
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 443");
//      B_s   -> Psi  X   
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 443");
//      Lambda_b -> Psi X   
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 443");
	break;
    case kBPsiPrimeDiElectron:
//      B0   -> Psi' X   
	fPythia8->ReadString("511:onMode = off");
	fPythia8->ReadString("511:onIfAny = 100443");
//      B+/-   -> Psi' X   
	fPythia8->ReadString("521:onMode = off");
	fPythia8->ReadString("521:onIfAny = 100443");
//      B_s   -> Psi'  X   
	fPythia8->ReadString("531:onMode = off");
	fPythia8->ReadString("531:onIfAny = 100443");
//      Lambda_b -> Psi' X   
	fPythia8->ReadString("5122:onMode = off");
	fPythia8->ReadString("5122:onIfAny = 100443");
//
//      Psi' -> mu+ mu-
	fPythia8->ReadString("100443:onMode = off");
	fPythia8->ReadString("100443:onIfAll = 11 11");
	break;
    case kPiToMu:
//      pi -> mu nu
	fPythia8->ReadString("211:onMode = off");
	fPythia8->ReadString("211:onIfAny = 13");
	break;
    case kKaToMu:
//      K -> mu nu
	fPythia8->ReadString("321:onMode = off");
	fPythia8->ReadString("321:onIfAny = 13");
	break;
    case kAllMuonic:
//      pi/K -> mu
	fPythia8->ReadString("211:onMode = off");
	fPythia8->ReadString("211:onIfAny = 13");
	fPythia8->ReadString("321:onMode = off");
	fPythia8->ReadString("321:onIfAny = 13");
	break;
    case kWToMuon:
//      W -> mu X
	fPythia8->ReadString("24:onMode = off");
	fPythia8->ReadString("24:onIfAny = 13");
	break;
    case kWToCharm:
//      W -> c X
	fPythia8->ReadString("24:onMode = off");
	fPythia8->ReadString("24:onIfAny = 4");
	break;
    case kWToCharmToMuon:
//      W -> c X
	fPythia8->ReadString("24:onMode = off");
	fPythia8->ReadString("24:onIfAny = 4");
//      D+- -> mu X
	fPythia8->ReadString("411:onMode = off");
	fPythia8->ReadString("411:onIfAll = 13");
//      D0 -> mu X
	fPythia8->ReadString("421:onMode = off");
	fPythia8->ReadString("421:onIfAll = 13");
//      D_s -> mu X
	fPythia8->ReadString("431:onMode = off");
	fPythia8->ReadString("431:onIfAll = 13");
//      Lambda_c -> mu X
	fPythia8->ReadString("4122:onMode = off");
	fPythia8->ReadString("4122:onIfAll = 13");
//      Sigma_c -> mu X
	fPythia8->ReadString("4132:onMode = off");
	fPythia8->ReadString("4132:onIfAll = 13");
//      Sigma_c -> mu X
	fPythia8->ReadString("4232:onMode = off");
	fPythia8->ReadString("4232:onIfAll = 13");
//      Omega_c -> mu X
	fPythia8->ReadString("4332:onMode = off");
	fPythia8->ReadString("4332:onIfAll = 13");
	break;
    case kZDiMuon:
//      Z -> mu+ mu-
	fPythia8->ReadString("23:onMode = off");
	fPythia8->ReadString("23:onIfAll = 13 13");
	break;
    case kZDiElectron:
//      Z -> e+ e-
	fPythia8->ReadString("23:onMode = off");
	fPythia8->ReadString("23:onIfAll = 11 11");
	break;
    case kHadronicD:
      ForceHadronicD(1,0,0);
	break;
    case kHadronicDWithV0:
        ForceHadronicD(1,1,0);
        break;
    case kHadronicDWithout4Bodies:
        ForceHadronicD(0,0,0);
        break;
    case kHadronicDWithout4BodiesWithV0:
        ForceHadronicD(0,1,0);
        break;
    case kLcpKpi:
      ForceHadronicD(0,0,1);  // Lc -> p K pi
      break;
    case kLcpK0S:
      ForceHadronicD(0,0,2);  // Lc -> p K0S
      break;
    case kPhiKK:
	// Phi-> K+ K-
	fPythia8->ReadString("333:onMode = off");
	fPythia8->ReadString("333:onIfAll = 321 321");
	break;
    case kOmega:
	// Omega -> Lambda K
	fPythia8->ReadString("3334:onMode = off");
	fPythia8->ReadString("3334:onIfAll = 3122 321 ");
	break;
    case kLambda:
	// Lambda -> p pi-
	fPythia8->ReadString("3122:onMode = off");
	fPythia8->ReadString("3122:onIfAll = 2212 211 ");
	break;
    case kBeautyUpgrade:
      ForceBeautyUpgrade();
	break;
    case kAll:
	break;
    case kNoDecay:
	fPythia8->ReadString("HadronLevel:Decay = off");
	break;
    case kNoDecayHeavy:
    case kNoDecayBeauty:
    case kNeutralPion:
    case kPsiPrimeJpsiDiElectron: 
    case kElectronEM:
    case kGammaEM:
    case kAllEM:
    case kDiElectronEM:
	break;
    }
}

Float_t AliDecayerPythia8::GetPartialBranchingRatio(Int_t ipart)
{
    // Get the partial branching ration for the forced decay channels
    
    Pythia8::Pythia* thePythia       = fPythia8->Pythia8();
    Pythia8::ParticleData & table    = thePythia->particleData;
    Pythia8::ParticleDataEntry* pd   = table.particleDataEntryPtr(ipart);

    Int_t nc = pd->sizeChannels();
    Float_t br = 0.;
//
//  Loop over decay channels
    for (Int_t ic = 0; ic < nc; ic++) {
      Pythia8::DecayChannel& decCh = pd->channel(ic);
	for (Int_t i = 0; i < decCh.multiplicity(); i++) {
	    br += decCh.bRatio();
	}
    }
    return (br);
}


Float_t AliDecayerPythia8::GetLifetime(Int_t kf)
{
    // Return lifetime of particle
    Pythia8::Pythia* thePythia       = fPythia8->Pythia8();
    Pythia8::ParticleData& table     = thePythia->particleData;
    Float_t tau = table.tau0(kf);
    return (tau);
}

void  AliDecayerPythia8::SwitchOffHeavyFlavour()
{
    // Switch off heavy flavour production
    //
// Maximum number of quark flavours used in pdf 
    fPythia8->ReadString("PDFinProcess:nQuarkIn = 3");
// Maximum number of flavors that can be used in showers
    fPythia8->ReadString("SpaceShower:nQuarkIn = 3");
    fPythia8->ReadString("TimeShower:nGammaToQuark = 3");
    fPythia8->ReadString("TimeShower:nGluonToQuark = 3");
}

void AliDecayerPythia8::ForceBeautyUpgrade(){
  // what about particle signs?

  // Bs -> Ds- pi+
  fPythia8->ReadString("531:onMode = off");
  fPythia8->ReadString("531:onIfMatch = 431 211");

  //Lb: 50% to Lc any, 50% to Lc pion
  fPythia8->ReadString("5122:onMode = off");
  if(gRandom->Rndm()<0.50) {
    fPythia8->ReadString("5122:onIfAll = 4122");
  }
  else {
    fPythia8->ReadString("5122:onIfMatch = 4122 211");
  }
  
  // B0 -> D*-pi+
  fPythia8->ReadString("511:onMode = off");
  fPythia8->ReadString("511:onIfMatch = 413 211");

  // B+ -> D0bar pi-
  fPythia8->ReadString("521:onMode = off");
  fPythia8->ReadString("521:onIfMatch = 421 211");

  ForceHadronicD(0,0,0);

}


void AliDecayerPythia8::ForceHadronicD(Int_t optUse4Bodies, Int_t optUseDtoV0, Int_t optForceLcChannel)
{
//
// Force golden D decay modes
//

    //add D+ decays absent in PYTHIA8 decay table and set BRs from PDG for other
    fPythia8->ReadString("411:oneChannel = 1 0.0752 0 -321 211 211");
    fPythia8->ReadString("411:addChannel = 1 0.0104 0 -313 211");
    fPythia8->ReadString("411:addChannel = 1 0.0156 0 311 211");
    //add Lc decays absent in PYTHIA8 decay table and set BRs from PDG for other
    fPythia8->ReadString("4122:oneChannel = 1 0.0196 100 2212 -313");
    fPythia8->ReadString("4122:addChannel = 1 0.0108 100 2224 -321");
    fPythia8->ReadString("4122:addChannel = 1 0.022 100 3124 211");
    fPythia8->ReadString("4122:addChannel = 1 0.035 0 2212 -321 211");
    fPythia8->ReadString("4122:addChannel = 1 0.0159 0 2212 311");
    fPythia8->ReadString("4122:addChannel = 1 0.0130 0 3122 211");
    //add Xic+ decays absent in PYTHIA8 decay table
    fPythia8->ReadString("4232:addChannel = 1 0.2 0 2212 313");
    fPythia8->ReadString("4232:addChannel = 1 0.2 0 2212 321 211");
    fPythia8->ReadString("4232:addChannel = 1 0.2 0 3324 211");
    fPythia8->ReadString("4232:addChannel = 1 0.2 0 3312 211 211");
    //add Xic0 decays absent in PYTHIA8 decay table
    fPythia8->ReadString("4132:addChannel = 1 0.2 0 3312 211");

    // K* -> K pi
    fPythia8->ReadString("313:onMode = off");
    fPythia8->ReadString("313:onIfAll = 321 211");
    // for Ds -> Phi pi+
    fPythia8->ReadString("333:onMode = off");
    fPythia8->ReadString("333:onIfAll = 321 321");
    // for D0 -> rho0 pi+ k-
    fPythia8->ReadString("113:onMode = off");
    fPythia8->ReadString("113:onIfAll = 211 211");
    // for Lambda_c -> Delta++ K-
    fPythia8->ReadString("2224:onMode = off");
    fPythia8->ReadString("2224:onIfAll = 2212 211");
    // for Lambda_c -> Lambda(1520) K-
    fPythia8->ReadString("3124:onMode = off");
    fPythia8->ReadString("3124:onIfAll = 2212 321");


    // Omega_c->Omega pi
    fPythia8->ReadString("4332:onMode = off");
    fPythia8->ReadString("4332:onIfMatch = 3334 211");


    fPythia8->ReadString("411:onMode = off");
    fPythia8->ReadString("421:onMode = off");
    fPythia8->ReadString("431:onMode = off");
    fPythia8->ReadString("4112:onMode = off");
    fPythia8->ReadString("4122:onMode = off");
    fPythia8->ReadString("4232:onMode = off");
    fPythia8->ReadString("4132:onMode = off");


    // D+/- -> K pi pi 
    fPythia8->ReadString("411:onIfMatch = 321 211 211");
    // D+/- -> K* pi
    fPythia8->ReadString("411:onIfMatch = 313 211");
    // D0 -> K pi
    fPythia8->ReadString("421:onIfMatch = 321 211");

    if (optUse4Bodies) {
      // D0 -> K pi pi pi
      fPythia8->ReadString("421:onIfMatch = 321 211 211 211");
      // D0 -> K pi rho
      fPythia8->ReadString("421:onIfMatch = 321 211 113");
      // D0 -> K*0 pi pi
      fPythia8->ReadString("421:onIfMatch = 313 211 211");
    }
    if(optUseDtoV0==1){
      // Ds -> K0K
      fPythia8->ReadString("431:onIfMatch = 311 321");
      // Ds -> K0pi
      fPythia8->ReadString("411:onIfMatch = 311 211");
    }


    // D_s -> K K*
    fPythia8->ReadString("431:onIfMatch = 321 313");
    // D_s -> Phi pi
    fPythia8->ReadString("431:onIfMatch = 333 211");

    if (optForceLcChannel == 0) { // all Lc channels
    // Lambda_c -> p K*
    fPythia8->ReadString("4122:onIfMatch = 2212 313");
    // Lambda_c -> Delta K
    fPythia8->ReadString("4122:onIfMatch = 2224 321");
    // Lambda_c -> Lambda(1520) pi
    fPythia8->ReadString("4122:onIfMatch = 3124 211");
    // Lambda_c -> p K pi
    fPythia8->ReadString("4122:onIfMatch = 2212 321 211");
    // Lambda_c -> Lambda pi 
    fPythia8->ReadString("4122:onIfMatch = 3122 211");
    // Lambda_c -> p K0
    fPythia8->ReadString("4122:onIfMatch = 2212 311");
    }

    if (optForceLcChannel == 1) { // force only Lc -> p K pi
    fPythia8->ReadString("4122:onIfMatch = 2212 321 211"); 
    }

    if (optForceLcChannel == 2) { // force only Lc -> p K0S
    fPythia8->ReadString("4122:onIfMatch = 2212 311");
    }

    // Xic+ -> pK*0
    fPythia8->ReadString("4232:onIfMatch = 2212 313");
    // Xic+ -> p K- pi+
    fPythia8->ReadString("4232:onIfMatch = 2212 321 211");
    // Xic+ -> Xi*0 pi+, Xi*->Xi- pi+
    fPythia8->ReadString("4232:onIfMatch = 3324 211");
    // Xic+ -> Xi- pi+ pi+
    fPythia8->ReadString("4232:onIfMatch = 3312 211 211");
    
    // Xic0 -> Xi- pi+
    fPythia8->ReadString("4132:onIfMatch = 3312 211");
       
    
}

//___________________________________________________________________________
void    AliDecayerPythia8::ReadDecayTable()
{
   //to read a decay table (not yet implemented)
}


//___________________________________________________________________________
void AliDecayerPythia8::AppendParticle(Int_t pdg, TLorentzVector* p)
{
   // Append a particle to the stack
   fPythia8->Pythia8()->event.append(pdg, 11, 0, 0, p->Px(), p->Py(), p->Pz(), p->E(), p->M());   
}


//___________________________________________________________________________
void AliDecayerPythia8::ClearEvent()
{
   // Clear the event stack
   fPythia8->Pythia8()->event.clear();
}

