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
Revision 1.25  2000/10/18 19:11:27  hristov
Division by zero fixed

Revision 1.24  2000/09/18 10:41:35  morsch
Add possibility to use nuclear structure functions from PDF library V8.

Revision 1.23  2000/09/14 14:05:40  morsch
dito

Revision 1.22  2000/09/14 14:02:22  morsch
- Correct conversion from mm to cm when passing particle vertex to MC.
- Correct handling of fForceDecay == all.

Revision 1.21  2000/09/12 14:14:55  morsch
Call fDecayer->ForceDecay() at the beginning of Generate().

Revision 1.20  2000/09/06 14:29:33  morsch
Use AliPythia for event generation an AliDecayPythia for decays.
Correct handling of "nodecay" option

Revision 1.19  2000/07/11 18:24:56  fca
Coding convention corrections + few minor bug fixes

Revision 1.18  2000/06/30 12:40:34  morsch
Pythia takes care of vertex smearing. Correct conversion from Pythia units (mm) to
Geant units (cm).

Revision 1.17  2000/06/09 20:34:07  morsch
All coding rule violations except RS3 corrected

Revision 1.16  2000/05/15 15:04:20  morsch
The full event is written for fNtrack = -1
Coding rule violations corrected.

Revision 1.15  2000/04/26 10:14:24  morsch
Particles array has one entry more than pythia particle list. Upper bound of
particle loop changed to np-1 (R. Guernane, AM)

Revision 1.14  2000/04/05 08:36:13  morsch
Check status code of particles in Pythia event
to avoid double counting as partonic state and final state particle.

Revision 1.13  1999/11/09 07:38:48  fca
Changes for compatibility with version 2.23 of ROOT

Revision 1.12  1999/11/03 17:43:20  fca
New version from G.Martinez & A.Morsch

Revision 1.11  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log
*/

#include "AliGenPythia.h"
#include "AliDecayerPythia.h"
#include "AliRun.h"
#include "AliPythia.h"
#include "AliPDG.h"
#include <TParticle.h>
#include <TSystem.h>

 ClassImp(AliGenPythia)

AliGenPythia::AliGenPythia()
                 :AliGenerator()
{
// Default Constructor
  fDecayer = new AliDecayerPythia();
}

AliGenPythia::AliGenPythia(Int_t npart)
                 :AliGenerator(npart)
{
// default charm production at 5. 5 TeV
// semimuonic decay
// structure function GRVHO
//
    fXsection  = 0.;
    fNucA1=0;
    fNucA2=0;
    fParentSelect.Set(5);
    fChildSelect.Set(5);
    for (Int_t i=0; i<5; i++) fParentSelect[i]=fChildSelect[i]=0;
    SetProcess();
    SetStrucFunc();
    SetForceDecay();
    SetPtHard();
    SetEnergyCMS();
    fDecayer = new AliDecayerPythia();
}

AliGenPythia::AliGenPythia(const AliGenPythia & Pythia)
{
// copy constructor
}

AliGenPythia::~AliGenPythia()
{
// Destructor
}

void AliGenPythia::Init()
{
// Initialisation
  SetMC(AliPythia::Instance());
    fPythia=(AliPythia*) fgMCEvGen;
//
    fParentWeight=1./Float_t(fNpart);
//
//  Forward Paramters to the AliPythia object
    //    gSystem->Exec("ln -s $ALICE_ROOT/data/Decay.table fort.1");
    //    fPythia->Pyupda(2,1);    
    //    gSystem->Exec("rm fort.1");

    fDecayer->SetForceDecay(fForceDecay);    
    fDecayer->Init();


    fPythia->SetCKIN(3,fPtHardMin);
    fPythia->SetCKIN(4,fPtHardMax);    
    if (fNucA1 > 0 && fNucA2 > 0) fPythia->SetNuclei(fNucA1, fNucA2);  
    fPythia->ProcInit(fProcess,fEnergyCMS,fStrucFunc);

    //    fPythia->Pylist(0);
    //    fPythia->Pystat(2);
//  Parent and Children Selection
    switch (fProcess) 
    {
    case charm:

	fParentSelect[0]=411;
	fParentSelect[1]=421;
	fParentSelect[2]=431;
	fParentSelect[3]=4122;	
	break;
    case charm_unforced:

	fParentSelect[0]=411;
	fParentSelect[1]=421;
	fParentSelect[2]=431;
	fParentSelect[3]=4122;	
	break;
    case beauty:
	fParentSelect[0]=511;
	fParentSelect[1]=521;
	fParentSelect[2]=531;
	fParentSelect[3]=5122;	
	break;
    case beauty_unforced:
	fParentSelect[0]=511;
	fParentSelect[1]=521;
	fParentSelect[2]=531;
	fParentSelect[3]=5122;	
	break;
    case jpsi_chi:
    case jpsi:
	fParentSelect[0]=443;
	break;
    case mb:
	break;
    }

    switch (fForceDecay) 
    {
    case semielectronic:
    case dielectron:
    case b_jpsi_dielectron:
    case b_psip_dielectron:
	fChildSelect[0]=kElectron;	
	break;
    case semimuonic:
    case dimuon:
    case b_jpsi_dimuon:
    case b_psip_dimuon:
    case pitomu:
    case katomu:
	fChildSelect[0]=kMuonMinus;
	break;
    case hadronicD:
      fChildSelect[0]=kPiPlus;
      fChildSelect[1]=kKPlus;
      break;
    case all:
    case nodecay:
      break;
    }
}

void AliGenPythia::Generate()
{
// Generate one event
    fDecayer->ForceDecay();

    Float_t polar[3] =   {0,0,0};
    Float_t origin[3]=   {0,0,0};
    Float_t originP[3]= {0,0,0};
    Float_t origin0[3]=  {0,0,0};
    Float_t p[3], pP[4];
//    Float_t random[6];
    static TClonesArray *particles;
//  converts from mm/c to s
    const Float_t kconv=0.001/2.999792458e8;
    
    
//
    Int_t nt=0;
    Int_t ntP=0;
    Int_t jev=0;
    Int_t j, kf;

    if(!particles) particles=new TClonesArray("TParticle",1000);
    
    fTrials=0;
    for (j=0;j<3;j++) origin0[j]=fOrigin[j];
    if(fVertexSmear==kPerEvent) {
	fPythia->SetMSTP(151,1);
	for (j=0;j<3;j++) {
	    fPythia->SetPARP(151+j, fOsigma[j]/10.);
	}
    } else if (fVertexSmear==kPerTrack) {
	fPythia->SetMSTP(151,0);
    }
    
    while(1)
    {
	fPythia->Pyevnt();
//	fPythia->Pylist(1);
	fTrials++;
	fPythia->ImportParticles(particles,"All");
	Int_t np = particles->GetEntriesFast();
	printf("\n **************************************************%d\n",np);
	Int_t nc=0;
	if (np == 0 ) continue;
	if (fProcess != mb) {
	    for (Int_t i = 0; i<np-1; i++) {
		TParticle *  iparticle = (TParticle *) particles->At(i);
		Int_t ks = iparticle->GetStatusCode();
		kf = CheckPDGCode(iparticle->GetPdgCode());
		if (ks==21) continue;

		fChildWeight=(fDecayer->GetPartialBranchingRatio(kf))*fParentWeight;	  
//
// Parent
		if (ParentSelected(TMath::Abs(kf))) {
		    if (KinematicSelection(iparticle)) {
			if (nc==0) {
//
// Store information concerning the hard scattering process
//
			    Float_t massP  = fPythia->GetPARI(13);
			    Float_t   ptP  = fPythia->GetPARI(17);
			    Float_t    yP  = fPythia->GetPARI(37);
			    Float_t  xmtP  = sqrt(ptP*ptP+massP*massP);
			    Float_t    ty  = Float_t(TMath::TanH(yP));
			    pP[0] = ptP;
			    pP[1] = 0;
			    pP[2] = xmtP*ty/sqrt(1.-ty*ty);
			    pP[3] = massP;
			    gAlice->SetTrack(0,-1,-1,
					     pP,originP,polar,
					     0,kPPrimary,ntP,fParentWeight);
//					     0,"Hard Scat.",ntP,fParentWeight);
			    gAlice->KeepTrack(ntP);
			}
			nc++;
//
// store parent track information
			p[0]=iparticle->Px();
			p[1]=iparticle->Py();
			p[2]=iparticle->Pz();
			origin[0]=origin0[0]+iparticle->Vx()/10.;
			origin[1]=origin0[1]+iparticle->Vy()/10.;
			origin[2]=origin0[2]+iparticle->Vz()/10.;

			Int_t ifch=iparticle->GetFirstDaughter();
			Int_t ilch=iparticle->GetLastDaughter();	

			if ((ifch !=0 && ilch !=0) || fForceDecay == nodecay) {
			  Int_t trackit=0;
			  if (fForceDecay == nodecay) trackit = 1;
			    gAlice->SetTrack(trackit,ntP,kf,
					     p,origin,polar,
					     0,kPPrimary,nt,fParentWeight);
			    gAlice->KeepTrack(nt);
			    Int_t iparent = nt;
//
// Children	    
			    if (fForceDecay != nodecay) {
			      for (j=ifch; j<=ilch; j++)
				{
				  TParticle *  ichild = 
				    (TParticle *) particles->At(j-1);
				  kf = CheckPDGCode(ichild->GetPdgCode());
//
// 
				  if (ChildSelected(TMath::Abs(kf))) {
				    origin[0]=origin0[0]+ichild->Vx()/10.;
				    origin[1]=origin0[1]+ichild->Vy()/10.;
				    origin[2]=origin0[2]+ichild->Vz()/10.;		
				    p[0]=ichild->Px();
				    p[1]=ichild->Py();
				    p[2]=ichild->Pz();
				    Float_t tof=kconv*ichild->T();
				    gAlice->SetTrack(fTrackIt, iparent, kf,
						     p,origin,polar,
						     tof,kPDecay,nt,fChildWeight);
				    gAlice->KeepTrack(nt);
				  } // select child
				} // child loop
			    } 
			}
		    } // kinematic selection
		} // select particle
	    } // particle loop
	} else {
	    for (Int_t i = 0; i<np-1; i++) {
		TParticle *  iparticle = (TParticle *) particles->At(i);
		kf = CheckPDGCode(iparticle->GetPdgCode());
		Int_t ks = iparticle->GetStatusCode();

		if (ks==1 && kf!=0 && KinematicSelection(iparticle)) {
			nc++;
//
// store track information
			p[0]=iparticle->Px();
			p[1]=iparticle->Py();
			p[2]=iparticle->Pz();
			origin[0]=origin0[0]+iparticle->Vx()/10.;
			origin[1]=origin0[1]+iparticle->Vy()/10.;
			origin[2]=origin0[2]+iparticle->Vz()/10.;
			Float_t tof=kconv*iparticle->T();
			gAlice->SetTrack(fTrackIt,-1,kf,p,origin,polar,
					 tof,kPPrimary,nt);
			gAlice->KeepTrack(nt);
		} // select particle
	    } // particle loop 
	    printf("\n I've put %i particles on the stack \n",nc);
	} // mb ?
	if (nc > 0) {
	    jev+=nc;
	    if (jev >= fNpart || fNpart == -1) {
		fKineBias=Float_t(fNpart)/Float_t(fTrials);
		printf("\n Trials: %i %i %i\n",fTrials, fNpart, jev);
// Print x-section summary
		fPythia->Pystat(1);
		break;
	    }
	}
    } // event loop
//  adjust weight due to kinematic selection
    AdjustWeights();
//  get cross-section
    fXsection=fPythia->GetPARI(1);
}

Bool_t AliGenPythia::ParentSelected(Int_t ip)
{
// True if particle is in list of parent particles to be selected
    for (Int_t i=0; i<5; i++)
    {
	if (fParentSelect[i]==ip) return kTRUE;
    }
    return kFALSE;
}

Bool_t AliGenPythia::ChildSelected(Int_t ip)
{
// True if particle is in list of decay products to be selected
    if (fForceDecay == all) return kTRUE;
    
    for (Int_t i=0; i<5; i++)
    {
	if (fChildSelect[i]==ip) return kTRUE;
    }
    return kFALSE;
}

Bool_t AliGenPythia::KinematicSelection(TParticle *particle)
{
// Perform kinematic selection
    Float_t px=particle->Px();
    Float_t py=particle->Py();
    Float_t pz=particle->Pz();
    Float_t  e=particle->Energy();

//
//  transverse momentum cut    
    Float_t pt=TMath::Sqrt(px*px+py*py);
    if (pt > fPtMax || pt < fPtMin) 
    {
//	printf("\n failed pt cut %f %f %f \n",pt,fPtMin,fPtMax);
	return kFALSE;
    }
//
// momentum cut
    Float_t p=TMath::Sqrt(px*px+py*py+pz*pz);
    if (p > fPMax || p < fPMin) 
    {
//	printf("\n failed p cut %f %f %f \n",p,fPMin,fPMax);
	return kFALSE;
    }
    
//
// theta cut
    Float_t  theta = Float_t(TMath::ATan2(Double_t(pt),Double_t(pz)));
    if (theta > fThetaMax || theta < fThetaMin) 
    {
//	printf("\n failed theta cut %f %f %f \n",theta,fThetaMin,fThetaMax);
	return kFALSE;
    }

//
// rapidity cut
    if (e==pz) {
      return kFALSE;
    }
    else {
      Float_t y = 0.5*TMath::Log((e+pz)/(e-pz));
      if (y > fYMax || y < fYMin)
        {
//	printf("\n failed y cut %f %f %f \n",y,fYMin,fYMax);
          return kFALSE;
        }
    }

//
// phi cut
    Float_t phi=Float_t(TMath::ATan2(Double_t(py),Double_t(px)));
    if (phi > fPhiMax || phi < fPhiMin)
    {
//	printf("\n failed phi cut %f %f %f \n",phi,fPhiMin,fPhiMax);
	return kFALSE;
    }

    return kTRUE;
}
void AliGenPythia::AdjustWeights()
{
// Adjust the weights after generation of all events
//
    TClonesArray *partArray = gAlice->Particles();
    TParticle *part;
    Int_t ntrack=gAlice->GetNtrack();
    for (Int_t i=0; i<ntrack; i++) {
	part= (TParticle*) partArray->UncheckedAt(i);
	part->SetWeight(part->GetWeight()*fKineBias);
    }
}

Int_t AliGenPythia::CheckPDGCode(Int_t pdgcode)
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

    
void AliGenPythia::SetNuclei(Int_t a1, Int_t a2)
{
// Treat protons as inside nuclei with mass numbers a1 and a2  
    fNucA1 = a1;
    fNucA2 = a2;
}
	
	  
AliGenPythia& AliGenPythia::operator=(const  AliGenPythia& rhs)
{
// Assignment operator
    return *this;
}



void AliGenPythia::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliGenPythia.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliGenerator::Streamer(R__b);
      R__b >> (Int_t&)fProcess;
      R__b >> (Int_t&)fStrucFunc;
      R__b >> (Int_t&)fForceDecay;
      R__b >> fEnergyCMS;
      R__b >> fKineBias;
      R__b >> fTrials;
      fParentSelect.Streamer(R__b);
      fChildSelect.Streamer(R__b);
      R__b >> fXsection;
//      (AliPythia::Instance())->Streamer(R__b);
      R__b >> fPtHardMin;
      R__b >> fPtHardMax;
//      if (fDecayer) fDecayer->Streamer(R__b);
   } else {
      R__b.WriteVersion(AliGenPythia::IsA());
      AliGenerator::Streamer(R__b);
      R__b << (Int_t)fProcess;
      R__b << (Int_t)fStrucFunc;
      R__b << (Int_t)fForceDecay;
      R__b << fEnergyCMS;
      R__b << fKineBias;
      R__b << fTrials;
      fParentSelect.Streamer(R__b);
      fChildSelect.Streamer(R__b);
      R__b << fXsection;
//      R__b << fPythia;
      R__b << fPtHardMin;
      R__b << fPtHardMax;
      //     fDecayer->Streamer(R__b);
   }
}



#ifndef WIN32
#define pyr    pyr_
#define pyrset pyrset_
#define pyrget pyrget_
#else
#define pyr    PYR
#define pyrset PYRSET
#define pyrget PYRGET
#endif

extern "C" {
  Double_t pyr(Int_t*) {return sRandom->Rndm();}
  void pyrset(Int_t*,Int_t*) {}
  void pyrget(Int_t*,Int_t*) {}
}
