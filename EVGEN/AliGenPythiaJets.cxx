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
Revision 1.1  2003/01/17 04:10:31  morsch
First commit.
*/

//
// Generator using the TPythia interface (via AliPythia)
// to generate jets in pp collisions.
// Using SetNuclei() also nuclear modifications to the structure functions
// can be taken into account. 
// Using SetQuenchingFactor(f) quenched jets can be modelled by superimposing 
// two jets with energies e * f and e * (1-f)
//
// andreas.morsch@cern.ch
//


#include "AliGenPythiaJets.h"
#include "AliRun.h"
ClassImp(AliGenPythiaJets)

AliGenPythiaJets::AliGenPythiaJets()
                 :AliGenPythia()
{
// Default Constructor
}

AliGenPythiaJets::AliGenPythiaJets(Int_t npart)
                 :AliGenPythia(npart)
{
    fName = "PythiaJets";
    fTitle= "Jet Generator using PYTHIA";
}

AliGenPythiaJets::AliGenPythiaJets(const AliGenPythiaJets & Pythia)
{
// copy constructor
    Pythia.Copy(*this);
}

AliGenPythiaJets::~AliGenPythiaJets()
{
// Destructor
}

void AliGenPythiaJets::Init()
{
// Initialization
//
    printf("AliGenPythiaJets::Init() \n");
    
    AliGenPythia::Init();
    
    if (fQuench > 0.) {
	 fEtMinJetQ[0]  = fEtMinJet  * fQuench;  
	 fEtMaxJetQ[0]  = fEtMaxJet  * fQuench;  
	 fEtMinJetQ[1]  = fEtMinJet  * (1. - fQuench);  
	 fEtMaxJetQ[1]  = fEtMaxJet  * (1. - fQuench);
	 fPtHardMinQ[0] = fPtHardMin * fQuench;  
	 fPtHardMaxQ[0] = fPtHardMax * fQuench;  
	 fPtHardMinQ[1] = fPtHardMin * (1. - fQuench);  
	 fPtHardMaxQ[1] = fPtHardMax * (1. - fQuench);  
    }
}

void AliGenPythiaJets::Generate()
{
// Generate one event
    fDecayer->ForceDecay();

    Float_t polar[3]   =   {0,0,0};
    Float_t origin[3]  =   {0,0,0};
    Float_t p[3];
//  converts from mm/c to s
    const Float_t kconv=0.001/2.999792458e8;
//
    Int_t nt  = 0;
    Int_t nc  = 0;
    Int_t jev = 0;
    Int_t j, kf, iparent;
    fTrials=0;
//
//  Set collision vertex position 
    if(fVertexSmear==kPerEvent) {
	fPythia->SetMSTP(151,1);
	for (j=0;j<3;j++) {
	    fPythia->SetPARP(151+j, fOsigma[j]*10.);
	}
    } else if (fVertexSmear==kPerTrack) {
	fPythia->SetMSTP(151,0);
    }
//  Event loop    
    while(1)
    {
	if (fQuench > 0.) {
	    fPythia->SetCKIN(3,fPtHardMinQ[jev]);
	    fPythia->SetCKIN(4,fPtHardMaxQ[jev]);
	    fEtMinJet = fEtMinJetQ[jev];
	    fEtMaxJet = fEtMaxJetQ[jev];	    
	}
	
	fPythia->Pyevnt();
	if (gAlice->GetEvNumber()>=fDebugEventFirst &&
	    gAlice->GetEvNumber()<=fDebugEventLast) fPythia->Pylist(1);
	fTrials++;
        //
	// Has this jet triggered
        //
	if ((fEtMinJet != -1) && ! CheckTrigger()) continue;
//
	fPythia->ImportParticles(fParticles,"All");
	Int_t i;
	Int_t np = fParticles->GetEntriesFast();
	if (np == 0 ) continue;
// Get event vertex and discard the event if the z coord. is too big	
	TParticle *iparticle = (TParticle *) fParticles->At(0);
	Float_t distz = iparticle->Vz()/10.;
	if(TMath::Abs(distz)>fCutVertexZ*fOsigma[2]) continue;
//
//
	fEventVertex[0] = iparticle->Vx()/10.+fOrigin.At(0);
	fEventVertex[1] = iparticle->Vy()/10.+fOrigin.At(1);
	fEventVertex[2] = iparticle->Vz()/10.+fOrigin.At(2);

	Int_t* pParent = new Int_t[np];
	for (i=0; i< np; i++) pParent[i] = -1;
	

	//
	for (i = 0; i<np; i++) {
	    Int_t trackIt = 0;
	    TParticle *  iparticle = (TParticle *) fParticles->At(i);
	    kf = CheckPDGCode(iparticle->GetPdgCode());
	    Int_t ks = iparticle->GetStatusCode();
	    Int_t km = iparticle->GetFirstMother();
	    if ((ks == 1  && kf !=0 && KinematicSelection(iparticle, 0)) ||
		(ks != 1) ||
		(fProcess == kPyJets && ks == 21 && km == 0 && i > 1)) {
		nc++;
		if (ks == 1) trackIt = 1;
		Int_t ipa = iparticle->GetFirstMother() - 1;
		iparent = (ipa > -1) ? pParent[ipa] : -1;
//
// Store track information
//
		p[0] = iparticle->Px();
		p[1] = iparticle->Py();
		p[2] = iparticle->Pz();
		origin[0] = fOrigin[0]+iparticle->Vx()/10.;
		origin[1] = fOrigin[1]+iparticle->Vy()/10.;
		origin[2] = fOrigin[2]+iparticle->Vz()/10.;
		Float_t tof=kconv*iparticle->T();
		SetTrack(fTrackIt*trackIt, iparent, kf, p, origin, polar,
			 tof, kPPrimary, nt, 1., ks);
		KeepTrack(nt);
		pParent[i] = nt;
	    } // select particle
	} // particle loop 
	
	if (pParent) delete[] pParent;
	printf("\n AliGenPythiaJets: I've put %i particles on the stack \n",nc);
	if (nc > 0) {
	    jev += 1;
	    if ((fQuench <= 0.) || (fQuench > 0. && jev == 2)) {
		fKineBias=Float_t(fNpart)/Float_t(fTrials);
		printf("\n Trials: %i %i %i\n",fTrials, fNpart, jev);
		fNev++;
		MakeHeader();
		break;
	    }
	}
    }
    SetHighWaterMark(nt);
//  Get cross-section
    fXsection=fPythia->GetPARI(1);
}

Bool_t AliGenPythiaJets::CheckTrigger()
{
// Check the kinematic trigger condition
//
    Bool_t   triggered = kFALSE;
    Int_t njets = 0;
    Int_t ntrig = 0;
    Float_t jets[4][10];
//
// Use Pythia clustering on parton level to determine jet axis
//
    GetJets(njets, ntrig, jets);
    
    if (ntrig) {
	triggered = kTRUE;
	Float_t px   = jets[0][0];
	Float_t py   = jets[1][0];
	Float_t pz   = jets[2][0];
	Float_t e    = jets[3][0];
	Float_t beta = pz/e;
	Float_t phi  =  TMath::Pi()+TMath::ATan2(-py,-px);
	TransformEvent(beta, -2. * TMath::Pi() / 3. + phi);
    }
    return triggered;
}
	  
AliGenPythiaJets& AliGenPythiaJets::operator=(const  AliGenPythiaJets& rhs)
{
// Assignment operator
    return *this;
}

void  AliGenPythiaJets::TransformEvent(Float_t beta, Float_t phi)
{
//
// Perform Lorentz Transformation and Rotation
//
    Float_t gamma = 1./TMath::Sqrt(1. - beta * beta);
    Int_t npart = (fPythia->GetPyjets())->N;

    for (Int_t part = 0; part < npart; part++) {
	Float_t px   = 	(fPythia->GetPyjets())->P[0][part];
	Float_t py   =  (fPythia->GetPyjets())->P[1][part];
	Float_t pz   =  (fPythia->GetPyjets())->P[2][part];
	Float_t e    =  (fPythia->GetPyjets())->P[3][part];
	//
	// Lorentz Transform
	//
	Float_t pzt =  gamma * pz        - gamma * beta * e;
	Float_t et  = -gamma * beta * pz + gamma        * e;
	//
	// Rotation
	//
	Float_t pxt =   TMath::Cos(phi) * px +  TMath::Sin(phi) * py;
	Float_t pyt = - TMath::Sin(phi) * px +  TMath::Cos(phi) * py;
	//
	//
	(fPythia->GetPyjets())->P[0][part] = pxt;
	(fPythia->GetPyjets())->P[1][part] = pyt;
	(fPythia->GetPyjets())->P[2][part] = pzt;
	(fPythia->GetPyjets())->P[3][part] = et;
    }
}

