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
Revision 1.17  2000/12/04 11:22:03  morsch
Init of sRandom as in 1.15

Revision 1.16  2000/12/02 11:41:39  morsch
Use SetRandom() to initialize random number generator in constructor.

Revision 1.15  2000/11/30 20:29:02  morsch
Initialise static variable sRandom in constructor: sRandom = fRandom;

Revision 1.14  2000/11/30 07:12:50  alibrary
Introducing new Rndm and QA classes

Revision 1.13  2000/11/09 17:40:27  morsch
Possibility to select/unselect spectator protons and neutrons.
Method SetSpectators(Int_t spect) added. (FCA, Ch. Oppedisano)

Revision 1.12  2000/10/20 13:38:38  morsch
Debug printouts commented.

Revision 1.11  2000/10/20 13:22:26  morsch
- skip particle type 92 (string)
- Charmed and beauty baryions (5122, 4122) are considered as stable consistent with
  mesons.

Revision 1.10  2000/10/17 15:10:20  morsch
Write first all the parent particles to the stack and then the final state particles.

Revision 1.9  2000/10/17 13:38:59  morsch
Protection against division by zero in EvaluateCrossSection() and KinematicSelection(..)     (FCA)

Revision 1.8  2000/10/17 12:46:31  morsch
Protect EvaluateCrossSections() against division by zero.

Revision 1.7  2000/10/02 21:28:06  fca
Removal of useless dependecies via forward declarations

Revision 1.6  2000/09/11 13:23:37  morsch
Write last seed to file (fortran lun 50) and reed back from same lun using calls to
luget_hijing and luset_hijing.

Revision 1.5  2000/09/07 16:55:40  morsch
fHijing->Initialize(); after change of parameters. (Dmitri Yurevitch Peressounko)

Revision 1.4  2000/07/11 18:24:56  fca
Coding convention corrections + few minor bug fixes

Revision 1.3  2000/06/30 12:08:36  morsch
In member data: char* replaced by TString, Init takes care of resizing the strings to
8 characters required by Hijing.

Revision 1.2  2000/06/15 14:15:05  morsch
Add possibility for heavy flavor selection: charm and beauty.

Revision 1.1  2000/06/09 20:47:27  morsch
AliGenerator interface class to HIJING using THijing (test version)

*/

#include "AliGenHijing.h"
#include "AliGenHijingEventHeader.h"
#include "AliRun.h"

#include <TArrayI.h>
#include <TParticle.h>
#include <THijing.h>


 ClassImp(AliGenHijing)

AliGenHijing::AliGenHijing()
                 :AliGenerator()
{
// Constructor
}

AliGenHijing::AliGenHijing(Int_t npart)
    :AliGenerator(npart)
{
// Default PbPb collisions at 5. 5 TeV
//
    SetEnergyCMS();
    SetImpactParameterRange();
    SetTarget();
    SetProjectile();
    fKeep=0;
    fQuench=1;
    fShadowing=1;
    fTrigger=0;
    fDecaysOff=1;
    fEvaluate=0;
    fSelectAll=0;
    fFlavor=0;
    fSpectators=1;
//
// Set random number generator   
    sRandom = fRandom;
}

AliGenHijing::AliGenHijing(const AliGenHijing & Hijing)
{
// copy constructor
}


AliGenHijing::~AliGenHijing()
{
// Destructor
}

void AliGenHijing::Init()
{
// Initialisation
    fFrame.Resize(8);
    fTarget.Resize(8);
    fProjectile.Resize(8);
    
    SetMC(new THijing(fEnergyCMS, fFrame, fProjectile, fTarget, 
		      fAProjectile, fZProjectile, fATarget, fZTarget, 
		      fMinImpactParam, fMaxImpactParam));

    fHijing=(THijing*) fgMCEvGen;

    fHijing->SetIHPR2(3,  fTrigger);
    fHijing->SetIHPR2(4,  fQuench);
    fHijing->SetIHPR2(6,  fShadowing);
    fHijing->SetIHPR2(12, fDecaysOff);    
    fHijing->SetIHPR2(21, fKeep);
    fHijing->Rluset(50,0);
    fHijing->Initialize();

    
//
    if (fEvaluate) EvaluateCrossSections();
//
//
//  Initialize random generator
}

void AliGenHijing::Generate()
{
// Generate one event

    Float_t polar[3] =   {0,0,0};
    Float_t origin[3]=   {0,0,0};
    Float_t origin0[3]=  {0,0,0};
    Float_t p[3], random[6];
    Float_t tof;

    static TClonesArray *particles;
//  converts from mm/c to s
    const Float_t kconv=0.001/2.999792458e8;
//
    Int_t nt=0;
    Int_t jev=0;
    Int_t j, kf, ks, imo;
    kf=0;
    
    if(!particles) particles=new TClonesArray("TParticle",10000);
    
    fTrials=0;
    for (j=0;j<3;j++) origin0[j]=fOrigin[j];
    if(fVertexSmear==kPerEvent) {
	Rndm(random,6);
	for (j=0;j<3;j++) {
	    origin0[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
//	    fHijing->SetMSTP(151,0);
	}
    } else if (fVertexSmear==kPerTrack) {
//	fHijing->SetMSTP(151,0);
	for (j=0;j<3;j++) {
//	    fHijing->SetPARP(151+j, fOsigma[j]*10.);
	}
    }
    while(1)
    {

	fHijing->GenerateEvent();
	fTrials++;
	fHijing->ImportParticles(particles,"All");
	Int_t np = particles->GetEntriesFast();
	printf("\n **************************************************%d\n",np);
	Int_t nc=0;
	if (np == 0 ) continue;
	Int_t i;
	Int_t * newPos = new Int_t[np];

	for (i = 0; i<np; i++) *(newPos+i)=i;
//
//      First write parent particles
//

	for (i = 0; i<np; i++) {
	    TParticle *  iparticle       = (TParticle *) particles->At(i);
// Is this a parent particle ?
	    if (Stable(iparticle)) continue;
//
	    Bool_t  hasMother            =  (iparticle->GetFirstMother()   >=0);
	    Bool_t  selected             =  kTRUE;
	    Bool_t  hasSelectedDaughters =  kFALSE;
	    
	    
	    kf        = iparticle->GetPdgCode();
	    ks        = iparticle->GetStatusCode();
	    if (kf == 92) continue;
	    
	    if (!fSelectAll) selected = KinematicSelection(iparticle)&&SelectFlavor(kf);
	    hasSelectedDaughters = DaughtersSelection(iparticle, particles);
//
// Put particle on the stack if it is either selected or it is the mother of at least one seleted particle
//
	    if (selected || hasSelectedDaughters) {
		nc++;
		p[0]=iparticle->Px();
		p[1]=iparticle->Py();
		p[2]=iparticle->Pz();
		origin[0]=origin0[0]+iparticle->Vx()/10;
		origin[1]=origin0[1]+iparticle->Vy()/10;
		origin[2]=origin0[2]+iparticle->Vz()/10;
		tof=kconv*iparticle->T();
		imo=-1;
		if (hasMother) {
		    imo=iparticle->GetFirstMother();
		    TParticle* mother= (TParticle *) particles->At(imo);
		    imo = (mother->GetPdgCode() != 92) ? imo=*(newPos+imo) : -1;
		}
// Put particle on the stack ... 
//		printf("\n set track mother: %d %d %d %d %d %d ",i,imo, kf, nt+1, selected, hasSelectedDaughters);

		gAlice->SetTrack(0,imo,kf,p,origin,polar,
				 tof,kPPrimary,nt);
// ... and keep it there
		gAlice->KeepTrack(nt);
//
		*(newPos+i)=nt;
	    } // selected
	} // particle loop parents
//
// Now write the final state particles
//

	for (i = 0; i<np; i++) {
	    TParticle *  iparticle       = (TParticle *) particles->At(i);
// Is this a final state particle ?
	    if (!Stable(iparticle)) continue;
//	    
	    Bool_t  hasMother            =  (iparticle->GetFirstMother()   >=0);
	    Bool_t  selected             =  kTRUE;
	    kf        = iparticle->GetPdgCode();
	    ks        = iparticle->GetStatusCode();
	    if (!fSelectAll) {
	      selected = KinematicSelection(iparticle)&&SelectFlavor(kf);
	      if (!fSpectators && selected) selected = (ks != 0 && ks != 10);
	    }
//
// Put particle on the stack if selected
//
	    if (selected) {
		nc++;
		p[0]=iparticle->Px();
		p[1]=iparticle->Py();
		p[2]=iparticle->Pz();
		origin[0]=origin0[0]+iparticle->Vx()/10;
		origin[1]=origin0[1]+iparticle->Vy()/10;
		origin[2]=origin0[2]+iparticle->Vz()/10;
		tof=kconv*iparticle->T();
		imo=-1;

		if (hasMother) {
		    imo=iparticle->GetFirstMother();
		    TParticle* mother= (TParticle *) particles->At(imo);
		    imo = (mother->GetPdgCode() != 92) ? imo=*(newPos+imo) : -1;
		}
// Put particle on the stack
		gAlice->SetTrack(fTrackIt,imo,kf,p,origin,polar,
				 tof,kPNoProcess,nt);
//				 tof,"Secondary",nt);

//		printf("\n set track final: %d %d %d",imo, kf, nt);
		gAlice->KeepTrack(nt);
		*(newPos+i)=nt;
	    } // selected
	} // particle loop final state
 
	delete newPos;

	printf("\n I've put %i particles on the stack \n",nc);
	if (nc > 0) {
	    jev+=nc;
	    if (jev >= fNpart || fNpart == -1) {
		fKineBias=Float_t(fNpart)/Float_t(fTrials);
		printf("\n Trials: %i %i %i\n",fTrials, fNpart, jev);
		break;
	    }
	}
    } // event loop
    fHijing->Rluget(50,-1);
}

Bool_t AliGenHijing::KinematicSelection(TParticle *particle)
{
// Perform kinematic selection
    Double_t px=particle->Px();
    Double_t py=particle->Py();
    Double_t pz=particle->Pz();
    Double_t  e=particle->Energy();

//
//  transverse momentum cut    
    Double_t pt=TMath::Sqrt(px*px+py*py);
    if (pt > fPtMax || pt < fPtMin) 
    {
//	printf("\n failed pt cut %f %f %f \n",pt,fPtMin,fPtMax);
	return kFALSE;
    }
//
// momentum cut
    Double_t p=TMath::Sqrt(px*px+py*py+pz*pz);
    if (p > fPMax || p < fPMin) 
    {
//	printf("\n failed p cut %f %f %f \n",p,fPMin,fPMax);
	return kFALSE;
    }
    
//
// theta cut
    Double_t  theta = Double_t(TMath::ATan2(Double_t(pt),Double_t(pz)));
    if (theta > fThetaMax || theta < fThetaMin) 
    {
	
//    	printf("\n failed theta cut %f %f %f \n",theta,fThetaMin,fThetaMax);
	return kFALSE;
    }

//
// rapidity cut
    Double_t y;
    if(e<=pz) y = 99;
    else if (e<=-pz)  y = -99;
    else y = 0.5*TMath::Log((e+pz)/(e-pz));
    if (y > fYMax || y < fYMin)
    {
//	printf("\n failed y cut %f %f %f \n",y,fYMin,fYMax);
	return kFALSE;
    }

//
// phi cut
    Double_t phi=Double_t(TMath::ATan2(Double_t(py),Double_t(px)));
    if (phi > fPhiMax || phi < fPhiMin)
    {
//	printf("\n failed phi cut %f %f %f \n",phi,fPhiMin,fPhiMax);
	return kFALSE;
    }

    return kTRUE;
}

void AliGenHijing::KeepFullEvent()
{
    fKeep=1;
}

void AliGenHijing::EvaluateCrossSections()
{
//     Glauber Calculation of geometrical x-section
//
    Float_t xTot=0.;          // barn
    Float_t xTotHard=0.;      // barn 
    Float_t xPart=0.;         // barn
    Float_t xPartHard=0.;     // barn 
    Float_t sigmaHard=0.1;    // mbarn
    Float_t bMin=0.;
    Float_t bMax=fHijing->GetHIPR1(34)+fHijing->GetHIPR1(35);
    const Float_t kdib=0.2;
    Int_t   kMax=Int_t((bMax-bMin)/kdib)+1;


    printf("\n Projectile Radius (fm): %f \n",fHijing->GetHIPR1(34));
    printf("\n Target     Radius (fm): %f \n",fHijing->GetHIPR1(35));    
    Int_t i;
    Float_t oldvalue=0.;
    
    for (i=0; i<kMax; i++)
    {
	Float_t xb=bMin+i*kdib;
	Float_t ov;
	ov=fHijing->Profile(xb);
	Float_t gb =  2.*0.01*fHijing->GetHIPR1(40)*kdib*xb*(1.-TMath::Exp(-fHijing->GetHINT1(12)*ov));
	Float_t gbh = 2.*0.01*fHijing->GetHIPR1(40)*kdib*xb*sigmaHard*ov;
	xTot+=gb;
	xTotHard+=gbh;
	if (xb > fMinImpactParam && xb < fMaxImpactParam)
	{
	    xPart+=gb;
	    xPartHard+=gbh;
	}
	
	if(oldvalue) if ((xTot-oldvalue)/oldvalue<0.0001) break;
	oldvalue=xTot;
	printf("\n Total cross section (barn): %d %f %f \n",i, xb, xTot);
	printf("\n Hard  cross section (barn): %d %f %f \n\n",i, xb, xTotHard);
    }
    printf("\n Total cross section (barn): %f \n",xTot);
    printf("\n Hard  cross section (barn): %f \n \n",xTotHard);
    printf("\n Partial       cross section (barn): %f %f \n",xPart, xPart/xTot*100.);
    printf("\n Partial  hard cross section (barn): %f %f \n",xPartHard, xPartHard/xTotHard*100.);
}

Bool_t AliGenHijing::DaughtersSelection(TParticle* iparticle, TClonesArray* particles)
{
//
// Looks recursively if one of the daughters has been selected
//
//    printf("\n Consider daughters %d:",iparticle->GetPdgCode());
    Int_t imin=-1;
    Int_t imax=-1;
    Int_t i;
    Bool_t hasDaughters= (iparticle->GetFirstDaughter() >=0);
    Bool_t selected=kFALSE;
    if (hasDaughters) {
	imin=iparticle->GetFirstDaughter();
	imax=iparticle->GetLastDaughter();       
	for (i=imin; i<= imax; i++){
	    TParticle *  jparticle       = (TParticle *) particles->At(i);	
	    Int_t ip=jparticle->GetPdgCode();
	    if (KinematicSelection(jparticle)&&SelectFlavor(ip)) {
		selected=kTRUE; break;
	    }
	    if (DaughtersSelection(jparticle, particles)) {selected=kTRUE; break; }
	}
    } else {
	return kFALSE;
    }

    return selected;
}


Bool_t AliGenHijing::SelectFlavor(Int_t pid)
{
// Select flavor of particle
// 0: all
// 4: charm and beauty
// 5: beauty
    if (fFlavor == 0) return kTRUE;
    
    Int_t ifl=TMath::Abs(pid/100);
    if (ifl > 10) ifl/=10;
    return (fFlavor == ifl);
}

Bool_t AliGenHijing::Stable(TParticle*  particle)
{
    Int_t kf = TMath::Abs(particle->GetPdgCode());
    
    if ( (particle->GetFirstDaughter() < 0 ) || (kf == 1000*fFlavor+122))
	 
    {
	return kTRUE;
    } else {
	return kFALSE;
    }
}

void AliGenHijing::MakeHeader()
{
// Builds the event header, to be called after each event
    AliGenHijingEventHeader* header = new AliGenHijingEventHeader("Hijing");
//    header->SetDate(date);
//    header->SetRunNumber(run);
//    header->SetEventNumber(event);
    header->SetNProduced(fHijing->GetNATT());
    header->SetImpactParameter(fHijing->GetHINT1(19));
    header->SetTotalEnergy(fHijing->GetEATT());
    header->SetHardScatters(fHijing->GetJATT());
    header->SetParticipants(fHijing->GetNP(), fHijing->GetNT());
    header->SetCollisions(fHijing->GetN0(),
			  fHijing->GetN01(),
			  fHijing->GetN10(),
			  fHijing->GetN11());
}

AliGenHijing& AliGenHijing::operator=(const  AliGenHijing& rhs)
{
// Assignment operator
    return *this;
}

#ifndef WIN32
# define rluget_hijing rluget_hijing_
# define rluset_hijing rluset_hijing_
# define rlu_hijing    rlu_hijing_
# define type_of_call
#else
# define rluget_hijing RLUGET_HIJING
# define rluset_hijing RLUSET_HIJING
# define rlu_hijing    RLU_HIJING
# define type_of_call _stdcall
#endif


extern "C" {
  void type_of_call rluget_hijing(Int_t & /*lfn*/, Int_t & /*move*/)
  {printf("Dummy version of rluget_hijing reached\n");}

  void type_of_call rluset_hijing(Int_t & /*lfn*/, Int_t & /*move*/)
  {printf("Dummy version of rluset_hijing reached\n");}

  Double_t type_of_call rlu_hijing(Int_t & /*idum*/) 
  {
      Float_t r;
      do r=sRandom->Rndm(); while(0 >= r || r >= 1);
      return r;
  }
}
