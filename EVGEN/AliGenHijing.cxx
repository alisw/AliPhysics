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
Revision 1.46  2003/01/07 14:12:33  morsch
Provides collision geometry.

Revision 1.45  2002/12/16 09:44:49  morsch
Default for fRadiation is 3.

Revision 1.44  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.42.4.1  2002/08/28 15:06:50  alibrary
Updating to v3-09-01

Revision 1.43  2002/08/09 12:09:52  morsch
Direct gamma trigger correctly included.

Revision 1.42  2002/03/12 11:07:08  morsch
Add status code of particle to SetTrack call.

Revision 1.41  2002/03/05 11:25:33  morsch
- New quenching options
- Correction in CheckTrigger()

Revision 1.40  2002/02/12 11:05:53  morsch
Get daughter indices right.

Revision 1.39  2002/02/12 09:16:39  morsch
Correction in SelectFlavor()

Revision 1.38  2002/02/12 08:53:21  morsch
SetNoGammas can be used to inhibit writing of gammas and pi0.

Revision 1.37  2002/02/08 16:50:50  morsch
Add name and title in constructor.

Revision 1.36  2002/01/31 20:17:55  morsch
Allow for triggered jets with simplified topology: Exact pT, back-to-back

Revision 1.35  2001/12/13 07:56:25  hristov
Set pointers to zero in the default constructor

Revision 1.34  2001/12/11 16:55:42  morsch
Correct initialization for jet phi-range.

Revision 1.33  2001/12/05 10:18:51  morsch
Possibility of kinematic biasing of jet phi range. (J. Klay)

Revision 1.32  2001/11/28 13:51:11  morsch
Introduce kinematic biasing (etamin, etamax) of jet trigger. Bookkeeping
(number of trials) done in AliGenHijingEventHeader.

Revision 1.31  2001/11/06 12:30:34  morsch
Add Boost() method to boost all particles to LHC lab frame. Needed for asymmetric collision systems.

Revision 1.30  2001/10/21 18:35:56  hristov
Several pointers were set to zero in the default constructors to avoid memory management problems

Revision 1.29  2001/10/15 08:12:24  morsch
- Vertex smearing with truncated gaussian.
- Store triggered jet info before and after final state radiation into mc-heade

Revision 1.28  2001/10/08 11:55:25  morsch
Store 4-momenta of trigegred jets in event header.
Possibility to switch of initial and final state radiation.

Revision 1.27  2001/10/08 07:13:14  morsch
Add setter for minimum transverse momentum of triggered jet.

Revision 1.26  2001/10/04 08:12:24  morsch
Redefinition of stable condition.

Revision 1.25  2001/07/27 17:09:36  morsch
Use local SetTrack, KeepTrack and SetHighWaterMark methods
to delegate either to local stack or to stack owned by AliRun.
(Piotr Skowronski, A.M.)

Revision 1.24  2001/07/20 09:34:56  morsch
Count the number of spectator neutrons and protons and add information
to the event header. (Chiara Oppedisano)

Revision 1.23  2001/07/13 17:30:22  morsch
Derive from AliGenMC.

Revision 1.22  2001/06/11 13:09:23  morsch
- Store cross-Section and number of binary collisions as a function of impact parameter
- Pass AliGenHijingEventHeader to gAlice.

Revision 1.21  2001/02/28 17:35:24  morsch
Consider elastic interactions (ks = 1 and ks = 11) as spectator (Chiara Oppedisano)

Revision 1.20  2001/02/14 15:50:40  hristov
The last particle in event marked using SetHighWaterMark

Revision 1.19  2000/12/21 16:24:06  morsch
Coding convention clean-up

Revision 1.18  2000/12/06 17:46:30  morsch
Avoid random numbers 1 and 0.

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



// Generator using HIJING as an external generator
// The main HIJING options are accessable for the user through this interface.
// Uses the THijing implementation of TGenerator.
//
// andreas.morsch@cern.ch

#include <TArrayI.h>
#include <TGraph.h>
#include <THijing.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TParticle.h>

#include "AliGenHijing.h"
#include "AliGenHijingEventHeader.h"
#include "AliRun.h"


  ClassImp(AliGenHijing)

AliGenHijing::AliGenHijing()
                 :AliGenMC()
{
// Constructor
    fParticles = 0;
    fHijing    = 0;
    fDsigmaDb  = 0;
    fDnDb      = 0;
}

AliGenHijing::AliGenHijing(Int_t npart)
    :AliGenMC(npart)
{
// Default PbPb collisions at 5. 5 TeV
//
    fName = "Hijing";
    fTitle= "Particle Generator using HIJING";

    SetEnergyCMS();
    SetImpactParameterRange();
    SetTarget();
    SetProjectile();
    SetBoostLHC();
    SetJetEtaRange();
    SetJetPhiRange();
    
    fKeep       =  0;
    fQuench     =  1;
    fShadowing  =  1;
    fTrigger    =  0;
    fDecaysOff  =  1;
    fEvaluate   =  0;
    fSelectAll  =  0;
    fFlavor     =  0;
    fSpectators =  1;
    fDsigmaDb   =  0;
    fDnDb       =  0;
    fPtMinJet   = -2.5; 	
    fRadiation  =  3;
    fEventVertex.Set(3);
//
    SetSimpleJets();
    SetNoGammas();
//
    fParticles = new TClonesArray("TParticle",10000);    
//
// Set random number generator   
    sRandom = fRandom;
    fHijing = 0;

}

AliGenHijing::AliGenHijing(const AliGenHijing & Hijing)
{
// copy constructor
}


AliGenHijing::~AliGenHijing()
{
// Destructor
    if ( fDsigmaDb) delete  fDsigmaDb;  
    if ( fDnDb)     delete  fDnDb;  
    delete fParticles;
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
    fHijing->SetIHPR2(2,  fRadiation);
    fHijing->SetIHPR2(3,  fTrigger);
    fHijing->SetIHPR2(6,  fShadowing);
    fHijing->SetIHPR2(12, fDecaysOff);    
    fHijing->SetIHPR2(21, fKeep);
    fHijing->SetHIPR1(10, fPtMinJet); 	
    fHijing->SetHIPR1(50, fSimpleJet);
//
//  Quenching
//
//
//  fQuench = 0:  no quenching
//  fQuench = 1:  hijing default
//  fQuench = 2:  new LHC  parameters for HIPR1(11) and HIPR1(14)
//  fQuench = 3:  new RHIC parameters for HIPR1(11) and HIPR1(14)
//  fQuench = 4:  new LHC  parameters with log(e) dependence
//  fQuench = 5:  new RHIC parameters with log(e) dependence
    fHijing->SetIHPR2(50, 0);
    if (fQuench > 0) 
	fHijing->SetIHPR2(4,  1);
    else
	fHijing->SetIHPR2(4,  0);
// New LHC parameters from Xin-Nian Wang
    if (fQuench == 2) {
	fHijing->SetHIPR1(14, 1.1);
	fHijing->SetHIPR1(11, 3.7);
    } else if (fQuench == 3) {
	fHijing->SetHIPR1(14, 0.20);
	fHijing->SetHIPR1(11, 2.5);
    } else if (fQuench == 4) {
	fHijing->SetIHPR2(50, 1);
	fHijing->SetHIPR1(14, 4.*0.34);
	fHijing->SetHIPR1(11, 3.7);
    } else if (fQuench == 5) {
	fHijing->SetIHPR2(50, 1);
	fHijing->SetHIPR1(14, 0.34);
	fHijing->SetHIPR1(11, 2.5);
    }
    
    
    
//
//  Initialize Hijing  
//    
    fHijing->Initialize();
//
    if (fEvaluate) EvaluateCrossSections();
//
}

void AliGenHijing::Generate()
{
// Generate one event

  Float_t polar[3]    =   {0,0,0};
  Float_t origin[3]   =   {0,0,0};
  Float_t origin0[3]  =   {0,0,0};
  Float_t p[3], random[6];
  Float_t tof;

//  converts from mm/c to s
  const Float_t kconv = 0.001/2.999792458e8;
//
  Int_t nt  = 0;
  Int_t jev = 0;
  Int_t j, kf, ks, imo;
  kf = 0;
    

    
  fTrials = 0;
  for (j = 0;j < 3; j++) origin0[j] = fOrigin[j];
  if(fVertexSmear == kPerEvent) {
      Float_t dv[3];
      dv[2] = 1.e10;
      while(TMath::Abs(dv[2]) > fCutVertexZ*fOsigma[2]) {
	  Rndm(random,6);
	  for (j=0; j < 3; j++) {
	      dv[j] = fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		  TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	  }
      }
      for (j=0; j < 3; j++) origin0[j] += dv[j];
  } else if (fVertexSmear == kPerTrack) {
//	    fHijing->SetMSTP(151,0);
      for (j = 0; j < 3; j++) {
//	      fHijing->SetPARP(151+j, fOsigma[j]*10.);
      }
  }
  while(1)
  {
//    Generate one event
// --------------------------------------------------------------------------
      fSpecn	= 0;  
      fSpecp	= 0;
// --------------------------------------------------------------------------
      fHijing->GenerateEvent();
      fTrials++;
      fHijing->ImportParticles(fParticles,"All");
      if (fTrigger != kNoTrigger) {
	  if (!CheckTrigger()) continue;
      }
      if (fLHC) Boost();
      
      
      Int_t np = fParticles->GetEntriesFast();
      printf("\n **************************************************%d\n",np);
      Int_t nc = 0;
      if (np == 0 ) continue;
      Int_t i;
      Int_t* newPos     = new Int_t[np];
      Int_t* pSelected  = new Int_t[np];

      for (i = 0; i < np; i++) {
	  newPos[i]    = i;
	  pSelected[i] = 0;
      }
      
//      Get event vertex
//
      TParticle *  iparticle = (TParticle *) fParticles->At(0);
      fEventVertex[0] = origin0[0];
      fEventVertex[1] = origin0[1];	
      fEventVertex[2] = origin0[2];
      
//
//      First select parent particles
//

      for (i = 0; i < np; i++) {
	  iparticle = (TParticle *) fParticles->At(i);

// Is this a parent particle ?
	  if (Stable(iparticle)) continue;
//
	  Bool_t  selected             =  kTRUE;
	  Bool_t  hasSelectedDaughters =  kFALSE;
	  
	  
	  kf        = iparticle->GetPdgCode();
	  ks        = iparticle->GetStatusCode();
	  if (kf == 92) continue;
	    
	  if (!fSelectAll) selected = KinematicSelection(iparticle, 0) && 
			       SelectFlavor(kf);
	  hasSelectedDaughters = DaughtersSelection(iparticle);
//
// Put particle on the stack if it is either selected or 
// it is the mother of at least one seleted particle
//
	  if (selected || hasSelectedDaughters) {
	      nc++;
	      pSelected[i] = 1;
	  } // selected
      } // particle loop parents
//
// Now select the final state particles
//

      for (i = 0; i<np; i++) {
	  TParticle *  iparticle = (TParticle *) fParticles->At(i);
// Is this a final state particle ?
	  if (!Stable(iparticle)) continue;
      
	  Bool_t  selected             =  kTRUE;
	  kf        = iparticle->GetPdgCode();
	  ks        = iparticle->GetStatusCode();
	  
// --------------------------------------------------------------------------
// Count spectator neutrons and protons
	  if(ks == 0 || ks == 1 || ks == 10 || ks == 11){
	      if(kf == kNeutron) fSpecn += 1;
	      if(kf == kProton)  fSpecp += 1;
	  }
// --------------------------------------------------------------------------
//	    
	  if (!fSelectAll) {
	      selected = KinematicSelection(iparticle,0)&&SelectFlavor(kf);
	      if (!fSpectators && selected) selected = (ks != 0 && ks != 1 && ks != 10
							&& ks != 11);
	  }
//
// Put particle on the stack if selected
//
	  if (selected) {
	      nc++;
	      pSelected[i] = 1;
	  } // selected
      } // particle loop final state
//
// Write particles to stack
//
      for (i = 0; i<np; i++) {
	  TParticle *  iparticle = (TParticle *) fParticles->At(i);
	  Bool_t  hasMother   = (iparticle->GetFirstMother()     >=0);
	  Bool_t  hasDaughter = (iparticle->GetFirstDaughter()   >=0);

	  if (pSelected[i]) {
	      kf   = iparticle->GetPdgCode();
	      ks   = iparticle->GetStatusCode();
	      p[0] = iparticle->Px();
	      p[1] = iparticle->Py();
	      p[2] = iparticle->Pz();
	      origin[0] = origin0[0]+iparticle->Vx()/10;
	      origin[1] = origin0[1]+iparticle->Vy()/10;
	      origin[2] = origin0[2]+iparticle->Vz()/10;
	      tof = kconv*iparticle->T();
	      imo = -1;
	      TParticle* mother = 0;
	      if (hasMother) {
		  imo = iparticle->GetFirstMother();
		  mother = (TParticle *) fParticles->At(imo);
		  imo = (mother->GetPdgCode() != 92) ? imo = newPos[imo] : -1;
	      } // if has mother   
	      Bool_t tFlag = (fTrackIt && !hasDaughter);
	      SetTrack(tFlag,imo,kf,p,origin,polar,
		       tof,kPNoProcess,nt, 1., ks);
	      KeepTrack(nt);
	      newPos[i] = nt;
	  } // if selected
      } // particle loop
      delete[] newPos;
      delete[] pSelected;
      
      printf("\n I've put %i particles on the stack \n",nc);
      if (nc > 0) {
	  jev += nc;
	  if (jev >= fNpart || fNpart == -1) {
	      fKineBias = Float_t(fNpart)/Float_t(fTrials);
	      printf("\n Trials: %i %i %i\n",fTrials, fNpart, jev);
	      break;
	  }
      }
  } // event loop
  MakeHeader();
  SetHighWaterMark(nt);
}

void AliGenHijing::KeepFullEvent()
{
    fKeep=1;
}

void AliGenHijing::EvaluateCrossSections()
{
//     Glauber Calculation of geometrical x-section
//
    Float_t xTot       = 0.;          // barn
    Float_t xTotHard   = 0.;          // barn 
    Float_t xPart      = 0.;          // barn
    Float_t xPartHard  = 0.;          // barn 
    Float_t sigmaHard  = 0.1;         // mbarn
    Float_t bMin       = 0.;
    Float_t bMax       = fHijing->GetHIPR1(34)+fHijing->GetHIPR1(35);
    const Float_t kdib = 0.2;
    Int_t   kMax       = Int_t((bMax-bMin)/kdib)+1;


    printf("\n Projectile Radius (fm): %f \n",fHijing->GetHIPR1(34));
    printf("\n Target     Radius (fm): %f \n",fHijing->GetHIPR1(35));    
    Int_t i;
    Float_t oldvalue= 0.;

    Float_t* b   = new Float_t[kMax];
    Float_t* si1 = new Float_t[kMax];    
    Float_t* si2 = new Float_t[kMax];    
    
    for (i = 0; i < kMax; i++)
    {
	Float_t xb  = bMin+i*kdib;
	Float_t ov;
	ov=fHijing->Profile(xb);
	Float_t gb  =  2.*0.01*fHijing->GetHIPR1(40)*kdib*xb*(1.-TMath::Exp(-fHijing->GetHINT1(12)*ov));
	Float_t gbh =  2.*0.01*fHijing->GetHIPR1(40)*kdib*xb*sigmaHard*ov;
	xTot+=gb;
	xTotHard += gbh;
	if (xb > fMinImpactParam && xb < fMaxImpactParam)
	{
	    xPart += gb;
	    xPartHard += gbh;
	}
	
	if(oldvalue) if ((xTot-oldvalue)/oldvalue<0.0001) break;
	oldvalue = xTot;
	printf("\n Total cross section (barn): %d %f %f \n",i, xb, xTot);
	printf("\n Hard  cross section (barn): %d %f %f \n\n",i, xb, xTotHard);
	if (i>0) {
	    si1[i] = gb/kdib;
	    si2[i] = gbh/gb;
	    b[i]  = xb;
	}
    }

    printf("\n Total cross section (barn): %f \n",xTot);
    printf("\n Hard  cross section (barn): %f \n \n",xTotHard);
    printf("\n Partial       cross section (barn): %f %f \n",xPart, xPart/xTot*100.);
    printf("\n Partial  hard cross section (barn): %f %f \n",xPartHard, xPartHard/xTotHard*100.);

//  Store result as a graph
    b[0] = 0;
    si1[0] = 0;
    si2[0]=si2[1];
    
    fDsigmaDb  = new TGraph(i, b, si1);
    fDnDb      = new TGraph(i, b, si2);
}

Bool_t AliGenHijing::DaughtersSelection(TParticle* iparticle)
{
//
// Looks recursively if one of the daughters has been selected
//
//    printf("\n Consider daughters %d:",iparticle->GetPdgCode());
    Int_t imin = -1;
    Int_t imax = -1;
    Int_t i;
    Bool_t hasDaughters = (iparticle->GetFirstDaughter() >=0);
    Bool_t selected = kFALSE;
    if (hasDaughters) {
	imin = iparticle->GetFirstDaughter();
	imax = iparticle->GetLastDaughter();       
	for (i = imin; i <= imax; i++){
	    TParticle *  jparticle = (TParticle *) fParticles->At(i);	
	    Int_t ip = jparticle->GetPdgCode();
	    if (KinematicSelection(jparticle,0)&&SelectFlavor(ip)) {
		selected=kTRUE; break;
	    }
	    if (DaughtersSelection(jparticle)) {selected=kTRUE; break; }
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
    Bool_t res = 0;
    
    if (fFlavor == 0) {
	res = kTRUE;
    } else {
	Int_t ifl = TMath::Abs(pid/100);
	if (ifl > 10) ifl/=10;
	res = (fFlavor == ifl);
    }
//
//  This part if gamma writing is inhibited
    if (fNoGammas) 
	res = res && (pid != kGamma && pid != kPi0);
//
    return res;
}

Bool_t AliGenHijing::Stable(TParticle*  particle)
{
// Return true for a stable particle
//
    
    if (particle->GetFirstDaughter() < 0 )
    {
	return kTRUE;
    } else {
	return kFALSE;
    }
}


void AliGenHijing::Boost()
{
//
// Boost cms into LHC lab frame
//
    Double_t dy    = - 0.5 * TMath::Log(Double_t(fZProjectile) * Double_t(fATarget) / 
				      (Double_t(fZTarget)    * Double_t(fAProjectile)));
    Double_t beta  = TMath::TanH(dy);
    Double_t gamma = 1./TMath::Sqrt(1.-beta*beta);
    Double_t gb    = gamma * beta;

    printf("\n Boosting particles to lab frame %f %f %f", dy, beta, gamma);
    
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


void AliGenHijing::MakeHeader()
{
// Builds the event header, to be called after each event
    AliGenEventHeader* header = new AliGenHijingEventHeader("Hijing");
    ((AliGenHijingEventHeader*) header)->SetNProduced(fHijing->GetNATT());
    ((AliGenHijingEventHeader*) header)->SetImpactParameter(fHijing->GetHINT1(19));
    ((AliGenHijingEventHeader*) header)->SetTotalEnergy(fHijing->GetEATT());
    ((AliGenHijingEventHeader*) header)->SetHardScatters(fHijing->GetJATT());
    ((AliGenHijingEventHeader*) header)->SetParticipants(fHijing->GetNP(), fHijing->GetNT());
    ((AliGenHijingEventHeader*) header)->SetCollisions(fHijing->GetN0(),
						       fHijing->GetN01(),
						       fHijing->GetN10(),
						       fHijing->GetN11());
    ((AliGenHijingEventHeader*) header)->SetSpectators(fSpecn, fSpecp);

// 4-momentum vectors of the triggered jets.
//
// Before final state gluon radiation.
    TLorentzVector* jet1 = new TLorentzVector(fHijing->GetHINT1(21), 
					      fHijing->GetHINT1(22),
					      fHijing->GetHINT1(23),
					      fHijing->GetHINT1(24));

    TLorentzVector* jet2 = new TLorentzVector(fHijing->GetHINT1(31), 
					      fHijing->GetHINT1(32),
					      fHijing->GetHINT1(33),
					      fHijing->GetHINT1(34));
// After final state gluon radiation.
    TLorentzVector* jet3 = new TLorentzVector(fHijing->GetHINT1(26), 
					      fHijing->GetHINT1(27),
					      fHijing->GetHINT1(28),
					      fHijing->GetHINT1(29));

    TLorentzVector* jet4 = new TLorentzVector(fHijing->GetHINT1(36), 
					      fHijing->GetHINT1(37),
					      fHijing->GetHINT1(38),
					      fHijing->GetHINT1(39));
    ((AliGenHijingEventHeader*) header)->SetJets(jet1, jet2, jet3, jet4);
// Bookkeeping for kinematic bias
    ((AliGenHijingEventHeader*) header)->SetTrials(fTrials);
// Event Vertex
    header->SetPrimaryVertex(fEventVertex);
    gAlice->SetGenEventHeader(header);   
    fCollisionGeometry = (AliGenHijingEventHeader*)  header;
}

Bool_t AliGenHijing::CheckTrigger()
{
// Check the kinematic trigger condition
//
    Bool_t   triggered = kFALSE;
 
    if (fTrigger == 1) {
//
//  jet-jet Trigger	
	
	TLorentzVector* jet1 = new TLorentzVector(fHijing->GetHINT1(26), 
						  fHijing->GetHINT1(27),
						  fHijing->GetHINT1(28),
						  fHijing->GetHINT1(29));
	
	TLorentzVector* jet2 = new TLorentzVector(fHijing->GetHINT1(36), 
						  fHijing->GetHINT1(37),
						  fHijing->GetHINT1(38),
						  fHijing->GetHINT1(39));
	Double_t eta1      = jet1->Eta();
	Double_t eta2      = jet2->Eta();
	Double_t phi1      = jet1->Phi();
	Double_t phi2      = jet2->Phi();
//    printf("\n Trigger: %f %f %f %f",
//	   fEtaMinJet, fEtaMaxJet, fPhiMinJet, fPhiMaxJet);
	if (
	    (eta1 < fEtaMaxJet && eta1 > fEtaMinJet &&  
	     phi1 < fPhiMaxJet && phi1 > fPhiMinJet) 
	    ||
	    (eta2 < fEtaMaxJet && eta2 > fEtaMinJet &&  
	     phi2 < fPhiMaxJet && phi2 > fPhiMinJet)
	    ) 
	    triggered = kTRUE;
    } else if (fTrigger == 2) {
//  Gamma Jet
//
	Int_t np = fParticles->GetEntriesFast();
	for (Int_t i = 0; i < np; i++) {
	    TParticle* part = (TParticle*) fParticles->At(i);
	    Int_t kf = part->GetPdgCode();
	    Int_t ks = part->GetStatusCode();
	    if (kf == 22 && ks == 40) {
		Float_t phi = part->Phi();
		Float_t eta = part->Eta();
		if  (eta < fEtaMaxJet && 
		     eta > fEtaMinJet &&
		     phi < fPhiMaxJet && 
		     phi > fPhiMinJet) {
		    triggered = 1;
		    break;
		} // check phi,eta within limits
	    } // direct gamma ? 
	} // particle loop
    } // fTrigger == 2
    return triggered;
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
