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


// Generator using DPMJET as an external generator
// The main DPMJET options are accessable for the user through this interface.
// Uses the TDPMjet implementation of TGenerator.

#include <TDPMjet.h>
#include <TRandom.h>
#include <TArrayI.h>
#include <TParticle.h>
#include <TGraph.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TParticleClassPDG.h>
#include <TPDGCode.h>
#include <TLorentzVector.h>

#include "AliGenDPMjet.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliPythia.h"
#include "AliRun.h"
#include "AliDpmJetRndm.h"

 ClassImp(AliGenDPMjet)


//______________________________________________________________________________
AliGenDPMjet::AliGenDPMjet()
                 :AliGenMC()
{
// Constructor
    fParticles = 0;
    fDPMjet = 0;
    AliDpmJetRndm::SetDpmJetRandom(GetRandom());
}


//______________________________________________________________________________
AliGenDPMjet::AliGenDPMjet(Int_t npart)
    :AliGenMC(npart)
{
// Default PbPb collisions at 5. 5 TeV
//
    fName = "DPMJET";
    fTitle= "Particle Generator using DPMJET";

    SetBeamEnergy();
    SetEnergyCMS();
    SetTarget();
    SetProjectile();
    SetCentral();
    SetImpactParameterRange();
    SetBoostLHC();
    
    fKeep       =  0;
    fDecaysOff  =  1;
    fEvaluate   =  0;
    fSelectAll  =  0;
    fFlavor     =  0;
    fSpectators =  1;
    fVertex.Set(3);
        
    fParticles = new TClonesArray("TParticle",10000);    

    // Set random number generator   
    fDPMjet = 0;
    // Instance AliPythia
    AliPythia::Instance(); 
    AliDpmJetRndm::SetDpmJetRandom(GetRandom());
}

AliGenDPMjet::AliGenDPMjet(const AliGenDPMjet &/*Dpmjet*/)
:AliGenMC()
{
}

//______________________________________________________________________________
AliGenDPMjet::~AliGenDPMjet()
{
// Destructor
    delete fParticles;
}


//______________________________________________________________________________
void AliGenDPMjet::Init()
{
// Initialization
    
    SetMC(new TDPMjet(fAProjectile, fZProjectile, fATarget, fZTarget, 
		      fBeamEn,fEnergyCMS));

    fDPMjet=(TDPMjet*) fMCEvGen;
    //
    // **** Flag to force central production
    // fICentr=1. central production forced 
    // fICentr<0 && fICentr>-100 -> bmin = fMinImpactParam, bmax = fMaxImpactParam	  
    // fICentr<-99 -> fraction of x-sec. = XSFRAC		  
    // fICentr=-1. -> evaporation/fzc suppressed		  
    // fICentr<-1. -> evaporation/fzc suppressed		  
    if (fAProjectile == 1 && fZProjectile == 1) fDPMjet->SetfIdp(1);
    
    fDPMjet->SetfFCentr(fICentr);  
    fDPMjet->SetbRange(fMinImpactParam,fMaxImpactParam); 
    
//
//  Initialize DPMjet  
//    
    fDPMjet->Initialize();
    //if (fEvaluate) EvaluateCrossSections();

    //  Initialize AliIonPDGCodes to define ion PDG codes
    /*AliIonPDGCodes *PDGcodes = new AliIonPDGCodes();
    PDGcodes->AddParticlesToPdgDataBase();
    PDGcodes->MapPDGGEant3Codes();*/
}


//______________________________________________________________________________
void AliGenDPMjet::Generate()
{
// Generate one event

  Float_t polar[3]    =   {0,0,0};
  Float_t origin[3]   =   {0,0,0};
  Float_t origin0[3]  =   {0,0,0};
  Float_t p[3], random[6];
  Float_t tof;

//  converts from mm/c to s
  const Float_t kconv = 0.001/2.999792458e8;
  Int_t nt  = 0;
  Int_t jev = 0;
  Int_t j, kf, ks, imo;
  kf = 0;
        
  fTrials = 0;
  for (j = 0;j < 3; j++) origin0[j] = fOrigin[j];
  if(fVertexSmear == kPerEvent) {
      Float_t dv[3];
      dv[2] = 1.e10;
      while(TMath::Abs(dv[2])>fCutVertexZ*fOsigma[2]) {
	  Rndm(random,6);
	  for (j=0; j < 3; j++) {
	      dv[j] = fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		  TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	  }
      }
      for (j=0; j < 3; j++) origin0[j] += dv[j];
  }
  
  while(1)
  {
//    Generate one event
// --------------------------------------------------------------------------
      fSpecn = 0;  
      fSpecp = 0;
// --------------------------------------------------------------------------
      fDPMjet->GenerateEvent();
      fTrials++;

      fDPMjet->ImportParticles(fParticles,"All");      
      if (fLHC) Boost();

      // Temporaneo
      fGenImpPar = fDPMjet->GetBImpac();
      
      Int_t np = fParticles->GetEntriesFast();
      printf("\n **************************************************%d\n",np);
      Int_t nc=0;
      if (np==0) continue;
      Int_t i;
      Int_t* newPos     = new Int_t[np];
      Int_t* pSelected  = new Int_t[np];

      for (i = 0; i<np; i++) {
	  newPos[i]    = i;
	  pSelected[i] = 0;
      }
      
//      Get event vertex
      
      //TParticle *iparticle = (TParticle*) fParticles->At(0);
      fVertex[0] = origin0[0];
      fVertex[1] = origin0[1];	
      fVertex[2] = origin0[2];
      //for(int jj=0; jj<3; jj++) printf("	fEventVertex[%d] = %f\n",jj,fEventVertex[jj]);
      
//      First select parent particles

      for (i = 0; i<np; i++) {
	  TParticle *iparticle = (TParticle *) fParticles->At(i);
          //if(!iparticle) iparticle->Print("");

// Is this a parent particle ?

	  if (Stable(iparticle)) continue;

	  Bool_t  selected             =  kTRUE;
	  Bool_t  hasSelectedDaughters =  kFALSE;
	  
	  kf = iparticle->GetPdgCode();
	  if (kf == 92) continue;
	  ks = iparticle->GetStatusCode();
// No initial state partons
          if (ks==21) continue;
	    
	  if (!fSelectAll) selected = KinematicSelection(iparticle, 0) && 
			       SelectFlavor(kf);
	  hasSelectedDaughters = DaughtersSelection(iparticle);

// Put particle on the stack if it is either selected or 
// it is the mother of at least one seleted particle

	  if (selected || hasSelectedDaughters) {
	      nc++;
	      pSelected[i] = 1;
	  } // selected
      } // particle loop parents

// Now select the final state particles


      for (i=0; i<np; i++) {
	  TParticle *iparticle = (TParticle *) fParticles->At(i);
	  //if(!iparticle) printf("	i = %d -> iparticle==0\n",i);

// Is this a final state particle ?

	  if (!Stable(iparticle)) continue;
      
	  Bool_t  selected =  kTRUE;
	  kf = iparticle->GetPdgCode();
	  ks = iparticle->GetStatusCode();

// --------------------------------------------------------------------------
// Count spectator neutrons and protons (ks == 13, 14)
	  if(ks == 13 || ks == 14){
	      if(kf == kNeutron) fSpecn += 1;
	      if(kf == kProton)  fSpecp += 1;
	  }
// --------------------------------------------------------------------------

	  if (!fSelectAll) {
	      selected = KinematicSelection(iparticle,0)&&SelectFlavor(kf);
	      if (!fSpectators && selected) selected = (ks == 13 || ks == 14);
	  }

// Put particle on the stack if selected

	  if (selected) {
	      nc++;
	      pSelected[i] = 1;
	  } // selected
      } // particle loop final state

// Write particles to stack

      for (i = 0; i<np; i++) {
	  TParticle *  iparticle = (TParticle *) fParticles->At(i);
	  Bool_t  hasMother   = (iparticle->GetFirstMother()>=0);
//	  Bool_t  hasDaughter = (iparticle->GetFirstDaughter()>=0);

	  if (pSelected[i]) {
	      
	      kf   = iparticle->GetPdgCode();	      
	      if(kf==99999) continue;
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

	      // --- Nucleons and photons from nucleus evaporation
	      //if(ks==-1 && (kf==2212 || kf==2112 || kf==22)) imo=-1;
	      
	      // --- Offset added to nuclei, alpha particles, deuteron 
	      // --- and triton PDG codes (according to TGeant3 PDG codes)
	      //if(ks==1000 || ks==1001) kf += 10000000;
	      if(kf>10000) kf += 10000000;
	      
	      if(kf>10000 || ks==1000 || ks==1001) printf(" kf = %d, ks = %d imo = %d \n", kf,ks,imo);
	     

	      Bool_t tFlag = (fTrackIt && (ks == 1));
	      
/*	      // --- Particle NOT to be tracked:
	      // --- (1) chains from valence quark system (idhkk=99999)
	      // --- (2) quark, diquarks and gluons (idhkk=1->8,...,idhkk=21)
	      if(kf==99999 || kf==1 || kf==2 || kf==3 || kf==4 || kf==5 ||
	         kf==6 || kf==7 || kf==8 || kf==21 || 
		 kf==1103 || kf==2101 || kf==2103 || kf==2203 || kf==3101 || 
		 kf==3103 || kf==3201 || kf==3203 || kf==3303 || kf==4101 || 
		 kf==4103 || kf==4201 || kf==4203 || kf==4301 || kf==4303 || 
		 kf==4403 || kf==5101 || kf==5103 || kf==5201 || kf==5203 || 
		 kf==5301 || kf==5303 || kf==5401 || kf==5403 || kf==5503) 
		 tFlag=kFALSE;
*/	      
	      PushTrack(tFlag,imo,kf,p,origin,polar,tof,kPNoProcess,nt, 1., ks);
	      KeepTrack(nt);
	      newPos[i] = nt;
	  } // if selected
      } // particle loop
      delete[] newPos;
      delete[] pSelected;
      
      printf("\n I've put %i particles on the stack \n",nc);
      if (nc>0) {
	  jev += nc;
	  if (jev >= fNpart || fNpart == -1) {
	      printf("\n Trials: %i %i %i\n",fTrials, fNpart, jev);
	      break;
	  }
	  printf("\n Not breaking ### Trials: %i %i %i\n",fTrials, fNpart, jev);
      }
  } // event loop
  MakeHeader();
  SetHighWaterMark(nt);
}


//______________________________________________________________________________
void AliGenDPMjet::KeepFullEvent()
{
    fKeep=1;
}


//______________________________________________________________________________
/*void AliGenDPMjet::EvaluateCrossSections()
{
//     Glauber Calculation of geometrical x-section
//
    Float_t xTot       = 0.;          // barn
    Float_t xTotHard   = 0.;          // barn 
    Float_t xPart      = 0.;          // barn
    Float_t xPartHard  = 0.;          // barn 
    Float_t sigmaHard  = 0.1;         // mbarn
    Float_t bMin       = 0.;
    Float_t bMax       = fDPMjet->GetProjRadius()+fDPMjet->GetTargRadius();
    const Float_t kdib = 0.2;
    Int_t   kMax       = Int_t((bMax-bMin)/kdib)+1;


    printf("\n Projectile Radius (fm): %f \n",fDPMjet->GetProjRadius());
    printf("\n Target     Radius (fm): %f \n",fDPMjet->GetTargRadius());    
    Int_t i;
    Float_t oldvalue= 0.;

    Float_t* b   = new Float_t[kMax];
    Float_t* si1 = new Float_t[kMax];    
    Float_t* si2 = new Float_t[kMax];    
    
    for (i = 0; i < kMax; i++)
    {
	Float_t xb  = bMin+i*kdib;
	Float_t ov;
	ov=fDPMjet->Profile(xb);
	// ATT!->Manca la x-sec anel. nucleone-nucleone
	Float_t gb  =  2.*0.01*fDPMjet->TMath::Pi()*kdib*xb*(1.-TMath::Exp(-fDPMjet->GetXSFrac()*ov));
	Float_t gbh =  2.*0.01*fDPMjet->TMath::Pi()*kdib*xb*sigmaHard*ov;
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
}*/


//______________________________________________________________________________
Bool_t AliGenDPMjet::DaughtersSelection(TParticle* iparticle)
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



//______________________________________________________________________________
Bool_t AliGenDPMjet::SelectFlavor(Int_t pid)
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

//______________________________________________________________________________
Bool_t AliGenDPMjet::Stable(TParticle*  particle)
{
// Return true for a stable particle
//
    
//    if (particle->GetFirstDaughter() < 0 ) return kTRUE;
    if (particle->GetStatusCode() == 1) return kTRUE;
    else return kFALSE;

}



//______________________________________________________________________________
/*void AliGenDPMjet::Boost(TClonesArray* particles)
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
    Int_t np = particles->GetEntriesFast();
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
}*/



//______________________________________________________________________________
void AliGenDPMjet::MakeHeader()
{
// Builds the event header, to be called after each event
    AliGenEventHeader* header = new AliGenDPMjetEventHeader("DPMJET");
    ((AliGenDPMjetEventHeader*) header)->SetNProduced(fDPMjet->GetNumStablePc());
    ((AliGenDPMjetEventHeader*) header)->SetImpactParameter(fDPMjet->GetBImpac());
    ((AliGenDPMjetEventHeader*) header)->SetTotalEnergy(fDPMjet->GetTotEnergy());
//    ((AliGenDPMjetEventHeader*) header)->SetWounded(fDPMjet->GetProjWounded(),
//    						    fDPMjet->GetTargWounded());
    ((AliGenDPMjetEventHeader*) header)->SetParticipants(fDPMjet->GetfIp(), 
    							 fDPMjet->GetfIt());
    /*((AliGenDPMjetEventHeader*) header)->SetCollisions(fDPMjet->GetN0(),
						       fDPMjet->GetN01(),
						       fDPMjet->GetN10(),
						       fDPMjet->GetN11());*/
    ((AliGenDPMjetEventHeader*) header)->SetSpectators(fSpecn, fSpecp);

// Bookkeeping for kinematic bias
    ((AliGenDPMjetEventHeader*) header)->SetTrials(fTrials);
// Event Vertex
    header->SetPrimaryVertex(fVertex);
    gAlice->SetGenEventHeader(header);    
}



//______________________________________________________________________________
AliGenDPMjet& AliGenDPMjet::operator=(const  AliGenDPMjet& /*rhs*/)
{
// Assignment operator
    return *this;
}



//______________________________________________________________________________
