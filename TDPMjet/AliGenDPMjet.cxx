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
#include "AliRunLoader.h"
#include "AliGenDPMjet.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliRun.h"
#include "AliDpmJetRndm.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliMC.h"

ClassImp(AliGenDPMjet)

//______________________________________________________________________________
AliGenDPMjet::AliGenDPMjet()
    :AliGenMC(), 
     fBeamEn(2750.),
     fEnergyCMS(5500.),
     fMinImpactParam(0.),
     fMaxImpactParam(5.),
     fICentr(0),
     fSelectAll(0),
     fFlavor(0),
     fTrials(0),
     fSpectators(1),
     fSpecn(0),
     fSpecp(0),
     fDPMjet(0),
     fParticles(0),
     fNoGammas(0),
     fLHC(0),
     fPi0Decay(0),
     fGenImpPar(0.),
     fProcess(kDpmMb)
{
// Constructor
    AliDpmJetRndm::SetDpmJetRandom(GetRandom());
}


//______________________________________________________________________________
AliGenDPMjet::AliGenDPMjet(Int_t npart)
    :AliGenMC(npart),
     fBeamEn(2750.),
     fEnergyCMS(5500.),
     fMinImpactParam(0.),
     fMaxImpactParam(5.),
     fICentr(0),
     fSelectAll(0),
     fFlavor(0),
     fTrials(0),
     fSpectators(1),
     fSpecn(0),
     fSpecp(0),
     fDPMjet(0),
     fParticles(new TClonesArray("TParticle",10000)),
     fNoGammas(0),
     fLHC(0),
     fPi0Decay(0),
     fGenImpPar(0.),
     fProcess(kDpmMb)
{
// Default PbPb collisions at 5. 5 TeV
//
    fName = "DPMJET";
    fTitle= "Particle Generator using DPMJET";
    SetTarget();
    SetProjectile();
    fVertex.Set(3);
    AliDpmJetRndm::SetDpmJetRandom(GetRandom());
}

AliGenDPMjet::AliGenDPMjet(const AliGenDPMjet &/*Dpmjet*/)
    :AliGenMC(),
     fBeamEn(2750.),
     fEnergyCMS(5500.),
     fMinImpactParam(0.),
     fMaxImpactParam(5.),
     fICentr(0),
     fSelectAll(0),
     fFlavor(0),
     fTrials(0),
     fSpectators(1),
     fSpecn(0),
     fSpecp(0),
     fDPMjet(0),
     fParticles(0),
     fNoGammas(0),
     fLHC(0),
     fPi0Decay(0),
     fGenImpPar(0.),
     fProcess(kDpmMb)
{
    // Dummy copy constructor
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
    
    SetMC(new TDPMjet(fProcess, fAProjectile, fZProjectile, fATarget, fZTarget, 
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
    fDPMjet->SetbRange(fMinImpactParam, fMaxImpactParam); 
    fDPMjet->SetPi0Decay(fPi0Decay);
//
//  Initialize DPMjet  
//    
    fDPMjet->Initialize();
}


//______________________________________________________________________________
void AliGenDPMjet::Generate()
{
// Generate one event

  Float_t polar[3]    =   {0,0,0};
  Float_t origin[3]   =   {0,0,0};
  Float_t origin0[3]  =   {0,0,0};
  Float_t p[3];
  Float_t tof;

//  converts from mm/c to s
  const Float_t kconv = 0.001/2.999792458e8;
  Int_t nt  = 0;
  Int_t jev = 0;
  Int_t kf, ks, imo;
  kf = 0;
  fTrials = 0;
  //  Set collision vertex position 
  if (fVertexSmear == kPerEvent) Vertex();
  
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
      
      fVertex[0] = origin0[0];
      fVertex[1] = origin0[1];	
      fVertex[2] = origin0[2];
      
//      First select parent particles

      for (i = 0; i<np; i++) {
	  TParticle *iparticle = (TParticle *) fParticles->At(i);

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
	  if (pSelected[i]) {
	      
	      kf   = iparticle->GetPdgCode();	      
	      ks   = iparticle->GetStatusCode();	      
	      
	      p[0] = iparticle->Px();
	      p[1] = iparticle->Py();
	      p[2] = iparticle->Pz();
	      origin[0] = fVertex[0]+iparticle->Vx()/10; // [cm]
	      origin[1] = fVertex[1]+iparticle->Vy()/10; // [cm]
	      origin[2] = fVertex[2]+iparticle->Vz()/10; // [cm]
		    
	      tof = kconv*iparticle->T();
	      
	      imo = -1;
	      TParticle* mother = 0;
	      if (hasMother) {
		  imo = iparticle->GetFirstMother();
		  mother = (TParticle *) fParticles->At(imo);
		  imo = (mother->GetPdgCode() != 92) ? imo = newPos[imo] : -1;
	      } // if has mother   

	      Bool_t tFlag = (fTrackIt && (ks == 1));
	      
	      PushTrack(tFlag,imo,kf,p,origin,polar,tof,kPNoProcess,nt, 1., ks);
	      KeepTrack(nt);
	      newPos[i] = nt;
	  } // if selected
      } // particle loop
      delete[] newPos;
      delete[] pSelected;
      if (nc>0) {
	  jev += nc;
	  if (jev >= fNpart || fNpart == -1) {
	      break;
	  }
      }
  } // event loop
  MakeHeader();
  SetHighWaterMark(nt);
}

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
void AliGenDPMjet::MakeHeader()
{
// Builds the event header, to be called after each event
    AliGenEventHeader* header = new AliGenDPMjetEventHeader("DPMJET");
    ((AliGenDPMjetEventHeader*) header)->SetNProduced(fDPMjet->GetNumStablePc());
    ((AliGenDPMjetEventHeader*) header)->SetImpactParameter(fDPMjet->GetBImpac());
    ((AliGenDPMjetEventHeader*) header)->SetTotalEnergy(fDPMjet->GetTotEnergy());
    ((AliGenDPMjetEventHeader*) header)->SetParticipants(fDPMjet->GetfIp(), 
    							 fDPMjet->GetfIt());
 ((AliGenDPMjetEventHeader*) header)->SetProcessType(fDPMjet->GetProcessCode());
// Bookkeeping for kinematic bias
    ((AliGenDPMjetEventHeader*) header)->SetTrials(fTrials);
// Event Vertex
    header->SetPrimaryVertex(fVertex);
    gAlice->SetGenEventHeader(header);    
 AddHeader(header);
}

void AliGenDPMjet::AddHeader(AliGenEventHeader* header)
{
    // Add header to container or runloader
    if (fContainer) {
        fContainer->AddHeader(header);
    } else {
        AliRunLoader::GetRunLoader()->GetHeader()->SetGenEventHeader(header);
    }
}


//______________________________________________________________________________
AliGenDPMjet& AliGenDPMjet::operator=(const  AliGenDPMjet& /*rhs*/)
{
// Assignment operator
    return *this;
}



//______________________________________________________________________________
