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

// Generator using Hydjet as an external generator
// The main Hydjet options are accessable for the user through this interface.
// Uses the THydjet implementation of TGenerator.
// Author:
// Rafael Diaz Valdes    (rafael.diaz.valdes@cern.ch)
//
#include <Riostream.h>
#include <THydjet.h>
#include <TParticle.h>
#include <TClonesArray.h>

#include "AliGenHydjet.h"
#include "AliLog.h"
#include "AliGenHydjetEventHeader.h"
#include "AliRun.h"
#include "AliPythiaRndm.h"

ClassImp(AliGenHydjet)

AliGenHydjet::AliGenHydjet(Int_t npart) :
  AliGenMC(npart),
    //initial parameters
  fFrame("CMS"),      // Reference frame
  fAtomicWeigth(207),        // Projectile-Target atomic weight
  fIfbtype(0),          // centrality type
  fFixImpactParam(0),          // fixed impact parameter
  fMinImpactParam(0.),         // minimum impact parameter
  fMaxImpactParam(1.),         // maximum impact parameter
  fSoftMult(20000),      // mean soft multiplicity
  //hydro parameters
  fJetProd(2),          // flag Jet production (nhsel)
  fYflow(5.),         // max longitudinal flow
  fTflow(1.),         // max transverse flow
  fSoftFract(1.),         // Soft multiplicity fraction
  // PYTHIA parameters
  fDijetProd(1),          // flag dijet production
  fMinPtHard(10.),        // min pt hard
  fStructFunction(7),          // Structure Function (default CTEQ5M)
  fMultipleInt(0),          // flag multiple interaction

  fHydjet(0),

  fSelectAll(0),
  fFlavor(0)
    
{
// Default PbPb collisions at 5. 5 TeV
//
  fEnergyCMS = 5500.;      //Energy cms
  fName = "Hydjet";
  fTitle = "Particle Generator using Hydjet";
  // Set random number generator 
  if (!AliPythiaRndm::GetPythiaRandom()) 
    AliPythiaRndm::SetPythiaRandom(GetRandom());
}

AliGenHydjet::~AliGenHydjet()
{
// Destructor
}

void AliGenHydjet::Init()
{
// Initialisation

    SetMC(new THydjet(fEnergyCMS, fFrame,fAtomicWeigth,
                      fIfbtype,fMinImpactParam, fMaxImpactParam,
                      fFixImpactParam,fSoftMult));

    fHydjet=(THydjet*) fMCEvGen;

    // init hydro common blocks
    fHydjet->SetNHSEL(fJetProd);
    fHydjet->SetYLFL(fYflow);
    fHydjet->SetYTFL(fTflow);
    fHydjet->SetFPART(fSoftFract);
    // init PYTHIA common block
    fHydjet->SetMSEL(fDijetProd);
    fHydjet->SetPTMIN(fMinPtHard);
    fHydjet->SetCKIN(3,fMinPtHard);
    fHydjet->SetMSTP(51,fStructFunction);
    fHydjet->SetMSTP(81,0);
    fHydjet->SetMSTU(21,1); 
    fHydjet->SetPARU(14,1.);
    
//  Initialize Hydjet  
    fHydjet->Initialize();

}

void AliGenHydjet::Generate()
{

///////////////////////////////////////////////////////////////////////////////////////////
  Float_t polar[3]    =   {0,0,0};
  Float_t origin[3]   =   {0,0,0};
  Float_t origin0[3]  =   {0,0,0};
  Float_t time0 = 0.;
  Float_t p[3];
  Float_t tof;

//  converts from mm/c to s
  const Float_t kconv = 0.001/2.999792458e8;
//
  Int_t nt  = 0;
  Int_t j, kf, ks, imo;
  kf = 0;

  for (j = 0;j < 3; j++) origin0[j] = fOrigin[j];
  time0 = fTimeOrigin;
  if(fVertexSmear == kPerEvent) {
      Vertex();
      for (j=0; j < 3; j++) origin0[j] = fVertex[j];
      time0 = fTime;
  }

  while(1)
  {
//    Generate one event
      fHydjet->GenerateEvent();
      fHydjet->ImportParticles(&fParticles,"All");

      Int_t np = fParticles.GetEntriesFast();
      Int_t nc = 0;
      if (np == 0 ) continue;
      Int_t i;
      Int_t* pSelected  = new Int_t[np];

      for (i = 0; i < np; i++) {
           pSelected[i] = 0;
      }

//      Get event vertex
      fVertex[0] = origin0[0];
      fVertex[1] = origin0[1];
      fVertex[2] = origin0[2];
      fTime = time0;

//
// Now select the final state particles
//

      for(i = 0; i<np; i++){
         TParticle *  iparticle = (TParticle *) fParticles.At(i);
         // Is this a final state particle ?
         if (!Stable(iparticle)) continue;
         Bool_t selected = kTRUE;
         kf = iparticle->GetPdgCode();
         ks = iparticle->GetStatusCode();
         if(!fSelectAll) selected = KinematicSelection(iparticle,0)&&SelectFlavor(kf);
         // Put particle on the stack if selected
         if(selected) {
           nc++;
           pSelected[i] = 1;
         } // selected
      } // particle loop final state


//
// Write particles to stack
//

      for(i = 0; i<np; i++) {
         TParticle *  iparticle = (TParticle *) fParticles.At(i);
         Bool_t  hasDaughter = (iparticle->GetFirstDaughter() > 0);
         if(pSelected[i]){
           kf = iparticle->GetPdgCode();
           ks = iparticle->GetStatusCode();
           p[0] = iparticle->Px();
           p[1] = iparticle->Py();
           p[2] = iparticle->Pz();
           origin[0] = fVertex[0]+iparticle->Vx()/10;
           origin[1] = fVertex[1]+iparticle->Vy()/10;
           origin[2] = fVertex[2]+iparticle->Vz()/10;
           tof   = fTime + kconv*iparticle->T();

           imo = -1;
           //imo = iparticle->GetFirstMother();
           //if(imo == 0) imo =-1;
           Bool_t tFlag = (fTrackIt && !hasDaughter);
           PushTrack(tFlag,imo,kf,p,origin,polar,tof,kPNoProcess,nt, 1., ks);
           KeepTrack(nt);
         } // if selected
      } // particle loop
      delete[] pSelected;

      AliInfo(Form("\n I've put %i particles on the stack \n",nc));
      if (nc > 0) break;
  } // event loop
  MakeHeader();
}


Bool_t AliGenHydjet::SelectFlavor(Int_t pid)
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
    return res;
}

Bool_t AliGenHydjet::Stable(TParticle*  particle) const
{
// Return true for a stable particle
    if (particle->GetFirstDaughter() == 0 )
    {
	return kTRUE;
    } else {
	return kFALSE;
    }
}


///////////////////////////////////
void AliGenHydjet::MakeHeader()
{
// Builds the event header, to be called after each event
    AliGenEventHeader* header = new AliGenHydjetEventHeader("Hydjet");
    ((AliGenHydjetEventHeader*) header)->SetNProduced(fHydjet->GetN());
    ((AliGenHydjetEventHeader*) header)->SetImpactParameter(fHydjet->GetBGEN());
    ((AliGenHydjetEventHeader*) header)->SetNbcol(fHydjet->GetNBCOL());
    ((AliGenHydjetEventHeader*) header)->SetNpart(fHydjet->GetNPART());
    ((AliGenHydjetEventHeader*) header)->SetNpyt(fHydjet->GetNPYT());
    ((AliGenHydjetEventHeader*) header)->SetNhyd(fHydjet->GetNHYD());
  // draw generation values
  cout << "*********Hyflow ********** " << endl;
  cout << " YTFL " << fHydjet->GetYTFL() << endl;
  cout << " YLFL " << fHydjet->GetYLFL() << endl;
  cout << " FPART " << fHydjet->GetFPART() << endl;
  cout << "*********Hyjpar ********* " << endl;
  cout << " GetNHSEL() " << fHydjet->GetNHSEL() << endl;
  cout << " GetPTMIN() " << fHydjet->GetPTMIN() << endl;
  cout << " GetNJET() " << fHydjet->GetNJET() << endl;
  cout << "*********Hyfpar ********* " << endl;
  cout << " GetBGEN() " << fHydjet->GetBGEN() << endl;
  cout << " GetNBCOL() " << fHydjet->GetNBCOL() << endl;
  cout << " GetNPART() " << fHydjet->GetNPART() << endl;
  cout << " GetNPYT() " << fHydjet->GetNPYT() << endl;
  cout << " GetNHYD() " << fHydjet->GetNHYD() << endl;
  

// Event Vertex
    header->SetPrimaryVertex(fVertex);
    header->SetInteractionTime(fTime);
    AddHeader(header);
    fCollisionGeometry = (AliGenHydjetEventHeader*)  header;
}

void AliGenHydjet::AddHeader(AliGenEventHeader* header)
{
    // Passes header either to the container or to gAlice
    if (fContainer) {
	fContainer->AddHeader(header);
    } else {
	gAlice->SetGenEventHeader(header);
    }
}


void AliGenHydjet::Copy(TObject &) const
{
  Fatal("Copy","Not implemented!\n");
}

AliGenHydjet& AliGenHydjet::operator=(const  AliGenHydjet& rhs)
{
    rhs.Copy(*this); 
    return (*this);
}

