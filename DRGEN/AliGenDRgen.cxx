/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
 *                                                                        *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////
// Event generator of two-pomeron processes
// in HE pp collisions. This is COMPASS generator adopted for ALICE exp.
// 
//%
// The example of the generator initialization for the process
// pp->pXp, X=f2(1270)
//%
//    AliGenDRgen *gener = new AliGenDRgen();
//    gener->SetProcess(4,4);
//    gener->SetBeamEnergy(3500.);
//    gener->Init();
//%
// Authors: Sergey.Evdokimov@cern.ch, Serguei.Sadovski@cern.ch
// 9 Oct 2013

#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>

//#include "AliRun.h"
#include "AliGenEventHeader.h"
#include "AliGenDRgen.h"
#include "TDRgen.h"
#include "DRGENcommon.h"

ClassImp(AliGenDRgen)

//------------------------------------------------------------

AliGenDRgen::AliGenDRgen() :
  AliGenMC(),
  fTDRgen(0x0),
  fParticles(0x0),
  fEvent(-1),
  fDebug(0),
  fDebugEventFirst(-1),
  fDebugEventLast(-1),
  kMassRange(kFALSE)
{
  // Constructor: create generator instance,
  // create particle array
  // Set DRgen parameters to default values:
  // f2(1270) production in PP collisions@7TeV

  SetMC(new TDRgen());
  fTDRgen = (TDRgen*) fMCEvGen;
  fParticles = new TClonesArray("TParticle",100);

  SetProcess   ();
  SetBeamEnergy();
  Double_t polarization[8]={0.,0.,0.,0.,0.,0.,0.,0.};
  SetF2Polarization(polarization);
}

//____________________________________________________________
AliGenDRgen::~AliGenDRgen()
{
  // Destroys the object, deletes and disposes all TParticles currently on list.
  if (fParticles) {
    fParticles->Delete();
    delete fParticles;
    fParticles = 0;
  }
}

//____________________________________________________________
void AliGenDRgen::Init()
{
  // Initialize the generator DRgen

  fTDRgen->Initialize();
  fEvent = 0;
}

//____________________________________________________________
void AliGenDRgen::Generate()
{
  // Generate one event of two-pomeron process.
  // Gaussian smearing on the vertex is done if selected. 
  // All particles from the DRgen/PYTHIA event listing
  // are stored in TreeK, and only final-state particles are tracked.
  // The event differectial cross section is assigned as a weight
  // to each track of the event.

  Float_t polar[3]   = {0,0,0};
  Float_t origin0[3] = {0,0,0};
  Float_t time0 = 0.;
  Float_t origin[3]  = {0,0,0};
  Float_t p[3], tof;
  Double_t weight;

  Int_t    ks,kf,iparent,nt, trackIt;
  Bool_t isWithinCuts;
  Float_t  random[6];
  const Float_t kconv=0.001/2.999792458e8;

  isWithinCuts = kFALSE;

  while (!isWithinCuts) {

    fTDRgen->GenerateEvent();
    weight = 1.;
    isWithinCuts=kTRUE;
    //  if (gAlice->GetEvNumber()>=fDebugEventFirst &&
    //      gAlice->GetEvNumber()<=fDebugEventLast) fPythia->Pylist(1);
    fTDRgen->ImportParticles(fParticles);

    if (fDebug == 1)
      Info("Generate()","one event before cuts is produced");
    Int_t np = fParticles->GetEntriesFast();
    TParticle *iparticle;
    //   Int_t* pParent   = new Int_t[np];
    for (Int_t ip=0; ip<np; ip++) {
      iparticle = (TParticle *) fParticles->At(ip);
      Int_t st = iparticle->GetStatusCode();
      if(st == 99){//this particle is central system
	if(kYRange)//cut on rapidity of central system  
	  if(fYMin > iparticle->Y() || fYMax <iparticle->Y() ) isWithinCuts = kFALSE;
	if(kPtRange)//cut on Pt of central system
	  if(fPtMin > iparticle->Pt() || fPtMax < iparticle->Pt() ) isWithinCuts = kFALSE;
	if(kMassRange){//cut on mass of central system
	  Float_t mass = TMath::Sqrt(TMath::Abs(iparticle->Energy()*iparticle->Energy() - iparticle->P()*iparticle->P()));
	  if(fMassMin > mass || fPtMax < mass ) isWithinCuts = kFALSE;
	}
      }

    }
  }

  Int_t j;
  for (j=0;j<3;j++) origin0[j]=fOrigin[j];
  time0 = fTimeOrigin;
  if(fVertexSmear==kPerEvent) {
    Rndm(random,6);
    for (j=0;j<3;j++) {
      origin0[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
    Rndm(random,2);
    time0 += fOsigma[2]/TMath::Ccgs()*
      TMath::Cos(2*random[0]*TMath::Pi())*
      TMath::Sqrt(-2*TMath::Log(random[1]));
  }

  Int_t ip;
  Int_t np = fParticles->GetEntriesFast();
  TParticle *iparticle;
//   Int_t* pParent   = new Int_t[np];
  for (ip=0; ip<np; ip++) {
    iparticle = (TParticle *) fParticles->At(ip);
    ks = iparticle->GetStatusCode();
//     // No initial state partons
//     if (ks==21) continue;
    p[0] = iparticle->Px();
    p[1] = iparticle->Py();
    p[2] = iparticle->Pz();
    origin[0] = origin0[0]+iparticle->Vx()/10.;
    origin[1] = origin0[1]+iparticle->Vy()/10.;
    origin[2] = origin0[2]+iparticle->Vz()/10.;
    kf = CheckPDGCode(iparticle->GetPdgCode());
    iparent = -1;
    tof     = time0 + kconv*iparticle->T();
    if (ks == 1) trackIt = 1;
    else         trackIt = 0;
    PushTrack(fTrackIt*trackIt,iparent,kf,p,origin,polar,tof,kPPrimary,nt,weight,ks);
    KeepTrack(nt); 

    if (fDebug == 2)
      printf("id=%+4d, parent=%3d, ks=%d, p = (%+11.4e,%+11.4e,%+11.4e) GeV\n",
	     kf,iparent,fTrackIt*trackIt,p[0],p[1],p[2]);
  }
  fEvent++;
  // fTDRgen->SetNEVENT(fEvent);
  if (fDebug == 1 && fEvent%100 == 0) {
    Info("Generate","Event %d\n",fEvent);
    fTDRgen->Finish();
  }

  // Make header
  AliGenEventHeader* header = new AliGenEventHeader("DRgen");
 
  // Pass header
  AddHeader(header);
}

//____________________________________________________________
void AliGenDRgen::SetEventListRange(Int_t eventFirst, Int_t eventLast)
{
  // Set a range of event numbers, for which a table
  // of generated particle will be printed
  fDebugEventFirst = eventFirst;
  fDebugEventLast  = eventLast;
  if (fDebugEventLast==-1) fDebugEventLast=fDebugEventFirst;
}

//____________________________________________________________
void AliGenDRgen::SetProcess       (Int_t   proc1, Int_t proc2  )
{
  // Set process number:

  fTDRgen->SetIPROC(proc1,proc2);
}
//____________________________________________________________
  void AliGenDRgen::SetBeamEnergy  (Double_t energy)
{
  // Set energy of the beam ion per nucleon in GeV
  fTDRgen->SetSqrtS(2.*energy);
}
//____________________________________________________________
  void AliGenDRgen::SetF2Polarization  (Double_t *polarization)
{
  // Set polarization of f2(1270) meson
  //f2 polarization array: |D0|^2, |D-|^2, |D+|^2, |D--|^2, |D++|^2, phase(D-,D0), phase(D--,D0), phase(D++,D+)
  if(polarization!=0x0)
    fTDRgen->SetF2Polarization(polarization);
}
//____________________________________________________________
  void AliGenDRgen::SetMassRange(Float_t min, Float_t max)
{
  if(min<max && max > 0){
    kMassRange = kTRUE;
    fMassMin = min;
    fMassMax = max;
  }
}
//____________________________________________________________
  TClonesArray*  AliGenDRgen::GetParticleList ()
{
  // Get array of particles of the event listing
  return fParticles;
}
//____________________________________________________________
