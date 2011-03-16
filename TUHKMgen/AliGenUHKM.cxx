/////////////////////////////////////////////////////////////////////////////
// Generator using UHKM 3.0 as an external generator.                      //
// ( only the HYDJET++ part is implemented for a moment)                   //
// temporary link:                                                         //
// http://lav01.sinp.msu.ru/~igor/hydjet++/hydjet++.txt                    //
// The main UHKM options are accessable  through this interface.           //
// Uses the TUHKMgen implementation of TGenerator.                         //
// Author of the first implementation: Sergey Zaporozhets                  //
// (zaporozh@sunhe.jinr.ru)                                                //
// Futhers modifications were made by                                      //
// Ionut Cristian Arsene (i.c.arsene@fys.uio.no)                           //
// & Malinina Liudmila(malinina@lav01.sinp.msu.ru) using as an example     //
//  AliGenTherminator.cxx created by Adam Kisiel                           //
//                                                                         //
////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include "TUHKMgen.h"
#ifndef DATABASE_PDG
#include "DatabasePDG.h"
#endif
#ifndef PARTICLE_PDG
#include "ParticlePDG.h"
#endif
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TMCProcess.h>
#include <TDatabasePDG.h>
#include <TSystem.h>

#include "AliGenUHKM.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliDecayer.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliLog.h"

using namespace std;


ClassImp(AliGenUHKM)

//_______________________________________
AliGenUHKM::AliGenUHKM()
  :AliGenMC(),
   fTrials(0),
   fUHKMgen(0),
   fHydjetParams(),
   fStableFlagged(0)
{
  // Default constructor setting up default reasonable parameter values
  // for central Pb+Pb collisions at 5.5TeV 
  
  // LHC
  fHydjetParams.fSqrtS=5500; //LHC
  fHydjetParams.fAw=207;//Pb-Pb
  fHydjetParams.fBmin=0.;
  fHydjetParams.fBmax=0.5; //0-5% centrality
  fHydjetParams.fT = 0.170;
  fHydjetParams.fMuB = 0.0;
  fHydjetParams.fMuS = 0.0;
  fHydjetParams.fMuI3 = 0.0;
  fHydjetParams.fThFO = 0.130;
  fHydjetParams.fMu_th_pip = 0.0;
  fHydjetParams.fSeed=0;
  fHydjetParams.fTau=10.;
  fHydjetParams.fSigmaTau=3.;
  fHydjetParams.fR=11.;
  fHydjetParams.fYlmax=4.0;
  fHydjetParams.fUmax=1.1;
  fHydjetParams.fDelta=0.;
  fHydjetParams.fEpsilon=0.;
  fHydjetParams.fWeakDecay=0; //>=0 on ,-1 off
  fHydjetParams.fEtaType=1;//gaus
  fHydjetParams.fCorrS=1.;
  fHydjetParams.fNhsel=2;
  fHydjetParams.fIshad=1;
  fHydjetParams.fPtmin=7.0;
  fHydjetParams.fT0=0.8;
  fHydjetParams.fTau0=0.1;
  fHydjetParams.fNf=0;
  fHydjetParams.fIenglu=0;
  fHydjetParams.fIanglu=0;


/* RHIC
  fHydjetParams.fSqrtS=200; //RHIC
  fHydjetParams.fAw=197;//Au-Au
  fHydjetParams.fBmin=0.;
  fHydjetParams.fBmax=0.5; //0-5% centrality
  fHydjetParams.fT = 0.165;
  fHydjetParams.fMuB = 0.0285;
  fHydjetParams.fMuS = 0.007;
  fHydjetParams.fMuI3 = -0.001;
  fHydjetParams.fThFO = 0.100;
  fHydjetParams.fMu_th_pip = 0.053;
  fHydjetParams.fSeed=0;
  fHydjetParams.fTau=8.;
  fHydjetParams.fSigmaTau=2.;
  fHydjetParams.fR=10.;
  fHydjetParams.fYlmax=3.3;
  fHydjetParams.fUmax=1.1;
  fHydjetParams.fDelta=0.;
  fHydjetParams.fEpsilon=0.;
  fHydjetParams.fWeakDecay=0; //>=0 on ,-1 off
  fHydjetParams.fEtaType=1;//gaus
  fHydjetParams.fCorrS=1.;
  fHydjetParams.fNhsel=2;
  fHydjetParams.fIshad=0;
  fHydjetParams.fPtmin=3.4;
  fHydjetParams.fT0=0.3;
  fHydjetParams.fTau0=0.4;
  fHydjetParams.fNf=2;
  fHydjetParams.fIenglu=0;
  fHydjetParams.fIanglu=0;
*/
  if(strlen(Form("%s/TUHKMgen/UHKM/particles.data", gSystem->Getenv("ALICE_ROOT")))<255)
    strncpy(fParticleFilename, Form("%s/TUHKMgen/UHKM/particles.data", gSystem->Getenv("ALICE_ROOT")), 256);
  if(strlen(Form("%s/TUHKMgen/UHKM/tabledecay.txt",gSystem->Getenv("ALICE_ROOT")))<255)
    strncpy(fDecayFilename, Form("%s/TUHKMgen/UHKM/tabledecay.txt", gSystem->Getenv("ALICE_ROOT")), 256);
  for(Int_t i=0; i<500; i++) {
    fStableFlagPDG[i] = 0;
    fStableFlagStatus[i] = kFALSE;
  }
  fStableFlagged = 0;

}

//_______________________________________
AliGenUHKM::AliGenUHKM(Int_t npart)
  :AliGenMC(npart),
   fTrials(0),
   fUHKMgen(0),
   fHydjetParams(),
   fStableFlagged(0)
{
  // Constructor specifying the size of the particle table
  // and setting up default reasonable parameter values
  // for central Pb+Pb collisions at 5.5TeV

  fName = "UHKM";
  fTitle= "Particle Generator using UHKM 3.0";

  fNprimaries = 0;   

  //LHC
  fHydjetParams.fSqrtS=5500; //LHC
  fHydjetParams.fAw=207;//Pb-Pb
  fHydjetParams.fBmin=0.;
  fHydjetParams.fBmax=0.5; //0-5% centrality
  fHydjetParams.fT = 0.170;
  fHydjetParams.fMuB = 0.0;
  fHydjetParams.fMuS = 0.0;
  fHydjetParams.fMuI3 = 0.0;
  fHydjetParams.fThFO = 0.130;
  fHydjetParams.fMu_th_pip = 0.0;
  fHydjetParams.fSeed=0;
  fHydjetParams.fTau=10.;
  fHydjetParams.fSigmaTau=3.;
  fHydjetParams.fR=11.;
  fHydjetParams.fYlmax=4.0;
  fHydjetParams.fUmax=1.1;
  fHydjetParams.fDelta=0.;
  fHydjetParams.fEpsilon=0.;
  fHydjetParams.fWeakDecay=0; //>=0 on ,-1 off
  fHydjetParams.fEtaType=1;//gaus
  fHydjetParams.fCorrS=1.;
  fHydjetParams.fNhsel=2;
  fHydjetParams.fIshad=1;
  fHydjetParams.fPtmin=7.0;
  fHydjetParams.fT0=0.8;
  fHydjetParams.fTau0=0.1;
  fHydjetParams.fNf=0;
  fHydjetParams.fIenglu=0;
  fHydjetParams.fIanglu=0;

/*RHIC
  fHydjetParams.fSqrtS=200; //RHIC
  fHydjetParams.fAw=197;//Au-Au
  fHydjetParams.fBmin=0.;
  fHydjetParams.fBmax=0.5; //0-5% centrality
  fHydjetParams.fT = 0.165;
  fHydjetParams.fMuB = 0.0285;
  fHydjetParams.fMuS = 0.007;
  fHydjetParams.fMuI3 = -0.001;
  fHydjetParams.fThFO = 0.100;
  fHydjetParams.fMu_th_pip = 0.053;
  fHydjetParams.fSeed=0;
  fHydjetParams.fTau=8.;
  fHydjetParams.fSigmaTau=2.;
  fHydjetParams.fR=10.;
  fHydjetParams.fYlmax=3.3;
  fHydjetParams.fUmax=1.1;
  fHydjetParams.fDelta=0.;
  fHydjetParams.fEpsilon=0.;
  fHydjetParams.fWeakDecay=0;//>=0 on ,-1 off
  fHydjetParams.fEtaType=1;//gaus
  fHydjetParams.fCorrS=1.;
  fHydjetParams.fNhsel=2;
  fHydjetParams.fIshad=1;
  fHydjetParams.fPtmin=3.4;
  fHydjetParams.fT0=0.3;
  fHydjetParams.fTau0=0.4;
  fHydjetParams.fNf=2;
  fHydjetParams.fIenglu=0;
  fHydjetParams.fIanglu=0;
*/

  if(strlen(Form("%s/TUHKMgen/UHKM/particles.data", gSystem->Getenv("ALICE_ROOT")))<255)
    strncpy(fParticleFilename, Form("%s/TUHKMgen/UHKM/particles.data", gSystem->Getenv("ALICE_ROOT")), 256);
  if(strlen(Form("%s/TUHKMgen/UHKM/tabledecay.txt", gSystem->Getenv("ALICE_ROOT")))<255)
    strncpy(fDecayFilename, Form("%s/TUHKMgen/UHKM/tabledecay.txt", gSystem->Getenv("ALICE_ROOT")), 256);
  for(Int_t i=0; i<500; i++) {
    fStableFlagPDG[i] = 0;
    fStableFlagStatus[i] = kFALSE;
  }
  fStableFlagged = 0;  

}

//__________________________________________
AliGenUHKM::~AliGenUHKM()
{
  // Destructor, do nothing
  //  delete fParticles;
}

void AliGenUHKM::SetAllParametersRHIC() 
{
  // Set reasonable default parameters for 0-5% central Au+Au collisions
  // at 200 GeV at RHIC
  SetEcms(200.0);                // RHIC top energy
  SetAw(197);                    // Au+Au
  SetBmin(0.0);                  // 0%
  SetBmax(0.5);                  // 5%
  SetChFrzTemperature(0.165);    // T_ch = 165 MeV
  SetMuB(0.0285);                // mu_B = 28.5 MeV
  SetMuS(0.007);                 // mu_S = 7 MeV
  SetMuQ(-0.001);                // mu_Q = -1 MeV
  SetThFrzTemperature(0.100);    // T_th = 100 MeV
  SetMuPionThermal(0.053);       // mu_th_pion = 53 MeV 
  SetSeed(0);                    // use UNIX time
  SetTauB(8.0);                  // tau = 8 fm/c
  SetSigmaTau(2.0);              // sigma_tau = 2 fm/c
  SetRmaxB(10.0);                // fR = 10 fm
  SetYlMax(3.3);                 // fYmax = 3.3
  SetEtaRMax(1.1);               // Umax = 1.1
  SetMomAsymmPar(0.0);           // delta = 0.0
  SetCoordAsymmPar(0.0);         // epsilon = 0.0
  //  SetFlagWeakDecay(0);           // weak decay on (<0 off !!!)
  SetEtaType(1);                 // gaus distributed with fYmax dispersion (0 means boost invariant)
  SetGammaS(1.0);                // gammaS = 1.0 (no strangeness canonical suppresion)
  SetPyquenNhsel(2);             // hydro on, jets on, jet quenching on
  SetPyquenShad(1);              // shadowing on (0 off)
  SetPyquenPtmin(3.4);           // ptmin = 3.4 GeV/c
  SetPyquenT0(0.3);              // T0 = 300 MeV
  SetPyquenTau0(0.4);            // tau0 = 0.4 fm/c
  SetPyquenNf(2);                // 2 flavours
  SetPyquenIenglu(0);            // radiative and collisional energy loss
  SetPyquenIanglu(0);            // small gluon angular distribution
}

void AliGenUHKM::SetAllParametersLHC()
{
  // Set reasonable default parameters for 0-5% central Pb+Pb collisions
  // at 5.5 TeV at LHC
  SetEcms(5500.0);               // LHC
  SetAw(207);                    // Pb+Pb
  SetBmin(0.0);                  // 0%
  SetBmax(0.5);                  // 5%
  SetChFrzTemperature(0.170);    // T_ch = 170 MeV
  SetMuB(0.0);                   // mu_B = 0 MeV
  SetMuS(0.0);                   // mu_S = 0 MeV
  SetMuQ(0.0);                   // mu_Q = 0 MeV
  SetThFrzTemperature(0.130);    // T_th = 130 MeV
  SetMuPionThermal(0.0);         // mu_th_pion = 0 MeV
  SetSeed(0);                    // use UNIX time
  SetTauB(10.0);                 // tau = 10 fm/c
  SetSigmaTau(3.0);              // sigma_tau = 3 fm/c
  SetRmaxB(11.0);                // fR = 11 fm
  SetYlMax(4.0);                 // fYmax = 4.0
  SetEtaRMax(1.1);               // Umax = 1.1
  SetMomAsymmPar(0.0);           // delta = 0.0
  SetCoordAsymmPar(0.0);         // epsilon = 0.0
  //  SetFlagWeakDecay(0);           // weak decay on (<0 off !!!)
  SetEtaType(1);                 // gaus distributed with fYmax dispersion (0 means boost invariant)
  SetGammaS(1.0);                // gammaS = 1.0 (no strangeness canonical suppresion)
  SetPyquenNhsel(2);             // hydro on, jets on, jet quenching on
  SetPyquenShad(1);              // shadowing on (0 off)
  SetPyquenPtmin(7.0);           // ptmin = 7.0 GeV/c
  SetPyquenT0(0.8);              // T0 = 800 MeV
  SetPyquenTau0(0.1);            // tau0 = 0.4 fm/c
  SetPyquenNf(0);                // 0 flavours
  SetPyquenIenglu(0);            // radiative and collisional energy loss
  SetPyquenIanglu(0);            // small gluon angular distribution
}

//_________________________________________
void AliGenUHKM::Init()
{
  // Initialization of the TGenerator::TUHKMgen interface object. 
  // Model input parameters are transmited to the TUHKMgen object which forwards them
  // further to the model.
  // HYDJET++ is initialized (average multiplicities are calculated, particle species definitions and decay
  // channels are loaded, etc.)

  SetMC(new TUHKMgen());
  fUHKMgen = (TUHKMgen*) fMCEvGen;
  SetAllParameters();
  
  AliGenMC::Init();
  fUHKMgen->Initialize();
  CheckPDGTable();    
  fUHKMgen->Print();  
}

//________________________________________
void AliGenUHKM::Generate()
{
  // Generate one HYDJET++ event, get the output and push particles further 
  // to AliRoot's stack 

  Float_t polar[3] = {0,0,0};
  Float_t origin[3]   = {0,0,0};
  Float_t origin0[3]  = {0,0,0};
  Float_t p[3];
  Float_t v[3];
  Float_t mass=0.0, energy=0.0;

  Vertex();
  for(Int_t j=0; j<3; j++) origin0[j] = fVertex[j];

  // Generate the event and import particles
  fUHKMgen->GenerateEvent();
  fUHKMgen->ImportParticles(&fParticles,"All");
  Int_t np = fParticles.GetEntriesFast();
  Int_t nt  =  0;


  // Handle the IDs of particles on the stack  
  Int_t* idsOnStack = new Int_t[np];
  Int_t* newPos     = new Int_t[np];
  for(Int_t i=0; i<np; i++) {
    newPos[i] = i;
    idsOnStack[i] = -1;
  }
  

  // Generate a random phi used to rotate the whole event
  Double_t eventRotation = gRandom->Rndm()*TMath::Pi();

  TParticle *iparticle;
  Double_t partMomPhi=0.0, partPt=0.0;
  Double_t partVtxPhi=0.0, partVtxR=0.0;
  //_________ Loop for particles in the stack
  for(Int_t i=0; i<np; i++) {
    iparticle = (TParticle*)fParticles.At(i);
    Int_t kf = iparticle->GetPdgCode();
    Bool_t hasMother = (iparticle->GetFirstMother() >= 0);
    Bool_t hasDaughter = (iparticle->GetNDaughters() > 0);
    
    if(hasDaughter) {
      // This particle has decayed
      // It will not be tracked
      // Add it only once with coordinates not
      // smeared with primary vertex position
      
      // rotate the direction of the particle
      partMomPhi = TMath::ATan2(iparticle->Py(), iparticle->Px());
      partPt = TMath::Hypot(iparticle->Px(), iparticle->Py());
      p[0] = partPt*TMath::Cos(partMomPhi+eventRotation);
      p[1] = partPt*TMath::Sin(partMomPhi+eventRotation);
      p[2] = iparticle->Pz();
      
      mass = TDatabasePDG::Instance()->GetParticle(kf)->Mass();
      energy = sqrt(mass*mass + p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

      // rotate the freezeout point
      partVtxPhi = TMath::ATan2(iparticle->Vy(), iparticle->Vx());
      partVtxR = TMath::Hypot(iparticle->Vx(), iparticle->Vy());
      v[0] = partVtxR*TMath::Cos(partVtxPhi + eventRotation);
      v[1] = partVtxR*TMath::Cos(partVtxPhi + eventRotation);
      v[2] = iparticle->Vz();

      Float_t time = iparticle->T();
      
      Int_t imo = -1;
      if(hasMother) {
	imo = iparticle->GetFirstMother(); //index of mother particle in fParticles
      } // if has mother
      Bool_t trackFlag = kFALSE;   // tFlag is kFALSE --> do not track the particle
      
      PushTrack(trackFlag, (imo>=0 ? idsOnStack[imo] : imo), kf,
		p[0], p[1], p[2], energy,
		v[0], v[1], v[2], time,
		polar[0], polar[1], polar[2],
		(hasMother ? kPDecay : kPNoProcess), nt);
      idsOnStack[i] = nt;
      
      fNprimaries++;
      KeepTrack(nt);
    }
    else {
      // This is a final state particle
      // It will be tracked
      // Add it TWICE to the stack !!!
      // First time with event-wide coordinates (for femtoscopy) -
      //   this one will not be tracked
      // Second time with event-wide c0ordinates and vertex smearing
      //   this one will be tracked

      // rotate the direction of the particle      
      partMomPhi = TMath::ATan2(iparticle->Py(), iparticle->Px());
      partPt = TMath::Hypot(iparticle->Px(), iparticle->Py());
      p[0] = partPt*TMath::Cos(partMomPhi+eventRotation);
      p[1] = partPt*TMath::Sin(partMomPhi+eventRotation);
      p[2] = iparticle->Pz();      

      energy = sqrt(mass*mass + p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

      // rotate the freezeout point
      partVtxPhi = TMath::ATan2(iparticle->Vy(), iparticle->Vx());
      partVtxR = TMath::Hypot(iparticle->Vx(), iparticle->Vy());
      v[0] = partVtxR*TMath::Cos(partVtxPhi + eventRotation);
      v[1] = partVtxR*TMath::Cos(partVtxPhi + eventRotation);
      v[2] = iparticle->Vz();      
      
      Int_t type    = iparticle->GetStatusCode(); // 1-from jet / 0-from hydro 
      Int_t coeffT=1;
      if(type==1) coeffT=-1; //to separate particles from jets
      
      
      Int_t imo = -1;
      
      if(hasMother) {
	imo = iparticle->GetFirstMother();
      } // if has mother
      
      Bool_t trackFlag = kFALSE;  // tFlag = kFALSE --> do not track this one, its for femtoscopy
      PushTrack(trackFlag, (imo>=0 ? idsOnStack[imo] : imo), kf,
		p[0], p[1], p[2], energy,
		v[0], v[1], v[2], (iparticle->T())*coeffT,    // freeze-out time is negative if the particle comes from jet
		polar[0], polar[1], polar[2],
		hasMother ? kPDecay:kPNoProcess, nt);
      
      idsOnStack[i] = nt;
      fNprimaries++;
      KeepTrack(nt);
      
      origin[0] = origin0[0]+v[0];
      origin[1] = origin0[1]+v[1];
      origin[2] = origin0[2]+v[2];
      imo = nt;
      
      trackFlag = fTrackIt;    // Track this particle, unless otherwise specified by fTrackIt
      
      PushTrack(trackFlag, imo, kf,
		p[0], p[1], p[2], energy,
		origin[0], origin[1], origin[2], iparticle->T(),
		polar[0], polar[1], polar[2],
		hasMother ? kPDecay:kPNoProcess, nt);
      
      fNprimaries++;
      KeepTrack(nt);
    }
  }
  
  SetHighWaterMark(fNprimaries);

  TArrayF eventVertex;
  eventVertex.Set(3);
  eventVertex[0] = origin0[0];
  eventVertex[1] = origin0[1];
  eventVertex[2] = origin0[2];

  // Builds the event header, to be called after each event
  AliGenEventHeader* header = new AliGenHijingEventHeader("UHKM");
  Double_t b       = 0.;
  Double_t npart   = 0; 
  Double_t nbin    = 0;
  fUHKMgen->GetCentrality(b, npart, nbin);
  printf("********** %13.3f %13.3f %13.3f \n", b, npart, nbin);
  


  ((AliGenHijingEventHeader*) header)->SetNProduced(fNprimaries);
  ((AliGenHijingEventHeader*) header)->SetPrimaryVertex(eventVertex);
  ((AliGenHijingEventHeader*) header)->SetImpactParameter(b);
  ((AliGenHijingEventHeader*) header)->SetTotalEnergy(0.0);
  ((AliGenHijingEventHeader*) header)->SetHardScatters(0);
  ((AliGenHijingEventHeader*) header)->SetParticipants(Int_t(npart), 0);
  ((AliGenHijingEventHeader*) header)->SetCollisions(Int_t(nbin), 0, 0, 0);
  ((AliGenHijingEventHeader*) header)->SetSpectators(0, 0, 0, 0);
  ((AliGenHijingEventHeader*) header)->SetReactionPlaneAngle(0);//evrot);

  header->SetPrimaryVertex(fVertex);
  AddHeader(header);
  fCollisionGeometry = (AliGenHijingEventHeader*)  header;

  delete [] idsOnStack;

}

void AliGenUHKM::Copy(TObject &) const
{
  Fatal("Copy","Not implemented!\n");
}

void AliGenUHKM::SetAllParameters() {
  // Forward all input parameters to the TGenerator::TUHKMgen object

  fUHKMgen->SetEcms(fHydjetParams.fSqrtS);
  fUHKMgen->SetBmin(fHydjetParams.fBmin);
  fUHKMgen->SetBmax(fHydjetParams.fBmax);
  fUHKMgen->SetAw(fHydjetParams.fAw);
  fUHKMgen->SetSeed(fHydjetParams.fSeed);

  fUHKMgen->SetChFrzTemperature(fHydjetParams.fT);
  fUHKMgen->SetMuB(fHydjetParams.fMuB);
  fUHKMgen->SetMuS(fHydjetParams.fMuS);
  fUHKMgen->SetMuQ(fHydjetParams.fMuI3);
  fUHKMgen->SetTauB(fHydjetParams.fTau);
  fUHKMgen->SetThFrzTemperature(fHydjetParams.fThFO);
  fUHKMgen->SetMuPionThermal(fHydjetParams.fMu_th_pip);

  fUHKMgen->SetSigmaTau(fHydjetParams.fSigmaTau);
  fUHKMgen->SetRmaxB(fHydjetParams.fR);
  fUHKMgen->SetYlMax(fHydjetParams.fYlmax);
  fUHKMgen->SetEtaRMax(fHydjetParams.fUmax);
  fUHKMgen->SetMomAsymmPar(fHydjetParams.fDelta);
  fUHKMgen->SetCoordAsymmPar(fHydjetParams.fEpsilon);

  fUHKMgen->SetGammaS(fHydjetParams.fCorrS);
  fUHKMgen->SetEtaType(fHydjetParams.fEtaType);
  fUHKMgen->SetFlagWeakDecay(fHydjetParams.fWeakDecay);

  //PYQUEN parameters

  fUHKMgen->SetPyquenNhsel(fHydjetParams.fNhsel);
  fUHKMgen->SetPyquenShad(fHydjetParams.fIshad);
  fUHKMgen->SetPyquenPtmin(fHydjetParams.fPtmin);
  fUHKMgen->SetPyquenT0(fHydjetParams.fT0);
  fUHKMgen->SetPyquenTau0(fHydjetParams.fTau0);
  fUHKMgen->SetPyquenNf(fHydjetParams.fNf);
  fUHKMgen->SetPyquenIenglu(fHydjetParams.fIenglu);
  fUHKMgen->SetPyquenIanglu(fHydjetParams.fIanglu);

  fUHKMgen->SetPDGParticleFile(fParticleFilename);
  fUHKMgen->SetPDGDecayFile(fDecayFilename);
  //  fUHKMgen->SetUseCharmParticles(fUseCharmParticles);
  //  fUHKMgen->SetMinimumWidth(fMinWidth);
  //  fUHKMgen->SetMaximumWidth(fMaxWidth);
  //  fUHKMgen->SetMinimumMass(fMinMass);
  //  fUHKMgen->SetMaximumMass(fMaxMass);

  //  cout << "AliGenUHKM::Init() no. stable flagged particles = " << fStableFlagged << endl;
  for(Int_t i=0; i<fStableFlagged; i++) {
    //    cout << "AliGenUHKM::Init() flag no. " << i
    //         << " PDG = " << fStableFlagPDG[i]
    //         << " flag = " << fStableFlagStatus[i] << endl;
    fUHKMgen->SetPDGParticleStable(fStableFlagPDG[i], fStableFlagStatus[i]);
  }

   cout<<" Print all parameters "<<endl;
   cout<<" SqrtS = "<<fHydjetParams.fSqrtS<<endl;
   cout<<" Bmin  = "<< fHydjetParams.fBmin<<endl;
 cout<<" Bmax= "<<fHydjetParams.fBmax<<endl;
 cout<<" Aw= "<<fHydjetParams.fAw<<endl;
 cout<<" Seed= "<<fHydjetParams.fSeed<<endl;

 cout<<" ---Stat-model parameters----------- "<<endl;

 cout<<" ChFrzTemperature= "<<fHydjetParams.fT<<endl;
 cout<<" MuB= "<<fHydjetParams.fMuB<<endl;
 cout<<" MuS= "<<fHydjetParams.fMuS<<endl;
 cout<<" MuQ= "<<fHydjetParams.fMuI3<<endl;
 cout<<" TauB= "<<fHydjetParams.fTau<<endl;
 cout<<" ThFrzTemperature= "<<fHydjetParams.fThFO<<endl;
 cout<<" MuPionThermal= "<<fHydjetParams.fMu_th_pip<<endl;

cout<<"-----Volume parameters -------------- "<<endl;

 cout<<" SigmaTau= "<<fHydjetParams.fSigmaTau<<endl;
 cout<<" RmaxB= "<<fHydjetParams.fR<<endl;
 cout<<" YlMax= "<<fHydjetParams.fYlmax<<endl;
 cout<<" EtaRMax= "<<fHydjetParams.fUmax<<endl;
 cout<<" MomAsymmPar= "<<fHydjetParams.fDelta<<endl;
 cout<<" CoordAsymmPar= "<<fHydjetParams.fEpsilon<<endl;

cout<<" --------Flags------ "<<endl;

 cout<<" GammaS= "<<fHydjetParams.fCorrS<<endl;
 cout<<" EtaType= "<<fHydjetParams.fEtaType<<endl;
 cout<<" FlagWeakDecay= "<<fHydjetParams.fWeakDecay<<endl;

  cout<<"----PYQUEN parameters---"<<endl;

  cout<<" Nhsel= "<<fHydjetParams.fNhsel<<endl;
  cout<<" Shad= "<<fHydjetParams.fIshad<<endl;
  cout<<" Ptmin= "<<fHydjetParams.fPtmin<<endl;
  cout<<" T0= "<<fHydjetParams.fT0<<endl;
  cout<<" Tau0= "<<fHydjetParams.fTau0<<endl;
  cout<<" Nf= "<<fHydjetParams.fNf<<endl;
  cout<<" Ienglu= "<<fHydjetParams.fIenglu<<endl;
  cout<<" Ianglu= "<<fHydjetParams.fIanglu<<endl;

  //  cout<<"----PDG table parameters---"<<endl;
  
  //  cout<<" UseCharmParticles= "<<fUseCharmParticles<<endl;
  //  cout<<" MinimumWidth= "<<fMinWidth<<endl;
  //  cout<<" MaximumWidth= "<<fMaxWidth<<endl;
  //  cout<<" MinimumMass= "<<fMinMass<<endl;
  //  cout<<" MaximumMass= "<<fMaxMass<<endl;

  //  cout << "AliGenUHKM::SetAllParameters() OUT" << endl;
}

// add the additional PDG codes from UHKM(SHARE table) to ROOT's table
void AliGenUHKM::CheckPDGTable() {
  // Add temporarely all particle definitions from HYDJET++ which miss in the ROOT's PDG tables
  // to the TDatabasePDG table.

  DatabasePDG *uhkmPDG = fUHKMgen->PDGInfo();              // UHKM's PDG table
  TParticlePDG *rootTestParticle;
  ParticlePDG *uhkmTestParticle;
  
  // loop over all particles in the SHARE table
  for(Int_t i=0; i<uhkmPDG->GetNParticles(); i++) {
    // get a particle specie
    uhkmTestParticle = uhkmPDG->GetPDGParticleByIndex(i);
    // check if this code exists in ROOT's table
    rootTestParticle = TDatabasePDG::Instance()->GetParticle(uhkmTestParticle->GetPDG());
    if(!rootTestParticle) {    // if not then add it to the ROOT's PDG database
      
      TDatabasePDG::Instance()->AddParticle(uhkmTestParticle->GetName(), uhkmTestParticle->GetName(), 
					    uhkmTestParticle->GetMass(), uhkmTestParticle->GetElectricCharge(),
					    (uhkmTestParticle->GetWidth()<1e-10 ? kTRUE : kFALSE),
					    uhkmTestParticle->GetWidth(), 
					    (Int_t(uhkmTestParticle->GetBaryonNumber())==0 ? "meson" : "baryon"),
					    uhkmTestParticle->GetPDG());    
      if(uhkmTestParticle->GetWidth()<1e-10) 
	cout << uhkmTestParticle->GetPDG() << " with mass "
	     << TDatabasePDG::Instance()->GetParticle(uhkmTestParticle->GetPDG())->Mass()
	     << TDatabasePDG::Instance()->GetParticle(uhkmTestParticle->GetPDG())->Width() << endl;  
    }
  }  // end for
}
