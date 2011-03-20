/////////////////////////////////////////////////////////////////////////////
// TUHKM is an interface to 
// UHKM 3.0                                                                //
// ( only the HYDJET++ part is implemented)                                // 
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

#ifndef TUHKMGEN_H
#include "TUHKMgen.h"
#endif
#include "TObjArray.h"
#include "TParticle.h"
#include "Particle.h"
#include "TROOT.h"
#include "TError.h"
#include "TClonesArray.h"
#include "TSystem.h"

#include <string>
#include <iostream>
using namespace std;

//class TObjArray;

ClassImp(TUHKMgen)

TUHKMgen::TUHKMgen() : 
  TGenerator("UHKM","UHKM"),
  fInitialState(0x0),
  fAllocator(),
  fSecondariesList(),
  fNPprim(0),
  fNPsec(0),
  fHydjetParams(),
  fStableFlagged(0)
{
  // default constructor setting reasonable defaults for initial parameters (central Pb+Pb at 5.5 TeV)

  //  ParticleAllocator fAllocator;
  //  List_t fSecondariesList;

  // Set reasonable default values for LHC
  
  fHydjetParams.fSqrtS=5500; //LHC
  fHydjetParams.fAw=207;//Au-Au
  fHydjetParams.fIfb    = 1;      // flag of type of centrality generation (=0 is fixed by fBfix, not 0 
                                   //impact parameter is generated in each event between fBfmin 
                                   //and fBmax according with Glauber model (f-la 30)
  fHydjetParams.fBfix=0.;
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


  // Set reasonable default values for RHIC
/*
  fHydjetParams.fSqrtS=200; //RHIC          
  fHydjetParams.fAw=197;//Au-Au        
  fHydjetParams.fIfb    = 1;      
  fHydjetParams.fBfix=0.;       
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

  strncpy(fParticleFilename, Form("%s/TUHKMgen/UHKM/particles.data", gSystem->Getenv("ALICE_ROOT")), 255);
  strncpy(fDecayFilename, Form("%s/TUHKMgen/UHKM/tabledecay.txt", gSystem->Getenv("ALICE_ROOT")), 255);
  for(Int_t i=0; i<500; i++) {
    fStableFlagPDG[i] = 0;
    fStableFlagStatus[i] = kFALSE;
  }
  fStableFlagged = 0;
  //  cout << "TUHKMgen::TUHKMgen() OUT" << endl;
}

//______________________________________________________________________________
TUHKMgen::~TUHKMgen()
{
  // destructor, deletes the InitialStateHydjet object
  if(fInitialState)
    delete fInitialState;
}

void TUHKMgen::SetAllParametersRHIC()
{
  // Set reasonable input parameters for 0-5% central Au+Au collisions
  // at 200 GeV at RHIC
  SetEcms(200.0);                // RHIC top energy
  SetAw(197);                    // Au+Au
  SetIfb(1);
  SetBfix(0.);
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
  SetFlagWeakDecay(0.0);           // weak decay on (<0 off !!!)
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

void TUHKMgen::SetAllParametersLHC()
{
  // Set reasonable input parameters for 0-5% Pb+Pb collisions at 5.5 TeV at LHC
  SetEcms(5500.0);               // LHC
  SetAw(207);                    // Pb+Pb
  SetIfb(1);
  SetBfix(0.);                  // 0
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
  SetFlagWeakDecay(0);           // weak decay on (<0 off !!!)
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

TObjArray* TUHKMgen::ImportParticles(const Option_t *)
{
  // Function overloading the TGenerator::ImportParticles() member function.
  // The particles from the local particle list (fSecondariesList) are
  // forwarded to the TGenerator::fParticles 
  fParticles->Clear();
  Int_t nump = 0;
  LPIT_t it,e;
   
  for(it = fSecondariesList.begin(), e = fSecondariesList.end(); it != e; it++) {
    TVector3 pos(it->Pos().Vect());
    TVector3 mom(it->Mom().Vect());
    Float_t m1 = it->TableMass();
    Int_t im1 = it->GetMother();
    Int_t im2 = -1;
    Int_t id1 = -1;
    Int_t id2 = -1;

    Int_t type = it->GetType();  // 0-hydro, 1-jets

    if (im1> -1) {
      TParticle *mother = (TParticle*) (fParticles->UncheckedAt(im1+1));	   
      mother->SetLastDaughter(nump);
      if (mother->GetFirstDaughter()==-1)
 	mother->SetFirstDaughter(nump+1);
    }

    nump++;

    TParticle* p = new TParticle(it->Encoding(), type,	                        //pdg,stat
				 im1, im2, id1, id2,				                //m1,m2,d1,d2
				 mom[0], mom[1], mom[2], TMath::Sqrt(mom.Mag2() + m1 * m1),	//px,py,pz,e
				 pos[0]*1.e-13, pos[1]*1.e-13, pos[2]*1.e-13, 
				 it->T()*1.e-13/3e10);    		        //x,y,z,t
  
     p->SetUniqueID(nump);
     fParticles->Add(p);
  } //end for
 
  fAllocator.FreeList(fSecondariesList);

  return fParticles;
}

Int_t TUHKMgen::ImportParticles(TClonesArray *particles, const Option_t* option)
{
  // Function overloading the TGenerator::ImportParticles() member function.
  // The particles from the local particle list (fSecondariesList) are
  // forwarded to the TGenerator::fParticles
  option = option;   // just to avoid the warning
  

  if(particles==0) return 0;
  TClonesArray &particlesR=*particles;
  particlesR.Clear();
 
  Int_t numprim,numsec;  numprim=numsec=0;
  Int_t nump = 0;
  LPIT_t it,e;
  
  for(it = fSecondariesList.begin(), e = fSecondariesList.end(); it != e; it++) {
    TVector3 pos(it->Pos().Vect());  
    TVector3 mom(it->Mom().Vect());  
    Float_t m1 = it->TableMass();  
    Int_t im1 = it->GetMother();  
    Int_t im2 = -1;
    Int_t id1 = -1;
    Int_t id2 = -1;
    
    Int_t type = it->GetType();  
      
    if (im1> -1) {
      // particle not a primary -> set the daughter indexes for the mother particle"<< endl;
      TParticle *mother = (TParticle*) (particlesR.UncheckedAt(im1));
      mother->SetLastDaughter(nump);
      if(mother->GetFirstDaughter()==-1)
	mother->SetFirstDaughter(nump);
    } else 
      ++numprim;
    
    new (particlesR[nump]) TParticle(it->Encoding(), type,			                        //pdg,stat
				     im1, im2, id1, id2,				                //m1,m2,d1,d2
				     mom[0], mom[1], mom[2], TMath::Sqrt(mom.Mag2() + m1 * m1),	//px,py,pz,e
				     pos[0]*1.e-13, pos[1]*1.e-13, pos[2]*1.e-13, 
				     it->T()*1.e-13/3e10);    		        //x,y,z,t
    
    particlesR[nump]->SetUniqueID(nump);
    nump++;
    numsec++;
  }//end for
  
  fSecondariesList.clear();
  printf("Scan and add prim %d sec %d and all %d particles\n",
  	 numprim,numsec,nump);
  return nump;
}

//______________________________________________________________________________
void TUHKMgen::Initialize()
{
  // Function overloading the TGenerator::Initialize() member function.
  // The Monte-Carlo model is initialized (input parameters are transmitted, 
  // particle list and decay channels are loaded, average multiplicities are calculated, etc.)
 
  fInitialState = new InitialStateHydjet();  
  SetAllParameters();  
  fInitialState->LoadPDGInfo();  
  // set the stable flags
  for(Int_t i=0; i<fStableFlagged; i++) 
    fInitialState->SetPDGParticleStable(fStableFlagPDG[i], fStableFlagStatus[i]); 

  if(!fInitialState->MultIni()) 
    Error("TUHKMgen::Initialize()", "Bad status return from MultIni(). Check it out!! \n");  //
}

void TUHKMgen::Print(const Option_t*) const
{
  cout << "**********************************************************************************" << endl;
  cout << "* UHKM Generator interface to ROOT::TGenerator                                   *" << endl;
  cout << "*  Documentation:                                                                *" << endl;
  cout << "* I.P.Lokhtin, L.V.Malinina, S.V.Petrushanko, A.M.Snigirev, I.Arsene, K.Tywoniuk *" << endl;
  cout << "*    Comput.Phys.Commun.180:779-799, 2009                                        *" << endl;
  cout << "**********************************************************************************" << endl;
}

void TUHKMgen::GenerateEvent()
{
  // Member function overloading the TGenerator::GenerateEvent()
  // The HYDJET++ model is run and the particle lists (fSourceList and fSecondariesList) are filled
  
  fInitialState->Initialize(fSecondariesList, fAllocator);  
 
  if(fSecondariesList.empty())Error("TUHKMgen::GenerateEvent()", "Source particle list empty after fireball initialization!! \n");  //

  // Run the decays

  if(fInitialState->RunDecays())  
    fInitialState->Evolve(fSecondariesList, fAllocator, fInitialState->GetWeakDecayLimit());  
  
}

void TUHKMgen::SetAllParameters() {
  // forward all model input parameters to the InitialStateHydjet object
  // which will handle all the Monte-Carlo simulation using HYDJET++ model

  fInitialState->fParams.fSqrtS = fHydjetParams.fSqrtS;
  fInitialState->fParams.fAw = fHydjetParams.fAw;
  fInitialState->fParams.fIfb = fHydjetParams.fIfb ;
  fInitialState->fParams.fBfix = fHydjetParams.fBfix;
  fInitialState->fParams.fBmin = fHydjetParams.fBmin;
  fInitialState->fParams.fBmax = fHydjetParams.fBmax;
  fInitialState->fParams.fSeed = fHydjetParams.fSeed;
  fInitialState->fParams.fT = fHydjetParams.fT;
  fInitialState->fParams.fMuB = fHydjetParams.fMuB;
  fInitialState->fParams.fMuS = fHydjetParams.fMuS;
  fInitialState->fParams.fMuI3 = fHydjetParams.fMuI3;
  fInitialState->fParams.fThFO = fHydjetParams.fThFO;
  fInitialState->fParams.fMu_th_pip = fHydjetParams.fMu_th_pip;

  fInitialState->fParams.fTau = fHydjetParams.fTau;
  fInitialState->fParams.fSigmaTau = fHydjetParams.fSigmaTau;
  fInitialState->fParams.fR = fHydjetParams.fR;
  fInitialState->fParams.fYlmax = fHydjetParams.fYlmax;
  fInitialState->fParams.fUmax = fHydjetParams.fUmax;
  fInitialState->fParams.fDelta = fHydjetParams.fDelta;
  fInitialState->fParams.fEpsilon = fHydjetParams.fEpsilon;

  fInitialState->fParams.fWeakDecay = fHydjetParams.fWeakDecay;
  fInitialState->fParams.fDecay = 1;                   // run decays
  fInitialState->fParams.fCorrS = fHydjetParams.fCorrS;
  fInitialState->fParams.fEtaType = fHydjetParams.fEtaType;
 
  fInitialState->fParams.fPtmin = fHydjetParams.fPtmin;
  fInitialState->fParams.fNhsel = fHydjetParams.fNhsel;
  fInitialState->fParams.fIshad = fHydjetParams.fIshad;
  fInitialState->fParams.fT0 = fHydjetParams.fT0;
  fInitialState->fParams.fTau0 = fHydjetParams.fTau0;
  fInitialState->fParams.fNf = fHydjetParams.fNf;      
  fInitialState->fParams.fIenglu = fHydjetParams.fIenglu;     
  fInitialState->fParams.fIanglu = fHydjetParams.fIanglu;  

  fInitialState->SetPDGParticleFilename(fParticleFilename);
  fInitialState->SetPDGDecayFilename(fDecayFilename);
  //  fInitialState->SetUseCharmParticles(fUseCharmParticles);
  //  fInitialState->SetWidthRange(fMinWidth, fMaxWidth);
  //  fInitialState->SetMassRange(fMinMass, fMaxMass);
  //  for(Int_t i=0; i<fStableFlagged; i++) 
  //    fInitialState->SetPDGParticleStable(fStableFlagPDG[i], fStableFlagStatus[i]);
  //  cout << "TUHKMgen::SetAllParameters() OUT" << endl;
}
