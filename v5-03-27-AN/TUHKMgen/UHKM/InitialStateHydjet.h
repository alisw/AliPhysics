// InitialStateHydjet is the class which controls the entire 
// event simulation (input parameters, physics, output)
// InitialStateHydjet inherits from class InitialState
// MultIni() member function initializes PYQUEN and calculates the average soft multiplicities 
// Initialize() member function calls the PYQUEN soubroutine for the hard part of the event
// The resonance decay is performed by the function InitialState::Evolve()

/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

#ifndef INITIALSTATEHYDJET_H
#define INITIALSTATEHYDJET_H

#include "InitialState.h"

class ParticleAllocator;

const Int_t kNPartTypes = 1000;

struct InitialParamsHydjet_t {

  Int_t fNevnt; //number of events
  Double_t fSqrtS;        //cms energy per nucleon
  Double_t fAw;        // atomic number of colliding nuclei
  Int_t fIfb;      // flag of type of centrality generation (=0 is fixed by fBfix, not 0 
                                   //impact parameter is generated in each event between fBfmin 
                                   //and fBmax according with Glauber model (f-la 30)
   Double_t fBmin;         //minimum impact parameter in units of nuclear radius RA 
   Double_t fBmax;         //maximum impact parameter in units of nuclear radius RA
   Double_t fBfix;         //fix impact parameter in units of nuclear radius RA

   Int_t fSeed;         //parameter to set the random nuber seed (=0 the current time is used
                                   //to set the random generator seed, !=0 the value fSeed is 
                                   //used to set the random generator seed and then the state of random
                                   //number generator in PYTHIA MRPY(1)=fSeed
       
   Double_t fT;     //chemical freeze-out temperature in GeV    
   Double_t fMuB;     //baryon potential 
   Double_t fMuS;    //strangeness potential 
   Double_t fMuI3;    //isospin potential   
   Double_t fThFO;    //thermal freeze-out temperature T^th in GeV
   Double_t fMu_th_pip;    // effective chemical potential of positivly charged pions at thermal in GeV 

       
   Double_t fTau;     //proper time value
   Double_t fSigmaTau;     //its standart deviation (emission duration)
   Double_t fR;     //maximal transverse radius 
   Double_t fYlmax;     //maximal longitudinal rapidity 
   Double_t fUmax;     //maximal transverse velocity multiplaed on \gamma_r 
   Double_t fDelta;     //momentum asymmetry parameter
   Double_t fEpsilon;     //coordinate asymmetry parameter
  
   Int_t fDecay;    // flag to switch on/off hadron decays<0: decays off,>=0: decays on, (default: 0)
   Double_t fWeakDecay;    //flag to switch on/off weak hadron decays <0: decays off, >0: decays on, (default: 0)
   Int_t fPythDecay;    //Flag to choose how to decay resonances in high-pt part, fPythDecay: 0 by PYTHIA decayer, 
                                     //1 by FASTMC decayer(mstj(21)=0)  
  
   Int_t fEtaType;     // flag to choose rapidity distribution, if fEtaType<=0, 
                                     //then uniform rapidity distribution in [-fYlmax,fYlmax] if fEtaType>0,
                                     //then Gaussian with dispertion = fYlmax 
  
   Int_t fTMuType;     // flag to use calculated chemical freeze-out temperature,
                                     //baryon potential and strangeness potential as a function of fSqrtS 

   Double_t fCorrS;     // flag and value to include strangeness supression factor    
   Int_t fNhsel;         //flag to switch on/off jet and hydro-state production (0: jet
                                     // production off and hydro on, 1: jet production on and jet quenching
                                     // off and hydro on, 2: jet production on and jet quenching on and
                                     // hydro on, 3: jet production on and jet quenching off and hydro
                                     // off, 4: jet production on and jet quenching on and hydro off
 
   Int_t fIshad;         //flag to switch on/off impact parameter dependent nuclear
                                    // shadowing for gluons and light sea quarks (u,d,s) (0: shadowing off,
                                    // 1: shadowing on for fAw=207, 197, 110, 40, default: 1
  
   Double_t fPtmin;       //minimal transverse momentum transfer p_T of hard
                                   // parton-parton scatterings in GeV (the PYTHIA parameter ckin(3)=fPtmin)
   
//  PYQUEN energy loss model parameters:
 
   Double_t fT0;          // initial temperature (in GeV) of QGP for
                                   //central Pb+Pb collisions at mid-rapidity (initial temperature for other
                                  //centralities and atomic numbers will be calculated automatically) (allowed range is 0.2<fT0<2) 
  
   Double_t fTau0;        //proper QGP formation time in fm/c (0.01<fTau0<10)
   Int_t fNf;          //number of active quark flavours N_f in QGP fNf=0, 1,2 or 3 
   Int_t fIenglu;      // flag to fix type of in-medium partonic energy loss 
                                  //(0: radiative and collisional loss, 1: radiative loss only, 2:
                                  //collisional loss only) (default: 0);
   Int_t fIanglu;      //flag to fix type of angular distribution of in-medium emitted
                                  // gluons (0: small-angular, 1: wide-angular, 2:collinear) (default: 0).
  
  
  Int_t    fNPartTypes;                   //counter of hadron species  
  Int_t    fPartEnc[kNPartTypes];                 //Hadron encodings. Maximal number of hadron species is 100!!!
  Double_t fPartMult[2*kNPartTypes];                //Multiplicities of hadron species
  Double_t fPartMu[2*kNPartTypes];                //Chemical potentials of hadron species

  Double_t fMuTh[kNPartTypes];                    //Chemical potentials at thermal freezeout of hadron species

};


class InitialStateHydjet : public InitialState {
 public:
   
  InitialStateHydjet() : fParams(), fVolEff(0), fBgen(0), fNpart(0), fNcoll(0) {};
  ~InitialStateHydjet() {};
  
  void SetVolEff(Double_t value) {fVolEff = value;}
  Double_t GetVolEff() const {return fVolEff;}
  //  virtual Double_t GetTime() {return fParams.fDecay;}
  virtual Bool_t RunDecays() {return (fParams.fDecay>0 ? kTRUE : kFALSE);}
  virtual Int_t GetNev() {return fParams.fNevnt;}
  virtual Double_t GetWeakDecayLimit() {return fParams.fWeakDecay;}  
  
  virtual void Initialize(List_t &source, ParticleAllocator &allocator);
  virtual Bool_t ReadParams();
  virtual Bool_t MultIni();
  virtual void GetCentrality(Double_t& b, Double_t & npart, Double_t & nbin)
      {b = fBgen; npart = fNpart; nbin = fNcoll;}
  
  Bool_t IniOfThFreezeoutParameters();
  
  InitialParamsHydjet_t fParams;             // the list of initial state parameters
  
  private:
   Double_t fVolEff;                     // the effective volume
   // Collision geometry
   Double_t    fBgen;                       // Generated impact parameter
   Double_t    fNpart;                      // Number of participants
   Double_t    fNcoll;                      // Number of collisions 
      
   Double_t F2(Double_t x, Double_t y);
   
   Double_t SimpsonIntegrator(Double_t a, Double_t b, Double_t phi);
   Double_t SimpsonIntegrator2(Double_t a, Double_t b);
   Double_t MidpointIntegrator2(Double_t a, Double_t b);
};

#endif
