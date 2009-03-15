// This class implements the FASTMC generator model into AliROOT.
// For the aliroot implementation we make use only of the Bjorken-like
// initial state.
// Please look into the FASTMC articles for more documentation on the 
// needed parameters.

#ifndef TUHKMgen_H
#define TUHKMgen_H

#ifndef ROOT_TGenerator
#include "TGenerator.h"
#endif

#ifndef INITIALSTATEHYDJET_INCLUDED
#include "InitialStateHydjet.h"
#endif

#ifndef DATABASE_PDG
#include "DatabasePDG.h"
#endif

#include <string>
using namespace std;

class TUHKMgen : public TGenerator {
 protected:
  InitialStateHydjet *fInitialState;
  ParticleAllocator fAllocator;
  List_t fSourceList;
  List_t fSecondariesList;
  Int_t  fNPprim;
  Int_t  fNPsec;
  InitialParamsHydjet_t fHydjetParams;       // list of parameters for the initial state
  // details for the PDG database
  Char_t fParticleFilename[256];               // particle list filename
  Char_t fDecayFilename[256];                  // decay table filename
  Int_t fStableFlagPDG[500];                   // array of PDG codes flagged to be stable
  Bool_t fStableFlagStatus[500];               // array of decay flag status
  Int_t fStableFlagged;                        // number of toggled decay flags
  Bool_t fUseCharmParticles;                   // flag to turn on/off the use of charm particles
  Double_t fMinWidth;                          // minimum decay width for the particles to be used from the PDG database
  Double_t fMaxWidth;                          // maximum ----
  Double_t fMinMass;                           // minimum mass for the particles to be used from the PDG database
  Double_t fMaxMass;                           // maximum ----

  void SetAllParameters();

 public:   
  TUHKMgen();
  virtual      ~TUHKMgen();
  virtual void  Initialize();
  virtual void  GenerateEvent();
  virtual Int_t ImportParticles(TClonesArray *particles,Option_t* option="prim");
  virtual TObjArray* ImportParticles(Option_t* option="prim");
  // this function makes available the PDG info in our database 
  virtual DatabasePDG* PDGInfo() const {return fInitialState->PDGInfo();}

  // Setters
  // set reasonable default parameters suited for central Au+Au collisions at RHIC(200GeV)
  void SetAllParametersRHIC();
  // set reasonable default parameters suited for central Pb+Pb collisions at LHC(5.5TeV)
  void SetAllParametersLHC();

  void SetEcms(Double_t value)             {fHydjetParams.fSqrtS = value;}          // CMS energy per nucleon [GeV] (<2.24 given temperature and ch pot are used)
  void SetAw(Double_t value)               {fHydjetParams.fAw = value;}             // nuclei mass number
  void SetIfb(Int_t value)                 {fHydjetParams.fIfb = value;}            //b-simulation: 0-fix,1-distributed
  void SetBfix(Double_t value)             {fHydjetParams.fBfix = value;}           // fix impact parameter
  void SetBmin(Double_t value)             {fHydjetParams.fBmin = value;}           // Minimum impact parameter
  void SetBmax(Double_t value)             {fHydjetParams.fBmax = value;}           // Maximum impact parameter
  void SetChFrzTemperature(Double_t value) {fHydjetParams.fT = value;}              // Temperature for the chemical freezeout [GeV]
  void SetMuB(Double_t value)              {fHydjetParams.fMuB = value;}            // Baryonic chemical potential [GeV]
  void SetMuS(Double_t value)              {fHydjetParams.fMuS = value;}            // Strangeness chemical potential [GeV]
  void SetMuQ(Double_t value)              {fHydjetParams.fMuI3 = value;}  // Isospin chemical potential [GeV]
  void SetThFrzTemperature(Double_t value) {fHydjetParams.fThFO = value;}           // Temperature for the thermal freezeout [GeV]
  void SetMuPionThermal(Double_t value)    {fHydjetParams.fMu_th_pip = value;}      // Chemical potential for pi+ at thermal freezeout [GeV]
  
  
  void SetSeed(Int_t value) {fHydjetParams.fSeed = value;}  //parameter to set the random nuber seed (=0 the current time is used
  //to set the random generator seed, !=0 the value fSeed is 
  //used to set the random generator seed and then the state of random
  //number generator in PYTHIA MRPY(1)=fSeed
  
  
  
  void SetTauB(Double_t value)             {fHydjetParams.fTau = value;}            // Proper time for the freeze-out hypersurface [fm/c]
  void SetSigmaTau(Double_t value)         {fHydjetParams.fSigmaTau = value;}       // Standard deviation for the proper time (emission duration) [fm/c]
  void SetRmaxB(Double_t value)            {fHydjetParams.fR = value;}              // Maximal transverse radius [fm]
  void SetYlMax(Double_t value)            {fHydjetParams.fYlmax = value;}          // Maximal fireball longitudinal rapidity
  void SetEtaRMax(Double_t value)          {fHydjetParams.fUmax = value;}           // Maximal transverse velocity
  void SetMomAsymmPar(Double_t value)      {fHydjetParams.fDelta = value;}          // Momentum asymmetry parameter
  void SetCoordAsymmPar(Double_t value)    {fHydjetParams.fEpsilon = value;}        // Coordinate asymmetry parameter


  void SetFlagWeakDecay(Int_t value)    {fHydjetParams.fWeakDecay = value;} //flag to switch on/off weak hadron decays <0: decays off, >0: decays on, (default: 0)


  void SetEtaType(Int_t value)    {fHydjetParams.fEtaType = value;}     // flag to choose rapidity distribution, if fEtaType<=0, 
                                     //then uniform rapidity distribution in [-fYlmax,fYlmax] if fEtaType>0,
                                     //then Gaussian with dispertion = fYlmax 
  

  void SetGammaS(Double_t value)           {fHydjetParams.fCorrS = value;}          // Strangeness suppression parameter (if gamma_s<=0 then it will be calculated)
 //not needed now  void SetHdec(Double_t value)             {fHydjetParams.fTime = value;}           // Enable/disable hadronic decays (<0 no decays, >=0 decays) 

//PYQUEN parameters 
  void SetPyquenNhsel(Int_t value)         {fHydjetParams.fNhsel = value;}          // Flag to choose the type of event to be generated
 
  void SetPyquenShad(Int_t value)         {fHydjetParams.fIshad = value;}         //flag to switch on/off impact parameter dependent nuclear
                                    // shadowing for gluons and light sea quarks (u,d,s) (0: shadowing off,
                                    // 1: shadowing on for fAw=207, 197, 110, 40, default: 1
 
  void SetPyquenPtmin(Double_t value)      {fHydjetParams.fPtmin = value;}          // Pyquen input parameter for minimum Pt of parton-parton scattering (5GeV<pt<500GeV)
                                                                                      // fNhsel = 0 --> UHKM fireball, no jets
                                                                                      // fNhsel = 1 --> UHKM fireball, jets with no quenching
                                                                                      // fNhsel = 2 --> UHKM fireball, jets with quenching 
                                                                                      // fNhsel = 3 --> no UHKM fireball, jets with no quenching
                                                                                      // fNhsel = 4 --> no UHKM fireball, jets with quenching

  void SetPyquenT0(Double_t value)      {fHydjetParams.fT0 = value;}        //proper QGP formation tempereture
  void SetPyquenTau0(Double_t value)      {fHydjetParams.fTau0 = value;}        //proper QGP formation time in fm/c (0.01<fTau0<10)
  void SetPyquenNf(Int_t value)      {fHydjetParams.fNf = value;}           //number of active quark flavours N_f in QGP fNf=0, 1,2 or 3  
  void SetPyquenIenglu(Int_t value)      {fHydjetParams.fIenglu = value;}      // flag to fix type of in-medium partonic energy loss 
                                  //(0: radiative and collisional loss, 1: radiative loss only, 2:
                                  //collisional loss only) (default: 0);  
  void SetPyquenIanglu(Int_t value)      {fHydjetParams.fIanglu = value;}  //flag to fix type of angular distribution of in-medium emitted
                                  // gluons (0: small-angular, 1: wide-angular, 2:collinear) (default: 0).
  
 
  void SetPDGParticleFile(Char_t *name)              {strcpy(fParticleFilename, name);}     // Set the filename containing the particle PDG info
  void SetPDGDecayFile(Char_t *name)                 {strcpy(fDecayFilename, name);}        // Set the filename containing the PDG decay channels info
  void SetPDGParticleStable(Int_t pdg, Bool_t value) {                                      // Turn on/off the decay flag for a PDG particle
    fStableFlagPDG[fStableFlagged] = pdg;
    fStableFlagStatus[fStableFlagged++] = value;
  } 
  void SetUseCharmParticles(Bool_t flag) {fUseCharmParticles = flag;}
  void SetMinimumWidth(Double_t value) {fMinWidth = value;}
  void SetMaximumWidth(Double_t value) {fMaxWidth = value;}
  void SetMinimumMass(Double_t value) {fMinMass = value;}
  void SetMaximumMass(Double_t value) {fMaxMass = value;}

  //Getters
   
  Double_t GetEcms()             {return fHydjetParams.fSqrtS;}          
  Double_t GetAw()               {return fHydjetParams.fAw;}        
  Double_t GetIfb()              {return fHydjetParams.fIfb;}
  Double_t GetBfix()             {return fHydjetParams.fBfix;}       
  Double_t GetBmin()             {return fHydjetParams.fBmin;}       
  Double_t GetBmax()             {return fHydjetParams.fBmax;}       
  Double_t GetChFrzTemperature() {return fHydjetParams.fT;} 
  Double_t GetMuB()              {return fHydjetParams.fMuB;}
  Double_t GetMuS()              {return fHydjetParams.fMuS;}
  Double_t GetMuQ()              {return fHydjetParams.fMuI3;}
  Double_t GetThFrzTemperature() {return fHydjetParams.fThFO;}
  Double_t GetMuPionThermal()    {return fHydjetParams.fMu_th_pip;}
  Int_t    GetSeed()             {return fHydjetParams.fSeed;}  
  Double_t GetTauB()             {return fHydjetParams.fTau;}
  Double_t GetSigmaTau()         {return fHydjetParams.fSigmaTau;}
  Double_t GetRmaxB()            {return fHydjetParams.fR;}
  Double_t GetYlMax()            {return fHydjetParams.fYlmax;}
  Double_t GetEtaRMax()          {return fHydjetParams.fUmax;}
  Double_t GetMomAsymmPar()      {return fHydjetParams.fDelta;}
  Double_t GetCoordAsymmPar()    {return fHydjetParams.fEpsilon;}
  Int_t GetFlagWeakDecay()       {return fHydjetParams.fWeakDecay;} 
  Int_t GetEtaType()          {return fHydjetParams.fEtaType;}
  Double_t GetGammaS()           {return fHydjetParams.fCorrS;}  
  Int_t    GetPyquenNhsel()      {return fHydjetParams.fNhsel;}
  Int_t    GetPyquenShad()       {return fHydjetParams.fIshad;}
  Double_t GetPyquenPtmin()      {return fHydjetParams.fPtmin;}
  Double_t GetPyquenT0()         {return fHydjetParams.fT0;}
  Double_t GetPyquenTau0()       {return fHydjetParams.fTau0;}
  Double_t GetPyquenNf()         {return fHydjetParams.fNf;}      
  Double_t GetPyquenIenglu()     {return fHydjetParams.fIenglu;}     
  Double_t GetPyquenIanglu()     {return fHydjetParams.fIanglu;}  
  
  Char_t*  GetPDGParticleFile()  {return fParticleFilename;}
  Char_t*  GetPDGDecayFile()     {return fDecayFilename;}
  Bool_t   GetUseCharmParticles(){return fUseCharmParticles;}
  Double_t GetMinimumWidth()     {return fMinWidth;}
  Double_t GetMaximumWidth()     {return fMaxWidth;}
  Double_t GetMinimumMass()      {return fMinMass;}
  Double_t GetMaximumMass()      {return fMaxMass;}
  
  void Print();

  ClassDef(TUHKMgen, 3)  //Interface to FASTMC Event Generator
};
#endif







