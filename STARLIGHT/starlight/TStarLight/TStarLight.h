////////////////////////////////////////////////////////////////////////
//
// Copyright 2013
//
////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                         $: revision of last commit
// $Authro::                      $: Author of last commit
// $Date::                        $: Date of last commit
//
// Description:
//     TStarLight.h is the include file defining data members and
// functions specifically needed to implement an interface of STARlight
// to ROOT's TGenerator (a part of ROOT's Virtual Monte Carlo).
//
// Based on work done by Bjoern Nielson
////////////////////////////////////////////////////////////////////////
#ifndef STARLIGHT_TSTARLIGHT_H
#define STARLIGHT_TSTARLIGHT_H
#include <cstring>
#include <iostream>

#include <TGenerator.h>
#include <TString.h>
#include <TMath.h>

#ifdef __CINT__
#  undef __GNUC__
#  define _SYS__SELECT_H_
struct pthread_cond_t;
struct pthread_mutex_t;
#endif

// Star Light includes
#include "starlight.h"
#include "upcevent.h"
#include "inputParameters.h"
#include "randomgenerator.h"

using std::ostream;
//
class TObjArray;
class TClonesArray;
class TMap;

class TStarLight : public TGenerator {
 public:
  TStarLight();

  TStarLight(const char* name,
	     const char* title,
	     const char* slConfigFile="./slight.in");

  virtual ~TStarLight();

  // TGenerator overloaded methods
  virtual void       GenerateEvent();
  virtual Int_t      ImportParticles(TClonesArray *particles, Option_t *opt="");
  virtual TObjArray* ImportParticles(Option_t *opt="");
  virtual void       SetParameter(const char* key, Double_t val);
  virtual Double_t   GetParameter(const char* name) const;

  void ImportEventInfo(TMap *) const;

  // read configuration from a file
  void ImportConfigurationFromFile(const char* filename){
    const std::string sf(filename);
    fInputParameters.configureFromFile(sf);
  }

  void PrintInputs(ostream& out) const {
    fInputParameters.print(out);
  }

  Bool_t InitStarLight() {
    if (not fInputParameters.init()) {
      Fatal("InitStarLight", "parameter initialization has failed");
      fErrorStatus = -1;
      return false;
    }
    fStarLight->setInputParameters(&fInputParameters);
    return fStarLight->init();
  }

  void SetParameter(const char* line);

  void SetInput(const inputParameters &in) {
    fInputParameters = in;
    fStarLight->setInputParameters(&fInputParameters);
    fStarLight->init();
  }
  inputParameters& GetInputParameters() { return fInputParameters; }
  const inputParameters& GetInputParameters() const { return fInputParameters; }

  Int_t GetErrorStatusFlag() const { return fErrorStatus; }

  // boost event to the experiment CM frame
  void BoostEvent() {
    fEvent.boost(0.5*(TMath::ACosH(fInputParameters.beam1LorentzGamma()) -
		      TMath::ACosH(fInputParameters.beam2LorentzGamma())));
  }

 private:
  bool Stable(Int_t pdgCode) const {
    switch(TMath::Abs(pdgCode)) {
    case 11: // electon
    case 12: // e_neutreno
    case 13: // Muon
    case 14: // Muon_neutreno
    case 16: // tau_neutreno
    case 22: // Photon
    case 211:// Charge Pion
    case 130:// K0 long
    case 310:// K0 short
    case 311:// K0
    case 321:// Charged K
    case 2212:// Proton
    case 2112:// neutron
    case 3122:// Lamda
    case 3222:// Sigma+
    case 3112:// Sigma-
    case 3322:// Exi0
    case 3312:// Exi-
    case 3334:// Omega-
      return kTRUE; // Stable enough.
    default:
      return kFALSE; // Not Stable.
    } // end switch
    return kFALSE;
  }
  Int_t            fErrorStatus;      //   Error status flag. 0=OK
  TString          fConfigFileName;   //   Input Configuration file name
  starlight       *fStarLight;        //!  Simulation Class
  inputParameters  fInputParameters;  //   simulation input information.
  upcEvent         fEvent;            //!  object holding STARlight simulated event.

  ClassDef(TStarLight,1); // STARlight interface to ROOT's Virtual Monte Carlo
} ;

#endif
