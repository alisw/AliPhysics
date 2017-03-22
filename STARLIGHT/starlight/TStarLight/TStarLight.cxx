// -*- C++ -*-
////////////////////////////////////////////////////////////////////////
//
// Copyright 2013
//
////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: -1                      $: revision of last commit
// $Authro::                      $: Author of last commit
// $Date::                        $: Date of last commit
//
// Description:
//     TStarLight.h is the include file defining data members and
// functions specifically needed to implement an interface of STARlight
// to ROOT's TGenerator (a part of ROOT's Virtual Monte Carlo.
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>

#include <Riostream.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TParticle.h>

#include "inputParameters.h"
#include "starlight.h"
#include "TStarLight.h"
#include "vector3.h"

ClassImp(TStarLight);

//----------------------------------------------------------------------
TStarLight::TStarLight()
  : TGenerator()           // Default initlization of base class
  , fErrorStatus(0)        // Error status flag 0=OK
  , fConfigFileName("")    // Confiuration file name
  , fStarLight(NULL)       // STARlight simulation class.
  , fInputParameters()     // Input to simulation class
  , fEvent() {}            // object holdng STARlight simulated event

//----------------------------------------------------------------------
TStarLight::TStarLight(const char* name,         // The name of this object in the root name tables
		       const char* title,        // A title for this object
		       const char* slConfigFile) // file used to configure STARlight.
  : TGenerator(name, title)       // Default initlization of base class
  , fErrorStatus(0)               // Error status flag 0=OK
  , fConfigFileName(slConfigFile) // Confiuration file name
  , fStarLight(NULL)              // STARlight simulation class.
  , fInputParameters()            // Input to simulation class
  , fEvent()                      // object holding STARlight simulated event.
{
  // if (NULL == fInputParameters) {
  //   fErrorStatus = -5; // Init failed. Creating inputParamtere class    
  //   Error("TStarLight", "creating inputParameters class failed");
  //   return;
  // } // end if
  if (fConfigFileName != "" && gSystem->ExpandPathName(fConfigFileName)) { // if error
    fErrorStatus = -4; // Init failed. error in path name
    Error("TStarLight", "Error expanding path name='%s'", fConfigFileName.Data());
    return;
  } // end if
  fStarLight = new starlight;
  if (!fStarLight) {
    fErrorStatus = -2; // Init failed. no simulation
    Error("TStarLight", "Failed to create simulator STARlight");
    return;
  } // end if
  if (fConfigFileName != "")
    fInputParameters.configureFromFile(slConfigFile);
}
//----------------------------------------------------------------------
TStarLight::~TStarLight()
{
  if (fStarLight)
    delete fStarLight;
  fStarLight = NULL;
}

//----------------------------------------------------------------------
void TStarLight::GenerateEvent() {

  if (NULL == fStarLight) {
    fErrorStatus = -1; // generate failed. No generator.
    Fatal("GenerateEvent", "TStarLight class/object not properly constructed");
    return;
  }
  fEvent = fStarLight->produceEvent();
}

// in some cases starlight does not fill the time(energy) component
TParticle* SetEnergyIfNeeded(TParticle *p) {
  if (p->Energy())
    return p;

  p->SetMomentum(p->Px(),
		 p->Py(),
		 p->Pz(),
		 TMath::Sqrt(p->P()*p->P() + p->GetMass()*p->GetMass()));
  return p;
}

//----------------------------------------------------------------------
Int_t TStarLight::ImportParticles(TClonesArray *part, // Pointer to array of particles to be filled.
				  Option_t *opt) {    // A character array of options.
  // Return:
  //   The number of particles added to the TClonesArray *part.

  if (NULL == part)
    return 0;

  TClonesArray &clonesParticles = *part;
  clonesParticles.Clear();

  Int_t nVtx(0);
  Double_t vtx(0), vty(0), vtz(0), vtt(0);
  const std::vector<vector3>* slVtx = fEvent.getVertices();
  if (NULL == slVtx) { // not vertex assume 0,0,0,0;
    vtx = vty = vtz = vtt = 0.0;
  } else { // a vertex exits
    slVtx = fEvent.getVertices();
    nVtx  = slVtx->size();
  }
  const std::vector<starlightParticle>* slPartArr = fEvent.getParticles();
  const Int_t npart(slPartArr->size());
  if (!strcmp(opt,"") || !strcmp(opt,"Final")) {
    for (Int_t ipart=0; ipart<npart; ++ipart) {
      const starlightParticle* slPart = &(slPartArr->at(ipart));
      if (Stable(slPart->getPdgCode())) {
	//slPart->setParent(-1);
	if (nVtx<1) { // No verticies
	  vtx = vty = vtz = vtt = 0.0;
	} else {
	  vtx = (slVtx->at(ipart <nVtx ? ipart : 0)).X();
	  vty = (slVtx->at(ipart <nVtx ? ipart : 0)).Y();
	  vtz = (slVtx->at(ipart <nVtx ? ipart : 0)).Z();
	  vtt = 0.0; // no time given.
	} // end if
 	new(clonesParticles[ipart]) TParticle(slPart->getPdgCode(),
 					      slPart->getStatus(),
 					      -1,/*slPart->getParent(),*/
 					      -1,
 					      slPart->getFirstDaughter(),
 					      slPart->getLastDaughter(),
 					      slPart->GetPx(),
 					      slPart->GetPy(),
 					      slPart->GetPz(),
 					      slPart->GetE(),
 					      vtx,vty,vtz,vtt);
      }
    }
  } else if (!strcmp(opt,"ALL")) { // store all particles
    for (Int_t ipart=0; ipart<npart; ++ipart) {
      const starlightParticle* slPart = &(slPartArr->at(ipart));
      //slPart->setParent(-1);
      if (nVtx < 1) { // No verticies
	vtx = vty = vtz = vtt = 0.0;
      } else {
	vtx = (slVtx->at((ipart<nVtx?ipart:0))).X();
	vty = (slVtx->at((ipart<nVtx?ipart:0))).Y();
	vtz = (slVtx->at((ipart<nVtx?ipart:0))).Z();
	vtt = 0.0; // no time given.
      }
      new(clonesParticles[ipart]) TParticle(slPart->getPdgCode(),
					    slPart->getStatus(),
					    -1,/*slPart->getParent(),*/
					    -1,
					    slPart->getFirstDaughter(),
					    slPart->getLastDaughter(),
					    slPart->GetPx(),
					    slPart->GetPy(),
					    slPart->GetPz(),
					    slPart->GetE(),
					    vtx,vty,vtz,vtt);
    } // end for ipart
  } else {
    fErrorStatus = -1; // Import particles failed unknown option
    Error("ImportParticles", "Unknown option '%s'", opt);
  } // end if opt

  for (Int_t ipart=0; ipart<npart; ++ipart)
    SetEnergyIfNeeded(dynamic_cast<TParticle*>(clonesParticles[ipart]));

  return npart;
}

//----------------------------------------------------------------------
TObjArray* TStarLight::ImportParticles(Option_t *opt) { // A character array of options
  Int_t nVtx(0);
  Double_t vtx(0), vty(0), vtz(0), vtt(0);
  const std::vector<vector3>* slVtx(fEvent.getVertices());
  if (slVtx == 0) { // not vertex assume 0,0,0,0;
    vtx = vty = vtz = vtt = 0.0;
  } else { // a vertex exits
    slVtx = fEvent.getVertices();
    nVtx = slVtx->size();
  } // end if
  const std::vector<starlightParticle>* slPartArr(fEvent.getParticles());
  const Int_t npart(fEvent.getParticles()->size());
  if (!strcmp(opt, "") || !strcmp(opt, "Final")) {
    for (Int_t ipart=0; ipart<npart; ipart++) {
      const starlightParticle* slPart(&(slPartArr->at(ipart)));
      if (Stable(slPart->getPdgCode())) {
	//slPart->setParent(-1);
	if (nVtx < 1) { // No verticies
	  vtx = vty = vtz = vtt = 0.0;
	} else {
	  vtx = (slVtx->at((ipart < nVtx ? ipart : 0))).X();
	  vty = (slVtx->at((ipart < nVtx ? ipart : 0))).Y();
	  vtz = (slVtx->at((ipart < nVtx ? ipart : 0))).Z();
	  vtt = 0.0; // no time given.
	} // end if
	TParticle *p = new TParticle(slPart->getPdgCode(),
				     slPart->getStatus(),
				     -1,/*slPart->getParent(),*/
				     -1,
				     slPart->getFirstDaughter(),
				     slPart->getLastDaughter(),
				     slPart->GetPx(),
				     slPart->GetPy(),
				     slPart->GetPz(),
				     slPart->GetE(),
				     vtx,vty,vtz,vtt);
	fParticles->Add(SetEnergyIfNeeded(p));
      }
    } // end for ipart
  } else if (!strcmp(opt,"ALL")) {// store all particles
    for(Int_t ipart=0;ipart<npart;ipart++) {
      const starlightParticle* slPart(&(slPartArr->at(ipart)));
      //slPart->setParent(-1);
      if (nVtx < 1) { // No verticies
	vtx = vty = vtz = vtt = 0.0;
      } else {
	vtx = (slVtx->at((ipart < nVtx ? ipart : 0))).X();
	vty = (slVtx->at((ipart < nVtx ? ipart : 0))).Y();
	vtz = (slVtx->at((ipart < nVtx ? ipart : 0))).Z();
	vtt = 0.0; // no time given.
      } // end if
      TParticle* p = new TParticle(slPart->getPdgCode(),
				   slPart->getStatus(),
				   -1,/*slPart->getParent(),*/
				   -1,
				   slPart->getFirstDaughter(),
				   slPart->getLastDaughter(),
				   slPart->GetPx(),
				   slPart->GetPy(),
				   slPart->GetPz(),
				   slPart->GetE(),
				   vtx,vty,vtz,vtt);
      fParticles->Add(SetEnergyIfNeeded(p));
    } // end for ipart
  }else{
    fErrorStatus = -1; // Import particles failed unknown option
    Error("ImportParticles", "Unknown option '%s'", opt);
  } // end if opt
  return fParticles;
}
//----------------------------------------------------------------------
void TStarLight::SetParameter(const char* line) {
  const std::string sl(line);
  if (not fInputParameters.setParameter(sl))
    Fatal("SetParameter", "cannot set parameter: '%s'", line);
}
void TStarLight::SetParameter(const char* name,
			      Double_t val) {
  const TString line(TString::Format("%s = %e", name, val));
  const std::string sl(line.Data());
  if (not fInputParameters.setParameter(sl))
    Fatal("SetParameter", "cannot set parameter: '%s'", name);
}

//----------------------------------------------------------------------
Double_t TStarLight::GetParameter(const char* name) const {
  const std::string sn(name);
  std::stringstream ioss;
  fInputParameters.write(ioss);
  std::string line;
  while (getline(ioss, line)) {
    const size_t index_colon(line.find(':'));
    if (index_colon == std::string::npos)
      continue;
    std::string key(line.substr(0,index_colon));
    if (key != sn) continue;
    std::string val(line.substr(index_colon+1, std::string::npos));
    std::istringstream iss(val);
    Double_t result(0);
    iss >> result;
    return result;
  }
  Fatal("GetParameter", "parameter '%s' not found", name);
  return 0.0;
}
