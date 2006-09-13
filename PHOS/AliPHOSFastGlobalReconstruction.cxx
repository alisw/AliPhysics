/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

// --- AliRoot header files ---

// Fast global reconstruction class.
// It performes fast reconstruction for charged particles only,
// assuming that they were detected by all central ALICE detectors but PHOS.
// This class acts as a filter for primary particles, selects them
// and deteriorates their 4-momenta.
// The filter can recognize the primary particle list from different
// event generators: Pythia, Hijing, HIJINGPara.
//!
// Usage:
// rec = new AliPHOSGlobalReconstruction("galice.root");
// rec->FastReconstruction(ievent);
// TClonesArray *recList = rec->GetRecParticles();
//!
// See also $ALICE_ROOT/macros/TestGlobalReconstruction.C
//!
// Author: Yuri Kharlov. 17 April 2003

#include "AliGenerator.h"
#include "AliPHOSGetter.h"
#include "AliPHOSFastGlobalReconstruction.h"

ClassImp(AliPHOSFastGlobalReconstruction)

//____________________________________________________________________________
AliPHOSFastGlobalReconstruction::AliPHOSFastGlobalReconstruction(): 
  fgime(0),
  fGenerator(0),
  fParticles(0),
  fNParticles(0)
{
  //Def ctor.
}

//____________________________________________________________________________
AliPHOSFastGlobalReconstruction::AliPHOSFastGlobalReconstruction(const char* headerFile):
  fgime(AliPHOSGetter::Instance(headerFile)),
  fGenerator(gAlice->Generator()),
  fParticles(new TClonesArray("TParticle",100)),
  fNParticles(0)
{
  // Constructor of fast global reconstruction:
  // create an instance of the PHOS getter,
  // create an array or reconstructed particles.
}

//____________________________________________________________________________
AliPHOSFastGlobalReconstruction::AliPHOSFastGlobalReconstruction(const AliPHOSFastGlobalReconstruction & rhs):
  TObject(rhs),
  fgime(0),
  fGenerator(0),
  fParticles(0),
  fNParticles(0)
{
  //Required by effc++, but not clear for me how to do correct copy. To be fixed.
  Fatal("AliPHOSFastGlobalReconstruction", "copy ctor not implemented");
}

//____________________________________________________________________________
AliPHOSFastGlobalReconstruction & AliPHOSFastGlobalReconstruction::operator = (const AliPHOSFastGlobalReconstruction &)
{
  //Required by effc++, but not clear for me how to do correct copy. To be fixed.
  Fatal("operator = ", "copy ctor not implemented");
  return *this;
}

//____________________________________________________________________________
AliPHOSFastGlobalReconstruction::~AliPHOSFastGlobalReconstruction()
{
  // Destructor of fast global reconstruction:
  // delete the array of reconstructed particles
  if (fParticles != 0) {
    delete fParticles;
    fParticles  = 0;
    fNParticles = 0;
  }
}

//____________________________________________________________________________
void AliPHOSFastGlobalReconstruction::FastReconstruction(Int_t event)
{
  // Perform a fast global reconstruction of event numbered "event".
  // Reconstructed particles will be stored into array recParticles

  TParticle *primary;
  TLorentzVector p,v;
  Int_t kf,ks,imom1,imom2,idaug1,idaug2;

  fgime->Event(event,"X") ;
  fParticles  ->Clear();
  fNParticles = 0;
  Int_t        nPrimaries = fgime->NPrimaries();
  TClonesArray *primaries = fgime->Primaries();

  for (Int_t iprim=0; iprim<nPrimaries; iprim++) {
    primary = (TParticle*)primaries->At(iprim);
    if ((strcmp(fGenerator->GetName(),"Pythia")    ==0 && primary->GetStatusCode() == 1) ||
	(strcmp(fGenerator->GetName(),"HIJINGpara")==0 && primary->GetFirstMother()==-1) ||
	(strcmp(fGenerator->GetName(),"Hijing")    ==0 && primary->GetStatusCode() == 3)) {
      if (Detected(primary)) {
	primary->Momentum(p);
	primary->ProductionVertex(v);
	kf     = primary->GetPdgCode();
	ks     = primary->GetStatusCode();
	imom1  = primary->GetFirstMother();
	imom2  = primary->GetSecondMother();
	idaug1 = primary->GetFirstDaughter();
	idaug2 = primary->GetLastDaughter();
	SmearMomentum(p);
	new((*fParticles)[fNParticles]) TParticle(kf,ks,imom1,imom2,idaug1,idaug2,p,v);
	fNParticles++;
      }
    }
  }
}

//____________________________________________________________________________
Bool_t AliPHOSFastGlobalReconstruction::Detected(TParticle *particle)
{
  // Returns kTRUE is a particle is reconstructed, kFALSE otherwise.
  // A particle is reconstructed if it is charged and accepted with the
  // probability Efficiency(pt,eta) depending on pt and eta.

  Bool_t detected = kFALSE;
  if (particle->GetPDG()->Charge() != 0) {
    Float_t pt  = particle->Pt();
    Float_t eta = particle->Eta();
    if (gRandom->Rndm() < Efficiency(pt,eta)) detected = kTRUE;
  }
  return detected;
}

//____________________________________________________________________________
Float_t AliPHOSFastGlobalReconstruction::Efficiency(Float_t pt, Float_t eta)
{
  // Detection probability vs. pt and eta, i.e. a probability to detect
  // a particle with transverse momentum pt and speudorapidity eta.
  // For the moment assume that charged particles are detected with
  // 80% efficiency within |eta|<0.9 and pt>0.15 GeV, and with 0% efficiency
  // beyond that acceptance

  const Float_t kEtaLimit = 0.9;
  const Float_t kPtLimit  = 0.15;
  Float_t efficiency = 0.0;
  if (TMath::Abs(eta) < kEtaLimit && pt > kPtLimit) efficiency = 0.80;
  return efficiency;
}

//____________________________________________________________________________
void AliPHOSFastGlobalReconstruction::SmearMomentum(TLorentzVector &p)
{
  // Smear 4-momentum according to known resolution (2% for the moment)

  const Float_t kAngleResolution    = 0.02;
  const Float_t kMomentumResolution = 0.02;
  Double_t mass = p.M();
  for (Int_t i=0; i<3; i++) {
    p[i] *= gRandom->Gaus(1.,kAngleResolution   );
    p[i] *= gRandom->Gaus(1.,kMomentumResolution);
  }
  p[3] = TMath::Sqrt(p.P()*p.P() + mass*mass);
}
