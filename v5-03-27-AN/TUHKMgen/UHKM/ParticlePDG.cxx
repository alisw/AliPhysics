//
//  Copyright   : The FASTMC and SPHMC Collaboration
//  Author      : Ionut Cristian Arsene 
//  Affiliation : Oslo University, Norway & Institute for Space Sciences, Bucharest, Romania
//  e-mail      : i.c.arsene@fys.uio.no
//  Date        : 2007/05/30
//
//  This class is using the particle and decay lists provided by the 
//  THERMINATOR (Computer Physics Communications 174 669 (2006)) and
//  SHARE (Computer Physics Communications 167 229 (2005)) collaborations.
//

#include <iostream>
#include <string>
using namespace std;
#include "ParticlePDG.h"

//________________________________________________________________________
ParticlePDG::ParticlePDG() :
  fPDG(kNonsensePDG),
  fMass(-1.0),
  fWidth(0.0),
  fSpin(0.0),
  fIsospin(0.0),
  fIsospinZ(0.0),
  fLightQuarkNumber(0.0),
  fAntiLightQuarkNumber(0.0),
  fStrangeQuarkNumber(0.0),
  fAntiStrangeQuarkNumber(0.0),
  fCharmQuarkNumber(0.0),
  fAntiCharmQuarkNumber(0.0),
  fNDecayChannels(0),
  fStable(0.0)
{
//
// default constructor
//
  memset(fName,'a',9);
  for(Int_t i=0; i<kMaxDecayChannels; i++)
    fDecayChannels[i] = new DecayChannel();
}

//________________________________________________________________________
ParticlePDG::ParticlePDG(const Char_t * const name, Int_t pdg, Double_t mass, Double_t width) :
  fPDG(pdg),
  fMass(mass),
  fWidth(width),
  fSpin(0.0),
  fIsospin(0.0),
  fIsospinZ(0.0),
  fLightQuarkNumber(0.0),
  fAntiLightQuarkNumber(0.0),
  fStrangeQuarkNumber(0.0),
  fAntiStrangeQuarkNumber(0.0),
  fCharmQuarkNumber(0.0),
  fAntiCharmQuarkNumber(0.0),
  fNDecayChannels(0),
  fStable(0.0)
{
  //
  // constructor
  //
  for(Int_t i=0; i<9; i++)
    if(*(name+i) != '\0') fName[i] = *(name+i);
    else break;
  for(Int_t i=0; i<kMaxDecayChannels; i++)
    fDecayChannels[i] = new DecayChannel();
}

//________________________________________________________________________
ParticlePDG::~ParticlePDG() {
  //
  // destructor
  //
  for(Int_t i=0; i<kMaxDecayChannels; i++)
    delete fDecayChannels[i];
}

//________________________________________________________________________
Double_t ParticlePDG::GetFullBranching() {
  //
  // calculate the sum of branching ratios from all decay channels (should add up to 1)
  //
  Double_t fullBranching = 0.0;
  for(Int_t i=0; i<fNDecayChannels; i++)
    fullBranching += fDecayChannels[i]->GetBranching();
  return fullBranching;
}

//________________________________________________________________________
void ParticlePDG::AddChannel(DecayChannel* channel) {
  //
  // add a decay channel
  //
  if(channel->GetMotherPDG() != fPDG) {
    cout << "ERROR in ParticlePDG::AddChannel() : You try to add a channel which has a different mother PDG" << endl;
    return;
  }
  fDecayChannels[fNDecayChannels]->SetMotherPDG(channel->GetMotherPDG());
  fDecayChannels[fNDecayChannels]->SetBranching(channel->GetBranching());
  fDecayChannels[fNDecayChannels]->SetDaughters(channel->GetDaughters(), channel->GetNDaughters());
  fNDecayChannels++;
}
