/*
  Copyright   : The FASTMC and SPHMC Collaboration
  Author      : Ionut Cristian Arsene 
  Affiliation : Oslo University, Norway & Institute for Space Sciences, Bucharest, Romania
  e-mail      : i.c.arsene@fys.uio.no
  Date        : 2007/05/30

  This class is using the particle and decay lists provided by the 
  THERMINATOR (Computer Physics Communications 174 669 (2006)) and
  SHARE (Computer Physics Communications 167 229 (2005)) collaborations.
*/

#ifndef PARTICLE_PDG
#include "ParticlePDG.h"
#endif

#include <iostream>
using std::cout;
using std::endl;

ParticlePDG::ParticlePDG() {
  fPDG   = kNonsensePDG;
  fMass  = -1.0;
  fWidth = 0.0;
  fNDecayChannels = 0;
  for(Int_t i=0; i<kMaxDecayChannels; i++)
    fDecayChannels[i] = new DecayChannel();
}

ParticlePDG::ParticlePDG(Char_t *name, Int_t pdg, Double_t mass, Double_t width) {
  for(Int_t i=0; i<9; i++)
    if(*(name+i) != '\0') fName[i] = *(name+i);
    else break;
  fPDG   = pdg;
  fMass  = mass;
  fWidth = width;
  fNDecayChannels = 0;
  for(Int_t i=0; i<kMaxDecayChannels; i++)
    fDecayChannels[i] = new DecayChannel();
}

ParticlePDG::~ParticlePDG() {
  for(Int_t i=0; i<kMaxDecayChannels; i++)
    delete fDecayChannels[i];
}

Double_t ParticlePDG::GetFullBranching() {
  Double_t fullBranching = 0.0;
  for(Int_t i=0; i<fNDecayChannels; i++)
    fullBranching += fDecayChannels[i]->GetBranching();
  return fullBranching;
}

void ParticlePDG::AddChannel(DecayChannel* channel) {
  if(channel->GetMotherPDG() != fPDG) {
    cout << "ERROR in ParticlePDG::AddChannel() : You try to add a channel which has a different mother PDG" << endl;
    return;
  }
  fDecayChannels[fNDecayChannels]->SetMotherPDG(channel->GetMotherPDG());
  fDecayChannels[fNDecayChannels]->SetBranching(channel->GetBranching());
  fDecayChannels[fNDecayChannels]->SetDaughters(channel->GetDaughters(), channel->GetNDaughters());
  fNDecayChannels++;
}
