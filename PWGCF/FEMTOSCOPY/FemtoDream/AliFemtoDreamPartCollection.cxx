/*
 * AliFemtoPPbpbLamPartCollection.cxx

 *
 *  Created on: Aug 29, 2017
 *      Author: gu74req
 */
#include <iostream>
#include "AliFemtoDreamPartCollection.h"
#include "AliLog.h"
ClassImp(AliFemtoDreamPartCollection)
AliFemtoDreamPartCollection::AliFemtoDreamPartCollection()
    : fHigherMath(),
      fNSpecies(0),
      fZVtxMultBuffer(),
      fValuesZVtxBins(),
      fValuesMultBins() {

}

AliFemtoDreamPartCollection::AliFemtoDreamPartCollection(
    const AliFemtoDreamPartCollection& coll)
    : fHigherMath(coll.fHigherMath),
      fNSpecies(coll.fNSpecies),
      fZVtxMultBuffer(coll.fZVtxMultBuffer),
      fValuesZVtxBins(coll.fValuesZVtxBins),
      fValuesMultBins(coll.fValuesMultBins) {

}
AliFemtoDreamPartCollection::AliFemtoDreamPartCollection(
    AliFemtoDreamCollConfig *conf, bool MinimalBooking)
    : fHigherMath(new AliFemtoDreamHigherPairMath(conf, MinimalBooking)),
      fNSpecies(conf->GetNParticles()),
      fZVtxMultBuffer(
          conf->GetNZVtxBins(),
          std::vector<AliFemtoDreamZVtxMultContainer>(
              conf->GetNMultBins(), AliFemtoDreamZVtxMultContainer(conf))),
      fValuesZVtxBins(conf->GetZVtxBins()),
      fValuesMultBins(conf->GetMultBins()) {
}

AliFemtoDreamPartCollection& AliFemtoDreamPartCollection::operator=(
    const AliFemtoDreamPartCollection& coll) {
  if (this != &coll) {
    this->fHigherMath = coll.fHigherMath;
    this->fNSpecies = coll.fNSpecies;
    this->fZVtxMultBuffer = coll.fZVtxMultBuffer;
    this->fValuesZVtxBins = coll.fValuesZVtxBins;
    this->fValuesMultBins = coll.fValuesMultBins;
  }

  return *this;
}

AliFemtoDreamPartCollection::~AliFemtoDreamPartCollection() {
}

void AliFemtoDreamPartCollection::SetEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles, float ZVtx,
    float Mult, float cent) {
  if (Particles.size() != fNSpecies) {
    TString fatalOut = Form("Too few Species %d for %d", (int) Particles.size(),
                            (int) fNSpecies);
    AliFatal(fatalOut.Data());
  }
  int bins[2] = { 0, 0 };
  FindBin(ZVtx, Mult, bins);
  if (!(bins[0] == -99 || bins[1] == -99)) {
    auto itZVtx = fZVtxMultBuffer.begin();
    itZVtx += bins[0];
    auto itMult = itZVtx->begin();
    itMult += bins[1];
    itMult->PairParticlesSE(Particles, fHigherMath, bins[1], cent);
    itMult->PairParticlesME(Particles, fHigherMath, bins[1], cent);
    itMult->SetEvent(Particles);
  }
  return;
}

void AliFemtoDreamPartCollection::SetEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamEvent* evt) {
  if (Particles.size() != fNSpecies) {
    TString fatalOut = Form("Too few Species %d for %d", (int) Particles.size(),
                            (int) fNSpecies);
    AliFatal(fatalOut.Data());
  }
  int bins[2] = { 0, 0 };
  float ZVtx = evt->GetZVertex();
  float Mult = evt->GetMultiplicity();
  float cent = evt->GetV0MCentrality();
  fHigherMath->SetBField(evt->GetBField());
  FindBin(ZVtx, Mult, bins);
  if (!(bins[0] == -99 || bins[1] == -99)) {
    auto itZVtx = fZVtxMultBuffer.begin();
    itZVtx += bins[0];
    auto itMult = itZVtx->begin();
    itMult += bins[1];
    itMult->PairParticlesSE(Particles, fHigherMath, bins[1], cent);
    itMult->PairParticlesME(Particles, fHigherMath, bins[1], cent);
    itMult->SetEvent(Particles);
  }
  return;
}

void AliFemtoDreamPartCollection::PrintEvent(int ZVtx, int Mult) {
  auto itZVtx = fZVtxMultBuffer.begin();
  itZVtx += ZVtx;
  auto itMult = itZVtx->begin();
  itMult += Mult;
  return;
}

void AliFemtoDreamPartCollection::FindBin(float ZVtxPos, float Multiplicity,
                                          int *returnBins) {
  returnBins[0] = -99;
  returnBins[1] = -99;
  for (auto itBin = fValuesZVtxBins.begin(); itBin != fValuesZVtxBins.end() - 1;
      ++itBin) {
    if (*itBin < ZVtxPos && ZVtxPos <= *(itBin + 1)) {
      returnBins[0] = itBin - fValuesZVtxBins.begin();
      break;
    }
  }
  int binCounter = fValuesMultBins.size();
  for (std::vector<int>::reverse_iterator itBin = fValuesMultBins.rbegin();
      itBin != fValuesMultBins.rend(); ++itBin) {
    binCounter--;
    if (Multiplicity >= *itBin) {
      returnBins[1] = binCounter;
      break;
    }
  }
  return;
}
