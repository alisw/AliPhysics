/*
 * AliFemtoPPbpbLamPartContainer.cxx
 *
 *  Created on: Aug 30, 2017
 *      Author: gu74req
 */

#include <iostream>
#include "AliFemtoDreamPartContainer.h"
#include "TLorentzVector.h"
#include "TVector3.h"
ClassImp(AliFemtoDreamPartContainer)
AliFemtoDreamPartContainer::AliFemtoDreamPartContainer()
    : fPartBuffer(),
      fMixingDepth(0) {

}

AliFemtoDreamPartContainer::AliFemtoDreamPartContainer(int MixingDepth)
    : fPartBuffer(),
      fMixingDepth(MixingDepth) {

}

AliFemtoDreamPartContainer& AliFemtoDreamPartContainer::operator=(
    const AliFemtoDreamPartContainer&obj) {
  if (this == &obj) {
    return *this;
  }
//  std::deque<std::deque<AliFemtoDreamBasePart>>::iterator copyIter;
//  std::deque<std::deque<AliFemtoDreamBasePart>>::iterator objIter;
//  for (objIter=obj.fPartBuffer.begin();objIter!=obj.fPartBuffer.end();
//      ++objIter)
//  {
//
//  }
  this->fMixingDepth = obj.fMixingDepth;
  this->fPartBuffer = obj.fPartBuffer;
  return (*this);
}

AliFemtoDreamPartContainer::~AliFemtoDreamPartContainer() {
}

void AliFemtoDreamPartContainer::SetEvent(
    std::vector<AliFemtoDreamBasePart> &Particles) {
  if (!(fPartBuffer.size() < fMixingDepth)) {
//    std::cout << "Popping Front" << std::endl;
    fPartBuffer.pop_front();
  }
  fPartBuffer.push_back(Particles);
//  std::cout << "PartBuffer Size: "<<fPartBuffer.size()<<'\t'<<"Input Size: "
//      << Particles.size() << '\n';
  return;
}

void AliFemtoDreamPartContainer::PrintLastEvent() {
  for (std::deque<std::vector<AliFemtoDreamBasePart>>::iterator itEvt =
      fPartBuffer.begin(); itEvt != fPartBuffer.end(); ++itEvt) {
    std::cout << "Printing Last Event with size: " << itEvt->size() << '\n';
    for (std::vector<AliFemtoDreamBasePart>::iterator itPart = itEvt->begin();
        itPart != itEvt->end(); ++itPart) {
      TVector3 P(itPart->GetMomentum());
      std::cout << "Px: " << P.X() << '\t' << "Py: " << P.Y() << '\t' << "Pz: "
                << P.Z() << std::endl;
    }
  }
}
std::vector<AliFemtoDreamBasePart> &AliFemtoDreamPartContainer::GetEvent(
    int Depth) {
  std::deque<std::vector<AliFemtoDreamBasePart>>::iterator itEvt = fPartBuffer
      .begin() + Depth;
  return *itEvt;
}
