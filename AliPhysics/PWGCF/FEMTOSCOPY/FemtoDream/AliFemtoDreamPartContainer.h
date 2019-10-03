/*
 * AliFemtoPPbpbLamPartContainer.h
 *
 *  Created on: Aug 30, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMPARTCONTAINER_H_
#define ALIFEMTODREAMPARTCONTAINER_H_
#include <deque>
#include <vector>
#include "Rtypes.h"

#include "AliFemtoDreamBasePart.h"

//Class Containing the Particles from previous Events up to a certain mixing
//depth for one Particle Species and Mult/ZVtx Bin
//ZVtx bin.
class AliFemtoDreamPartContainer {
 public:
  AliFemtoDreamPartContainer();
  AliFemtoDreamPartContainer(int MixingDepth);
  AliFemtoDreamPartContainer& operator=(const AliFemtoDreamPartContainer& obj);
  virtual ~AliFemtoDreamPartContainer();
  void PrintLastEvent();
  void SetEvent(std::vector<AliFemtoDreamBasePart> &Particles);
  std::deque<std::vector<AliFemtoDreamBasePart>> GetEventBuffer() const {
    return fPartBuffer;
  }
  ;
  std::vector<AliFemtoDreamBasePart> &GetEvent(int Depth);
  unsigned int GetMixingDepth() const {
    return fPartBuffer.size();
  }
  ;
 private:
  std::deque<std::vector<AliFemtoDreamBasePart>> fPartBuffer;
  unsigned int fMixingDepth;ClassDef(AliFemtoDreamPartContainer,2)
  ;
};

#endif /* ALIFEMTODREAMPARTCONTAINER_H_ */
