/*
 * AliFemtoPPbpbLamZVtxMultContainer.h
 *
 *  Created on: Aug 30, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMZVTXMULTCONTAINER_H_
#define ALIFEMTODREAMZVTXMULTCONTAINER_H_
#include <vector>
#include "Rtypes.h"

#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamCorrHists.h"
#include "AliFemtoDreamPartContainer.h"
//Class containing the array buffer of the different particle species for one
//Multiplicity bin
class AliFemtoDreamZVtxMultContainer {
 public:
  AliFemtoDreamZVtxMultContainer();
  AliFemtoDreamZVtxMultContainer(AliFemtoDreamCollConfig *conf);
  virtual ~AliFemtoDreamZVtxMultContainer();
  void PairParticlesSE(
      std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
      AliFemtoDreamCorrHists *ResultsHist,int iMult);
  void PairParticlesME(
      std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
      AliFemtoDreamCorrHists *ResultsHist,int iMult);
  void SetEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &Particles);
  TString ClassName() {return "zVtxMult Container";};
 private:
  double RelativePairMomentum(TVector3 Part1Momentum,int PDGPart1,
                              TVector3 Part2Momentum,int PDGPart2);
  std::vector<AliFemtoDreamPartContainer> fPartContainer;
  std::vector<int> fPDGParticleSpecies;
  ClassDef(AliFemtoDreamZVtxMultContainer,1);
};

#endif /* ALIFEMTODREAMZVTXMULTCONTAINER_H_ */
