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
#include "AliFemtoDreamHigherPairMath.h"

//Class containing the array buffer of the different particle species for one
//Multiplicity bin
class AliFemtoDreamZVtxMultContainer {
 public:
  AliFemtoDreamZVtxMultContainer();
  AliFemtoDreamZVtxMultContainer(AliFemtoDreamCollConfig *conf);
  virtual ~AliFemtoDreamZVtxMultContainer();
  void PairParticlesSE(
      std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
      AliFemtoDreamHigherPairMath *HigherMath, int iMult, float cent);
  void PairParticlesME(
      std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
      AliFemtoDreamHigherPairMath *HigherMath, int iMult, float cent);
  void DeltaEtaDeltaPhi(int Hist, AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2, bool SEorME,
                        AliFemtoDreamCorrHists *ResultsHist, float relk);
  float ComputeDeltaEta(AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2);
  float ComputeDeltaPhi(AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2);
  void SetEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &Particles);
  TString ClassName() {
    return "zVtxMult Container";
  }
  ;
 private:
  std::vector<AliFemtoDreamPartContainer> fPartContainer;
  std::vector<int> fPDGParticleSpecies;
  std::vector<unsigned int> fWhichPairs;
//  std::vector<bool> fRejPairs;
//  bool fDoDeltaEtaDeltaPhiCut;
//  float fDeltaEtaMax;
//  float fDeltaPhiMax;
//  float fDeltaPhiEtaMax;

ClassDef(AliFemtoDreamZVtxMultContainer, 4)
  ;
};

#endif /* ALIFEMTODREAMZVTXMULTCONTAINER_H_ */
