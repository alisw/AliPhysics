/*
 * AliFemtoPPbpbLamPartCollection.h
 *
 *  Created on: Aug 29, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMPARTCOLLECTION_H_
#define ALIFEMTODREAMPARTCOLLECTION_H_
#include <deque>
#include <vector>
#include "Rtypes.h"
#include "TList.h"

#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamCorrHists.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliFemtoDreamZVtxMultContainer.h"
//Class containing all the different multiplicity containers for all the
//particles
class AliFemtoDreamPartCollection {
 public:
  AliFemtoDreamPartCollection();
  AliFemtoDreamPartCollection(const AliFemtoDreamPartCollection& coll);
  AliFemtoDreamPartCollection(AliFemtoDreamCollConfig *conf,
                              bool MinimalBooking);
  AliFemtoDreamPartCollection& operator=(
      const AliFemtoDreamPartCollection& coll);
  virtual ~AliFemtoDreamPartCollection();
  void SetEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
                float ZVtx, float Mult, float cent);
  void SetEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
                AliFemtoDreamEvent* evt);
  void PrintEvent(int ZVtx, int Mult);
  TList* GetHistList() {
    return fHigherMath->GetHistList();
  }
  ;
  TList* GetQAList() {
    return fHigherMath->GetQAHists();
  }
  ;
  TString ClassName() {
    return "ParticleCollection";
  }
  ;
  void FindBin(float ZVtxPos, float Multiplicity, int *returnBins);
 private:
  AliFemtoDreamHigherPairMath* fHigherMath;
  unsigned int fNSpecies;
  std::vector<std::vector<AliFemtoDreamZVtxMultContainer>> fZVtxMultBuffer;
  std::vector<float> fValuesZVtxBins;
  std::vector<int> fValuesMultBins;
  ClassDef(AliFemtoDreamPartCollection,2);
};

#endif /* ALIFEMTODREAMPARTCOLLECTION_H_ */
