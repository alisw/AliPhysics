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
#include "AliFemtoDreamZVtxMultContainer.h"
//Class containing all the different multiplicity containers for all the
//particles
class AliFemtoDreamPartCollection {
 public:
  AliFemtoDreamPartCollection();
  AliFemtoDreamPartCollection(AliFemtoDreamCollConfig *conf);
  virtual ~AliFemtoDreamPartCollection();
  void SetEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
                double ZVtx,double Mult);
  void PrintEvent(int ZVtx,int Mult);
  TList* GetHistList(){return fResults->GetHistList();};
  TList* GetQAList(){return fResults->GetQAHists();};
 private:
  void FindBin(double ZVtxPos,double Multiplicity,int *returnBins);
  AliFemtoDreamCorrHists *fResults;
  unsigned int fNSpecies;
  std::vector<std::vector<AliFemtoDreamZVtxMultContainer>> fZVtxMultBuffer;
  std::vector<double> fValuesZVtxBins;
  std::vector<int> fValuesMultBins;
  ClassDef(AliFemtoDreamPartCollection,1);
};

#endif /* ALIFEMTODREAMPARTCOLLECTION_H_ */
