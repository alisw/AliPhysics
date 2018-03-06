/*
 * AliFemtoDreamCollConfig.h
 *
 *  Created on: Sep 13, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMCOLLCONFIG_H_
#define ALIFEMTODREAMCOLLCONFIG_H_
#include <vector>
#include "Rtypes.h"
#include "TNamed.h"
#include "TNtuple.h"

class AliFemtoDreamCollConfig : public TNamed {
 public:
  AliFemtoDreamCollConfig();
  AliFemtoDreamCollConfig(const char *name, const char *title);
  virtual ~AliFemtoDreamCollConfig();
  void SetMultBinning(bool doIt){fMultBinning=doIt;};
  void SetMomentumResolution(bool doIt){fMomentumResolution=doIt;};
  void SetZBins(std::vector<float> ZBins);
  void SetMultBins(std::vector<int> MultBins);
  void SetPDGCodes(std::vector<int> PDGCodes);
  void SetNBinsHist(std::vector<int> NBins);
  void SetMinKRel(std::vector<float> minKRel);
  void SetMaxKRel(std::vector<float> maxKRel);
  void SetMixingDepth(int MixingDepth){fMixingDepth=MixingDepth;};

  bool GetDoMultBinning() {return fMultBinning;};
  bool GetDoMomResolution() {return fMomentumResolution;};
  std::vector<float> GetZVtxBins();
  int GetNZVtxBins(){return (fZVtxBins->GetEntries()-1);};
  std::vector<int> GetMultBins();
  int GetNMultBins(){return fMultBins->GetEntries();};
  std::vector<int> GetPDGCodes();
  int GetNParticles() {return fPDGParticleSpecies->GetEntries();};
  int GetNParticleCombinations();
  std::vector<int> GetNBinsHist();
  std::vector<float> GetMinKRel();
  std::vector<float> GetMaxKRel();
  int GetMixingDepth(){return fMixingDepth;};
 private:
  bool fMultBinning;            //
  bool fMomentumResolution;     //
  TNtuple *fZVtxBins;           //
  TNtuple *fMultBins;           //
  TNtuple *fPDGParticleSpecies; //
  TNtuple *fNBinsHists;         //
  TNtuple *fMinK_rel;           //
  TNtuple *fMaxK_rel;           //
  int fMixingDepth;             //
  ClassDef(AliFemtoDreamCollConfig,2);
};

#endif /* ALIFEMTODREAMCOLLCONFIG_H_ */
