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
  void SetZBins(std::vector<double> ZBins);
  void SetMultBins(std::vector<int> MultBins);
  void SetPDGCodes(std::vector<int> PDGCodes);
  void SetNBinsHist(std::vector<int> NBins);
  void SetMinKRel(std::vector<double> minKRel);
  void SetMaxKRel(std::vector<double> maxKRel);
  void SetMixingDepth(int MixingDepth){fMixingDepth=MixingDepth;};

  bool GetDoMultBinning() {return fMultBinning;};
  std::vector<double> GetZVtxBins();
  int GetNZVtxBins(){return (fZVtxBins->GetEntries()-1);};
  std::vector<int> GetMultBins();
  int GetNMultBins(){return fMultBins->GetEntries();};
  std::vector<int> GetPDGCodes();
  int GetNParticles() {return fPDGParticleSpecies->GetEntries();};
  int GetNParticleCombinations();
  std::vector<int> GetNBinsHist();
  std::vector<double> GetMinKRel();
  std::vector<double> GetMaxKRel();
  int GetMixingDepth(){return fMixingDepth;};
 private:
  bool fMultBinning;            //
  TNtuple *fZVtxBins;           //
  TNtuple *fMultBins;           //
  TNtuple *fPDGParticleSpecies; //
  TNtuple *fNBinsHists;         //
  TNtuple *fMinK_rel;           //
  TNtuple *fMaxK_rel;           //
  int fMixingDepth;             //
  ClassDef(AliFemtoDreamCollConfig,1);
};

#endif /* ALIFEMTODREAMCOLLCONFIG_H_ */
