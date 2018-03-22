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
  void SetkTBinning(bool doIt){fkTBinning=doIt;};
  void SetmTBinning(bool doIt){fmTBinning=doIt;};
  void SetMomentumResolution(bool doIt){fMomentumResolution=doIt;};
  void SetkTCentralityBinning(bool doIt){fkTCentrality=doIt;};
  void SetPhiEtaBinnign(bool doIt){
    fPhiEtaBinning=doIt;fNumberRadii=9;
  };
  void SetZBins(std::vector<float> ZBins);
  void SetMultBins(std::vector<int> MultBins);
  void SetPDGCodes(std::vector<int> PDGCodes);
  void SetNBinsHist(std::vector<int> NBins);
  void SetMinKRel(std::vector<float> minKRel);
  void SetMaxKRel(std::vector<float> maxKRel);
  void SetCentBins(std::vector<float> CentBins);
  void SetMixingDepth(int MixingDepth){fMixingDepth=MixingDepth;};
  void SetSpinningDepth(int SpinningDepth){fSpinningDepth=SpinningDepth;};
  void SetSECommonAncestor(bool doit) {fMCCommonAncestor=doit;};

  bool GetDoMultBinning() {return fMultBinning;};
  bool GetDokTBinning() {return fkTBinning;};
  bool GetDomTBinning() {return fmTBinning;};
  bool GetDoMomResolution() {return fMomentumResolution;};
  bool GetDoPhiEtaBinning() {return fPhiEtaBinning;};
  bool GetDokTCentralityBinning() {return fkTCentrality;};
  bool GetDoSECommonAncestor() {return fMCCommonAncestor;};

  int GetNRadii() {return fNumberRadii;};
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
  std::vector<float> GetCentBins();
  int GetMixingDepth(){return fMixingDepth;};
  int GetSpinningDepth(){return fSpinningDepth;};
 private:
  bool fMultBinning;            //
  bool fkTBinning;            //
  bool fmTBinning;            //
  bool fMomentumResolution;     //
  bool fPhiEtaBinning;          //
  int fNumberRadii;             //
  TNtuple *fZVtxBins;           //
  TNtuple *fMultBins;           //
  TNtuple *fPDGParticleSpecies; //
  TNtuple *fNBinsHists;         //
  TNtuple *fMinK_rel;           //
  TNtuple *fMaxK_rel;           //
  TNtuple *fCentBins;           //
  int fMixingDepth;             //
  int fSpinningDepth;			//
  bool fkTCentrality;           //
  bool fMCCommonAncestor;       // Setter used in MC Only to obtain the SE distribution for common ancestor and non common ancestor
  ClassDef(AliFemtoDreamCollConfig,3);
};

#endif /* ALIFEMTODREAMCOLLCONFIG_H_ */
