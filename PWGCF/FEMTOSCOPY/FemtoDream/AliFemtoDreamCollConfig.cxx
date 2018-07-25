/*
 * AliFemtoDreamCollConfig.cxx
 *
 *  Created on: Sep 13, 2017
 *      Author: gu74req
 */
#include "TMath.h"
#include "AliFemtoDreamCollConfig.h"
ClassImp(AliFemtoDreamCollConfig)
AliFemtoDreamCollConfig::AliFemtoDreamCollConfig()
:TNamed()
,fMultBinning(false)
,fCentBinning(false)
,fkTBinning(false)
,fmTBinning(false)
,fMomentumResolution(false)
,fPhiEtaBinning(false)
,fdPhidEtaPlots(false)
,fMixedEventStatistics(true)
,fGetTheControlSampel(false)
,fInvMassPairs(false)
,fMinimalBookingME(false)
,fMinimalBookingSample(false)
,fNumberRadii(0)
,fZVtxBins(0)
,fMultBins(0)
,fPDGParticleSpecies(0)
,fNBinsHists(0)
,fMinK_rel(0)
,fMaxK_rel(0)
,fCentBins(0)
,fMixingDepth(0)
,fSpinningDepth(0)
,fkTCentrality(false)
,fMCCommonAncestor(false)
,fEst(AliFemtoDreamEvent::kSPD)
{
  //should not be used, since we need a name to deal with root objects
}

AliFemtoDreamCollConfig::AliFemtoDreamCollConfig(const char *name,
                                                 const char *title)
:TNamed(name,title)
,fMultBinning(false)
,fCentBinning(false)
,fkTBinning(false)
,fmTBinning(false)
,fMomentumResolution(false)
,fPhiEtaBinning(false)
,fdPhidEtaPlots(false)
,fMixedEventStatistics(true)
,fGetTheControlSampel(false)
,fInvMassPairs(false)
,fMinimalBookingME(false)
,fMinimalBookingSample(false)
,fNumberRadii(0)
,fMixingDepth(0)
,fSpinningDepth(0)
,fkTCentrality(false)
,fMCCommonAncestor(false)
,fEst(AliFemtoDreamEvent::kSPD)
{
  fZVtxBins=new TNtuple("ZBins","ZBins","zvtx");
  fMultBins=new TNtuple("MultBins","MultBins","mult");
  fPDGParticleSpecies=new TNtuple("PDGCodes","PDGCodes","PDGCodes");
  fNBinsHists=new TNtuple("NmbBins","NmbBins","NmbBins");
  fMinK_rel=new TNtuple("MinK_rel","MinK_rel","minkRel");
  fMaxK_rel=new TNtuple("MaxK_rel","MaxK_rel","maxkRel");
  fCentBins=new TNtuple("CentBins","CentBins","centBin");
}

AliFemtoDreamCollConfig::~AliFemtoDreamCollConfig() {
  delete fZVtxBins;
  delete fMultBins;
  delete fPDGParticleSpecies;
  delete fNBinsHists;
  delete fMinK_rel;
  delete fMaxK_rel;

}

void AliFemtoDreamCollConfig::SetZBins(std::vector<float> ZBins) {
  //Make sure to set the entries in ascending order!
  //Todo: maybe build in a check for this
  for (std::vector<float>::iterator it=ZBins.begin();it!=ZBins.end();++it) {
    fZVtxBins->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetZVtxBins() {
  //Make sure to set the entries in ascending order!
  std::vector<float> ZBins;
  float out=0;
  fZVtxBins->SetBranchAddress("zvtx",&out);
  for (int iBins=0;iBins<fZVtxBins->GetEntries();++iBins) {
    fZVtxBins->GetEntry(iBins);
    ZBins.push_back(out);
  }
  return ZBins;
}
void AliFemtoDreamCollConfig::SetMultBins(std::vector<int> MultBins) {
  //Make sure to set the entries in ascending order! The last bin to infinite
  //is implicit
  //Todo: maybe build in a check for this
  for (std::vector<int>::iterator it=MultBins.begin();it!=MultBins.end();++it)
  {
    fMultBins->Fill(*it);
  }
}
std::vector<int> AliFemtoDreamCollConfig::GetMultBins() {
  std::vector<int> MultBins;
  float out=0;
  fMultBins->SetBranchAddress("mult",&out);
  for (int iBins=0;iBins<fMultBins->GetEntries();++iBins) {
    fMultBins->GetEntry(iBins);
    MultBins.push_back(out);
  }
  return MultBins;
}
void AliFemtoDreamCollConfig::SetPDGCodes(std::vector<int> PDGCodes) {
  //the order needs to correspond the first particle array in your vector that
  //you hand over in the AliFemtoDreamPartCollection::SetEvent Method!
  for (std::vector<int>::iterator it=PDGCodes.begin();it!=PDGCodes.end();++it)
  {
    fPDGParticleSpecies->Fill(*it);
  }
}
std::vector<int> AliFemtoDreamCollConfig::GetPDGCodes() {
  std::vector<int> PDGCodes;
  float out=0;
  fPDGParticleSpecies->SetBranchAddress("PDGCodes",&out);
  for (int iBins=0;iBins<fPDGParticleSpecies->GetEntries();++iBins) {
    fPDGParticleSpecies->GetEntry(iBins);
    PDGCodes.push_back(out);
  }
  return PDGCodes;
}
int AliFemtoDreamCollConfig::GetNParticleCombinations(){
  //The possible number of combinations for pairing two particles species
  //with itself and all other species is for n species given by:
  //-Combinations within the same species n
  //-Combinations with all other species Binominal(n,2)
  int comb=fPDGParticleSpecies->GetEntries();
  if(comb>1){
    comb+=TMath::Binomial(comb,2);
  }
  return comb;
}

void AliFemtoDreamCollConfig::SetNBinsHist(std::vector<int> NBins) {
  //The way the histograms are assigned later is going to be for example for
  //4 different particle species X1,X2,X3,X4:
  //    X1  X2  X3  X4
  //X1  1   2   3   4
  //X2      5   6   7
  //X3          8   9
  //X4              10<-----Number of the Histogram=Position in input vector
  //Assign your binnig accordingly. X1 corresponds the first particle array
  //in your vector that you hand over in the
  //AliFemtoDreamPartCollection::SetEvent Method, X2 to the second and so on.
  //Same binning and ranges for Same Event and Mixed Event Distribution
  for (std::vector<int>::iterator it=NBins.begin();it!=NBins.end();++it)
  {
    fNBinsHists->Fill(*it);
  }
}
std::vector<int> AliFemtoDreamCollConfig::GetNBinsHist() {
  std::vector<int> NBinsHist;
  float out=0;
  fNBinsHists->SetBranchAddress("NmbBins",&out);
  for (int iBins=0;iBins<fNBinsHists->GetEntries();++iBins) {
    fNBinsHists->GetEntry(iBins);
    NBinsHist.push_back(out);
  }
  return NBinsHist;
}
void AliFemtoDreamCollConfig::SetMinKRel(std::vector<float> minKRel) {
  //See SetNBinsHist
  for (std::vector<float>::iterator it=minKRel.begin();it!=minKRel.end();++it)
  {
    fMinK_rel->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetMinKRel() {
  std::vector<float> MinKRel;
  float out=0;
  fMinK_rel->SetBranchAddress("minkRel",&out);
  for (int iBins=0;iBins<fMinK_rel->GetEntries();++iBins) {
    fMinK_rel->GetEntry(iBins);
    MinKRel.push_back(out);
  }
  return MinKRel;
}
void AliFemtoDreamCollConfig::SetMaxKRel(std::vector<float> maxKRel) {
  //See SetNBinsHist
  for (std::vector<float>::iterator it=maxKRel.begin();it!=maxKRel.end();++it)
  {
    fMaxK_rel->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetMaxKRel() {
  std::vector<float> MaxKRel;
  float out=0;
  fMaxK_rel->SetBranchAddress("maxkRel",&out);
  for (int iBins=0;iBins<fMaxK_rel->GetEntries();++iBins) {
    fMaxK_rel->GetEntry(iBins);
    MaxKRel.push_back(out);
  }
  return MaxKRel;
}
void AliFemtoDreamCollConfig::SetCentBins(std::vector<float> CentBins) {
  //Set Centrality Bins for the kT Centrality Binning
  for (std::vector<float>::iterator it=CentBins.begin();it!=CentBins.end();++it)
  {
    fCentBins->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetCentBins() {
  std::vector<float> CentBins;
  float out=0;
  fCentBins->SetBranchAddress("centBin",&out);
  for (int iBins=0;iBins<fCentBins->GetEntries();++iBins) {
    fCentBins->GetEntry(iBins);
    CentBins.push_back(out);
  }
  return CentBins;
}
