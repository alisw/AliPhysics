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
    : TNamed(),
      fMultBinning(false),
      fCentBinning(false),
      fkTBinning(false),
      fmTBinning(false),
      fPtQA(false),
      fMomentumResolution(false),
      fPhiEtaBinning(false),
      fdPhidEtaPlots(false),
      fdPhidEtaPlotsSmallK(false),
      fMixedEventStatistics(true),
      fGetTheControlSampel(false),
      fStravinsky(false),
      fInvMassPairs(false),
      fMinimalBookingME(false),
      fMinimalBookingSample(false),
      fNumberRadii(0),
      fZVtxBins(0),
      fMultBins(0),
      fPDGParticleSpecies(0),
      fNBinsHists(0),
      fMinK_rel(0),
      fMaxK_rel(0),
      fCentBins(0),
      fmTBins(0),
      fWhichQAPairs(0),
      fClosePairRej(0),
      fMixingDepth(0),
      fSpinningDepth(0),
      fkTCentrality(false),
      fmTdEtadPhi(false),
      fEst(AliFemtoDreamEvent::kSPD),
      fDeltaEtaMax(0.f),
      fDeltaPhiMax(0.f),
      fDoDeltaEtaDeltaPhiCut(false) {
  //should not be used, since we need a name to deal with root objects
}

AliFemtoDreamCollConfig::AliFemtoDreamCollConfig(
    const AliFemtoDreamCollConfig& config)
    : TNamed(config),
      fMultBinning(config.fMultBinning),
      fCentBinning(config.fCentBinning),
      fkTBinning(config.fkTBinning),
      fmTBinning(config.fmTBinning),
      fPtQA(config.fPtQA),
      fMomentumResolution(config.fMomentumResolution),
      fPhiEtaBinning(config.fPhiEtaBinning),
      fdPhidEtaPlots(config.fdPhidEtaPlots),
      fdPhidEtaPlotsSmallK(config.fdPhidEtaPlotsSmallK),
      fMixedEventStatistics(config.fMixedEventStatistics),
      fGetTheControlSampel(config.fGetTheControlSampel),
      fStravinsky(config.fStravinsky),
      fInvMassPairs(config.fInvMassPairs),
      fMinimalBookingME(config.fMinimalBookingME),
      fMinimalBookingSample(config.fMinimalBookingSample),
      fNumberRadii(config.fNumberRadii),
      fZVtxBins(config.fZVtxBins),
      fMultBins(config.fMultBins),
      fPDGParticleSpecies(config.fPDGParticleSpecies),
      fNBinsHists(config.fNBinsHists),
      fMinK_rel(config.fMinK_rel),
      fMaxK_rel(config.fMaxK_rel),
      fCentBins(config.fCentBins),
      fmTBins(config.fmTBins),
      fWhichQAPairs(config.fWhichQAPairs),
      fClosePairRej(config.fClosePairRej),
      fMixingDepth(config.fMixingDepth),
      fSpinningDepth(config.fSpinningDepth),
      fkTCentrality(config.fkTCentrality),
      fmTdEtadPhi(config.fmTdEtadPhi),
      fEst(config.fEst),
      fDeltaEtaMax(config.fDeltaEtaMax),
      fDeltaPhiMax(config.fDeltaPhiMax),
      fDoDeltaEtaDeltaPhiCut(config.fDoDeltaEtaDeltaPhiCut) {
}

AliFemtoDreamCollConfig::AliFemtoDreamCollConfig(const char *name,
                                                 const char *title)
    : TNamed(name, title),
      fMultBinning(false),
      fCentBinning(false),
      fkTBinning(false),
      fmTBinning(false),
      fPtQA(false),
      fMomentumResolution(false),
      fPhiEtaBinning(false),
      fdPhidEtaPlots(false),
      fdPhidEtaPlotsSmallK(false),
      fMixedEventStatistics(true),
      fGetTheControlSampel(false),
      fStravinsky(false),
      fInvMassPairs(false),
      fMinimalBookingME(false),
      fMinimalBookingSample(false),
      fNumberRadii(0),
      fZVtxBins(nullptr),
      fMultBins(nullptr),
      fPDGParticleSpecies(nullptr),
      fNBinsHists(nullptr),
      fMinK_rel(nullptr),
      fMaxK_rel(nullptr),
      fCentBins(nullptr),
      fmTBins(nullptr),
      fWhichQAPairs(nullptr),
      fClosePairRej(nullptr),
      fMixingDepth(0),
      fSpinningDepth(0),
      fkTCentrality(false),
      fmTdEtadPhi(false),
      fEst(AliFemtoDreamEvent::kSPD),
      fDeltaEtaMax(0.f),
      fDeltaPhiMax(0.f),
      fDoDeltaEtaDeltaPhiCut(false) {
  fZVtxBins = new TNtuple("ZBins", "ZBins", "zvtx");
  fMultBins = new TNtuple("MultBins", "MultBins", "mult");
  fPDGParticleSpecies = new TNtuple("PDGCodes", "PDGCodes", "PDGCodes");
  fNBinsHists = new TNtuple("NmbBins", "NmbBins", "NmbBins");
  fMinK_rel = new TNtuple("MinK_rel", "MinK_rel", "minkRel");
  fMaxK_rel = new TNtuple("MaxK_rel", "MaxK_rel", "maxkRel");
  fCentBins = new TNtuple("CentBins", "CentBins", "centBin");
  fmTBins = new TNtuple("mTBins", "mTBins", "mTBin");
  fWhichQAPairs = new TNtuple("DoPairs", "DoPairs", "DoPair");
  fClosePairRej = new TNtuple("PairRej", "RejPairs", "RejPair");
}
AliFemtoDreamCollConfig& AliFemtoDreamCollConfig::operator=(
    const AliFemtoDreamCollConfig& config) {
  if (this != &config) {
    TNamed::operator=(config);
    this->fMultBinning = config.fMultBinning;
    this->fCentBinning = config.fCentBinning;
    this->fkTBinning = config.fkTBinning;
    this->fmTBinning = config.fmTBinning;
    this->fPtQA = config.fPtQA;
    this->fMomentumResolution = config.fMomentumResolution;
    this->fPhiEtaBinning = config.fPhiEtaBinning;
    this->fdPhidEtaPlots = config.fdPhidEtaPlots;
    this->fdPhidEtaPlotsSmallK = config.fdPhidEtaPlotsSmallK;
    this->fMixedEventStatistics = config.fMixedEventStatistics;
    this->fGetTheControlSampel = config.fGetTheControlSampel;
    this->fInvMassPairs = config.fInvMassPairs;
    this->fMinimalBookingME = config.fMinimalBookingME;
    this->fMinimalBookingSample = config.fMinimalBookingSample;
    this->fNumberRadii = config.fNumberRadii;
    this->fZVtxBins = config.fZVtxBins;
    this->fMultBins = config.fMultBins;
    this->fPDGParticleSpecies = config.fPDGParticleSpecies;
    this->fNBinsHists = config.fNBinsHists;
    this->fMinK_rel = config.fMinK_rel;
    this->fMaxK_rel = config.fMaxK_rel;
    this->fCentBins = config.fCentBins;
    this->fmTBins = config.fmTBins;
    this->fWhichQAPairs = config.fWhichQAPairs;
    this->fClosePairRej = config.fClosePairRej;
    this->fMixingDepth = config.fMixingDepth;
    this->fSpinningDepth = config.fSpinningDepth;
    this->fkTCentrality = config.fkTCentrality;
    this->fmTdEtadPhi = config.fmTdEtadPhi;
    this->fEst = config.fEst;
    this->fDeltaEtaMax = config.fDeltaEtaMax;
    this->fDeltaPhiMax = config.fDeltaPhiMax;
    this->fDoDeltaEtaDeltaPhiCut = config.fDoDeltaEtaDeltaPhiCut;
  }
  return *this;
}

AliFemtoDreamCollConfig::~AliFemtoDreamCollConfig() {
  delete fZVtxBins;
  delete fMultBins;
  delete fPDGParticleSpecies;
  delete fNBinsHists;
  delete fMinK_rel;
  delete fMaxK_rel;
  delete fCentBins;
  delete fmTBins;
  delete fWhichQAPairs;
  delete fClosePairRej;

}

void AliFemtoDreamCollConfig::SetZBins(std::vector<float> ZBins) {
  //Make sure to set the entries in ascending order!
  //Todo: maybe build in a check for this
  for (std::vector<float>::iterator it = ZBins.begin(); it != ZBins.end();
      ++it) {
    fZVtxBins->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetZVtxBins() {
  //Make sure to set the entries in ascending order!
  std::vector<float> ZBins;
  float out = 0;
  fZVtxBins->SetBranchAddress("zvtx", &out);
  for (int iBins = 0; iBins < fZVtxBins->GetEntries(); ++iBins) {
    fZVtxBins->GetEntry(iBins);
    ZBins.push_back(out);
  }
  return ZBins;
}
void AliFemtoDreamCollConfig::SetMultBins(std::vector<int> MultBins) {
  //Make sure to set the entries in ascending order! The last bin to infinite
  //is implicit
  //Todo: maybe build in a check for this
  for (std::vector<int>::iterator it = MultBins.begin(); it != MultBins.end();
      ++it) {
    fMultBins->Fill(*it);
  }
}
std::vector<int> AliFemtoDreamCollConfig::GetMultBins() {
  std::vector<int> MultBins;
  float out = 0;
  fMultBins->SetBranchAddress("mult", &out);
  for (int iBins = 0; iBins < fMultBins->GetEntries(); ++iBins) {
    fMultBins->GetEntry(iBins);
    MultBins.push_back(out);
  }
  return MultBins;
}
void AliFemtoDreamCollConfig::SetPDGCodes(std::vector<int> PDGCodes) {
  //the order needs to correspond the first particle array in your vector that
  //you hand over in the AliFemtoDreamPartCollection::SetEvent Method!
  for (std::vector<int>::iterator it = PDGCodes.begin(); it != PDGCodes.end();
      ++it) {
    fPDGParticleSpecies->Fill(*it);
  }
}
std::vector<int> AliFemtoDreamCollConfig::GetPDGCodes() {
  std::vector<int> PDGCodes;
  float out = 0;
  fPDGParticleSpecies->SetBranchAddress("PDGCodes", &out);
  for (int iBins = 0; iBins < fPDGParticleSpecies->GetEntries(); ++iBins) {
    fPDGParticleSpecies->GetEntry(iBins);
    PDGCodes.push_back(out);
  }
  return PDGCodes;
}
int AliFemtoDreamCollConfig::GetNParticleCombinations() {
  //The possible number of combinations for pairing two particles species
  //with itself and all other species is for n species given by:
  //-Combinations within the same species n
  //-Combinations with all other species Binominal(n,2)
  int comb = fPDGParticleSpecies->GetEntries();
  if (comb > 1) {
    comb += TMath::Binomial(comb, 2);
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
  for (std::vector<int>::iterator it = NBins.begin(); it != NBins.end(); ++it) {
    fNBinsHists->Fill(*it);
  }
}
std::vector<int> AliFemtoDreamCollConfig::GetNBinsHist() {
  std::vector<int> NBinsHist;
  float out = 0;
  fNBinsHists->SetBranchAddress("NmbBins", &out);
  for (int iBins = 0; iBins < fNBinsHists->GetEntries(); ++iBins) {
    fNBinsHists->GetEntry(iBins);
    NBinsHist.push_back(out);
  }
  return NBinsHist;
}
void AliFemtoDreamCollConfig::SetMinKRel(std::vector<float> minKRel) {
  //See SetNBinsHist
  for (std::vector<float>::iterator it = minKRel.begin(); it != minKRel.end();
      ++it) {
    fMinK_rel->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetMinKRel() {
  std::vector<float> MinKRel;
  float out = 0;
  fMinK_rel->SetBranchAddress("minkRel", &out);
  for (int iBins = 0; iBins < fMinK_rel->GetEntries(); ++iBins) {
    fMinK_rel->GetEntry(iBins);
    MinKRel.push_back(out);
  }
  return MinKRel;
}
void AliFemtoDreamCollConfig::SetMaxKRel(std::vector<float> maxKRel) {
  //See SetNBinsHist
  for (std::vector<float>::iterator it = maxKRel.begin(); it != maxKRel.end();
      ++it) {
    fMaxK_rel->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetMaxKRel() {
  std::vector<float> MaxKRel;
  float out = 0;
  fMaxK_rel->SetBranchAddress("maxkRel", &out);
  for (int iBins = 0; iBins < fMaxK_rel->GetEntries(); ++iBins) {
    fMaxK_rel->GetEntry(iBins);
    MaxKRel.push_back(out);
  }
  return MaxKRel;
}
void AliFemtoDreamCollConfig::SetCentBins(std::vector<float> CentBins) {
  //Set Centrality Bins for the kT Centrality Binning
  for (std::vector<float>::iterator it = CentBins.begin(); it != CentBins.end();
      ++it) {
    fCentBins->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetCentBins() {
  std::vector<float> CentBins;
  float out = 0;
  fCentBins->SetBranchAddress("centBin", &out);
  for (int iBins = 0; iBins < fCentBins->GetEntries(); ++iBins) {
    fCentBins->GetEntry(iBins);
    CentBins.push_back(out);
  }
  return CentBins;
}
void AliFemtoDreamCollConfig::SetmTdEtadPhiBins(std::vector<float> mTBins) {
  //Set Bins for the deta dphi mT Binning
  fmTdEtadPhi = true;
  for (std::vector<float>::iterator it = mTBins.begin(); it != mTBins.end();
      ++it) {
    fmTBins->Fill(*it);
  }
}
std::vector<float> AliFemtoDreamCollConfig::GetmTBins() {
  std::vector<float> mTBins;
  float out = 0;
  fmTBins->SetBranchAddress("mTBin", &out);
  for (int iBins = 0; iBins < fmTBins->GetEntries(); ++iBins) {
    fmTBins->GetEntry(iBins);
    mTBins.push_back(out);
  }
  return mTBins;
}
void AliFemtoDreamCollConfig::SetExtendedQAPairs(std::vector<int> whichPairs) {
  // Decider for if one wants plots like mT, kT, TrackSplitting etc. for a pair
  // 0 means no, > 0 means yes. In particular for track splitting one can steer
  // which combinations to check:
  // 12 for example means: take only the first track from particle 1 of the pair
  // and check it against the 2 tracks of particle 2 ( if it is for example a v0
  // Lambda candidate).
  // Indices follow the scheme explained in SetNBinsHist.
  for (auto it : whichPairs) {
    fWhichQAPairs->Fill(it);
  }
  //how to select which pairs?
}

std::vector<unsigned int> AliFemtoDreamCollConfig::GetWhichPairs() {
  std::vector<unsigned int> Pairs;
  float out = 0;
  if (fWhichQAPairs->GetEntries() == 0) {
    AliWarning("===========================================");
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("No Pair QA Specified, setting all to false ");
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("===========================================");
    for (int iQA = 0; iQA < this->GetNParticleCombinations(); ++iQA) {
      Pairs.push_back(0);
    }
  } else if (fWhichQAPairs->GetEntries() != this->GetNParticleCombinations()) {
    AliFatal("Not all Pairs have a specified QA Behaviour, terminating \n");
  } else {
    fWhichQAPairs->SetBranchAddress("DoPair", &out);
    for (int iBins = 0; iBins < fWhichQAPairs->GetEntries(); ++iBins) {
      fWhichQAPairs->GetEntry(iBins);
      Pairs.push_back(TMath::Abs(out));
    }
  }
  return Pairs;
}

std::vector<float> AliFemtoDreamCollConfig::GetStandardmTBins() {
  std::vector<float> mTBins;
  mTBins.push_back(1.2);
  mTBins.push_back(2.0);
  mTBins.push_back(4.5);
  mTBins.push_back(999.);
  return mTBins;
}

std::vector<int> AliFemtoDreamCollConfig::GetStandardPairs() {
  std::vector<int> pairs;
  pairs.push_back(11);        // p p
  pairs.push_back(0);         // p barp
  pairs.push_back(12);        // p Lambda
  pairs.push_back(0);         // p barLambda
  pairs.push_back(13);         // p Xi
  pairs.push_back(0);         // p barXi
  pairs.push_back(11);        // barp barp
  pairs.push_back(0);         // barp Lambda
  pairs.push_back(12);        // barp barLambda
  pairs.push_back(0);         // barp Xi
  pairs.push_back(13);         // barp barXi
  pairs.push_back(0);         // Lambda Lambda
  pairs.push_back(0);         // Lambda barLambda
  pairs.push_back(0);         // Lambda Xi
  pairs.push_back(0);         // Lambda barXi
  pairs.push_back(0);         // barLambda barLamb
  pairs.push_back(0);         // barLambda Xi
  pairs.push_back(0);         // barLambda barXi
  pairs.push_back(0);         // Xi Xi
  pairs.push_back(0);         // Xi barXi
  pairs.push_back(0);         // barXi barXi
  return pairs;
}

void AliFemtoDreamCollConfig::SetClosePairRejection(
    std::vector<bool> whichPairs) {
  // Decider for if one wants plots like mT, kT, TrackSplitting etc. for a pair
  // 0 means no, > 0 means yes. In particular for track splitting one can steer
  // which combinations to check:
  // 12 for example means: take only the first track from particle 1 of the pair
  // and check it against the 2 tracks of particle 2 ( if it is for example a v0
  // Lambda candidate).
  // Indices follow the scheme explained in SetNBinsHist.
  for (auto it : whichPairs) {
    fClosePairRej->Fill(it ? 1 : 0);
  }
  //how to select which pairs?
}

std::vector<bool> AliFemtoDreamCollConfig::GetClosePairRej() {
  std::vector<bool> Pairs;
  float out = 0;
  if (fClosePairRej->GetEntries() == 0) {
    AliWarning("=========================================================");
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("!No Close Pair Rejection Specified, setting all to false!");
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("=========================================================");
    for (int iRej = 0; iRej < this->GetNParticleCombinations(); ++iRej) {
      Pairs.push_back(false);
    }
  } else if (fClosePairRej->GetEntries() != this->GetNParticleCombinations()) {
    AliFatal("Not all Pairs have a specified QA Behaviour, terminating \n");
  } else {
    fClosePairRej->SetBranchAddress("RejPair", &out);
    for (int iRej = 0; iRej < fClosePairRej->GetEntries(); ++iRej) {
      fClosePairRej->GetEntry(iRej);
      Pairs.push_back(TMath::Abs(out) < 1e-6 ? false : true);
      std::cout << "Close Pair Rejection for Pair " << iRej << " is "
                << (TMath::Abs(out) < 1e-6 ? "deactivated" : "activated")
                << std::endl;
    }
  }
  return Pairs;
}

std::vector<bool> AliFemtoDreamCollConfig::GetStandardPairRejection() {
  std::vector<bool> pairs;
  pairs.push_back(true);        // p p
  pairs.push_back(true);         // p barp
  pairs.push_back(false);        // p Lambda
  pairs.push_back(false);         // p barLambda
  pairs.push_back(false);         // p Xi
  pairs.push_back(false);         // p barXi
  pairs.push_back(true);        // barp barp
  pairs.push_back(false);         // barp Lambda
  pairs.push_back(false);        // barp barLambda
  pairs.push_back(false);         // barp Xi
  pairs.push_back(false);         // barp barXi
  pairs.push_back(false);         // Lambda Lambda
  pairs.push_back(false);         // Lambda barLambda
  pairs.push_back(false);         // Lambda Xi
  pairs.push_back(false);         // Lambda barXi
  pairs.push_back(false);         // barLambda barLamb
  pairs.push_back(false);         // barLambda Xi
  pairs.push_back(false);         // barLambda barXi
  pairs.push_back(false);         // Xi Xi
  pairs.push_back(false);         // Xi barXi
  pairs.push_back(false);         // barXi barXi
  return pairs;
}

std::vector<bool> AliFemtoDreamCollConfig::GetAllPairRejection() {
  std::vector<bool> pairs;
  pairs.push_back(true);        // p p
  pairs.push_back(true);         // p barp
  pairs.push_back(true);        // p Lambda
  pairs.push_back(true);         // p barLambda
  pairs.push_back(true);         // p Xi
  pairs.push_back(true);         // p barXi
  pairs.push_back(true);        // barp barp
  pairs.push_back(true);         // barp Lambda
  pairs.push_back(true);        // barp barLambda
  pairs.push_back(true);         // barp Xi
  pairs.push_back(true);         // barp barXi
  pairs.push_back(true);         // Lambda Lambda
  pairs.push_back(true);         // Lambda barLambda
  pairs.push_back(true);         // Lambda Xi
  pairs.push_back(true);         // Lambda barXi
  pairs.push_back(true);         // barLambda barLamb
  pairs.push_back(true);         // barLambda Xi
  pairs.push_back(true);         // barLambda barXi
  pairs.push_back(true);         // Xi Xi
  pairs.push_back(true);         // Xi barXi
  pairs.push_back(true);         // barXi barXi
  return pairs;
}
