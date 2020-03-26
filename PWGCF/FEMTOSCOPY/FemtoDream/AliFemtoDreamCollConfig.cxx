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
      fkTandMultBinning(false),
      fPtQA(false),
      fMassQA(false),
      fMomentumResolution(false),
      fPhiEtaBinning(false),
      fdPhidEtaPlots(false),
      fdPhidEtaPlotsSmallK(false),
      fMixedEventStatistics(true),
      fGetTheControlSampel(false),
      fMode(AliFemtoDreamCollConfig::kNone),
      fMinimalBookingME(false),
      fMinimalBookingSample(false),
      fNumberRadii(0),
      fZVtxBins(),
      fMultBins(),
      fPDGParticleSpecies(),
      fNBinsHists(),
      fMinK_rel(),
      fMaxK_rel(),
      fCentBins(),
      fmTBins(),
      fWhichQAPairs(),
      fClosePairRej(),
      fMixingDepth(0),
      fSpinningDepth(0),
      fCorrelationRange(0.),
      fkTCentrality(false),
      fmTdEtadPhi(false),
      fmTMultBinning(false),
      fEst(AliFemtoDreamEvent::kSPD),
      fAncestors(false),
      fDeltaEtaMax(0.f),
      fDeltaPhiMax(0.f),
      fDoDeltaEtaDeltaPhiCut(false),
      fCoutVariables(false) {
  //should not be used, since we need a name to deal with root objects
}

AliFemtoDreamCollConfig::AliFemtoDreamCollConfig(
    const AliFemtoDreamCollConfig& config)
    : TNamed(config),
      fMultBinning(config.fMultBinning),
      fCentBinning(config.fCentBinning),
      fkTBinning(config.fkTBinning),
      fmTBinning(config.fmTBinning),
      fkTandMultBinning(config.fkTandMultBinning),
      fPtQA(config.fPtQA),
      fMassQA(config.fMassQA),
      fMomentumResolution(config.fMomentumResolution),
      fPhiEtaBinning(config.fPhiEtaBinning),
      fdPhidEtaPlots(config.fdPhidEtaPlots),
      fdPhidEtaPlotsSmallK(config.fdPhidEtaPlotsSmallK),
      fMixedEventStatistics(config.fMixedEventStatistics),
      fGetTheControlSampel(config.fGetTheControlSampel),
      fMode(config.fMode),
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
      fCorrelationRange(config.fCorrelationRange),
      fkTCentrality(config.fkTCentrality),
      fmTdEtadPhi(config.fmTdEtadPhi),
      fmTMultBinning(config.fmTMultBinning),
      fEst(config.fEst),
      fAncestors(config.fAncestors),
      fDeltaEtaMax(config.fDeltaEtaMax),
      fDeltaPhiMax(config.fDeltaPhiMax),
      fDoDeltaEtaDeltaPhiCut(config.fDoDeltaEtaDeltaPhiCut),
      fCoutVariables(config.fCoutVariables) {
}

AliFemtoDreamCollConfig::AliFemtoDreamCollConfig(const char *name,
                                                 const char *title,
                                                 bool QACouts)
    : TNamed(name, title),
      fMultBinning(false),
      fCentBinning(false),
      fkTBinning(false),
      fmTBinning(false),
      fkTandMultBinning(false),
      fPtQA(false),
      fMassQA(false),
      fMomentumResolution(false),
      fPhiEtaBinning(false),
      fdPhidEtaPlots(false),
      fdPhidEtaPlotsSmallK(false),
      fMixedEventStatistics(true),
      fGetTheControlSampel(false),
      fMode(AliFemtoDreamCollConfig::kNone),
      fMinimalBookingME(false),
      fMinimalBookingSample(false),
      fNumberRadii(0),
      fZVtxBins(),
      fMultBins(),
      fPDGParticleSpecies(),
      fNBinsHists(),
      fMinK_rel(),
      fMaxK_rel(),
      fCentBins(),
      fmTBins(),
      fWhichQAPairs(),
      fClosePairRej(),
      fMixingDepth(0),
      fSpinningDepth(0),
      fCorrelationRange(0.),
      fkTCentrality(false),
      fmTdEtadPhi(false),
      fmTMultBinning(false),
      fEst(AliFemtoDreamEvent::kSPD),
      fAncestors(false),
      fDeltaEtaMax(0.f),
      fDeltaPhiMax(0.f),
      fDoDeltaEtaDeltaPhiCut(false),
      fCoutVariables(QACouts) {
}
AliFemtoDreamCollConfig& AliFemtoDreamCollConfig::operator=(
    const AliFemtoDreamCollConfig& config) {
  if (this != &config) {
    TNamed::operator=(config);
    this->fMultBinning = config.fMultBinning;
    this->fCentBinning = config.fCentBinning;
    this->fkTBinning = config.fkTBinning;
    this->fmTBinning = config.fmTBinning;
    this->fkTandMultBinning = config.fkTandMultBinning;
    this->fPtQA = config.fPtQA;
    this->fMassQA = config.fMassQA;
    this->fMomentumResolution = config.fMomentumResolution;
    this->fPhiEtaBinning = config.fPhiEtaBinning;
    this->fdPhidEtaPlots = config.fdPhidEtaPlots;
    this->fdPhidEtaPlotsSmallK = config.fdPhidEtaPlotsSmallK;
    this->fMixedEventStatistics = config.fMixedEventStatistics;
    this->fGetTheControlSampel = config.fGetTheControlSampel;
    this->fMode = config.fMode;
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
    this->fCorrelationRange = config.fCorrelationRange;
    this->fkTCentrality = config.fkTCentrality;
    this->fmTdEtadPhi = config.fmTdEtadPhi;
    this->fmTMultBinning = config.fmTMultBinning; 
    this->fEst = config.fEst;
    this->fAncestors = config.fAncestors;
    this->fDeltaEtaMax = config.fDeltaEtaMax;
    this->fDeltaPhiMax = config.fDeltaPhiMax;
    this->fDoDeltaEtaDeltaPhiCut = config.fDoDeltaEtaDeltaPhiCut;
    this->fCoutVariables = config.fCoutVariables;
  }
  return *this;
}

AliFemtoDreamCollConfig::~AliFemtoDreamCollConfig() {
}

void AliFemtoDreamCollConfig::SetZBins(std::vector<float> ZBins) {
  //Make sure to set the entries in ascending order!
  //Todo: maybe build in a check for this
  fZVtxBins = ZBins;
}
std::vector<float> AliFemtoDreamCollConfig::GetZVtxBins() {
  //Make sure to set the entries in ascending order!
  if (fCoutVariables) {
    for (auto it : fZVtxBins) {
      std::cout << "Stored z-Vtx bins: " << it << std::endl;
    }
  }
  return fZVtxBins;
}
void AliFemtoDreamCollConfig::SetMultBins(std::vector<int> MultBins) {
  //Make sure to set the entries in ascending order! The last bin to infinite
  //is implicit
  //Todo: maybe build in a check for this
  fMultBins = MultBins;
}
std::vector<int> AliFemtoDreamCollConfig::GetMultBins() {
  if (fCoutVariables) {
    for (auto it : fMultBins) {
      std::cout << "Stored Multiplicity bins: " << it << std::endl;
    }
  }
  return fMultBins;
}
void AliFemtoDreamCollConfig::SetPDGCodes(std::vector<int> PDGCodes) {
  //the order needs to correspond the first particle array in your vector that
  //you hand over in the AliFemtoDreamPartCollection::SetEvent Method!
  fPDGParticleSpecies = PDGCodes;
}
std::vector<int> AliFemtoDreamCollConfig::GetPDGCodes() {
  if (fCoutVariables) {
    for (auto it : fPDGParticleSpecies) {
      std::cout << "Stored PDG Codes for your Analysis " << it << std::endl;
    }
  }
  return fPDGParticleSpecies;
}
int AliFemtoDreamCollConfig::GetNParticleCombinations() {
  //The possible number of combinations for pairing two particles species
  //with itself and all other species is for n species given by:
  //-Combinations within the same species n
  //-Combinations with all other species Binominal(n,2)
  int comb = (int) fPDGParticleSpecies.size();
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
  fNBinsHists = NBins;
}
std::vector<int> AliFemtoDreamCollConfig::GetNBinsHist() {
  if (fCoutVariables) {
    for (auto it : fNBinsHists) {
      std::cout << "Stored  NBins for your Analysis " << it << std::endl;
    }
  }
  return fNBinsHists;
}
void AliFemtoDreamCollConfig::SetMinKRel(std::vector<float> minKRel) {
  //See SetNBinsHist
  fMinK_rel = minKRel;
}
std::vector<float> AliFemtoDreamCollConfig::GetMinKRel() {
  if (fCoutVariables) {
    for (auto it : fMinK_rel) {
      std::cout << "Stored kMin for your Analysis " << it << std::endl;
    }
  }
  return fMinK_rel;
}
void AliFemtoDreamCollConfig::SetMaxKRel(std::vector<float> maxKRel) {
  fMaxK_rel = maxKRel;
}
std::vector<float> AliFemtoDreamCollConfig::GetMaxKRel() {
  if (fCoutVariables) {
    for (auto it : fMaxK_rel) {
      std::cout << "Stored kMax for your Analysis " << it << std::endl;
    }
  }
  return fMaxK_rel;
}
void AliFemtoDreamCollConfig::SetCentBins(std::vector<int> CentBins) {
  fCentBins = CentBins;
}
std::vector<int> AliFemtoDreamCollConfig::GetCentBins() {
  if (fCoutVariables) {
    for (auto it : fCentBins) {
      std::cout << "Stored Centrality Ranges for your Analysis " << it
                << std::endl;
    }
  }
  return fCentBins;
}
void AliFemtoDreamCollConfig::SetmTBins(std::vector<float> mTBins) {
  //Set Bins for the deta dphi and multiplicity mT Binning 
  fmTBins = mTBins;
}
std::vector<float> AliFemtoDreamCollConfig::GetmTBins() {
  if (fCoutVariables) {
    for (auto it : fmTBins) {
      std::cout << "Stored mTbins for your Analysis " << it << std::endl;
    }
  }
  return fmTBins;
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
    fWhichQAPairs.push_back((unsigned int) it);
  }
}

std::vector<unsigned int> AliFemtoDreamCollConfig::GetWhichPairs() {
  if ((int) fWhichQAPairs.size() == 0) {
    AliWarning("===========================================");
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("No Pair QA Specified, setting all to false ");
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("===========================================");
    for (int iQA = 0; iQA < this->GetNParticleCombinations(); ++iQA) {
      fWhichQAPairs.push_back(0);
    }
  } else if ((int) fWhichQAPairs.size() != this->GetNParticleCombinations()) {
    AliFatal("Not all Pairs have a specified QA Behaviour, terminating \n");
  } else {
    for (auto it : fWhichQAPairs) {
      std::cout << "Stored WhichQAPairs for your Analysis " << it << std::endl;
    }
  }
  return fWhichQAPairs;
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
  fClosePairRej = whichPairs;
  //how to select which pairs?
}

std::vector<bool> AliFemtoDreamCollConfig::GetClosePairRej() {
  if ((int) fClosePairRej.size() == 0) {
    AliWarning("=========================================================");
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("!No Close Pair Rejection Specified, setting all to false!");
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("=========================================================");
    for (int iRej = 0; iRej < this->GetNParticleCombinations(); ++iRej) {
      fClosePairRej.push_back(false);
    }
  } else if ((int) fClosePairRej.size() != this->GetNParticleCombinations()) {
    AliFatal("Not all Pairs have a specified QA Behaviour, terminating \n");
  } else {
    int counter = 0;
    std::cout << "=========================================================\n";
    for (auto it : fClosePairRej) {
      std::cout << "Stored CPR for your Analysis "
                << (it ? "active" : "inactive") << std::endl;
      if (it) {
        std::cout
            << "Using the number of pairs set from SetExtendedQAPairs() \n"
            << "to determine the number of combinations for this pair: \n"
            << "Checking " << (unsigned int) fWhichQAPairs.at(counter) / 10
            << " Tracks against "
            << (unsigned int) fWhichQAPairs.at(counter) % 10 << std::endl;
      }
      counter++;
    }
    std::cout << "=========================================================\n";
  }
  return fClosePairRej;
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
