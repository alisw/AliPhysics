// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.)
Class steers the initialization and calculation of n-particle correlations. Uses recursive function, all terms are calculated only once.
Latest version includes the calculation of any number of gaps and any combination of harmonics (including eg symmetric cumulants, etc.)
If used, modified, or distributed, please aknowledge the author of this code.
*/
#ifndef PWGCF_FLOW_GF_AliGFW_H_
#define PWGCF_FLOW_GF_AliGFW_H_

#include "AliGFWCumulant.h"
#include "AliGFWPowerArray.h"
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <complex>

class AliGFW
{
 public:
  struct Region {
    int Nhar, NpT;
    std::vector<int> NparVec{};
    double EtaMin = -999;
    double EtaMax = -999;
    int BitMask = 1;
    std::string rName = "";
    bool powsDefined = false;
    bool operator<(const Region& a) const
    {
      return EtaMin < a.EtaMin;
    };
    void PrintStructure() { printf("%s: eta [%f.. %f].", rName.c_str(), EtaMin, EtaMax); }
  };
  struct CorrConfig {
    std::vector<std::vector<int>> Regs{};
    std::vector<std::vector<int>> Hars{};
    std::vector<int> Overlap;
    std::vector<int> ptInd;
    bool pTDif = false;
    std::string Head = "";
  };
  AliGFW();
  ~AliGFW();
  std::vector<Region> fRegions;
  std::vector<AliGFWCumulant> fCumulants;
  void AddRegion(std::string refName, double lEtaMin, double lEtaMax, int lNpT, int BitMask);
  void AddRegion(std::string refName, std::vector<int> lNparVec, double lEtaMin, double lEtaMax, int lNpT, int BitMask); // Legacy
  void AddRegion(std::string refName, int lNhar, int lNpar, double lEtaMin, double lEtaMax, int lNpT, int BitMask);      // Legacy support, all powers are the same
  void AddRegion(std::string refName, int lNhar, int* lNparVec, double lEtaMin, double lEtaMax, int lNpT, int BitMask);  // Legacy support, array instead of a vector
  int CreateRegions();
  void Fill(double eta, int ptin, double phi, double weight, int mask, double secondWeight = -1);
  void Clear();
  AliGFWCumulant GetCumulant(int index) { return fCumulants.at(index); }
  CorrConfig GetCorrelatorConfig(std::string config, std::string head = "", bool ptdif = false);
  std::complex<double> Calculate(CorrConfig corconf, int ptbin, bool SetHarmsToZero);
  void InitializePowerArrays();

 protected:
  bool fInitialized;
  std::vector<CorrConfig> fListOfCFGs;
  std::complex<double> TwoRec(int n1, int n2, int p1, int p2, int ptbin, AliGFWCumulant*, AliGFWCumulant*, AliGFWCumulant*);
  std::complex<double> RecursiveCorr(AliGFWCumulant* qpoi, AliGFWCumulant* qref, AliGFWCumulant* qol, int ptbin, std::vector<int>& hars, std::vector<int>& pows); // POI, Ref. flow, overlapping region
  std::complex<double> RecursiveCorr(AliGFWCumulant* qpoi, AliGFWCumulant* qref, AliGFWCumulant* qol, int ptbin, std::vector<int>& hars);                         // POI, Ref. flow, overlapping region
  void AddRegion(Region inreg) { fRegions.push_back(inreg); }
  Region GetRegion(int index) { return fRegions.at(index); }
  int FindRegionByName(std::string refName);
  std::vector<std::pair<int, std::vector<int>>> GetHarmonicsSingleConfig(const CorrConfig&);
  // Calculating functions:
  std::complex<double> Calculate(int poi, int ref, std::vector<int> hars, int ptbin = 0); // For differential, need POI and reference
  std::complex<double> Calculate(int poi, std::vector<int> hars);                         // For integrated case
  // Operations on strings. Equivalent to TString operations, but one to rid of root dependence
  int s_index(std::string& instr, const std::string& pattern, const int& spos = 0);
  bool s_contains(std::string& instr, const std::string& pattern);
  void s_replace(std::string& instr, const std::string& pattern1, const std::string& pattern2, const int& spos = 0);
  void s_replace_all(std::string& instr, const std::string& pattern1, const std::string& pattern2);
  bool s_tokenize(std::string& instr, std::string& substr, int& spos, const std::string& delim);
};
#endif // PWGCF_FLOW_GF_AliGFW_H_
