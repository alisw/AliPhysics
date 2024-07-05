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

#include "AliGFW.h"

using std::complex;
using std::pair;
using std::string;
using std::vector;

AliGFW::AliGFW() : fInitialized(false) {}

AliGFW::~AliGFW()
{
  for (auto pItr = fCumulants.begin(); pItr != fCumulants.end(); ++pItr)
    pItr->DestroyComplexVectorArray();
};
void AliGFW::AddRegion(string refName, double lEtaMin, double lEtaMax, int lNpT, int BitMask)
{
  if (lNpT < 1) {
    printf("Number of pT bins cannot be less than 1! Not adding anything.\n");
    return;
  }
  if (lEtaMin >= lEtaMax) {
    printf("Eta min. cannot be more than eta max! Not adding...\n");
    return;
  }
  if (refName == "") {
    printf("Region must have a name!\n");
    return;
  }
  Region lOneRegion;
  lOneRegion.Nhar = 0;            // Empty for now
  lOneRegion.powsDefined = false; // If vector with powers defined, set this to zero
  lOneRegion.NparVec = {};        // Empty for now
  lOneRegion.EtaMin = lEtaMin;    // Min. eta
  lOneRegion.EtaMax = lEtaMax;    // Max. eta
  lOneRegion.NpT = lNpT;          // Number of pT bins
  lOneRegion.rName = refName;     // Name of the region
  lOneRegion.BitMask = BitMask;   // Bit mask
  AddRegion(lOneRegion);
};
void AliGFW::AddRegion(string refName, vector<int> lNparVec, double lEtaMin, double lEtaMax, int lNpT, int BitMask)
{
  AddRegion(refName, lEtaMin, lEtaMax, lNpT, BitMask);
  (fRegions.end() - 1)->Nhar = static_cast<int>(lNparVec.size());
  (fRegions.end() - 1)->NparVec = lNparVec;
  (fRegions.end() - 1)->powsDefined = true;
};
void AliGFW::AddRegion(string refName, int lNhar, int lNpar, double lEtaMin, double lEtaMax, int lNpT, int BitMask)
{
  vector<int> tVec = {};
  for (int i = 0; i < lNhar; i++)
    tVec.push_back(lNpar);
  AddRegion(refName, tVec, lEtaMin, lEtaMax, lNpT, BitMask);
};
void AliGFW::AddRegion(string refName, int lNhar, int* lNparVec, double lEtaMin, double lEtaMax, int lNpT, int BitMask)
{
  vector<int> tVec = {};
  for (int i = 0; i < lNhar; i++)
    tVec.push_back(lNparVec[i]);
  AddRegion(refName, tVec, lEtaMin, lEtaMax, lNpT, BitMask);
};
int AliGFW::CreateRegions()
{
  for (auto pItr = fCumulants.begin(); pItr != fCumulants.end(); ++pItr)
    pItr->DestroyComplexVectorArray();
  fCumulants.clear();
  InitializePowerArrays();
  if (fRegions.size() < 1) {
    printf("No regions set. Skipping...\n");
    return 0;
  }
  int nRegions = 0;
  for (auto pItr = fRegions.begin(); pItr != fRegions.end(); pItr++) {
    AliGFWCumulant* lCumulant = new AliGFWCumulant();
    lCumulant->CreateComplexVectorArrayVarPower(pItr->Nhar, pItr->NparVec, pItr->NpT);
    fCumulants.push_back(*lCumulant);
    ++nRegions;
  }
  if (nRegions)
    fInitialized = true;
  return nRegions;
};
void AliGFW::Fill(double eta, int ptin, double phi, double weight, int mask, double SecondWeight)
{
  // if(!fInitialized) return;
  for (int i = 0; i < static_cast<int>(fRegions.size()); ++i) {
    if (fRegions.at(i).EtaMin < eta && fRegions.at(i).EtaMax > eta && (fRegions.at(i).BitMask & mask))
      fCumulants.at(i).FillArray(ptin, phi, weight, SecondWeight);
  }
};
complex<double> AliGFW::TwoRec(int n1, int n2, int p1, int p2, int ptbin, AliGFWCumulant* r1, AliGFWCumulant* r2, AliGFWCumulant* r3)
{
  complex<double> part1 = r1->Vec(n1, p1, ptbin);
  complex<double> part2 = r2->Vec(n2, p2, ptbin);
  complex<double> part3 = r3 ? r3->Vec(n1 + n2, p1 + p2, ptbin) : complex<double>(0., 0.);
  complex<double> formula = part1 * part2 - part3;
  return formula;
};
complex<double> AliGFW::RecursiveCorr(AliGFWCumulant* qpoi, AliGFWCumulant* qref, AliGFWCumulant* qol, int ptbin, vector<int>& hars)
{
  vector<int> pows;
  for (int i = 0; i < static_cast<int>(hars.size()); i++)
    pows.push_back(1);
  return RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows);
};

complex<double> AliGFW::RecursiveCorr(AliGFWCumulant* qpoi, AliGFWCumulant* qref, AliGFWCumulant* qol, int ptbin, vector<int>& hars, vector<int>& pows)
{
  if ((pows.at(0) != 1) && qol)
    qpoi = qol; // if the power of POI is not unity, then always use overlap (if defined).
  // Only valid for 1 particle of interest though!
  if (hars.size() < 2)
    return qpoi->Vec(hars.at(0), pows.at(0), ptbin);
  if (hars.size() < 3)
    return TwoRec(hars.at(0), hars.at(1), pows.at(0), pows.at(1), ptbin, qpoi, qref, qol);
  int harlast = hars.at(hars.size() - 1);
  int powlast = pows.at(pows.size() - 1);
  hars.erase(hars.end() - 1);
  pows.erase(pows.end() - 1);
  complex<double> formula = RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows) * qref->Vec(harlast, powlast);
  int lDegeneracy = 1;
  int harSize = static_cast<int>(hars.size());
  for (int i = harSize - 1; i >= 0; i--) {
    // checking if current configuration is a permutation of the next one.
    // Need to have more than 2 harmonics though, otherwise it doesn't make sense.
    if (i > 2) {                                                          // only makes sense when we have more than two harmonics remaining
      if (hars.at(i) == hars.at(i - 1) && pows.at(i) == pows.at(i - 1)) { // if it is a permutation, then increase degeneracy and continue;
        lDegeneracy++;
        continue;
      }
    }
    hars.at(i) += harlast;
    pows.at(i) += powlast;
    complex<double> subtractVal = RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows);
    if (lDegeneracy > 1) {
      subtractVal *= lDegeneracy;
      lDegeneracy = 1;
    }
    formula -= subtractVal;
    hars.at(i) -= harlast;
    pows.at(i) -= powlast;
  }
  hars.push_back(harlast);
  pows.push_back(powlast);
  return formula;
};
void AliGFW::Clear()
{
  if (!fInitialized)
    CreateRegions();
  for (auto ptr = fCumulants.begin(); ptr != fCumulants.end(); ++ptr)
    ptr->ResetQs();
};
AliGFW::CorrConfig AliGFW::GetCorrelatorConfig(string config, string head, bool ptdif)
{
  // First remove all ; and ,:
  s_replace_all(config, ",", " ");
  s_replace_all(config, ";", " ");
  s_replace_all(config, "| ", "|");
  // If pT-bin is provided, then look for & remove space before "(" (so that it's clean afterwards)
  while (s_index(config, " (") > -1)
    s_replace_all(config, " (", "(");
  // Then make sure we don't have any double-spaces:
  while (s_index(config, "  ") > -1)
    s_replace_all(config, "  ", " ");
  vector<int> regs;
  vector<int> hars;
  int sz1 = 0;
  int szend = 0;
  string ts, ts2;
  CorrConfig ReturnConfig;
  // Fetch region descriptor
  if (!s_tokenize(config, ts, szend, "{")) {
    printf("Could not find any harmonics!\n");
    return ReturnConfig;
  }
  szend = 0;
  int counter = 0;
  while (s_tokenize(config, ts, szend, "{")) {
    counter++;
    ReturnConfig.Regs.push_back(vector<int>{});
    ReturnConfig.Hars.push_back(vector<int>{});
    ReturnConfig.Overlap.push_back(-1); // initially, assume no overlap
    // Check if there's a particular pT bin I should be using here. If so, store it (otherwise, it's bin 0)
    int ptbin = -1;
    int sz2 = 0;
    if (s_contains(ts, "(")) {
      if (!s_contains(ts, ")")) {
        printf("Missing \")\" in the configurator. Returning...\n");
        return ReturnConfig;
      }
      sz2 = s_index(ts, "(");
      sz1 = sz2 + 1;
      s_tokenize(ts, ts2, sz1, ")");
      ptbin = stoi(ts2);
      ts.erase(sz2, (sz1 - sz2 + 1));
      szend -= (sz1 - sz2); // szend also becomes shorter
      // also need to remove this from config now:
      sz2 = s_index(config, "(");
      sz1 = s_index(config, ")");
      config.erase(sz2, sz1 - sz2 + 1);
    }
    ReturnConfig.ptInd.push_back(ptbin);
    sz1 = 0;
    // Fetch regions
    while (s_tokenize(ts, ts2, sz1, " ")) {
      if (sz1 >= szend)
        break;
      bool isOverlap = s_contains(ts2, "|");
      if (isOverlap)
        ts2.erase(0, 1); // If overlap, remove the delimiter |
      int ind = FindRegionByName(ts2);
      if (ts2 == " " || ts2 == "")
        continue;
      if (ind < 0) {
        printf("Could not find region named %s!\n", ts2.c_str());
        break;
      }
      if (!isOverlap)
        ReturnConfig.Regs.at(counter - 1).push_back(ind);
      else
        ReturnConfig.Overlap.at(static_cast<int>(ReturnConfig.Overlap.size()) - 1) = ind;
    }
    string harstr;
    s_tokenize(config, harstr, szend, "}");
    int dummys = 0;
    // Fetch harmonics
    while (s_tokenize(harstr, ts, dummys, " "))
      ReturnConfig.Hars.at(counter - 1).push_back(stoi(ts));
  }
  ReturnConfig.Head = head;
  ReturnConfig.pTDif = ptdif;
  // ReturnConfig.pTbin = ptbin;
  fListOfCFGs.push_back(ReturnConfig);
  return ReturnConfig;
};

complex<double> AliGFW::Calculate(int poi, int ref, vector<int> hars, int ptbin)
{
  AliGFWCumulant* qref = &fCumulants.at(ref);
  AliGFWCumulant* qpoi = &fCumulants.at(poi);
  AliGFWCumulant* qovl = qpoi;
  return RecursiveCorr(qpoi, qref, qovl, ptbin, hars);
};
complex<double> AliGFW::Calculate(CorrConfig corconf, int ptbin, bool SetHarmsToZero)
{
  // if(!fInitialized) return complex<double>(0,0); //First check if initialised, if not -- initialize, and if it fails, return
  if (corconf.Regs.size() == 0)
    return complex<double>(0, 0); // Check if we have any regions at all
  complex<double> retval(1, 0);
  int ptInd;
  for (int i = 0; i < static_cast<int>(corconf.Regs.size()); i++) { // looping over all regions
    if (corconf.Regs.at(i).size() == 0)
      return complex<double>(0, 0); // again, if no regions in the current subevent, then quit immediatelly
    ptInd = corconf.ptInd.at(i);    // for i=0 (potentially, POI)
    if (ptInd < 0)
      ptInd = ptbin;
    // picking up the indecies of regions...
    int poi = corconf.Regs.at(i).at(0);
    int ref = (corconf.Regs.at(i).size() > 1) ? corconf.Regs.at(i).at(1) : corconf.Regs.at(i).at(0);
    int ovl = corconf.Overlap.at(i);
    // and regions themselves
    AliGFWCumulant* qref = &fCumulants.at(ref);
    AliGFWCumulant* qpoi = &fCumulants.at(poi);
    if (!qref->IsPtBinFilled(ptInd))
      return complex<double>(0, 0); // if REF is not filled, don't even continue. Could be redundant, but should save little CPU time
    if (!qpoi->IsPtBinFilled(ptInd))
      return complex<double>(0, 0); // if POI is not filled, don't even continue. Could be redundant, but should save little CPU time
    AliGFWCumulant* qovl = 0;
    // Check if in the ref. region we have enough particles (no. of particles in the region >= no of harmonics for subevent)
    int sz1 = corconf.Hars.at(i).size();
    if (poi != ref)
      sz1--;
    if (qref->GetN() < sz1)
      return complex<double>(0, 0);
    // Then, figure the overlap
    if (ovl > -1) // if overlap is defined, then (unless it's explicitly disabled)
      qovl = &fCumulants.at(ovl);
    else if (ref == poi)
      qovl = qref; // If ref and poi are the same, then the same is for overlap. Only, when OL not explicitly defined
    if (SetHarmsToZero) {
      for (int j = 0; j < static_cast<int>(corconf.Hars.at(i).size()); j++) {
        corconf.Hars.at(i).at(j) = 0;
      }
    }
    retval *= RecursiveCorr(qpoi, qref, qovl, ptInd, corconf.Hars.at(i));
  }
  return retval;
};
vector<pair<int, vector<int>>> AliGFW::GetHarmonicsSingleConfig(const CorrConfig& incfg)
{
  vector<pair<int, vector<int>>> retPair;
  for (int iR = 0; iR < static_cast<int>(incfg.Regs.size()); iR++) {
    if (static_cast<int>(incfg.Regs[iR].size()) > 1) {
      retPair.push_back(make_pair(incfg.Regs[iR][0], vector<int>{incfg.Hars[iR][0]})); // If we have a PoI, then it comes with the first harmonic
      retPair.push_back(make_pair(incfg.Regs[iR][1], incfg.Hars[iR]));                 // Then the second is ref. with full harmonics
    } else {
      retPair.push_back(make_pair(incfg.Regs[iR][0], incfg.Hars[iR])); // Otherwise, it's only ref with all harmonics
    }
    if (incfg.Overlap[iR] > -1)
      retPair.push_back(make_pair(incfg.Overlap[iR], incfg.Hars[iR])); // if overlap provided, then also fetch its harmonics
  }
  return retPair;
};
void AliGFW::InitializePowerArrays()
{
  vector<vector<vector<int>>> harSets(static_cast<int>(fRegions.size()));
  for (const CorrConfig& lConf : fListOfCFGs) {
    auto HarPerReg = GetHarmonicsSingleConfig(lConf);
    for (auto oneHar : HarPerReg)
      harSets[oneHar.first].push_back(oneHar.second);
  }
  // Now, loop through all combinations of different harmonics for each region and calculate power arrays
  for (int i = 0; i < static_cast<int>(harSets.size()); i++) {
    if (fRegions[i].powsDefined)
      continue; // Only do if powers have not been externally defined
    vector<int> powerArray = AliGFWPowerArray::GetPowerArray(harSets[i]);
    fRegions[i].Nhar = static_cast<int>(powerArray.size());
    fRegions[i].NparVec = powerArray;
    fRegions[i].powsDefined = true;
  }
};
complex<double> AliGFW::Calculate(int poi, vector<int> hars)
{
  AliGFWCumulant* qpoi = &fCumulants.at(poi);
  return RecursiveCorr(qpoi, qpoi, qpoi, 0, hars);
};
int AliGFW::FindRegionByName(string refName)
{
  for (int i = 0; i < static_cast<int>(fRegions.size()); i++)
    if (fRegions.at(i).rName == refName)
      return i;
  return -1;
};
// String processing:
int AliGFW::s_index(string& instr, const string& pattern, const int& spos)
{
  return instr.find(pattern, spos);
};
bool AliGFW::s_contains(string& instr, const string& pattern)
{
  return (s_index(instr, pattern) > -1);
};
void AliGFW::s_replace(string& instr, const string& pattern1, const string& pattern2, const int& spos)
{
  int lpos = s_index(instr, pattern1, spos);
  if (lpos < 0)
    return;
  instr.replace(lpos, pattern1.size(), pattern2);
};
void AliGFW::s_replace_all(string& instr, const string& pattern1, const string& pattern2)
{
  int lpos = s_index(instr, pattern1);
  while (lpos > -1) {
    s_replace(instr, pattern1, pattern2, lpos);
    lpos = s_index(instr, pattern1, lpos);
  }
};
bool AliGFW::s_tokenize(string& instr, string& subs, int& spos, const string& delim)
{
  if (spos < 0 || spos >= static_cast<int>(instr.size())) {
    spos = -1;
    subs = "";
    return false;
  }
  int lpos = s_index(instr, delim, spos);
  if (lpos < 0)
    lpos = instr.size();
  subs = instr.substr(spos, lpos - spos);
  spos = lpos + 1;
  return true;
}
