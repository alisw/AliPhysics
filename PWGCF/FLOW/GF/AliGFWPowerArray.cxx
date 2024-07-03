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

#include "AliGFWPowerArray.h"

using std::string;
using std::vector;

int AliGFWPowerArray::getHighestHarmonic(const HarSet& inhar)
{
  // Highest possible harmonic: sum of same-sign harmonics
  int maxPos = 0, maxNeg = 0;
  for (int val : inhar)
    if (val > 0)
      maxPos += val;
    else
      maxNeg += abs(val);
  return maxPos > maxNeg ? maxPos : maxNeg;
};
HarSet AliGFWPowerArray::TrimVec(HarSet hars, int ind)
{
  HarSet retVec = hars;
  retVec.erase(retVec.begin() + ind);
  return retVec;
};
HarSet AliGFWPowerArray::AddConstant(HarSet hars, int offset)
{
  HarSet retVec = hars;
  for (int& val : retVec)
    val += offset;
  return retVec;
};
void AliGFWPowerArray::FlushVectorToMaster(HarSet& masterVector, HarSet& comVec, const int& MaxPower)
{
  int nPartLoc = MaxPower - comVec.size() + 1;
  for (auto& val : comVec) {
    int absVal = abs(val);
    if (masterVector.at(absVal) < nPartLoc) {
      masterVector.at(absVal) = nPartLoc;
    }
  }
};
void AliGFWPowerArray::RecursiveFunction(HarSet& masterVector, HarSet hars, int offset, const int& MaxPower)
{
  HarSet compVec = AddConstant(hars, offset);
  FlushVectorToMaster(masterVector, compVec, MaxPower);
  for (int i = 0; i < hars.size(); i++)
    RecursiveFunction(masterVector, TrimVec(hars, i), offset + hars.at(i), MaxPower);
  ;
};
void AliGFWPowerArray::PrintVector(const HarSet& singleSet)
{
  int vcSize = static_cast<int>(singleSet.size());
  if (!vcSize)
    printf("Vector is empty!\n");
  printf("{%i", singleSet[0]);
  for (int i = 1; i < vcSize; i++)
    printf(", %i", singleSet[i]);
  printf("}\n");
}
HarSet AliGFWPowerArray::GetPowerArray(vector<HarSet> inHarmonics)
{
  // First, find maximum number of particle correlations ( = max power) and maximum (sum of) harmonics
  int MaxHar = 0;
  int nMaxPart = 0;
  for (HarSet singleSet : inHarmonics) {
    int harSum = getHighestHarmonic(singleSet);
    MaxHar = harSum > MaxHar ? harSum : MaxHar;
  }
  // Make a vector with MaxHar+1 entries (entry 0 for sum=0)
  HarSet retVec = HarSet(MaxHar + 1);
  // Then loop over all combinations and calculate max powers
  for (HarSet singleSet : inHarmonics) {
    int lNPart = static_cast<int>(singleSet.size()); // Total number of particles correlated
    RecursiveFunction(retVec, singleSet, 0, lNPart);
    // Harmonic sum = 0 is a special case. In principle all 0 cases with non-zero harmonics are captured by the function above, but to calculate normalization, we set all harmonics to 0. This means that sum=0 power is the max number of harmonics/particles being correlated
    nMaxPart = (lNPart > nMaxPart) ? lNPart : nMaxPart;
  }
  // Override the sum=0 power with the number of correlated particles
  if (retVec[0] < nMaxPart)
    retVec[0] = nMaxPart;
  // Need an extra power ( = 0) for all non-zero powers
  for (int& val : retVec)
    if (val != 0)
      val++;
  return retVec;
};
void AliGFWPowerArray::PowerArrayTest()
{
  vector<HarSet> AllHars = {
    HarSet{2},
    HarSet{3},
    HarSet{2, 2},
    HarSet{3, 3}};
  printf("Input harmonics are:\n");
  for (HarSet inSet : AllHars)
    PrintVector(inSet);
  printf("The configuration of powers must then be:\n");
  auto vc = GetPowerArray(AllHars);
  PrintVector(vc);
};
