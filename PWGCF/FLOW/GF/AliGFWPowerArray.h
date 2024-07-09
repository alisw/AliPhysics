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

#ifndef PWGCF_FLOW_GF_AliGFWPOWERARRAY_H_
#define PWGCF_FLOW_GF_AliGFWPOWERARRAY_H_

#include <vector>
#include <cmath>
#include <string>

typedef std::vector<int> HarSet;
class AliGFWPowerArray
{
 public:
  static HarSet GetPowerArray(std::vector<HarSet> inHarmonics);
  static void PowerArrayTest();

 private:
  static int getHighestHarmonic(const HarSet& inhar);
  static HarSet TrimVec(HarSet hars, int ind);
  static HarSet AddConstant(HarSet hars, int offset);
  static void FlushVectorToMaster(HarSet& masterVector, HarSet& comVec, const int& MaxPower);
  static void RecursiveFunction(HarSet& masterVector, HarSet hars, int offset, const int& MaxPower);
  static void PrintVector(const HarSet& singleSet);
};
#endif // PWGCF_FLOW_GF_AliGFWPOWERARRAY_H_
