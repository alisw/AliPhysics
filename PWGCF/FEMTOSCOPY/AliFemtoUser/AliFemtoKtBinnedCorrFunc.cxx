///
/// \file AliFemtoUser/AliFemtoKtBinnedCorrFunc.cxx
///

#include "AliFemtoKtBinnedCorrFunc.h"

#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"

#include <TObjArray.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <tuple>



AliFemtoKtBinnedCorrFunc::AliFemtoKtBinnedCorrFunc(const TString& name, AliFemtoCorrFctn *cf):
  fName(name),
  fPrototypeCF(cf)
{ // no-op
}

AliFemtoKtBinnedCorrFunc::~AliFemtoKtBinnedCorrFunc()
{
  delete fPrototypeCF;
}

UInt_t AliFemtoKtBinnedCorrFunc::AddKtRange(float low, float high)
{
  assert(low < high);

  AliFemtoCorrFctn *clone = fPrototypeCF->Clone();
  if (!clone) {
    return NPos;
  }

  const std::pair<Float_t, Float_t> range(low, high);

  const auto insert_location = std::lower_bound(fRanges.begin(), fRanges.end(), range);
  fRanges.insert(insert_location, range);

  auto offset = std::distance(fRanges.begin(), insert_location);
  const auto buffer_location = std::next(fCFBuffer.begin(), offset);

  fCFBuffer.insert(buffer_location, clone);

  return static_cast<Int_t>(offset);
}

UInt_t AliFemtoKtBinnedCorrFunc::AddKtRange(const TString& name, float low, float high)
{
  UInt_t index = AddKtRange(low, high);
  // fCFBuffer[index]->GetName() += name;
  return index;
}

void AliFemtoKtBinnedCorrFunc::AddKtRanges(const float data[], float stop_at)
{
  UInt_t i = 1;
  do {
    AddKtRange(data[i-1], data[i]);
  } while (data[++i] != stop_at);
}

inline
std::tuple<Int_t, Int_t>
load_integers(
  Float_t kt,
  const std::vector<std::pair<Float_t, Float_t>> &ranges)
{
  auto lower = std::find_if(ranges.begin(), ranges.end(),
                            [=](std::pair<Float_t, Float_t> r) {
                              return r.first <= kt && kt < r.second; });

  auto higher = std::find_if(lower, ranges.end(),
                             [=](std::pair<Float_t, Float_t> r) {
                               return kt < r.first; });

  return std::make_tuple(std::distance(ranges.begin(), lower),
                         std::distance(ranges.begin(), higher));
}

void AliFemtoKtBinnedCorrFunc::AddRealPair(AliFemtoPair *pair)
{
  const Float_t kt = pair->KT();

  Int_t idx, stop;
  std::tie(idx, stop) = load_integers(kt, fRanges);

  for (; idx < stop; ++idx) {
    const auto r = fRanges[idx];
    if (r.first <= kt && kt < r.second) {
      fCFBuffer[idx]->AddRealPair(pair);
    }
  }
}

void AliFemtoKtBinnedCorrFunc::AddMixedPair(AliFemtoPair *pair)
{
  const Float_t kt = pair->KT();

  Int_t idx, stop;
  std::tie(idx, stop) = load_integers(kt, fRanges);

  for (; idx < stop; ++idx) {
    const auto r = fRanges[idx];
    if (r.first <= kt && kt < r.second) {
      fCFBuffer[idx]->AddMixedPair(pair);
    }
  }
}

TList* AliFemtoKtBinnedCorrFunc::GetOutputList()
{
  TList *olist = new TList();

  TObjArray *output = new TObjArray();
  output->SetName(fName);

  for (UInt_t i = 0; i < fRanges.size(); ++i) {
    const std::pair<Float_t, Float_t> &range = fRanges[i];
    AliFemtoCorrFctn *cf = fCFBuffer[i];

    auto stringify = [] (const std::pair<float, float> pair) {
      auto stringify_float = [] (float n) {
        const auto res = (std::floor(n) == n)
                       ? Form("%0.1f", n)
                       : Form("%g", n);
        return TString(res);
      };

      TString lo = stringify_float(pair.first),
              hi = stringify_float(pair.second);

      return lo + "_" + hi;
    };

    const TString range_name = stringify(range);
    TObjArray *cf_output = new TObjArray();
    cf_output->SetName(range_name);

    TList *cf_list = cf->GetOutputList();
    cf_output->AddAll(cf_list);
    delete cf_list;

    output->Add(cf_output);
  }

  olist->Add(output);
  return olist;
}
