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

  const auto r_insert_loc = std::lower_bound(fRanges.begin(), fRanges.end(), range);
  const auto loc = fRanges.insert(r_insert_loc, range);

  const auto offset = std::distance(fRanges.begin(), loc);
  const auto cf_insert_loc = std::next(fCFBuffer.begin(), offset);

  // std::cout << "Adding "
  //           << Form("{%g, %g}", low, high)
  //           << " to location: " << offset << "/" << fCFBuffer.size()
  //           << "\n";

  fCFBuffer.insert(cf_insert_loc, clone);

  return offset;
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

void AliFemtoKtBinnedCorrFunc::AddPair(AliFemtoPair *pair, bool is_same_event)
{
  auto add_pair_to = [is_same_event, pair] (AliFemtoCorrFctn *cf)
    {
      (is_same_event) ? cf->AddRealPair(pair)
                      : cf->AddMixedPair(pair);
    };

  // auto add_pair_to = mixed
  //                  ? [pair] (AliFemtoCorrFctn *cf) { cf->AddMixedPair(pair); }
  //                  : [pair] (AliFemtoCorrFctn *cf) { cf->AddRealPair(pair); };

  const Float_t kt = pair->KT();

  // auto lower = std::lower_bound(fRanges.begin(), fRanges.end(), kt,
  //                               [](std::pair<Float_t, Float_t> r, float kt) {
  //                                 return r.first <= kt; });
  // std::size_t idx = std::distance(fRanges.begin(), lower);

  std::size_t idx = 0;
  for (; idx < fRanges.size(); ++idx) {
    const auto r = fRanges[idx];
    if (r.first <= kt && kt < r.second) {
      break;
    }
  }

  for (; idx < fRanges.size(); ++idx) {
    const auto r = fRanges[idx];

    // no remaining ranges
    if (kt < r.first) {
      break;
    }

    // we are within this range
    if (kt < r.second) {
      add_pair_to(fCFBuffer[idx]);
    }
  }
}

void AliFemtoKtBinnedCorrFunc::AddRealPair(AliFemtoPair *pair)
{
  AddPair(pair, true);
}

void AliFemtoKtBinnedCorrFunc::AddMixedPair(AliFemtoPair *pair)
{
  AddPair(pair, false);
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
