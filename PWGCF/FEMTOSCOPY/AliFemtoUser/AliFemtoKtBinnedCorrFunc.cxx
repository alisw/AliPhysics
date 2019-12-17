///
/// \file AliFemtoUser/AliFemtoKtBinnedCorrFunc.cxx
///

#include "AliFemtoKtBinnedCorrFunc.h"

#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"

#include <TObjArray.h>
#include <TH1D.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <tuple>


AliFemtoKtBinnedCorrFunc::AliFemtoKtBinnedCorrFunc(const TString& name, AliFemtoCorrFctn *cf)
  : AliFemtoCorrFctn()
  , fName(name)
  , fPrototypeCF(cf)
  , fCFBuffer()
  , fRanges()
  , fKtMonitor(nullptr)
{
}

AliFemtoKtBinnedCorrFunc::~AliFemtoKtBinnedCorrFunc()
{
  delete fPrototypeCF;
  delete fKtMonitor;
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

  const Float_t kt = pair->KT();

  std::size_t idx = 0;
  for (; idx < fRanges.size(); ++idx) {
    const auto r = fRanges[idx];
    if (r.first <= kt && kt < r.second) {
      break;
    }
  }

  bool ktmonitor_needs_update = false;
  for (; idx < fRanges.size(); ++idx) {
    const auto r = fRanges[idx];

    // no remaining ranges
    if (kt < r.first) {
      break;
    }

    // we are within this range
    if (kt < r.second) {
      add_pair_to(fCFBuffer[idx]);
      ktmonitor_needs_update = true;
    }
  }

  // update kt-monitor for real events
  if (is_same_event and ktmonitor_needs_update) {
    fKtMonitor->Fill(kt);
  }
}

inline
TH1D* build_kt_monitor(const std::vector<std::pair<float, float>> ranges)
{
  // determinte high and low points
  auto low = ranges.front().first,
       high = -1.0f;

  for (auto r : ranges) {
    high = std::max(high, r.second);
  }

  UInt_t max_bins = 12 * 15 * 7;  // <- easy rebinning, not too big
  UInt_t small_bins = std::lrint((high - low) / 0.005);  // 5MeV bins

  // determine best limit
  auto nbins = std::min(small_bins, max_bins);

  auto hist = new TH1D("kTDist", "k_{T} Distribution", nbins, low, high);
  return hist;
}

void AliFemtoKtBinnedCorrFunc::AddRealPair(AliFemtoPair *pair)
{
  if (__builtin_expect(fKtMonitor == nullptr, 0)) {
    fKtMonitor = build_kt_monitor(fRanges);
  }
  AddPair(pair, true);
}

void AliFemtoKtBinnedCorrFunc::AddMixedPair(AliFemtoPair *pair)
{
  AddPair(pair, false);
}

TList* AliFemtoKtBinnedCorrFunc::GetOutputList()
{
  TList *olist = new TList();
  AddOutputObjectsTo(*olist);
  return olist;
}

void AliFemtoKtBinnedCorrFunc::AddOutputObjectsTo(TCollection &olist)
{
  TObjArray *output = new TObjArray();
  output->SetName(fName);

  if (fKtMonitor == nullptr) {
    fKtMonitor = build_kt_monitor(fRanges);
  }

  output->Add(fKtMonitor);

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

    cf->AddOutputObjectsTo(*cf_output);

    output->Add(cf_output);
  }

  olist.Add(output);
}

void
AliFemtoKtBinnedCorrFunc::EventBegin(const AliFemtoEvent* ev)
{
  for (auto *cf : fCFBuffer) {
    cf->EventBegin(ev);
  }
}

void
AliFemtoKtBinnedCorrFunc::EventEnd(const AliFemtoEvent* ev)
{
  for (auto *cf : fCFBuffer) {
    cf->EventEnd(ev);
  }
}
