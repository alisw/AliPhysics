///
/// \file AliFemtoUser/AliFemtoKtBinnedCorrFunc.cxx
///

#include "AliFemtoKtBinnedCorrFunc.h"

#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"

#include <TObjArray.h>
#include <cassert>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoKtBinnedCorrFunc);
  /// \endcond
#endif


AliFemtoKtBinnedCorrFunc::AliFemtoKtBinnedCorrFunc(const TString& name, AliFemtoCorrFctn *cf):
  fName(name),
  fPrototypeCF(cf)
{ // no-op
}

UInt_t AliFemtoKtBinnedCorrFunc::AddKtRange(float low, float high)
{
  assert(low < high);

  // 'index' is the location of the correlation function in the buffer
  const UInt_t index = fCFBuffer.size();

  AliFemtoCorrFctn *clone = fPrototypeCF->Clone();
  if (clone == NULL) {
    AliFemtoQinvCorrFctn *qcorr = dynamic_cast<AliFemtoQinvCorrFctn*>(fPrototypeCF);
    if (qcorr != NULL) {
      clone = new AliFemtoQinvCorrFctn(*qcorr);
    } else {
      AliFemtoCorrFctn3DLCMSSym *q3dcorr = dynamic_cast<AliFemtoCorrFctn3DLCMSSym*>(fPrototypeCF);
      if (q3dcorr != NULL) {
        clone = new AliFemtoCorrFctn3DLCMSSym(*q3dcorr);
      } else {
        return NPos;
      }
    }
  }
  fCFBuffer.push_back(clone);

  std::vector<std::pair<Float_t, Float_t> >::iterator it;

  // Find first element which has a higher lower bound than this range.
  // We will insert at this location.
  for (it = fRanges.begin(); it != fRanges.end(); ++it) {
    if (low < it->first) {
      break;
    }
  }

  // get index of where the range will be stored
  const UInt_t mapped_index = std::distance(fRanges.begin(), it);

  // get corresponding iterator for the map
  std::vector<UInt_t>::iterator map_it = fIndexMap.begin() + mapped_index;

  fRanges.insert(it, std::pair<Float_t, Float_t>(low, high));
  fIndexMap.insert(map_it, index);

  return index;
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

void AliFemtoKtBinnedCorrFunc::AddRealPair(AliFemtoPair *pair)
{
  UInt_t index = FindKtBin(pair);
  if (index == NPos) {
    return;
  }
  fCFBuffer[index]->AddRealPair(pair);
}

void AliFemtoKtBinnedCorrFunc::AddMixedPair(AliFemtoPair *pair)
{
  UInt_t index = FindKtBin(pair);
  if (index == NPos) {
    return;
  }
  fCFBuffer[index]->AddMixedPair(pair);
}

TList* AliFemtoKtBinnedCorrFunc::GetOutputList()
{
  TList *olist = new TList();

  TObjArray *output = new TObjArray();
  output->SetName(fName);

  for (UInt_t i = 0; i < fRanges.size(); ++i) {
    const std::pair<Float_t, Float_t> &range = fRanges[i];
    AliFemtoCorrFctn *cf = fCFBuffer[i];

    const TString range_name = Form("%0.3f_%0.3f", range.first, range.second);
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
