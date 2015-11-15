///
/// \file AliFemtoModelCorrFctnKStar.cxx
/// \author Andrew Kubera
///

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelCorrFctnKStar);
  /// \endcond
#endif

#include "AliFemtoModelCorrFctnKStar.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelAllHiddenInfo.h"
#include "AliFemtoPair.h"

#include <TH1F.h>
#include <TH2F.h>


AliFemtoModelCorrFctnKStar::AliFemtoModelCorrFctnKStar():
  AliFemtoModelCorrFctn()
  , fResNum(NULL)
  , fResDen(NULL)
  , fTrueNum(NULL)
  , fTrueDen(NULL)
{
  fResNum = new TH1F("kstar_num", "KStar - Numerator; k*(GeV);", 200, 0.0, 1.0);
  fResNum->Sumw2();

  fResDen = new TH1F("kstar_num", "KStar - Numerator; k*(GeV);", 200, 0.0, 1.0);
  fResDen->Sumw2();

  fTrueNum = new TH1F("true_kstar_num", "KStar - Numerator (True Pairs); k*(GeV);", 200, 0.0, 1.0);
  fTrueNum->Sumw2();

  fTrueDen = new TH1F("true_kstar_den", "KStar - Denominator (True Pairs); k*(GeV);", 200, 0.0, 1.0);
  fTrueDen->Sumw2();
}

AliFemtoModelCorrFctnKStar::AliFemtoModelCorrFctnKStar(const char *name,
                                                       const int nbins,
                                                       const float KStarLo,
                                                       const float KStarHi):
  AliFemtoModelCorrFctn()
  , fResNum(NULL)
  , fResDen(NULL)
  , fTrueNum(NULL)
  , fTrueDen(NULL)
{
  fResNum = new TH1F(
    TString::Format("%s_Num", name),
    "KStar - Numerator; k*(GeV);",
    nbins, KStarLo, KStarHi);
  fResNum->Sumw2();

  fResDen = new TH1F(
    TString::Format("%s_Den", name),
    "KStar - Denominator; k*(GeV);",
    nbins, KStarLo, KStarHi);
  fResDen->Sumw2();

  fTrueNum = new TH1F(
    TString::Format("true_%s_Num", name),
     "KStar - Numerator (True Pairs); k*(GeV);",
    nbins, KStarLo, KStarHi);
  fTrueNum->Sumw2();

  fTrueDen = new TH1F(
    TString::Format("true_%s_Den", name),
    "KStar - Denominator; k*(GeV);",
    nbins, KStarLo, KStarHi);
  fTrueDen->Sumw2();
}

AliFemtoModelCorrFctnKStar::AliFemtoModelCorrFctnKStar(const AliFemtoModelCorrFctnKStar &orig):
  AliFemtoModelCorrFctn(orig)
  , fPairType(orig.fPairType)
  , fExpectedTrack1Code(orig.fExpectedTrack1Code)
  , fExpectedTrack2Code(orig.fExpectedTrack2Code)
  , fResNum(new TH1F(*orig.fResNum))
  , fResDen(new TH1F(*orig.fResDen))
  , fTrueNum(new TH1F(*orig.fTrueNum))
  , fTrueDen(new TH1F(*orig.fTrueDen))
{
}

AliFemtoModelCorrFctnKStar::~AliFemtoModelCorrFctnKStar()
{
  delete fResNum;
  delete fResDen;
  delete fTrueNum;
  delete fTrueDen;
}

AliFemtoString
AliFemtoModelCorrFctnKStar::Report()
{
  TString report;
  return AliFemtoString(report);
}


TList* AliFemtoModelCorrFctnKStar::GetOutputList()
{
  return AppendOutputList(new TList());
}

TList* AliFemtoModelCorrFctnKStar::AppendOutputList(TList *output_list)
{
  output_list->Add(fResNum);
  output_list->Add(fResDen);
  output_list->Add(fTrueNum);
  output_list->Add(fTrueDen);
  return output_list;
}

bool AliFemtoModelCorrFctnKStar::PairContainsExpectedTypes(const AliFemtoPair *pair)
{
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelAllHiddenInfo*>(pair->Track1()->HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelAllHiddenInfo*>(pair->Track2()->HiddenInfo());

  return PairContainsExpectedTypes(info1, info2);
}

bool AliFemtoModelCorrFctnKStar::PairContainsExpectedTypes(
  const AliFemtoModelHiddenInfo *info1,
  const AliFemtoModelHiddenInfo *info2)
{
  const bool track_1_is_expected = fExpectedTrack1Code == info1->GetPDGPid(),
             track_2_is_expected = fExpectedTrack2Code == info2->GetPDGPid();

  return track_1_is_expected && track_2_is_expected;
}


void AliFemtoModelCorrFctnKStar::AddRealPair(AliFemtoPair* aPair)
{
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track1()->HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track2()->HiddenInfo());

  if (info1 == NULL || info2 == NULL) {
    // cout << "NULL ("
    //      << aPair->Track1()->HiddenInfo() << " "
    //      << aPair->Track2()->HiddenInfo() << ")\n";
    return;
  }

  const Float_t kstar = CalcTrueKStar(aPair);
  fResNum->Fill(kstar);

  // track type not expected - skip the 'true' histogram
  if (!PairContainsExpectedTypes(info1, info2)) {
    return;
  }

  fTrueNum->Fill(kstar);
}

void AliFemtoModelCorrFctnKStar::AddMixedPair(AliFemtoPair* aPair)
{
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track1()->HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track2()->HiddenInfo());

  if (info1 == NULL || info2 == NULL) {
    // cout << "NULL ("
    //      << aPair->Track1()->HiddenInfo() << " "
    //      << aPair->Track2()->HiddenInfo() << ")\n";
    return;
  }

  const Float_t kstar = CalcTrueKStar(aPair);
  fResDen->Fill(kstar);
  // track type not expected - skip the 'true' histogram
  if (!PairContainsExpectedTypes(info1, info2)) {
    return;
  }

  fTrueDen->Fill(kstar);
}
