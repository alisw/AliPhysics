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

#include <vector>
#include <string>

static const char* true_types[] =
{ "Other"     // 0
, "Primary"   // 1
, "Xi"        // 2
, "Xi0"       // 3
, "Sigma0"    // 4
, "Omega"     // 5
, "Lamba_{c}" // 6
, "Sigma_{c}" // 7
, NULL
};

static const int xi_bin = 2
        , xi0_bin = 3
        , sigma0_bin = 4
        , omega_bin = 5
        , lambda_c_bin = 6
        , sigma_c_bin = 7
        ;

static const int xi_code = 3312
        , xi0_code = 3322
        , sigma0_code = 3212
        , omega_code = 3334
        , lambda_c_code = 4122
        , sigma0_c_code = 3214
        , sigma_c_code = 3224
        ;

const size_t true_type_count = sizeof(true_types) / sizeof(true_types[0]);

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

  fTrueNum = new TH2F("true_kstar_num", "KStar - Numerator (True Pairs); k*(GeV);", 200, 0.0, 1.0, true_type_count, -0.5, true_type_count - 0.5);
  fTrueNum->Sumw2();

  fTrueDen = new TH2F("true_kstar_den", "KStar - Denominator (True Pairs); k*(GeV);", 200, 0.0, 1.0, true_type_count, -0.5, true_type_count - 0.5);
  fTrueDen->Sumw2();

  for (size_t i = 0; i < true_type_count; ++i) {
    fTrueNum->GetYaxis()->SetBinLabel(i+1, true_types[i]);
    fTrueDen->GetYaxis()->SetBinLabel(i+1, true_types[i]);
  }
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
    nbins, KStarLo, KStarHi
  );
  fResNum->Sumw2();

  fResDen = new TH1F(
    TString::Format("%s_Den", name),
    "KStar - Denominator; k*(GeV);",
    nbins, KStarLo, KStarHi
  );
  fResDen->Sumw2();

  fTrueNum = new TH2F(
    TString::Format("true_%s_Num", name),
    "KStar - Numerator (True Pairs); k*(GeV);",
    nbins, KStarLo, KStarHi,
    true_type_count, -0.5, true_type_count - 0.5
  );
  fTrueNum->Sumw2();

  fTrueDen = new TH2F(
    TString::Format("true_%s_Den", name),
    "KStar - Denominator; k*(GeV);",
    nbins, KStarLo, KStarHi,
    true_type_count, -0.5, true_type_count - 0.5
  );
  fTrueDen->Sumw2();

  for (size_t i = 0; i < true_type_count; ++i) {
    fTrueNum->GetYaxis()->SetBinLabel(i+1, true_types[i]);
    fTrueDen->GetYaxis()->SetBinLabel(i+1, true_types[i]);
  }
}

AliFemtoModelCorrFctnKStar::AliFemtoModelCorrFctnKStar(const AliFemtoModelCorrFctnKStar &orig):
  AliFemtoModelCorrFctn(orig)
  , fPairType(orig.fPairType)
  , fExpectedTrack1Code(orig.fExpectedTrack1Code)
  , fExpectedTrack2Code(orig.fExpectedTrack2Code)
  , fResNum(new TH1F(*orig.fResNum))
  , fResDen(new TH1F(*orig.fResDen))
  , fTrueNum(new TH2F(*orig.fTrueNum))
  , fTrueDen(new TH2F(*orig.fTrueDen))
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


static inline int GetTruthBinFrom(const AliFemtoModelHiddenInfo *info)
{
  switch (abs(info->GetMotherPdgCode())) {
    default: return 0;
    case 0: return 1;
    case xi_code: return xi_bin;
    case xi0_code: return xi0_bin;
    case sigma0_code: return sigma0_bin;
    case omega_code: return omega_bin;
    case lambda_c_code: return lambda_c_bin;
    case sigma0_c_code:
    case sigma_c_code: return sigma_c_bin;
  }
}

void AliFemtoModelCorrFctnKStar::AddRealPair(AliFemtoPair* aPair)
{
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track1()->HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track2()->HiddenInfo());

  if (info1 == NULL || info2 == NULL) {
    return;
  }

  const Float_t kstar = CalcTrueKStar(aPair);
  fResNum->Fill(kstar);

  const int truth_bin = GetTruthBinFrom(info1);
  fTrueNum->Fill(kstar, truth_bin);
}

void AliFemtoModelCorrFctnKStar::AddMixedPair(AliFemtoPair* aPair)
{
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track1()->HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track2()->HiddenInfo());

  if (info1 == NULL || info2 == NULL) {
    return;
  }

  const Float_t kstar = CalcTrueKStar(aPair);
  fResDen->Fill(kstar);

  const int truth_bin = GetTruthBinFrom(info1);
  fTrueDen->Fill(kstar, truth_bin);
}
