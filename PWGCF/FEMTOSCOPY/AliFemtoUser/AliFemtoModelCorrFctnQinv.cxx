///
/// \file AliFemtoModelCorrFctnQinv.cxx
/// \author Andrew Kubera
///

#include "AliFemtoModelCorrFctnQinv.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelAllHiddenInfo.h"
#include "AliFemtoPair.h"

#include <TH1F.h>
#include <TH2F.h>

#include <vector>
#include <string>
#include <map>

const static int UNSPECIFIED_PARTICLE = -999999999;


static
const int xi_code = 3312
        , xi0_code = 3322
        , sigma0_code = 3212
        , omega_code = 3334
        , lambda_c_code = 4122
        , sigma0_c_code = 3214
        , sigma_c_code = 3224
        ;

const std::map<Int_t, std::string> code_to_label = {
  {UNSPECIFIED_PARTICLE, "other"},
  {0, "none"},
  {11, "e^{-}"},
  {-11, "e^{+}"},
  {13, "#mu^{-}"},
  {-13, "#mu^{+}"},
  {15, "#tau^{-}"},
  {-15, "#tau^{+}"},
  {-211, "#pi^{-}"},
  {211, "#pi^{+}"},
  {311, "K^{0}"},
  {321, "K^{+}"},
  {-321, "K^{-}"},
  {113, "#rho^{0}"},
  {213, "#rho^{+}"},
  {-213, "#rho^{-}"},
  {221, "#eta"},
  {223, "#omega"}, // is this right?
  {2212, "p^{+}"},
  {-2212, /* "\\bar{p}" */ "p^{-}"},
  {3122, "#Lambda"},
  {-3122, "$\\bar{\\Lambda}$"},
  {3222, "#Sigma^{+}"},
  {sigma0_code, "#Sigma^{0}"},
  {3112, "#Sigma^{-}"},
  {xi_code, "#xi"},
  {xi0_code, "#xi^{0}"},
  {omega_code, "#Omega"},
  {lambda_c_code, "#lambda_{c}"},
  {sigma0_c_code, "#sigma^{0}_{c}"},
  {sigma_c_code, "#sigma_{c}"},
};

const std::vector<int> true_type_codes = {
  UNSPECIFIED_PARTICLE,
  0,
  xi_code, xi0_code, sigma0_code, omega_code, lambda_c_code, sigma_c_code
};

// use default values
AliFemtoModelCorrFctnQinv::AliFemtoModelCorrFctnQinv():
  AliFemtoModelCorrFctnQinv("qinv", 200, 0.0, 1.0)
{
}

AliFemtoModelCorrFctnQinv::AliFemtoModelCorrFctnQinv(const char *name,
                                                     const int nbins,
                                                     const float KStarLo,
                                                     const float KStarHi):
  AliFemtoModelCorrFctn()
  , fResNum(nullptr)
  , fResDen(nullptr)
  , fTrueNum(nullptr)
  , fTrueDen(nullptr)
{
  fResNum = new TH1F(
    TString::Format("%s_Num", name),
    "q_{inv} - Numerator; q_{inv} (GeV);",
    nbins, KStarLo, KStarHi
  );
  fResNum->Sumw2();

  fResDen = new TH1F(
    TString::Format("%s_Den", name),
    "q_{inv} - Denominator; q_{inv} (GeV);",
    nbins, KStarLo, KStarHi
  );
  fResDen->Sumw2();

  const float binstart = -0.5,
              binstop = true_type_codes.size() - 0.5;

  fTrueNum = new TH2F(
    TString::Format("true_%s_Num", name),
    "q_{inv} - Numerator (True Pairs); q_{inv} (GeV);",
    nbins, KStarLo, KStarHi,
    true_type_codes.size(), binstart, binstop
  );
  fTrueNum->Sumw2();

  fTrueDen = new TH2F(
    TString::Format("true_%s_Den", name),
    "q_{inv} - Denominator; q_{inv} (GeV);",
    nbins, KStarLo, KStarHi,
    true_type_codes.size(), binstart, binstop
  );
  fTrueDen->Sumw2();

  for (size_t bin = 0; bin < true_type_codes.size(); ++bin) {
    Int_t code = true_type_codes[bin];
    auto label = code_to_label.find(code)->second.c_str();
    fTrueNum->GetYaxis()->SetBinLabel(bin, label);
    fTrueDen->GetYaxis()->SetBinLabel(bin, label);
  }
}

AliFemtoModelCorrFctnQinv::AliFemtoModelCorrFctnQinv(const AliFemtoModelCorrFctnQinv &orig):
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

AliFemtoModelCorrFctnQinv::~AliFemtoModelCorrFctnQinv()
{
  delete fResNum;
  delete fResDen;
  delete fTrueNum;
  delete fTrueDen;
}

AliFemtoString
AliFemtoModelCorrFctnQinv::Report()
{
  TString report;
  return AliFemtoString(report);
}


TList* AliFemtoModelCorrFctnQinv::GetOutputList()
{
  return AppendOutputList(new TList());
}

TList* AliFemtoModelCorrFctnQinv::AppendOutputList(TList *output_list) const
{
  output_list->Add(fResNum);
  output_list->Add(fResDen);
  output_list->Add(fTrueNum);
  output_list->Add(fTrueDen);
  return output_list;
}

bool AliFemtoModelCorrFctnQinv::PairContainsExpectedTypes(const AliFemtoPair *pair)
{
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelAllHiddenInfo*>(pair->Track1()->HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelAllHiddenInfo*>(pair->Track2()->HiddenInfo());

  return PairContainsExpectedTypes(info1, info2);
}

bool AliFemtoModelCorrFctnQinv::PairContainsExpectedTypes(
  const AliFemtoModelHiddenInfo *info1,
  const AliFemtoModelHiddenInfo *info2)
{
  const bool track_1_is_expected = fExpectedTrack1Code == info1->GetPDGPid(),
             track_2_is_expected = fExpectedTrack2Code == info2->GetPDGPid();

  return track_1_is_expected && track_2_is_expected;
}


static inline int GetTruthBinFrom(const AliFemtoModelHiddenInfo *info)
{
  const auto code = abs(info->GetMotherPdgCode());
  const auto type_bin = std::find(true_type_codes.begin(), true_type_codes.end(), code);
  if (type_bin == true_type_codes.end()) {
    return 0;
  } else {
    return std::distance(true_type_codes.begin(), type_bin);
  }
}


void AliFemtoModelCorrFctnQinv::AddPair(const AliFemtoPair *pair,
                                        TH1F* combined_hist,
                                        TH2F* tue_hist)
{
  // const AliFemtoModelHiddenInfo
  //   *info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track1()->HiddenInfo()),
  //   *info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(aPair->Track2()->HiddenInfo());
  const AliFemtoModelHiddenInfo
    *info1 = (const AliFemtoModelHiddenInfo*)pair->Track1()->HiddenInfo(),
    *info2 = (const AliFemtoModelHiddenInfo*)pair->Track2()->HiddenInfo();

  if (info1 == nullptr || info2 == nullptr) {
    return;
  }

  // const Float_t q = CalcTrueQinv(pair);

  const Float_t mass1 = info1->GetMass(),
                mass2 = info2->GetMass();

  // block all zero-mass particles from the correlation function
  if (mass1 == 0.0 || mass2 == 0.0) {
    return;
  }

  const AliFemtoThreeVector *momentum1 = info1->GetTrueMomentum(),
                            *momentum2 = info2->GetTrueMomentum();

  const Float_t e1 = sqrt(mass1 * mass1 + momentum1->Mag2()),
                e2 = sqrt(mass2 * mass2 + momentum2->Mag2());

  const AliFemtoLorentzVector p1 = AliFemtoLorentzVector(e1, *momentum1),
                              p2 = AliFemtoLorentzVector(e2, *momentum2);

  const Float_t q = CalcQinv(p1, p2);
  combined_hist->Fill(q);

  const int truth_bin = GetTruthBinFrom(info1);
  tue_hist->Fill(q, truth_bin);
}

void AliFemtoModelCorrFctnQinv::AddRealPair(AliFemtoPair* aPair)
{
  AddPair(aPair, fResNum, fTrueNum);
}

void AliFemtoModelCorrFctnQinv::AddMixedPair(AliFemtoPair* aPair)
{
  AddPair(aPair, fResDen, fTrueDen);
}
