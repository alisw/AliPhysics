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
  11,
  13,
  -211,
  xi_code, xi0_code, sigma0_code, omega_code, lambda_c_code, sigma_c_code
};

// use default values
AliFemtoModelCorrFctnQinv::AliFemtoModelCorrFctnQinv():
  AliFemtoModelCorrFctnQinv("qinv", 200, 0.0, 1.0)
{
}

AliFemtoModelCorrFctnQinv::AliFemtoModelCorrFctnQinv(const char *suffix,
                                                     const int nbins,
                                                     const float qinv_low_limit,
                                                     const float qinv_high_limit)
  : AliFemtoModelCorrFctn(suffix, nbins, qinv_low_limit, qinv_high_limit)
  , fPairType(AliFemtoAvgSepCorrFctn::kTracks)
  , fExpectedTrack1Code(0)
  , fExpectedTrack2Code(0)
  , fNumPid(nullptr)
  , fDenPid(nullptr)
{
  fNumeratorTrue->SetTitle("Reconstructed Numerator; q_{inv} (Gev)");
  fNumeratorFake->SetTitle("Simulated Reconstructed Numerator; q_{inv} (Gev)");
  fDenominator->SetTitle("Reconstructed Denominator; q_{inv} (Gev)");

  fNumeratorTrueIdeal->SetTitle("Generated Numerator; q_{inv} (Gev)");
  fNumeratorFakeIdeal->SetTitle("Simulated Generated Numerator; q_{inv} (Gev)");
  fDenominatorIdeal->SetTitle("Generated Denominator; q_{inv} (Gev)");

  fDenominatorIdeal->Sumw2(false);
  fDenominator->Sumw2(false);

  fQgenQrec->SetTitle("Q_{inv,reconstructed} vs Q_{inv,generated};"
                      "q_{inv,gen} (GeV);"
                      "q_{inv,rec} (GeV);"
                      "N_pairs");

  /*
  const float binstart = -0.5,
              binstop = true_type_codes.size() - 0.5;

  auto _build_1d_hist = [&](const TString &name, const TString &title)
    {
      return new TH1F(name + suffix, title,
                      nbins, qinv_low_limit, qinv_high_limit);
    };

  auto _build_2d_hist = [&](const TString &name, const TString &title)
    {
      return new TH2F(name + suffix, title,
                      nbins, qinv_low_limit, qinv_high_limit,
                      true_type_codes.size(), binstart, binstop);
    };

  fNumPid = _build_2d_hist("TrueNum", "True q_{inv} - Numerator; q_{inv} (GeV);");
  fDenPid = _build_2d_hist("TrueDen", "True q_{inv} - Denominator; q_{inv} (GeV);");

  for (size_t bin = 0; bin < true_type_codes.size(); ++bin) {
    Int_t code = true_type_codes[bin];
    auto label = code_to_label.find(code)->second.c_str();
    fNumPid->GetYaxis()->SetBinLabel(bin, label);
    fDenPid->GetYaxis()->SetBinLabel(bin, label);
  }
  */
}

AliFemtoModelCorrFctnQinv::AliFemtoModelCorrFctnQinv(const AliFemtoModelCorrFctnQinv &orig):
  AliFemtoModelCorrFctn(orig)
  , fPairType(orig.fPairType)
  , fExpectedTrack1Code(orig.fExpectedTrack1Code)
  , fExpectedTrack2Code(orig.fExpectedTrack2Code)
  , fNumPid(nullptr) // new TH2F(*orig.fNumPid))
  , fDenPid(nullptr) // new TH2F(*orig.fDenPid))
{
}

AliFemtoModelCorrFctn*
AliFemtoModelCorrFctnQinv::Clone() const
{
  AliFemtoModelCorrFctnQinv *result = new AliFemtoModelCorrFctnQinv(*this);
  // result->fNumPid->Reset();
  // result->fDenPid->Reset();
  return result;
};


AliFemtoModelCorrFctnQinv&
AliFemtoModelCorrFctnQinv::operator=(const AliFemtoModelCorrFctnQinv &rhs)
{
  if (this == &rhs) {
    return *this;
  }

  fPairType = rhs.fPairType;

  fExpectedTrack1Code = rhs.fExpectedTrack1Code;
  fExpectedTrack2Code = rhs.fExpectedTrack2Code;

  // *fNumPid = *rhs.fNumPid;
  // *fDenPid = *rhs.fDenPid;

  return *this;
}

AliFemtoModelCorrFctnQinv::~AliFemtoModelCorrFctnQinv()
{
  delete fNumPid;
  delete fDenPid;
}

AliFemtoString
AliFemtoModelCorrFctnQinv::Report()
{
  AliFemtoString report;
  return report;
}


TList* AliFemtoModelCorrFctnQinv::GetOutputList()
{
  TList *result = AliFemtoModelCorrFctn::GetOutputList();
  AppendOutputList(result);
  return result;
}

TList* AliFemtoModelCorrFctnQinv::AppendOutputList(TList *output_list) const
{
  // output_list->Add(fNumPid);
  // output_list->Add(fDenPid);
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


void AliFemtoModelCorrFctnQinv::AddPair(const AliFemtoPair &pair,
                                        TH1F* combined_hist,
                                        TH2F* pid_hist)
{
  const AliFemtoModelHiddenInfo
    *info1 = static_cast<const AliFemtoModelHiddenInfo*>(pair.Track1()->HiddenInfo()),
    *info2 = static_cast<const AliFemtoModelHiddenInfo*>(pair.Track2()->HiddenInfo());

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
  pid_hist->Fill(q, truth_bin);
}

void AliFemtoModelCorrFctnQinv::AddRealPair(AliFemtoPair* aPair)
{
  AliFemtoModelCorrFctn::AddRealPair(aPair);
  // AddPair(*aPair, nullptr, fNumPid);
}

void AliFemtoModelCorrFctnQinv::AddMixedPair(AliFemtoPair* aPair)
{
  AliFemtoModelCorrFctn::AddMixedPair(aPair);
  // AddPair(*aPair, nullptr, fDenPid);
}
