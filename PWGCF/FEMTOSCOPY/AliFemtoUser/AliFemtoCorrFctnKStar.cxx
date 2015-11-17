///
/// \file AliFemtoCorrFctnKStar.cxx
///

#include "AliFemtoCorrFctnKStar.h"

#include "TH1D.h"
#include "TH2F.h"
#include "AliFemtoCorrFctn.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctnKStar);
  /// \endcond
#endif

const double LambdaMass = 1.115683,
               PionMass = 0.13956995;


const float kstar_min = 0.0, kstar_max = 1.0,
               kt_min = 0.0,    kt_max = 3.5,
               mt_min = 1.5,    mt_max = 10.5; 

const float q_o_min = 0.0, q_o_max = 7.5,
            q_s_min = -5.0, q_s_max = 5.0;

AliFemtoCorrFctnKStar::AliFemtoCorrFctnKStar():
  AliFemtoCorrFctn()
  , fNumerator(new TH1D("kstar_num",
                        "KStar - Numerator; k*(GeV);",
                        200, kstar_min, kstar_max))
  , fDenominator(new TH1D("kstar_den",
                          "KStar - Denominator; k*(GeV);",
                          200, kstar_min, kstar_max))

  , fNumerator_kT(new TH2F("kstar_num_kt",
                           "KStar vs kT Numerator; k*(GeV); k_T (GeV);",
                           200, kstar_min, kstar_max,
                           128, kt_min, kt_max))
  , fDenominator_kT(new TH2F("kstar_den_kt",
                             "KStar vs kT Denominator; k* (GeV); k_T (GeV)",
                             200, kstar_min, kstar_max,
                             128, kt_min, kt_max))

  , fNumerator_mT(new TH2F("kstar_num_mt",
                           "KStar vs m_{T} Numerator; k*(GeV); m_{T} (GeV);",
                           200, kstar_min, kstar_max,
                           128, mt_min, mt_max))
  , fDenominator_mT(new TH2F("kstar_den_mt",
                             "KStar vs m_{T} Denominator; k* (GeV); m_{T} (GeV)",
                             200, kstar_min, kstar_max,
                             128, mt_min, mt_max))

  , fNumerator_qq(new TH2F("qo_qs_num",
                           "Q_{out} vs Q_{side} Numerator;"
                           "Q_{side} (GeV); Q_{out} (GeV);",
                           200, q_s_min, q_s_max,
                           200, q_o_min, q_o_max))
  , fDenominator_qq(new TH2F("qo_qs_den",
                           "Q_{out} vs Q_{side} Denominator;"
                           "Q_{side} (GeV); Q_{out} (GeV);",
                           200, q_s_min, q_s_max,
                           200, q_o_min, q_o_max))

{
  fNumerator->Sumw2();
  fDenominator->Sumw2();

  fNumerator_kT->Sumw2();
  fDenominator_kT->Sumw2();

  fNumerator_mT->Sumw2();
  fDenominator_mT->Sumw2();
  
  fNumerator_qq->Sumw2();
  fDenominator_qq->Sumw2();
}

AliFemtoCorrFctnKStar::AliFemtoCorrFctnKStar(const char *title,
                                             const int nbins,
                                             const float KStarLo,
                                             const float KStarHi):
  AliFemtoCorrFctn()
  , fNumerator(new TH1D(TString::Format("%s_Num", title),
                        "KStar - Numerator; k*(GeV);",
                         nbins, KStarLo, KStarHi))
  , fDenominator(new TH1D(TString::Format("%s_Den", title),
                        "KStar - Denominator; k*(GeV);",
                        nbins, KStarLo, KStarHi))

  , fNumerator_kT(new TH2F("kstar_num_kt",
                         "KStar vs kT Numerator; k*(GeV); k_T (GeV);",
                         nbins, KStarLo, KStarHi,
                         128, kt_min, kt_max))
  , fDenominator_kT(new TH2F("kstar_den_kt",
                         "KStar vs kT Denominator; k* (GeV); k_T (GeV)",
                         nbins, KStarLo, KStarHi,
                         128, kt_min, kt_max))

  , fNumerator_mT(new TH2F("kstar_num_mt",
                         "KStar vs m_{T} Numerator; k*(GeV); m_{T} (GeV);",
                         nbins, KStarLo, KStarHi,
                         128, mt_min, mt_max))
  , fDenominator_mT(new TH2F("kstar_den_mt",
                         "KStar vs m_{T} Denominator; k* (GeV); m_{T} (GeV)",
                         nbins, KStarLo, KStarHi,
                         128, mt_min, mt_max))

  , fNumerator_qq(new TH2F("qo_qs_num",
                           "Q_{out} vs Q_{side} Numerator;"
                           "Q_{side} (GeV); Q_{out} (GeV);",
                           200, q_s_min, q_s_max,
                           200, q_o_min, q_o_max))
  , fDenominator_qq(new TH2F("qo_qs_den",
                           "Q_{out} vs Q_{side} Denominator;"
                           "Q_{side} (GeV); Q_{out} (GeV);",
                           200, q_s_min, q_s_max,
                           200, q_o_min, q_o_max))

{
  fNumerator->Sumw2();
  fDenominator->Sumw2();

  fNumerator_kT->Sumw2();
  fDenominator_kT->Sumw2();

  fNumerator_mT->Sumw2();
  fDenominator_mT->Sumw2();

  fNumerator_qq->Sumw2();
  fDenominator_qq->Sumw2();
}


AliFemtoString
AliFemtoCorrFctnKStar::Report()
{
  TString report;
  return AliFemtoString(report);
}

void AliFemtoCorrFctnKStar::Finish()
{ // no-op
}

TList* AliFemtoCorrFctnKStar::GetOutputList()
{
  TList *output_list = new TList();
  output_list->Add(fNumerator);
  output_list->Add(fDenominator);
  output_list->Add(fNumerator_kT);
  output_list->Add(fDenominator_kT);
  output_list->Add(fNumerator_mT);
  output_list->Add(fDenominator_mT);
  output_list->Add(fNumerator_qq);
  output_list->Add(fDenominator_qq);
  return output_list;
}

void AliFemtoCorrFctnKStar::AddRealPair(AliFemtoPair* aPair)
{
  fNumerator->Fill(aPair->KStar());
  fNumerator_kT->Fill(aPair->KStar(), aPair->KT());
  fNumerator_mT->Fill(aPair->KStar(), CalcMt(aPair));
  fNumerator_qq->Fill(aPair->QSideCMS(), aPair->QOutCMS());
}

void AliFemtoCorrFctnKStar::AddMixedPair(AliFemtoPair* aPair)
{
  fDenominator->Fill(aPair->KStar());
  fDenominator_kT->Fill(aPair->KStar(), aPair->KT());
  fDenominator_mT->Fill(aPair->KStar(), CalcMt(aPair));
  fDenominator_qq->Fill(aPair->QSideCMS(), aPair->QOutCMS());
}


float AliFemtoCorrFctnKStar::CalcMt(const AliFemtoPair* aPair)
{
  const double mass_lam_2 = LambdaMass * LambdaMass,
                mass_pi_2 = PionMass * PionMass;
  
  AliFemtoThreeVector p1 = aPair->Track1()->V0()->MomV0(),
                      p2 = aPair->Track2()->Track()->P();
  
  const float et1 = TMath::Sqrt(mass_lam_2 + p1.x() * p1.x() + p1.y() * p1.y()),
              et2 = TMath::Sqrt(mass_pi_2 + p2.x() * p2.x() + p2.y() * p2.y()),
               pt = p1.x() * p2.x() + p1.y() * p2.y();
  
  const float mt = mass_lam_2 + mass_pi_2 + 2 * (et1 * et2 - pt);
  return mt;
}
