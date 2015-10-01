///
/// \file AliFemtoCorrFctnKStar.cxx
///

#include "AliFemtoCorrFctnKStar.h"

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"

/// \cond CLASSIMP
ClassImp(AliFemtoCorrFctnKStar);
/// \endcond

AliFemtoCorrFctnKStar::AliFemtoCorrFctnKStar():
  AliFemtoCorrFctn(),
  fNumerator(new TH1D("kstar_num", "KStar - Numerator; k*(GeV);", 200, 0.0, 1.0)),
  fDenominator(new TH1D("kstar_den", "KStar - Denominator; k*(GeV);", 200, 0.0, 1.0))
{
  fNumerator->Sumw2();
  fDenominator->Sumw2();
}

AliFemtoCorrFctnKStar::AliFemtoCorrFctnKStar(const char *title,
                                             const int nbins,
                                             const float KStarLo,
                                             const float KStarHi):
  AliFemtoCorrFctn(),
  fNumerator(new TH1D(TString::Format("%s_Num", title), "KStar - Numerator; k*(GeV);", nbins, KStarLo, KStarHi)),
  fDenominator(new TH1D(TString::Format("%s_Den", title), "KStar - Denominator; k*(GeV);", nbins, KStarLo, KStarHi))
{
  fNumerator->Sumw2();
  fDenominator->Sumw2();
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
  return output_list;
}

void AliFemtoCorrFctnKStar::AddRealPair(AliFemtoPair* aPair)
{
  fNumerator->Fill(aPair->KStar());
}

void AliFemtoCorrFctnKStar::AddMixedPair(AliFemtoPair* aPair)
{
  fDenominator->Fill(aPair->KStar());
}
