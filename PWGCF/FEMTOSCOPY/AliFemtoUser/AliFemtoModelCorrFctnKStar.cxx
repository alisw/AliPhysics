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
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"

#include <TH1D.h>
#include <TH2D.h>


AliFemtoModelCorrFctnKStar::AliFemtoModelCorrFctnKStar():
  AliFemtoCorrFctn(),
  fNumerator(new TH1D("kstar_num", "KStar - Numerator; k*(GeV);", 200, 0.0, 1.0)),
  fDenominator(new TH1D("kstar_den", "KStar - Denominator; k*(GeV);", 200, 0.0, 1.0))
{
  fNumerator->Sumw2();
  fDenominator->Sumw2();
}

AliFemtoModelCorrFctnKStar::AliFemtoModelCorrFctnKStar(const char *title,
                                                       const int nbins,
                                                       const float KStarLo,
                                                       const float KStarHi):
  AliFemtoModelCorrFctnKStar()
  , fNumerator(new TH1D(TString::Format("%s_Num", title), "KStar - Numerator; k*(GeV);", nbins, KStarLo, KStarHi)),
  , fDenominator(new TH1D(TString::Format("%s_Den", title), "KStar - Denominator; k*(GeV);", nbins, KStarLo, KStarHi))
{
  fNumerator->Sumw2();
  fDenominator->Sumw2();
}

AliFemtoString
AliFemtoModelCorrFctnKStar::Report()
{
  TString report;
  return AliFemtoString(report);
}

void AliFemtoModelCorrFctnKStar::Finish()
{ // no-op
}

TList* AliFemtoModelCorrFctnKStar::GetOutputList()
{
  TList *output_list = new TList();
  output_list->Add(fNumerator);
  output_list->Add(fDenominator);
  return output_list;
}

void AliFemtoModelCorrFctnKStar::AddRealPair(AliFemtoPair* aPair)
{
  fNumerator->Fill(aPair->KStar());
}

void AliFemtoModelCorrFctnKStar::AddMixedPair(AliFemtoPair* aPair)
{
  fDenominator->Fill(aPair->KStar());
}
