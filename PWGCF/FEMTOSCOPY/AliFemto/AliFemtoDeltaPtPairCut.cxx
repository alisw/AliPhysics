///
/// \file AliFemtoDeltaPtPairCut.cxx
///

#include "AliFemtoDeltaPtPairCut.h"
#include <string>
#include <cstdio>

#include <TMath.h>
#include <TString.h>
#include <TObjString.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoDeltaPtPairCut);
  /// \endcond
#endif

//__________________
AliFemtoDeltaPtPairCut::AliFemtoDeltaPtPairCut():
  AliFemtoPairCut(),
  fDeltaPtMin(0.0),
  fDeltaPtMax(0.0)
{

}
//__________________
AliFemtoDeltaPtPairCut::AliFemtoDeltaPtPairCut(double lo, double hi) :
  AliFemtoPairCut(),
  fDeltaPtMin(lo),
  fDeltaPtMax(hi)
{  /* no-op */
}
//__________________
AliFemtoDeltaPtPairCut::AliFemtoDeltaPtPairCut(const AliFemtoDeltaPtPairCut& c) :
  AliFemtoPairCut(c),
  fDeltaPtMin(c.fDeltaPtMin),
  fDeltaPtMax(c.fDeltaPtMax)
{  /* no-op */
}

//__________________
AliFemtoDeltaPtPairCut::~AliFemtoDeltaPtPairCut()
{  /* no-op */
}

AliFemtoDeltaPtPairCut& AliFemtoDeltaPtPairCut::operator=(const AliFemtoDeltaPtPairCut& c)
{
  if (this != &c) {
    AliFemtoPairCut::operator=(c);
    fDeltaPtMin = c.fDeltaPtMin;
    fDeltaPtMax = c.fDeltaPtMax;
  }

  return *this;
}
//__________________
AliFemtoString AliFemtoDeltaPtPairCut::Report()
{
  // Prepare a report from the execution
  TString report("AliFemtoKT Pair Cut\n");
  report += TString::Format("Accept pair with DeltaPt in range %f , %f",fDeltaPtMin,fDeltaPtMax);

  return AliFemtoString((const char *)report);
}
//__________________

TList *AliFemtoDeltaPtPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  TString next_setting;

  next_setting = TString::Format("AliFemtoDeltaPtPairCut.ktmax=%f", fDeltaPtMax);
  tListSetttings->Add(new TObjString(next_setting));

  next_setting = TString::Format("AliFemtoDeltaPtPairCut.ktmin=%f", fDeltaPtMin);
  tListSetttings->Add(new TObjString(next_setting));

  return tListSetttings;
}

void AliFemtoDeltaPtPairCut::SetDeltaPtRange(double ptmin, double ptmax)
{
  fDeltaPtMin = ptmin;
  fDeltaPtMax = ptmax;
}
//______________________________________________________
bool AliFemtoDeltaPtPairCut::Pass(const AliFemtoPair* pair)
{

  const double pT1 = pair->Track1()->Track()->Pt(),
               pT2 = pair->Track2()->Track()->Pt();

  const double DeltaPt = TMath::Abs(pT1 - pT2);

  bool passes = (fDeltaPtMin <= DeltaPt && DeltaPt <= fDeltaPtMax);

  return passes;
}
