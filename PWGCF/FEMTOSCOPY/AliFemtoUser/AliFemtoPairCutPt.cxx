///
/// \file AliFemtoPairCutPt.cxx
///

#include "AliFemtoPairCutPt.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoPairCutPt);
  /// \endcond
#endif

//__________________
AliFemtoPairCutPt::AliFemtoPairCutPt():
  AliFemtoPairCut(),
  fSumPtMin(0),
  fSumPtMax(10000),
  fNPairsFailed(0),
  fNPairsPassed(0)
{ /* no-op */
}
//__________________
AliFemtoPairCutPt::AliFemtoPairCutPt(double lo, double hi):
  AliFemtoPairCut(),
  fSumPtMin(lo),
  fSumPtMax(hi),
  fNPairsFailed(0),
  fNPairsPassed(0)
{ /* no-op */
}
//__________________
AliFemtoPairCutPt::AliFemtoPairCutPt(const AliFemtoPairCutPt& c):
  AliFemtoPairCut(c),
  fSumPtMin(c.fSumPtMin),
  fSumPtMax(c.fSumPtMax),
  fNPairsFailed(0),
  fNPairsPassed(0)
{ /* no-op */
}
AliFemtoPairCutPt& AliFemtoPairCutPt::operator=(const AliFemtoPairCutPt& c)
{
  if (this != &c) {
    AliFemtoPairCut::operator=(c);
    fSumPtMin = c.fSumPtMin;
    fSumPtMax = c.fSumPtMax;
  }

  return *this;
}

//__________________
AliFemtoPairCutPt::~AliFemtoPairCutPt()
{ /* no-op */
}
//__________________
bool AliFemtoPairCutPt::Pass(const AliFemtoPair* pair)
{
  const double pt1 = pair->Track1()->Track()->Pt(),
               pt2 = pair->Track2()->Track()->Pt();

  const double pt_sum = pt1 + pt2;

  // passes if pt_sum is between fSumPtMin and fSumPtMax
  const bool temp = (fSumPtMin <= pt_sum && pt_sum <= fSumPtMax);

  if (temp) {
    fNPairsPassed++;
  } else {
    fNPairsFailed++;
  }

  return temp;
}
//__________________
AliFemtoString AliFemtoPairCutPt::Report()
{
  // Prepare a report from the execution
  TString report("AliFemtoPairCutPt Pair Cut\n");
  report += TString::Format("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", (long int)fNPairsPassed, (long int)fNPairsFailed);

  return AliFemtoString((const char *)report);
}
//__________________

TList *AliFemtoPairCutPt::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  tListSetttings->Add(new TObjString(
    TString::Format("AliFemtoPairCutPt.sumptmin=%f", fSumPtMin)
  ));
  tListSetttings->Add(new TObjString(
    TString::Format("AliFemtoPairCutPr.sumptmax=%f", fSumPtMax)
  ));

  return tListSetttings;
}

void AliFemtoPairCutPt::SetMinSumPt(Double_t sumptmin)
{
  fSumPtMin = sumptmin;
}

void AliFemtoPairCutPt::SetMaxSumPt(Double_t sumptmax)
{
  fSumPtMax = sumptmax;
}
