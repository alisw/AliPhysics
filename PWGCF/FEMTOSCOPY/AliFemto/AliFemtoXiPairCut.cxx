#include "AliFemtoXiPairCut.h"

#include <string>
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiPairCut);
  /// \endcond
#endif

//__________________
AliFemtoXiPairCut::AliFemtoXiPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fDataType(kAOD)
{
  /* no-op */
}
//__________________
AliFemtoXiPairCut::~AliFemtoXiPairCut()
{
  /* no-op */
}

AliFemtoXiPairCut &AliFemtoXiPairCut::operator=(const AliFemtoXiPairCut &cut)
{
  if (this == &cut) {
    return *this;
  }

  AliFemtoPairCut::operator=(cut);
  fNPairsPassed = cut.fNPairsPassed;
  fNPairsFailed = cut.fNPairsFailed;
  fDataType = cut.fDataType;


  return *this;
}

//__________________
bool AliFemtoXiPairCut::Pass(const AliFemtoPair *pair)
{
  const AliFemtoXi *Xi_1 = pair->Track1()->Xi(),
                   *Xi_2 = pair->Track2()->Xi();

  // Assert pair is of two Xi particles
  if (Xi_1 == NULL || Xi_2 == NULL) {
    return false;
  }

  // Assert that both particles do not share daughters
  if (Xi_1->IdNeg() == Xi_2->IdNeg() || Xi_1->IdPos() == Xi_2->IdPos()) {
    return false;
  }

  // Assert that both Xi's do not share bachelor track
  if (Xi_1->IdBac() == Xi_2->IdBac()) {
    return false;
  }

  if (Xi_1->IdBac() == Xi_1->IdPos() || Xi_1->IdBac() == Xi_1->IdNeg()) {
    return false;
  }

  if (Xi_2->IdBac() == Xi_2->IdPos() || Xi_2->IdBac() == Xi_2->IdNeg()) {
    return false;
  }


  return true;
}
//__________________
AliFemtoString AliFemtoXiPairCut::Report()
{
  TString report = "AliFemtoXi Pair Cut - remove shared and split pairs\n";
  report += TString::Format("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);

  return AliFemtoString((const char *)report);
}
//__________________



TList *AliFemtoXiPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  // The TString format patterns (F is float, I is integer, L is long)
  const TString prefix = "AliFemtoXiPairCut.";

  tListSetttings->Add(new TObjString(prefix + Form("datatype=%d", fDataType)));
  tListSetttings->Add(new TObjString(prefix + Form("pairs_passed=%ld", fNPairsPassed)));
  tListSetttings->Add(new TObjString(prefix + Form("pairs_failed=%ld", fNPairsFailed)));

  return tListSetttings;
}


void AliFemtoXiPairCut::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}
