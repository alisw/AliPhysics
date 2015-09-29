#include "AliFemtoXiV0PairCut.h"

#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoXiV0PairCut)
#endif

//__________________
AliFemtoXiV0PairCut::AliFemtoXiV0PairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fDataType(kAOD)
{
  /* no-op */
}
//__________________
AliFemtoXiV0PairCut::~AliFemtoXiV0PairCut()
{
  /* no-op */
}

AliFemtoXiV0PairCut &AliFemtoXiV0PairCut::operator=(const AliFemtoXiV0PairCut &cut)
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
bool AliFemtoXiV0PairCut::Pass(const AliFemtoPair *pair) 
{
  const AliFemtoXi *Xi = pair->Track1()->Xi();
  const AliFemtoV0 *V0 = pair->Track2()->V0();

  // Assert pair is of Xi-V0 particles
  if (Xi == NULL || V0 == NULL) {
    return false;
  }

  // Assert that both particles do not share daughters
  if (Xi->IdNeg() == V0->IdNeg() || Xi->IdPos() == V0->IdPos()) {
    return false;
  }

  if (Xi->IdBac() == Xi->IdPos() || Xi->IdBac() == Xi->IdNeg()) {
    return false;
  }

  if (Xi->IdBac() == V0->IdPos() || Xi->IdBac() == V0->IdNeg()) {
    return false;
  }


  return true;
}
//__________________
AliFemtoString AliFemtoXiV0PairCut::Report()
{
  TString report = "AliFemtoXi Pair Cut - remove shared and split pairs\n";
  report += TString::Format("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);

  return AliFemtoString(report);
}
//__________________



TList *AliFemtoXiV0PairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  // The TString format patterns (F is float, I is integer, L is long)
  const char ptrnF[] = "AliFemtoXiV0PairCut.%s=%f",
             ptrnI[] = "AliFemtoXiV0PairCut.%s=%d",
             ptrnL[] = "AliFemtoXiV0PairCut.%s=%ld";


  tListSetttings->Add(new TObjString(TString::Format(ptrnI, "datatype", fDataType)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnL, "pairs_passed", fNPairsPassed)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnL, "pairs_failed", fNPairsFailed)));

  return tListSetttings;
}


void AliFemtoXiV0PairCut::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}
