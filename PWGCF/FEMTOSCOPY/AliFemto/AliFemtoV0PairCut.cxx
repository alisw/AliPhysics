///
/// \file AliFemtoV0PairCut.cxx
///

#include "AliFemtoV0PairCut.h"

#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoV0PairCut)
#endif

//__________________
AliFemtoV0PairCut::AliFemtoV0PairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fV0Max(1.0),
  fShareFractionMax(1.0),
  fRemoveSameLabel(0),
  fDataType(kAOD),
  fDTPCMin(0),
  fDTPCExitMin(0),
  fMinAvgSepPosPos(0),
  fMinAvgSepPosNeg(0),
  fMinAvgSepNegPos(0),
  fMinAvgSepNegNeg(0)
{
  /* no-op */
}
//__________________
AliFemtoV0PairCut::~AliFemtoV0PairCut()
{
  /* no-op */
}

AliFemtoV0PairCut &AliFemtoV0PairCut::operator=(const AliFemtoV0PairCut &cut)
{
  if (this == &cut) {
    return *this;
  }

  AliFemtoPairCut::operator=(cut);
  fNPairsPassed = cut.fNPairsPassed;
  fNPairsFailed = cut.fNPairsFailed;
  fV0Max = cut.fV0Max;
  fShareFractionMax = cut.fShareFractionMax;
  fRemoveSameLabel = cut.fRemoveSameLabel;
  fDataType = cut.fDataType;
  fDTPCMin = cut.fDTPCMin;
  fDTPCExitMin = cut.fDTPCExitMin;
  fMinAvgSepPosPos = cut.fMinAvgSepPosPos;
  fMinAvgSepPosNeg = cut.fMinAvgSepPosNeg;
  fMinAvgSepNegPos = cut.fMinAvgSepNegPos;
  fMinAvgSepNegNeg = cut.fMinAvgSepNegNeg;

  return *this;
}

//__________________
bool AliFemtoV0PairCut::Pass(const AliFemtoPair *pair)
{
  const AliFemtoV0 *V0_1 = pair->Track1()->V0(),
                   *V0_2 = pair->Track2()->V0();

  // Assert pair is of two V0 particles
  if (V0_1 == NULL || V0_2 == NULL) {
    return false;
  }

  // Assert that both particles do not share daughters
  if (V0_1->IdNeg() == V0_2->IdNeg() || V0_1->IdPos() == V0_2->IdPos()) {
    return false;
  }

  // Test separation between track daughters' entrance and exit points
  if (fDataType == kESD || fDataType == kAOD) {
    const AliFemtoThreeVector diffPosEntrance = V0_1->NominalTpcEntrancePointPos() - V0_2->NominalTpcEntrancePointPos(),
                                  diffPosExit = V0_1->NominalTpcExitPointPos() - V0_2->NominalTpcExitPointPos(),

                              diffNegEntrance = V0_1->NominalTpcEntrancePointNeg() - V0_2->NominalTpcEntrancePointNeg(),
                                  diffNegExit = V0_1->NominalTpcExitPointNeg() - V0_2->NominalTpcExitPointNeg();

    if (diffPosEntrance.Mag() < fDTPCMin ||
        diffNegEntrance.Mag() < fDTPCMin ||
        diffPosExit.Mag() < fDTPCExitMin ||
        diffNegExit.Mag() < fDTPCExitMin) {
      return false;
    }
  }

  // Find average separations between tracks

  // Remove const qualifier so we can use the non-const NominalTpcPoint functions (an oversight in AliFemtoV0 definition)
  AliFemtoV0 *V1 = const_cast<AliFemtoV0*>(V0_1),
             *V2 = const_cast<AliFemtoV0*>(V0_2);

  double avgSep = 0.0;

  for (int i = 0; i < 8; i++) {
    AliFemtoThreeVector delta = V1->NominalTpcPointPos(i) - V2->NominalTpcPointPos(i);
    avgSep += delta.Mag();
  }
  avgSep /= 8.0;

  if (avgSep < fMinAvgSepPosPos) {
    return false;
  }

  avgSep = 0.0;

  for (int i = 0; i < 8; i++) {
    AliFemtoThreeVector delta = V1->NominalTpcPointPos(i) - V2->NominalTpcPointNeg(i);
    avgSep += delta.Mag();
  }
  avgSep /= 8.0;
  if (avgSep < fMinAvgSepPosNeg) {
    return false;
  }

  avgSep = 0.0;

  for (int i = 0; i < 8 ; i++) {
    AliFemtoThreeVector delta = V1->NominalTpcPointNeg(i) - V2->NominalTpcPointPos(i);
    avgSep += delta.Mag();
  }
  avgSep /= 8.0;
  if (avgSep < fMinAvgSepNegPos) {
    return false;
  }

  avgSep = 0.0;

  for (int i = 0; i < 8 ; i++) {
    AliFemtoThreeVector delta = V1->NominalTpcPointNeg(i) - V2->NominalTpcPointNeg(i);
    avgSep += delta.Mag();
  }
  avgSep /= 8;
  if (avgSep < fMinAvgSepNegNeg) {
    return false;
  }

  return true;
}
//__________________
AliFemtoString AliFemtoV0PairCut::Report()
{
  TString report = "AliFemtoV0 Pair Cut - remove shared and split pairs\n";
  report += TString::Format("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);

  return AliFemtoString(report);
}
//__________________

void AliFemtoV0PairCut::SetV0Max(Double_t aV0Max)
{
  fV0Max = aV0Max;
}

Double_t AliFemtoV0PairCut::GetAliFemtoV0Max() const
{
  return fV0Max;
}

TList *AliFemtoV0PairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  // The TString format patterns (F is float, I is integer, L is long)
  const char ptrnF[] = "AliFemtoV0PairCut.%s=%f",
             ptrnI[] = "AliFemtoV0PairCut.%s=%d",
             ptrnL[] = "AliFemtoV0PairCut.%s=%ld";

  tListSetttings->Add(new TObjString(TString::Format(ptrnF, "V0max", fV0Max)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnF, "sharefractionmax", fShareFractionMax)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnF, "TPCmin", fDTPCMin)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnF, "TPCexitmin", fDTPCExitMin)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnI, "datatype", fDataType)));

  tListSetttings->Add(new TObjString(TString::Format(ptrnF, "minAvgSepPosPos", fMinAvgSepPosPos)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnF, "minAvgSepPosNeg", fMinAvgSepPosNeg)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnF, "minAvgSepNegPos", fMinAvgSepNegPos)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnF, "minAvgSepNegNeg", fMinAvgSepNegNeg)));

  tListSetttings->Add(new TObjString(TString::Format(ptrnL, "pairs_passed", fNPairsPassed)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnL, "pairs_failed", fNPairsFailed)));

  return tListSetttings;
}

void AliFemtoV0PairCut::SetRemoveSameLabel(Bool_t aRemove)
{
  fRemoveSameLabel = aRemove;
}

void AliFemtoV0PairCut::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}

void AliFemtoV0PairCut::SetTPCEntranceSepMinimum(double dtpc)
{
  fDTPCMin = dtpc;
}

void AliFemtoV0PairCut::SetTPCExitSepMinimum(double dtpc)
{
  fDTPCExitMin = dtpc;
}

void AliFemtoV0PairCut::SetMinAvgSeparation(int type, double minSep)
{
  if (type == 0) //Pos-Pos
    fMinAvgSepPosPos = minSep;
  else if (type == 1) //Pos-Neg
    fMinAvgSepPosNeg = minSep;
  else if (type == 2) //Neg-Pos
    fMinAvgSepNegPos = minSep;
  else if (type == 3) //Neg-Neg
    fMinAvgSepNegNeg = minSep;
}
