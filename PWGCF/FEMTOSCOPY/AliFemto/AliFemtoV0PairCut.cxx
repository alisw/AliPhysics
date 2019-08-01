///
/// \file AliFemtoV0PairCut.cxx
///

#include "AliFemtoV0PairCut.h"

#include <string>
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoV0PairCut);
  /// \endcond
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
  fMinAvgSepNegNeg(0),
  fNanoAODAnalysis(false)
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
  fNanoAODAnalysis = cut.fNanoAODAnalysis;
  return *this;
}

//__________________
bool AliFemtoV0PairCut::Pass(const AliFemtoPair *pair)
{
  const AliFemtoV0 *V0_1 = pair->Track1()->V0(),
                   *V0_2 = pair->Track2()->V0();

  // Assert pair is of two V0 particles
  if (V0_1 == nullptr || V0_2 == nullptr) {
    return false;
  }

  // Assert that both particles do not share daughters
  if (V0_1->IdNeg() == V0_2->IdNeg() || V0_1->IdPos() == V0_2->IdPos()) {
    return false;
  }

  // Test separation between track daughters' entrance and exit points
  if(!fNanoAODAnalysis){
  if (fDataType == kESD || fDataType == kAOD) {
    const AliFemtoThreeVector diffPosEntrance = V0_1->NominalTpcEntrancePointPos() - V0_2->NominalTpcEntrancePointPos(),
                                  diffPosExit = V0_1->NominalTpcExitPointPos() - V0_2->NominalTpcExitPointPos(),

                              diffNegEntrance = V0_1->NominalTpcEntrancePointNeg() - V0_2->NominalTpcEntrancePointNeg(),
                                  diffNegExit = V0_1->NominalTpcExitPointNeg() - V0_2->NominalTpcExitPointNeg();

    if (diffPosEntrance.Mag() < fDTPCMin
        || diffNegEntrance.Mag() < fDTPCMin
        || diffPosExit.Mag() < fDTPCExitMin
        || diffNegExit.Mag() < fDTPCExitMin) {
      return false;
    }
  }
}
  // Find average separations between tracks
  double avgSep_pp = 0.0,
         avgSep_pn = 0.0,
         avgSep_np = 0.0,
         avgSep_nn = 0.0;

  AliFemtoPair::CalcAvgSepV0V0(*V0_1, *V0_2,
                               avgSep_nn,
                               avgSep_np,
                               avgSep_pn,
                               avgSep_pp);
if(!fNanoAODAnalysis){
  if (avgSep_pp < fMinAvgSepPosPos) {
    return false;
  }

  if (avgSep_pn < fMinAvgSepPosNeg) {
    return false;
  }

  if (avgSep_np < fMinAvgSepNegPos) {
    return false;
  }

  if (avgSep_nn < fMinAvgSepNegNeg) {
    return false;
  }
}

  return true;
}
//__________________
AliFemtoString AliFemtoV0PairCut::Report()
{
  AliFemtoString report = "AliFemtoV0 Pair Cut - remove shared and split pairs\n";
  report += Form("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);

  return report;
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
  TString prefix = "AliFemtoV0PairCut.";

  tListSetttings->AddVector(
    new TObjString(prefix + Form("V0max=%g", fV0Max)),
    new TObjString(prefix + Form("sharefractionmax=%g", fShareFractionMax)),
    new TObjString(prefix + Form("TPCmin=%g", fDTPCMin)),
    new TObjString(prefix + Form("TPCexitmin=%g", fDTPCExitMin)),
    new TObjString(prefix + Form("datatype=%d", fDataType)),

    new TObjString(prefix + Form("minAvgSepPosPos=%g", fMinAvgSepPosPos)),
    new TObjString(prefix + Form("minAvgSepPosNeg=%g", fMinAvgSepPosNeg)),
    new TObjString(prefix + Form("minAvgSepNegPos=%g", fMinAvgSepNegPos)),
    new TObjString(prefix + Form("minAvgSepNegNeg=%g", fMinAvgSepNegNeg)),

    new TObjString(prefix + Form("pairs_passed=%ld", fNPairsPassed)),
    new TObjString(prefix + Form("pairs_failed=%ld", fNPairsFailed)),

    nullptr);

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
  switch (type) {
  case 0: // Pos-Pos
    fMinAvgSepPosPos = minSep;
    break;
  case 1: // Pos-Neg
    fMinAvgSepPosNeg = minSep;
    break;
  case 2: // Neg-Pos
    fMinAvgSepNegPos = minSep;
    break;
  case 3: // Neg-Neg
    fMinAvgSepNegNeg = minSep;
    break;
  }
}
