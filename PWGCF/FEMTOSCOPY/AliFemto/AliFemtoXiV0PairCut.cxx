#include "AliFemtoXiV0PairCut.h"

#include <string>
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiV0PairCut);
  /// \endcond
#endif

//__________________
AliFemtoXiV0PairCut::AliFemtoXiV0PairCut():
  fV0PairCut(nullptr),
  fNPairsPassed(0),
  fNPairsFailed(0),
  fDataType(kAOD),
  fMinAvgSepBacPos(0.),
  fMinAvgSepBacNeg(0.)
{
  fV0PairCut = new AliFemtoV0PairCut();
}
//__________________
AliFemtoXiV0PairCut::~AliFemtoXiV0PairCut()
{
  delete fV0PairCut;
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

  fMinAvgSepBacPos = cut.fMinAvgSepBacPos;
  fMinAvgSepBacNeg = cut.fMinAvgSepBacNeg;

  if(fV0PairCut) delete fV0PairCut;
  fV0PairCut = new AliFemtoV0PairCut(*cut.fV0PairCut);

  return *this;
}

//__________________
inline
bool is_unset_vector(const AliFemtoThreeVector& v)
{
  return v.x() <= -9999.0 || v.y() <= -9999.0 ||  v.z() <= -9999.0;
}

//__________________
inline
double calculate_avg_separation_BacV0Daughter(const AliFemtoXi* aXi,
                                              const AliFemtoV0* aV0,
                                              const bool use_pos_daughter)
{
  int counter = 0;
  double avgSep = 0.0;

  for (int i = 0; i < 8; i++) {
    const AliFemtoThreeVector &bac_p = aXi->NominalTpcPointBac(i);
    const AliFemtoThreeVector &V0Daught_p = use_pos_daughter
      ? aV0->NominalTpcPointPos(i)
      : aV0->NominalTpcPointNeg(i);


    if (is_unset_vector(bac_p) || is_unset_vector(V0Daught_p)) {
      continue;
    }
    avgSep += (bac_p - V0Daught_p).Mag();
    counter++;
  }
  if(counter==0) return 0;
  return avgSep / counter;
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

  //Calling tPassV0 = fV0PairCut->Pass(tPair) below will handle the average separation of
  //the Xi's V0 daughters to the V0 daughters.  So, all that needs to be checked here in the average separation
  //of the bachelor pion to the V0 daughters


  // Find average separations between bachelor pion and V0 daughters
  double avgSep = 0.0;

  avgSep = calculate_avg_separation_BacV0Daughter(Xi, V0, true);
  if (avgSep < fMinAvgSepBacPos) {
    return false;
  }

  avgSep = calculate_avg_separation_BacV0Daughter(Xi, V0, false);
  if (avgSep < fMinAvgSepBacNeg) {
    return false;
  }

  //Make sure it passes AliFemtoV0PairCut
  double tLambdaMass = 1.115683;
  double tK0ShortMass = 0.497611;
  AliFemtoPair *tPair = new AliFemtoPair();
  AliFemtoParticle* tXiV0 = new AliFemtoParticle((AliFemtoV0*)Xi, tLambdaMass);
  AliFemtoParticle* tV0 = new AliFemtoParticle(V0, tK0ShortMass);
  tPair->SetTrack1(tXiV0);
  tPair->SetTrack2(tV0);
  bool tPassV0 = false;
  tPassV0 = fV0PairCut->Pass(tPair);
  delete tPair;
  delete tXiV0;
  delete tV0;
  if(!tPassV0) return false;

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
