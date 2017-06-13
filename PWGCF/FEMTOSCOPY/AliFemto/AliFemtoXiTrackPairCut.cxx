#include "AliFemtoXiTrackPairCut.h"

#include <string>
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiTrackPairCut);
  /// \endcond
#endif

//__________________
AliFemtoXiTrackPairCut::AliFemtoXiTrackPairCut():
  fV0TrackPairCut(nullptr),
  fNPairsPassed(0),
  fNPairsFailed(0),
  fMinAvgSepTrackBacPion(0.),
  fDataType(kAOD),
  fTrackTPCOnly(0)
{
  fV0TrackPairCut = new AliFemtoV0TrackPairCut();
}
//__________________
AliFemtoXiTrackPairCut::~AliFemtoXiTrackPairCut()
{
  delete fV0TrackPairCut;
}

AliFemtoXiTrackPairCut &AliFemtoXiTrackPairCut::operator=(const AliFemtoXiTrackPairCut &cut)
{
  if (this == &cut) {
    return *this;
  }

  AliFemtoPairCut::operator=(cut);
  fNPairsPassed = cut.fNPairsPassed;
  fNPairsFailed = cut.fNPairsFailed;
  fMinAvgSepTrackBacPion = cut.fMinAvgSepTrackBacPion;
  fDataType = cut.fDataType;
  fTrackTPCOnly = cut.fTrackTPCOnly;

  if(fV0TrackPairCut) delete fV0TrackPairCut;
  fV0TrackPairCut = new AliFemtoV0TrackPairCut(*cut.fV0TrackPairCut);

  return *this; 
}

//__________________
bool AliFemtoXiTrackPairCut::Pass(const AliFemtoPair *pair) 
{
  //Track1 - Xi
  //Track2 - track
  const AliFemtoXi *Xi = pair->Track1()->Xi();
  const AliFemtoTrack *track = pair->Track2()->Track();

  // check that we have a Xi and a track
  if (Xi == NULL || track == NULL) {
    fNPairsFailed++;
    return false;
  }

  // Check if track is a daughter of the Xi or the Xi bachelor
  const int pos_id = (fTrackTPCOnly) ? -(Xi->IdPos() + 1) : Xi->IdPos(),
            neg_id = (fTrackTPCOnly) ? -(Xi->IdNeg() + 1) : Xi->IdNeg(),
            bac_id = (fTrackTPCOnly) ? -(Xi->IdBac() + 1) : Xi->IdBac();

  if (pos_id == track->TrackId() || neg_id == track->TrackId() || bac_id == track->TrackId()) {
    fNPairsFailed++;
    return false;
  }


  if (Xi->IdBac() == Xi->IdPos() || Xi->IdBac() == Xi->IdNeg()) {
    return false;
  }

  //Calling tPassV0 = fV0TrackPairCut->Pass(tPair) below will handle the average separation of
  //the V0 daughters to the track.  So, all that needs to be checked here in the average separation
  //of the bachelor pion to the track
  UInt_t bac_point_cnt = 0;
  Double_t bac_avgSep = 0.0;

  // loop through NominalTpcPoints of the track and bachelor pion
  for (int i = 0; i < 8; i++) {
    // Grab references to each of the i'th points
    const AliFemtoThreeVector &bac_p = Xi->NominalTpcPointBac(i),
                            &track_p = track->NominalTpcPoint(i);

    // if any track points are outside the boundary - skip
    if (track_p.x() < -9990.0 || track_p.y() < -9990.0 || track_p.z() < -9990.0) {
      continue;
    }

    // If the bachelor pion points are not bad, increment point count and
    // increase the cumulative average separation
    if (!(bac_p.x() < -9990.0 || bac_p.y() < -9990.0 || bac_p.z() < -9990.0)) {
      bac_avgSep += (bac_p - track_p).Mag();
      bac_point_cnt++;
    }
  }

  if (bac_point_cnt == 0 || bac_avgSep / bac_point_cnt < fMinAvgSepTrackBacPion) {
    fNPairsFailed++;
    return false;
  }




  //Make sure it passes AliFemtoV0TrackPairCuts
  double tLambdaMass = 1.115683;
  double tPionMass = 0.19357018;
  AliFemtoPair *tPair = new AliFemtoPair();
  AliFemtoParticle* tV0 = new AliFemtoParticle((AliFemtoV0*)Xi, tLambdaMass);
  AliFemtoParticle* tTrack = new AliFemtoParticle(track, tPionMass);
  tPair->SetTrack1(tV0);
  tPair->SetTrack2(tTrack);
  bool tPassV0 = false;
  tPassV0 = fV0TrackPairCut->Pass(tPair);
  delete tPair;
  delete tV0;
  delete tTrack;
  if(!tPassV0) return false;

  return true;
}
//__________________
AliFemtoString AliFemtoXiTrackPairCut::Report()
{
  TString report = "AliFemtoXi Pair Cut - remove shared and split pairs\n";
  report += TString::Format("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);

  return AliFemtoString(report);
}
//__________________



TList *AliFemtoXiTrackPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  // The TString format patterns (F is float, I is integer, L is long)
  const char ptrnF[] = "AliFemtoXiTrackPairCut.%s=%f",
             ptrnI[] = "AliFemtoXiTrackPairCut.%s=%d",
             ptrnL[] = "AliFemtoXiTrackPairCut.%s=%ld";


  tListSetttings->Add(new TObjString(TString::Format(ptrnI, "datatype", fDataType)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnL, "pairs_passed", fNPairsPassed)));
  tListSetttings->Add(new TObjString(TString::Format(ptrnL, "pairs_failed", fNPairsFailed)));

  return tListSetttings;
}


void AliFemtoXiTrackPairCut::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}

void AliFemtoXiTrackPairCut::SetTPCOnly(Bool_t tpconly)
{
  fTrackTPCOnly = tpconly;
}
