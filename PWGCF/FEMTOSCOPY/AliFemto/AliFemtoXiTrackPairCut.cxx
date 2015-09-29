#include "AliFemtoXiTrackPairCut.h"

#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoXiTrackPairCut)
#endif

//__________________
AliFemtoXiTrackPairCut::AliFemtoXiTrackPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fDataType(kAOD),
  fTrackTPCOnly(0)
{
  /* no-op */
}
//__________________
AliFemtoXiTrackPairCut::~AliFemtoXiTrackPairCut()
{
  /* no-op */
}

AliFemtoXiTrackPairCut &AliFemtoXiTrackPairCut::operator=(const AliFemtoXiTrackPairCut &cut)
{
  if (this == &cut) {
    return *this;
  }

  AliFemtoPairCut::operator=(cut);
  fNPairsPassed = cut.fNPairsPassed;
  fNPairsFailed = cut.fNPairsFailed;
  fDataType = cut.fDataType;
  fTrackTPCOnly = cut.fTrackTPCOnly;

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
