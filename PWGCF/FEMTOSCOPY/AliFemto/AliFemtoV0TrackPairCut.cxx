///
/// \file AliFemtoV0TrackPairCut.cxx
///

#include "AliFemtoV0TrackPairCut.h"
// #include "AliFemtoAvgSepCalculator.h"

#include <string>
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoV0TrackPairCut);
  /// \endcond
#endif

//__________________
AliFemtoV0TrackPairCut::AliFemtoV0TrackPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fV0Max(1.0),
  fShareQualityMax(1.0),
  fShareFractionMax(1.0),
  fRemoveSameLabel(kFALSE),
  fTrackTPCOnly(kFALSE),
  fDataType(kAOD),
  fDTPCMin(0.0),
  fDTPCExitMin(0.0),
  fKstarCut(0.0),
  fFirstParticleType(kLambda),
  fSecondParticleType(kProton),
  fMinAvgSepTrackPos(0.0),
  fMinAvgSepTrackNeg(0.0),
  fMinDEtaStarPos(0.0),
  fMinDEtaStarNeg(0.0),
  fMinDPhiStarPos(0.0),
  fMinDPhiStarNeg(0.0),
  fMinRad(0.0),
  fNanoAODAnalysis(false)
{
  /* no-op */
}
//__________________
AliFemtoV0TrackPairCut::~AliFemtoV0TrackPairCut()
{
  /* no-op */
}
AliFemtoV0TrackPairCut &AliFemtoV0TrackPairCut::operator=(const AliFemtoV0TrackPairCut &cut)
{
  if (this == &cut) {
    return *this;
  }

  AliFemtoPairCut::operator=(cut);
  fNPairsPassed = cut.fNPairsPassed;
  fNPairsFailed = cut.fNPairsFailed;
  fV0Max = cut.fV0Max;
  fShareQualityMax = cut.fShareQualityMax;
  fShareFractionMax = cut.fShareFractionMax;
  fRemoveSameLabel = cut.fRemoveSameLabel;
  fTrackTPCOnly = cut.fTrackTPCOnly;
  fDataType = cut.fDataType;
  fDTPCMin = cut.fDTPCMin;
  fDTPCExitMin = cut.fDTPCExitMin;

  fMinAvgSepTrackPos = cut.fMinAvgSepTrackPos;
  fMinAvgSepTrackNeg = cut.fMinAvgSepTrackNeg;

  fMinDEtaStarPos = cut.fMinDEtaStarPos;
  fMinDEtaStarNeg = cut.fMinDEtaStarNeg;
  fMinDPhiStarPos = cut.fMinDPhiStarPos;
  fMinDPhiStarNeg = cut.fMinDPhiStarNeg;

  fMinRad = cut.fMinRad;
  fNanoAODAnalysis = cut.fNanoAODAnalysis;

  return *this;
}
//__________________
bool AliFemtoV0TrackPairCut::Pass(const AliFemtoPair *pair)
{
  //Track1 - V0
  //Track2 - track
  const AliFemtoV0 *V0 = pair->Track1()->V0();
  const AliFemtoTrack *track = pair->Track2()->Track();

  // check that we have a V0 and a track
  if (V0 == NULL || track == NULL) {
    fNPairsFailed++;
    return false;
  }

  // Check if track is a daughter of the V0
  const int pos_id = (fTrackTPCOnly) ? -(V0->IdPos() + 1) : V0->IdPos(),
            neg_id = (fTrackTPCOnly) ? -(V0->IdNeg() + 1) : V0->IdNeg();

  if (pos_id == track->TrackId() || neg_id == track->TrackId()) {
    fNPairsFailed++;
    return false;
  }

  //
  // Test separation between the track and the V0 daughters' entrance and exit points
  //
if(!fNanoAODAnalysis){
  if (fDataType == kESD || fDataType == kAOD) {
    const AliFemtoThreeVector diffPosEntrance = V0->NominalTpcEntrancePointPos() - track->NominalTpcEntrancePoint(),
                                  diffPosExit = V0->NominalTpcExitPointPos() - track->NominalTpcExitPoint(),

                              diffNegEntrance = V0->NominalTpcEntrancePointNeg() - track->NominalTpcEntrancePoint(),
                                  diffNegExit = V0->NominalTpcExitPointNeg() - track->NominalTpcExitPoint();

    if (diffPosEntrance.Mag() < fDTPCMin || diffPosExit.Mag() < fDTPCExitMin ||
        diffNegEntrance.Mag() < fDTPCMin || diffNegExit.Mag() < fDTPCExitMin) {
      fNPairsFailed++;
      return false;
    }
  }
}

  //
  // Check for pairs that are possibly shared/double reconstruction
  //
  if ((fShareQualityMax < 1.0) || (fShareFractionMax < 1.0)) {

    Int_t number_of_cluster_hits = 0,
          an = 0,
          number_of_shared_hits = 0;

    const TBits &track_clusters = track->TPCclusters(),
                  &pos_clusters = V0->TPCclustersPos(),
                  &neg_clusters = V0->TPCclustersNeg(),

                 &track_sharing = track->TPCsharing(),
                   &pos_sharing = V0->TPCsharingPos(),
                   &neg_sharing = V0->TPCsharingNeg();

    for (unsigned int imap = 0; imap < pos_clusters.GetNbits(); imap++) {

      // If both have clusters in the same row
      if (pos_clusters.TestBitNumber(imap) && track_clusters.TestBitNumber(imap)) {

        number_of_cluster_hits += 2;

        // Do they share the cluster?
        if (pos_sharing.TestBitNumber(imap) && track_sharing.TestBitNumber(imap)) {
          an++;
          number_of_shared_hits += 2;

        // Different hits on the same padrow
        } else {
          an--;
        }

      // One track has a hit, the other does not
      } else if (pos_clusters.TestBitNumber(imap) || track_clusters.TestBitNumber(imap)) {
        an++;
        number_of_cluster_hits++;
      }
    }

    Float_t hsmval = (number_of_cluster_hits > 0) ? float(an) / number_of_cluster_hits : 0.0,
            hsfval = (number_of_cluster_hits > 0) ? float(number_of_shared_hits) / number_of_cluster_hits : 0.0;

    if (((fShareQualityMax < 1.0) && (hsmval > fShareQualityMax)) || ((fShareFractionMax < 1.0) && (hsfval > fShareFractionMax))) {
      fNPairsFailed++;
      return false;
    }

    number_of_cluster_hits = 0;
    an = 0;
    number_of_shared_hits = 0;

    for (unsigned int imap = 0; imap < neg_clusters.GetNbits(); imap++) {
      // If both have clusters in the same row
      if (neg_clusters.TestBitNumber(imap) && track_clusters.TestBitNumber(imap)) {

        number_of_cluster_hits += 2;

        // Do they share the cluster?
        if (neg_sharing.TestBitNumber(imap) && track_sharing.TestBitNumber(imap)) {
          an++;
          number_of_shared_hits += 2;

        // Different hits on the same padrow
        } else {
          an--;
        }

      } else if (neg_clusters.TestBitNumber(imap) || track_clusters.TestBitNumber(imap)) {
        // One track has a hit, the other does not
        an++;
        number_of_cluster_hits++;
      }
    }

    hsmval = (number_of_cluster_hits > 0) ? float(an) / number_of_cluster_hits : 0.0;
    hsfval = (number_of_cluster_hits > 0) ? float(number_of_shared_hits) / number_of_cluster_hits : 0.0;

    if (((fShareQualityMax < 1.0) && (hsmval > fShareQualityMax)) || ((fShareFractionMax < 1.0) && (hsfval > fShareFractionMax))) {
      fNPairsFailed++;
      return false;
    }
  }

  //
  // Test the average separation between the track and each daughter in TPC
  //
if(!fNanoAODAnalysis){
    Double_t pos_avgSep = pair->NominalTpcAverageSeparationTrackV0Pos(),
             neg_avgSep = pair->NominalTpcAverageSeparationTrackV0Neg();

    if (pos_avgSep < fMinAvgSepTrackPos) {
      fNPairsFailed++;
      return false;
    }

    if (neg_avgSep < fMinAvgSepTrackNeg) {
      fNPairsFailed++;
      return false;
    }

  }

  //
  // Qinv cut (we are trying to get rid of antiresidual correlation between primary protons)
  //
  if (fKstarCut > 0.0) {

    // daughter... cut?
    const AliFemtoThreeVector pos_p = V0->MomPos(),
                              neg_p = V0->MomNeg(),
                            track_p = track->P();

    //double PionMass = 0.13956995;
    //double KaonMass = 0.493677;
    const double ProtonMass = 0.938272;
    //double LambdaMass = 1.115683;

    AliFemtoThreeVector temp1;
    double ener1 = 0.0;

    if (fFirstParticleType == kLambda && fSecondParticleType == kProton) {
      temp1 = pos_p;
      ener1 = ::sqrt(temp1.Mag2() + ProtonMass * ProtonMass);
    } else if (fFirstParticleType == kAntiLambda && fSecondParticleType == kAntiProton) {
      temp1 = neg_p;
      ener1 = ::sqrt(temp1.Mag2() + ProtonMass * ProtonMass);
    }

    AliFemtoThreeVector temp2 = track_p;
    double ener2 = 0;

    if (fSecondParticleType == kProton || fSecondParticleType == kAntiProton) {
      ener2 = ::sqrt(temp2.Mag2() + ProtonMass * ProtonMass);
    }

    // Particle momentum
    const AliFemtoLorentzVector fourMomentum1(ener1, temp1);
    const AliFemtoLorentzVector fourMomentum2(ener2, temp2);

    // Calculate qInv
    AliFemtoLorentzVector tDiff = (fourMomentum1 - fourMomentum2);

    float tQinv = fabs(tDiff.m());   // note - qInv() will be negative for identical pairs...

    if (tQinv / 2.0 < fKstarCut) {
      fNPairsFailed++;
      return false;
    }
  }

  //
  // Delta Eta* Delta Phi*cut
  //
  if (fMinRad > 0.0) {

    const auto v0_shifted_tpc_point_pos = V0->NominalTpcPointPosShifted(),
               v0_shifted_tpc_point_neg = V0->NominalTpcPointNegShifted(),
               track_shifted_tpc_point = track->NominalTpcPointShifted();

    const auto diff_point_pos = v0_shifted_tpc_point_pos - track_shifted_tpc_point,
               diff_point_neg = v0_shifted_tpc_point_neg - track_shifted_tpc_point;

    const double thetas1_pos = TMath::Pi()/2. - TMath::ATan(v0_shifted_tpc_point_pos.z()/(fMinRad*1e2)),
                 thetas1_neg = TMath::Pi()/2. - TMath::ATan(v0_shifted_tpc_point_neg.z()/(fMinRad*1e2)),
                     thetas2 = TMath::Pi()/2. - TMath::ATan(track_shifted_tpc_point.z()/(fMinRad*1e2));

    const double etas1_pos = -TMath::Log( TMath::Tan(thetas1_pos/2.) ),
                 etas1_neg = -TMath::Log( TMath::Tan(thetas1_neg/2.) ),
                     etas2 = -TMath::Log( TMath::Tan(thetas2/2.) );

    const double detas_pos = etas1_pos - etas2,
                 detas_neg = etas1_neg - etas2;

    const double distSft_pos = diff_point_pos.Perp(),
                 distSft_neg = diff_point_neg.Perp();

    const double dPhiS_pos = 2.0 * TMath::ATan(distSft_pos/2./((fMinRad*1e2))),
                 dPhiS_neg = 2.0 * TMath::ATan(distSft_neg/2./((fMinRad*1e2)));

    if ( (TMath::Abs(detas_pos) < fMinDEtaStarPos &&
          TMath::Abs(dPhiS_pos) < fMinDPhiStarPos) ||
         (TMath::Abs(detas_neg) < fMinDEtaStarNeg &&
          TMath::Abs(dPhiS_neg) < fMinDPhiStarNeg) ) {
      fNPairsFailed++;
      return false;
    }
  }

  fNPairsPassed++;

  return true;
}
//__________________
AliFemtoString AliFemtoV0TrackPairCut::Report()
{
  TString report = "AliFemtoV0 Pair Cut - remove shared and split pairs\n";
  report += TString::Format("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);

  return AliFemtoString((const char *)report);
}
//__________________

void AliFemtoV0TrackPairCut::SetV0Max(Double_t aV0Max)
{
  fV0Max = aV0Max;
}

Double_t AliFemtoV0TrackPairCut::GetAliFemtoV0Max() const
{
  return fV0Max;
}


TList *AliFemtoV0TrackPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *list = new TList();

  // The TString format patterns (F is float, I is integer, L is long)
  const TString prefix = "AliFemtoV0TrackPairCut.";

  list->Add(new TObjString(prefix + Form("sharequalitymax=%g", fShareQualityMax)));
  list->Add(new TObjString(prefix + Form("sharefractionmax=%g", fShareFractionMax)));

  list->Add(new TObjString(prefix + Form("pairs_passed=%ld", fNPairsPassed)));
  list->Add(new TObjString(prefix + Form("pairs_failed=%ld", fNPairsFailed)));

  return list;
}

void AliFemtoV0TrackPairCut::SetRemoveSameLabel(Bool_t aRemove)
{
  fRemoveSameLabel = aRemove;
}

void AliFemtoV0TrackPairCut::SetTPCOnly(Bool_t tpconly)
{
  fTrackTPCOnly = tpconly;
}

void AliFemtoV0TrackPairCut::SetShareQualityMax(Double_t aShareQualityMax)
{
  fShareQualityMax = aShareQualityMax;
}

void AliFemtoV0TrackPairCut::SetShareFractionMax(Double_t aShareFractionMax)
{
  fShareFractionMax = aShareFractionMax;
}

void AliFemtoV0TrackPairCut::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}

void AliFemtoV0TrackPairCut::SetTPCEntranceSepMinimum(double dtpc)
{
  fDTPCMin = dtpc;
}

void AliFemtoV0TrackPairCut::SetTPCExitSepMinimum(double dtpc)
{
  fDTPCExitMin = dtpc;
}

void AliFemtoV0TrackPairCut::SetKstarCut(double kstar, AliFemtoParticleType firstParticle, AliFemtoParticleType secondParticle)
{
  fKstarCut = kstar;
  fFirstParticleType = firstParticle;  //for kstar - first particle type (V0 type)
  fSecondParticleType = secondParticle;
}

void AliFemtoV0TrackPairCut::SetMinAvgSeparation(int type, double minSep)
{
  if (type == 0) //Track-Pos
    fMinAvgSepTrackPos = minSep;
  else if (type == 1) //Track-Neg
    fMinAvgSepTrackNeg = minSep;
}

void AliFemtoV0TrackPairCut::SetMinDEtaStar(Int_t type, Double_t minDEtaStar)
{
  if (type == 0) //Track-Pos
    fMinDEtaStarPos = minDEtaStar;
  else if (type == 1) //Track-Neg
    fMinDEtaStarNeg = minDEtaStar;
}

void AliFemtoV0TrackPairCut::SetMinDPhiStar(Int_t type, Double_t minDPhiStar)
{
  if (type == 0) //Track-Pos
    fMinDPhiStarPos = minDPhiStar;
  else if (type == 1) //Track-Neg
    fMinDPhiStarNeg = minDPhiStar;
}

void AliFemtoV0TrackPairCut::SetShiftPosition(Double_t rad)
{
  fMinRad = rad;
}
