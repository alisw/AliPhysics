///
/// \file AliFemtoAvgSepCorrFctn.cxx
///

#include "AliFemtoAvgSepCorrFctn.h"
//#include "AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoAvgSepCorrFctn)
#endif

//____________________________
AliFemtoAvgSepCorrFctn::AliFemtoAvgSepCorrFctn(const char *title,
                                               const int &nbins,
                                               const float &Low,
                                               const float &High):
  fNumerator(0), //2 tracks
  fDenominator(0),
  fNumeratorPos(0), //track + V0
  fDenominatorPos(0),
  fNumeratorNeg(0),
  fDenominatorNeg(0),
  fNumeratorPosPos(0), //2 V0s
  fDenominatorPosPos(0),
  fNumeratorPosNeg(0),
  fDenominatorPosNeg(0),
  fNumeratorNegPos(0),
  fDenominatorNegPos(0),
  fNumeratorNegNeg(0),
  fDenominatorNegNeg(0),
  fNumeratorBacTrack(0),
  fDenominatorBacTrack(0),
  fNumeratorBacPos(0),
  fDenominatorBacPos(0),
  fNumeratorBacNeg(0),
  fDenominatorBacNeg(0),
  fRatio(0),
  fPairType(kTracks)
{
  // set up numerator
  //  title = "Num Qinv (MeV/c)";
  char tTitNum[101] = "Num";
  strncat(tTitNum, title, 100);
  fNumerator = new TH1D(tTitNum, title, nbins, Low, High);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTitDen[101] = "Den";
  strncat(tTitDen, title, 100);
  fDenominator = new TH1D(tTitDen, title, nbins, Low, High);
  // set up ratio
  //title = "Ratio Qinv (MeV/c)";
  char tTitRat[101] = "Rat";
  strncat(tTitRat, title, 100);
  fRatio = new TH1D(tTitRat, title, nbins, Low, High);


  char tTitNumPos[101] = "NumV0TrackPos";
  strncat(tTitNumPos, title, 100);
  fNumeratorPos = new TH1D(tTitNumPos, title, nbins, Low, High);
  char tTitDenPos[101] = "DenV0TrackPos";
  strncat(tTitDenPos, title, 100);
  fDenominatorPos = new TH1D(tTitDenPos, title, nbins, Low, High);
  char tTitNumNeg[101] = "NumV0TrackNeg";
  strncat(tTitNumNeg, title, 100);
  fNumeratorNeg = new TH1D(tTitNumNeg, title, nbins, Low, High);
  char tTitDenNeg[101] = "DenV0TrackNeg";
  strncat(tTitDenNeg, title, 100);
  fDenominatorNeg = new TH1D(tTitDenNeg, title, nbins, Low, High);


  char tTitNumPosPos[101] = "NumV0sPosPos";
  strncat(tTitNumPosPos, title, 100);
  fNumeratorPosPos = new TH1D(tTitNumPosPos, title, nbins, Low, High);
  char tTitDenPosPos[101] = "DenV0sPosPos";
  strncat(tTitDenPosPos, title, 100);
  fDenominatorPosPos = new TH1D(tTitDenPosPos, title, nbins, Low, High);
  char tTitNumPosNeg[101] = "NumV0sPosNeg";
  strncat(tTitNumPosNeg, title, 100);
  fNumeratorPosNeg = new TH1D(tTitNumPosNeg, title, nbins, Low, High);
  char tTitDenPosNeg[101] = "DenV0sPosNeg";
  strncat(tTitDenPosNeg, title, 100);
  fDenominatorPosNeg = new TH1D(tTitDenPosNeg, title, nbins, Low, High);
  char tTitNumNegPos[101] = "NumV0sNegPos";
  strncat(tTitNumNegPos, title, 100);
  fNumeratorNegPos = new TH1D(tTitNumNegPos, title, nbins, Low, High);
  char tTitDenNegPos[101] = "DenV0sNegPos";
  strncat(tTitDenNegPos, title, 100);
  fDenominatorNegPos = new TH1D(tTitDenNegPos, title, nbins, Low, High);
  char tTitNumNegNeg[101] = "NumV0sNegNeg";
  strncat(tTitNumNegNeg, title, 100);
  fNumeratorNegNeg = new TH1D(tTitNumNegNeg, title, nbins, Low, High);
  char tTitDenNegNeg[101] = "DenV0sNegNeg";
  strncat(tTitDenNegNeg, title, 100);
  fDenominatorNegNeg = new TH1D(tTitDenNegNeg, title, nbins, Low, High);

  char tTitNumBac[101] = "NumXiTrackBac";
  strncat(tTitNumBac, title, 100);
  fNumeratorBacTrack = new TH1D(tTitNumBac, title, nbins, Low, High);
  char tTitDenBac[101] = "DenXiTrackBac";
  strncat(tTitDenBac, title, 100);
  fDenominatorBacTrack = new TH1D(tTitDenBac, title, nbins, Low, High);

  char tTitNumBacPos[101] = "NumXiBacV0Pos";
  strncat(tTitNumBacPos, title, 100);
  fNumeratorBacPos = new TH1D(tTitNumBacPos, title, nbins, Low, High);
  char tTitDenBacPos[101] = "DenXiBacV0Pos";
  strncat(tTitDenBacPos, title, 100);
  fDenominatorBacPos = new TH1D(tTitDenBacPos, title, nbins, Low, High);

  char tTitNumBacNeg[101] = "NumXiBacV0Neg";
  strncat(tTitNumBacNeg, title, 100);
  fNumeratorBacNeg = new TH1D(tTitNumBacNeg, title, nbins, Low, High);
  char tTitDenBacNeg[101] = "DenXiBacV0Neg";
  strncat(tTitDenBacNeg, title, 100);
  fDenominatorBacNeg = new TH1D(tTitDenBacNeg, title, nbins, Low, High);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();

  fNumeratorPos->Sumw2();
  fDenominatorPos->Sumw2();
  fNumeratorNeg->Sumw2();
  fDenominatorNeg->Sumw2();

  fNumeratorPosPos->Sumw2();
  fDenominatorPosPos->Sumw2();
  fNumeratorPosNeg->Sumw2();
  fDenominatorPosNeg->Sumw2();
  fNumeratorNegPos->Sumw2();
  fDenominatorNegPos->Sumw2();
  fNumeratorNegNeg->Sumw2();
  fDenominatorNegNeg->Sumw2();

  fNumeratorBacTrack->Sumw2();
  fDenominatorBacTrack->Sumw2();

  fNumeratorBacPos->Sumw2();
  fDenominatorBacPos->Sumw2();
  fNumeratorBacNeg->Sumw2();
  fDenominatorBacNeg->Sumw2();

  fRatio->Sumw2();
}

//____________________________
AliFemtoAvgSepCorrFctn::AliFemtoAvgSepCorrFctn(const AliFemtoAvgSepCorrFctn &aCorrFctn) :
  AliFemtoCorrFctn(),
  fNumerator(0),
  fDenominator(0),
  fNumeratorPos(0),
  fDenominatorPos(0),
  fNumeratorNeg(0),
  fDenominatorNeg(0),
  fNumeratorPosPos(0),
  fDenominatorPosPos(0),
  fNumeratorPosNeg(0),
  fDenominatorPosNeg(0),
  fNumeratorNegPos(0),
  fDenominatorNegPos(0),
  fNumeratorNegNeg(0),
  fDenominatorNegNeg(0),
  fNumeratorBacTrack(0),
  fDenominatorBacTrack(0),
  fNumeratorBacPos(0),
  fDenominatorBacPos(0),
  fNumeratorBacNeg(0),
  fDenominatorBacNeg(0),
  fRatio(0),
  fPairType(kTracks)
{
  // copy constructor
  fNumerator = new TH1D(*aCorrFctn.fNumerator);
  fDenominator = new TH1D(*aCorrFctn.fDenominator);

  fNumeratorPos = new TH1D(*aCorrFctn.fNumeratorPos);
  fDenominatorPos = new TH1D(*aCorrFctn.fDenominatorPos);
  fNumeratorNeg = new TH1D(*aCorrFctn.fNumeratorNeg);
  fDenominatorNeg = new TH1D(*aCorrFctn.fDenominatorNeg);

  fNumeratorPosPos = new TH1D(*aCorrFctn.fNumeratorPosPos);
  fDenominatorPosPos = new TH1D(*aCorrFctn.fDenominatorPosPos);
  fNumeratorPosNeg = new TH1D(*aCorrFctn.fNumeratorPosNeg);
  fDenominatorPosNeg = new TH1D(*aCorrFctn.fDenominatorPosNeg);
  fNumeratorNegPos = new TH1D(*aCorrFctn.fNumeratorNegPos);
  fDenominatorNegPos = new TH1D(*aCorrFctn.fDenominatorNegPos);
  fNumeratorNegNeg = new TH1D(*aCorrFctn.fNumeratorNegNeg);
  fDenominatorNegNeg = new TH1D(*aCorrFctn.fDenominatorNegNeg);

  fNumeratorBacTrack = new TH1D(*aCorrFctn.fNumeratorBacTrack);
  fDenominatorBacTrack = new TH1D(*aCorrFctn.fDenominatorBacTrack);

  fNumeratorBacPos = new TH1D(*aCorrFctn.fNumeratorBacPos);
  fDenominatorBacPos = new TH1D(*aCorrFctn.fDenominatorBacPos);
  fNumeratorBacNeg = new TH1D(*aCorrFctn.fNumeratorBacNeg);
  fDenominatorBacNeg = new TH1D(*aCorrFctn.fDenominatorBacNeg);

  fRatio = new TH1D(*aCorrFctn.fRatio);
}
//____________________________
AliFemtoAvgSepCorrFctn::~AliFemtoAvgSepCorrFctn()
{
  // destructor
  delete fNumerator;
  delete fDenominator;

  delete fNumeratorPos;
  delete fDenominatorPos;
  delete fNumeratorNeg;
  delete fDenominatorNeg;

  delete fNumeratorPosPos;
  delete fDenominatorPosPos;
  delete fNumeratorPosNeg;
  delete fDenominatorPosNeg;
  delete fNumeratorNegPos;
  delete fDenominatorNegPos;
  delete fNumeratorNegNeg;
  delete fDenominatorNegNeg;

  delete fNumeratorBacTrack;
  delete fDenominatorBacTrack;

  delete fNumeratorBacPos;
  delete fDenominatorBacPos;
  delete fNumeratorBacNeg;
  delete fDenominatorBacNeg;

  delete fRatio;
}
//_________________________
AliFemtoAvgSepCorrFctn &AliFemtoAvgSepCorrFctn::operator=(const AliFemtoAvgSepCorrFctn &aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (fNumerator) delete fNumerator;
  fNumerator = new TH1D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH1D(*aCorrFctn.fDenominator);

  if (fNumeratorPos) delete fNumeratorPos;
  fNumeratorPos = new TH1D(*aCorrFctn.fNumeratorPos);
  if (fDenominatorPos) delete fDenominatorPos;
  fDenominatorPos = new TH1D(*aCorrFctn.fDenominatorPos);
  if (fNumeratorNeg) delete fNumeratorNeg;
  fNumeratorNeg = new TH1D(*aCorrFctn.fNumeratorNeg);
  if (fDenominatorNeg) delete fDenominatorNeg;
  fDenominatorNeg = new TH1D(*aCorrFctn.fDenominatorNeg);

  if (fNumeratorPosPos) delete fNumeratorPosPos;
  fNumeratorPosPos = new TH1D(*aCorrFctn.fNumeratorPosPos);
  if (fDenominatorPosPos) delete fDenominatorPosPos;
  fDenominatorPosPos = new TH1D(*aCorrFctn.fDenominatorPosPos);
  if (fNumeratorPosNeg) delete fNumeratorPosNeg;
  fNumeratorPosNeg = new TH1D(*aCorrFctn.fNumeratorPosNeg);
  if (fDenominatorPosNeg) delete fDenominatorPosNeg;
  fDenominatorPosNeg = new TH1D(*aCorrFctn.fDenominatorPosNeg);
  if (fNumeratorNegPos) delete fNumeratorNegPos;
  fNumeratorNegPos = new TH1D(*aCorrFctn.fNumeratorNegPos);
  if (fDenominatorNegPos) delete fDenominatorNegPos;
  fDenominatorNegPos = new TH1D(*aCorrFctn.fDenominatorNegPos);
  if (fNumeratorNegNeg) delete fNumeratorNegNeg;
  fNumeratorNegNeg = new TH1D(*aCorrFctn.fNumeratorNegNeg);
  if (fDenominatorNegNeg) delete fDenominatorNegNeg;
  fDenominatorNegNeg = new TH1D(*aCorrFctn.fDenominatorNegNeg);

  if (fNumeratorBacTrack) delete fNumeratorBacTrack;
  fNumeratorBacTrack = new TH1D(*aCorrFctn.fNumeratorBacTrack);
  if (fDenominatorBacTrack) delete fDenominatorBacTrack;
  fDenominatorBacTrack = new TH1D(*aCorrFctn.fDenominatorBacTrack);

  if (fNumeratorBacPos) delete fNumeratorBacPos;
  fNumeratorBacPos = new TH1D(*aCorrFctn.fNumeratorBacPos);
  if (fDenominatorBacPos) delete fDenominatorBacPos;
  fDenominatorBacPos = new TH1D(*aCorrFctn.fDenominatorBacPos);
  if (fNumeratorBacNeg) delete fNumeratorBacNeg;
  fNumeratorBacNeg = new TH1D(*aCorrFctn.fNumeratorBacNeg);
  if (fDenominatorBacNeg) delete fDenominatorBacNeg;
  fDenominatorBacNeg = new TH1D(*aCorrFctn.fDenominatorBacNeg);

  if (fRatio) delete fRatio;
  fRatio = new TH1D(*aCorrFctn.fRatio);

  return *this;
}

//_________________________
void AliFemtoAvgSepCorrFctn::Finish()
{
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  fNumerator->Draw();
  //fDenominator->Draw();
  //fRatio->Draw();
  fRatio->Divide(fNumerator, fDenominator, 1.0, 1.0);

}

//____________________________
AliFemtoString AliFemtoAvgSepCorrFctn::Report()
{
  // construct report
  string stemp = "Qinv Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n", fDenominator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in ratio:\t%E\n", fRatio->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}

// Function to check if the vector has a point near the default value (-9999)
// (i.e. it was never set)
static bool TpcPointIsUnset(const AliFemtoThreeVector& v) {
  return v.x() < -9000. ||
         v.y() < -9000. ||
         v.z() < -9000.;
}

static void StoreAvgSepBetweenTracks(const AliFemtoTrack *track_1,
                                     const AliFemtoTrack *track_2,
                                     TH1D *output)
{
  // sums the separation magnitude
  double avgSep = 0.0;
  int count = 0;

  // loop through the 8 points of the 'NominalTpcPoint' methods
  for (int i = 0; i < 8; i++) {
    const AliFemtoThreeVector &point_1 = track_1->NominalTpcPoint(i),
                              &point_2 = track_2->NominalTpcPoint(i);

    if (TpcPointIsUnset(point_1) || TpcPointIsUnset(point_2)) {
      break;
    }

    avgSep += (point_1 - point_2).Mag();
    count++;
  }
  if (count != 0) {
    output->Fill(avgSep / count);
  }
}

static void StoreAvgSepBetweenV0AndTrack(const AliFemtoV0 *V0,
                                         const AliFemtoTrack *track,
                                         TH1D *pos_output,
                                         TH1D *neg_output)
{
  // store number of successful points for pos an neg daughters
  int countPos = 0,
      countNeg = 0;

  // sums of the separation magnitude
  double avgSepPos = 0.0,
         avgSepNeg = 0.0;

  // Need non-const V0 for some reason.
  AliFemtoV0 *mutable_V0 = const_cast<AliFemtoV0*>(V0);

  // loop through the 8 points of the 'NominalTpcPoint' methods
  for (int i = 0; i < 8; i++) {
    const AliFemtoThreeVector &pos_vec = mutable_V0->NominalTpcPointPos(i),
                              &neg_vec = mutable_V0->NominalTpcPointNeg(i),
                              &trac_vec = track->NominalTpcPoint(i);

    if (TpcPointIsUnset(trac_vec)) {
      break;
    }

    const bool bad_pos = TpcPointIsUnset(pos_vec),
               bad_neg = TpcPointIsUnset(neg_vec);

    if (bad_pos && bad_neg) {
      break;
    }

    if (!bad_pos) {
      avgSepPos += (pos_vec - trac_vec).Mag();
      countPos++;
    }

    if (!bad_neg) {
      avgSepNeg += (neg_vec - trac_vec).Mag();
      countNeg++;
    }
  }
  if (countPos != 0) {
    pos_output->Fill(avgSepPos / countPos);
  }

  if (countNeg != 0) {
    neg_output->Fill(avgSepNeg / countNeg);
  }
}

static void StoreAvgSepBetweenV0s(const AliFemtoV0 *V0_1,
                                  const AliFemtoV0 *V0_2,
                                  TH1D *pospos_output,
                                  TH1D *posneg_output,
                                  TH1D *negpos_output,
                                  TH1D *negneg_output)
{
  // keep track of each combination of pos+neg tracks that are "good"
  int countPosPos = 0,
      countPosNeg = 0,
      countNegPos = 0,
      countNegNeg = 0;

  // sums of the separation magnitude
  double avgSepPosPos = 0.0,
         avgSepPosNeg = 0.0,
         avgSepNegPos = 0.0,
         avgSepNegNeg = 0.0;

  // Getting the TPC points requires a non-const AliFemtoV0, so we const_cast
  // for now, until this changes.
  AliFemtoV0 *mutable_v0_1 = const_cast<AliFemtoV0*>(V0_1),
             *mutable_v0_2 = const_cast<AliFemtoV0*>(V0_2);

  for (int i = 0; i < 8; i++) {

    const AliFemtoThreeVector &pos_1 = mutable_v0_1->NominalTpcPointPos(i),
                              &neg_1 = mutable_v0_1->NominalTpcPointNeg(i),
                              &pos_2 = mutable_v0_2->NominalTpcPointPos(i),
                              &neg_2 = mutable_v0_2->NominalTpcPointNeg(i);

    const bool bad_pos_1 = TpcPointIsUnset(pos_1),
               bad_neg_1 = TpcPointIsUnset(neg_1),
               bad_pos_2 = TpcPointIsUnset(pos_2),
               bad_neg_2 = TpcPointIsUnset(neg_2);

    if ((bad_pos_1 && bad_neg_1) || (bad_pos_2 && bad_neg_2)) {
      break;
    }

    if (!bad_pos_1 && !bad_pos_2) {
      avgSepPosPos += (pos_1 - pos_2).Mag();
      countPosPos++;
    }

    if (!bad_pos_1 && !bad_neg_2) {
      avgSepPosNeg += (pos_1 - neg_2).Mag();
      countPosNeg++;
    }

    if (!bad_neg_1 && !bad_pos_2) {
      avgSepNegPos += (neg_1 - pos_2).Mag();
      countNegPos++;
    }

    if (!bad_neg_1 && !bad_neg_2) {
      avgSepNegNeg += (neg_1 - neg_2).Mag();
      countNegNeg++;
    }
  }

  if (countPosPos != 0) {
    pospos_output->Fill(avgSepPosPos / countPosPos);
  }
  if (countPosNeg != 0) {
    posneg_output->Fill(avgSepPosNeg / countPosNeg);
  }
  if (countNegPos != 0) {
    negpos_output->Fill(avgSepNegPos / countNegPos);
  }
  if (countNegNeg != 0) {
    negneg_output->Fill(avgSepNegNeg / countNegNeg);
  }
}

static void StoreAvgSepBetweenXiBacAndTrack(const AliFemtoXi *Xi,
                                            const AliFemtoTrack *track,
                                            TH1D *bac_output)
{
  // store number of successful points for bachelor pion
  int countBac = 0;

  // sums of the separation magnitude
  double avgSepBac = 0.0;

  // Need non-const V0 for some reason.
  AliFemtoXi *mutable_Xi = const_cast<AliFemtoXi*>(Xi);

  // loop through the 8 points of the 'NominalTpcPoint' methods
  for (int i = 0; i < 8; i++) {
    const AliFemtoThreeVector &bac_vec = mutable_Xi->NominalTpcPointBac(i),
                              &trac_vec = track->NominalTpcPoint(i);

    if (TpcPointIsUnset(trac_vec)) {
      break;
    }

    const bool bad_bac = TpcPointIsUnset(bac_vec);

    if (bad_bac) {
      break;
    }

    if (!bad_bac) {
      avgSepBac += (bac_vec - trac_vec).Mag();
      countBac++;
    }
  }
  if (countBac != 0) {
    bac_output->Fill(avgSepBac / countBac);
  }
}

static void StoreAvgSepBetweenXiBacAndV0(const AliFemtoXi *Xi,
                                         const AliFemtoV0 *V0,
                                         TH1D *bacPos_output,
                                         TH1D *bacNeg_output)
{
  // store number of successful points for bachelor pion
  int countBacPos = 0, countBacNeg = 0;

  // sums of the separation magnitude
  double avgSepBacPos = 0.0, avgSepBacNeg = 0.0;

  // Need non-const V0 for some reason.
  AliFemtoXi *mutable_Xi = const_cast<AliFemtoXi*>(Xi);
  AliFemtoV0 *mutable_V0 = const_cast<AliFemtoV0*>(V0);

  // loop through the 8 points of the 'NominalTpcPoint' methods
  for (int i = 0; i < 8; i++) {
    const AliFemtoThreeVector &bac_vec = mutable_Xi->NominalTpcPointBac(i),
                              &pos_vec = mutable_V0->NominalTpcPointPos(i),
                              &neg_vec = mutable_V0->NominalTpcPointNeg(i);

    const bool bad_bac = TpcPointIsUnset(bac_vec),
               bad_pos = TpcPointIsUnset(pos_vec),
               bad_neg = TpcPointIsUnset(neg_vec);

    if(!bad_bac && !bad_pos)
    {
      avgSepBacPos += (bac_vec - pos_vec).Mag();
      countBacPos++;
    }

    if(!bad_bac && !bad_neg)
    {
      avgSepBacNeg += (bac_vec - neg_vec).Mag();
      countBacNeg++;
    }
  }

  if (countBacPos != 0) {
    bacPos_output->Fill(avgSepBacPos / countBacPos);
  }
  if (countBacNeg != 0) {
    bacNeg_output->Fill(avgSepBacNeg / countBacNeg);
  }

}

//____________________________
void AliFemtoAvgSepCorrFctn::AddRealPair(AliFemtoPair *pair)
{
  // only add passing pairs if cut is set
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  const AliFemtoParticle *track_1 = pair->Track1(),
                         *track_2 = pair->Track2();

  switch (fPairType) {
  // 2 tracks
  case kTracks:
    StoreAvgSepBetweenTracks(track_1->Track(),
                             track_2->Track(),
                             fNumerator);
    break;

  // track + V0
  case kTrackV0:
    StoreAvgSepBetweenV0AndTrack(track_1->V0(),
                                 track_2->Track(),
                                 fNumeratorPos,
                                 fNumeratorNeg);
    break;

  // 2 V0s
  case kV0s:
    StoreAvgSepBetweenV0s(track_1->V0(),
                          track_2->V0(),
                          fNumeratorPosPos,
                          fNumeratorPosNeg,
                          fNumeratorNegPos,
                          fNumeratorNegNeg);
    break;

  // track + Xi
  case kTrackXi:
    StoreAvgSepBetweenXiBacAndTrack(track_1->Xi(),
                                    track_2->Track(),
                                    fNumeratorBacTrack);
    StoreAvgSepBetweenV0AndTrack((AliFemtoV0*)track_1->Xi(),
                                              track_2->Track(),
                                              fNumeratorPos,
                                              fNumeratorNeg);
    break;

  // Xi + V0
  case kXiV0:
    StoreAvgSepBetweenV0s(track_1->Xi(),
                              track_2->V0(),
                              fNumeratorPosPos,
                              fNumeratorPosNeg,
                              fNumeratorNegPos,
                              fNumeratorNegNeg);
    StoreAvgSepBetweenXiBacAndV0(track_1->Xi(),
                                 track_2->V0(),
                                 fNumeratorBacPos,
                                 fNumeratorBacNeg);
    break;
  }
}
//____________________________
void AliFemtoAvgSepCorrFctn::AddMixedPair(AliFemtoPair *pair)
{
  // add mixed (background) pair

  // only add passing pairs if cut is set
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  const AliFemtoParticle *track_1 = pair->Track1(),
                         *track_2 = pair->Track2();

  switch (fPairType) {
  // 2 tracks
  case kTracks:
    StoreAvgSepBetweenTracks(track_1->Track(),
                             track_2->Track(),
                             fDenominator);
    break;

  // track + V0
  case kTrackV0:
    StoreAvgSepBetweenV0AndTrack(track_1->V0(),
                                 track_2->Track(),
                                 fDenominatorPos,
                                 fDenominatorNeg);
    break;

  // 2 V0s
  case kV0s:
    StoreAvgSepBetweenV0s(track_1->V0(),
                          track_2->V0(),
                          fDenominatorPosPos,
                          fDenominatorPosNeg,
                          fDenominatorNegPos,
                          fDenominatorNegNeg);
    break;

  // track + Xi
  case kTrackXi:
    StoreAvgSepBetweenXiBacAndTrack(track_1->Xi(),
                                    track_2->Track(),
                                    fDenominatorBacTrack);
    StoreAvgSepBetweenV0AndTrack((AliFemtoV0*)track_1->Xi(),
                                              track_2->Track(),
                                              fDenominatorPos,
                                              fDenominatorNeg);
    break;

  // Xi + V0
  case kXiV0:
    StoreAvgSepBetweenV0s(track_1->Xi(),
                              track_2->V0(),
                              fDenominatorPosPos,
                              fDenominatorPosNeg,
                              fDenominatorNegPos,
                              fDenominatorNegNeg);
    StoreAvgSepBetweenXiBacAndV0(track_1->Xi(),
                                 track_2->V0(),
                                 fDenominatorBacPos,
                                 fDenominatorBacNeg);
    break;
  }

}
//____________________________
void AliFemtoAvgSepCorrFctn::Write()
{
  // Write out neccessary objects
  if (fPairType == kTracks) {
    fNumerator->Write();
    fDenominator->Write();
    fRatio->Write();
  } else if (fPairType == kTrackV0) {
    fNumeratorPos->Write();
    fDenominatorPos->Write();
    fNumeratorNeg->Write();
    fDenominatorNeg->Write();
  } else if (fPairType == kV0s) {
    fNumeratorPosPos->Write();
    fDenominatorPosPos->Write();
    fNumeratorPosNeg->Write();
    fDenominatorPosNeg->Write();
    fNumeratorNegPos->Write();
    fDenominatorNegPos->Write();
    fNumeratorNegNeg->Write();
    fDenominatorNegNeg->Write();
  } else if (fPairType == kTrackXi) {
    fNumeratorPos->Write();
    fDenominatorPos->Write();
    fNumeratorNeg->Write();
    fDenominatorNeg->Write();
    fNumeratorBacTrack->Write();
    fDenominatorBacTrack->Write();
  }
    else if (fPairType == kXiV0) {
    fNumeratorPosPos->Write();
    fDenominatorPosPos->Write();
    fNumeratorPosNeg->Write();
    fDenominatorPosNeg->Write();
    fNumeratorNegPos->Write();
    fDenominatorNegPos->Write();
    fNumeratorNegNeg->Write();
    fDenominatorNegNeg->Write();
    fNumeratorBacPos->Write();
    fDenominatorBacPos->Write();
    fNumeratorBacNeg->Write();
    fDenominatorBacNeg->Write();
  }
}
//______________________________
TList *AliFemtoAvgSepCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  if (fPairType == kTracks) {
    tOutputList->Add(fNumerator);
    tOutputList->Add(fDenominator);
    tOutputList->Add(fRatio);
  } else if (fPairType == kTrackV0) {
    tOutputList->Add(fNumeratorPos);
    tOutputList->Add(fDenominatorPos);
    tOutputList->Add(fNumeratorNeg);
    tOutputList->Add(fDenominatorNeg);
  } else if (fPairType == kV0s) {
    tOutputList->Add(fNumeratorPosPos);
    tOutputList->Add(fDenominatorPosPos);
    tOutputList->Add(fNumeratorPosNeg);
    tOutputList->Add(fDenominatorPosNeg);
    tOutputList->Add(fNumeratorNegPos);
    tOutputList->Add(fDenominatorNegPos);
    tOutputList->Add(fNumeratorNegNeg);
    tOutputList->Add(fDenominatorNegNeg);
  } else if (fPairType == kTrackXi) {
    tOutputList->Add(fNumeratorPos);
    tOutputList->Add(fDenominatorPos);
    tOutputList->Add(fNumeratorNeg);
    tOutputList->Add(fDenominatorNeg);
    tOutputList->Add(fNumeratorBacTrack);
    tOutputList->Add(fDenominatorBacTrack);
  }
    else if (fPairType == kXiV0) {
    tOutputList->Add(fNumeratorPosPos);
    tOutputList->Add(fDenominatorPosPos);
    tOutputList->Add(fNumeratorPosNeg);
    tOutputList->Add(fDenominatorPosNeg);
    tOutputList->Add(fNumeratorNegPos);
    tOutputList->Add(fDenominatorNegPos);
    tOutputList->Add(fNumeratorNegNeg);
    tOutputList->Add(fDenominatorNegNeg);
    tOutputList->Add(fNumeratorBacPos);
    tOutputList->Add(fDenominatorBacPos);
    tOutputList->Add(fNumeratorBacNeg);
    tOutputList->Add(fDenominatorBacNeg);
  }
  return tOutputList;
}


void AliFemtoAvgSepCorrFctn::SetPairType(AliFemtoPairType pairtype)
{
  fPairType = pairtype;
}
