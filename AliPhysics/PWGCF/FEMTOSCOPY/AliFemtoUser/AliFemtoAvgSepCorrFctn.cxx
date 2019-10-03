///
/// \file AliFemtoAvgSepCorrFctn.cxx
///

#include "AliFemtoAvgSepCorrFctn.h"


//____________________________
AliFemtoAvgSepCorrFctn::AliFemtoAvgSepCorrFctn(const char *suffix,
                                               const int &nbins,
                                               const float &Low,
                                               const float &High):
  AliFemtoCorrFctn(),
  fPairType(kTracks),
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
  fRatio(0)
{
  TString suf = suffix;

  auto hist_title = [] (TString title) {
    return "Average Separation " + title + "; $\\left< \\Delta X \\right>$ (cm)";
  };

  fNumerator = new TH1D("Num" + suf, hist_title("Numerator"), nbins, Low, High);
  fDenominator = new TH1D("Den" + suf, hist_title("Denominator"), nbins, Low, High);
  fRatio = new TH1D("Rat" + suf, hist_title("Ratio"), nbins, Low, High);

  fNumeratorPos = new TH1D("NumV0TrackPos" + suf, "Num : Track & V0.Pos", nbins, Low, High);
  fDenominatorPos = new TH1D("DenV0TrackPos" + suf, "Den : Track & V0.Pos", nbins, Low, High);
  fNumeratorNeg = new TH1D("NumV0TrackNeg" + suf, "Num : Track & V0.Neg", nbins, Low, High);
  fDenominatorNeg = new TH1D("DenV0TrackNeg" + suf, "Den : Track & V0.Neg", nbins, Low, High);

  fNumeratorPosPos   = new TH1D("NumV0sPosPos" + suf, "Num : V0-1.Pos & V0-2.Pos", nbins, Low, High);
  fDenominatorPosPos = new TH1D("DenV0sPosPos" + suf, "Den : V0-1.Pos & V0-2.Pos", nbins, Low, High);
  fNumeratorPosNeg   = new TH1D("NumV0sPosNeg" + suf, "Num : V0-1.Pos & V0-2.Neg", nbins, Low, High);
  fDenominatorPosNeg = new TH1D("DenV0sPosNeg" + suf, "Den : V0-1.Pos & V0-2.Neg", nbins, Low, High);
  fNumeratorNegPos   = new TH1D("NumV0sNegPos" + suf, "Num : V0-1.Neg & V0-2.Pos", nbins, Low, High);
  fDenominatorNegPos = new TH1D("DenV0sNegPos" + suf, "Den : V0-1.Neg & V0-2.Pos", nbins, Low, High);
  fNumeratorNegNeg   = new TH1D("NumV0sNegNeg" + suf, "Num : V0-1.Neg & V0-2.Neg", nbins, Low, High);
  fDenominatorNegNeg = new TH1D("DenV0sNegNeg" + suf, "Den : V0-1.Neg & V0-2.Neg", nbins, Low, High);

  TString title = suffix;

  fNumeratorBacTrack = new TH1D("NumXiTrackBac" + suf, title, nbins, Low, High);
  fDenominatorBacTrack = new TH1D("DenXiTrackBac" + suf, title, nbins, Low, High);

  fNumeratorBacPos = new TH1D("NumXiBacV0Pos" + suf, title, nbins, Low, High);
  fDenominatorBacPos = new TH1D("DenXiBacV0Pos" + suf, title, nbins, Low, High);

  fNumeratorBacNeg = new TH1D("NumXiBacV0Neg" + suf, title, nbins, Low, High);
  fDenominatorBacNeg = new TH1D("DenXiBacV0Neg" + suf, title, nbins, Low, High);

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
  AliFemtoCorrFctn(aCorrFctn),
  fPairType(aCorrFctn.fPairType),
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
  fRatio(0)
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

  AliFemtoCorrFctn::operator=(aCorrFctn);

  *fNumerator = *aCorrFctn.fNumerator;
  *fDenominator = *aCorrFctn.fDenominator;

  *fNumeratorPos = *aCorrFctn.fNumeratorPos;
  *fDenominatorPos = *aCorrFctn.fDenominatorPos;
  *fNumeratorNeg = *aCorrFctn.fNumeratorNeg;
  *fDenominatorNeg = *aCorrFctn.fDenominatorNeg;

  *fNumeratorPosPos = *aCorrFctn.fNumeratorPosPos;
  *fDenominatorPosPos = *aCorrFctn.fDenominatorPosPos;
  *fNumeratorPosNeg = *aCorrFctn.fNumeratorPosNeg;
  *fDenominatorPosNeg = *aCorrFctn.fDenominatorPosNeg;
  *fNumeratorNegPos = *aCorrFctn.fNumeratorNegPos;
  *fDenominatorNegPos = *aCorrFctn.fDenominatorNegPos;
  *fNumeratorNegNeg = *aCorrFctn.fNumeratorNegNeg;
  *fDenominatorNegNeg = *aCorrFctn.fDenominatorNegNeg;

  *fNumeratorBacTrack = *aCorrFctn.fNumeratorBacTrack;
  *fDenominatorBacTrack = *aCorrFctn.fDenominatorBacTrack;

  *fNumeratorBacPos = *aCorrFctn.fNumeratorBacPos;
  *fDenominatorBacPos = *aCorrFctn.fDenominatorBacPos;
  *fNumeratorBacNeg = *aCorrFctn.fNumeratorBacNeg;
  *fDenominatorBacNeg = *aCorrFctn.fDenominatorBacNeg;

  *fRatio = *aCorrFctn.fRatio;

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
  AliFemtoString report = "Qinv Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());
  report += Form("Number of entries in ratio:\t%E\n", fRatio->GetEntries());
  return report;
}

// Function to check if the vector has a point near the default value (-9999)
// (i.e. it was never set)
static bool TpcPointIsUnset(const AliFemtoThreeVector& v) {
  return v.x() < -9000. ||
         v.y() < -9000. ||
         v.z() < -9000.;
}

static void StoreAvgSepBetweenV0AndTrack(const AliFemtoV0 &V0,
                                         const AliFemtoTrack &track,
                                         TH1D *pos_output,
                                         TH1D *neg_output)
{
  double neg = -1.0, pos = -1.0;
  AliFemtoPair::CalcAvgSepTrackV0(track, V0, neg, pos);

  neg_output->Fill(neg);
  pos_output->Fill(pos);
}

static void StoreAvgSepBetweenV0s(const AliFemtoV0 &V0_1,
                                  const AliFemtoV0 &V0_2,
                                  TH1D *pospos_output,
                                  TH1D *posneg_output,
                                  TH1D *negpos_output,
                                  TH1D *negneg_output)
{
  // sums of the separation magnitude
  double avgSepPosPos = 0.0,
         avgSepPosNeg = 0.0,
         avgSepNegPos = 0.0,
         avgSepNegNeg = 0.0;

  AliFemtoPair::CalcAvgSepV0V0(V0_1, V0_2, avgSepNegNeg, avgSepNegPos, avgSepPosNeg, avgSepPosPos);

  pospos_output->Fill(avgSepPosPos);
  posneg_output->Fill(avgSepPosNeg);
  negpos_output->Fill(avgSepNegPos);
  negneg_output->Fill(avgSepNegNeg);
}

static void StoreAvgSepBetweenXiBacAndTrack(const AliFemtoXi &Xi,
                                            const AliFemtoTrack &track,
                                            TH1D *bac_output)
{
  // store number of successful points for bachelor pion
  int countBac = 0;

  // sums of the separation magnitude
  double avgSepBac = 0.0;

  // loop through the 'NominalTpcPoints'
  for (int i = 0; i < 9; i++) {
    const AliFemtoThreeVector &bac_vec = Xi.NominalTpcPointBac(i),
                              &trac_vec = track.NominalTpcPoint(i);

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

static void StoreAvgSepBetweenXiBacAndV0(const AliFemtoXi &Xi,
                                         const AliFemtoV0 &V0,
                                         TH1D *bacPos_output,
                                         TH1D *bacNeg_output)
{
  // store number of successful points for bachelor pion
  int countBacPos = 0, countBacNeg = 0;

  // sums of the separation magnitude
  double avgSepBacPos = 0.0, avgSepBacNeg = 0.0;

  for (int i = 0; i < 9; i++) {
    const AliFemtoThreeVector &bac_vec = Xi.NominalTpcPointBac(i),
                              &pos_vec = V0.NominalTpcPointPos(i),
                              &neg_vec = V0.NominalTpcPointNeg(i);

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
    fNumerator->Fill(pair->NominalTpcAverageSeparationTracks());
    break;

  // track + V0
  case kTrackV0:
    StoreAvgSepBetweenV0AndTrack(*track_1->V0(),
                                 *track_2->Track(),
                                 fNumeratorPos,
                                 fNumeratorNeg);
    break;

  // 2 V0s
  case kV0s:
    StoreAvgSepBetweenV0s(*track_1->V0(),
                          *track_2->V0(),
                          fNumeratorPosPos,
                          fNumeratorPosNeg,
                          fNumeratorNegPos,
                          fNumeratorNegNeg);
    break;

  // track + Xi
  case kTrackXi:
    StoreAvgSepBetweenXiBacAndTrack(*track_1->Xi(),
                                    *track_2->Track(),
                                    fNumeratorBacTrack);
    StoreAvgSepBetweenV0AndTrack(*track_1->Xi(),
                                 *track_2->Track(),
                                 fNumeratorPos,
                                 fNumeratorNeg);
    break;

  // Xi + V0
  case kXiV0:
    StoreAvgSepBetweenV0s(*track_1->Xi(),
                          *track_2->V0(),
                          fNumeratorPosPos,
                          fNumeratorPosNeg,
                          fNumeratorNegPos,
                          fNumeratorNegNeg);
    StoreAvgSepBetweenXiBacAndV0(*track_1->Xi(),
                                 *track_2->V0(),
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
    fDenominator->Fill(pair->NominalTpcAverageSeparationTracks());
    break;

  // track + V0
  case kTrackV0:
    StoreAvgSepBetweenV0AndTrack(*track_1->V0(),
                                 *track_2->Track(),
                                 fDenominatorPos,
                                 fDenominatorNeg);
    break;

  // 2 V0s
  case kV0s:
    StoreAvgSepBetweenV0s(*track_1->V0(),
                          *track_2->V0(),
                          fDenominatorPosPos,
                          fDenominatorPosNeg,
                          fDenominatorNegPos,
                          fDenominatorNegNeg);
    break;

  // track + Xi
  case kTrackXi:
    StoreAvgSepBetweenXiBacAndTrack(*track_1->Xi(),
                                    *track_2->Track(),
                                    fDenominatorBacTrack);
    StoreAvgSepBetweenV0AndTrack(*track_1->Xi(),
                                 *track_2->Track(),
                                 fDenominatorPos,
                                 fDenominatorNeg);
    break;

  // Xi + V0
  case kXiV0:
    StoreAvgSepBetweenV0s(*track_1->Xi(),
                          *track_2->V0(),
                          fDenominatorPosPos,
                          fDenominatorPosNeg,
                          fDenominatorNegPos,
                          fDenominatorNegNeg);
    StoreAvgSepBetweenXiBacAndV0(*track_1->Xi(),
                                 *track_2->V0(),
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
