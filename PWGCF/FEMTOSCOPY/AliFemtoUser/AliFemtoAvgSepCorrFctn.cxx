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
AliFemtoAvgSepCorrFctn::AliFemtoAvgSepCorrFctn(char *title, const int &nbins, const float &Low, const float &High):
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
  fRatio(0),
  fPairType(kTracks)
{
  // copy constructor
  fNumerator = new TH1D(*aCorrFctn.fNumerator);
  fDenominator = new TH1D(*aCorrFctn.fDenominator);

  fNumeratorPos = new TH1D(*aCorrFctn.fNumerator);
  fDenominatorPos = new TH1D(*aCorrFctn.fDenominator);
  fNumeratorNeg = new TH1D(*aCorrFctn.fNumerator);
  fDenominatorNeg = new TH1D(*aCorrFctn.fDenominator);

  fNumeratorPosPos = new TH1D(*aCorrFctn.fNumerator);
  fDenominatorPosPos = new TH1D(*aCorrFctn.fDenominator);
  fNumeratorPosNeg = new TH1D(*aCorrFctn.fNumerator);
  fDenominatorPosNeg = new TH1D(*aCorrFctn.fDenominator);
  fNumeratorNegPos = new TH1D(*aCorrFctn.fNumerator);
  fDenominatorNegPos = new TH1D(*aCorrFctn.fDenominator);
  fNumeratorNegNeg = new TH1D(*aCorrFctn.fNumerator);
  fDenominatorNegNeg = new TH1D(*aCorrFctn.fDenominator);

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
//____________________________
void AliFemtoAvgSepCorrFctn::AddRealPair(AliFemtoPair *pair)
{
  // add true pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double avgSep = 0;
  AliFemtoThreeVector first, second, tmp;


  if (fPairType == 0) { //2 tracks
    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->Track()->NominalTpcPoint(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->Track()->NominalTpcPoint(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fNumerator->Fill(avgSep);
  } else if (fPairType == 1) { // track + V0
    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointPos(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());
      //cout<<"V0pos x: "<<tmp.x()<<", y: "<<tmp.y()<<", z: "<<tmp.z()<<endl;
      //cout<<"V0pos x: "<<first.x()<<", y: "<<first.y()<<", z: "<<first.z()<<endl;

      tmp = pair->Track2()->Track()->NominalTpcPoint(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());
      //cout<<"Track x: "<<tmp.x()<<", y: "<<tmp.y()<<", z: "<<tmp.z()<<endl;
      //cout<<"Track x: "<<second.x()<<", y: "<<second.y()<<", z: "<<second.z()<<endl;

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    //cout<<"****************************************************"<<endl;
    avgSep /= 8;
    fNumeratorPos->Fill(avgSep);
    //cout<<"Track + Pos V0 Avg Sep: "<<avgSep<<endl;
    avgSep = 0;

    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointNeg(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());
      //cout<<"V0 X: "<<tmp.x()<<endl;

      tmp = pair->Track2()->Track()->NominalTpcPoint(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());
      //cout<<"Track X: "<<tmp.x()<<endl;

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    fNumeratorNeg->Fill(avgSep);

  } else if (fPairType == 2) {
    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointPos(i);
      //cout<<"X pos: "<<tmp.x()<<endl;
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->V0()->NominalTpcPointPos(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fNumeratorPosPos->Fill(avgSep);
    //cout<<"PovV0 + PosV0 Avg Sep: "<<avgSep<<endl;

    avgSep = 0;

    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointPos(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->V0()->NominalTpcPointNeg(i);
      //cout<<"X neg: "<<tmp.x()<<endl;
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fNumeratorPosNeg->Fill(avgSep);

    avgSep = 0;

    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointNeg(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->V0()->NominalTpcPointPos(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fNumeratorNegPos->Fill(avgSep);


    avgSep = 0;

    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointNeg(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->V0()->NominalTpcPointNeg(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fNumeratorNegNeg->Fill(avgSep);


  }


  //2 daughters of first V0





}
//____________________________
void AliFemtoAvgSepCorrFctn::AddMixedPair(AliFemtoPair *pair)
{
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double avgSep = 0;
  AliFemtoThreeVector first, second, tmp;


  if (fPairType == 0) { //2 tracks
    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->Track()->NominalTpcPoint(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->Track()->NominalTpcPoint(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fDenominator->Fill(avgSep);
  } else if (fPairType == 1) { // track + V0
    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointPos(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->Track()->NominalTpcPoint(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fDenominatorPos->Fill(avgSep);

    avgSep = 0;

    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointNeg(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->Track()->NominalTpcPoint(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    fDenominatorNeg->Fill(avgSep);

  } else if (fPairType == 2) {
    for (int i = 0; i < 8 ; i++) {

      tmp = pair->Track1()->V0()->NominalTpcPointPos(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->V0()->NominalTpcPointPos(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fDenominatorPosPos->Fill(avgSep);

    avgSep = 0;

    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointPos(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->V0()->NominalTpcPointNeg(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fDenominatorPosNeg->Fill(avgSep);

    avgSep = 0;

    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointNeg(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->V0()->NominalTpcPointPos(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fDenominatorNegPos->Fill(avgSep);

    avgSep = 0;

    for (int i = 0; i < 8 ; i++) {
      tmp = pair->Track1()->V0()->NominalTpcPointNeg(i);
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());

      tmp = pair->Track2()->V0()->NominalTpcPointNeg(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z());

      avgSep += TMath::Sqrt(((double)first.x() - (double)second.x()) * ((double)first.x() - (double)second.x()) + ((double)first.y() - (double)second.y()) * ((double)first.y() - second.y()) + ((double)first.z() - (double)second.z()) * ((double)first.z() - (double)second.z()));
    }
    avgSep /= 8;
    fDenominatorNegNeg->Fill(avgSep);


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
  }
  return tOutputList;
}


void AliFemtoAvgSepCorrFctn::SetPairType(AliFemtoPairType pairtype)
{
  fPairType = pairtype;
}
