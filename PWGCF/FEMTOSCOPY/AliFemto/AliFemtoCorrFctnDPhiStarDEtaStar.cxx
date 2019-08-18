///
/// \file AliFemtoCorrFctnDPhiStarDEtaStar.cxx
///

#include "AliFemtoCorrFctnDPhiStarDEtaStar.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctnDPhiStarDEtaStar);
  /// \endcond
#endif

//____________________________
AliFemtoCorrFctnDPhiStarDEtaStar::AliFemtoCorrFctnDPhiStarDEtaStar(const char* title,
								   double radius=1.2,
								   const int& aEtaBins=50,
								   double aEtaRangeLow=0.0,
								   double aEtaRangeUp=0.2,
								   const int& aPhiStarBins=50,
								   double aPhiStarRangeLow=0.0,
								   double aPhiStarRangeUp=0.2):
AliFemtoCorrFctn(),
  fDPhiStarDEtaStarNumerator(0),
  fDPhiStarDEtaStarDenominator(0),
  fDPhiStarDEtaStarNumeratorTrackPos(0),
  fDPhiStarDEtaStarDenominatorTrackPos(0),
  fDPhiStarDEtaStarNumeratorTrackNeg(0),
  fDPhiStarDEtaStarDenominatorTrackNeg(0),
  fDPhiStarDEtaStarNumeratorPosPos(0),
  fDPhiStarDEtaStarDenominatorPosPos(0),
  fDPhiStarDEtaStarNumeratorPosNeg(0),
  fDPhiStarDEtaStarDenominatorPosNeg(0),
  fDPhiStarDEtaStarNumeratorNegPos(0),
  fDPhiStarDEtaStarDenominatorNegPos(0),
  fDPhiStarDEtaStarNumeratorNegNeg(0),
  fDPhiStarDEtaStarDenominatorNegNeg(0),
  fEtaStarRangeLow(0),
  fEtaStarRangeUp(0),
  fPhiStarRangeLow(0),
  fPhiStarRangeUp(0),
  fMinRad(1.2),
  fPairType(kTracks)
{

  // Set up lower and upper range of EtaStar and PhiStar:
  fEtaStarRangeLow = aEtaRangeLow;
  fEtaStarRangeUp = aEtaRangeUp;
  fPhiStarRangeLow = aPhiStarRangeLow;
  fPhiStarRangeUp = aPhiStarRangeUp;

  // Set up radial distance
  fMinRad = radius;

  // Set up numerator:
  char tTitNum[101] = "NumDPhiStarDEtaStar";
  strncat(tTitNum, title, 100);
  fDPhiStarDEtaStarNumerator = new TH2D(tTitNum, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);
  // Set up denominator:
  char tTitDen[101] = "DenDPhiStarDEtaStar";
  strncat(tTitDen, title, 100);
  fDPhiStarDEtaStarDenominator = new TH2D(tTitDen, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);

  // Set up numerator:
  char tTitNumTrackPos[101] = "NumDPhiStarDEtaStarTrackPos";
  strncat(tTitNumTrackPos, title, 100);
  fDPhiStarDEtaStarNumeratorTrackPos = new TH2D(tTitNumTrackPos, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);
  // Set up denominator:
  char tTitDenTrackPos[101] = "DenDPhiStarDEtaStarTrackPos";
  strncat(tTitDenTrackPos, title, 100);
  fDPhiStarDEtaStarDenominatorTrackPos = new TH2D(tTitDenTrackPos, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);

  // Set up numerator:
  char tTitNumTrackNeg[101] = "NumDPhiStarDEtaStarTrackNeg";
  strncat(tTitNumTrackNeg, title, 100);
  fDPhiStarDEtaStarNumeratorTrackNeg = new TH2D(tTitNumTrackNeg, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);
  // Set up denominator:
  char tTitDenTrackNeg[101] = "DenDPhiStarDEtaStarTrackNeg";
  strncat(tTitDenTrackNeg, title, 100);
  fDPhiStarDEtaStarDenominatorTrackNeg = new TH2D(tTitDenTrackNeg, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);

  // Set up numerator:
  char tTitNumPosPos[101] = "NumDPhiStarDEtaStarPosPos";
  strncat(tTitNumPosPos, title, 100);
  fDPhiStarDEtaStarNumeratorPosPos = new TH2D(tTitNumPosPos, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);
  // Set up denominator:
  char tTitDenPosPos[101] = "DenDPhiStarDEtaStarPosPos";
  strncat(tTitDenPosPos, title, 100);
  fDPhiStarDEtaStarDenominatorPosPos = new TH2D(tTitDenPosPos, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);

  // Set up numerator:
  char tTitNumPosNeg[101] = "NumDPhiStarDEtaStarPosNeg";
  strncat(tTitNumPosNeg, title, 100);
  fDPhiStarDEtaStarNumeratorPosNeg = new TH2D(tTitNumPosNeg, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);
  // Set up denominator:
  char tTitDenPosNeg[101] = "DenDPhiStarDEtaStarPosNeg";
  strncat(tTitDenPosNeg, title, 100);
  fDPhiStarDEtaStarDenominatorPosNeg = new TH2D(tTitDenPosNeg, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);

  // Set up numerator:
  char tTitNumNegPos[101] = "NumDPhiStarDEtaStarNegPos";
  strncat(tTitNumNegPos, title, 100);
  fDPhiStarDEtaStarNumeratorNegPos = new TH2D(tTitNumNegPos, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);
  // Set up denominator:
  char tTitDenNegPos[101] = "DenDPhiStarDEtaStarNegPos";
  strncat(tTitDenNegPos, title, 100);
  fDPhiStarDEtaStarDenominatorNegPos = new TH2D(tTitDenNegPos, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);

  // Set up numerator:
  char tTitNumNegNeg[101] = "NumDPhiStarDEtaStarNegNeg";
  strncat(tTitNumNegNeg, title, 100);
  fDPhiStarDEtaStarNumeratorNegNeg = new TH2D(tTitNumNegNeg, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);
  // Set up denominator:
  char tTitDenNegNeg[101] = "DenDPhiStarDEtaStarNegNeg";
  strncat(tTitDenNegNeg, title, 100);
  fDPhiStarDEtaStarDenominatorNegNeg = new TH2D(tTitDenNegNeg, title, aEtaBins, fEtaStarRangeLow, fEtaStarRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);

  // To enable error bars calculations:
  fDPhiStarDEtaStarNumerator->Sumw2();
  fDPhiStarDEtaStarDenominator->Sumw2();
  fDPhiStarDEtaStarNumeratorTrackPos->Sumw2();
  fDPhiStarDEtaStarDenominatorTrackPos->Sumw2();
  fDPhiStarDEtaStarNumeratorTrackNeg->Sumw2();
  fDPhiStarDEtaStarDenominatorTrackNeg->Sumw2();
  fDPhiStarDEtaStarNumeratorPosPos->Sumw2();
  fDPhiStarDEtaStarDenominatorPosPos->Sumw2();
  fDPhiStarDEtaStarNumeratorPosNeg->Sumw2();
  fDPhiStarDEtaStarDenominatorPosNeg->Sumw2();
  fDPhiStarDEtaStarNumeratorNegPos->Sumw2();
  fDPhiStarDEtaStarDenominatorNegPos->Sumw2();
  fDPhiStarDEtaStarNumeratorNegNeg->Sumw2();
  fDPhiStarDEtaStarDenominatorNegNeg->Sumw2();
}

//____________________________
AliFemtoCorrFctnDPhiStarDEtaStar::AliFemtoCorrFctnDPhiStarDEtaStar(const AliFemtoCorrFctnDPhiStarDEtaStar& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiStarDEtaStarNumerator(0),
  fDPhiStarDEtaStarDenominator(0),
  fDPhiStarDEtaStarNumeratorTrackPos(0),
  fDPhiStarDEtaStarDenominatorTrackPos(0),
  fDPhiStarDEtaStarNumeratorTrackNeg(0),
  fDPhiStarDEtaStarDenominatorTrackNeg(0),
  fDPhiStarDEtaStarNumeratorPosPos(0),
  fDPhiStarDEtaStarDenominatorPosPos(0),
  fDPhiStarDEtaStarNumeratorPosNeg(0),
  fDPhiStarDEtaStarDenominatorPosNeg(0),
  fDPhiStarDEtaStarNumeratorNegPos(0),
  fDPhiStarDEtaStarDenominatorNegPos(0),
  fDPhiStarDEtaStarNumeratorNegNeg(0),
  fDPhiStarDEtaStarDenominatorNegNeg(0),
  fEtaStarRangeLow(0),
  fEtaStarRangeUp(0),
  fPhiStarRangeLow(0),
  fPhiStarRangeUp(0),
  fMinRad(1.2),
  fPairType(kTracks)
{
  // Copy constructor
  if (aCorrFctn.fDPhiStarDEtaStarNumerator)
    fDPhiStarDEtaStarNumerator = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumerator);
  else
    fDPhiStarDEtaStarNumerator = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominator)
    fDPhiStarDEtaStarDenominator = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominator);
  else
    fDPhiStarDEtaStarDenominator = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorTrackPos)
    fDPhiStarDEtaStarNumeratorTrackPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorTrackPos);
  else
    fDPhiStarDEtaStarNumeratorTrackPos = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorTrackPos)
    fDPhiStarDEtaStarDenominatorTrackPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorTrackPos);
  else
    fDPhiStarDEtaStarDenominatorTrackPos = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorTrackNeg)
    fDPhiStarDEtaStarNumeratorTrackNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorTrackNeg);
  else
    fDPhiStarDEtaStarNumeratorTrackNeg = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorTrackNeg)
    fDPhiStarDEtaStarDenominatorTrackNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorTrackNeg);
  else
    fDPhiStarDEtaStarDenominatorTrackNeg = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorPosPos)
    fDPhiStarDEtaStarNumeratorPosPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorPosPos);
  else
    fDPhiStarDEtaStarNumeratorPosPos = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorPosPos)
    fDPhiStarDEtaStarDenominatorPosPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorPosPos);
  else
    fDPhiStarDEtaStarDenominatorPosPos = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorPosNeg)
    fDPhiStarDEtaStarNumeratorPosNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorPosNeg);
  else
    fDPhiStarDEtaStarNumeratorPosNeg = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorPosNeg)
    fDPhiStarDEtaStarDenominatorPosNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorPosNeg);
  else
    fDPhiStarDEtaStarDenominatorPosNeg = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorNegPos)
    fDPhiStarDEtaStarNumeratorNegPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorNegPos);
  else
    fDPhiStarDEtaStarNumeratorNegPos = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorNegPos)
    fDPhiStarDEtaStarDenominatorNegPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorNegPos);
  else
    fDPhiStarDEtaStarDenominatorNegPos = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorNegNeg)
    fDPhiStarDEtaStarNumeratorNegNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorNegNeg);
  else
    fDPhiStarDEtaStarNumeratorNegNeg = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorNegNeg)
    fDPhiStarDEtaStarDenominatorNegNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorNegNeg);
  else
    fDPhiStarDEtaStarDenominatorNegNeg = 0;

  fEtaStarRangeLow = aCorrFctn.fEtaStarRangeLow;
  fEtaStarRangeUp = aCorrFctn.fEtaStarRangeUp;
  fPhiStarRangeLow = aCorrFctn.fPhiStarRangeLow;
  fPhiStarRangeUp = aCorrFctn.fPhiStarRangeUp;
  fMinRad = aCorrFctn.fMinRad;
  fPairType = aCorrFctn.fPairType;
}

//____________________________
AliFemtoCorrFctnDPhiStarDEtaStar::~AliFemtoCorrFctnDPhiStarDEtaStar(){
  // Destructor
  delete fDPhiStarDEtaStarNumerator;
  delete fDPhiStarDEtaStarDenominator;
  delete fDPhiStarDEtaStarNumeratorTrackPos;
  delete fDPhiStarDEtaStarDenominatorTrackPos;
  delete fDPhiStarDEtaStarNumeratorTrackNeg;
  delete fDPhiStarDEtaStarDenominatorTrackNeg;
  delete fDPhiStarDEtaStarNumeratorPosPos;
  delete fDPhiStarDEtaStarDenominatorPosPos;
  delete fDPhiStarDEtaStarNumeratorPosNeg;
  delete fDPhiStarDEtaStarDenominatorPosNeg;
  delete fDPhiStarDEtaStarNumeratorNegPos;
  delete fDPhiStarDEtaStarDenominatorNegPos;
  delete fDPhiStarDEtaStarNumeratorNegNeg;
  delete fDPhiStarDEtaStarDenominatorNegNeg;
}

//_________________________
AliFemtoCorrFctnDPhiStarDEtaStar& AliFemtoCorrFctnDPhiStarDEtaStar::operator=(const AliFemtoCorrFctnDPhiStarDEtaStar& aCorrFctn)
{
  // Assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fDPhiStarDEtaStarNumerator)
    fDPhiStarDEtaStarNumerator = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumerator);
  else
    fDPhiStarDEtaStarNumerator = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominator)
    fDPhiStarDEtaStarDenominator = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominator);
  else
    fDPhiStarDEtaStarDenominator = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorTrackPos)
    fDPhiStarDEtaStarNumeratorTrackPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorTrackPos);
  else
    fDPhiStarDEtaStarNumeratorTrackPos = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorTrackPos)
    fDPhiStarDEtaStarDenominatorTrackPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorTrackPos);
  else
    fDPhiStarDEtaStarDenominatorTrackPos = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorTrackNeg)
    fDPhiStarDEtaStarNumeratorTrackNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorTrackNeg);
  else
    fDPhiStarDEtaStarNumeratorTrackNeg = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorTrackNeg)
    fDPhiStarDEtaStarDenominatorTrackNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorTrackNeg);
  else
    fDPhiStarDEtaStarDenominatorTrackNeg = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorPosPos)
    fDPhiStarDEtaStarNumeratorPosPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorPosPos);
  else
    fDPhiStarDEtaStarNumeratorPosPos = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorPosPos)
    fDPhiStarDEtaStarDenominatorPosPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorPosPos);
  else
    fDPhiStarDEtaStarDenominatorPosPos = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorPosNeg)
    fDPhiStarDEtaStarNumeratorPosNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorPosNeg);
  else
    fDPhiStarDEtaStarNumeratorPosNeg = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorPosNeg)
    fDPhiStarDEtaStarDenominatorPosNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorPosNeg);
  else
    fDPhiStarDEtaStarDenominatorPosNeg = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorNegPos)
    fDPhiStarDEtaStarNumeratorNegPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorNegPos);
  else
    fDPhiStarDEtaStarNumeratorNegPos = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorNegPos)
    fDPhiStarDEtaStarDenominatorNegPos = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorNegPos);
  else
    fDPhiStarDEtaStarDenominatorNegPos = 0;

  if (aCorrFctn.fDPhiStarDEtaStarNumeratorNegNeg)
    fDPhiStarDEtaStarNumeratorNegNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarNumeratorNegNeg);
  else
    fDPhiStarDEtaStarNumeratorNegNeg = 0;
  if (aCorrFctn.fDPhiStarDEtaStarDenominatorNegNeg)
    fDPhiStarDEtaStarDenominatorNegNeg = new TH2D(*aCorrFctn.fDPhiStarDEtaStarDenominatorNegNeg);
  else
    fDPhiStarDEtaStarDenominatorNegNeg = 0;

  fEtaStarRangeLow = aCorrFctn.fEtaStarRangeLow;
  fEtaStarRangeUp = aCorrFctn.fEtaStarRangeUp;
  fPhiStarRangeLow = aCorrFctn.fPhiStarRangeLow;
  fPhiStarRangeUp = aCorrFctn.fPhiStarRangeUp;
  fMinRad = aCorrFctn.fMinRad;
  fPairType = aCorrFctn.fPairType;

  return *this;
}

//_________________________
void AliFemtoCorrFctnDPhiStarDEtaStar::Finish(){
  // Here is where we should normalize, fit, etc...
  // We should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //  mShareDenominator->Draw();
  //  mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDPhiStarDEtaStar::Report()
{
  // Create report
  AliFemtoString report = "DPhi* DEta* Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n",fDPhiStarDEtaStarNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n",fDPhiStarDEtaStarDenominator->GetEntries());
  //  report += mCoulombWeight->Report();

  return report;
}

//____________________________
static void StoreDPhiStarDEtaStarBetweenTracks(const AliFemtoTrack *track1,
					       const AliFemtoTrack *track2,
					       TH2D *output,
					       Double_t minRad)
{
  double thetas1 = TMath::Pi()/2. - TMath::ATan(track1->NominalTpcPointShifted().z()/(minRad*1e2));
  double thetas2 = TMath::Pi()/2. - TMath::ATan(track2->NominalTpcPointShifted().z()/(minRad*1e2));
  double etas1 = -TMath::Log( TMath::Tan(thetas1/2.) );
  double etas2 = -TMath::Log( TMath::Tan(thetas2/2.) );
  double detas = TMath::Abs(etas1 - etas2);
  double distSft= TMath::Sqrt(TMath::Power(track1->NominalTpcPointShifted().x() -
					   track2->NominalTpcPointShifted().x(),2) +
			      TMath::Power(track1->NominalTpcPointShifted().y() -
					   track2->NominalTpcPointShifted().y(),2));
  double dPhiS = 2.0 * TMath::ATan(distSft/2./((minRad*1e2)));

  output->Fill(detas, dPhiS);
}

//____________________________
static void StoreDPhiStarDEtaStarBetweenV0AndTrack(const AliFemtoV0 *V0,
						   const AliFemtoTrack *track,
						   TH2D *output_pos,
						   TH2D *output_neg,
						   Double_t minRad)
{
  AliFemtoV0 *mut_V0 = const_cast<AliFemtoV0*>(V0);
  double thetas1_pos = TMath::Pi()/2. - TMath::ATan(mut_V0->NominalTpcPointPosShifted().z()/(minRad*1e2));
  double thetas2_pos = TMath::Pi()/2. - TMath::ATan(track->NominalTpcPointShifted().z()/(minRad*1e2));
  double etas1_pos = -TMath::Log( TMath::Tan(thetas1_pos/2.) );
  double etas2_pos = -TMath::Log( TMath::Tan(thetas2_pos/2.) );
  double detas_pos = TMath::Abs(etas1_pos - etas2_pos);
  double distSft_pos = TMath::Sqrt(TMath::Power(mut_V0->NominalTpcPointPosShifted().x() -
						track->NominalTpcPointShifted().x(),2) +
				   TMath::Power(mut_V0->NominalTpcPointPosShifted().y() -
						track->NominalTpcPointShifted().y(),2));
  double dPhiS_pos = 2.0 * TMath::ATan(distSft_pos/2./((minRad*1e2)));

  double thetas1_neg = TMath::Pi()/2. - TMath::ATan(mut_V0->NominalTpcPointNegShifted().z()/(minRad*1e2));
  double thetas2_neg = TMath::Pi()/2. - TMath::ATan(track->NominalTpcPointShifted().z()/(minRad*1e2));
  double etas1_neg = -TMath::Log( TMath::Tan(thetas1_neg/2.) );
  double etas2_neg = -TMath::Log( TMath::Tan(thetas2_neg/2.) );
  double detas_neg = TMath::Abs(etas1_neg - etas2_neg);
  double distSft_neg = TMath::Sqrt(TMath::Power(mut_V0->NominalTpcPointNegShifted().x() -
						track->NominalTpcPointShifted().x(),2) +
				   TMath::Power(mut_V0->NominalTpcPointNegShifted().y() -
						track->NominalTpcPointShifted().y(),2));
  double dPhiS_neg = 2.0 * TMath::ATan(distSft_neg/2./((minRad*1e2)));

  output_pos->Fill(detas_pos, dPhiS_pos);
  output_neg->Fill(detas_neg, dPhiS_neg);
}

//____________________________
static void StoreDPhiStarDEtaStarBetweenV0s(const AliFemtoV0 *V0_1,
					    const AliFemtoV0 *V0_2,
					    TH2D *output_pospos,
					    TH2D *output_posneg,
					    TH2D *output_negpos,
					    TH2D *output_negneg,
					    Double_t minRad)
{
  AliFemtoV0 *mut_V0_1 = const_cast<AliFemtoV0*>(V0_1);
  AliFemtoV0 *mut_V0_2 = const_cast<AliFemtoV0*>(V0_2);

  double thetas1_pospos = TMath::Pi()/2. - TMath::ATan(mut_V0_1->NominalTpcPointPosShifted().z()/(minRad*1e2));
  double thetas2_pospos = TMath::Pi()/2. - TMath::ATan(mut_V0_2->NominalTpcPointPosShifted().z()/(minRad*1e2));
  double etas1_pospos = -TMath::Log( TMath::Tan(thetas1_pospos/2.) );
  double etas2_pospos = -TMath::Log( TMath::Tan(thetas2_pospos/2.) );
  double detas_pospos = TMath::Abs(etas1_pospos - etas2_pospos);
  double distSft_pospos = TMath::Sqrt(TMath::Power(mut_V0_1->NominalTpcPointPosShifted().x() -
						   mut_V0_2->NominalTpcPointPosShifted().x(),2) +
				      TMath::Power(mut_V0_1->NominalTpcPointPosShifted().y() -
						   mut_V0_2->NominalTpcPointPosShifted().y(),2));
  double dPhiS_pospos = 2.0 * TMath::ATan(distSft_pospos/2./((minRad*1e2)));

  double thetas1_posneg = TMath::Pi()/2. - TMath::ATan(mut_V0_1->NominalTpcPointPosShifted().z()/(minRad*1e2));
  double thetas2_posneg = TMath::Pi()/2. - TMath::ATan(mut_V0_2->NominalTpcPointNegShifted().z()/(minRad*1e2));
  double etas1_posneg = -TMath::Log( TMath::Tan(thetas1_posneg/2.) );
  double etas2_posneg = -TMath::Log( TMath::Tan(thetas2_posneg/2.) );
  double detas_posneg = TMath::Abs(etas1_posneg - etas2_posneg);
  double distSft_posneg = TMath::Sqrt(TMath::Power(mut_V0_1->NominalTpcPointPosShifted().x() -
						   mut_V0_2->NominalTpcPointNegShifted().x(),2) +
				      TMath::Power(mut_V0_1->NominalTpcPointPosShifted().y() -
						   mut_V0_2->NominalTpcPointNegShifted().y(),2));
  double dPhiS_posneg = 2.0 * TMath::ATan(distSft_posneg/2./((minRad*1e2)));

  double thetas1_negpos = TMath::Pi()/2. - TMath::ATan(mut_V0_1->NominalTpcPointNegShifted().z()/(minRad*1e2));
  double thetas2_negpos = TMath::Pi()/2. - TMath::ATan(mut_V0_2->NominalTpcPointPosShifted().z()/(minRad*1e2));
  double etas1_negpos = -TMath::Log( TMath::Tan(thetas1_negpos/2.) );
  double etas2_negpos = -TMath::Log( TMath::Tan(thetas2_negpos/2.) );
  double detas_negpos = TMath::Abs(etas1_negpos - etas2_negpos);
  double distSft_negpos = TMath::Sqrt(TMath::Power(mut_V0_1->NominalTpcPointNegShifted().x() -
						   mut_V0_2->NominalTpcPointPosShifted().x(),2) +
				      TMath::Power(mut_V0_1->NominalTpcPointNegShifted().y() -
						   mut_V0_2->NominalTpcPointPosShifted().y(),2));
  double dPhiS_negpos = 2.0 * TMath::ATan(distSft_negpos/2./((minRad*1e2)));

  double thetas1_negneg = TMath::Pi()/2. - TMath::ATan(mut_V0_1->NominalTpcPointNegShifted().z()/(minRad*1e2));
  double thetas2_negneg = TMath::Pi()/2. - TMath::ATan(mut_V0_2->NominalTpcPointNegShifted().z()/(minRad*1e2));
  double etas1_negneg = -TMath::Log( TMath::Tan(thetas1_negneg/2.) );
  double etas2_negneg = -TMath::Log( TMath::Tan(thetas2_negneg/2.) );
  double detas_negneg = TMath::Abs(etas1_negneg - etas2_negneg);
  double distSft_negneg = TMath::Sqrt(TMath::Power(mut_V0_1->NominalTpcPointNegShifted().x() -
						   mut_V0_2->NominalTpcPointNegShifted().x(),2) +
				      TMath::Power(mut_V0_1->NominalTpcPointNegShifted().y() -
						   mut_V0_2->NominalTpcPointNegShifted().y(),2));
  double dPhiS_negneg = 2.0 * TMath::ATan(distSft_negneg/2./((minRad*1e2)));

  output_pospos->Fill(detas_pospos, dPhiS_pospos);
  output_posneg->Fill(detas_posneg, dPhiS_posneg);
  output_negpos->Fill(detas_negpos, dPhiS_negpos);
  output_negneg->Fill(detas_negneg, dPhiS_negneg);

}

//____________________________
void AliFemtoCorrFctnDPhiStarDEtaStar::AddRealPair( AliFemtoPair* pair){

  // Add real (effect) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  const AliFemtoParticle *track_1 = pair->Track1();
  const AliFemtoParticle *track_2 = pair->Track2();

  switch (fPairType) {
    // 2 tracks
  case kTracks:
    StoreDPhiStarDEtaStarBetweenTracks(track_1->Track(),
    				       track_2->Track(),
    				       fDPhiStarDEtaStarNumerator,
    				       fMinRad);
    break;
    // track + V0
  case kTrackV0:
    StoreDPhiStarDEtaStarBetweenV0AndTrack(track_1->V0(),
    					   track_2->Track(),
    					   fDPhiStarDEtaStarNumeratorTrackPos,
    					   fDPhiStarDEtaStarNumeratorTrackNeg,
    					   fMinRad);

    break;

    // 2 V0s
  case kV0s:
    StoreDPhiStarDEtaStarBetweenV0s(track_1->V0(),
    				    track_2->V0(),
    				    fDPhiStarDEtaStarNumeratorPosPos,
    				    fDPhiStarDEtaStarNumeratorPosNeg,
    				    fDPhiStarDEtaStarNumeratorNegPos,
    				    fDPhiStarDEtaStarNumeratorNegNeg,
    				    fMinRad);

    break;
  }
}

//____________________________
void AliFemtoCorrFctnDPhiStarDEtaStar::AddMixedPair( AliFemtoPair* pair){
  // Add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  const AliFemtoParticle *track_1 = pair->Track1();
  const AliFemtoParticle *track_2 = pair->Track2();

  switch (fPairType) {
    // 2 tracks
  case kTracks:
    StoreDPhiStarDEtaStarBetweenTracks(track_1->Track(),
    				       track_2->Track(),
    				       fDPhiStarDEtaStarDenominator,
    				       fMinRad);
    break;
    // track + V0
  case kTrackV0:
    StoreDPhiStarDEtaStarBetweenV0AndTrack(track_1->V0(),
    					   track_2->Track(),
    					   fDPhiStarDEtaStarDenominatorTrackPos,
    					   fDPhiStarDEtaStarDenominatorTrackNeg,
    					   fMinRad);

    break;

    // 2 V0s
  case kV0s:
    StoreDPhiStarDEtaStarBetweenV0s(track_1->V0(),
    				    track_2->V0(),
    				    fDPhiStarDEtaStarDenominatorPosPos,
    				    fDPhiStarDEtaStarDenominatorPosNeg,
    				    fDPhiStarDEtaStarDenominatorNegPos,
    				    fDPhiStarDEtaStarDenominatorNegNeg,
    				    fMinRad);

    break;
  }
}


void AliFemtoCorrFctnDPhiStarDEtaStar::WriteHistos()
{
  // Write out result histograms
  if (fPairType == kTracks) {
    fDPhiStarDEtaStarNumerator->Write();
    fDPhiStarDEtaStarDenominator->Write();
  }
  else  if (fPairType == kTrackV0) {
    fDPhiStarDEtaStarNumeratorTrackPos->Write();
    fDPhiStarDEtaStarDenominatorTrackPos->Write();
    fDPhiStarDEtaStarNumeratorTrackNeg->Write();
    fDPhiStarDEtaStarDenominatorTrackNeg->Write();
  }
  else  if (fPairType == kV0s) {
    fDPhiStarDEtaStarNumeratorPosPos->Write();
    fDPhiStarDEtaStarDenominatorPosPos->Write();
    fDPhiStarDEtaStarNumeratorPosNeg->Write();
    fDPhiStarDEtaStarDenominatorPosNeg->Write();
    fDPhiStarDEtaStarNumeratorNegPos->Write();
    fDPhiStarDEtaStarDenominatorNegPos->Write();
    fDPhiStarDEtaStarNumeratorNegNeg->Write();
    fDPhiStarDEtaStarDenominatorNegNeg->Write();
  }
}

TList* AliFemtoCorrFctnDPhiStarDEtaStar::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  if (fPairType == kTracks) {
    tOutputList->Add(fDPhiStarDEtaStarNumerator);
    tOutputList->Add(fDPhiStarDEtaStarDenominator);
  }
  else  if (fPairType == kTrackV0) {
    tOutputList->Add(fDPhiStarDEtaStarNumeratorTrackPos);
    tOutputList->Add(fDPhiStarDEtaStarDenominatorTrackPos);
    tOutputList->Add(fDPhiStarDEtaStarNumeratorTrackNeg);
    tOutputList->Add(fDPhiStarDEtaStarDenominatorTrackNeg);
  }
  else  if (fPairType == kV0s) {
    tOutputList->Add(fDPhiStarDEtaStarNumeratorPosPos);
    tOutputList->Add(fDPhiStarDEtaStarDenominatorPosPos);
    tOutputList->Add(fDPhiStarDEtaStarNumeratorPosNeg);
    tOutputList->Add(fDPhiStarDEtaStarDenominatorPosNeg);
    tOutputList->Add(fDPhiStarDEtaStarNumeratorNegPos);
    tOutputList->Add(fDPhiStarDEtaStarDenominatorNegPos);
    tOutputList->Add(fDPhiStarDEtaStarNumeratorNegNeg);
    tOutputList->Add(fDPhiStarDEtaStarDenominatorNegNeg);
  }
  return tOutputList;
}

void AliFemtoCorrFctnDPhiStarDEtaStar::SetRadius(double minrad)
{
  fMinRad = minrad;
}

void AliFemtoCorrFctnDPhiStarDEtaStar::SetPairType(AliFemtoPairType pairtype)
{
  fPairType = pairtype;
}
