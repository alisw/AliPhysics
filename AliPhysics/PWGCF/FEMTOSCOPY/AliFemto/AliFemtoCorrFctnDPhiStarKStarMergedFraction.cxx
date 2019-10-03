////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDPhiStarKStarMergedFraction - correlation function for two //
// particle correlations which uses dPhi* and k* as a function variables,     //
// calculates the fraction of "merged" points                                 //
//                                                                            //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDPhiStarKStarMergedFraction.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctnDPhiStarKStarMergedFraction);
  /// \endcond
#endif

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnDPhiStarKStarMergedFraction::AliFemtoCorrFctnDPhiStarKStarMergedFraction(const char* title, Double_t aRadiusMin, Double_t aRadiusMax, Double_t aDistanceMax, Double_t aMergedFractionLimit, Double_t aDEtaMax, const int& aKStarBins, Double_t aKStarRangeLow, Double_t aKStarRangeUp, const Int_t& aPhiStarBins, Double_t aPhiStarRangeLow, Double_t aPhiStarRangeUp):
AliFemtoCorrFctn(),
  fDPhiStarKStarMergedNumerator(0),
  fDPhiStarKStarTotalNumerator(0),
  fDPhiStarKStarMergedDenominator(0),
  fDPhiStarKStarTotalDenominator(0),
  fDPhiStarRangeLow(0),
  fDPhiStarRangeUp(0),
  fKStarRangeLow(0),
  fKStarRangeUp(0),
  fDistanceMax(0.0),
  fMergedFractionLimit(0.0),
  fDEtaMax(0.0),
  fRadiusMin(0.8),
  fRadiusMax(2.5),
  fMagSign(1)
{

  // Calculate lower and upper range of Eta and PhiStar:
  fKStarRangeLow = aKStarRangeLow;
  fKStarRangeUp = aKStarRangeUp;
  fDPhiStarRangeLow = aPhiStarRangeLow;
  fDPhiStarRangeUp = aPhiStarRangeUp;

  // Calculate parameters:
  fDistanceMax = aDistanceMax;
  fMergedFractionLimit = aMergedFractionLimit;
  fDEtaMax = aDEtaMax;

  // Calculate radii range
  fRadiusMin = aRadiusMin;
  fRadiusMax = aRadiusMax;

  // Set up numerator:
  char tTitMergedNum[101] = "NumDPhiStarKStarMerged";
  strncat(tTitMergedNum, title, 100);
  fDPhiStarKStarMergedNumerator = new TH2D(tTitMergedNum, title, aKStarBins, fKStarRangeLow, fKStarRangeUp, aPhiStarBins, fDPhiStarRangeLow, fDPhiStarRangeUp);
  char tTitTotalNum[101] = "NumDPhiStarKStarTotal";
  strncat(tTitTotalNum, title, 100);
  fDPhiStarKStarTotalNumerator = new TH2D(tTitTotalNum, title, aKStarBins, fKStarRangeLow, fKStarRangeUp, aPhiStarBins, fDPhiStarRangeLow, fDPhiStarRangeUp);
  // Set up denominator:
  char tTitMergedDen[101] = "DenDPhiStarKStarMerged";
  strncat(tTitMergedDen, title, 100);
  fDPhiStarKStarMergedDenominator = new TH2D(tTitMergedDen, title, aKStarBins, fKStarRangeLow, fKStarRangeUp, aPhiStarBins, fDPhiStarRangeLow, fDPhiStarRangeUp);
  char tTitTotalDen[101] = "DenDPhiStarKStarTotal";
  strncat(tTitTotalDen, title, 100);
  fDPhiStarKStarTotalDenominator = new TH2D(tTitTotalDen, title, aKStarBins, fKStarRangeLow, fKStarRangeUp, aPhiStarBins, fDPhiStarRangeLow, fDPhiStarRangeUp);

  // To enable error bars calculations:
  fDPhiStarKStarMergedNumerator->Sumw2();
  fDPhiStarKStarTotalNumerator->Sumw2();
  fDPhiStarKStarMergedDenominator->Sumw2();
  fDPhiStarKStarTotalDenominator->Sumw2();
}

//____________________________
AliFemtoCorrFctnDPhiStarKStarMergedFraction::AliFemtoCorrFctnDPhiStarKStarMergedFraction(const AliFemtoCorrFctnDPhiStarKStarMergedFraction& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiStarKStarMergedNumerator(0),
  fDPhiStarKStarTotalNumerator(0),
  fDPhiStarKStarMergedDenominator(0),
  fDPhiStarKStarTotalDenominator(0),
  fDPhiStarRangeLow(0),
  fDPhiStarRangeUp(0),
  fKStarRangeLow(0),
  fKStarRangeUp(0),
  fDistanceMax(0.0),
  fMergedFractionLimit(0.0),
  fDEtaMax(0.0),
  fRadiusMin(0.8),
  fRadiusMax(2.5),
  fMagSign(1)
{
  // Copy constructor
  if (aCorrFctn.fDPhiStarKStarMergedNumerator)
    fDPhiStarKStarMergedNumerator = new TH2D(*aCorrFctn.fDPhiStarKStarMergedNumerator);
  else
    fDPhiStarKStarMergedNumerator = 0;
  if (aCorrFctn.fDPhiStarKStarTotalNumerator)
    fDPhiStarKStarTotalNumerator = new TH2D(*aCorrFctn.fDPhiStarKStarTotalNumerator);
  else
    fDPhiStarKStarTotalNumerator = 0;
  if (aCorrFctn.fDPhiStarKStarMergedDenominator)
    fDPhiStarKStarMergedDenominator = new TH2D(*aCorrFctn.fDPhiStarKStarMergedDenominator);
  else
    fDPhiStarKStarMergedDenominator = 0;
  if (aCorrFctn.fDPhiStarKStarTotalDenominator)
    fDPhiStarKStarTotalDenominator = new TH2D(*aCorrFctn.fDPhiStarKStarTotalDenominator);
  else
    fDPhiStarKStarTotalDenominator = 0;

  fKStarRangeLow = aCorrFctn.fKStarRangeLow;
  fKStarRangeUp = aCorrFctn.fKStarRangeUp;
  fDPhiStarRangeLow = aCorrFctn.fDPhiStarRangeLow;
  fDPhiStarRangeUp = aCorrFctn.fDPhiStarRangeUp;
  fDistanceMax = aCorrFctn.fDistanceMax;
  fMergedFractionLimit = aCorrFctn.fMergedFractionLimit;
  fDEtaMax = aCorrFctn.fDEtaMax;
  fRadiusMin = aCorrFctn.fRadiusMin;
  fRadiusMax = aCorrFctn.fRadiusMax;
  fMagSign = aCorrFctn.fMagSign;

}

//____________________________
AliFemtoCorrFctnDPhiStarKStarMergedFraction::~AliFemtoCorrFctnDPhiStarKStarMergedFraction(){
  // Destructor
  delete fDPhiStarKStarMergedNumerator;
  delete fDPhiStarKStarMergedDenominator;
  delete fDPhiStarKStarTotalNumerator;
  delete fDPhiStarKStarTotalDenominator;
}

//_________________________
AliFemtoCorrFctnDPhiStarKStarMergedFraction& AliFemtoCorrFctnDPhiStarKStarMergedFraction::operator=(const AliFemtoCorrFctnDPhiStarKStarMergedFraction& aCorrFctn)
{
  // Assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fDPhiStarKStarMergedNumerator)
    fDPhiStarKStarMergedNumerator = new TH2D(*aCorrFctn.fDPhiStarKStarMergedNumerator);
  else
    fDPhiStarKStarMergedNumerator = 0;
  if (aCorrFctn.fDPhiStarKStarTotalNumerator)
    fDPhiStarKStarTotalNumerator = new TH2D(*aCorrFctn.fDPhiStarKStarTotalNumerator);
  else
    fDPhiStarKStarTotalNumerator = 0;
  if (aCorrFctn.fDPhiStarKStarMergedDenominator)
    fDPhiStarKStarMergedDenominator = new TH2D(*aCorrFctn.fDPhiStarKStarMergedDenominator);
  else
    fDPhiStarKStarMergedDenominator = 0;
  if (aCorrFctn.fDPhiStarKStarTotalDenominator)
    fDPhiStarKStarTotalDenominator = new TH2D(*aCorrFctn.fDPhiStarKStarTotalDenominator);
  else
    fDPhiStarKStarTotalDenominator = 0;

  fKStarRangeLow = aCorrFctn.fKStarRangeLow;
  fKStarRangeUp = aCorrFctn.fKStarRangeUp;
  fDPhiStarRangeLow = aCorrFctn.fDPhiStarRangeLow;
  fDPhiStarRangeUp = aCorrFctn.fDPhiStarRangeUp;
  fDistanceMax = aCorrFctn.fDistanceMax;
  fMergedFractionLimit = aCorrFctn.fMergedFractionLimit;
  fDEtaMax = aCorrFctn.fDEtaMax;
  fRadiusMin = aCorrFctn.fRadiusMin;
  fRadiusMax = aCorrFctn.fRadiusMax;
  fMagSign = aCorrFctn.fMagSign;

  return *this;
}

//_________________________
void AliFemtoCorrFctnDPhiStarKStarMergedFraction::Finish(){
  // Here is where we should normalize, fit, etc...
  // We should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //  mShareDenominator->Draw();
  //  mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDPhiStarKStarMergedFraction::Report()
{
  // Create report
  AliFemtoString report = "DPhiStarKStarMergedFraction Correlation Function Report:\n";
  report += Form("Number of entries in merged numerator:\t%E\n",fDPhiStarKStarMergedNumerator->GetEntries());
  report += Form("Number of entries in total numerator:\t%E\n",fDPhiStarKStarTotalNumerator->GetEntries());
  report += Form("Number of entries in merged denominator:\t%E\n",fDPhiStarKStarMergedDenominator->GetEntries());
  report += Form("Number of entries in total denominator:\t%E\n",fDPhiStarKStarTotalDenominator->GetEntries());
  //  report += mCoulombWeight->Report();

  return report;
}

//____________________________
void AliFemtoCorrFctnDPhiStarKStarMergedFraction::AddRealPair( AliFemtoPair* pair){
  // Add real (effect) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  // Prepare variables:
  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double chg1 = pair->Track1()->Track()->Charge();
  double chg2 = pair->Track2()->Track()->Charge();
  double pt1 = pair->Track1()->Track()->Pt();
  double pt2 = pair->Track2()->Track()->Pt();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();
  double kstar = pair->KStar();

  // Check magnetic field sign:
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  Double_t magsign = 0.0;
  if (!aodH) {
    return;
  }
  else {
    AliAODEvent *fAOD;
    fAOD = aodH->GetEvent();
    magsign = fAOD->GetMagneticField();
  }
  if (magsign > 1)
    fMagSign = 1;
  else if ( magsign < 1)
    fMagSign = -1;
  else
    fMagSign = magsign;

  // Calculate dEta:
  double deta = eta2 - eta1;

  if(TMath::Abs(deta) < TMath::Abs(fDEtaMax)) {

    Double_t badpoints = 0.0;
    Double_t allpoints = 0.0;

    // Iterate through all radii in range (fRadiusMin, fRadiusMax):
    for(double irad = fRadiusMin; irad < fRadiusMax; irad += 0.01) {

      // Calculate dPhiStar:
      double afsi0b = -0.07510020733*chg1*fMagSign*irad/pt1;
      double afsi1b = -0.07510020733*chg2*fMagSign*irad/pt2;
      Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
      dphistar = TVector2::Phi_mpi_pi(dphistar);

      // Calculate distance between two points in track:
      double distance = 2 * TMath::Sin(TMath::Abs(dphistar) * 0.5) * irad;

      // Check if pair parameters meet the requirements:
      if(distance < fDistanceMax) {
	badpoints += 1.0;
      }
      allpoints += 1.0;

    }

    if(allpoints != 0.0) {
      // Calculate fraction:
      Double_t fraction = badpoints / allpoints;

      // Add pair if the fraction is above limit:
      if(fraction > fMergedFractionLimit) {
	double rad = fRadiusMin;
	double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
	double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
	Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
	fDPhiStarKStarMergedNumerator->Fill(kstar, dphistar);
      }
    }
  }

  double rad = fRadiusMin;
  double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
  double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
  Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
  fDPhiStarKStarTotalNumerator->Fill(kstar, dphistar);
}

//____________________________
void AliFemtoCorrFctnDPhiStarKStarMergedFraction::AddMixedPair( AliFemtoPair* pair){
  // Add real (effect) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  // Prepare variables:
  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double chg1 = pair->Track1()->Track()->Charge();
  double chg2 = pair->Track2()->Track()->Charge();
  double pt1 = pair->Track1()->Track()->Pt();
  double pt2 = pair->Track2()->Track()->Pt();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();
  double kstar = pair->KStar();

  // Check magnetic field sign:
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  Double_t magsign = 0.0;
  if (!aodH) {
    return;
  }
  else {
    AliAODEvent *fAOD;
    fAOD = aodH->GetEvent();
    magsign = fAOD->GetMagneticField();
  }
  if (magsign > 1)
    fMagSign = 1;
  else if ( magsign < 1)
    fMagSign = -1;
  else
    fMagSign = magsign;

  // Calculate dEta:
  double deta = eta2 - eta1;

  if(TMath::Abs(deta) < TMath::Abs(fDEtaMax)) {

    Double_t badpoints = 0.0;
    Double_t allpoints = 0.0;

    // Iterate through all radii in range (fRadiusMin, fRadiusMax):
    for(double irad = fRadiusMin; irad < fRadiusMax; irad += 0.01) {

      // Calculate dPhiStar:
      double afsi0b = -0.07510020733*chg1*fMagSign*irad/pt1;
      double afsi1b = -0.07510020733*chg2*fMagSign*irad/pt2;
      Double_t dphistar = phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
      dphistar = TVector2::Phi_mpi_pi(dphistar);

      // Calculate distance between two points in track:
      double distance = 2 * TMath::Sin(TMath::Abs(dphistar) * 0.5) * irad;

      // Check if pair parameters meet the requirements:
      if(distance < fDistanceMax) {
	badpoints += 1.0;
      }
      allpoints += 1.0;

    }

    if(allpoints != 0.0) {
      // Calculate fraction:
      Double_t fraction = badpoints / allpoints;

      // Add pair if the fraction is above limit:
      if(fraction > fMergedFractionLimit) {
	double rad = fRadiusMin;
	double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
	double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
	Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
	fDPhiStarKStarMergedDenominator->Fill(kstar, dphistar);
      }
    }
  }

  double rad = fRadiusMin;
  double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
  double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
  Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
  fDPhiStarKStarTotalDenominator->Fill(kstar, dphistar);
}


void AliFemtoCorrFctnDPhiStarKStarMergedFraction::WriteHistos()
{
  // Write out result histograms
  fDPhiStarKStarMergedNumerator->Write();
  fDPhiStarKStarTotalNumerator->Write();
  fDPhiStarKStarMergedDenominator->Write();
  fDPhiStarKStarTotalDenominator->Write();
}

TList* AliFemtoCorrFctnDPhiStarKStarMergedFraction::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiStarKStarMergedNumerator);
  tOutputList->Add(fDPhiStarKStarTotalNumerator);
  tOutputList->Add(fDPhiStarKStarMergedDenominator);
  tOutputList->Add(fDPhiStarKStarTotalDenominator);

  return tOutputList;
}

void AliFemtoCorrFctnDPhiStarKStarMergedFraction::SetRadiusMax(double maxrad)
{
  fRadiusMax = maxrad;
}

void AliFemtoCorrFctnDPhiStarKStarMergedFraction::SetRadiusMin(double minrad)
{
  fRadiusMin = minrad;
}

void AliFemtoCorrFctnDPhiStarKStarMergedFraction::SetDistanceMax(double maxdist)
{
  fDistanceMax = maxdist;
}

void AliFemtoCorrFctnDPhiStarKStarMergedFraction::SetMergedFractionLimit(double frac)
{
  fMergedFractionLimit = frac;
}

void AliFemtoCorrFctnDPhiStarKStarMergedFraction::SetDEtaMax(double deta)
{
  fDEtaMax = deta;
}

void AliFemtoCorrFctnDPhiStarKStarMergedFraction::SetMagneticFieldSign(int magsign)
{
  if(magsign>1) fMagSign = 1;
  else if(magsign<1) fMagSign = -1;
  else fMagSign = magsign;
}
