////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction - correlation     //
// function for two particle correlations which uses dPhi* and k* as          //
// a function variables, calculates the average fraction of "merged" points   //
// in two tracks.                                                             //
//                                                                            //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction);
  /// \endcond
#endif
  
#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction(char* title, Double_t aRadiusMin, Double_t aRadiusMax, Double_t aDistanceMax, Double_t aDEtaMax, const int& aKStarBins, Double_t aKStarRangeLow, Double_t aKStarRangeUp, const Int_t& aPhiStarBins, Double_t aPhiStarRangeLow, Double_t aPhiStarRangeUp):
AliFemtoCorrFctn(),
  fDPhiStarKStarMergedNumerator(0),
  fDPhiStarKStarTotalNumerator(0),
  fDPhiStarKStarMergedDenominator(0),
  fDPhiStarKStarTotalDenominator(0),
  fKStarRangeLow(0),
  fKStarRangeUp(0),
  fDPhiStarRangeLow(0),
  fDPhiStarRangeUp(0),
  fDistanceMax(0.0),
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
  fDEtaMax = aDEtaMax;
  
  // Calculate radii range:
  fRadiusMin = aRadiusMin;
  fRadiusMax = aRadiusMax;

  // Set up numerator:
  char tTitMergedNum[101] = "NumDPhiStarKStarAverageMergedPoints";
  strncat(tTitMergedNum, title, 100);
  fDPhiStarKStarMergedNumerator = new TH2D(tTitMergedNum, title, aKStarBins, fKStarRangeLow, fKStarRangeUp, aPhiStarBins, fDPhiStarRangeLow, fDPhiStarRangeUp);
  char tTitTotalNum[101] = "NumDPhiStarKStarAverageTotalPoints";
  strncat(tTitTotalNum, title, 100);
  fDPhiStarKStarTotalNumerator = new TH2D(tTitTotalNum, title, aKStarBins, fKStarRangeLow, fKStarRangeUp, aPhiStarBins, fDPhiStarRangeLow, fDPhiStarRangeUp);
  // Set up denominator:
  char tTitMergedDen[101] = "DenDPhiStarKStarAverageMergedPoints";
  strncat(tTitMergedDen, title, 100);
  fDPhiStarKStarMergedDenominator = new TH2D(tTitMergedDen, title, aKStarBins, fKStarRangeLow, fKStarRangeUp, aPhiStarBins, fDPhiStarRangeLow, fDPhiStarRangeUp);
  char tTitTotalDen[101] = "DenDPhiStarKStarAverageTotalPoints";
  strncat(tTitTotalDen, title, 100);
  fDPhiStarKStarTotalDenominator = new TH2D(tTitTotalDen, title, aKStarBins, fKStarRangeLow, fKStarRangeUp, aPhiStarBins, fDPhiStarRangeLow, fDPhiStarRangeUp);

  // To enable error bars calculations:
  fDPhiStarKStarMergedNumerator->Sumw2();
  fDPhiStarKStarTotalNumerator->Sumw2();
  fDPhiStarKStarMergedDenominator->Sumw2();
  fDPhiStarKStarTotalDenominator->Sumw2();
}

//____________________________
AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction(const AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiStarKStarMergedNumerator(0),
  fDPhiStarKStarTotalNumerator(0),
  fDPhiStarKStarMergedDenominator(0),
  fDPhiStarKStarTotalDenominator(0),
  fKStarRangeLow(0),
  fKStarRangeUp(0),
  fDPhiStarRangeLow(0),
  fDPhiStarRangeUp(0),
  fDistanceMax(0.0),
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
  fDEtaMax = aCorrFctn.fDEtaMax;
  fRadiusMin = aCorrFctn.fRadiusMin;
  fRadiusMax = aCorrFctn.fRadiusMax;
  fMagSign = aCorrFctn.fMagSign;
  
}

//____________________________
AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::~AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction(){
  // Destructor
  delete fDPhiStarKStarMergedNumerator;
  delete fDPhiStarKStarMergedDenominator;
  delete fDPhiStarKStarTotalNumerator;
  delete fDPhiStarKStarTotalDenominator;
}

//_________________________
AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction& AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::operator=(const AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction& aCorrFctn)
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
  fDEtaMax = aCorrFctn.fDEtaMax;
  fRadiusMin = aCorrFctn.fRadiusMin;
  fRadiusMax = aCorrFctn.fRadiusMax;
  fMagSign = aCorrFctn.fMagSign;

  return *this;
}

//_________________________
void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::Finish(){
  // Here is where we should normalize, fit, etc...
  // We should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //  mShareDenominator->Draw();
  //  mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::Report(){
  // Create report
  string stemp = "TPC Ncls Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in merged numerator:\t%E\n",fDPhiStarKStarMergedNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in total numerator:\t%E\n",fDPhiStarKStarTotalNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in merged denominator:\t%E\n",fDPhiStarKStarMergedDenominator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in total denominator:\t%E\n",fDPhiStarKStarTotalDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}

//____________________________
void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::AddRealPair( AliFemtoPair* pair){
  // Add real (effect) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

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

  // Calculate dPhiStar for minimal radius:
  double afsi0b_radiusmin = -0.07510020733*chg1*fMagSign*fRadiusMin/pt1;
  double afsi1b_radiusmin = -0.07510020733*chg2*fMagSign*fRadiusMin/pt2;
  Double_t dphistar_radiusmin =  phi2 - phi1 + TMath::ASin(afsi1b_radiusmin) - TMath::ASin(afsi0b_radiusmin);
  dphistar_radiusmin = TVector2::Phi_mpi_pi(dphistar_radiusmin);
  
  if(TMath::Abs(deta) < TMath::Abs(fDEtaMax)) {
    
    // Iterate through all radii in range (fRadiusMin, fRadiusMax):
    for(double irad = fRadiusMin; irad < fRadiusMax; irad += 0.01) {
      
      // Calculate dPhiStar:
      double afsi0b = 0.07510020733*chg1*fMagSign*irad/pt1;
      double afsi1b = 0.07510020733*chg2*fMagSign*irad/pt2;
      Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
      dphistar = TVector2::Phi_mpi_pi(dphistar);

      // Calculate distance between two points in track:
      double distance = 2 * TMath::Sin(TMath::Abs(dphistar) * 0.5) * irad;

      // Check if pair parameters meet the requirements:
      if(distance < fDistanceMax) {
	fDPhiStarKStarMergedNumerator->Fill(kstar, dphistar_radiusmin);
      }
      fDPhiStarKStarTotalNumerator->Fill(kstar, dphistar_radiusmin);
    }
  }
}

//____________________________
void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::AddMixedPair( AliFemtoPair* pair){
  // Add real (effect) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

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

  // Calculate dPhiStar for minimal radius:
  double afsi0b_radiusmin = -0.07510020733*chg1*fMagSign*fRadiusMin/pt1;
  double afsi1b_radiusmin = -0.07510020733*chg2*fMagSign*fRadiusMin/pt2;
  Double_t dphistar_radiusmin =  phi2 - phi1 + TMath::ASin(afsi1b_radiusmin) - TMath::ASin(afsi0b_radiusmin);
  dphistar_radiusmin = TVector2::Phi_mpi_pi(dphistar_radiusmin);
  
  if(TMath::Abs(deta) < TMath::Abs(fDEtaMax)) {
    
    // Iterate through all radii in range (fRadiusMin, fRadiusMax):
    for(double irad = fRadiusMin; irad < fRadiusMax; irad += 0.01) {
      
      // Calculate dPhiStar:
      double afsi0b = 0.07510020733*chg1*fMagSign*irad/pt1;
      double afsi1b = 0.07510020733*chg2*fMagSign*irad/pt2;
      Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
      dphistar = TVector2::Phi_mpi_pi(dphistar);

      // Calculate distance between two points in track:
      double distance = 2 * TMath::Sin(TMath::Abs(dphistar) * 0.5) * irad;

      // Check if pair parameters meet the requirements:
      if(distance < fDistanceMax) {
	fDPhiStarKStarMergedDenominator->Fill(kstar, dphistar_radiusmin);
      }
      fDPhiStarKStarTotalDenominator->Fill(kstar, dphistar_radiusmin);
    }
  }
}


void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::WriteHistos()
{
  // Write out result histograms
  fDPhiStarKStarMergedNumerator->Write();
  fDPhiStarKStarTotalNumerator->Write();
  fDPhiStarKStarMergedDenominator->Write();
  fDPhiStarKStarTotalDenominator->Write();
}

TList* AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();
  
  tOutputList->Add(fDPhiStarKStarMergedNumerator);
  tOutputList->Add(fDPhiStarKStarTotalNumerator);
  tOutputList->Add(fDPhiStarKStarMergedDenominator);
  tOutputList->Add(fDPhiStarKStarTotalDenominator);
  
  return tOutputList;
}

void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::SetRadiusMax(double maxrad)
{
  fRadiusMax = maxrad;
}

void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::SetRadiusMin(double minrad)
{
  fRadiusMin = minrad;
}

void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::SetDistanceMax(double maxdist)
{
  fDistanceMax = maxdist;
}

void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::SetDEtaMax(double deta)
{
  fDEtaMax = deta;
}

void AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction::SetMagneticFieldSign(int magsign)
{
  if(magsign>1) fMagSign = 1;
  else if(magsign<1) fMagSign = -1;
  else fMagSign = magsign;
}
