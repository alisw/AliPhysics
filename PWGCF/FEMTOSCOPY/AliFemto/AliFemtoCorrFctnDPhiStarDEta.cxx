////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDPhiStarDEta - correlation function for two particle       //
// correlations which uses dPhi* and dEta as a function variables.            //
//                                                                            //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDPhiStarDEta.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__ 
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctnDPhiStarDEta);
  /// \endcond
#endif
  
#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnDPhiStarDEta::AliFemtoCorrFctnDPhiStarDEta(char* title, double radius=0.8, const int& aEtaBins=50, double aEtaRangeLow=-0.1, double aEtaRangeUp=0.1, const int& aPhiStarBins=50, double aPhiStarRangeLow=-0.1, double aPhiStarRangeUp=0.1):
  AliFemtoCorrFctn(),
  fDPhiStarDEtaNumerator(0),
  fDPhiStarDEtaDenominator(0),
  fEtaRangeLow(0),
  fEtaRangeUp(0),
  fPhiStarRangeLow(0),
  fPhiStarRangeUp(0),
  fMinRad(0.8),
  fMagSign(1)
{

  // Calculate lower and upper range of Eta and PhiStar:
  fEtaRangeLow = aEtaRangeLow;
  fEtaRangeUp = aEtaRangeUp;
  fPhiStarRangeLow = aPhiStarRangeLow;
  fPhiStarRangeUp = aPhiStarRangeUp;

  // Calculate radial distance
  fMinRad = radius;

  // Set up numerator:
  char tTitNum[101] = "NumDPhiStarDEta";
  strncat(tTitNum, title, 100);
  fDPhiStarDEtaNumerator = new TH2D(tTitNum, title, aEtaBins, fEtaRangeLow, fEtaRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);
  // Set up denominator:
  char tTitDen[101] = "DenDPhiStarDEta";
  strncat(tTitDen, title, 100);
  fDPhiStarDEtaDenominator = new TH2D(tTitDen, title, aEtaBins, fEtaRangeLow, fEtaRangeUp, aPhiStarBins, fPhiStarRangeLow, fPhiStarRangeUp);

  // To enable error bars calculations:
  fDPhiStarDEtaNumerator->Sumw2();
  fDPhiStarDEtaDenominator->Sumw2();
}

//____________________________
AliFemtoCorrFctnDPhiStarDEta::AliFemtoCorrFctnDPhiStarDEta(const AliFemtoCorrFctnDPhiStarDEta& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiStarDEtaNumerator(0),
  fDPhiStarDEtaDenominator(0),
  fEtaRangeLow(0),
  fEtaRangeUp(0),
  fPhiStarRangeLow(0),
  fPhiStarRangeUp(0),
  fMinRad(0.8),
  fMagSign(1)
{
  // Copy constructor
  if (aCorrFctn.fDPhiStarDEtaNumerator)
    fDPhiStarDEtaNumerator = new TH2D(*aCorrFctn.fDPhiStarDEtaNumerator);
  else
    fDPhiStarDEtaNumerator = 0;
  if (aCorrFctn.fDPhiStarDEtaDenominator)
    fDPhiStarDEtaDenominator = new TH2D(*aCorrFctn.fDPhiStarDEtaDenominator);
  else
    fDPhiStarDEtaDenominator = 0;

  fEtaRangeLow = aCorrFctn.fEtaRangeLow;
  fEtaRangeUp = aCorrFctn.fEtaRangeUp;
  fPhiStarRangeLow = aCorrFctn.fPhiStarRangeLow;
  fPhiStarRangeUp = aCorrFctn.fPhiStarRangeUp;
  fMinRad = aCorrFctn.fMinRad;
  fMagSign = aCorrFctn.fMagSign;
  
}

//____________________________
AliFemtoCorrFctnDPhiStarDEta::~AliFemtoCorrFctnDPhiStarDEta(){
  // Destructor
  delete fDPhiStarDEtaNumerator;
  delete fDPhiStarDEtaDenominator;
}

//_________________________
AliFemtoCorrFctnDPhiStarDEta& AliFemtoCorrFctnDPhiStarDEta::operator=(const AliFemtoCorrFctnDPhiStarDEta& aCorrFctn)
{
  // Assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fDPhiStarDEtaNumerator)
    fDPhiStarDEtaNumerator = new TH2D(*aCorrFctn.fDPhiStarDEtaNumerator);
  else
    fDPhiStarDEtaNumerator = 0;
  if (aCorrFctn.fDPhiStarDEtaDenominator)
    fDPhiStarDEtaDenominator = new TH2D(*aCorrFctn.fDPhiStarDEtaDenominator);
  else
    fDPhiStarDEtaDenominator = 0;

  fEtaRangeLow = aCorrFctn.fEtaRangeLow;
  fEtaRangeUp = aCorrFctn.fEtaRangeUp;
  fPhiStarRangeLow = aCorrFctn.fPhiStarRangeLow;
  fPhiStarRangeUp = aCorrFctn.fPhiStarRangeUp;
  fMinRad = aCorrFctn.fMinRad;
  fMagSign = aCorrFctn.fMagSign;

  return *this;
}

//_________________________
void AliFemtoCorrFctnDPhiStarDEta::Finish(){
  // Here is where we should normalize, fit, etc...
  // We should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //  mShareDenominator->Draw();
  //  mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDPhiStarDEta::Report(){
  // Create report
  string stemp = "TPC Ncls Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fDPhiStarDEtaNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDPhiStarDEtaDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}

//____________________________
void AliFemtoCorrFctnDPhiStarDEta::AddRealPair( AliFemtoPair* pair){
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

  // Calculate radial distance:
  Double_t rad;
  rad = fMinRad;

  // Calculate dPhiStar:
  double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
  double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
  Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
  
  //double dphistar = phistar1 - phistar2;
  //while (dphistar<fPhiStarRangeLow) dphistar += PIT;
  //while (dphistar>fPhiStarRangeUp) dphistar -= PIT;
 

  // Calculate dEta:
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();
  double deta = eta2 - eta1;

  // Fill numerator:
  fDPhiStarDEtaNumerator->Fill(deta, dphistar);
}

//____________________________
void AliFemtoCorrFctnDPhiStarDEta::AddMixedPair( AliFemtoPair* pair){
  // Add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  // Prepare variables:
  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double chg1 = pair->Track1()->Track()->Charge();
  double chg2 = pair->Track2()->Track()->Charge();
  double pt1 = pair->Track1()->Track()->Pt();
  double pt2 = pair->Track2()->Track()->Pt();

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

  // Calculate radial distance:
  Double_t rad;
  rad = fMinRad;

  // Calculate dPhiStar:
  double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
  double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
  Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
  
  //double dphistar = phistar1 - phistar2;
  //while (dphistar<fPhiStarRangeLow) dphistar += PIT;
  //while (dphistar>fPhiStarRangeUp) dphistar -= PIT;
 

  // Calculate dEta:
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();
  double deta = eta2 - eta1;

  // Fill denominator:
  fDPhiStarDEtaDenominator->Fill(deta, dphistar);
}


void AliFemtoCorrFctnDPhiStarDEta::WriteHistos()
{
  // Write out result histograms
  fDPhiStarDEtaNumerator->Write();
  fDPhiStarDEtaDenominator->Write();
}

TList* AliFemtoCorrFctnDPhiStarDEta::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();
  
  tOutputList->Add(fDPhiStarDEtaNumerator);
  tOutputList->Add(fDPhiStarDEtaDenominator);
  
  return tOutputList;
}

void AliFemtoCorrFctnDPhiStarDEta::SetMinimumRadius(double minrad)
{
  fMinRad = minrad;
}

void AliFemtoCorrFctnDPhiStarDEta::SetMagneticFieldSign(int magsign)
{
  if(magsign>1) fMagSign = 1;
  else if(magsign<1) fMagSign = -1;
  else fMagSign = magsign;
}
