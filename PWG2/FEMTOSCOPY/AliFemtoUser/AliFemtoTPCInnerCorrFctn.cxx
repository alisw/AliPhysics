////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoTPCInnerCorrFctn - A correlation function that saves the         ///
/// distance at the entrance to the TPC between two tracks as a function     ///
/// of qinv                                                                  ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoTPCInnerCorrFctn.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoTPCInnerCorrFctn)
#endif

//____________________________
AliFemtoTPCInnerCorrFctn::AliFemtoTPCInnerCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  fDTPCNumerator(0),
  fDTPCDenominator(0),
  fRadDNumerator(0),
  fRadDDenominator(0),
  fRadius(100)
{
  // set up numerator
  //  title = "Num Qinv (MeV/c)";
  char tTitNum[101] = "NumDTPC";
  strncat(tTitNum,title, 100);
  fDTPCNumerator = new TH2D(tTitNum,title,nbins,QinvLo,QinvHi,100,0.0,20.0);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTitDen[101] = "DenDTPC";
  strncat(tTitDen,title, 100);
  fDTPCDenominator = new TH2D(tTitDen,title,nbins,QinvLo,QinvHi,100,0.0,20.0);

  char tTitNumR[101] = "NumRadD";
  strncat(tTitNumR,title, 100);
  fRadDNumerator = new TH2D(tTitNumR,title,50,-0.1,0.1,50,-0.1,0.1);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTitDenR[101] = "DenRadD";
  strncat(tTitDenR,title, 100);
  fRadDDenominator = new TH2D(tTitDenR,title,50,-0.1,0.1,50,-0.1,0.1);

  // to enable error bar calculation...
  fDTPCNumerator->Sumw2();
  fDTPCDenominator->Sumw2();
  fRadDNumerator->Sumw2();
  fRadDDenominator->Sumw2();
}

//____________________________
AliFemtoTPCInnerCorrFctn::AliFemtoTPCInnerCorrFctn(const AliFemtoTPCInnerCorrFctn& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDTPCNumerator(0),
  fDTPCDenominator(0),
  fRadDNumerator(0),
  fRadDDenominator(0),
  fRadius(100)
{
  // copy constructor
  if (aCorrFctn.fDTPCNumerator)
    fDTPCNumerator = new TH2D(*aCorrFctn.fDTPCNumerator);
  if (aCorrFctn.fDTPCDenominator)
    fDTPCDenominator = new TH2D(*aCorrFctn.fDTPCDenominator);
  if (aCorrFctn.fRadDNumerator)
    fRadDNumerator = new TH2D(*aCorrFctn.fRadDNumerator);
  if (aCorrFctn.fRadDDenominator)
    fRadDDenominator = new TH2D(*aCorrFctn.fRadDDenominator);
}
//____________________________
AliFemtoTPCInnerCorrFctn::~AliFemtoTPCInnerCorrFctn(){
  // destructor
  delete fDTPCNumerator;
  delete fDTPCDenominator;
  delete fRadDNumerator;
  delete fRadDDenominator;
}
//_________________________
AliFemtoTPCInnerCorrFctn& AliFemtoTPCInnerCorrFctn::operator=(const AliFemtoTPCInnerCorrFctn& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fDTPCNumerator)
    fDTPCNumerator = new TH2D(*aCorrFctn.fDTPCNumerator);
  else
    fDTPCNumerator = 0;
  if (aCorrFctn.fDTPCDenominator)
    fDTPCDenominator = new TH2D(*aCorrFctn.fDTPCDenominator);
  else
    fDTPCDenominator = 0;

  if (aCorrFctn.fRadDNumerator)
    fRadDNumerator = new TH2D(*aCorrFctn.fRadDNumerator);
  else
    fRadDNumerator = 0;
  if (aCorrFctn.fRadDDenominator)
    fRadDDenominator = new TH2D(*aCorrFctn.fRadDDenominator);
  else
    fRadDDenominator = 0;

  fRadius = aCorrFctn.fRadius;

  return *this;
}
//_________________________
void AliFemtoTPCInnerCorrFctn::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoTPCInnerCorrFctn::Report(){
  // create report
  string stemp = "Entrace TPC distance Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fDTPCNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDTPCDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoTPCInnerCorrFctn::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double pih = TMath::Pi();
  double pit = TMath::Pi()*2;

  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
  double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
  double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
  double dist = sqrt(distx*distx + disty*disty + distz*distz);

  fDTPCNumerator->Fill(tQinv, dist);

  if (tQinv < 0.1) {
    double phi1 = pair->Track1()->Track()->P().Phi();
    double phi2 = pair->Track2()->Track()->P().Phi();
    double chg1 = pair->Track1()->Track()->Charge();
    double chg2 = pair->Track2()->Track()->Charge();
    double ptv1 = pair->Track1()->Track()->Pt();
    double ptv2 = pair->Track2()->Track()->Pt();
    double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
    double eta2 = pair->Track2()->Track()->P().PseudoRapidity();
    double arg1 = -0.3 * 0.5 * chg1 * fRadius/(2*ptv1);
    double arg2 = -0.3 * 0.5 * chg2 * fRadius/(2*ptv2);
    double phid = phi2 - phi1 + TMath::ASin(arg2) - TMath::ASin(arg1);

    while (phid>pih) phid -= pit;
    while (phid<-pih) phid += pit;
    //    dist = phi2 - phi1 + TMath::ASin(-0.3 * 0.5 * chg2 * fRadius/(2*ptv2)) - TMath::ASin(-0.3 * 0.5 * chg1 * fRadius/(2*ptv1));
    double etad = eta2 - eta1;
    fRadDNumerator->Fill(phid, etad);
  }
}
//____________________________
void AliFemtoTPCInnerCorrFctn::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double pih = TMath::Pi();
  double pit = TMath::Pi()*2;

  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
  double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
  double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
  double dist = sqrt(distx*distx + disty*disty + distz*distz);

  fDTPCDenominator->Fill(tQinv,dist);

  if (tQinv < 0.1) {
    double phi1 = pair->Track1()->Track()->P().Phi();
    double phi2 = pair->Track2()->Track()->P().Phi();
    double chg1 = pair->Track1()->Track()->Charge();
    double chg2 = pair->Track2()->Track()->Charge();
    double ptv1 = pair->Track1()->Track()->Pt();
    double ptv2 = pair->Track2()->Track()->Pt();
    double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
    double eta2 = pair->Track2()->Track()->P().PseudoRapidity();
    double arg1 = -0.3 * 0.5 * chg1 * fRadius/(2*ptv1);
    double arg2 = -0.3 * 0.5 * chg2 * fRadius/(2*ptv2);
    double phid = phi2 - phi1 + TMath::ASin(arg2) - TMath::ASin(arg1);

    while (phid>pih) phid -= pit;
    while (phid<-pih) phid += pit;
    //    dist = phi2 - phi1 + TMath::ASin(-0.3 * 0.5 * chg2 * fRadius/(2*ptv2)) - TMath::ASin(-0.3 * 0.5 * chg1 * fRadius/(2*ptv1));
    double etad = eta2 - eta1;
    fRadDDenominator->Fill(phid, etad);
  }
}


void AliFemtoTPCInnerCorrFctn::WriteHistos()
{
  // Write out result histograms
  fDTPCNumerator->Write();
  fDTPCDenominator->Write();
  fRadDNumerator->Write();
  fRadDDenominator->Write();
}
//______________________________
TList* AliFemtoTPCInnerCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDTPCNumerator); 
  tOutputList->Add(fDTPCDenominator);  
  tOutputList->Add(fRadDNumerator); 
  tOutputList->Add(fRadDDenominator);  

  return tOutputList;
}

void AliFemtoTPCInnerCorrFctn::SetRadius(double rad)
{
  fRadius = rad;
}
