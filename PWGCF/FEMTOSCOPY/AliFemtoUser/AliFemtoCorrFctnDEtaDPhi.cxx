////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDEtaDPhi.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctnDEtaDPhi)
#endif
  
#define PIH 1.57079632679489656
#define PIT 6.28318530717958623
#define PIQ 4.71238898038468967

//____________________________
AliFemtoCorrFctnDEtaDPhi::AliFemtoCorrFctnDEtaDPhi(char* title, const int& aPhiBins=20, const int& aEtaBins=20):
  AliFemtoCorrFctn(),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
  fDoPtAnalysis(0),
  fDPhiPtNumerator(0),
  fDPhiPtDenominator(0),
  fDCosPtNumerator(0),
  fDCosPtDenominator(0)
{
  // set up numerator
  char tTitNumD[101] = "NumDPhiDEta";
  strncat(tTitNumD,title, 100);
  fDPhiDEtaNumerator = new TH2D(tTitNumD,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.5,2.5);
  // set up denominator
  char tTitDenD[101] = "DenDPhiDEta";
  strncat(tTitDenD,title, 100);
  fDPhiDEtaDenominator = new TH2D(tTitDenD,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.5,2.5);

  // set up numerator
  char tTitNumDPhi[101] = "NumDPhi";
  strncat(tTitNumDPhi,title, 100);
  fDPhiNumerator = new TH1D(tTitNumDPhi,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());
  // set up denominator
  char tTitDenDPhi[101] = "DenDPhi";
  strncat(tTitDenDPhi,title, 100);
  fDPhiDenominator = new TH1D(tTitDenDPhi,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());

  // set up numerator
  char tTitNumDCos[101] = "NumDCos";
  strncat(tTitNumDCos,title, 100);
  fDCosNumerator = new TH1D(tTitNumDCos,title,aPhiBins*2,-1.0,1.0);
  // set up denominator
  char tTitDenDCos[101] = "DenDCos";
  strncat(tTitDenDCos,title, 100);
  fDCosDenominator = new TH1D(tTitDenDCos,title,aPhiBins*2,-1.0,1.0);

  // to enable error bar calculation...
  fDPhiDEtaNumerator->Sumw2();
  fDPhiDEtaDenominator->Sumw2();
  fDPhiNumerator->Sumw2();
  fDPhiDenominator->Sumw2();
  fDCosNumerator->Sumw2();
  fDCosDenominator->Sumw2();
}

//____________________________
AliFemtoCorrFctnDEtaDPhi::AliFemtoCorrFctnDEtaDPhi(const AliFemtoCorrFctnDEtaDPhi& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
  fDoPtAnalysis(0),
  fDPhiPtNumerator(0),
  fDPhiPtDenominator(0),
  fDCosPtNumerator(0),
  fDCosPtDenominator(0)
{
  // copy constructor
  if (aCorrFctn.fDPhiDEtaNumerator)
    fDPhiDEtaNumerator = new TH2D(*aCorrFctn.fDPhiDEtaNumerator);
  else
    fDPhiDEtaNumerator = 0;
  if (aCorrFctn.fDPhiDEtaDenominator)
    fDPhiDEtaDenominator = new TH2D(*aCorrFctn.fDPhiDEtaDenominator);
  else
    fDPhiDEtaDenominator = 0;

  if (aCorrFctn.fDPhiNumerator)
    fDPhiNumerator = new TH1D(*aCorrFctn.fDPhiNumerator);
  else
    fDPhiNumerator = 0;
  if (aCorrFctn.fDPhiDenominator)
    fDPhiDenominator = new TH1D(*aCorrFctn.fDPhiDenominator);
  else
    fDPhiDenominator = 0;

  if (aCorrFctn.fDCosNumerator)
    fDCosNumerator = new TH1D(*aCorrFctn.fDCosNumerator);
  else
    fDCosNumerator = 0;
  if (aCorrFctn.fDCosDenominator)
    fDCosDenominator = new TH1D(*aCorrFctn.fDCosDenominator);
  else
    fDCosDenominator = 0;

  if (aCorrFctn.fDPhiPtNumerator)
    fDPhiPtNumerator = new TH2D(*aCorrFctn.fDPhiPtNumerator);
  else
    fDPhiPtNumerator = 0;
  if (aCorrFctn.fDPhiPtDenominator)
    fDPhiPtDenominator = new TH2D(*aCorrFctn.fDPhiPtDenominator);
  else
    fDPhiPtDenominator = 0;

  if (aCorrFctn.fDCosPtNumerator)
    fDCosPtNumerator = new TH2D(*aCorrFctn.fDCosPtNumerator);
  else
    fDCosPtNumerator = 0;
  if (aCorrFctn.fDCosPtDenominator)
    fDCosPtDenominator = new TH2D(*aCorrFctn.fDCosPtDenominator);
  else
    fDCosPtDenominator = 0;

}
//____________________________
AliFemtoCorrFctnDEtaDPhi::~AliFemtoCorrFctnDEtaDPhi(){
  // destructor
  delete fDPhiDEtaNumerator;
  delete fDPhiDEtaDenominator;
  delete fDPhiNumerator;
  delete fDPhiDenominator;
  delete fDCosNumerator;
  delete fDCosDenominator;
  if (fDoPtAnalysis) {
    delete fDPhiPtNumerator;
    delete fDPhiPtDenominator;
    delete fDCosPtNumerator;
    delete fDCosPtDenominator;
  }
}
//_________________________
AliFemtoCorrFctnDEtaDPhi& AliFemtoCorrFctnDEtaDPhi::operator=(const AliFemtoCorrFctnDEtaDPhi& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fDPhiDEtaNumerator)
    fDPhiDEtaNumerator = new TH2D(*aCorrFctn.fDPhiDEtaNumerator);
  else
    fDPhiDEtaNumerator = 0;
  if (aCorrFctn.fDPhiDEtaDenominator)
    fDPhiDEtaDenominator = new TH2D(*aCorrFctn.fDPhiDEtaDenominator);
  else
    fDPhiDEtaDenominator = 0;

  if (aCorrFctn.fDPhiNumerator)
    fDPhiNumerator = new TH1D(*aCorrFctn.fDPhiNumerator);
  else
    fDPhiNumerator = 0;
  if (aCorrFctn.fDPhiDenominator)
    fDPhiDenominator = new TH1D(*aCorrFctn.fDPhiDenominator);
  else
    fDPhiDenominator = 0;

  if (aCorrFctn.fDCosNumerator)
    fDCosNumerator = new TH1D(*aCorrFctn.fDCosNumerator);
  else
    fDCosNumerator = 0;
  if (aCorrFctn.fDCosDenominator)
    fDCosDenominator = new TH1D(*aCorrFctn.fDCosDenominator);
  else
    fDCosDenominator = 0;

  if (aCorrFctn.fDPhiPtNumerator)
    fDPhiPtNumerator = new TH2D(*aCorrFctn.fDPhiPtNumerator);
  else
    fDPhiPtNumerator = 0;
  if (aCorrFctn.fDPhiPtDenominator)
    fDPhiPtDenominator = new TH2D(*aCorrFctn.fDPhiPtDenominator);
  else
    fDPhiPtDenominator = 0;

  if (aCorrFctn.fDCosPtNumerator)
    fDCosPtNumerator = new TH2D(*aCorrFctn.fDCosPtNumerator);
  else
    fDCosPtNumerator = 0;
  if (aCorrFctn.fDCosPtDenominator)
    fDCosPtDenominator = new TH2D(*aCorrFctn.fDCosPtDenominator);
  else
    fDCosPtDenominator = 0;

  return *this;
}
//_________________________
void AliFemtoCorrFctnDEtaDPhi::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDEtaDPhi::Report(){
  // create report
  string stemp = "TPC Ncls Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fDPhiDEtaNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDPhiDEtaDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnDEtaDPhi::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

  double dphi = phi1 - phi2;
  while (dphi<-PIH) dphi+=PIT;
  while (dphi>PIQ) dphi-=PIT;

  double deta = eta1 - eta2;

//   double px1 = pair->Track1()->Track()->P().x();
//   double py1 = pair->Track1()->Track()->P().y();
//   double pz1 = pair->Track1()->Track()->P().z();

//   double px2 = pair->Track2()->Track()->P().x();
//   double py2 = pair->Track2()->Track()->P().y();
//   double pz2 = pair->Track2()->Track()->P().z();

//   double pt1 = TMath::Hypot(px1, py1);
//   double pt2 = TMath::Hypot(px2, py2);
//   double ptmin = pt1>pt2 ? pt2 : pt1;

//   double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
//     sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));

  fDPhiDEtaNumerator->Fill(dphi, deta);

  fDPhiNumerator->Fill(dphi);
//   fDCosNumerator->Fill(cosphi);

  if (fDoPtAnalysis) {
//     fDPhiPtNumerator->Fill(dphi, ptmin);
//     fDCosPtNumerator->Fill(cosphi, ptmin);
  }

}
//____________________________
void AliFemtoCorrFctnDEtaDPhi::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

  double dphi = phi1 - phi2;
  while (dphi<-PIH) dphi+=PIT;
  while (dphi>PIQ) dphi-=PIT;

  double deta = eta1 - eta2;

//   double px1 = pair->Track1()->Track()->P().x();
//   double py1 = pair->Track1()->Track()->P().y();
//   double pz1 = pair->Track1()->Track()->P().z();

//   double px2 = pair->Track2()->Track()->P().x();
//   double py2 = pair->Track2()->Track()->P().y();
//   double pz2 = pair->Track2()->Track()->P().z();

//   double pt1 = TMath::Hypot(px1, py1);
//   double pt2 = TMath::Hypot(px2, py2);
//   double ptmin = pt1>pt2 ? pt2 : pt1;

//   double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
//     sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));

  fDPhiDEtaDenominator->Fill(dphi, deta);

  fDPhiDenominator->Fill(dphi);
//   fDCosDenominator->Fill(cosphi);

  if (fDoPtAnalysis) {
    //   fDPhiPtDenominator->Fill(dphi, ptmin);
    //   fDCosPtDenominator->Fill(cosphi, ptmin);
  }
}


void AliFemtoCorrFctnDEtaDPhi::WriteHistos()
{
  // Write out result histograms
  fDPhiDEtaNumerator->Write();
  fDPhiDEtaDenominator->Write();
  fDPhiNumerator->Write();
  fDPhiDenominator->Write();
  fDCosNumerator->Write();
  fDCosDenominator->Write();
  if (fDoPtAnalysis) {
    fDPhiPtNumerator->Write();
    fDPhiPtDenominator->Write();
    fDCosPtNumerator->Write();
    fDCosPtDenominator->Write();
  }
}

TList* AliFemtoCorrFctnDEtaDPhi::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDEtaNumerator);
  tOutputList->Add(fDPhiDEtaDenominator);
  tOutputList->Add(fDPhiNumerator);
  tOutputList->Add(fDPhiDenominator);
  tOutputList->Add(fDCosNumerator);
  tOutputList->Add(fDCosDenominator);
  if (fDoPtAnalysis) {
    tOutputList->Add(fDPhiPtNumerator);
    tOutputList->Add(fDPhiPtDenominator);
    tOutputList->Add(fDCosPtNumerator);
    tOutputList->Add(fDCosPtDenominator);
  }

  return tOutputList;

}

void AliFemtoCorrFctnDEtaDPhi::SetDoPtAnalysis(int do2d)
{
  fDoPtAnalysis = do2d;
  
  int aPhiBins = fDPhiDEtaNumerator->GetNbinsX();
  const char *title = fDPhiDEtaNumerator->GetTitle();

  // set up numerator
  char tTitNumDPhiPt[101] = "NumDPhiPt";
  strncat(tTitNumDPhiPt,title, 100);
  fDPhiPtNumerator = new TH2D(tTitNumDPhiPt,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi(), 30, 0.0, 3.0);
  // set up denominator
  char tTitDenDPhiPt[101] = "DenDPhiPt";
  strncat(tTitDenDPhiPt,title, 100);
  fDPhiPtDenominator = new TH2D(tTitDenDPhiPt,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi(), 30, 0.0, 3.0);

  // set up numerator
  char tTitNumDCosPt[101] = "NumDCosPt";
  strncat(tTitNumDCosPt,title, 100);
  fDCosPtNumerator = new TH2D(tTitNumDCosPt,title,aPhiBins*2,-1.0,1.0, 30, 0.0, 3.0);
  // set up denominator
  char tTitDenDCosPt[101] = "DenDCosPt";
  strncat(tTitDenDCosPt,title, 100);
  fDCosPtDenominator = new TH2D(tTitDenDCosPt,title,aPhiBins*2,-1.0,1.0, 30, 0.0, 3.0);

  fDPhiPtNumerator->Sumw2();
  fDPhiPtDenominator->Sumw2();
  fDCosPtNumerator->Sumw2();
  fDCosPtDenominator->Sumw2();
  
}
