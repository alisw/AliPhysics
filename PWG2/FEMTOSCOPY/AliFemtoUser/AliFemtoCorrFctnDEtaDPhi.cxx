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
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctnDEtaDPhi)
#endif

//____________________________
AliFemtoCorrFctnDEtaDPhi::AliFemtoCorrFctnDEtaDPhi(char* title, const int& aPhiBins=20, const int& aEtaBins=20):
  AliFemtoCorrFctn(),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiDEtaColNumerator(0),
  fDPhiDEtaColDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
  fDPhiPtNumerator(0),
  fDPhiPtDenominator(0),
  fDCosPtNumerator(0),
  fDCosPtDenominator(0)
{
  // set up numerator
  char tTitNumD[100] = "NumDPhiDEta";
  strcat(tTitNumD,title);
  fDPhiDEtaNumerator = new TH2D(tTitNumD,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);
  // set up denominator
  char tTitDenD[100] = "DenDPhiDEta";
  strcat(tTitDenD,title);
  fDPhiDEtaDenominator = new TH2D(tTitDenD,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);

  // set up numerator
  char tTitNumR[100] = "NumDPhiDEtaCol";
  strcat(tTitNumR,title);
  fDPhiDEtaColNumerator = new TH2D(tTitNumR,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);
  // set up denominator
  char tTitDenR[100] = "DenDPhiDEtaCol";
  strcat(tTitDenR,title);
  fDPhiDEtaColDenominator = new TH2D(tTitDenR,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);

  // set up numerator
  char tTitNumDPhi[100] = "NumDPhi";
  strcat(tTitNumDPhi,title);
  fDPhiNumerator = new TH1D(tTitNumDPhi,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());
  // set up denominator
  char tTitDenDPhi[100] = "DenDPhi";
  strcat(tTitDenDPhi,title);
  fDPhiDenominator = new TH1D(tTitDenDPhi,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());

  // set up numerator
  char tTitNumDCos[100] = "NumDCos";
  strcat(tTitNumDCos,title);
  fDCosNumerator = new TH1D(tTitNumDCos,title,aPhiBins*2,-1.0,1.0);
  // set up denominator
  char tTitDenDCos[100] = "DenDCos";
  strcat(tTitDenDCos,title);
  fDCosDenominator = new TH1D(tTitDenDCos,title,aPhiBins*2,-1.0,1.0);

  // set up numerator
  char tTitNumDPhiPt[100] = "NumDPhiPt";
  strcat(tTitNumDPhiPt,title);
  fDPhiPtNumerator = new TH2D(tTitNumDPhiPt,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi(), 30, 0.0, 3.0);
  // set up denominator
  char tTitDenDPhiPt[100] = "DenDPhiPt";
  strcat(tTitDenDPhiPt,title);
  fDPhiPtDenominator = new TH2D(tTitDenDPhiPt,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi(), 30, 0.0, 3.0);

  // set up numerator
  char tTitNumDCosPt[100] = "NumDCosPt";
  strcat(tTitNumDCosPt,title);
  fDCosPtNumerator = new TH2D(tTitNumDCosPt,title,aPhiBins*2,-1.0,1.0, 30, 0.0, 3.0);
  // set up denominator
  char tTitDenDCosPt[100] = "DenDCosPt";
  strcat(tTitDenDCosPt,title);
  fDCosPtDenominator = new TH2D(tTitDenDCosPt,title,aPhiBins*2,-1.0,1.0, 30, 0.0, 3.0);

  // to enable error bar calculation...
  fDPhiDEtaNumerator->Sumw2();
  fDPhiDEtaDenominator->Sumw2();
  fDPhiDEtaColNumerator->Sumw2();
  fDPhiDEtaColDenominator->Sumw2();
  fDPhiNumerator->Sumw2();
  fDPhiDenominator->Sumw2();
  fDCosNumerator->Sumw2();
  fDCosDenominator->Sumw2();
  fDPhiPtNumerator->Sumw2();
  fDPhiPtDenominator->Sumw2();
  fDCosPtNumerator->Sumw2();
  fDCosPtDenominator->Sumw2();

}

//____________________________
AliFemtoCorrFctnDEtaDPhi::AliFemtoCorrFctnDEtaDPhi(const AliFemtoCorrFctnDEtaDPhi& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiDEtaColNumerator(0),
  fDPhiDEtaColDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
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

  if (aCorrFctn.fDPhiDEtaColNumerator)
    fDPhiDEtaColNumerator = new TH2D(*aCorrFctn.fDPhiDEtaColNumerator);
  else
    fDPhiDEtaColNumerator = 0;
  if (aCorrFctn.fDPhiDEtaColDenominator)
    fDPhiDEtaColDenominator = new TH2D(*aCorrFctn.fDPhiDEtaColDenominator);
  else
    fDPhiDEtaColDenominator = 0;

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
  delete fDPhiDEtaColNumerator;
  delete fDPhiDEtaColDenominator;
  delete fDPhiNumerator;
  delete fDPhiDenominator;
  delete fDCosNumerator;
  delete fDCosDenominator;
  delete fDPhiPtNumerator;
  delete fDPhiPtDenominator;
  delete fDCosPtNumerator;
  delete fDCosPtDenominator;
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

  if (aCorrFctn.fDPhiDEtaColNumerator)
    fDPhiDEtaColNumerator = new TH2D(*aCorrFctn.fDPhiDEtaColNumerator);
  else
    fDPhiDEtaColNumerator = 0;
  if (aCorrFctn.fDPhiDEtaColDenominator)
    fDPhiDEtaColDenominator = new TH2D(*aCorrFctn.fDPhiDEtaColDenominator);
  else
    fDPhiDEtaColDenominator = 0;

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
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fDPhiDEtaNumerator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDPhiDEtaDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnDEtaDPhi::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  double phi1 = pair->Track1()->Track()->P().phi();
  double phi2 = pair->Track2()->Track()->P().phi();
  double eta1 = pair->Track1()->Track()->P().pseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().pseudoRapidity();

  double dphi = phi1 - phi2;
  while (dphi<-TMath::Pi()/2) dphi+=TMath::Pi()*2;
  while (dphi>3*TMath::Pi()/2) dphi-=TMath::Pi()*2;

  double deta = eta1 - eta2;

  double px1 = pair->Track1()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  double pz1 = pair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->Track()->P().x();
  double py2 = pair->Track2()->Track()->P().y();
  double pz2 = pair->Track2()->Track()->P().z();

  double pt1 = TMath::Hypot(px1, py1);
  double pt2 = TMath::Hypot(px2, py2);
  double ptmin = pt1>pt2 ? pt2 : pt1;

  double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
    sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));

  fDPhiDEtaNumerator->Fill(dphi, deta);

  if (cosphi > 0) {
    fDPhiDEtaColNumerator->Fill(dphi, deta);
  }
  else {
    fDPhiDEtaColNumerator->Fill(dphi, -eta1-eta2);
  }

  fDPhiNumerator->Fill(dphi);
  fDCosNumerator->Fill(cosphi);

  fDPhiPtNumerator->Fill(dphi, ptmin);
  fDCosPtNumerator->Fill(cosphi, ptmin);

}
//____________________________
void AliFemtoCorrFctnDEtaDPhi::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  double phi1 = pair->Track1()->Track()->P().phi();
  double phi2 = pair->Track2()->Track()->P().phi();
  double eta1 = pair->Track1()->Track()->P().pseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().pseudoRapidity();

  double dphi = phi1 - phi2;
  while (dphi<-TMath::Pi()/2) dphi+=TMath::Pi()*2;
  while (dphi>3*TMath::Pi()/2) dphi-=TMath::Pi()*2;

  double deta = eta1 - eta2;

  double px1 = pair->Track1()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  double pz1 = pair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->Track()->P().x();
  double py2 = pair->Track2()->Track()->P().y();
  double pz2 = pair->Track2()->Track()->P().z();

  double pt1 = TMath::Hypot(px1, py1);
  double pt2 = TMath::Hypot(px2, py2);
  double ptmin = pt1>pt2 ? pt2 : pt1;

  double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
    sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));

  fDPhiDEtaDenominator->Fill(dphi, deta);

  if (cosphi > 0) {
    fDPhiDEtaColDenominator->Fill(dphi, deta);
  }
  else {
    fDPhiDEtaColDenominator->Fill(dphi, -eta1-eta2);
  }

  fDPhiDenominator->Fill(dphi);
  fDCosDenominator->Fill(cosphi);

  fDPhiPtDenominator->Fill(dphi, ptmin);
  fDCosPtDenominator->Fill(cosphi, ptmin);
}


void AliFemtoCorrFctnDEtaDPhi::WriteHistos()
{
  // Write out result histograms
  fDPhiDEtaNumerator->Write();
  fDPhiDEtaDenominator->Write();
  fDPhiDEtaColNumerator->Write();
  fDPhiDEtaColDenominator->Write();
  fDPhiNumerator->Write();
  fDPhiDenominator->Write();
  fDCosNumerator->Write();
  fDCosDenominator->Write();
  fDPhiPtNumerator->Write();
  fDPhiPtDenominator->Write();
  fDCosPtNumerator->Write();
  fDCosPtDenominator->Write();
}

TList* AliFemtoCorrFctnDEtaDPhi::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDEtaNumerator);
  tOutputList->Add(fDPhiDEtaDenominator);
  tOutputList->Add(fDPhiDEtaColNumerator);
  tOutputList->Add(fDPhiDEtaColDenominator);
  tOutputList->Add(fDPhiNumerator);
  tOutputList->Add(fDPhiDenominator);
  tOutputList->Add(fDCosNumerator);
  tOutputList->Add(fDCosDenominator);
  tOutputList->Add(fDPhiPtNumerator);
  tOutputList->Add(fDPhiPtDenominator);
  tOutputList->Add(fDCosPtNumerator);
  tOutputList->Add(fDCosPtDenominator);

  return tOutputList;

}
