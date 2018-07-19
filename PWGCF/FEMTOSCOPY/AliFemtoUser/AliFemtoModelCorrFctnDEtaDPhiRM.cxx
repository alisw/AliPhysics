////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoModelCorrFctnDEtaDPhiRM - A correlation function that analyzes       //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference, contains histograms for MInv distr                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch, rmaselek@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoModelCorrFctnDEtaDPhiRM.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoModelCorrFctnDEtaDPhiRM)
#endif

/*And some Model libraries..*/
//1
//#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
//#include "AliFemtoModelHiddenInfo.h"
//2
#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
//

//____________________________
AliFemtoModelCorrFctnDEtaDPhiRM::AliFemtoModelCorrFctnDEtaDPhiRM(const char* title, const int& aPhiBins=20, const int& aEtaBins=20, const double m1=0.13956995, const double m2=0.13956995):
  AliFemtoModelCorrFctn(),
  fDPhiDEtaNumeratorTrue(0),
  fDPhiDEtaNumeratorFake(0),
  fDPhiDEtaDenominator(0),
  fDPhiDEtaColNumerator(0),
  fDPhiDEtaColDenominator(0),
  fDPhiNumeratorTrue(0),
  fDPhiNumeratorFake(0),
  fDPhiDenominator(0),
  fDCosNumeratorTrue(0),
  fDCosNumeratorFake(0),
  fDCosDenominator(0),
  fDPhiPtNumerator(0),
  fDPhiPtDenominator(0),
  fDCosPtNumerator(0),
  fDCosPtDenominator(0),
  fPtSumDist(0),
  fInvMassDist(0),
  fM1(m1),
  fM2(m2)
{
  // set up numerator
  char tTitNumDT[101] = "NumDPhiDEtaTrue";
  strncat(tTitNumDT,title, 100);
  fDPhiDEtaNumeratorTrue = new TH2D(tTitNumDT,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);

  char tTitNumDF[101] = "NumDPhiDEtaFake";
  strncat(tTitNumDF,title, 100);
  fDPhiDEtaNumeratorFake = new TH2D(tTitNumDF,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);


  // set up denominator
  char tTitDenD[101] = "DenDPhiDEta";
  strncat(tTitDenD,title, 100);
  fDPhiDEtaDenominator = new TH2D(tTitDenD,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);

  // set up numerator
  char tTitNumR[101] = "NumDPhiDEtaCol";
  strncat(tTitNumR,title, 100);
  fDPhiDEtaColNumerator = new TH2D(tTitNumR,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);
  // set up denominator
  char tTitDenR[101] = "DenDPhiDEtaCol";
  strncat(tTitDenR,title, 100);
  fDPhiDEtaColDenominator = new TH2D(tTitDenR,title,aPhiBins,-0.5*TMath::Pi(),1.5*TMath::Pi(),aEtaBins,-2.0,2.0);

  // set up numerator
  char tTitNumDPhiT[101] = "NumDPhiTrue";
  strncat(tTitNumDPhiT,title, 100);
  fDPhiNumeratorTrue = new TH1D(tTitNumDPhiT,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());
  // set up numerator
  char tTitNumDPhiF[101] = "NumDPhiFake";
  strncat(tTitNumDPhiF,title, 100);
  fDPhiNumeratorFake = new TH1D(tTitNumDPhiF,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());

  // set up denominator
  char tTitDenDPhi[101] = "DenDPhi";
  strncat(tTitDenDPhi,title, 100);
  fDPhiDenominator = new TH1D(tTitDenDPhi,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());

  // set up numerator
  char tTitNumDCosT[101] = "NumDCosTrue";
  strncat(tTitNumDCosT,title, 100);
  fDCosNumeratorTrue = new TH1D(tTitNumDCosT,title,aPhiBins*2,-1.0,1.0);
  // set up numerator
  char tTitNumDCosF[101] = "NumDCosFake";
  strncat(tTitNumDCosF,title, 100);
  fDCosNumeratorFake = new TH1D(tTitNumDCosF,title,aPhiBins*2,-1.0,1.0);

  // set up denominator
  char tTitDenDCos[101] = "DenDCos";
  strncat(tTitDenDCos,title, 100);
  fDCosDenominator = new TH1D(tTitDenDCos,title,aPhiBins*2,-1.0,1.0);

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

  char tTitNum[101] = "PtSumDist";
  strncat(tTitNum,title, 100);
  fPtSumDist = new TH1D(tTitNum,title,200,0,10);
  fPtSumDist->Sumw2();

  char tTitInvMass[101] = "InvariantMassDist";
  strncat(tTitInvMass,title, 100);
  fInvMassDist = new TH1D(tTitInvMass,title,2000,0,4);

  // to enable error bar calculation...
  fDPhiDEtaNumeratorTrue->Sumw2();
  fDPhiDEtaNumeratorFake->Sumw2();
  fDPhiDEtaDenominator->Sumw2();
  fDPhiDEtaColNumerator->Sumw2();
  fDPhiDEtaColDenominator->Sumw2();
  fDPhiNumeratorTrue->Sumw2();
  fDPhiNumeratorFake->Sumw2();
  fDPhiDenominator->Sumw2();
  fDCosNumeratorTrue->Sumw2();
  fDCosNumeratorFake->Sumw2();
  fDCosDenominator->Sumw2();
  fDPhiPtNumerator->Sumw2();
  fDPhiPtDenominator->Sumw2();
  fDCosPtNumerator->Sumw2();
  fDCosPtDenominator->Sumw2();
  fPtSumDist->Sumw2();
  fInvMassDist->Sumw2();

}

//____________________________
AliFemtoModelCorrFctnDEtaDPhiRM::AliFemtoModelCorrFctnDEtaDPhiRM(const AliFemtoModelCorrFctnDEtaDPhiRM& aCorrFctn) :
  AliFemtoModelCorrFctn(),
  fDPhiDEtaNumeratorTrue(0),
  fDPhiDEtaNumeratorFake(0),
  fDPhiDEtaDenominator(0),
  fDPhiDEtaColNumerator(0),
  fDPhiDEtaColDenominator(0),
  fDPhiNumeratorTrue(0),
  fDPhiNumeratorFake(0),
  fDPhiDenominator(0),
  fDCosNumeratorTrue(0),
  fDCosNumeratorFake(0),
  fDCosDenominator(0),
  fDPhiPtNumerator(0),
  fDPhiPtDenominator(0),
  fDCosPtNumerator(0),
  fDCosPtDenominator(0),
  fPtSumDist(0),
  fInvMassDist(0),
  fM1(0),
  fM2(0)
{
  // copy constructor
  if (aCorrFctn.fDPhiDEtaNumeratorTrue)
    fDPhiDEtaNumeratorTrue = new TH2D(*aCorrFctn.fDPhiDEtaNumeratorTrue);
  else
    fDPhiDEtaNumeratorTrue = 0;
  if (aCorrFctn.fDPhiDEtaNumeratorFake)
    fDPhiDEtaNumeratorFake = new TH2D(*aCorrFctn.fDPhiDEtaNumeratorFake);
  else
    fDPhiDEtaNumeratorFake = 0;

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

  if (aCorrFctn.fDPhiNumeratorTrue)
    fDPhiNumeratorTrue = new TH1D(*aCorrFctn.fDPhiNumeratorTrue);
  else
    fDPhiNumeratorTrue = 0;

  if (aCorrFctn.fDPhiNumeratorFake)
    fDPhiNumeratorFake = new TH1D(*aCorrFctn.fDPhiNumeratorFake);
  else
    fDPhiNumeratorFake = 0;

  if (aCorrFctn.fDPhiDenominator)
    fDPhiDenominator = new TH1D(*aCorrFctn.fDPhiDenominator);
  else
    fDPhiDenominator = 0;

  if (aCorrFctn.fDCosNumeratorTrue)
    fDCosNumeratorTrue = new TH1D(*aCorrFctn.fDCosNumeratorTrue);
  else
    fDCosNumeratorTrue = 0;
  if (aCorrFctn.fDCosNumeratorFake)
    fDCosNumeratorFake = new TH1D(*aCorrFctn.fDCosNumeratorFake);
  else
    fDCosNumeratorFake = 0;

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
  if (aCorrFctn.fPtSumDist)
   fPtSumDist = new TH1D(*aCorrFctn.fPtSumDist);
 else
   fPtSumDist = 0;

 if (aCorrFctn.fInvMassDist)
   fInvMassDist = new TH1D(*aCorrFctn.fInvMassDist);
 else
   fInvMassDist = 0;

 fM1 = aCorrFctn.fM1;
 fM2 = aCorrFctn.fM2;
}
//____________________________
AliFemtoModelCorrFctnDEtaDPhiRM::~AliFemtoModelCorrFctnDEtaDPhiRM(){
  // destructor
  if(fDPhiDEtaNumeratorTrue)      delete fDPhiDEtaNumeratorTrue;
  if(fDPhiDEtaNumeratorFake)      delete fDPhiDEtaNumeratorFake;
  if(fDPhiDEtaDenominator)    delete fDPhiDEtaDenominator;
  if(fDPhiDEtaColNumerator)   delete fDPhiDEtaColNumerator;
  if(fDPhiDEtaColDenominator) delete fDPhiDEtaColDenominator;
  if(fDPhiNumeratorTrue)       delete fDPhiNumeratorTrue;
  if(fDPhiNumeratorFake)       delete fDPhiNumeratorFake;
  if(fDPhiDenominator)     delete fDPhiDenominator;
  if(fDCosNumeratorTrue)       delete fDCosNumeratorTrue;
  if(fDCosNumeratorFake)       delete fDCosNumeratorFake;
  if(fDCosDenominator)     delete fDCosDenominator;
  if(fDPhiPtNumerator)     delete fDPhiPtNumerator;
  if(fDPhiPtDenominator)   delete fDPhiPtDenominator;
  if(fDCosPtNumerator)     delete fDCosPtNumerator;
  if(fDCosPtDenominator)   delete fDCosPtDenominator;
  if(fPtSumDist)  delete fPtSumDist;
  if(fInvMassDist) delete fInvMassDist;
}
//_________________________
AliFemtoModelCorrFctnDEtaDPhiRM& AliFemtoModelCorrFctnDEtaDPhiRM::operator=(const AliFemtoModelCorrFctnDEtaDPhiRM& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fDPhiDEtaNumeratorTrue)
    fDPhiDEtaNumeratorTrue = new TH2D(*aCorrFctn.fDPhiDEtaNumeratorTrue);
  else
    fDPhiDEtaNumeratorTrue = 0;

  if (aCorrFctn.fDPhiDEtaNumeratorFake)
    fDPhiDEtaNumeratorFake = new TH2D(*aCorrFctn.fDPhiDEtaNumeratorFake);
  else
    fDPhiDEtaNumeratorFake = 0;

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

  if (aCorrFctn.fDPhiNumeratorTrue)
    fDPhiNumeratorTrue = new TH1D(*aCorrFctn.fDPhiNumeratorTrue);
  else
    fDPhiNumeratorTrue = 0;
  if (aCorrFctn.fDPhiNumeratorFake)
    fDPhiNumeratorFake = new TH1D(*aCorrFctn.fDPhiNumeratorFake);
  else
    fDPhiNumeratorFake = 0;

  if (aCorrFctn.fDPhiDenominator)
    fDPhiDenominator = new TH1D(*aCorrFctn.fDPhiDenominator);
  else
    fDPhiDenominator = 0;

  if (aCorrFctn.fDCosNumeratorTrue)
    fDCosNumeratorTrue = new TH1D(*aCorrFctn.fDCosNumeratorTrue);
  else
    fDCosNumeratorTrue = 0;
  if (aCorrFctn.fDCosNumeratorFake)
    fDCosNumeratorFake = new TH1D(*aCorrFctn.fDCosNumeratorFake);
  else
    fDCosNumeratorFake = 0;

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
  if (aCorrFctn.fPtSumDist)
   fPtSumDist = new TH1D(*aCorrFctn.fPtSumDist);
 else
   fPtSumDist = 0;

 if (aCorrFctn.fInvMassDist)
   fInvMassDist = new TH1D(*aCorrFctn.fInvMassDist);
 else
   fInvMassDist = 0;

 fM1 = aCorrFctn.fM1;
 fM2 = aCorrFctn.fM2;

  return *this;
}
//_________________________
void AliFemtoModelCorrFctnDEtaDPhiRM::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoModelCorrFctnDEtaDPhiRM::Report(){
  // create report
  string stemp = "TPC Ncls Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp,100,"Number of entries in numerator true:\t%E\n",fDPhiDEtaNumeratorTrue->GetEntries());
  snprintf(ctemp,100,"Number of entries in numerator fake:\t%E\n",fDPhiDEtaNumeratorFake->GetEntries());
  stemp += ctemp;
  snprintf(ctemp,100,"Number of entries in denominator:\t%E\n",fDPhiDEtaDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoModelCorrFctnDEtaDPhiRM::AddRealPair( AliFemtoPair* pair)
{
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  // add real (effect) pair
  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

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
  fPtSumDist->Fill(pt1+pt2);

  double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
    sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));

/*weights*/
  Double_t weight = fManager->GetWeight(pair);
  fDPhiDEtaNumeratorTrue->Fill(dphi, deta,weight);

  if (cosphi > 0) {
    fDPhiDEtaColNumerator->Fill(dphi, deta);
  }
  else {
    fDPhiDEtaColNumerator->Fill(dphi, -eta1-eta2);
  }

  fDPhiNumeratorTrue->Fill(dphi,weight);
  fDCosNumeratorTrue->Fill(cosphi,weight);

  fDPhiPtNumerator->Fill(dphi, ptmin);
  fDCosPtNumerator->Fill(cosphi, ptmin);

  // double minv = abs(pair->Track1()->FourMomentum() + pair->Track2()->FourMomentum());

  //Invariant Mass
  double Invpx1=pt1*cos(phi1);
  double Invpy1=pt1*sin(phi1);
  double Invpz1=pt1*TMath::SinH(eta1);

  double Invpx2=pt2*cos(phi2);
  double Invpy2=pt2*sin(phi2);
  double Invpz2=pt2*TMath::SinH(eta2);

  double p21 = Invpx1*Invpx1+Invpy1*Invpy1+Invpz1*Invpz1;
  double p22 = Invpx2*Invpx2+Invpy2*Invpy2+Invpz2*Invpz2;

  //double KaonMass = 0.493677;

  double e1 = TMath::Sqrt(fM1*fM1 + p21);
  double e2 = TMath::Sqrt(fM2*fM2 + p22);

  double minv = TMath::Sqrt(fM1*fM1 + fM2*fM2 + 2*(e1*e2 - Invpx1*Invpx2 - Invpy1*Invpy2 - Invpz1*Invpz2));
  fInvMassDist->Fill(minv);
  // cout<<"Corr: "<<minv<<" masses "<<fM1<<" "<<fM2<<endl;
}
//____________________________
void AliFemtoModelCorrFctnDEtaDPhiRM::AddMixedPair( AliFemtoPair* pair)
{
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  // add mixed (background) pair
  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

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


  Double_t weight = fManager->GetWeight(pair);
  fDPhiDEtaNumeratorFake->Fill(dphi, deta,weight);

  fDPhiDEtaDenominator->Fill(dphi, deta,1.0);

  if (cosphi > 0) {
    fDPhiDEtaColDenominator->Fill(dphi, deta);
  }
  else {
    fDPhiDEtaColDenominator->Fill(dphi, -eta1-eta2);
  }

  fDPhiNumeratorFake->Fill(dphi,weight);
  fDCosNumeratorFake->Fill(cosphi,weight);

  fDPhiDenominator->Fill(dphi,1.0);
  fDCosDenominator->Fill(cosphi,1.0);

  fDPhiPtDenominator->Fill(dphi, ptmin);
  fDCosPtDenominator->Fill(cosphi, ptmin);
}


void AliFemtoModelCorrFctnDEtaDPhiRM::WriteHistos()
{
  // Write out result histograms
  fDPhiDEtaNumeratorTrue->Write();
  fDPhiDEtaNumeratorFake->Write();
  fDPhiDEtaDenominator->Write();
  // fDPhiDEtaColNumerator->Write();
  // fDPhiDEtaColDenominator->Write();
  // fDPhiNumeratorTrue->Write();
  // fDPhiNumeratorFake->Write();
  // fDPhiDenominator->Write();
  // fDCosNumeratorTrue->Write();
  // fDCosNumeratorFake->Write();
  // fDCosDenominator->Write();
  // fDPhiPtNumerator->Write();
  // fDPhiPtDenominator->Write();
  // fDCosPtNumerator->Write();
  // fDCosPtDenominator->Write();
  fPtSumDist->Write();
  fInvMassDist->Write();
}

TList* AliFemtoModelCorrFctnDEtaDPhiRM::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDEtaNumeratorTrue);
  tOutputList->Add(fDPhiDEtaNumeratorFake);
  tOutputList->Add(fDPhiDEtaDenominator);
  // tOutputList->Add(fDPhiDEtaColNumerator);
  // tOutputList->Add(fDPhiDEtaColDenominator);
  // tOutputList->Add(fDPhiNumeratorTrue);
  // tOutputList->Add(fDPhiNumeratorFake);
  // tOutputList->Add(fDPhiDenominator);
  // tOutputList->Add(fDCosNumeratorTrue);
  // tOutputList->Add(fDCosNumeratorFake);
  // tOutputList->Add(fDCosDenominator);
  // tOutputList->Add(fDPhiPtNumerator);
  // tOutputList->Add(fDPhiPtDenominator);
  // tOutputList->Add(fDCosPtNumerator);
  // tOutputList->Add(fDCosPtDenominator);
  tOutputList->Add(fPtSumDist);
  tOutputList->Add(fInvMassDist);

  return tOutputList;

}
