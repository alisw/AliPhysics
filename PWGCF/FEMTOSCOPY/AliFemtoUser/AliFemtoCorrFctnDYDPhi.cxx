////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDYDPhi - A correlation function that analyzes              //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and rapidity (y) difference                                                //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch,                                  //
//          Piotr Modzelewski Piotr.Mateusz.Modzelewski@cern.ch               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDYDPhi.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnDYDPhi)
#endif

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnDYDPhi::AliFemtoCorrFctnDYDPhi(const char* title, const int& aPhiBins=20, const int& aYBins=20, const double& mass=0):
  AliFemtoCorrFctn(),
  fDPhiDYNumerator(0),
  fDPhiDYDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
  fDoPtAnalysis(0),
  fDPhiPtNumerator(0),
  fDPhiPtDenominator(0),
  fDCosPtNumerator(0),
  fDCosPtDenominator(0),
  fPhi(0),
  fY(0),
  fPtSumDist(0),
  fYtYtNumerator(0),
  fYtYtDenominator(0),
  fIfCorrectionHist(kNone),
  fPtCorrectionsNum(0),
  fPtCorrectionsDen(0),
  fYCorrectionsNum(0),
  fYCorrectionsDen(0),
  fphiL(0),
  fphiT(0),
  fMass(0)
{

  fphiL = (-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  fphiT = 2*TMath::Pi()+(-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  fMass = mass;

  // set up numerator
  char tTitNumD[101] = "NumDPhiDY";
  strncat(tTitNumD,title, 100);
  fDPhiDYNumerator = new TH2D(tTitNumD,title,aPhiBins,fphiL,fphiT,aYBins,-2.0,2.0);
  // set up denominator
  char tTitDenD[101] = "DenDPhiDY";
  strncat(tTitDenD,title, 100);
  fDPhiDYDenominator = new TH2D(tTitDenD,title,aPhiBins,fphiL,fphiT,aYBins,-2.0,2.0);

  char tTitNum[101] = "PtSumDist";
  strncat(tTitNum,title, 100);
  fPtSumDist = new TH1D(tTitNum,title,200,0,10);
  fPtSumDist->Sumw2();

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

  char tTitPhi[101] = "Phi";
  strncat(tTitPhi,title, 100);
  fPhi = new TH1D(tTitPhi,title,90,-TMath::Pi(),TMath::Pi());

  char tTitY[101] = "Y";
  strncat(tTitY,title, 100);
  fY = new TH1D(tTitY,title,90,-1.2,1.2);

  // set up numerator
  char tTitYtNum[101] = "NumYtYt";
  strncat(tTitYtNum,title, 100);
  fYtYtNumerator = new TH2D(tTitYtNum,title,aPhiBins,1,5,aYBins,1,5);
  // set up denominator
  char tTitYtYtDen[101] = "DenYtYt";
  strncat(tTitYtYtDen,title, 100);
  fYtYtDenominator = new TH2D(tTitYtYtDen,title,aPhiBins,1,5,aYBins,1,5);


  char tTitPtCorrectionsNum[101] = "NumpT1pT2YPhi";
  strncat(tTitPtCorrectionsNum,title, 100);
  char tTitPtCorrectionsDen[101] = "DenpT1pT2YPhi";
  strncat(tTitPtCorrectionsDen,title, 100);

  Int_t nbins[4] = {20,20,aPhiBins,aYBins};
  Double_t xmin[4] = {0,0,-0.5*TMath::Pi(),-2.0};
  Double_t xmax[4] = {4,4,1.5*TMath::Pi(),2.0};


  fPtCorrectionsNum = new THnSparseF(tTitPtCorrectionsNum,title,4,nbins,xmin,xmax);
  fPtCorrectionsDen = new THnSparseF(tTitPtCorrectionsDen,title,4,nbins,xmin,xmax);

  char tTitYCorrectionsNum[101] = "NumY1Y2YPhi";
  strncat(tTitYCorrectionsNum,title, 100);
  char tTitYCorrectionsDen[101] = "DenY1Y2YPhi";
  strncat(tTitYCorrectionsDen,title, 100);

  Double_t xminy[4] = {-1,1,-0.5*TMath::Pi(),-2.0};
  Double_t xmaxy[4] = {-1,1,1.5*TMath::Pi(),2.0};

  fYCorrectionsNum = new THnSparseF(tTitYCorrectionsNum,title,4,nbins,xminy,xmaxy);
  fYCorrectionsDen = new THnSparseF(tTitYCorrectionsDen,title,4,nbins,xminy,xmaxy);

  // THnSparse(const char* name, const char* title, Int_t dim,
  //           const Int_t* nbins, const Double_t* xmin, const Double_t* xmax,
  //           Int_t chunksize);

  // to enable error bar calculation...
  fDPhiDYNumerator->Sumw2();
  fDPhiDYDenominator->Sumw2();
  fDPhiNumerator->Sumw2();
  fDPhiDenominator->Sumw2();
  fDCosNumerator->Sumw2();
  fDCosDenominator->Sumw2();
  fPhi->Sumw2();
  fY->Sumw2();
  fYtYtNumerator->Sumw2();
  fYtYtDenominator->Sumw2();
  fPtCorrectionsNum->Sumw2();
  fPtCorrectionsDen->Sumw2();
}

//____________________________
AliFemtoCorrFctnDYDPhi::AliFemtoCorrFctnDYDPhi(const AliFemtoCorrFctnDYDPhi& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiDYNumerator(0),
  fDPhiDYDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
  fDoPtAnalysis(0),
  fDPhiPtNumerator(0),
  fDPhiPtDenominator(0),
  fDCosPtNumerator(0),
  fDCosPtDenominator(0),
  fPhi(0),
  fY(0),
  fPtSumDist(0),
  fYtYtNumerator(0),
  fYtYtDenominator(0),
  fIfCorrectionHist(kNone),
  fPtCorrectionsNum(0),
  fPtCorrectionsDen(0),
  fYCorrectionsNum(0),
  fYCorrectionsDen(0),
  fphiL(0),
  fphiT(0),
  fMass(0)
{
  // copy constructor
  if (aCorrFctn.fDPhiDYNumerator)
    fDPhiDYNumerator = new TH2D(*aCorrFctn.fDPhiDYNumerator);
  else
    fDPhiDYNumerator = 0;
  if (aCorrFctn.fDPhiDYDenominator)
    fDPhiDYDenominator = new TH2D(*aCorrFctn.fDPhiDYDenominator);
  else
    fDPhiDYDenominator = 0;

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
 if (aCorrFctn.fPhi)
    fPhi = new TH1D(*aCorrFctn.fPhi);
  else
    fPhi = 0;
 if (aCorrFctn.fY)
    fY = new TH1D(*aCorrFctn.fY);
  else
    fY = 0;

 if (aCorrFctn.fPtSumDist)
   fPtSumDist = new TH1D(*aCorrFctn.fPtSumDist);
 else
   fPtSumDist = 0;

 if (aCorrFctn.fYtYtNumerator)
   fYtYtNumerator = new TH2D(*aCorrFctn.fDPhiDYDenominator);
 else
   fYtYtNumerator = 0;

 if (aCorrFctn.fYtYtDenominator)
   fYtYtDenominator = new TH2D(*aCorrFctn.fDPhiDYDenominator);
 else
   fYtYtDenominator = 0;

  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;
  fMass = aCorrFctn.fMass;

//  if (aCorrFctn.fPtCorrectionsNum)
//    fPtCorrectionsNum = new THnSparseF(*aCorrFctn.fPtCorrectionsNum);
//    else
//    fPtCorrectionsNum = 0;

// if (aCorrFctn.fPtCorrectionsDen)
//    fPtCorrectionsDen = new THnSparseF(*aCorrFctn.fPtCorrectionsDen);
//  else
//    fPtCorrectionsDen = 0;



}
//____________________________
AliFemtoCorrFctnDYDPhi::~AliFemtoCorrFctnDYDPhi()
{
  // destructor
  delete fDPhiDYNumerator;
  delete fDPhiDYDenominator;
  delete fPtSumDist;


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
  delete fPhi;
  delete fY;

  delete fYtYtNumerator;
  delete fYtYtDenominator;

  delete fPtCorrectionsNum;
  delete fPtCorrectionsDen;
  delete fYCorrectionsNum;
  delete fYCorrectionsDen;
}
//_________________________
AliFemtoCorrFctnDYDPhi& AliFemtoCorrFctnDYDPhi::operator=(const AliFemtoCorrFctnDYDPhi& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fDPhiDYNumerator)
    fDPhiDYNumerator = new TH2D(*aCorrFctn.fDPhiDYNumerator);
  else
    fDPhiDYNumerator = 0;
  if (aCorrFctn.fDPhiDYDenominator)
    fDPhiDYDenominator = new TH2D(*aCorrFctn.fDPhiDYDenominator);
  else
    fDPhiDYDenominator = 0;

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
  if (aCorrFctn.fPhi)
    fPhi = new TH1D(*aCorrFctn.fPhi);
  else
    fPhi = 0;
  if (aCorrFctn.fY)
    fY = new TH1D(*aCorrFctn.fY);
  else
    fY = 0;

 if (aCorrFctn.fPtSumDist)
   fPtSumDist = new TH1D(*aCorrFctn.fPtSumDist);
 else
   fPtSumDist = 0;

 if (aCorrFctn.fYtYtNumerator)
   fYtYtNumerator = new TH2D(*aCorrFctn.fDPhiDYDenominator);
 else
   fYtYtNumerator = 0;

 if (aCorrFctn.fYtYtDenominator)
   fYtYtDenominator = new TH2D(*aCorrFctn.fDPhiDYDenominator);
 else
   fYtYtDenominator = 0;

 fIfCorrectionHist = kNone;

  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;
  fMass = aCorrFctn.fMass;

// if (aCorrFctn.fPtCorrectionsNum)
//    fPtCorrectionsNum = new THnSparseF(*aCorrFctn.fPtCorrectionsNum);
//  else
//    fPtCorrectionsNum = 0;

// if (aCorrFctn.fPtCorrectionsDen)
//    fPtCorrectionsDen = new THnSparseF(*aCorrFctn.fPtCorrectionsDen);
//  else
//    fPtCorrectionsDen = 0;



  return *this;
}
//_________________________
void AliFemtoCorrFctnDYDPhi::Finish()
{
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDYDPhi::Report()
{
  // create report
  AliFemtoString report = "TPC Ncls Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n", fDPhiDYNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDPhiDYDenominator->GetEntries());
  //  report += mCoulombWeight->Report();

  return report;
}
//____________________________
void AliFemtoCorrFctnDYDPhi::AddRealPair(AliFemtoPair* pair)
{
  // add real (effect) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  double px1 = pair->Track1()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  double pz1 = pair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->Track()->P().x();
  double py2 = pair->Track2()->Track()->P().y();
  double pz2 = pair->Track2()->Track()->P().z();

  double p1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
  double p2 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2);

  double E1 = TMath::Sqrt(p1*p1+fMass*fMass);
  double E2 = TMath::Sqrt(p2*p2+fMass*fMass);

  double y1 = 0.5*TMath::Log((E1+pz1)/(E1-pz1));
  double y2 = 0.5*TMath::Log((E2+pz2)/(E2-pz2));

  double phi1 = pair->Track1()->FourMomentum().Phi();
  double phi2 = pair->Track2()->FourMomentum().Phi();


  double dphi = phi1 - phi2;
  while (dphi<fphiL) dphi+=PIT;
  while (dphi>fphiT) dphi-=PIT;

  double dy = y1 - y2;


   double pt1 = TMath::Hypot(px1, py1);
   double pt2 = TMath::Hypot(px2, py2);
//   double ptmin = pt1>pt2 ? pt2 : pt1;
  fPtSumDist->Fill(pt1+pt2);

//   double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
//     sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));

  fDPhiDYNumerator->Fill(dphi, dy);

  fDPhiNumerator->Fill(dphi);
//   fDCosNumerator->Fill(cosphi);

  if (fDoPtAnalysis) {
//     fDPhiPtNumerator->Fill(dphi, ptmin);
//     fDCosPtNumerator->Fill(cosphi, ptmin);
  }

  fPhi->Fill(phi1);
  fY->Fill(y1);

  double PionMass = 0.13956995;
  double yt1 = TMath::Log(sqrt(1+(pt1/PionMass)*(pt1/PionMass))+(pt1/PionMass));
  double yt2 = TMath::Log(sqrt(1+(pt2/PionMass)*(pt2/PionMass))+(pt2/PionMass));
  fYtYtNumerator->Fill(yt1,yt2);

  if(fIfCorrectionHist)
    {
      if(fIfCorrectionHist == kPt){
	Double_t val[] = {pt1,pt2,dphi,dy};
	fPtCorrectionsNum->Fill(val);
      }
      if(fIfCorrectionHist == kY){
	Double_t val[] = {y1,y2,dphi,dy};
	fYCorrectionsNum->Fill(val);
      }
    }
}
//____________________________
void AliFemtoCorrFctnDYDPhi::AddMixedPair(AliFemtoPair* pair)
{
  // add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  double px1 = pair->Track1()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  double pz1 = pair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->Track()->P().x();
  double py2 = pair->Track2()->Track()->P().y();
  double pz2 = pair->Track2()->Track()->P().z();

  double p1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
  double p2 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2);

  double E1 = TMath::Sqrt(p1*p1+fMass*fMass);
  double E2 = TMath::Sqrt(p2*p2+fMass*fMass);

  double y1 = 0.5*TMath::Log((E1+pz1)/(E1-pz1));
  double y2 = 0.5*TMath::Log((E2+pz2)/(E2-pz2));

  double phi1 = pair->Track1()->FourMomentum().Phi();
  double phi2 = pair->Track2()->FourMomentum().Phi();


  double dphi = phi1 - phi2;
  while (dphi<fphiL) dphi+=PIT;
  while (dphi>fphiT) dphi-=PIT;

  double dy = y1 - y2;

   double pt1 = TMath::Hypot(px1, py1);
   double pt2 = TMath::Hypot(px2, py2);
//   double ptmin = pt1>pt2 ? pt2 : pt1;

//   double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
//     sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));

  fDPhiDYDenominator->Fill(dphi, dy);

  fDPhiDenominator->Fill(dphi);
//   fDCosDenominator->Fill(cosphi);

  //if (fDoPtAnalysis) {
    //   fDPhiPtDenominator->Fill(dphi, ptmin);
    //   fDCosPtDenominator->Fill(cosphi, ptmin);
  //}

  double PionMass = 0.13956995;
    double yt1 = TMath::Log(sqrt(1+(pt1/PionMass)*(pt1/PionMass))+(pt1/PionMass));
    double yt2 = TMath::Log(sqrt(1+(pt2/PionMass)*(pt2/PionMass))+(pt2/PionMass));
    fYtYtDenominator->Fill(yt1,yt2);

  if(fIfCorrectionHist)
    {
      if(fIfCorrectionHist == kPt){
	Double_t val[] = {pt1,pt2,dphi,dy};
	fPtCorrectionsDen->Fill(val);
      }
      if(fIfCorrectionHist == kY){
	Double_t val[] = {y1,y2,dphi,dy};
	fYCorrectionsDen->Fill(val);
      }
    }

}


void AliFemtoCorrFctnDYDPhi::WriteHistos()
{
  // Write out result histograms
  fDPhiDYNumerator->Write();
  fDPhiDYDenominator->Write();
  fPtSumDist->Write();
  /*fDPhiNumerator->Write();
  fDPhiDenominator->Write();
  fDCosNumerator->Write();
  fDCosDenominator->Write();*/
  if (fDoPtAnalysis) {
    fDPhiPtNumerator->Write();
    fDPhiPtDenominator->Write();
    fDCosPtNumerator->Write();
    fDCosPtDenominator->Write();
  }
  fPhi->Write();
  fY->Write();

  if(fIfCorrectionHist){
    if(fIfCorrectionHist==kPt){
      fPtCorrectionsNum->Write();
      fPtCorrectionsDen->Write();
    }
    if(fIfCorrectionHist==kY){
      fYCorrectionsNum->Write();
      fYCorrectionsDen->Write();
    }
  }
}

TList* AliFemtoCorrFctnDYDPhi::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDYNumerator);
  tOutputList->Add(fDPhiDYDenominator);
  tOutputList->Add(fPtSumDist);
  /*tOutputList->Add(fDPhiNumerator);
  tOutputList->Add(fDPhiDenominator);
  tOutputList->Add(fDCosNumerator);
  tOutputList->Add(fDCosDenominator);
  if (fDoPtAnalysis) {
    tOutputList->Add(fDPhiPtNumerator);
    tOutputList->Add(fDPhiPtDenominator);
    tOutputList->Add(fDCosPtNumerator);
    tOutputList->Add(fDCosPtDenominator);
    }*/
  tOutputList->Add(fPhi);
  tOutputList->Add(fY);
  tOutputList->Add(fYtYtNumerator);
  tOutputList->Add(fYtYtDenominator);

  if(fIfCorrectionHist){
    if(fIfCorrectionHist==kPt){
      tOutputList->Add(fPtCorrectionsNum);
      tOutputList->Add(fPtCorrectionsDen);
    }
    if(fIfCorrectionHist==kY){
      tOutputList->Add(fYCorrectionsNum);
      tOutputList->Add(fYCorrectionsDen);
    }
  }

  return tOutputList;
}

void AliFemtoCorrFctnDYDPhi::SetDoPtAnalysis(int do2d)
{
  fDoPtAnalysis = do2d;

  int aPhiBins = fDPhiDYNumerator->GetNbinsX();
  const char *title = fDPhiDYNumerator->GetTitle();

  // set up numerator
  char tTitNumDPhiPt[101] = "NumDPhiPt";
  strncat(tTitNumDPhiPt,title, 100);
  fDPhiPtNumerator = new TH2D(tTitNumDPhiPt,title,aPhiBins*2,-0.5*TMath::Pi(),3./2.*TMath::Pi(), 30, 0.0, 3.0);
  // set up denominator
  char tTitDenDPhiPt[101] = "DenDPhiPt";
  strncat(tTitDenDPhiPt,title, 100);
  fDPhiPtDenominator = new TH2D(tTitDenDPhiPt,title,aPhiBins*2,-0.5*TMath::Pi(),3./2.*TMath::Pi(), 30, 0.0, 3.0);

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

void AliFemtoCorrFctnDYDPhi::SetDo4DCorrectionHist(CorrectionType doCorr)
{
  fIfCorrectionHist = doCorr;
}
