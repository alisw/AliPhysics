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

//____________________________
AliFemtoCorrFctnDEtaDPhi::AliFemtoCorrFctnDEtaDPhi(const char* title, const int& aPhiBins=20, const int& aEtaBins=20):
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
  fDCosPtDenominator(0),
  fPhi(0),
  fEta(0),
  fPtSumDist(0),
  fYtYtNumerator(0),
  fYtYtDenominator(0),
  fIfCorrectionHist(kNone),
  fPtCorrectionsNum(0),
  fPtCorrectionsDen(0),
  fEtaCorrectionsNum(0),
  fEtaCorrectionsDen(0),
  fphiL(0),
  fphiT(0)
{

  fphiL = (-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  fphiT = 2*TMath::Pi()+(-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;

  TString suffix = title;

  auto tTitNumD = "NumDPhiDEta" + suffix;
  fDPhiDEtaNumerator = new TH2D(tTitNumD,title,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  auto tTitDenD = "DenDPhiDEta" + suffix;
  fDPhiDEtaDenominator = new TH2D(tTitDenD,title,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);



  
  // set up numerator
  //char tTitNumDPhi[101] = "NumDPhi";
  //strncat(tTitNumDPhi,title, 100);

  //auto tTitNum = "PtSumDist" + suffix;
  //fPtSumDist = new TH1D(tTitNum,title,200,0,10);
  //fPtSumDist->Sumw2();

  auto tTitNumDPhi = "NumDPhi" + suffix;
  fDPhiNumerator = new TH1D(tTitNumDPhi,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());
  auto tTitDenDPhi = "DenDPhi" + suffix;
  fDPhiDenominator = new TH1D(tTitDenDPhi,title,aPhiBins*2,-0.5*TMath::Pi(),1.5*TMath::Pi());

  auto tTitNumDCos = "NumDCos" + suffix;
  fDCosNumerator = new TH1D(tTitNumDCos,title,aPhiBins*2,-1.0,1.0);
  auto tTitDenDCos = "DenDCos" + suffix;
  fDCosDenominator = new TH1D(tTitDenDCos,title,aPhiBins*2,-1.0,1.0);

  auto tTitPhi = "Phi" + suffix;
  fPhi = new TH1D(tTitPhi,title,90,-TMath::Pi(),TMath::Pi());

  auto tTitEta = "Eta" + suffix;
  fEta = new TH1D(tTitEta,title,90,-1.2,1.2);


  // THnSparse(const char* name, const char* title, Int_t dim,
  //           const Int_t* nbins, const Double_t* xmin, const Double_t* xmax,
  //           Int_t chunksize);

  // to enable error bar calculation...
  fDPhiDEtaNumerator->Sumw2();
  fDPhiDEtaDenominator->Sumw2();
  fDPhiNumerator->Sumw2();
  fDPhiDenominator->Sumw2();
  fDCosNumerator->Sumw2();
  fDCosDenominator->Sumw2();
  fPhi->Sumw2();
  fEta->Sumw2();

}

//____________________________
AliFemtoCorrFctnDEtaDPhi::AliFemtoCorrFctnDEtaDPhi(const AliFemtoCorrFctnDEtaDPhi& aCorrFctn) :
  AliFemtoCorrFctn(aCorrFctn),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
  fDoPtAnalysis(aCorrFctn.fDoPtAnalysis),
  fDPhiPtNumerator(0),
  fDPhiPtDenominator(0),
  fDCosPtNumerator(0),
  fDCosPtDenominator(0),
  fPhi(0),
  fEta(0),
  fPtSumDist(0),
  fYtYtNumerator(0),
  fYtYtDenominator(0),
  fIfCorrectionHist(kNone),
  fPtCorrectionsNum(0),
  fPtCorrectionsDen(0),
  fEtaCorrectionsNum(0),
  fEtaCorrectionsDen(0),
  fphiL(aCorrFctn.fphiL),
  fphiT(aCorrFctn.fphiT)
{
  // copy constructor
  fDPhiDEtaNumerator = new TH2D(*aCorrFctn.fDPhiDEtaNumerator);
  fDPhiDEtaDenominator = new TH2D(*aCorrFctn.fDPhiDEtaDenominator);

  fDPhiNumerator = new TH1D(*aCorrFctn.fDPhiNumerator);
  fDPhiDenominator = new TH1D(*aCorrFctn.fDPhiDenominator);

  fDCosNumerator = new TH1D(*aCorrFctn.fDCosNumerator);
  fDCosDenominator = new TH1D(*aCorrFctn.fDCosDenominator);

  fPhi = new TH1D(*aCorrFctn.fPhi);
  fEta = new TH1D(*aCorrFctn.fEta);
  fPtSumDist = new TH1D(*aCorrFctn.fPtSumDist);

  if (fDoPtAnalysis) {

    fDPhiPtNumerator = new TH2D(*aCorrFctn.fDPhiPtNumerator);
    fDPhiPtDenominator = new TH2D(*aCorrFctn.fDPhiPtDenominator);

    fDCosPtNumerator = new TH2D(*aCorrFctn.fDCosPtNumerator);
    fDCosPtDenominator = new TH2D(*aCorrFctn.fDCosPtDenominator);

    fYtYtNumerator = new TH2D(*aCorrFctn.fYtYtNumerator);
    fYtYtDenominator = new TH2D(*aCorrFctn.fYtYtDenominator);
  }

  if (fIfCorrectionHist) {
    fPtCorrectionsNum = static_cast<THnSparseF*>(aCorrFctn.fPtCorrectionsNum->Clone());
    fPtCorrectionsDen = static_cast<THnSparseF*>(aCorrFctn.fPtCorrectionsDen->Clone());
    fEtaCorrectionsNum = static_cast<THnSparseF*>(aCorrFctn.fEtaCorrectionsNum->Clone());
    fEtaCorrectionsDen = static_cast<THnSparseF*>(aCorrFctn.fEtaCorrectionsDen->Clone());
  }

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
    delete fPtSumDist;
    delete fDPhiPtNumerator;
    delete fDPhiPtDenominator;
    delete fDCosPtNumerator;
    delete fDCosPtDenominator;
    delete fYtYtNumerator;
    delete fYtYtDenominator;

  }
  delete fPhi;
  delete fEta;


  delete fPtCorrectionsNum;
  delete fPtCorrectionsDen;
  delete fEtaCorrectionsNum;
  delete fEtaCorrectionsDen;
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
  if (aCorrFctn.fPhi)
    fPhi = new TH1D(*aCorrFctn.fPhi);
  else
    fPhi = 0;
  if (aCorrFctn.fEta)
    fEta = new TH1D(*aCorrFctn.fEta);
  else
    fEta = 0;

 if (aCorrFctn.fPtSumDist)
   fPtSumDist = new TH1D(*aCorrFctn.fPtSumDist);
 else
   fPtSumDist = 0;

 if (aCorrFctn.fYtYtNumerator)
   fYtYtNumerator = new TH2D(*aCorrFctn.fDPhiDEtaDenominator);
 else
   fYtYtNumerator = 0;

 if (aCorrFctn.fYtYtDenominator)
   fYtYtDenominator = new TH2D(*aCorrFctn.fDPhiDEtaDenominator);
 else
   fYtYtDenominator = 0;

 fIfCorrectionHist = kNone;

  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;

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
  AliFemtoString report = "TPC Ncls Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n",fDPhiDEtaNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n",fDPhiDEtaDenominator->GetEntries());
  //  stemp += mCoulombWeight->Report();
  return report;
}
//____________________________
void AliFemtoCorrFctnDEtaDPhi::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  /*double phi1 = pair->Track1()->Track()->P().Phi();
    double phi2 = pair->Track2()->Track()->P().Phi();
    double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
    double eta2 = pair->Track2()->Track()->P().PseudoRapidity();*/

  double phi1 = pair->Track1()->FourMomentum().Phi();
  double phi2 = pair->Track2()->FourMomentum().Phi();
  double eta1 = pair->Track1()->FourMomentum().PseudoRapidity();
  double eta2 = pair->Track2()->FourMomentum().PseudoRapidity();

  double dphi = phi1 - phi2;
  while (dphi<fphiL) dphi+=PIT;
  while (dphi>fphiT) dphi-=PIT;

  double deta = eta1 - eta2;

  double px1 = pair->Track1()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  //double pz1 = pair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->Track()->P().x();
  double py2 = pair->Track2()->Track()->P().y();
  //double pz2 = pair->Track2()->Track()->P().z();

  double pt1 = TMath::Hypot(px1, py1);
  double pt2 = TMath::Hypot(px2, py2);
  //   double ptmin = pt1>pt2 ? pt2 : pt1;
  

  //   double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
  //   sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));

  fDPhiDEtaNumerator->Fill(dphi, deta);
  fDPhiNumerator->Fill(dphi);
  //   fDCosNumerator->Fill(cosphi);

  if (fDoPtAnalysis) {
    //   fDPhiPtNumerator->Fill(dphi, ptmin);
    //   fDCosPtNumerator->Fill(cosphi, ptmin);

    double PionMass = 0.13956995;
    double yt1 = TMath::Log(sqrt(1+(pt1/PionMass)*(pt1/PionMass))+(pt1/PionMass));
    double yt2 = TMath::Log(sqrt(1+(pt2/PionMass)*(pt2/PionMass))+(pt2/PionMass));
    fYtYtNumerator->Fill(yt1,yt2);
	fPtSumDist->Fill(pt1+pt2);

  }

  fPhi->Fill(phi1);
  fEta->Fill(eta1);


  if(fIfCorrectionHist)
    {
      if(fIfCorrectionHist == kPt){
        Double_t val[] = {pt1,pt2,dphi,deta};
        fPtCorrectionsNum->Fill(val);
      }
      if(fIfCorrectionHist == kEta){
        Double_t val[] = {eta1,eta2,dphi,deta};
        fEtaCorrectionsNum->Fill(val);
      }

    }

}
//____________________________
void AliFemtoCorrFctnDEtaDPhi::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  /*double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();*/

  double phi1 = pair->Track1()->FourMomentum().Phi();
  double phi2 = pair->Track2()->FourMomentum().Phi();
  double eta1 = pair->Track1()->FourMomentum().PseudoRapidity();
  double eta2 = pair->Track2()->FourMomentum().PseudoRapidity();

  double dphi = phi1 - phi2;
  while (dphi<fphiL) dphi+=PIT;
  while (dphi>fphiT) dphi-=PIT;

  double deta = eta1 - eta2;



  fDPhiDEtaDenominator->Fill(dphi, deta);

    double px1 = pair->Track1()->Track()->P().x();
    double py1 = pair->Track1()->Track()->P().y();
    //double pz1 = pair->Track1()->Track()->P().z();

    double px2 = pair->Track2()->Track()->P().x();
    double py2 = pair->Track2()->Track()->P().y();
    //double pz2 = pair->Track2()->Track()->P().z();

    double pt1 = TMath::Hypot(px1, py1);
    double pt2 = TMath::Hypot(px2, py2);
//   double ptmin = pt1>pt2 ? pt2 : pt1;

//   double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
//     sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));


  fDPhiDenominator->Fill(dphi);
//   fDCosDenominator->Fill(cosphi);

  if (fDoPtAnalysis) {




    //   fDPhiPtDenominator->Fill(dphi, ptmin);
    //   fDCosPtDenominator->Fill(cosphi, ptmin);

    double PionMass = 0.13956995;
    double yt1 = TMath::Log(sqrt(1+(pt1/PionMass)*(pt1/PionMass))+(pt1/PionMass));
    double yt2 = TMath::Log(sqrt(1+(pt2/PionMass)*(pt2/PionMass))+(pt2/PionMass));
    fYtYtDenominator->Fill(yt1,yt2);
  }


  if(fIfCorrectionHist)
    {
      if(fIfCorrectionHist == kPt){
        Double_t val[] = {pt1,pt2,dphi,deta};
        fPtCorrectionsDen->Fill(val);
      }
      if(fIfCorrectionHist == kEta){
        Double_t val[] = {eta1,eta2,dphi,deta};
        fEtaCorrectionsDen->Fill(val);
      }
    }

}


void AliFemtoCorrFctnDEtaDPhi::WriteHistos()
{
  // Write out result histograms
  fDPhiDEtaNumerator->Write();
  fDPhiDEtaDenominator->Write();
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
  fEta->Write();

  if(fIfCorrectionHist){
    if(fIfCorrectionHist==kPt){
    fPtCorrectionsNum->Write();
    fPtCorrectionsDen->Write();}
    if(fIfCorrectionHist==kEta){
    fEtaCorrectionsNum->Write();
    fEtaCorrectionsDen->Write();}
  }
}

TList* AliFemtoCorrFctnDEtaDPhi::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDEtaNumerator);
  tOutputList->Add(fDPhiDEtaDenominator);
  tOutputList->Add(fPtSumDist);
  /*tOutputList->Add(fDPhiNumerator);
    tOutputList->Add(fDPhiDenominator);
    tOutputList->Add(fDCosNumerator);
    tOutputList->Add(fDCosDenominator);*/
  if (fDoPtAnalysis) {
    /*tOutputList->Add(fDPhiPtNumerator);
      tOutputList->Add(fDPhiPtDenominator);
      tOutputList->Add(fDCosPtNumerator);
      tOutputList->Add(fDCosPtDenominator);*/
    tOutputList->Add(fYtYtNumerator);
    tOutputList->Add(fYtYtDenominator);
  }
  tOutputList->Add(fPhi);
  tOutputList->Add(fEta);


  if(fIfCorrectionHist){
    if(fIfCorrectionHist==kPt){
      tOutputList->Add(fPtCorrectionsNum);
      tOutputList->Add(fPtCorrectionsDen);
    }
    if(fIfCorrectionHist==kEta){
      tOutputList->Add(fEtaCorrectionsNum);
      tOutputList->Add(fEtaCorrectionsDen);
    }
  }
  return tOutputList;

}

void AliFemtoCorrFctnDEtaDPhi::SetDoPtAnalysis(int do2d)
{
  fDoPtAnalysis = do2d;

  int aPhiBins = fDPhiDEtaNumerator->GetNbinsX();
  int aEtaBins = fDPhiDEtaNumerator->GetNbinsY();
  const char *title = fDPhiDEtaNumerator->GetTitle();

  TString suffix = title;
  
  char tTitNum[101] = "PtSumDist";
  strncat(tTitNum,title, 100);
  fPtSumDist = new TH1D(tTitNum,title,200,0,10);
  fPtSumDist->Sumw2();

  // set up numerator
  //char tTitNumDPhiPt[101] = "NumDPhiPt";
  //strncat(tTitNumDPhiPt,title, 100);
  
  

  auto tTitNumDPhiPt = "NumDPhiPt" + suffix;
  ////
  fDPhiPtNumerator = new TH2D(tTitNumDPhiPt,title,aPhiBins*2,-0.5*TMath::Pi(),3./2.*TMath::Pi(), 30, 0.0, 3.0);
  auto tTitDenDPhiPt = "DenDPhiPt" + suffix;
  fDPhiPtDenominator = new TH2D(tTitDenDPhiPt,title,aPhiBins*2,-0.5*TMath::Pi(),3./2.*TMath::Pi(), 30, 0.0, 3.0);

  auto tTitNumDCosPt = "NumDCosPt" + suffix;
  fDCosPtNumerator = new TH2D(tTitNumDCosPt,title,aPhiBins*2,-1.0,1.0, 30, 0.0, 3.0);
  auto tTitDenDCosPt = "DenDCosPt" + suffix;
  fDCosPtDenominator = new TH2D(tTitDenDCosPt,title,aPhiBins*2,-1.0,1.0, 30, 0.0, 3.0);

  // set up numerator
  auto tTitYtNum = "NumYtYt" + suffix;
  fYtYtNumerator = new TH2D(tTitYtNum,title,aPhiBins,1,5,aEtaBins,1,5);

  // set up denominator
  auto tTitYtYtDen = "DenYtYt" + suffix;
  fYtYtDenominator = new TH2D(tTitYtYtDen,title,aPhiBins,1,5,aEtaBins,1,5);
  fYtYtNumerator->Sumw2();
  fYtYtDenominator->Sumw2();


  fDPhiPtNumerator->Sumw2();
  fDPhiPtDenominator->Sumw2();
  fDCosPtNumerator->Sumw2();
  fDCosPtDenominator->Sumw2();

}

void AliFemtoCorrFctnDEtaDPhi::SetDo4DCorrectionHist(CorrectionType doCorr)
{
  fIfCorrectionHist = doCorr;

  const char *title = fDPhiDEtaNumerator->GetTitle();
  int aPhiBins = fDPhiDEtaNumerator->GetNbinsX();
  int aEtaBins = fDPhiDEtaNumerator->GetNbinsY();

  char tTitPtCorrectionsNum[101] = "NumpT1pT2EtaPhi";
  strncat(tTitPtCorrectionsNum,title, 100);
  char tTitPtCorrectionsDen[101] = "DenpT1pT2EtaPhi";
  strncat(tTitPtCorrectionsDen,title, 100);

  Int_t nbins[4] = {20,20,aPhiBins,aEtaBins};
  Double_t xmin[4] = {0,0,-0.5*TMath::Pi(),-2.0};
  Double_t xmax[4] = {4,4,1.5*TMath::Pi(),2.0};


  fPtCorrectionsNum = new THnSparseF(tTitPtCorrectionsNum,title,4,nbins,xmin,xmax);
  fPtCorrectionsDen = new THnSparseF(tTitPtCorrectionsDen,title,4,nbins,xmin,xmax);

  char tTitEtaCorrectionsNum[101] = "NumEta1Eta2EtaPhi";
  strncat(tTitEtaCorrectionsNum,title, 100);
  char tTitEtaCorrectionsDen[101] = "DenEta1Eta2EtaPhi";
  strncat(tTitEtaCorrectionsDen,title, 100);

  Double_t xmineta[4] = {-1,1,-0.5*TMath::Pi(),-2.0};
  Double_t xmaxeta[4] = {-1,1,1.5*TMath::Pi(),2.0};

  fEtaCorrectionsNum = new THnSparseF(tTitEtaCorrectionsNum,title,4,nbins,xmineta,xmaxeta);
  fEtaCorrectionsDen = new THnSparseF(tTitEtaCorrectionsDen,title,4,nbins,xmineta,xmaxeta);

  fPtCorrectionsNum->Sumw2();
  fPtCorrectionsDen->Sumw2();
}
