////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhiCorrections - A correlation function that analyzes //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDEtaDPhiCorrections.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>
#include "THn.h"

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnDEtaDPhiCorrections)
#endif

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623




//____________________________
AliFemtoCorrFctnDEtaDPhiCorrections::AliFemtoCorrFctnDEtaDPhiCorrections(const char* title, const int& aPhiBins=20, const int& aEtaBins=20):
  AliFemtoCorrFctn(),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
  fDoFullAnalysis(kFALSE),
  fCalculatePairPurity(kFALSE),
  fPhi(0),
  fEta(0),
  fPtSumDist(0),
  fYtYtNumerator(0),
  fYtYtDenominator(0),
  fPairPurity(0),
  fDPhiDEtaNumeratorNoCorr(0),
  fIfCorrection(0),
  fPtCorrectionsNum(0),
  fPtCorrectionsDen(0),
  fEtaCorrectionsNum(0),
  fEtaCorrectionsDen(0),
  fCorrFactorTab(0),
  fpTab(0),
  fPartType(kNoCorrection),
  fphiL(0),
  fphiT(0),
  ifileCorrTab(0),
  fdoPtCorr(0),
  fdoEtaCorr(0),
  fdoPhiCorr(0),
  fdoZVertCorr(0),
  fpartType1(0),
  fpartType2(0),
  fhntReco1(0),
  fhntReco2(0),
  fh1Reco1(0),
  fh1Reco2(0),
  fh2Reco1(0),
  fh2Reco2(0),
  fh3Reco1(0),
  fh3Reco2(0),
  fhCont1(0),
  fhCont2(0),
  fSinglePurity1(0),
  fSinglePurity2(0),
  fCorr1D(kFALSE)
{

  fphiL = (-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  fphiT = 2*TMath::Pi()+(-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;

  // set up numerator
  char tTitNumD[101] = "NumDPhiDEta";
  strncat(tTitNumD,title, 100);
  fDPhiDEtaNumerator = new TH2D(tTitNumD,title,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  char tTitDenD[101] = "DenDPhiDEta";
  strncat(tTitDenD,title, 100);
  fDPhiDEtaDenominator = new TH2D(tTitDenD,title,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);

    // set up numerator
  char tTitNum[101] = "PtSumDist";
  strncat(tTitNum,title, 100);
  fPtSumDist = new TH1D(tTitNum,title,200,0,10);

  fDPhiDEtaNumerator->Sumw2();
  fDPhiDEtaDenominator->Sumw2();
  fPtSumDist->Sumw2();



  // to enable error bar calculation...




}

//____________________________
AliFemtoCorrFctnDEtaDPhiCorrections::AliFemtoCorrFctnDEtaDPhiCorrections(const AliFemtoCorrFctnDEtaDPhiCorrections& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiNumerator(0),
  fDPhiDenominator(0),
  fDCosNumerator(0),
  fDCosDenominator(0),
  fDoFullAnalysis(kFALSE),
  fCalculatePairPurity(kFALSE),
  fPhi(0),
  fEta(0),
  fPtSumDist(0),
  fYtYtNumerator(0),
  fYtYtDenominator(0),
  fPairPurity(0),
  fDPhiDEtaNumeratorNoCorr(0),
  fIfCorrection(0),
  fPtCorrectionsNum(0),
  fPtCorrectionsDen(0),
  fEtaCorrectionsNum(0),
  fEtaCorrectionsDen(0),
  fCorrFactorTab(0),
  fpTab(0),
  fPartType(kNoCorrection),
  fphiL(0),
  fphiT(0),
  ifileCorrTab(0),
  fdoPtCorr(0),
  fdoEtaCorr(0),
  fdoPhiCorr(0),
  fdoZVertCorr(0),
  fpartType1(0),
  fpartType2(0),
  fhntReco1(0),
  fhntReco2(0),
  fh1Reco1(0),
  fh1Reco2(0),
  fh2Reco1(0),
  fh2Reco2(0),
  fh3Reco1(0),
  fh3Reco2(0),
  fhCont1(0),
  fhCont2(0),
  fSinglePurity1(0),
  fSinglePurity2(0),
  fCorr1D(kFALSE)
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


 if (aCorrFctn.fPairPurity)
   fPairPurity = new TH2F(*aCorrFctn.fPairPurity);
 else
    fPairPurity = 0;


 if (aCorrFctn.fDPhiDEtaNumeratorNoCorr)
   fDPhiDEtaNumeratorNoCorr = new TH2F(*aCorrFctn.fDPhiDEtaNumeratorNoCorr);
 else
   fDPhiDEtaNumeratorNoCorr = 0;

  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;

  fPartType = aCorrFctn.fPartType;

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
AliFemtoCorrFctnDEtaDPhiCorrections::~AliFemtoCorrFctnDEtaDPhiCorrections(){
  // destructor
  delete fDPhiDEtaNumerator;
  delete fDPhiDEtaDenominator;
  delete fPtSumDist;

  if(fCalculatePairPurity){
      delete fPairPurity;
      delete fDPhiDEtaNumeratorNoCorr;
  }

  if (fDoFullAnalysis) {
    delete fDPhiNumerator;
    delete fDPhiDenominator;
    delete fDCosNumerator;
    delete fDCosDenominator;

    delete fYtYtNumerator;
    delete fYtYtDenominator;
    delete fPhi;
    delete fEta;
    delete fPtCorrectionsNum;
    delete fPtCorrectionsDen;
    delete fEtaCorrectionsNum;
    delete fEtaCorrectionsDen;
  }

  //corrctions
  if(fhntReco1){
    delete fhntReco1;
    delete fhntReco2;
    delete fhCont1;
    delete fhCont2;
  }
  if(fh1Reco1) delete fh1Reco1;
  if(fh1Reco2) delete fh1Reco2;
  if(fh2Reco1) delete fh2Reco1;
  if(fh2Reco2) delete fh2Reco2;
  if(fh3Reco1) delete fh3Reco1;
  if(fh3Reco2) delete fh3Reco2;

  delete fSinglePurity1;
  delete fSinglePurity2;

}
//_________________________
AliFemtoCorrFctnDEtaDPhiCorrections& AliFemtoCorrFctnDEtaDPhiCorrections::operator=(const AliFemtoCorrFctnDEtaDPhiCorrections& aCorrFctn)
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


 if (aCorrFctn.fPairPurity)
   fPairPurity = new TH2F(*aCorrFctn.fPairPurity);
 else
    fPairPurity = 0;


 if (aCorrFctn.fDPhiDEtaNumeratorNoCorr)
   fDPhiDEtaNumeratorNoCorr = new TH2F(*aCorrFctn.fDPhiDEtaNumeratorNoCorr);
 else
   fDPhiDEtaNumeratorNoCorr = 0;

  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;

  fPartType = aCorrFctn.fPartType;

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
void AliFemtoCorrFctnDEtaDPhiCorrections::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDEtaDPhiCorrections::Report(){
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
void AliFemtoCorrFctnDEtaDPhiCorrections::AddRealPair( AliFemtoPair* pair){
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

  double vert1[3];
  pair->Track1()->Track()->GetPrimaryVertex(vert1);
  double vert2[3];
  pair->Track2()->Track()->GetPrimaryVertex(vert2);

  double corrweight=0; //double corrweightpT1=1;  double corrweightpT2=1;
  //if (fIfCorrection) corrweight = CalculateCorrectionWeight(pt1, pt2);
  if (fIfCorrection)
    {
      corrweight = CalculateCorrectionWeight(pt1, pt2, eta1, eta2, phi1, phi2, vert1[2], vert2[2]);
    }
  else if(fCorr1D)
    {
      corrweight = CalculateCorrectionWeight(pt1, pt2);
      //corrweightpT1 = CalculateCorrectionWeight(pt1);
      //corrweightpT2 = CalculateCorrectionWeight(pt2);
    }

  fPtSumDist->Fill(pt1+pt2,corrweight);
  /*   double ptmin = pt1>pt2 ? pt2 : pt1;

       double cosphi = (px1*px2 + py1*py2 + pz1*pz2)/
       sqrt((px1*px1 + py1*py1 + pz1*pz1)*(px2*px2 + py2*py2 + pz2*pz2));
  */
  if (fIfCorrection || fCorr1D)
    {
      fDPhiDEtaNumerator->Fill(dphi, deta, corrweight);
    }
  else{
    fDPhiDEtaNumerator->Fill(dphi, deta);
  }

  if(fPairPurity){
    double purityweight = GetPurity(pt1,1)*GetPurity(pt2,2);
    //cout<<"Pair purity: "<<purityweight<<" 1: "<<GetPurity(pt1,1)<<" 2: "<<GetPurity(pt2,2)<<endl;
    fPairPurity->Fill(dphi,deta,purityweight);
    fDPhiDEtaNumeratorNoCorr->Fill(dphi,deta,1);
  }

  if (fDoFullAnalysis) {
    //fDPhiPtNumerator->Fill(dphi, ptmin);
    //fDCosPtNumerator->Fill(cosphi, ptmin);

    fDPhiNumerator->Fill(dphi);
    //fDCosNumerator->Fill(cosphi);
    double PionMass = 0.13956995;
    double yt1 = TMath::Log(sqrt(1+(pt1/PionMass)*(pt1/PionMass))+(pt1/PionMass));
    double yt2 = TMath::Log(sqrt(1+(pt2/PionMass)*(pt2/PionMass))+(pt2/PionMass));
    fYtYtNumerator->Fill(yt1,yt2);


    fPhi->Fill(phi1);
    fEta->Fill(eta1);
  }

}
//____________________________
void AliFemtoCorrFctnDEtaDPhiCorrections::AddMixedPair( AliFemtoPair* pair){
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


  double vert1[3];
  pair->Track1()->Track()->GetPrimaryVertex(vert1);
  double vert2[3];
  pair->Track2()->Track()->GetPrimaryVertex(vert2);

  double corrweight=-999;
  //if (fIfCorrection) corrweight = CalculateCorrectionWeight(pt1, pt2);
  if (fIfCorrection)
    {
      corrweight = CalculateCorrectionWeight(pt1, pt2, eta1, eta2, phi1, phi2, vert1[2], vert2[2]);
    }
  else if(fCorr1D)
    {
      corrweight = CalculateCorrectionWeight(pt1, pt2);
    }


  if(fIfCorrection || fCorr1D)
    fDPhiDEtaDenominator->Fill(dphi, deta, corrweight);
  else
    fDPhiDEtaDenominator->Fill(dphi, deta);

  if (fDoFullAnalysis) {
    //fDPhiPtDenominator->Fill(dphi, ptmin);
    //fDCosPtDenominator->Fill(cosphi, ptmin);
    fDPhiDenominator->Fill(dphi);

    double PionMass = 0.13956995;
    double yt1 = TMath::Log(sqrt(1+(pt1/PionMass)*(pt1/PionMass))+(pt1/PionMass));
    double yt2 = TMath::Log(sqrt(1+(pt2/PionMass)*(pt2/PionMass))+(pt2/PionMass));
    fYtYtDenominator->Fill(yt1,yt2);

  }
}


void AliFemtoCorrFctnDEtaDPhiCorrections::WriteHistos()
{
  // Write out result histograms
  fDPhiDEtaNumerator->Write();
  fDPhiDEtaDenominator->Write();
  fPtSumDist->Write();
  if(fCalculatePairPurity){
    fPairPurity->Write();
    fDPhiDEtaNumeratorNoCorr->Write();
  }
  /*fDPhiNumerator->Write();
  fDPhiDenominator->Write();
  fDCosNumerator->Write();
  fDCosDenominator->Write();
  if (fDoFullAnalysis) {
    fDPhiPtNumerator->Write();
    fDPhiPtDenominator->Write();
    fDCosPtNumerator->Write();
    fDCosPtDenominator->Write();
    }*/
  // fPhi->Write();
  // fEta->Write();

}

TList* AliFemtoCorrFctnDEtaDPhiCorrections::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDEtaNumerator);
  tOutputList->Add(fDPhiDEtaDenominator);
  tOutputList->Add(fPtSumDist);

  if(fCalculatePairPurity){
    tOutputList->Add(fPairPurity);
    tOutputList->Add(fDPhiDEtaNumeratorNoCorr);
  }

  if (fDoFullAnalysis) {
    // tOutputList->Add(fDPhiPtNumerator);
    // tOutputList->Add(fDPhiPtDenominator);
    //tOutputList->Add(fDCosPtNumerator);
    //tOutputList->Add(fDCosPtDenominator);
    tOutputList->Add(fDPhiNumerator);
    tOutputList->Add(fDPhiDenominator);
    tOutputList->Add(fDCosNumerator);
    tOutputList->Add(fDCosDenominator);

    tOutputList->Add(fYtYtNumerator);
    tOutputList->Add(fYtYtDenominator);
    tOutputList->Add(fPhi);
    tOutputList->Add(fEta);
  }




  return tOutputList;

}

void AliFemtoCorrFctnDEtaDPhiCorrections::SetCalculatePairPurity(Bool_t dopp)
{
  fCalculatePairPurity = dopp;

  if(fCalculatePairPurity){
    int aPhiBins = fDPhiDEtaNumerator->GetNbinsX();
    int aEtaBins = fDPhiDEtaNumerator->GetNbinsY();
    const char *title = fDPhiDEtaNumerator->GetTitle();

    char tTitNumDPhi[101] = "PairPurity";
    strncat(tTitNumDPhi,title, 100);
    fPairPurity = new TH2F(tTitNumDPhi,title,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
    fPairPurity->Sumw2();

    char tTitNumDPhiNoCorr[101] = "NumDPhiDEtaNoCorr";
    strncat(tTitNumDPhiNoCorr,title, 100);
    fDPhiDEtaNumeratorNoCorr = new TH2F(tTitNumDPhiNoCorr,title,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
    fDPhiDEtaNumeratorNoCorr->Sumw2();
  }
}


void AliFemtoCorrFctnDEtaDPhiCorrections::SetDoFullAnalysis(Bool_t do2d)
{
  fDoFullAnalysis = do2d;

  if(fDoFullAnalysis)
    {

      int aPhiBins = fDPhiDEtaNumerator->GetNbinsX();
      int aEtaBins = fDPhiDEtaNumerator->GetNbinsY();
      const char *title = fDPhiDEtaNumerator->GetTitle();

      // set up numerator
      char tTitNumDPhi[101] = "NumDPhi";
      strncat(tTitNumDPhi,title, 100);
      fDPhiNumerator = new TH1D(tTitNumDPhi,title,aPhiBins*2,fphiL, fphiT);
      // set up denominator
      char tTitDenDPhi[101] = "DenDPhi";
      strncat(tTitDenDPhi,title, 100);
      fDPhiDenominator = new TH1D(tTitDenDPhi,title,aPhiBins*2,fphiL, fphiT);

      // set up numerator
      char tTitNumDCos[101] = "NumDCos";
      strncat(tTitNumDCos,title, 100);
      fDCosNumerator = new TH1D(tTitNumDCos,title,aPhiBins*2,-1.0,1.0);
      // set up denominator
      char tTitDenDCos[101] = "DenDCos";
      strncat(tTitDenDCos,title, 100);
      fDCosDenominator = new TH1D(tTitDenDCos,title,aPhiBins*2,-1.0,1.0);

      // set up numerator
      char tTitYtNum[101] = "NumYtYt";
      strncat(tTitYtNum,title, 100);
      fYtYtNumerator = new TH2D(tTitYtNum,title,aPhiBins,1,5,aEtaBins,1,5);
      // set up denominator
      char tTitYtYtDen[101] = "DenYtYt";
      strncat(tTitYtYtDen,title, 100);
      fYtYtDenominator = new TH2D(tTitYtYtDen,title,aPhiBins,1,5,aEtaBins,1,5);


      char tTitPhi[101] = "Phi";
      strncat(tTitPhi,title, 100);
      fPhi = new TH1D(tTitPhi,title,90,-TMath::Pi(),TMath::Pi());

      char tTitEta[101] = "Eta";
      strncat(tTitEta,title, 100);
      fEta = new TH1D(tTitEta,title,90,-1.2,1.2);

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
      // THnSparse(const char* name, const char* title, Int_t dim,
      //           const Int_t* nbins, const Double_t* xmin, const Double_t* xmax,
      //           Int_t chunksize);

      // to enable error bar calculation...
      fDPhiNumerator->Sumw2();
      fDPhiDenominator->Sumw2();
      fDCosNumerator->Sumw2();
      fDCosDenominator->Sumw2();
      fYtYtNumerator->Sumw2();
      fPhi->Sumw2();
      fEta->Sumw2();
      fYtYtDenominator->Sumw2();
      fPtCorrectionsNum->Sumw2();
      fPtCorrectionsDen->Sumw2();
    }


}



void AliFemtoCorrFctnDEtaDPhiCorrections::LoadCorrectionTabFromROOTFile(const char *file, ParticleType partType1, ParticleType partType2, bool doPtCorr, bool doEtaCorr, bool doPhiCorr, bool doZVertCorr)
{
  fIfCorrection = kTRUE;

  ifileCorrTab = TFile::Open(file);
  fdoPtCorr = doPtCorr;
  fdoEtaCorr = doEtaCorr;
  fdoPhiCorr = doPhiCorr;
  fdoZVertCorr = doZVertCorr;
  fpartType1 = partType1;
  fpartType2 = partType2;

  char* type1 = new char[12];
  char* type2 = new char[12];


  if(fpartType1==kPion) strcpy(type1,"Pion");
  else if(fpartType1==kKaon) strcpy(type1,"Kaon");
  else if (fpartType1==kProton)strcpy(type1,"Proton");
  else if (fpartType1==kAll) strcpy(type1,"All");
  else if(fpartType1==kPionMinus) strcpy(type1,"PionMinus");
  else if(fpartType1==kKaonMinus) strcpy(type1,"KaonMinus");
  else if (fpartType1==kProtonMinus)strcpy(type1,"ProtonMinus");
  else strcpy(type1,"");

  if(fpartType2==kPion) strcpy(type2,"Pion");
  else if(fpartType2==kKaon) strcpy(type2,"Kaon");
  else if (fpartType2==kProton) strcpy(type2,"Proton");
  else if (fpartType2==kAll) strcpy(type2,"All");
  else if(fpartType2==kPionMinus) strcpy(type1,"PionMinus");
  else if(fpartType2==kKaonMinus) strcpy(type1,"KaonMinus");
  else if (fpartType2==kProtonMinus)strcpy(type1,"ProtonMinus");
  else strcpy(type1,"");



  fhntReco1 = (THnT<float>*)(ifileCorrTab->Get(Form("fCorrectionMapData%s",type1)))->Clone();
  fhntReco2 = (THnT<float>*)(ifileCorrTab->Get(Form("fCorrectionMapData%s",type2)))->Clone();
  fhCont1 = (TH1D*)(ifileCorrTab->Get(Form("SecondariesContamination%s",type1)))->Clone();
  fhCont2 = (TH1D*)(ifileCorrTab->Get(Form("SecondariesContamination%s",type2)))->Clone();

  cout<<"!!!!!!!!!!!!!"<<endl;
  if(fCalculatePairPurity){
    fSinglePurity1 = (TH1F*)(ifileCorrTab->Get(Form("hPurity%s",type1)))->Clone();
    fSinglePurity2 = (TH1F*)(ifileCorrTab->Get(Form("hPurity%s",type2)))->Clone();

  cout<<"!!! load first single purity:" <<  fSinglePurity1<<endl;
  cout<< "!!! load seconf single purity: "<< fSinglePurity2<<endl;

  }

  delete[] type1;
  delete[] type2;

  double fhntReco1_nbins = fhntReco1->GetNbins();
  double fhntReco2_nbins = fhntReco2->GetNbins();

  int boolSum = fdoPtCorr+fdoEtaCorr+fdoPhiCorr+fdoZVertCorr;
  /*if(boolSum == 0)
    {
      return 1;
      }*/
  if(boolSum == 1)
    {

      if(fdoPtCorr == 1)
	{
	  fh1Reco1 = (TH1F*)(fhntReco1->Projection(0))->Clone();
	  fh1Reco2 = (TH1F*)(fhntReco2->Projection(0))->Clone();
	  fh1Reco1->Scale(1./fhntReco1_nbins*fh1Reco1->GetNbinsX());
	  fh1Reco2->Scale(1./fhntReco2_nbins*fh1Reco2->GetNbinsX());

	}

      else if(fdoEtaCorr == 1)
	{
	  fh1Reco1 = (TH1F*)(fhntReco1->Projection(1))->Clone();
	  fh1Reco2 = (TH1F*)(fhntReco2->Projection(1))->Clone();
	  fh1Reco1->Scale(1./fhntReco1_nbins*fh1Reco1->GetNbinsX());
	  fh1Reco2->Scale(1./fhntReco2_nbins*fh1Reco2->GetNbinsX());
	}

      else if(fdoPhiCorr == 1)
	{
	  fh1Reco1 = (TH1F*)(fhntReco1->Projection(2))->Clone();
	  fh1Reco2 = (TH1F*)(fhntReco2->Projection(2))->Clone();
	  fh1Reco1->Scale(1./fhntReco1_nbins*fh1Reco1->GetNbinsX());
	  fh1Reco2->Scale(1./fhntReco2_nbins*fh1Reco2->GetNbinsX());
	}

      else if(fdoZVertCorr == 1)
	{
	  fh1Reco1 = (TH1F*)(fhntReco1->Projection(3))->Clone();
	  fh1Reco2 = (TH1F*)(fhntReco2->Projection(3))->Clone();
	  fh1Reco1->Scale(1./fhntReco1_nbins*fh1Reco1->GetNbinsX());
	  fh1Reco2->Scale(1./fhntReco2_nbins*fh1Reco2->GetNbinsX());
	}

    }

  else if(boolSum == 2)
    {
      if(fdoPtCorr == 1 && fdoEtaCorr == 1)
	{
	  fh2Reco1 = (TH2F*)(fhntReco1->Projection(1,0))->Clone();
	  fh2Reco2 = (TH2F*)(fhntReco2->Projection(1,0))->Clone();
	  fh2Reco1->Scale(1./fhntReco1_nbins*fh2Reco1->GetNbinsX()*fh2Reco1->GetNbinsY());
	  fh2Reco2->Scale(1./fhntReco2_nbins*fh2Reco2->GetNbinsX()*fh2Reco2->GetNbinsY());
	}

      else if(fdoPtCorr == 1 && fdoPhiCorr == 1)
	{
	  fh2Reco1 = (TH2F*)(fhntReco1->Projection(2,0))->Clone();
	  fh2Reco2 = (TH2F*)(fhntReco2->Projection(2,0))->Clone();
	  fh2Reco1->Scale(1./fhntReco1_nbins*fh2Reco1->GetNbinsX()*fh2Reco1->GetNbinsY());
	  fh2Reco2->Scale(1./fhntReco2_nbins*fh2Reco2->GetNbinsX()*fh2Reco2->GetNbinsY());
	}

      else if(fdoPtCorr == 1 && fdoZVertCorr == 1)
	{
	  fh2Reco1 = (TH2F*)(fhntReco1->Projection(3,0))->Clone();
	  fh2Reco2 = (TH2F*)(fhntReco2->Projection(3,0))->Clone();
	  fh2Reco1->Scale(1./fhntReco1_nbins*fh2Reco1->GetNbinsX()*fh2Reco1->GetNbinsY());
	  fh2Reco2->Scale(1./fhntReco2_nbins*fh2Reco2->GetNbinsX()*fh2Reco2->GetNbinsY());
	}
      else if(fdoEtaCorr == 1 && fdoPhiCorr == 1)
	{
	  fh2Reco1 = (TH2F*)(fhntReco1->Projection(2,1))->Clone();
	  fh2Reco2 = (TH2F*)(fhntReco2->Projection(2,1))->Clone();
	  fh2Reco1->Scale(1./fhntReco1_nbins*fh2Reco1->GetNbinsX()*fh2Reco1->GetNbinsY());
	  fh2Reco2->Scale(1./fhntReco2_nbins*fh2Reco2->GetNbinsX()*fh2Reco2->GetNbinsY());
	}
      else if(fdoEtaCorr == 1 && fdoZVertCorr == 1)
	{
	  fh2Reco1 = (TH2F*)(fhntReco1->Projection(3,1))->Clone();
	  fh2Reco2 = (TH2F*)(fhntReco2->Projection(3,1))->Clone();
	  fh2Reco1->Scale(1./fhntReco1_nbins*fh2Reco1->GetNbinsX()*fh2Reco1->GetNbinsY());
	  fh2Reco2->Scale(1./fhntReco2_nbins*fh2Reco2->GetNbinsX()*fh2Reco2->GetNbinsY());
	}
      else if(fdoPhiCorr == 1 && fdoZVertCorr == 1)
	{
	  fh2Reco1 = (TH2F*)(fhntReco1->Projection(3,2))->Clone();
	  fh2Reco2 = (TH2F*)(fhntReco2->Projection(3,2))->Clone();
	  fh2Reco1->Scale(1./fhntReco1_nbins*fh2Reco1->GetNbinsX()*fh2Reco1->GetNbinsY());
	  fh2Reco2->Scale(1./fhntReco2_nbins*fh2Reco2->GetNbinsX()*fh2Reco2->GetNbinsY());
	}
    }


  else if(boolSum == 3)
    {
      if(fdoPtCorr == 1 && fdoEtaCorr == 1 && fdoPhiCorr == 1)
	{
	  fh3Reco1 = (TH3F*)(fhntReco1->Projection(0,1,2))->Clone();
	  fh3Reco2 = (TH3F*)(fhntReco2->Projection(0,1,2))->Clone();
	  fh3Reco1->Scale(1./fhntReco1_nbins*fh3Reco1->GetNbinsX()*fh3Reco1->GetNbinsY()*fh3Reco1->GetNbinsZ());
	  fh3Reco2->Scale(1./fhntReco2_nbins*fh3Reco2->GetNbinsX()*fh3Reco2->GetNbinsY()*fh3Reco2->GetNbinsZ());

	}

      else if(fdoPtCorr == 1 && fdoEtaCorr == 1 && fdoZVertCorr == 1)
	{
	  fh3Reco1 = (TH3F*)(fhntReco1->Projection(0,1,3))->Clone();
	  fh3Reco2 = (TH3F*)(fhntReco2->Projection(0,1,3))->Clone();
	  fh3Reco1->Scale(1./fhntReco1_nbins*fh3Reco1->GetNbinsX()*fh3Reco1->GetNbinsY()*fh3Reco1->GetNbinsZ());
	  fh3Reco2->Scale(1./fhntReco2_nbins*fh3Reco2->GetNbinsX()*fh3Reco2->GetNbinsY()*fh3Reco2->GetNbinsZ());
	}

      else if(fdoPtCorr == 1 && fdoPhiCorr == 1 && fdoZVertCorr == 1)
	{
	  fh3Reco1 = (TH3F*)(fhntReco1->Projection(0,2,3))->Clone();
	  fh3Reco2 = (TH3F*)(fhntReco2->Projection(0,2,3))->Clone();
	  fh3Reco1->Scale(1./fhntReco1_nbins*fh3Reco1->GetNbinsX()*fh3Reco1->GetNbinsY()*fh3Reco1->GetNbinsZ());
	  fh3Reco2->Scale(1./fhntReco2_nbins*fh3Reco2->GetNbinsX()*fh3Reco2->GetNbinsY()*fh3Reco2->GetNbinsZ());
	}

      else if(fdoEtaCorr == 1 && fdoPhiCorr == 1 && fdoZVertCorr == 1)
	{
	  fh3Reco1 = (TH3F*)(fhntReco1->Projection(1,2,3))->Clone();
	  fh3Reco2 = (TH3F*)(fhntReco2->Projection(1,2,3))->Clone();
	  fh3Reco1->Scale(1./fhntReco1_nbins*fh3Reco1->GetNbinsX()*fh3Reco1->GetNbinsY()*fh3Reco1->GetNbinsZ());
	  fh3Reco2->Scale(1./fhntReco2_nbins*fh3Reco2->GetNbinsX()*fh3Reco2->GetNbinsY()*fh3Reco2->GetNbinsZ());
	}
    }

  /*else if(boolSum == 4)
    {
    }*/

  ifileCorrTab->Close();

}

void AliFemtoCorrFctnDEtaDPhiCorrections::LoadCorrectionTabFromROOTFile1D(const char *file, ParticleType partType1, ParticleType partType2)
{
  fCorr1D = kTRUE;

  ifileCorrTab = TFile::Open(file);

  fpartType1 = partType1;
  fpartType2 = partType2;


  char type1[12];
  char type2[12];


  if(fpartType1==kPion) strcpy(type1,"Pion");
  else if(fpartType1==kKaon) strcpy(type1,"Kaon");
  else if (fpartType1==kProton)strcpy(type1,"Proton");
  else if (fpartType1==kAll) strcpy(type1,"All");
  else if(fpartType1==kPionMinus) strcpy(type1,"PionMinus");
  else if (fpartType1==kKaonMinus) strcpy(type1,"KaonMinus");
  else if (fpartType1==kProtonMinus)strcpy(type1,"ProtonMinus");
  else strcpy(type1,"");

  if(fpartType2==kPion) strcpy(type2,"Pion");
  else if(fpartType2==kKaon) strcpy(type2,"Kaon");
  else if (fpartType2==kProton) strcpy(type2,"Proton");
  else if (fpartType2==kAll) strcpy(type2,"All");
  else if(fpartType2==kPionMinus) strcpy(type1,"PionMinus");
  else if(fpartType2==kKaonMinus) strcpy(type1,"KaonMinus");
  else if (fpartType2==kProtonMinus)strcpy(type1,"ProtonMinus");
  else strcpy(type1,"");

  fhCont1 = (TH1D*)(ifileCorrTab->Get(Form("CorrectionFactorPtEffandCont%s",type1)));//->Clone();
  fhCont2 = (TH1D*)(ifileCorrTab->Get(Form("CorrectionFactorPtEffandCont%s",type2)));//->Clone();


  if(fCalculatePairPurity){
    fSinglePurity1 = (TH1F*)(ifileCorrTab->Get(Form("hPurity%s",type1)));
    fSinglePurity2 = (TH1F*)(ifileCorrTab->Get(Form("hPurity%s",type2)));
  }

  ifileCorrTab->Close();

}

void AliFemtoCorrFctnDEtaDPhiCorrections::SetCorrectionTab(ParticleType partType)
{

  double pttab[] = {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2, 2.225, 2.25, 2.275, 2.3, 2.325, 2.35, 2.375, 2.4, 2.425, 2.45, 2.475, 2.5, 2.525, 2.55, 2.575, 2.6, 2.625, 2.65, 2.675, 2.7, 2.725, 2.75, 2.775, 2.8, 2.825, 2.85, 2.875, 2.9, 2.925, 2.95, 2.975, 3, 3.025, 3.05, 3.075, 3.1, 3.125, 3.15, 3.175, 3.2, 3.225, 3.25, 3.275, 3.3, 3.325, 3.35, 3.375, 3.4, 3.425, 3.45, 3.475, 3.5, 3.525, 3.55, 3.575, 3.6, 3.625, 3.65, 3.675, 3.7, 3.725, 3.75, 3.775, 3.8, 3.825, 3.85, 3.875, 3.9, 3.925, 3.95, 3.975, 4, 4.025, 4.05, 4.075, 4.1, 4.125, 4.15, 4.175, 4.2, 4.225, 4.25, 4.275, 4.3, 4.325, 4.35, 4.375, 4.4, 4.425, 4.45, 4.475, 4.5, 4.525, 4.55, 4.575, 4.6, 4.625, 4.65, 4.675, 4.7, 4.725, 4.75};

  double pioncorrtab[] = {0, 0, 0, 0, 0, 0, 0, 1.40089, 1.40089, 1.29482, 1.29482, 1.25595, 1.22529, 1.22529, 1.23099, 1.32027, 1.32027, 1.44774, 1.44774, 1.74645, 1.8619, 1.8619, 1.82089, 1.78506, 1.78506, 1.75918, 1.75918, 1.74951, 1.74614, 1.74614, 1.74006, 1.73229, 1.73229, 1.72844, 1.72844, 1.72306, 1.71906, 1.71906, 1.71375, 1.71301, 1.71301, 1.70381, 1.70381, 1.69975, 1.69242, 1.69242, 1.69013, 1.67698, 1.67698, 1.6772, 1.6772, 1.67118, 1.66607, 1.66607, 1.66131, 1.67228, 1.67228, 1.66834, 1.66834, 1.66031, 1.6588, 1.6588, 1.6555, 1.64923, 1.64923, 1.6467, 1.6467, 1.63894, 1.63682, 1.63682, 1.6297, 1.62904, 1.62904, 1.63007, 1.63007, 1.62832, 1.62557, 1.62557, 1.62687, 1.62928, 1.62928, 1.62767, 1.62767, 1.62767, 1.62767, 1.62767, 1.62767, 1.63415, 1.63415, 1.63415, 1.63415, 1.63415, 1.63415, 1.63415, 1.64141, 1.64141, 1.64141, 1.64141, 1.64141, 1.64141, 1.65191, 1.65191, 1.65191, 1.65191, 1.65191, 1.65191, 1.65191, 1.66838, 1.66838, 1.66838, 1.66838, 1.66838, 1.66838, 1.6839, 1.6839, 1.6839, 1.6839, 1.6839, 1.6839, 1.69601, 1.69601, 1.69601, 1.69601, 1.69601, 1.69601, 1.69601, 1.70062, 1.70062, 1.70062, 1.70062, 1.70062, 1.70062, 1.68668, 1.68668, 1.68668, 1.68668, 1.68668, 1.68668, 1.68668, 1.68182, 1.68182, 1.68182, 1.68182, 1.68182, 1.68182, 1.681, 1.681, 1.681, 1.681, 1.681, 1.681, 1.67749, 1.67749, 1.67749, 1.67749, 1.67749, 1.67749, 1.67749, 1.66558, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.67223, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.66872, 1.64419};

  double protoncorrtab[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 357.585, 357.585, 8.66944, 2.10995, 2.10995, 1.50443, 1.73168, 1.73168, 2.3605, 2.3605, 4.7726, 4.40359, 4.40359, 3.0307, 2.49649, 2.49649, 2.2231, 2.2231, 2.11247, 2.05862, 2.05862, 2.00703, 1.9623, 1.9623, 1.93393, 1.93393, 1.9101, 1.89334, 1.89334, 1.87734, 1.86342, 1.86342, 1.85075, 1.85075, 1.83985, 1.83684, 1.83684, 1.82915, 1.81832, 1.81832, 1.81215, 1.81215, 1.7998, 1.79524, 1.79524, 1.78568, 1.79989, 1.79989, 1.7973, 1.7973, 1.79591, 1.78468, 1.78468, 1.78037, 1.77394, 1.77394, 1.77198, 1.77198, 1.76736, 1.76875, 1.76875, 1.76221, 1.75729, 1.75729, 1.75397, 1.75397, 1.75229, 1.74918, 1.74918, 1.75064, 1.75643, 1.75643, 1.75765, 1.75765, 1.75765, 1.75765, 1.75765, 1.75765, 1.76345, 1.76345, 1.76345, 1.76345, 1.76345, 1.76345, 1.76345, 1.76901, 1.76901, 1.76901, 1.76901, 1.76901, 1.76901, 1.78291, 1.78291, 1.78291, 1.78291, 1.78291, 1.78291, 1.78291, 1.80009, 1.80009, 1.80009, 1.80009, 1.80009, 1.80009, 1.81064, 1.81064, 1.81064, 1.81064, 1.81064, 1.81064, 1.81765, 1.81765, 1.81765, 1.81765, 1.81765, 1.81765, 1.81765, 1.79549, 1.79549, 1.79549, 1.79549, 1.79549, 1.79549, 1.80455, 1.80455, 1.80455, 1.80455, 1.80455, 1.80455, 1.80455, 1.78912, 1.78912, 1.78912, 1.78912, 1.78912, 1.78912, 1.78501, 1.78501, 1.78501, 1.78501, 1.78501, 1.78501, 1.79512, 1.79512, 1.79512, 1.79512, 1.79512, 1.79512, 1.79512, 1.77138, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.784, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152, 1.75152};

  double kaoncorrtab[] = {0, 0, 0, 0, 0, 0, 0, 8.43268, 8.43268, 3.30657, 3.30657, 2.5102, 2.16256, 2.16256, 2.03757, 2.27166, 2.27166, 2.70432, 2.70432, 4.06234, 4.69199, 4.69199, 4.13074, 3.75139, 3.75139, 3.48381, 3.48381, 3.29762, 3.15261, 3.15261, 3.03022, 2.91874, 2.91874, 2.82421, 2.82421, 2.7388, 2.65961, 2.65961, 2.58426, 2.5174, 2.5174, 2.45378, 2.45378, 2.39687, 2.34699, 2.34699, 2.30247, 2.25299, 2.25299, 2.22443, 2.22443, 2.18303, 2.16012, 2.16012, 2.13083, 2.12806, 2.12806, 2.11376, 2.11376, 2.09566, 2.07526, 2.07526, 2.05378, 2.03252, 2.03252, 2.02466, 2.02466, 2.00531, 1.98945, 1.98945, 1.97877, 1.97226, 1.97226, 1.95475, 1.95475, 1.94838, 1.9314, 1.9314, 1.92571, 1.96346, 1.96346, 1.92849, 1.92849, 1.92849, 1.92849, 1.92849, 1.92849, 1.90949, 1.90949, 1.90949, 1.90949, 1.90949, 1.90949, 1.90949, 1.88743, 1.88743, 1.88743, 1.88743, 1.88743, 1.88743, 1.87486, 1.87486, 1.87486, 1.87486, 1.87486, 1.87486, 1.87486, 1.87785, 1.87785, 1.87785, 1.87785, 1.87785, 1.87785, 1.8757, 1.8757, 1.8757, 1.8757, 1.8757, 1.8757, 1.87948, 1.87948, 1.87948, 1.87948, 1.87948, 1.87948, 1.87948, 1.86148, 1.86148, 1.86148, 1.86148, 1.86148, 1.86148, 1.84329, 1.84329, 1.84329, 1.84329, 1.84329, 1.84329, 1.84329, 1.83105, 1.83105, 1.83105, 1.83105, 1.83105, 1.83105, 1.81955, 1.81955, 1.81955, 1.81955, 1.81955, 1.81955, 1.79944, 1.79944, 1.79944, 1.79944, 1.79944, 1.79944, 1.79944, 1.79345, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.80077, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.78333, 1.74958};

  double allcorrtab[] = {0, 0, 0, 0, 0, 0, 0, 1.46883, 1.46883, 1.3528, 1.3528, 1.30939, 1.26936, 1.26936, 1.23645, 1.21359, 1.21359, 1.19759, 1.19759, 1.18565, 1.17772, 1.17772, 1.17203, 1.16739, 1.16739, 1.16398, 1.16398, 1.16201, 1.16065, 1.16065, 1.16012, 1.16009, 1.16009, 1.16044, 1.16044, 1.16104, 1.16139, 1.16139, 1.16134, 1.16278, 1.16278, 1.1631, 1.1631, 1.16227, 1.16152, 1.16152, 1.16066, 1.15984, 1.15984, 1.15932, 1.15932, 1.15912, 1.15818, 1.15818, 1.15877, 1.16754, 1.16754, 1.17075, 1.17075, 1.17047, 1.16995, 1.16995, 1.16885, 1.16845, 1.16845, 1.16824, 1.16824, 1.16771, 1.16704, 1.16704, 1.16681, 1.16723, 1.16723, 1.16819, 1.16819, 1.16811, 1.16974, 1.16974, 1.17217, 1.16759, 1.16759, 1.17376, 1.17376, 1.17376, 1.17376, 1.17376, 1.17376, 1.18247, 1.18247, 1.18247, 1.18247, 1.18247, 1.18247, 1.18247, 1.18916, 1.18916, 1.18916, 1.18916, 1.18916, 1.18916, 1.19649, 1.19649, 1.19649, 1.19649, 1.19649, 1.19649, 1.19649, 1.20315, 1.20315, 1.20315, 1.20315, 1.20315, 1.20315, 1.20984, 1.20984, 1.20984, 1.20984, 1.20984, 1.20984, 1.21236, 1.21236, 1.21236, 1.21236, 1.21236, 1.21236, 1.21236, 1.21272, 1.21272, 1.21272, 1.21272, 1.21272, 1.21272, 1.21416, 1.21416, 1.21416, 1.21416, 1.21416, 1.21416, 1.21416, 1.21308, 1.21308, 1.21308, 1.21308, 1.21308, 1.21308, 1.21332, 1.21332, 1.21332, 1.21332, 1.21332, 1.21332, 1.21204, 1.21204, 1.21204, 1.21204, 1.21204, 1.21204, 1.21204, 1.21006, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.21141, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2092, 1.2066};

  fpTab = new double[190];
  for(int i=0;i<190;i++)
    fpTab[i]=pttab[i];

  if(partType==kPion || partType==kPionMinus)
    {
      fCorrFactorTab = new double[190];
      for(int i=0;i<190;i++)
	fCorrFactorTab[i] = pioncorrtab[i];
    }
  else if(partType==kKaon || partType==kKaonMinus)
    {
      fCorrFactorTab = new double[190];
      for(int i=0;i<190;i++)
	fCorrFactorTab[i] = kaoncorrtab[i];
    }
  else if(partType==kProton||partType==kProtonMinus)
    {
      fCorrFactorTab = new double[190];
      for(int i=0;i<190;i++)
	fCorrFactorTab[i] = protoncorrtab[i];
    }
  else if(partType==kAll)
    {
      fCorrFactorTab = new double[190];
      for(int i=0;i<190;i++)
	fCorrFactorTab[i] = allcorrtab[i];
    }
}

double AliFemtoCorrFctnDEtaDPhiCorrections::CalculateCorrectionWeight(double pT1, double pT2)
{
   double w1=0., w2=0.;
   if(pT1 > fhCont1->GetXaxis()->GetXmin() && pT1 < fhCont1->GetXaxis()->GetXmax() && pT2 > fhCont2->GetXaxis()->GetXmin() && pT2 < fhCont2->GetXaxis()->GetXmax())
     {
       w1 = fhCont1->GetBinContent(fhCont1->FindFixBin(pT1));
       w2 = fhCont2->GetBinContent(fhCont2->FindFixBin(pT2));

       return w1*w2;
     }
   else
     return 0;
}


double AliFemtoCorrFctnDEtaDPhiCorrections::CalculateCorrectionWeight(double pT1)
{
   double w1=0.;
   if(pT1 > fhCont1->GetXaxis()->GetXmin() && pT1 < fhCont1->GetXaxis()->GetXmax())
     {
       w1 = fhCont1->GetBinContent(fhCont1->FindFixBin(pT1));
       return w1;
     }
   else
     return 0;
}

double AliFemtoCorrFctnDEtaDPhiCorrections::CalculateCorrectionWeight(double pT1, double pT2, double eta1, double eta2, double phi1, double phi2, double zvert1, double zvert2)
{

    double w1=0., w2=0.;
    double eps1=0., eps2=0;
    double cont1=0., cont2=0; //w=(1-cont)/eps
    phi1 += TMath::Pi();
    phi2 += TMath::Pi();

    if(pT1 > fhCont1->GetXaxis()->GetXmin() && pT1 < fhCont1->GetXaxis()->GetXmax() && pT2 > fhCont2->GetXaxis()->GetXmin() && pT2 < fhCont2->GetXaxis()->GetXmax())
      {
	cont1 = fhCont1->GetBinContent(fhCont1->FindFixBin(pT1));
	cont2 = fhCont1->GetBinContent(fhCont2->FindFixBin(pT2));
      }
    else
      return 0;

    int boolSum = fdoPtCorr+fdoEtaCorr+fdoPhiCorr+fdoZVertCorr;
    if(boolSum == 0)
      {
	return 1;
      }
    else if(boolSum == 1)
      {

	if(fdoPtCorr == 1)
	  {
	    if(pT1 > fh1Reco1->GetXaxis()->GetXmin() && pT1 < fh1Reco1->GetXaxis()->GetXmax() && pT2 > fh1Reco2->GetXaxis()->GetXmin() && pT2 < fh1Reco2->GetXaxis()->GetXmax())
	      {
		eps1 = fh1Reco1->GetBinContent(fh1Reco1->FindFixBin(pT1));
		eps2 = fh1Reco2->GetBinContent(fh1Reco2->FindFixBin(pT2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;
		return w1*w2;
	      }
	    else
	      return 0;
	  }

	else if(fdoEtaCorr == 1)
	  {
	    if(eta1 > fh1Reco1->GetXaxis()->GetXmin() && eta1 < fh1Reco1->GetXaxis()->GetXmax() && eta2 > fh1Reco2->GetXaxis()->GetXmin() && eta2 < fh1Reco2->GetXaxis()->GetXmax())
	      {
		eps1 = fh1Reco1->GetBinContent(fh1Reco1->FindFixBin(eta1));
		eps2 = fh1Reco2->GetBinContent(fh1Reco2->FindFixBin(eta2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;
	  }

	else if(fdoPhiCorr == 1)
	  {
	    if(phi1 > fh1Reco1->GetXaxis()->GetXmin() && phi1 < fh1Reco1->GetXaxis()->GetXmax() && phi2 > fh1Reco2->GetXaxis()->GetXmin() && phi2 < fh1Reco2->GetXaxis()->GetXmax())
	      {
		eps1 = fh1Reco1->GetBinContent(fh1Reco1->FindFixBin(phi1));
		eps2 = fh1Reco2->GetBinContent(fh1Reco2->FindFixBin(phi2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;

	  }

	else if(fdoZVertCorr == 1)
	  {
	    if(zvert1 > fh1Reco1->GetXaxis()->GetXmin() && zvert1 < fh1Reco1->GetXaxis()->GetXmax() && zvert2 > fh1Reco2->GetXaxis()->GetXmin() && zvert2 < fh1Reco2->GetXaxis()->GetXmax())
	      {
		eps1 = fh1Reco1->GetBinContent(fh1Reco1->FindFixBin(zvert1));
		eps2 = fh1Reco2->GetBinContent(fh1Reco2->FindFixBin(zvert2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;
	  }

      }

    else if(boolSum == 2)
      {
	if(fdoPtCorr == 1 && fdoEtaCorr == 1)
	  {
	    if(pT1 > fh2Reco1->GetXaxis()->GetXmin() && pT1 < fh2Reco1->GetXaxis()->GetXmax() && pT2 > fh2Reco2->GetXaxis()->GetXmin() && pT2 < fh2Reco2->GetXaxis()->GetXmax() && eta1 > fh2Reco1->GetYaxis()->GetXmin() && eta1 < fh2Reco1->GetYaxis()->GetXmax() && eta2 > fh2Reco2->GetYaxis()->GetXmin() && eta2 < fh2Reco2->GetYaxis()->GetXmax())
	      {
		eps1 = fh2Reco1->GetBinContent(fh2Reco1->GetXaxis()->FindFixBin(pT1),fh2Reco1->GetYaxis()->FindFixBin(eta1));
		eps2 = fh2Reco2->GetBinContent(fh2Reco2->GetXaxis()->FindFixBin(pT2),fh2Reco2->GetYaxis()->FindFixBin(eta2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;

	  }

	if(fdoPtCorr == 1 && fdoPhiCorr == 1)
	  {

	    if(pT1 > fh2Reco1->GetXaxis()->GetXmin() && pT1 < fh2Reco1->GetXaxis()->GetXmax() && pT2 > fh2Reco2->GetXaxis()->GetXmin() && pT2 < fh2Reco2->GetXaxis()->GetXmax() && phi1 > fh2Reco1->GetYaxis()->GetXmin() && phi1 < fh2Reco1->GetYaxis()->GetXmax() && phi2 > fh2Reco2->GetYaxis()->GetXmin() && phi2 < fh2Reco2->GetYaxis()->GetXmax())
	      {
		eps1 = fh2Reco1->GetBinContent(fh2Reco1->GetXaxis()->FindFixBin(pT1),fh2Reco1->GetYaxis()->FindFixBin(phi1));
		eps2 = fh2Reco2->GetBinContent(fh2Reco2->GetXaxis()->FindFixBin(pT2),fh2Reco2->GetYaxis()->FindFixBin(phi2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;

	  }

	else if(fdoPtCorr == 1 && fdoZVertCorr == 1)
	  {

	    if(pT1 > fh2Reco1->GetXaxis()->GetXmin() && pT1 < fh2Reco1->GetXaxis()->GetXmax() && pT2 > fh2Reco2->GetXaxis()->GetXmin() && pT2 < fh2Reco2->GetXaxis()->GetXmax() && zvert1 > fh2Reco1->GetYaxis()->GetXmin() && zvert1 < fh2Reco1->GetYaxis()->GetXmax() && zvert2 > fh2Reco2->GetYaxis()->GetXmin() && zvert2 < fh2Reco2->GetYaxis()->GetXmax())
	      {
		eps1 = fh2Reco1->GetBinContent(fh2Reco1->GetXaxis()->FindFixBin(pT1),fh2Reco1->GetYaxis()->FindFixBin(zvert1));
		eps2 = fh2Reco2->GetBinContent(fh2Reco2->GetXaxis()->FindFixBin(pT2),fh2Reco2->GetYaxis()->FindFixBin(zvert2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;
	  }
	else if(fdoEtaCorr == 1 && fdoPhiCorr == 1)
	  {

	    if(eta1 > fh2Reco1->GetXaxis()->GetXmin() && eta1 < fh2Reco1->GetXaxis()->GetXmax() && eta2 > fh2Reco2->GetXaxis()->GetXmin() && eta2 < fh2Reco2->GetXaxis()->GetXmax() && phi1 > fh2Reco1->GetYaxis()->GetXmin() && phi1 < fh2Reco1->GetYaxis()->GetXmax() && phi2 > fh2Reco2->GetYaxis()->GetXmin() && phi2 < fh2Reco2->GetYaxis()->GetXmax())
	      {
		eps1 = fh2Reco1->GetBinContent(fh2Reco1->GetXaxis()->FindFixBin(eta1),fh2Reco1->GetYaxis()->FindFixBin(phi1));
		eps2 = fh2Reco2->GetBinContent(fh2Reco2->GetXaxis()->FindFixBin(eta2),fh2Reco2->GetYaxis()->FindFixBin(phi2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;


	  }
	else if(fdoEtaCorr == 1 && fdoZVertCorr == 1)
	  {

	    if(eta1 > fh2Reco1->GetXaxis()->GetXmin() && eta1 < fh2Reco1->GetXaxis()->GetXmax() && eta2 > fh2Reco2->GetXaxis()->GetXmin() && eta2 < fh2Reco2->GetXaxis()->GetXmax() && zvert1 > fh2Reco1->GetYaxis()->GetXmin() && zvert1 < fh2Reco1->GetYaxis()->GetXmax() && zvert2 > fh2Reco2->GetXaxis()->GetXmin() && zvert2 < fh2Reco2->GetXaxis()->GetXmax())
	      {
		eps1 = fh2Reco1->GetBinContent(fh2Reco1->GetXaxis()->FindFixBin(eta1),fh2Reco1->GetYaxis()->FindFixBin(zvert1));
		eps1 = fh2Reco2->GetBinContent(fh2Reco2->GetXaxis()->FindFixBin(eta2),fh2Reco2->GetYaxis()->FindFixBin(zvert2));


		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;

	  }
	else if(fdoPhiCorr == 1 && fdoZVertCorr == 1)
	  {

	    if(phi1 > fh2Reco1->GetXaxis()->GetXmin() && phi1 < fh2Reco1->GetXaxis()->GetXmax() && phi2 > fh2Reco2->GetXaxis()->GetXmin() && phi2 < fh2Reco2->GetXaxis()->GetXmax() && zvert1 > fh2Reco1->GetYaxis()->GetXmin() && zvert1 < fh2Reco1->GetYaxis()->GetXmax() && zvert2 > fh2Reco2->GetYaxis()->GetXmin() && zvert2 < fh2Reco2->GetYaxis()->GetXmax())
	      {
		eps1 = fh2Reco1->GetBinContent(fh2Reco1->GetXaxis()->FindFixBin(phi1),fh2Reco1->GetYaxis()->FindFixBin(zvert1));
		eps2 = fh2Reco2->GetBinContent(fh2Reco2->GetXaxis()->FindFixBin(phi2),fh2Reco2->GetYaxis()->FindFixBin(zvert2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;

	  }
      }


    else if(boolSum == 3)
      {
	if(fdoPtCorr == 1 && fdoEtaCorr == 1 && fdoPhiCorr == 1)
	  {
	    if(pT1 >fh3Reco1->GetXaxis()->GetXmin() && pT1 <fh3Reco1->GetXaxis()->GetXmax() &&
               pT2 > fh3Reco2->GetXaxis()->GetXmin() && pT2 <fh3Reco2->GetXaxis()->GetXmax() &&
               eta1 > fh3Reco1->GetYaxis()->GetXmin() && eta1 <fh3Reco1->GetYaxis()->GetXmax() &&
               eta2 > fh3Reco2->GetYaxis()->GetXmin() && eta2 <fh3Reco2->GetYaxis()->GetXmax() &&
               phi1 > fh3Reco1->GetZaxis()->GetXmin() && phi1 < fh3Reco1->GetZaxis()->GetXmax() &&
               phi2 > fh3Reco2->GetZaxis()->GetXmin() && phi2 < fh3Reco2->GetZaxis()->GetXmax())
	      {
		eps1 = fh3Reco1->GetBinContent(fh3Reco1->GetXaxis()->FindFixBin(pT1),fh3Reco1->GetYaxis()->FindFixBin(eta1),fh3Reco1->GetZaxis()->FindFixBin(phi1));
		eps2 = fh3Reco2->GetBinContent(fh3Reco2->GetXaxis()->FindFixBin(pT2),fh3Reco2->GetYaxis()->FindFixBin(eta2),fh3Reco2->GetZaxis()->FindFixBin(phi2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;


	}

	else if(fdoPtCorr == 1 && fdoEtaCorr == 1 && fdoZVertCorr == 1)
	  {

	    if(pT1 >fh3Reco1->GetXaxis()->GetXmin() && pT1 <fh3Reco1->GetXaxis()->GetXmax() && pT2 > fh3Reco2->GetXaxis()->GetXmin() && pT2 <fh3Reco2->GetXaxis()->GetXmax() && eta1 > fh3Reco1->GetYaxis()->GetXmin() && eta1 <fh3Reco1->GetYaxis()->GetXmax() &&  eta2 > fh3Reco2->GetYaxis()->GetXmin() && eta2 <fh3Reco2->GetYaxis()->GetXmax() &&  zvert1 > fh3Reco1->GetZaxis()->GetXmin() && zvert1 < fh3Reco1->GetZaxis()->GetXmax() &&  zvert2 > fh3Reco2->GetZaxis()->GetXmin() && zvert2 < fh3Reco2->GetZaxis()->GetXmax())
	      {
		eps1 = fh3Reco1->GetBinContent(fh3Reco1->GetXaxis()->FindFixBin(pT1),fh3Reco1->GetYaxis()->FindFixBin(eta1),fh3Reco1->GetZaxis()->FindFixBin(zvert1));
		eps2 = fh3Reco2->GetBinContent(fh3Reco2->GetXaxis()->FindFixBin(pT2),fh3Reco2->GetYaxis()->FindFixBin(eta2),fh3Reco2->GetZaxis()->FindFixBin(zvert2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;
	  }

	else if(fdoPtCorr == 1 && fdoPhiCorr == 1 && fdoZVertCorr == 1)
	  {

	    if(pT1 >fh3Reco1->GetXaxis()->GetXmin() && pT1 <fh3Reco1->GetXaxis()->GetXmax() && pT2 > fh3Reco2->GetXaxis()->GetXmin() && pT2 <fh3Reco2->GetXaxis()->GetXmax() && phi1 > fh3Reco1->GetYaxis()->GetXmin() && phi1 <fh3Reco1->GetYaxis()->GetXmax() &&  phi2 > fh3Reco2->GetYaxis()->GetXmin() && phi2 <fh3Reco2->GetYaxis()->GetXmax() &&  zvert1 > fh3Reco1->GetZaxis()->GetXmin() && zvert1 < fh3Reco1->GetZaxis()->GetXmax() &&  zvert2 > fh3Reco2->GetZaxis()->GetXmin() && zvert2 < fh3Reco2->GetZaxis()->GetXmax())
	      {
		eps1 = fh3Reco1->GetBinContent(fh3Reco1->GetXaxis()->FindFixBin(pT1),fh3Reco1->GetYaxis()->FindFixBin(phi1),fh3Reco1->GetZaxis()->FindFixBin(zvert1));
		eps2 = fh3Reco2->GetBinContent(fh3Reco2->GetXaxis()->FindFixBin(pT2),fh3Reco2->GetYaxis()->FindFixBin(phi2),fh3Reco2->GetZaxis()->FindFixBin(zvert2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;

	  }

	else if(fdoEtaCorr == 1 && fdoPhiCorr == 1 && fdoZVertCorr == 1)
	  {

	    if(eta1 >fh3Reco1->GetXaxis()->GetXmin() && eta1 <fh3Reco1->GetXaxis()->GetXmax() && eta2 > fh3Reco2->GetXaxis()->GetXmin() && eta2 <fh3Reco2->GetXaxis()->GetXmax() && phi1 > fh3Reco1->GetYaxis()->GetXmin() && phi1 <fh3Reco1->GetYaxis()->GetXmax() &&  phi2 > fh3Reco2->GetYaxis()->GetXmin() && phi2 <fh3Reco2->GetYaxis()->GetXmax() &&  zvert1 > fh3Reco1->GetZaxis()->GetXmin() && zvert1 < fh3Reco1->GetZaxis()->GetXmax() &&  zvert2 > fh3Reco2->GetZaxis()->GetXmin() && zvert2 < fh3Reco2->GetZaxis()->GetXmax())
	      {
		eps1 = fh3Reco1->GetBinContent(fh3Reco1->GetXaxis()->FindFixBin(eta1),fh3Reco1->GetYaxis()->FindFixBin(phi1),fh3Reco1->GetZaxis()->FindFixBin(zvert1));
		eps2 = fh3Reco2->GetBinContent(fh3Reco2->GetXaxis()->FindFixBin(eta2),fh3Reco2->GetYaxis()->FindFixBin(phi2),fh3Reco2->GetZaxis()->FindFixBin(zvert2));

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;

		return w1*w2;
	      }
	    else
	      return 0;

	  }
      }

    else if(boolSum == 4)
      {

	if(pT1 > fhntReco1->GetAxis(0)->GetXmin() && pT1 < fhntReco1->GetAxis(0)->GetXmax() && pT2 > fhntReco2->GetAxis(0)->GetXmin() && pT2 < fhntReco2->GetAxis(0)->GetXmax() && eta1 > fhntReco1->GetAxis(1)->GetXmin() && eta1 <fhntReco1->GetAxis(1)->GetXmax() && eta2 > fhntReco2->GetAxis(1)->GetXmin() && eta2 < fhntReco2->GetAxis(1)->GetXmax() && phi1 > fhntReco1->GetAxis(2)->GetXmin() && phi2 < fhntReco2->GetAxis(2)->GetXmax() && phi2 > fhntReco2->GetAxis(2)->GetXmin() && phi2 < fhntReco2->GetAxis(2)->GetXmax() && zvert1 > fhntReco1->GetAxis(3)->GetXmin() && zvert1 < fhntReco1->GetAxis(3)->GetXmax() && zvert2 > fhntReco2->GetAxis(3)->GetXmin() && zvert2 < fhntReco2->GetAxis(3)->GetXmax())
	      {

		int tab1[] = {fhntReco1->GetAxis(0)->FindFixBin(pT1),fhntReco1->GetAxis(1)->FindFixBin(eta1),fhntReco1->GetAxis(2)->FindFixBin(phi1),fhntReco1->GetAxis(3)->FindFixBin(zvert1)};
		int tab2[] = {fhntReco2->GetAxis(0)->FindFixBin(pT2),fhntReco2->GetAxis(1)->FindFixBin(eta2),fhntReco2->GetAxis(2)->FindFixBin(phi2),fhntReco2->GetAxis(3)->FindFixBin(zvert2)};

		eps1 = fhntReco1->GetBinContent(tab1);
		eps2 = fhntReco2->GetBinContent(tab2);

		w1 = (1-cont1)/eps1;
		w2 = (1-cont2)/eps2;
		return w1*w2;

	      }
	    else
	      return 0;

      }

    return 0;

}



double AliFemtoCorrFctnDEtaDPhiCorrections::GetPurity(double pT1, int n)
{
  double w1=0.;
  if(n==1){
    if(pT1 > fSinglePurity1->GetXaxis()->GetXmin() && pT1 < fSinglePurity1->GetXaxis()->GetXmax())
      {
	w1 = fSinglePurity1->GetBinContent(fSinglePurity1->FindFixBin(pT1));
	return w1;
      }
    else
      return 0;
  }
  if(n==2){
    if(pT1 > fSinglePurity2->GetXaxis()->GetXmin() && pT1 < fSinglePurity2->GetXaxis()->GetXmax())
      {
	w1 = fSinglePurity2->GetBinContent(fSinglePurity2->FindFixBin(pT1));
	return w1;
      }
    else
      return 0;
  }
  return 0;
}
