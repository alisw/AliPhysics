////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhiTHn - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Malgorzata Janik majanik@cern.ch                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDEtaDPhiTHn.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctnDEtaDPhiTHn)
#endif
  
#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnDEtaDPhiTHn::AliFemtoCorrFctnDEtaDPhiTHn(char* title, const int& aPhiBins=20, const int& aEtaBins=20, const int &pT1Bins=1, const double& pT1min=0, const double& pT1max=4, const int &pT2Bins=1, const double& pT2min=0, const double& pT2max=4, const int &zvtxBins=10, const double& zvtxmin=-10, const double& zvtxmax=10, const int &multBins=5, const int& multmin=0, const int& multmax=100):
AliFemtoCorrFctn(),
  fDPhiDEtaNum(0),
  fDPhiDEtaDen(0),
  fphiL(0),
  fphiT(0),
  fPt1Min(pT1min),
  fPt1Max(pT1max),
  fPt2Min(pT2min),
  fPt2Max(pT2max),
  fZvtxMin(zvtxmin),
  fZvtxMax(zvtxmax),
  fMultMin(multmin),
  fMultMax(multmax)  
{
  fphiL = (-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  fphiT = 2*TMath::Pi()+(-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;

 
  // THnSparse(const char* name, const char* title, Int_t dim,
  //           const Int_t* nbins, const Double_t* xmin, const Double_t* xmax,
  //           Int_t chunksize);

 
  const Int_t nbins[] = {aPhiBins,aEtaBins,pT1Bins,pT2Bins,multBins,zvtxBins};
  const Double_t xmin[] = {fphiL,-2.0 ,pT1min,pT2min,fMultMin,zvtxmin};
  const Double_t xmax[] = {fphiT, 2.0 ,pT1max,pT2max,fMultMax,zvtxmax};
  
  // set up numerator
  char tTitNumD[101] = "NumDPhiDEta";
  strncat(tTitNumD,title, 100);
  fDPhiDEtaNum = new THnSparseF(tTitNumD,title,6,nbins,xmin,xmax);
  
  // set up denominator
  char tTitDenD[101] = "DenDPhiDEta";
  strncat(tTitDenD,title, 100);
  fDPhiDEtaDen = new THnSparseF(tTitDenD,title,6,nbins,xmin,xmax);

  // to enable error bar calculation...
  fDPhiDEtaNum->Sumw2();
  fDPhiDEtaDen->Sumw2();
}

//____________________________
AliFemtoCorrFctnDEtaDPhiTHn::AliFemtoCorrFctnDEtaDPhiTHn(const AliFemtoCorrFctnDEtaDPhiTHn& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiDEtaNum(0),
  fDPhiDEtaDen(0),
  fphiL(0),
  fphiT(0),
  fPt1Min(aCorrFctn.fPt1Min),
  fPt1Max(aCorrFctn.fPt1Max),
  fPt2Min(aCorrFctn.fPt2Min),
  fPt2Max(aCorrFctn.fPt2Max),
  fZvtxMin(aCorrFctn.fZvtxMin),
  fZvtxMax(aCorrFctn.fZvtxMax),
  fMultMin(aCorrFctn.fMultMin),
  fMultMax(aCorrFctn.fMultMax)  
{
  // copy constructor
  /*
  if (aCorrFctn.fDPhiDEtaNum)
    fDPhiDEtaNum = new THnSparseF(*aCorrFctn.fDPhiDEtaNum);
  else
    fDPhiDEtaNum = 0;
  
  if (aCorrFctn.fDPhiDEtaDen)
    fDPhiDEtaDen = new THnSparseF(*aCorrFctn.fDPhiDEtaDen);
  else
    fDPhiDEtaDen = 0;
  */
    char title[]={"bad constructor"};
    Int_t aPhiBins=20;
    Int_t   aEtaBins=20;
    Int_t   pT1Bins=1;
    Int_t   pT2Bins=1;
    Int_t   multBins=5;
    Int_t   zvtxBins=10;
   
    const Int_t nbins[] = {aPhiBins,aEtaBins,pT1Bins,pT2Bins,multBins,zvtxBins};
    const Double_t xmin[] = {fphiL,-2,0,0,0,-10};
    const Double_t xmax[] = {fphiT,2,4,4,100,10};
  
    // set up numerator
    char tTitNumD[101] = "NumDPhiDEta";
    strncat(tTitNumD,title, 100);
    fDPhiDEtaNum = new THnSparseF(tTitNumD,title,6,nbins,xmin,xmax);
  
    // set up denominator
    char tTitDenD[101] = "DenDPhiDEta";
    strncat(tTitDenD,title, 100);
    fDPhiDEtaDen = new THnSparseF(tTitDenD,title,6,nbins,xmin,xmax);
     
    
  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;


}
//____________________________
AliFemtoCorrFctnDEtaDPhiTHn::~AliFemtoCorrFctnDEtaDPhiTHn(){
  // destructor

  delete fDPhiDEtaNum;
  delete fDPhiDEtaDen;
}
//_________________________
AliFemtoCorrFctnDEtaDPhiTHn& AliFemtoCorrFctnDEtaDPhiTHn::operator=(const AliFemtoCorrFctnDEtaDPhiTHn& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;
    char title[]={"bad constructor"};
    Int_t aPhiBins=20;
    Int_t   aEtaBins=20;
    Int_t   pT1Bins=1;
    Int_t   pT2Bins=1;
    Int_t   multBins=5;
    Int_t   zvtxBins=10;
   
    const Int_t nbins[] = {aPhiBins,aEtaBins,pT1Bins,pT2Bins,multBins,zvtxBins};
    const Double_t xmin[] = {fphiL,-2,0,0,0,-10};
    const Double_t xmax[] = {fphiT,2,4,4,100,10};
  
    // set up numerator
    char tTitNumD[101] = "NumDPhiDEta";
    strncat(tTitNumD,title, 100);
    fDPhiDEtaNum = new THnSparseF(tTitNumD,title,6,nbins,xmin,xmax);
  
    // set up denominator
    char tTitDenD[101] = "DenDPhiDEta";
    strncat(tTitDenD,title, 100);
    fDPhiDEtaDen = new THnSparseF(tTitDenD,title,6,nbins,xmin,xmax);


      /*
  if (aCorrFctn.fDPhiDEtaNum)
    fDPhiDEtaNum = new THnSparseF(*aCorrFctn.fDPhiDEtaNum);
  else
    fDPhiDEtaNum = 0;

  if (aCorrFctn.fDPhiDEtaDen)
    fDPhiDEtaDen = new THnSparseF(*aCorrFctn.fDPhiDEtaDen);
  else
    fDPhiDEtaDen = 0;
  */

  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;

  return *this;
}
//_________________________
void AliFemtoCorrFctnDEtaDPhiTHn::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDEtaDPhiTHn::Report(){
  // create report
  string stemp = "TPC Ncls Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fDPhiDEtaNum->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDPhiDEtaDen->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}

//____________________________
void AliFemtoCorrFctnDEtaDPhiTHn::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

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

  double px1 = pair->Track1()->FourMomentum()[0];
  double py1 = pair->Track1()->FourMomentum()[1];
  //double pz1 = pair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->FourMomentum()[0];
  double py2 = pair->Track2()->FourMomentum()[1];
  //double pz2 = pair->Track2()->Track()->P().z();

  double pt1 = TMath::Hypot(px1, py1);
  double pt2 = TMath::Hypot(px2, py2);
  //   double ptmin = pt1>pt2 ? pt2 : pt1;

  double mult = -1;
  double zvtx = -11;
  if(pair->Track1()->Track()){
    mult = pair->Track1()->Track()->Multiplicity();
    zvtx = pair->Track1()->Track()->Zvtx();
  }

  else if(pair->Track1()->V0()){
    mult = pair->Track1()->V0()->Multiplicity();
    zvtx = pair->Track1()->V0()->Zvtx();
  }

   
  //fDPhiDEtaNumerator->Fill(dphi, deta);
  //cout<<dphi<<" "<<deta<<" "<<pt1<<" "<<pt2<<" "<<mult<<" "<<zvtx<<endl;
  Double_t value[] = {dphi, deta, pt1, pt2, mult, zvtx};
  fDPhiDEtaNum->Fill(value);

}

//____________________________
void AliFemtoCorrFctnDEtaDPhiTHn::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

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

  double px1 = pair->Track1()->FourMomentum()[0];
  double py1 = pair->Track1()->FourMomentum()[1];
  //double pz1 = pair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->FourMomentum()[0];
  double py2 = pair->Track2()->FourMomentum()[1];
  //double pz2 = pair->Track2()->Track()->P().z();

  double pt1 = TMath::Hypot(px1, py1);
  double pt2 = TMath::Hypot(px2, py2);
  //   double ptmin = pt1>pt2 ? pt2 : pt1;

  double mult = -1;
  double zvtx = -11;
  if(pair->Track1()->Track()){
    mult = pair->Track1()->Track()->Multiplicity();
    zvtx = pair->Track1()->Track()->Zvtx();
  }

  else if(pair->Track1()->V0()){
    mult = pair->Track1()->V0()->Multiplicity();
    zvtx = pair->Track1()->V0()->Zvtx();
  }


  //fDPhiDEtaNumerator->Fill(dphi, deta);
  Double_t value[] = {dphi, deta, pt1, pt2, mult, zvtx};
  fDPhiDEtaDen->Fill(value);
}


void AliFemtoCorrFctnDEtaDPhiTHn::WriteHistos()
{
  // Write out result histograms
  fDPhiDEtaNum->Write();
  fDPhiDEtaDen->Write();
}

TList* AliFemtoCorrFctnDEtaDPhiTHn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();
  tOutputList->Add(fDPhiDEtaNum);
  tOutputList->Add(fDPhiDEtaDen);
  return tOutputList;

}
