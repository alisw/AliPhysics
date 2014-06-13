////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnPairFractions - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Malgorzata Janik, majanik@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnPairFractions.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctnPairFractions)
#endif
  
#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnPairFractions::AliFemtoCorrFctnPairFractions(char* title):
AliFemtoCorrFctn(),
  fPairFractions(0),
  fphiL(0),
  fphiT(0)
{

  //fphiL = (-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  //fphiT = 2*TMath::Pi()+(-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;

  TString  hname  = "hPairFraction"; hname+= title;
  TString  htitle = "Pair Fraction "; htitle+= title;
  fPairFractions = new TH1F(hname.Data(),htitle.Data(), 7, 0, 7);
  fPairFractions->GetXaxis()->SetBinLabel(1,"#pi#pi, MC");
  fPairFractions->GetXaxis()->SetBinLabel(2,"KK, MC");
  fPairFractions->GetXaxis()->SetBinLabel(3,"pp, MC");
  fPairFractions->GetXaxis()->SetBinLabel(4,"#pi K, MC");
  fPairFractions->GetXaxis()->SetBinLabel(5,"#pi p, MC");
  fPairFractions->GetXaxis()->SetBinLabel(6,"Kp, MC");
  fPairFractions->GetXaxis()->SetBinLabel(7,"Other, MC");


  // to enable error bar calculation...

  fPairFractions->Sumw2();
  fPairFractions->Sumw2();
}

//____________________________
AliFemtoCorrFctnPairFractions::AliFemtoCorrFctnPairFractions(const AliFemtoCorrFctnPairFractions& aCorrFctn) :
  AliFemtoCorrFctn(),
  fPairFractions(0),
  fphiL(0),
  fphiT(0)
{
  // copy constructor
  if (aCorrFctn.fPairFractions)
    fPairFractions = new TH1F(*aCorrFctn.fPairFractions);
  else
    fPairFractions = 0;

 if (aCorrFctn.fPairFractions)
    fPairFractions = new TH1F(*aCorrFctn.fPairFractions);
  else
    fPairFractions = 0;

  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;


}
//____________________________
AliFemtoCorrFctnPairFractions::~AliFemtoCorrFctnPairFractions(){
  // destructor
  delete fPairFractions;
  delete fPairFractions;
}
//_________________________
AliFemtoCorrFctnPairFractions& AliFemtoCorrFctnPairFractions::operator=(const AliFemtoCorrFctnPairFractions& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fPairFractions)
    fPairFractions = new TH1F(*aCorrFctn.fPairFractions);
  else
    fPairFractions = 0;

  
  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;

  return *this;
}
//_________________________
void AliFemtoCorrFctnPairFractions::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  // mShareDenominator->Draw();
  // mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnPairFractions::Report(){
  // create report
  string stemp = "Pair Fractions Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fPairFractions->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fPairFractions->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnPairFractions::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair

  //Applying pair cuts
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;



  Int_t pdg1=0;
  AliFemtoModelHiddenInfo *info1 = ( AliFemtoModelHiddenInfo *) pair->Track1()->GetHiddenInfo();
  if(info1)pdg1 = info1->GetPDGPid();

  Int_t pdg2=0;
  AliFemtoModelHiddenInfo *info2 = ( AliFemtoModelHiddenInfo *) pair->Track2()->GetHiddenInfo();
  if(info2)pdg2 = info2->GetPDGPid();

  if(abs(pdg1)==211 && abs(pdg2)==211) //pi pi
      fPairFractions->Fill(0.5);
  else if(abs(pdg1)==321 && abs(pdg2)==321)// K K
      fPairFractions->Fill(1.5);
  else if(abs(pdg1)==2212 && abs(pdg2)==2212)// p p
      fPairFractions->Fill(2.5);
  else if(abs(pdg1)==211 && abs(pdg2)==321)// pi K
      fPairFractions->Fill(3.5);
  else if(abs(pdg1)==211 && abs(pdg2)==2212)// pi p
      fPairFractions->Fill(4.5);
  else if(abs(pdg1)==321 && abs(pdg2)==2212)//K p
      fPairFractions->Fill(5.5);
  else //other
      fPairFractions->Fill(6.5);

  /*double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

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
   double px2 = pair->Track2()->Track()->P().x();
   double py2 = pair->Track2()->Track()->P().y();
   double pt1 = TMath::Hypot(px1, py1);
   double pt2 = TMath::Hypot(px2, py2);


   double PionMass = 0.13956995;*/
 
}
//____________________________
void AliFemtoCorrFctnPairFractions::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

}


void AliFemtoCorrFctnPairFractions::WriteHistos()
{
  // Write out result histograms
  fPairFractions->Write();
}

TList* AliFemtoCorrFctnPairFractions::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fPairFractions);


  return tOutputList;

}
