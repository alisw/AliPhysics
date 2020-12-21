////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnPairsForCorrFit - A correlation function that analyzes        //
// two particle mass minvariant with various mass assumptions                 //
//                                                                            //
// Authors: Małgorzata Janik majanik@cern.ch
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnPairsForCorrFit.h"

#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnPairsForCorrFit)
#endif

//____________________________
AliFemtoCorrFctnPairsForCorrFit::AliFemtoCorrFctnPairsForCorrFit(const char* title):
AliFemtoCorrFctn(),
  mNtuple(0),
  hKstar(0)
{
  mNtuple = new TNtuple(Form("pair%s",title), Form("pair%s",title), "px1:py1:pz1:e1:px2:py2:pz2:e2");
  hKstar = new TH1F(Form("kstar%s",title), Form("kstar%s",title),400,0,2);
}

//____________________________
AliFemtoCorrFctnPairsForCorrFit::AliFemtoCorrFctnPairsForCorrFit(const AliFemtoCorrFctnPairsForCorrFit& aCorrFctn) :
  AliFemtoCorrFctn(),  
  mNtuple(0),
  hKstar(0)
{
  // copy constructor
  if (mNtuple) delete mNtuple;
  mNtuple = dynamic_cast<TNtuple*>(aCorrFctn.mNtuple->Clone());

  if (hKstar) delete hKstar;
  hKstar = dynamic_cast<TH1F*>(aCorrFctn.hKstar->Clone());
}
//____________________________
AliFemtoCorrFctnPairsForCorrFit::~AliFemtoCorrFctnPairsForCorrFit(){
  // destructor
  delete  mNtuple;
  delete hKstar;
}
//_________________________
AliFemtoCorrFctnPairsForCorrFit& AliFemtoCorrFctnPairsForCorrFit::operator=(const AliFemtoCorrFctnPairsForCorrFit& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;
  
  if (mNtuple) delete mNtuple;
   mNtuple = dynamic_cast<TNtuple*>(aCorrFctn.mNtuple->Clone());
   
  if (hKstar) delete hKstar;
  hKstar = dynamic_cast<TH1F*>(aCorrFctn.hKstar->Clone());
  
  return *this;
}
//_________________________
void AliFemtoCorrFctnPairsForCorrFit::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnPairsForCorrFit::Report(){
  // create report
  string stemp = "Kstar vs Pt Monitor Report\n";
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnPairsForCorrFit::AddRealPair( AliFemtoPair* pair){

 
 
}
//____________________________
void AliFemtoCorrFctnPairsForCorrFit::AddMixedPair( AliFemtoPair* pair){

  double px1 = pair->Track1()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  double pz1 = pair->Track1()->Track()->P().z();
  Double_t e1 = pair->Track1()->FourMomentum().e();
 
  double px2 = pair->Track2()->Track()->P().x();
  double py2 = pair->Track2()->Track()->P().y();
  double pz2 = pair->Track2()->Track()->P().z();
  Double_t e2 = pair->Track2()->FourMomentum().e();

  double tKStar = fabs(pair->KStar());
  hKstar->Fill(tKStar);
  int bin = hKstar->FindBin(tKStar);
  if(tKStar<0.21)
    mNtuple->Fill(px1, py1, pz1, e1, px2, py2, pz2, e2); 
 

  
}
//____________________________

void AliFemtoCorrFctnPairsForCorrFit::WriteHistos()
{
  // Write out result histograms
  mNtuple->Write();
  hKstar->Write();
}

TList* AliFemtoCorrFctnPairsForCorrFit::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(mNtuple);
  tOutputList->Add(hKstar);
  return tOutputList;
}
