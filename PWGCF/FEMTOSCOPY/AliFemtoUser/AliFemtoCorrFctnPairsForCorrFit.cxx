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
  mNtuple(0)
{
  mNtuple = new TNtuple("pair", "pair", "px1:py1:pz1:e1:px2:py2:pz2:e2");
}

//____________________________
AliFemtoCorrFctnPairsForCorrFit::AliFemtoCorrFctnPairsForCorrFit(const AliFemtoCorrFctnPairsForCorrFit& aCorrFctn) :
  AliFemtoCorrFctn(),  
  mNtuple(0)
{
  // copy constructor
  if (mNtuple) delete mNtuple;
  mNtuple = dynamic_cast<TNtuple*>(aCorrFctn.mNtuple->Clone());
}
//____________________________
AliFemtoCorrFctnPairsForCorrFit::~AliFemtoCorrFctnPairsForCorrFit(){
  // destructor
  delete  mNtuple;
}
//_________________________
AliFemtoCorrFctnPairsForCorrFit& AliFemtoCorrFctnPairsForCorrFit::operator=(const AliFemtoCorrFctnPairsForCorrFit& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;
  
  if (mNtuple) delete mNtuple;
   mNtuple = dynamic_cast<TNtuple*>(aCorrFctn.mNtuple->Clone());
  
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

  double PionMass = 0.13957018;//0.13956995;
  double KaonMass = 0.493677;
 
 double px1 = pair->Track1()->Track()->P().x();
 double py1 = pair->Track1()->Track()->P().y();
 double pz1 = pair->Track1()->Track()->P().z();
 
 double px2 = pair->Track2()->Track()->P().x();
 double py2 = pair->Track2()->Track()->P().y();
 double pz2 = pair->Track2()->Track()->P().z();
 
 const AliFemtoThreeVector p1 =  pair->Track1()->Track()->P();
 const AliFemtoThreeVector p2 =  pair->Track2()->Track()->P();
 double e1 = TMath::Sqrt(PionMass*PionMass + p1.Mag2());
 double e2 = TMath::Sqrt(KaonMass*KaonMass + p2.Mag2());

 double tKStar = fabs(pair->KStar());
 
 mNtuple->Fill(px1, py1, pz1, e1, px2, py2, pz2, e2); 
 

}
//____________________________
void AliFemtoCorrFctnPairsForCorrFit::AddMixedPair( AliFemtoPair* pair){

}
//____________________________

void AliFemtoCorrFctnPairsForCorrFit::WriteHistos()
{
  // Write out result histograms
  mNtuple->Write();
}

TList* AliFemtoCorrFctnPairsForCorrFit::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(mNtuple);

  return tOutputList;
}
