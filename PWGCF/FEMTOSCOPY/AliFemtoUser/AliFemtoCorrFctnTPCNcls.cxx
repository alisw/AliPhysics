////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnTPCNcls - A correlation function that saves the correlation//
// function as a function of number of TPC clusters of the track              //
//                                                                            //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnTPCNcls.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctnTPCNcls)
#endif

//____________________________
AliFemtoCorrFctnTPCNcls::AliFemtoCorrFctnTPCNcls(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  AliFemtoCorrFctn(),
  fNclsTPCMinNumerator(0),
  fNclsTPCMinDenominator(0)
{
  // set up numerator
  char tTitNum[101] = "NumNclsTPCMin";
  strncat(tTitNum,title, 100);
  fNclsTPCMinNumerator = new TH2D(tTitNum,title,nbins,QinvLo,QinvHi,159,0.5,159.5);
  // set up denominator
  char tTitDen[101] = "DenNclsTPCMin";
  strncat(tTitDen,title, 100);
  fNclsTPCMinDenominator = new TH2D(tTitDen,title,nbins,QinvLo,QinvHi,159,0.5,159.5);

  // to enable error bar calculation...
  fNclsTPCMinNumerator->Sumw2();
  fNclsTPCMinDenominator->Sumw2();
}

//____________________________
AliFemtoCorrFctnTPCNcls::AliFemtoCorrFctnTPCNcls(const AliFemtoCorrFctnTPCNcls& aCorrFctn) :
  AliFemtoCorrFctn(),
  fNclsTPCMinNumerator(0),
  fNclsTPCMinDenominator(0)
{
  // copy constructor
  if (aCorrFctn.fNclsTPCMinNumerator)
    fNclsTPCMinNumerator = new TH2D(*aCorrFctn.fNclsTPCMinNumerator);
  if (aCorrFctn.fNclsTPCMinDenominator)
    fNclsTPCMinDenominator = new TH2D(*aCorrFctn.fNclsTPCMinDenominator);
}
//____________________________
AliFemtoCorrFctnTPCNcls::~AliFemtoCorrFctnTPCNcls(){
  // destructor
  delete fNclsTPCMinNumerator;
  delete fNclsTPCMinDenominator;
}
//_________________________
AliFemtoCorrFctnTPCNcls& AliFemtoCorrFctnTPCNcls::operator=(const AliFemtoCorrFctnTPCNcls& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fNclsTPCMinNumerator)
    fNclsTPCMinNumerator = new TH2D(*aCorrFctn.fNclsTPCMinNumerator);
  else
    fNclsTPCMinNumerator = 0;
  if (aCorrFctn.fNclsTPCMinDenominator)
    fNclsTPCMinDenominator = new TH2D(*aCorrFctn.fNclsTPCMinDenominator);
  else
    fNclsTPCMinDenominator = 0;

  return *this;
}
//_________________________
void AliFemtoCorrFctnTPCNcls::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnTPCNcls::Report(){
  // create report
  string stemp = "TPC Ncls Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fNclsTPCMinNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fNclsTPCMinDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnTPCNcls::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...

  if (pair->Track1()->Track()->TPCncls()>pair->Track2()->Track()->TPCncls())
    fNclsTPCMinNumerator->Fill(tQinv, pair->Track2()->Track()->TPCncls());
  else
    fNclsTPCMinNumerator->Fill(tQinv, pair->Track1()->Track()->TPCncls());
}
//____________________________
void AliFemtoCorrFctnTPCNcls::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...

  if (pair->Track1()->Track()->TPCncls()>pair->Track2()->Track()->TPCncls())
    fNclsTPCMinDenominator->Fill(tQinv, pair->Track2()->Track()->TPCncls());
  else
    fNclsTPCMinDenominator->Fill(tQinv, pair->Track1()->Track()->TPCncls());
}


void AliFemtoCorrFctnTPCNcls::WriteHistos()
{
  // Write out result histograms
  fNclsTPCMinNumerator->Write();
  fNclsTPCMinDenominator->Write();
}

TList* AliFemtoCorrFctnTPCNcls::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNclsTPCMinNumerator);
  tOutputList->Add(fNclsTPCMinDenominator);
  
  return tOutputList;

}
