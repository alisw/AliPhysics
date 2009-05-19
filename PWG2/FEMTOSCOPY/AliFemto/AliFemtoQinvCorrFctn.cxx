///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoQinvCorrFctn:                                                 //
// a simple Q-invariant correlation function                             // 
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoQinvCorrFctn.h"
//#include "AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoQinvCorrFctn)
#endif

//____________________________
AliFemtoQinvCorrFctn::AliFemtoQinvCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  fNumerator(0),
  fDenominator(0),
  fRatio(0)
{
  // set up numerator
  //  title = "Num Qinv (MeV/c)";
  char tTitNum[100] = "Num";
  strcat(tTitNum,title);
  fNumerator = new TH1D(tTitNum,title,nbins,QinvLo,QinvHi);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTitDen[100] = "Den";
  strcat(tTitDen,title);
  fDenominator = new TH1D(tTitDen,title,nbins,QinvLo,QinvHi);
  // set up ratio
  //title = "Ratio Qinv (MeV/c)";
  char tTitRat[100] = "Rat";
  strcat(tTitRat,title);
  fRatio = new TH1D(tTitRat,title,nbins,QinvLo,QinvHi);
  // this next bit is unfortunately needed so that we can have many histos of same "title"
  // it is neccessary if we typedef TH1D to TH1d (which we do)
  //fNumerator->SetDirectory(0);
  //fDenominator->SetDirectory(0);
  //fRatio->SetDirectory(0);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fRatio->Sumw2();

}

//____________________________
AliFemtoQinvCorrFctn::AliFemtoQinvCorrFctn(const AliFemtoQinvCorrFctn& aCorrFctn) :
  AliFemtoCorrFctn(),
  fNumerator(0),
  fDenominator(0),
  fRatio(0)
{
  // copy constructor
  fNumerator = new TH1D(*aCorrFctn.fNumerator);
  fDenominator = new TH1D(*aCorrFctn.fDenominator);
  fRatio = new TH1D(*aCorrFctn.fRatio);
}
//____________________________
AliFemtoQinvCorrFctn::~AliFemtoQinvCorrFctn(){
  // destructor
  delete fNumerator;
  delete fDenominator;
  delete fRatio;
}
//_________________________
AliFemtoQinvCorrFctn& AliFemtoQinvCorrFctn::operator=(const AliFemtoQinvCorrFctn& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (fNumerator) delete fNumerator;
  fNumerator = new TH1D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH1D(*aCorrFctn.fDenominator);
  if (fRatio) delete fRatio;
  fRatio = new TH1D(*aCorrFctn.fRatio);

  return *this;
}

//_________________________
void AliFemtoQinvCorrFctn::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  fNumerator->Draw();
  //fDenominator->Draw();
  //fRatio->Draw();
  fRatio->Divide(fNumerator,fDenominator,1.0,1.0);

}

//____________________________
AliFemtoString AliFemtoQinvCorrFctn::Report(){
  // construct report
  string stemp = "Qinv Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fNumerator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in ratio:\t%E\n",fRatio->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoQinvCorrFctn::AddRealPair(AliFemtoPair* pair){
  // add true pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;
  
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  fNumerator->Fill(tQinv);
  //  cout << "AliFemtoQinvCorrFctn::AddRealPair : " << pair->qInv() << " " << tQinv <<
  //" " << pair->track1().FourMomentum() << " " << pair->track2().FourMomentum() << endl;
}
//____________________________
void AliFemtoQinvCorrFctn::AddMixedPair(AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;
  
  double weight = 1.0;
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  fDenominator->Fill(tQinv,weight);
}
//____________________________
void AliFemtoQinvCorrFctn::Write(){
  // Write out neccessary objects
  fNumerator->Write(); 
  fDenominator->Write();  
}
//______________________________
TList* AliFemtoQinvCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator); 
  tOutputList->Add(fDenominator);  

  return tOutputList;
}


