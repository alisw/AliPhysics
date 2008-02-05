////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoTPCInnerCorrFctn - A correlation function that saves the         ///
/// distance at the entrance to the TPC between two tracks as a function     ///
/// of qinv                                                                  ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoTPCInnerCorrFctn.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoTPCInnerCorrFctn)
#endif

//____________________________
AliFemtoTPCInnerCorrFctn::AliFemtoTPCInnerCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  fDTPCNumerator(0),
  fDTPCDenominator(0)
{
  // set up numerator
  //  title = "Num Qinv (MeV/c)";
  char tTitNum[100] = "NumDTPC";
  strcat(tTitNum,title);
  fDTPCNumerator = new TH2D(tTitNum,title,nbins,QinvLo,QinvHi,100,0.0,20.0);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTitDen[100] = "DenDTPC";
  strcat(tTitDen,title);
  fDTPCDenominator = new TH2D(tTitDen,title,nbins,QinvLo,QinvHi,100,0.0,20.0);

  // to enable error bar calculation...
  fDTPCNumerator->Sumw2();
  fDTPCDenominator->Sumw2();
}

//____________________________
AliFemtoTPCInnerCorrFctn::AliFemtoTPCInnerCorrFctn(const AliFemtoTPCInnerCorrFctn& aCorrFctn) :
  fDTPCNumerator(0),
  fDTPCDenominator(0)
{
  // copy constructor
  if (aCorrFctn.fDTPCNumerator)
    fDTPCNumerator = new TH2D(*aCorrFctn.fDTPCNumerator);
  if (aCorrFctn.fDTPCDenominator)
    fDTPCDenominator = new TH2D(*aCorrFctn.fDTPCDenominator);
}
//____________________________
AliFemtoTPCInnerCorrFctn::~AliFemtoTPCInnerCorrFctn(){
  // destructor
  delete fDTPCNumerator;
  delete fDTPCDenominator;
}
//_________________________
AliFemtoTPCInnerCorrFctn& AliFemtoTPCInnerCorrFctn::operator=(const AliFemtoTPCInnerCorrFctn& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fDTPCNumerator)
    fDTPCNumerator = new TH2D(*aCorrFctn.fDTPCNumerator);
  else
    fDTPCNumerator = 0;
  if (aCorrFctn.fDTPCDenominator)
    fDTPCDenominator = new TH2D(*aCorrFctn.fDTPCDenominator);
  else
    fDTPCDenominator = 0;

  return *this;
}
//_________________________
void AliFemtoTPCInnerCorrFctn::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoTPCInnerCorrFctn::Report(){
  // create report
  string stemp = "Entrace TPC distance Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fDTPCNumerator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDTPCDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoTPCInnerCorrFctn::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
  double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
  double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
  double dist = sqrt(distx*distx + disty*disty + distz*distz);

  fDTPCNumerator->Fill(tQinv, dist);
//   cout << "AliFemtoTPCInnerCorrFctn::AddRealPair : " << tQinv << " " << dist << endl;
//   cout << distx << " " << disty << " " << distz << endl;
}
//____________________________
void AliFemtoTPCInnerCorrFctn::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
  double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
  double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
  double dist = sqrt(distx*distx + disty*disty + distz*distz);

  fDTPCDenominator->Fill(tQinv,dist);
}


void AliFemtoTPCInnerCorrFctn::WriteHistos()
{
  // Write out result histograms
  fDTPCNumerator->Write();
  fDTPCDenominator->Write();
}
//______________________________
TList* AliFemtoTPCInnerCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDTPCNumerator); 
  tOutputList->Add(fDTPCDenominator);  

  return tOutputList;
}
