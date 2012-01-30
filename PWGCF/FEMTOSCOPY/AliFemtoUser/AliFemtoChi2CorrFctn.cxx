////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoChi2CorrFctn - A correlation function that saves the correlation ///
/// function as a function of single track quality (chi2/ndof) for its and   ///
/// tpc                                                                      ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoChi2CorrFctn.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoChi2CorrFctn)
#endif

//____________________________
AliFemtoChi2CorrFctn::AliFemtoChi2CorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  AliFemtoCorrFctn(),
  fChi2ITSSUMNumerator(0),
  fChi2ITSSUMDenominator(0),
  fChi2TPCSUMNumerator(0),
  fChi2TPCSUMDenominator(0),
  fChi2ITSONENumerator(0),
  fChi2ITSONEDenominator(0),
  fChi2TPCONENumerator(0),
  fChi2TPCONEDenominator(0),
  fSigmaToVertexNumerator(0),
  fSigmaToVertexDenominator(0)
{
  // set up numerator
  char tTitNum[101] = "NumChi2ITSSUM";
  strncat(tTitNum,title, 100);
  fChi2ITSSUMNumerator = new TH2D(tTitNum,title,nbins,QinvLo,QinvHi,200,0.0,6.0);
  // set up denominator
  char tTitDen[101] = "DenChi2ITSSUM";
  strncat(tTitDen,title, 100);
  fChi2ITSSUMDenominator = new TH2D(tTitDen,title,nbins,QinvLo,QinvHi,200,0.0,6.0);

  // set up numerator
  char tTit2Num[101] = "NumChi2TPCSUM";
  strncat(tTit2Num,title, 100);
  fChi2TPCSUMNumerator = new TH2D(tTit2Num,title,nbins,QinvLo,QinvHi,200,0.0,6.0);
  // set up denominator
  char tTit2Den[101] = "DenChi2TPCSUM";
  strncat(tTit2Den,title, 100);
  fChi2TPCSUMDenominator = new TH2D(tTit2Den,title,nbins,QinvLo,QinvHi,200,0.0,6.0);

  // to enable error bar calculation...
  fChi2ITSSUMNumerator->Sumw2();
  fChi2ITSSUMDenominator->Sumw2();

  fChi2TPCSUMNumerator->Sumw2();
  fChi2TPCSUMDenominator->Sumw2();
  // set up numerator
  snprintf(tTitNum , 100, "%s%s","NumChi2ITSONE",title);
  fChi2ITSONENumerator = new TH2D(tTitNum,title,nbins,QinvLo,QinvHi,200,0.0,6.0);
  // set up denominator
  snprintf(tTitDen , 100, "%s%s", "DenChi2ITSONE", title);
  fChi2ITSONEDenominator = new TH2D(tTitDen,title,nbins,QinvLo,QinvHi,200,0.0,6.0);

  // set up numerator
  snprintf(tTit2Num , 100, "%s%s","NumChi2TPCONE",title);
  fChi2TPCONENumerator = new TH2D(tTit2Num,title,nbins,QinvLo,QinvHi,200,0.0,6.0);
  // set up denominator
  snprintf(tTit2Den , 100, "%s%s", "DenChi2TPCONE", title);
  fChi2TPCONEDenominator = new TH2D(tTit2Den,title,nbins,QinvLo,QinvHi,200,0.0,6.0);

  // set up numerator
  snprintf(tTit2Num , 100, "%s%s","NumSigmaToVertex",title);
  fSigmaToVertexNumerator = new TH2D(tTit2Num,title,nbins,QinvLo,QinvHi,200,0.0,6.0);
  // set up denominator
  snprintf(tTit2Den , 100, "%s%s", "DenSigmaToVertex", title);
  fSigmaToVertexDenominator = new TH2D(tTit2Den,title,nbins,QinvLo,QinvHi,200,0.0,6.0);

  // to enable error bar calculation...
  fChi2ITSONENumerator->Sumw2();
  fChi2ITSONEDenominator->Sumw2();

  fChi2TPCONENumerator->Sumw2();
  fChi2TPCONEDenominator->Sumw2();
}

//____________________________
AliFemtoChi2CorrFctn::AliFemtoChi2CorrFctn(const AliFemtoChi2CorrFctn& aCorrFctn) :
  AliFemtoCorrFctn(),
  fChi2ITSSUMNumerator(0),
  fChi2ITSSUMDenominator(0),
  fChi2TPCSUMNumerator(0),
  fChi2TPCSUMDenominator(0),
  fChi2ITSONENumerator(0),
  fChi2ITSONEDenominator(0),
  fChi2TPCONENumerator(0),
  fChi2TPCONEDenominator(0),
  fSigmaToVertexNumerator(0),
  fSigmaToVertexDenominator(0)
{
  // copy constructor
  if (aCorrFctn.fChi2ITSSUMNumerator)
    fChi2ITSSUMNumerator = new TH2D(*aCorrFctn.fChi2ITSSUMNumerator);
  if (aCorrFctn.fChi2ITSSUMDenominator)
    fChi2ITSSUMDenominator = new TH2D(*aCorrFctn.fChi2ITSSUMDenominator);
  if (aCorrFctn.fChi2TPCSUMNumerator)
    fChi2TPCSUMNumerator = new TH2D(*aCorrFctn.fChi2TPCSUMNumerator);
  if (aCorrFctn.fChi2TPCSUMDenominator)
    fChi2TPCSUMDenominator = new TH2D(*aCorrFctn.fChi2TPCSUMDenominator);
  if (aCorrFctn.fChi2ITSONENumerator)
    fChi2ITSONENumerator = new TH2D(*aCorrFctn.fChi2ITSONENumerator);
  if (aCorrFctn.fChi2ITSONEDenominator)
    fChi2ITSONEDenominator = new TH2D(*aCorrFctn.fChi2ITSONEDenominator);
  if (aCorrFctn.fChi2TPCONENumerator)
    fChi2TPCONENumerator = new TH2D(*aCorrFctn.fChi2TPCONENumerator);
  if (aCorrFctn.fChi2TPCONEDenominator)
    fChi2TPCONEDenominator = new TH2D(*aCorrFctn.fChi2TPCONEDenominator);
  if (aCorrFctn.fSigmaToVertexNumerator)
    fSigmaToVertexNumerator = new TH2D(*aCorrFctn.fSigmaToVertexNumerator);
  if (aCorrFctn.fSigmaToVertexDenominator)
    fSigmaToVertexDenominator = new TH2D(*aCorrFctn.fSigmaToVertexDenominator);
}
//____________________________
AliFemtoChi2CorrFctn::~AliFemtoChi2CorrFctn(){
  // destructor
  delete fChi2ITSSUMNumerator;
  delete fChi2ITSSUMDenominator;
  delete fChi2TPCSUMNumerator;
  delete fChi2TPCSUMDenominator;
  delete fChi2ITSONENumerator;
  delete fChi2ITSONEDenominator;
  delete fChi2TPCONENumerator;
  delete fChi2TPCONEDenominator;
  delete fSigmaToVertexNumerator;
  delete fSigmaToVertexDenominator;
}
//_________________________
AliFemtoChi2CorrFctn& AliFemtoChi2CorrFctn::operator=(const AliFemtoChi2CorrFctn& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fChi2ITSSUMNumerator)
    fChi2ITSSUMNumerator = new TH2D(*aCorrFctn.fChi2ITSSUMNumerator);
  else
    fChi2ITSSUMNumerator = 0;
  if (aCorrFctn.fChi2ITSSUMDenominator)
    fChi2ITSSUMDenominator = new TH2D(*aCorrFctn.fChi2ITSSUMDenominator);
  else
    fChi2ITSSUMDenominator = 0;
  if (aCorrFctn.fChi2TPCSUMNumerator)
    fChi2TPCSUMNumerator = new TH2D(*aCorrFctn.fChi2TPCSUMNumerator);
  else
    fChi2TPCSUMNumerator = 0;
  if (aCorrFctn.fChi2TPCSUMDenominator)
    fChi2TPCSUMDenominator = new TH2D(*aCorrFctn.fChi2TPCSUMDenominator);
  else
    fChi2TPCSUMDenominator = 0;
  if (aCorrFctn.fChi2ITSONENumerator)
    fChi2ITSONENumerator = new TH2D(*aCorrFctn.fChi2ITSONENumerator);
  else
    fChi2ITSONENumerator = 0;
  if (aCorrFctn.fChi2ITSONEDenominator)
    fChi2ITSONEDenominator = new TH2D(*aCorrFctn.fChi2ITSONEDenominator);
  else
    fChi2ITSONEDenominator = 0;
  if (aCorrFctn.fChi2TPCONENumerator)
    fChi2TPCONENumerator = new TH2D(*aCorrFctn.fChi2TPCONENumerator);
  else
    fChi2TPCONENumerator = 0;
  if (aCorrFctn.fChi2TPCONEDenominator)
    fChi2TPCONEDenominator = new TH2D(*aCorrFctn.fChi2TPCONEDenominator);
  else
    fChi2TPCONEDenominator = 0;
  if (aCorrFctn.fSigmaToVertexNumerator)
    fSigmaToVertexNumerator = new TH2D(*aCorrFctn.fSigmaToVertexNumerator);
  else
    fSigmaToVertexNumerator = 0;
  if (aCorrFctn.fSigmaToVertexDenominator)
    fSigmaToVertexDenominator = new TH2D(*aCorrFctn.fSigmaToVertexDenominator);
  else
    fSigmaToVertexDenominator = 0;

  return *this;
}
//_________________________
void AliFemtoChi2CorrFctn::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoChi2CorrFctn::Report(){
  // create report
  string stemp = "ITS and TPC quality Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fChi2ITSSUMNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fChi2ITSSUMDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoChi2CorrFctn::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...

  if ((pair->Track1()->Track()->ITSncls() == 0) && (pair->Track2()->Track()->ITSncls() == 0))
    fChi2ITSSUMNumerator->Fill(tQinv, 1000.0);
  else 
    fChi2ITSSUMNumerator->Fill(tQinv, 
			       (pair->Track1()->Track()->ITSchi2() + 
				pair->Track2()->Track()->ITSchi2())/
			       (pair->Track1()->Track()->ITSncls() +
				pair->Track2()->Track()->ITSncls()));
  if ((pair->Track1()->Track()->TPCncls() == 0) && (pair->Track2()->Track()->TPCncls() == 0))
    fChi2TPCSUMNumerator->Fill(tQinv, 1000.0);
  else
    fChi2TPCSUMNumerator->Fill(tQinv, 
			       (pair->Track1()->Track()->TPCchi2() +
				pair->Track2()->Track()->TPCchi2())/
			       (pair->Track1()->Track()->TPCncls() +
				pair->Track2()->Track()->TPCncls()));
  double chi2perpointITS1, chi2perpointITS2;
  if (pair->Track1()->Track()->ITSncls() == 0)
    chi2perpointITS1 = 1000.0;
  else
    chi2perpointITS1 = pair->Track1()->Track()->ITSchi2()/pair->Track1()->Track()->ITSncls();

  if (pair->Track2()->Track()->ITSncls() == 0)
    chi2perpointITS2 = 1000.0;
  else
    chi2perpointITS2 = pair->Track2()->Track()->ITSchi2()/pair->Track2()->Track()->ITSncls();


  if (chi2perpointITS1 > chi2perpointITS2) {
    fChi2ITSONENumerator->Fill(tQinv, chi2perpointITS1);
  }
  else {
    fChi2ITSONENumerator->Fill(tQinv, chi2perpointITS2);
  }

  double chi2perpointTPC1, chi2perpointTPC2;
  if (pair->Track1()->Track()->TPCncls() == 0)
    chi2perpointTPC1 = 1000.0;
  else
    chi2perpointTPC1 = pair->Track1()->Track()->TPCchi2()/pair->Track1()->Track()->TPCncls();

  if (pair->Track2()->Track()->TPCncls() == 0)
    chi2perpointTPC2 = 1000.0;
  else
    chi2perpointTPC2 = pair->Track2()->Track()->TPCchi2()/pair->Track2()->Track()->TPCncls();


  if (chi2perpointTPC1 > chi2perpointTPC2) {
    fChi2TPCONENumerator->Fill(tQinv, chi2perpointTPC1);
  }
  else {
    fChi2TPCONENumerator->Fill(tQinv, chi2perpointTPC2);
  }

  if (pair->Track1()->Track()->SigmaToVertex() > pair->Track2()->Track()->SigmaToVertex()) {
    fSigmaToVertexNumerator->Fill(tQinv, 
				  pair->Track1()->Track()->SigmaToVertex());
  }
  else {
    fSigmaToVertexNumerator->Fill(tQinv, 
				  pair->Track2()->Track()->SigmaToVertex());
  }
}
//____________________________
void AliFemtoChi2CorrFctn::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...

  if ((pair->Track1()->Track()->ITSncls() == 0) && (pair->Track2()->Track()->ITSncls() == 0))
    fChi2ITSSUMDenominator->Fill(tQinv, 1000.0);
  else 
    fChi2ITSSUMDenominator->Fill(tQinv, 
				 (pair->Track1()->Track()->ITSchi2() + 
				  pair->Track2()->Track()->ITSchi2())/
				 (pair->Track1()->Track()->ITSncls() +
				  pair->Track2()->Track()->ITSncls()));
  if ((pair->Track1()->Track()->TPCncls() == 0) && (pair->Track2()->Track()->TPCncls() == 0))
    fChi2TPCSUMDenominator->Fill(tQinv, 1000.0);
  else
    fChi2TPCSUMDenominator->Fill(tQinv, 
				 (pair->Track1()->Track()->TPCchi2() +
				  pair->Track2()->Track()->TPCchi2())/
				 (pair->Track1()->Track()->TPCncls() +
				  pair->Track2()->Track()->TPCncls()));
  double chi2perpointITS1, chi2perpointITS2;
  if (pair->Track1()->Track()->ITSncls() == 0)
    chi2perpointITS1 = 1000.0;
  else
    chi2perpointITS1 = pair->Track1()->Track()->ITSchi2()/pair->Track1()->Track()->ITSncls();

  if (pair->Track2()->Track()->ITSncls() == 0)
    chi2perpointITS2 = 1000.0;
  else
    chi2perpointITS2 = pair->Track2()->Track()->ITSchi2()/pair->Track2()->Track()->ITSncls();


  if (chi2perpointITS1 > chi2perpointITS2) {
    fChi2ITSONEDenominator->Fill(tQinv, chi2perpointITS1);
  }
  else {
    fChi2ITSONEDenominator->Fill(tQinv, chi2perpointITS2);
  }

  double chi2perpointTPC1, chi2perpointTPC2;
  if (pair->Track1()->Track()->TPCncls() == 0)
    chi2perpointTPC1 = 1000.0;
  else
    chi2perpointTPC1 = pair->Track1()->Track()->TPCchi2()/pair->Track1()->Track()->TPCncls();

  if (pair->Track2()->Track()->TPCncls() == 0)
    chi2perpointTPC2 = 1000.0;
  else
    chi2perpointTPC2 = pair->Track2()->Track()->TPCchi2()/pair->Track2()->Track()->TPCncls();


  if (chi2perpointTPC1 > chi2perpointTPC2) {
    fChi2TPCONEDenominator->Fill(tQinv, chi2perpointTPC1);
  }
  else {
    fChi2TPCONEDenominator->Fill(tQinv, chi2perpointTPC2);
  }
  if (pair->Track1()->Track()->SigmaToVertex() > pair->Track2()->Track()->SigmaToVertex()) {
    fSigmaToVertexDenominator->Fill(tQinv, 
				  pair->Track1()->Track()->SigmaToVertex());
  }
  else {
    fSigmaToVertexDenominator->Fill(tQinv, 
				  pair->Track2()->Track()->SigmaToVertex());
  }
}


void AliFemtoChi2CorrFctn::WriteHistos()
{
  // Write out result histograms
  fChi2ITSSUMNumerator->Write();
  fChi2ITSSUMDenominator->Write();
  fChi2TPCSUMNumerator->Write();
  fChi2TPCSUMDenominator->Write();
  fChi2ITSONENumerator->Write();
  fChi2ITSONEDenominator->Write();
  fChi2TPCONENumerator->Write();
  fChi2TPCONEDenominator->Write();
  fSigmaToVertexNumerator->Write();
  fSigmaToVertexDenominator->Write();
  
}

TList* AliFemtoChi2CorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fChi2ITSSUMNumerator);
  tOutputList->Add(fChi2ITSSUMDenominator);
  tOutputList->Add(fChi2TPCSUMNumerator);
  tOutputList->Add(fChi2TPCSUMDenominator);
  tOutputList->Add(fChi2ITSONENumerator);
  tOutputList->Add(fChi2ITSONEDenominator);
  tOutputList->Add(fChi2TPCONENumerator);
  tOutputList->Add(fChi2TPCONEDenominator);
  tOutputList->Add(fSigmaToVertexNumerator);
  tOutputList->Add(fSigmaToVertexDenominator);
  
  return tOutputList;

}
