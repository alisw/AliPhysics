////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoShareQualityCorrFctn - A correlation function that saves the     ///
/// amount of sharing and splitting hits per pair as a function of qinv      ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoShareQualityCorrFctn.h"
//#include "Infrastructure/AliFemtoHisto.hh"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoShareQualityCorrFctn)
#endif

//____________________________
AliFemtoShareQualityCorrFctn::AliFemtoShareQualityCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  fShareNumerator(0),
  fShareDenominator(0),
  fQualityNumerator(0),
  fQualityDenominator(0)
{
  // set up numerator
  //  title = "Num Qinv (MeV/c)";
  char TitNum[100] = "NumShare";
  strcat(TitNum,title);
  fShareNumerator = new TH2D(TitNum,title,nbins,QinvLo,QinvHi,50,0.0,1.00001);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char TitDen[100] = "DenShare";
  strcat(TitDen,title);
  fShareDenominator = new TH2D(TitDen,title,nbins,QinvLo,QinvHi,50,0.0,1.00001);

  char Tit2Num[100] = "NumQuality";
  strcat(Tit2Num,title);
  fQualityNumerator = new TH2D(Tit2Num,title,nbins,QinvLo,QinvHi,75,-0.500001,1.000001);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char Tit2Den[100] = "DenQuality";
  strcat(Tit2Den,title);
  fQualityDenominator = new TH2D(Tit2Den,title,nbins,QinvLo,QinvHi,75,-0.500001,1.000001);
  // set up ratio
  //title = "Ratio Qinv (MeV/c)";
  // this next bit is unfortunately needed so that we can have many histos of same "title"
  // it is neccessary if we typedef TH2D to TH1d (which we do)
  //mShareNumerator->SetDirectory(0);
  //mShareDenominator->SetDirectory(0);
  //mRatio->SetDirectory(0);

  // to enable error bar calculation...
  fShareNumerator->Sumw2();
  fShareDenominator->Sumw2();

  fQualityNumerator->Sumw2();
  fQualityDenominator->Sumw2();
}

//____________________________
AliFemtoShareQualityCorrFctn::AliFemtoShareQualityCorrFctn(const AliFemtoShareQualityCorrFctn& aCorrFctn) :
  fShareNumerator(0),
  fShareDenominator(0),
  fQualityNumerator(0),
  fQualityDenominator(0)
{
  if (aCorrFctn.fShareNumerator)
    fShareNumerator = new TH2D(*aCorrFctn.fShareNumerator);
  if (aCorrFctn.fShareDenominator)
    fShareDenominator = new TH2D(*aCorrFctn.fShareDenominator);
  if (aCorrFctn.fQualityNumerator)
    fQualityNumerator = new TH2D(*aCorrFctn.fQualityNumerator);
  if (aCorrFctn.fQualityDenominator)
    fQualityDenominator = new TH2D(*aCorrFctn.fQualityDenominator);
}
//____________________________
AliFemtoShareQualityCorrFctn::~AliFemtoShareQualityCorrFctn(){
  delete fShareNumerator;
  delete fShareDenominator;
  delete fQualityNumerator;
  delete fQualityDenominator;
}
//_________________________
AliFemtoShareQualityCorrFctn& AliFemtoShareQualityCorrFctn::operator=(const AliFemtoShareQualityCorrFctn& aCorrFctn)
{
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fShareNumerator)
    fShareNumerator = new TH2D(*aCorrFctn.fShareNumerator);
  else
    fShareNumerator = 0;
  if (aCorrFctn.fShareDenominator)
    fShareDenominator = new TH2D(*aCorrFctn.fShareDenominator);
  else
    fShareDenominator = 0;
  if (aCorrFctn.fQualityNumerator)
    fQualityNumerator = new TH2D(*aCorrFctn.fQualityNumerator);
  else
    fQualityNumerator = 0;
  if (aCorrFctn.fQualityDenominator)
    fQualityDenominator = new TH2D(*aCorrFctn.fQualityDenominator);
  else
    fQualityDenominator = 0;

  return *this;
}
//_________________________
void AliFemtoShareQualityCorrFctn::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoShareQualityCorrFctn::Report(){
  string stemp = "Qinv Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fShareNumerator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fShareDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoShareQualityCorrFctn::AddRealPair(const AliFemtoPair* pair){
  double Qinv = fabs(pair->qInv());   // note - qInv() will be negative for identical pairs...
  Int_t nh = 0;
  Int_t an = 0;
  Int_t ns = 0;
  
  for (unsigned int imap=0; imap<pair->track1()->Track()->TPCclusters().GetNbits(); imap++) {
    // If both have clusters in the same row
    if (pair->track1()->Track()->TPCclusters().TestBitNumber(imap) && 
	pair->track2()->Track()->TPCclusters().TestBitNumber(imap)) {
      // Do they share it ?
      if (pair->track1()->Track()->TPCsharing().TestBitNumber(imap) ||
	  pair->track2()->Track()->TPCsharing().TestBitNumber(imap))
	{
	  if (Qinv < 0.01) {
	    cout << "Shared cluster in row " << imap << endl; 
	  }
	  an++;
	  nh+=2;
	  ns+=2;
	}
      
      // Different hits on the same padrow
      else {
	an--;
	nh+=2;
      }
    }
    else if (pair->track1()->Track()->TPCclusters().TestBitNumber(imap) ||
	     pair->track2()->Track()->TPCclusters().TestBitNumber(imap)) {
      // One track has a hit, the other does not
      an++;
      nh++;
    }
  }
  if (Qinv < 0.01) {
    cout << "Qinv of the pair is " << Qinv << endl;
    cout << "Clusters: " << endl;
    for (unsigned int imap=0; imap<pair->track1()->Track()->TPCclusters().GetNbits(); imap++) {
      cout << imap ;
      if (pair->track1()->Track()->TPCclusters().TestBitNumber(imap)) cout << " 1 ";
      else cout << " 0 " ;
      if (pair->track2()->Track()->TPCclusters().TestBitNumber(imap)) cout << " 1 ";
      else cout << " 0 " ;
      cout << "     ";
      if (pair->track1()->Track()->TPCsharing().TestBitNumber(imap)) cout << " S ";
      else cout << " X ";
      if (pair->track2()->Track()->TPCsharing().TestBitNumber(imap)) cout << " S ";
      else cout << " X ";
      cout << endl;
    }
  }

  Float_t hsmval = 0.0;
  Float_t hsfval = 0.0;

  if (nh >0) {
    hsmval = an*1.0/nh;
    hsfval = ns*1.0/nh;
  }

  if (Qinv < 0.01) {
    cout << "Quality  Sharity " << hsmval << " " << hsfval << " " << pair->track1()->Track() << " " << pair->track2()->Track() << endl;
  }

  fShareNumerator->Fill(Qinv, hsfval);
  fQualityNumerator->Fill(Qinv, hsmval);
  //  cout << "AliFemtoShareQualityCorrFctn::AddRealPair : " << pair->qInv() << " " << Qinv <<
  //" " << pair->track1().FourMomentum() << " " << pair->track2().FourMomentum() << endl;
}
//____________________________
void AliFemtoShareQualityCorrFctn::AddMixedPair(const AliFemtoPair* pair){
  double weight = 1.0;
  double Qinv = fabs(pair->qInv());   // note - qInv() will be negative for identical pairs...
  Int_t nh = 0;
  Int_t an = 0;
  Int_t ns = 0;
  
  for (unsigned int imap=0; imap<pair->track1()->Track()->TPCclusters().GetNbits(); imap++) {
    // If both have clusters in the same row
    if (pair->track1()->Track()->TPCclusters().TestBitNumber(imap) && 
	pair->track2()->Track()->TPCclusters().TestBitNumber(imap)) {
      // Do they share it ?
      if (pair->track1()->Track()->TPCsharing().TestBitNumber(imap) ||
	  pair->track2()->Track()->TPCsharing().TestBitNumber(imap))
	{
	  //	  cout << "A shared cluster !!!" << endl;
	  //	cout << "imap idx1 idx2 " 
	  //	     << imap << " "
	  //	     << tP1idx[imap] << " " << tP2idx[imap] << endl;
	  an++;
	  nh+=2;
	  ns+=2;
	}
      
      // Different hits on the same padrow
      else {
	an--;
	nh+=2;
      }
    }
    else if (pair->track1()->Track()->TPCclusters().TestBitNumber(imap) ||
	     pair->track2()->Track()->TPCclusters().TestBitNumber(imap)) {
      // One track has a hit, the other does not
      an++;
      nh++;
    }
  }
  
  Float_t hsmval = 0.0;
  Float_t hsfval = 0.0;

  if (nh >0) {
    hsmval = an*1.0/nh;
    hsfval = ns*1.0/nh;
  }

  fShareDenominator->Fill(Qinv,hsfval,weight);
  fQualityDenominator->Fill(Qinv,hsmval,weight);
}


void AliFemtoShareQualityCorrFctn::WriteHistos()
{
  fShareNumerator->Write();
  fShareDenominator->Write();
  fQualityNumerator->Write();
  fQualityDenominator->Write();
  
}
