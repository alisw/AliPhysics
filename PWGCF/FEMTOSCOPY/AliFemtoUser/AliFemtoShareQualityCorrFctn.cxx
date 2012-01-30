////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoShareQualityCorrFctn - A correlation function that saves the     ///
/// amount of sharing and splitting hits per pair as a function of qinv      ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoShareQualityCorrFctn.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoShareQualityCorrFctn)
#endif

//____________________________
AliFemtoShareQualityCorrFctn::AliFemtoShareQualityCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  AliFemtoCorrFctn(),
  fShareNumerator(0),
  fShareDenominator(0),
  fQualityNumerator(0),
  fQualityDenominator(0),
  fTPCSepNumerator(0),
  fTPCSepDenominator(0)
{
  // set up numerator
  //  title = "Num Qinv (MeV/c)";
  char tTitNum[101] = "NumShare";
  strncat(tTitNum,title, 100);
  fShareNumerator = new TH2D(tTitNum,title,nbins,QinvLo,QinvHi,100,0.0,1.00001);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTitDen[101] = "DenShare";
  strncat(tTitDen,title, 100);
  fShareDenominator = new TH2D(tTitDen,title,nbins,QinvLo,QinvHi,100,0.0,1.00001);

  char tTit2Num[101] = "NumQuality";
  strncat(tTit2Num,title, 100);
  fQualityNumerator = new TH2D(tTit2Num,title,nbins,QinvLo,QinvHi,150,-0.500001,1.000001);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTit2Den[101] = "DenQuality";
  strncat(tTit2Den,title, 100);
  fQualityDenominator = new TH2D(tTit2Den,title,nbins,QinvLo,QinvHi,150,-0.500001,1.000001);
  // set up ratio
  //title = "Ratio Qinv (MeV/c)";
  // this next bit is unfortunately needed so that we can have many histos of same "title"
  // it is neccessary if we typedef TH2D to TH1d (which we do)
  //mShareNumerator->SetDirectory(0);
  //mShareDenominator->SetDirectory(0);
  //mRatio->SetDirectory(0);

  char tTit3Num[101] = "NumTPCSep";
  strncat(tTit3Num,title, 100);
  fTPCSepNumerator = new TH2D(tTit3Num,title,nbins,QinvLo,QinvHi,150,0.0,100.0);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTit3Den[101] = "DenTPCSep";
  strncat(tTit3Den,title, 100);
  fTPCSepDenominator = new TH2D(tTit3Den,title,nbins,QinvLo,QinvHi,150,0.0,100.0);

  // to enable error bar calculation...
  fShareNumerator->Sumw2();
  fShareDenominator->Sumw2();

  fQualityNumerator->Sumw2();
  fQualityDenominator->Sumw2();
  
  fTPCSepNumerator->Sumw2();
  fTPCSepDenominator->Sumw2();
}

//____________________________
AliFemtoShareQualityCorrFctn::AliFemtoShareQualityCorrFctn(const AliFemtoShareQualityCorrFctn& aCorrFctn) :
  AliFemtoCorrFctn(),
  fShareNumerator(0),
  fShareDenominator(0),
  fQualityNumerator(0),
  fQualityDenominator(0),
  fTPCSepNumerator(0),
  fTPCSepDenominator(0)
{
  // copy constructor
  if (aCorrFctn.fShareNumerator)
    fShareNumerator = new TH2D(*aCorrFctn.fShareNumerator);
  if (aCorrFctn.fShareDenominator)
    fShareDenominator = new TH2D(*aCorrFctn.fShareDenominator);
  if (aCorrFctn.fQualityNumerator)
    fQualityNumerator = new TH2D(*aCorrFctn.fQualityNumerator);
  if (aCorrFctn.fQualityDenominator)
    fQualityDenominator = new TH2D(*aCorrFctn.fQualityDenominator);
  if (aCorrFctn.fTPCSepNumerator)
    fTPCSepNumerator = new TH2D(*aCorrFctn.fTPCSepNumerator);
  if (aCorrFctn.fTPCSepDenominator)
    fTPCSepDenominator = new TH2D(*aCorrFctn.fTPCSepDenominator);
}
//____________________________
AliFemtoShareQualityCorrFctn::~AliFemtoShareQualityCorrFctn(){
  // destructor
  delete fShareNumerator;
  delete fShareDenominator;
  delete fQualityNumerator;
  delete fQualityDenominator;
  delete fTPCSepNumerator;
  delete fTPCSepDenominator;
}
//_________________________
AliFemtoShareQualityCorrFctn& AliFemtoShareQualityCorrFctn::operator=(const AliFemtoShareQualityCorrFctn& aCorrFctn)
{
  // assignment operator
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
  if (aCorrFctn.fTPCSepNumerator)
    fTPCSepNumerator = new TH2D(*aCorrFctn.fTPCSepNumerator);
  else
    fTPCSepNumerator = 0;
  if (aCorrFctn.fTPCSepDenominator)
    fTPCSepDenominator = new TH2D(*aCorrFctn.fTPCSepDenominator);
  else
    fTPCSepDenominator = 0;

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
  // create report
  string stemp = "Qinv Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fShareNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fShareDenominator->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoShareQualityCorrFctn::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  Int_t nh = 0;
  Int_t an = 0;
  Int_t ns = 0;
  
  for (unsigned int imap=0; imap<pair->Track1()->Track()->TPCclusters().GetNbits(); imap++) {
    // If both have clusters in the same row
    if (pair->Track1()->Track()->TPCclusters().TestBitNumber(imap) && 
	pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
      // Do they share it ?
      if (pair->Track1()->Track()->TPCsharing().TestBitNumber(imap) &&
	  pair->Track2()->Track()->TPCsharing().TestBitNumber(imap))
	{
// 	  if (tQinv < 0.01) {
// 	    cout << "Shared cluster in row " << imap << endl; 
// 	  }
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
    else if (pair->Track1()->Track()->TPCclusters().TestBitNumber(imap) ||
	     pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
      // One track has a hit, the other does not
      an++;
      nh++;
    }
  }
//    if (tQinv < 0.01) {
//     cout << "Qinv of the pair is " << tQinv << endl;
//     cout << "Clusters: " << endl;
//     for (unsigned int imap=0; imap<pair->Track1()->Track()->TPCclusters().GetNbits(); imap++) {
//       cout << imap ;
//       if (pair->Track1()->Track()->TPCclusters().TestBitNumber(imap)) cout << " 1 ";
//       else cout << " 0 " ;
//       if (pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) cout << " 1 ";
//       else cout << " 0 " ;
//       cout << "     ";
//       if (pair->Track1()->Track()->TPCsharing().TestBitNumber(imap)) cout << " S ";
//       else cout << " X ";
//       if (pair->Track2()->Track()->TPCsharing().TestBitNumber(imap)) cout << " S ";
//       else cout << " X ";
//       cout << endl;
//     }
//   }

  Float_t hsmval = 0.0;
  Float_t hsfval = 0.0;

  if (nh >0) {
    hsmval = an*1.0/nh;
    hsfval = ns*1.0/nh;
  }

//   if ((tQinv < 0.005) && (hsmval<-0.0)) {
//     cout << "Quality  Sharity " << hsmval << " " << hsfval << " " << pair->Track1()->Track() << " " << pair->Track2()->Track() << endl;
//     cout << "Qinv of the pair is " << tQinv << endl;
//     cout << "Clusters: " << endl;
//     for (unsigned int imap=0; imap<pair->Track1()->Track()->TPCclusters().GetNbits(); imap++) {
//       cout << imap ;
//       if (pair->Track1()->Track()->TPCclusters().TestBitNumber(imap)) cout << " 1 ";
//       else cout << " 0 " ;
//       if (pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) cout << " 1 ";
//       else cout << " 0 " ;
//       cout << "     ";
//       if (pair->Track1()->Track()->TPCsharing().TestBitNumber(imap)) cout << " S ";
//       else cout << " X ";
//       if (pair->Track2()->Track()->TPCsharing().TestBitNumber(imap)) cout << " S ";
//       else cout << " X ";
//       cout << endl;
//     }
//     cout << "Momentum1 " 
// 	 << pair->Track1()->Track()->P().x() << " " 
// 	 << pair->Track1()->Track()->P().y() << " "  
// 	 << pair->Track1()->Track()->P().z() << " "  
// 	 << pair->Track1()->Track()->Label() << " "  
// 	 << pair->Track1()->Track()->TrackId() << " "  
// 	 << pair->Track1()->Track()->Flags() << " "
// 	 << pair->Track1()->Track()->KinkIndex(0) << " "
// 	 << pair->Track1()->Track()->KinkIndex(1) << " "
// 	 << pair->Track1()->Track()->KinkIndex(2) << " "
// 	 << pair->Track1()->Track()->ITSchi2() << " "
// 	 << pair->Track1()->Track()->ITSncls() << " "
// 	 << pair->Track1()->Track()->TPCchi2() << " "
// 	 << pair->Track1()->Track()->TPCncls() << " "
// 	 << endl;
//     cout << "Momentum2 " 
// 	 << pair->Track2()->Track()->P().x() << " "  
// 	 << pair->Track2()->Track()->P().y() << " "  
// 	 << pair->Track2()->Track()->P().z() << " "  
// 	 << pair->Track2()->Track()->Label() << " "  
// 	 << pair->Track2()->Track()->TrackId() << " "  
// 	 << pair->Track2()->Track()->Flags() << " " 
// 	 << pair->Track2()->Track()->KinkIndex(0) << " "
// 	 << pair->Track2()->Track()->KinkIndex(1) << " "
// 	 << pair->Track2()->Track()->KinkIndex(2) << " "
// 	 << pair->Track2()->Track()->ITSchi2() << " "
// 	 << pair->Track2()->Track()->ITSncls() << " "
// 	 << pair->Track2()->Track()->TPCchi2() << " "
// 	 << pair->Track2()->Track()->TPCncls() << " "
// 	 << endl;
//   }

  double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
  double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
  double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
  double dist = sqrt(distx*distx + disty*disty + distz*distz);

  fShareNumerator->Fill(tQinv, hsfval);
  fQualityNumerator->Fill(tQinv, hsmval);
  fTPCSepNumerator->Fill(tQinv, dist);

  //  cout << "AliFemtoShareQualityCorrFctn::AddRealPair : " << pair->qInv() << " " << tQinv <<
  //" " << pair->Track1().FourMomentum() << " " << pair->Track2().FourMomentum() << endl;
}
//____________________________
void AliFemtoShareQualityCorrFctn::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double weight = 1.0;
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  Int_t nh = 0;
  Int_t an = 0;
  Int_t ns = 0;
  
  for (unsigned int imap=0; imap<pair->Track1()->Track()->TPCclusters().GetNbits(); imap++) {
    // If both have clusters in the same row
    if (pair->Track1()->Track()->TPCclusters().TestBitNumber(imap) && 
	pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
      // Do they share it ?
      if (pair->Track1()->Track()->TPCsharing().TestBitNumber(imap) &&
	  pair->Track2()->Track()->TPCsharing().TestBitNumber(imap))
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
    else if (pair->Track1()->Track()->TPCclusters().TestBitNumber(imap) ||
	     pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
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

  double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
  double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
  double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
  double dist = sqrt(distx*distx + disty*disty + distz*distz);

  fShareDenominator->Fill(tQinv,hsfval,weight);
  fQualityDenominator->Fill(tQinv,hsmval,weight);
  fTPCSepDenominator->Fill(tQinv, dist, weight);

}


void AliFemtoShareQualityCorrFctn::WriteHistos()
{
  // Write out result histograms
  fShareNumerator->Write();
  fShareDenominator->Write();
  fQualityNumerator->Write();
  fQualityDenominator->Write();
  fTPCSepNumerator->Write();
  fTPCSepDenominator->Write();
  
}
//______________________________
TList* AliFemtoShareQualityCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fShareNumerator); 
  tOutputList->Add(fShareDenominator);  
  tOutputList->Add(fQualityNumerator);  
  tOutputList->Add(fQualityDenominator);  
  tOutputList->Add(fTPCSepNumerator);  
  tOutputList->Add(fTPCSepDenominator);  

  return tOutputList;
}
