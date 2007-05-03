/***************************************************************************
 *
 * $Id$
 *
 * Author: Adam Kisiel, Ohio State, kisiel@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a cut to remove "shared" and "split" pairs
 *
 ***************************************************************************
 *
 *
 **************************************************************************/

#include "AliFemtoShareQualityPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoShareQualityPairCut)
#endif

//__________________
AliFemtoShareQualityPairCut::AliFemtoShareQualityPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fShareQualityMax(1.0)
{
}
//__________________
AliFemtoShareQualityPairCut::~AliFemtoShareQualityPairCut(){
  /* no-op */
}
//__________________
bool AliFemtoShareQualityPairCut::Pass(const AliFemtoPair* pair){
  bool temp;
  
  Int_t nh = 0;
  Int_t an = 0;
  Int_t ns = 0;
  
  for (unsigned int imap=0; imap<pair->track1()->Track()->TPCclusters().GetNbits(); imap++) {
    // If both have clusters in the same row
    if (pair->track1()->Track()->TPCclusters().TestBitNumber(imap) && 
	pair->track2()->Track()->TPCclusters().TestBitNumber(imap)) {
      // Do they share it ?
      if (pair->track1()->Track()->TPCsharing().TestBitNumber(imap) &&
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
  //  if (hsmval > -0.4) {
  cout << "Pair quality: " << hsmval << " " << an << " " << nh << " " 
       << (pair->track1()->Track()) << " " 
       << (pair->track2()->Track()) << endl;
  cout << "Bits: " << pair->track1()->Track()->TPCclusters().GetNbits() << endl;
    //  }
  if (hsfval > 0.0) {
    cout << "Pair sharity: " << hsfval << " " << ns << " " << nh << "    " << hsmval << " " << an << " " << nh << endl;
  }

  temp = hsmval < fShareQualityMax;

  temp ? fNPairsPassed++ : fNPairsFailed++;
  return temp;
}
//__________________
AliFemtoString AliFemtoShareQualityPairCut::Report(){
  string Stemp = "AliFemtoShareQuality Pair Cut - remove shared and split pairs\n";  char Ctemp[100];
  sprintf(Ctemp,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  Stemp += Ctemp;
  AliFemtoString returnThis = Stemp;
  return returnThis;}
//__________________

void AliFemtoShareQualityPairCut::SetShareQualityMax(Double_t aShareQualityMax) {
  fShareQualityMax = aShareQualityMax;
}

Double_t AliFemtoShareQualityPairCut::GetAliFemtoShareQualityMax() {
  return fShareQualityMax;
}
