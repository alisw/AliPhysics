/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityQAPairCut - a pair cut which checks for some pair     //
// qualities that attempt to identify slit/doubly reconstructed tracks     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityQAPairCut.cxx 24360 2008-03-10 09:48:27Z akisiel $
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

#include "AliFemtoShareQualityQAPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoShareQualityQAPairCut)
#endif

//__________________
AliFemtoShareQualityQAPairCut::AliFemtoShareQualityQAPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fShareQualityMax(1.0),
  fShareQualitymin(-0.5),
  fShareFractionMax(1.0),
  fShareFractionmin(0.0),
  fRemoveSameLabel(0),
  fShareQualityQASwitch(0),
  fShareFractionQASwitch(0)
{ 
  fShareQualityQASwitch  = false;
  fShareQualityQAExclusionZone[0] = -0.5;
  fShareQualityQAExclusionZone[1] = 1.0;
  fShareFractionQASwitch = false; 
  fShareFractionQAExclusionZone[0] = 0.0;
  fShareFractionQAExclusionZone[1] = 1.0;
}

//__________________
AliFemtoShareQualityQAPairCut::~AliFemtoShareQualityQAPairCut(){
  /* no-op */
}
//__________________
bool AliFemtoShareQualityQAPairCut::Pass(const AliFemtoPair* pair){
  // Check for pairs that are possibly shared/double reconstruction
  bool pass;
  
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
  //  if (hsmval > -0.4) {
//   cout << "Pair quality: " << hsmval << " " << an << " " << nh << " " 
//        << (pair->Track1()->Track()) << " " 
//        << (pair->Track2()->Track()) << endl;
//   cout << "Bits: " << pair->Track1()->Track()->TPCclusters().GetNbits() << endl;
    //  }
//   if (hsfval > 0.0) {
//     cout << "Pair sharity: " << hsfval << " " << ns << " " << nh << "    " << hsmval << " " << an << " " << nh << endl;
//   }

  // Determine if pair pass/fail cuts:
  if (fShareQualityQASwitch) {
  pass = ((hsmval >= fShareQualitymin)  && (hsmval < fShareQualityQAExclusionZone[0])  || 
          (hsmval >= fShareQualityQAExclusionZone[1])  && (hsmval < fShareQualityMax)) &&
         (hsfval >= fShareFractionmin) && (hsfval < fShareFractionMax);
  }
  else if (fShareFractionQASwitch) {
  pass = ((hsfval >= fShareFractionmin)  && (hsfval < fShareFractionQAExclusionZone[0])  || 
          (hsfval >= fShareFractionQAExclusionZone[1])  && (hsfval < fShareFractionMax)) &&
         (hsmval >= fShareQualitymin) && (hsmval < fShareQualityMax);
  }
  else {
  pass = (hsmval >= fShareQualitymin)  && (hsmval < fShareQualityMax) && 
         (hsfval >= fShareFractionmin) && (hsfval < fShareFractionMax);
  }

  if (fRemoveSameLabel) {
    if (abs(pair->Track1()->Track()->Label()) == abs(pair->Track2()->Track()->Label())) {
      cout << "Found a pair with same label " << pair->Track1()->Track()->Label() << endl;
      cout << "Quality Sharity Passed " << hsmval << " " << hsfval << " " << pair->QInv() << " " << pass << endl;
      pass = kFALSE;
    }
  }

  pass ? fNPairsPassed++ : fNPairsFailed++;
  return pass;
}
//__________________
AliFemtoString AliFemtoShareQualityQAPairCut::Report(){
  // Prepare the report from the execution
  string stemp = "AliFemtoShareQuality Pair Cut - remove shared and split pairs\n";  char ctemp[100];
  sprintf(ctemp,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

void AliFemtoShareQualityQAPairCut::SetShareQualityMax(Double_t aShareQualityMax) {
  fShareQualityMax = aShareQualityMax;
}

void AliFemtoShareQualityQAPairCut::SetShareQualitymin(Double_t aShareQualitymin) {
  fShareQualitymin = aShareQualitymin;
}

void AliFemtoShareQualityQAPairCut::SetShareQualityQASwitch(bool aSwitch) {
  fShareQualityQASwitch = aSwitch;
}

void AliFemtoShareQualityQAPairCut::SetShareQualityQAExclusionZone(Double_t lo, Double_t hi) {
  fShareQualityQAExclusionZone[0] = lo;
  fShareQualityQAExclusionZone[1] = hi;
}

Double_t AliFemtoShareQualityQAPairCut::GetAliFemtoShareQualityMax() const {
  return fShareQualityMax;
}

void AliFemtoShareQualityQAPairCut::SetShareFractionMax(Double_t aShareFractionMax) {
  fShareFractionMax = aShareFractionMax;
}

void AliFemtoShareQualityQAPairCut::SetShareFractionmin(Double_t aShareFractionmin) {
  fShareFractionmin = aShareFractionmin;
}

void AliFemtoShareQualityQAPairCut::SetShareFractionQASwitch(bool aSwitch) {
  fShareFractionQASwitch = aSwitch;
}

void AliFemtoShareQualityQAPairCut::SetShareFractionQAExclusionZone(Double_t lo, Double_t hi) {
  fShareFractionQAExclusionZone[0] = lo;
  fShareFractionQAExclusionZone[1] = hi;
}

Double_t AliFemtoShareQualityQAPairCut::GetAliFemtoShareFractionMax() const {
  return fShareFractionMax;
}

TList *AliFemtoShareQualityQAPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoShareQualityQAPairCut.sharequalitymax=%f", fShareQualityMax);
  snprintf(buf, 200, "AliFemtoShareQualityQAPairCut.sharefractionmax=%f", fShareFractionMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void     AliFemtoShareQualityQAPairCut::SetRemoveSameLabel(Bool_t aRemove)
{
  fRemoveSameLabel = aRemove;
}
