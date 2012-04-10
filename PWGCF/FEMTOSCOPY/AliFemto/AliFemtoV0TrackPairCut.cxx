/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityPairCut - a pair cut which checks for some pair     //
// qualities that attempt to identify slit/doubly reconstructed tracks     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityPairCut.cxx 50722 2011-07-21 15:18:38Z akisiel $
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

#include "AliFemtoV0TrackPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoV0TrackPairCut)
#endif

//__________________
AliFemtoV0TrackPairCut::AliFemtoV0TrackPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fV0Max(1.0),
  fShareQualityMax(1.0),
  fShareFractionMax(1.0),
  fRemoveSameLabel(0),
  fTrackTPCOnly(0)
{
  // Default constructor
  // Nothing to do
}
//__________________
AliFemtoV0TrackPairCut::~AliFemtoV0TrackPairCut(){
  /* no-op */
}
//__________________
AliFemtoV0TrackPairCut& AliFemtoV0TrackPairCut::operator=(const AliFemtoV0TrackPairCut& cut) 
{
  if (this != &cut) {
   
    AliFemtoPairCut::operator=(cut);
    fNPairsPassed = 0;
    fNPairsFailed =0;
    fV0Max = 1.0;
    fShareQualityMax = 1.0;
    fShareFractionMax = 1.0;
    fRemoveSameLabel = 0;
    fTrackTPCOnly =0;
  }

  return *this;
}

//__________________
bool AliFemtoV0TrackPairCut::Pass(const AliFemtoPair* pair){
  // Check for pairs that are possibly shared/double reconstruction

  bool temp = true;
  //Track1 - V0
  //Track2 - track

  //!!!!!USUN MNIE!!!!
  /*if(pair->KStar()>0.15)
    {
      return false;
      }*/
  //!!!!!!!!


  if(!(pair->Track1()->V0() && pair->Track2()->Track()))
    {
      return false;
    }
  if(fTrackTPCOnly)
    {
      if( (-(pair->Track1()->V0()->IdNeg()+1)) ==pair->Track2()->TrackId() || (-(pair->Track1()->V0()->IdPos()+1)) ==pair->Track2()->TrackId())
	{
	  return false;
	}
    }
  else
    {
      if(pair->Track1()->V0()->IdNeg()==pair->Track2()->TrackId() || pair->Track1()->V0()->IdPos()==pair->Track2()->TrackId())
	{
	  return false;
	}
    }
  


  //reject merged trakcs in TPC
  

  //temp = dist > fDTPCMin;
  //koniec kopii
  //if(!temp)
  //return false;

  //kopia z AliFemtoShareQualityPairCut.cxx
  Int_t nh = 0;
  Int_t an = 0;
  Int_t ns = 0;

 
  
  if ((fShareFractionMax < 1.0) && ( fShareQualityMax < 1.0)) {
    for (unsigned int imap=0; imap<pair->Track1()->V0()->TPCclustersPos().GetNbits(); imap++) {
      // If both have clusters in the same row
      if (pair->Track1()->V0()->TPCclustersPos().TestBitNumber(imap) && 
	  pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
	// Do they share it ?
	if (pair->Track1()->V0()->TPCsharingPos().TestBitNumber(imap) &&
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
      else if (pair->Track1()->V0()->TPCclustersPos().TestBitNumber(imap) ||
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
    
    temp = (hsmval < fShareQualityMax) && (hsfval < fShareFractionMax);
    if(!temp) return false;

    nh = 0;
    an = 0;
    ns = 0;

    for (unsigned int imap=0; imap<pair->Track1()->V0()->TPCclustersNeg().GetNbits(); imap++) {
      // If both have clusters in the same row
      if (pair->Track1()->V0()->TPCclustersNeg().TestBitNumber(imap) && 
	  pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
	// Do they share it ?
	if (pair->Track1()->V0()->TPCsharingNeg().TestBitNumber(imap) &&
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
      else if (pair->Track1()->V0()->TPCclustersNeg().TestBitNumber(imap) ||
	       pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
	// One track has a hit, the other does not
	an++;
	nh++;
      }
    }
    
    hsmval = 0.0;
    hsfval = 0.0;
    
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
    
    temp = (hsmval < fShareQualityMax) && (hsfval < fShareFractionMax);


  }
  else
    temp = true;
  //koniec kopii





  return temp;
}
//__________________
AliFemtoString AliFemtoV0TrackPairCut::Report(){
  // Prepare the report from the execution
  string stemp = "AliFemtoV0 Pair Cut - remove shared and split pairs\n";  char ctemp[100];
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

void AliFemtoV0TrackPairCut::SetV0Max(Double_t aV0Max) {
  fV0Max = aV0Max;
}

Double_t AliFemtoV0TrackPairCut::GetAliFemtoV0Max() const {
  return fV0Max;
}


TList *AliFemtoV0TrackPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoV0TrackPairCut.sharequalitymax=%f", fV0Max);
  snprintf(buf, 200, "AliFemtoV0TrackPairCut.sharefractionmax=%f", fShareFractionMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void     AliFemtoV0TrackPairCut::SetRemoveSameLabel(Bool_t aRemove)
{
  fRemoveSameLabel = aRemove;
}

void AliFemtoV0TrackPairCut::SetTPCOnly(Bool_t tpconly)
{
  fTrackTPCOnly = tpconly;
}

void AliFemtoV0TrackPairCut::SetShareQualityMax(Double_t aShareQualityMax) {
  fShareQualityMax = aShareQualityMax;
}

void AliFemtoV0TrackPairCut::SetShareFractionMax(Double_t aShareFractionMax) {
  fShareFractionMax = aShareFractionMax;
}
