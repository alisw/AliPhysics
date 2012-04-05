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

#include "AliFemtoV0PairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoV0PairCut)
#endif

//__________________
AliFemtoV0PairCut::AliFemtoV0PairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fV0Max(1.0),
  fShareFractionMax(1.0),
  fRemoveSameLabel(0)
{
  // Default constructor
  // Nothing to do
}
//__________________
AliFemtoV0PairCut::~AliFemtoV0PairCut(){
  /* no-op */
}
//__________________
bool AliFemtoV0PairCut::Pass(const AliFemtoPair* pair){
  // Check for pairs that are possibly shared/double reconstruction

  bool temp = true;

  /*cout<<"pair->Track1(): "<<pair->Track1()<<endl;
  cout<<"pair->Track2(): "<<pair->Track2()<<endl;
  cout<<"pair->Track1()->V0(): "<<pair->Track1()->V0()<<endl;
  cout<<"pair->Track2()->V0(): "<<pair->Track2()->V0()<<endl;
  cout<<"pair->Track1()->V0()->IdNeg(): "<<pair->Track1()->V0()->IdNeg()<<endl;
  cout<<"pair->Track2()->V0()->IdNeg(): "<<pair->Track2()->V0()->IdNeg()<<endl;
  cout<<"pair->Track1()->V0()->IdPos(): "<<pair->Track1()->V0()->IdPos()<<endl;
  cout<<"pair->Track2()->V0()->IdPos(): "<<pair->Track2()->V0()->IdPos()<<endl;*/

  if(!(pair->Track1()->V0() && pair->Track2()->V0()))
    {
      return false;
    }
  if(pair->Track1()->V0()->IdNeg()==pair->Track2()->V0()->IdNeg() || pair->Track1()->V0()->IdPos()==pair->Track2()->V0()->IdPos())
    {

      return false;
    }
  

  return temp;
}
//__________________
AliFemtoString AliFemtoV0PairCut::Report(){
  // Prepare the report from the execution
  string stemp = "AliFemtoV0 Pair Cut - remove shared and split pairs\n";  char ctemp[100];
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

void AliFemtoV0PairCut::SetV0Max(Double_t aV0Max) {
  fV0Max = aV0Max;
}

Double_t AliFemtoV0PairCut::GetAliFemtoV0Max() const {
  return fV0Max;
}


TList *AliFemtoV0PairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoV0PairCut.sharequalitymax=%f", fV0Max);
  snprintf(buf, 200, "AliFemtoV0PairCut.sharefractionmax=%f", fShareFractionMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void     AliFemtoV0PairCut::SetRemoveSameLabel(Bool_t aRemove)
{
  fRemoveSameLabel = aRemove;
}
