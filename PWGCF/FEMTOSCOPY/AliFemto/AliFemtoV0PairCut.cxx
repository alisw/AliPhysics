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
  fRemoveSameLabel(0),
  fDataType(kAOD),
  fDTPCMin(0),
  fDTPCExitMin(0)
{
  // Default constructor
  // Nothing to do
}
//__________________
AliFemtoV0PairCut::~AliFemtoV0PairCut(){
  /* no-op */
}

AliFemtoV0PairCut& AliFemtoV0PairCut::operator=(const AliFemtoV0PairCut& cut) 
{
  if (this != &cut) {
   
    AliFemtoPairCut::operator=(cut);
    fNPairsPassed = 0;
    fNPairsFailed =0;
    fV0Max = 1.0;
    fShareFractionMax = 1.0;
    fRemoveSameLabel = 0;
    fDataType = kAOD;
    fDTPCMin = 0;
    fDTPCExitMin = 0;
  }

  return *this;
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

  bool tempTPCEntrancePos = true;
  bool tempTPCEntranceNeg = true;
  bool tempTPCExitPos = true;
  bool tempTPCExitNeg = true;
  if(fDataType==kESD || fDataType==kAOD)
    {
      double distx = pair->Track1()->V0()->NominalTpcEntrancePointPos().x() - pair->Track2()->V0()->NominalTpcEntrancePointPos().x();
      double disty = pair->Track1()->V0()->NominalTpcEntrancePointPos().y() - pair->Track2()->V0()->NominalTpcEntrancePointPos().y();
      double distz = pair->Track1()->V0()->NominalTpcEntrancePointPos().z() - pair->Track2()->V0()->NominalTpcEntrancePointPos().z();
      double distPos = sqrt(distx*distx + disty*disty + distz*distz);

      distx = pair->Track1()->V0()->NominalTpcEntrancePointNeg().x() - pair->Track2()->V0()->NominalTpcEntrancePointNeg().x();
      disty = pair->Track1()->V0()->NominalTpcEntrancePointNeg().y() - pair->Track2()->V0()->NominalTpcEntrancePointNeg().y();
      distz = pair->Track1()->V0()->NominalTpcEntrancePointNeg().z() - pair->Track2()->V0()->NominalTpcEntrancePointNeg().z();
      double distNeg = sqrt(distx*distx + disty*disty + distz*distz);

      double distExitX = pair->Track1()->V0()->NominalTpcExitPointPos().x() - pair->Track2()->V0()->NominalTpcExitPointPos().x();
      double distExitY = pair->Track1()->V0()->NominalTpcExitPointPos().y() - pair->Track2()->V0()->NominalTpcExitPointPos().y();
      double distExitZ = pair->Track1()->V0()->NominalTpcExitPointPos().z() - pair->Track2()->V0()->NominalTpcExitPointPos().z();
      double distExitPos = sqrt(distExitX*distExitX + distExitY*distExitY + distExitZ*distExitZ);

      distExitX = pair->Track1()->V0()->NominalTpcExitPointNeg().x() - pair->Track2()->V0()->NominalTpcExitPointNeg().x();
      distExitY = pair->Track1()->V0()->NominalTpcExitPointNeg().y() - pair->Track2()->V0()->NominalTpcExitPointNeg().y();
      distExitZ = pair->Track1()->V0()->NominalTpcExitPointNeg().z() - pair->Track2()->V0()->NominalTpcExitPointNeg().z();
      double distExitNeg = sqrt(distExitX*distExitX + distExitY*distExitY + distExitZ*distExitZ);

      tempTPCEntrancePos = distPos > fDTPCMin;
      tempTPCEntranceNeg = distNeg > fDTPCMin;

      tempTPCExitPos = distExitPos > fDTPCExitMin;
      tempTPCExitNeg = distExitNeg > fDTPCExitMin;
    }
 
  if(!tempTPCEntrancePos || !tempTPCEntranceNeg || !tempTPCExitPos || !tempTPCExitNeg) return false;


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

void AliFemtoV0PairCut::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}

void AliFemtoV0PairCut::SetTPCEntranceSepMinimum(double dtpc)
{
  fDTPCMin = dtpc;
}

void AliFemtoV0PairCut::SetTPCExitSepMinimum(double dtpc)
{
  fDTPCExitMin = dtpc;
}
