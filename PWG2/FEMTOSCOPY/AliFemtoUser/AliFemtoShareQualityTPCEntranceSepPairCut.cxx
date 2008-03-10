/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityTPCEntranceSepPairCut - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityTPCEntranceSepPairCut.cxx,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
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

#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoShareQualityTPCEntranceSepPairCut)
#endif

//__________________
AliFemtoShareQualityTPCEntranceSepPairCut::AliFemtoShareQualityTPCEntranceSepPairCut():
  AliFemtoShareQualityPairCut(),
  fDTPCMin(0)
{
}
//__________________
AliFemtoShareQualityTPCEntranceSepPairCut::AliFemtoShareQualityTPCEntranceSepPairCut(const AliFemtoShareQualityTPCEntranceSepPairCut& c) : 
  AliFemtoShareQualityPairCut(c),
  fDTPCMin(0)
{ 
  fDTPCMin = c.fDTPCMin;
}

//__________________
AliFemtoShareQualityTPCEntranceSepPairCut::~AliFemtoShareQualityTPCEntranceSepPairCut(){
  /* no-op */
}
//__________________
bool AliFemtoShareQualityTPCEntranceSepPairCut::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC entrance separation and
  // quality and sharity
  bool temp = true;
  
  double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
  double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
  double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
  double dist = sqrt(distx*distx + disty*disty + distz*distz);

  temp = dist > fDTPCMin;

  if (temp) {
    temp = AliFemtoShareQualityPairCut::Pass(pair);
  }
  else
    fNPairsFailed++;

  return temp;
}
//__________________
AliFemtoString AliFemtoShareQualityTPCEntranceSepPairCut::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoShareQualityTPCEntranceSep Pair Cut - remove shared and split pairs and pairs with small separation at the entrance to the TPC\n";  char ctemp[100];
  sprintf(ctemp,"Accept pair with TPC entrance separation more that %f",fDTPCMin);
  stemp += ctemp;
  sprintf(ctemp,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoShareQualityTPCEntranceSepPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoShareQualityPairCut::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoShareQualityTPCEntranceSepPairCut.tpcentsepmin=%f", fDTPCMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoShareQualityTPCEntranceSepPairCut::SetTPCEntranceSepMinimum(double dtpc)
{
  fDTPCMin = dtpc;
}
