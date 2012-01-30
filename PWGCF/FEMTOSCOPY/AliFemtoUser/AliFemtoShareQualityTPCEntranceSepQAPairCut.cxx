/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityTPCEntranceSepQAPairCut - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityTPCEntranceSepQAPairCut.cxx,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
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

#include "AliFemtoShareQualityTPCEntranceSepQAPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoShareQualityTPCEntranceSepQAPairCut)
#endif

//__________________
AliFemtoShareQualityTPCEntranceSepQAPairCut::AliFemtoShareQualityTPCEntranceSepQAPairCut():
  AliFemtoShareQualityQAPairCut(),
  fDTPCMin(0.0),
  fDTPCMax(1000.0),
  fDTPCQASwitch(0)
{
  fDTPCQASwitch = false;
  fDTPCQAExclusionZone[0] = 0.0;
  fDTPCQAExclusionZone[1] = 1000.0;
}
//__________________
AliFemtoShareQualityTPCEntranceSepQAPairCut::AliFemtoShareQualityTPCEntranceSepQAPairCut(const AliFemtoShareQualityTPCEntranceSepQAPairCut& c) : 
  AliFemtoShareQualityQAPairCut(c),
  fDTPCMin(0),
  fDTPCMax(1000.0),
  fDTPCQASwitch(0)
{ 
  fDTPCMin = c.fDTPCMin;
  fDTPCMax = c.fDTPCMax;
  fDTPCQASwitch = c.fDTPCQASwitch;
  fDTPCQAExclusionZone[0] = c.fDTPCQAExclusionZone[0];
  fDTPCQAExclusionZone[1] = c.fDTPCQAExclusionZone[1];
}

//__________________
AliFemtoShareQualityTPCEntranceSepQAPairCut& AliFemtoShareQualityTPCEntranceSepQAPairCut::operator=(const AliFemtoShareQualityTPCEntranceSepQAPairCut& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  fDTPCMin = aCut.fDTPCMin;
  fDTPCMax = aCut.fDTPCMax;
  fDTPCQASwitch = aCut.fDTPCQASwitch;
  fDTPCQAExclusionZone[0] = aCut.fDTPCQAExclusionZone[0];
  fDTPCQAExclusionZone[1] = aCut.fDTPCQAExclusionZone[1];

  return *this;
}

//__________________
AliFemtoShareQualityTPCEntranceSepQAPairCut::~AliFemtoShareQualityTPCEntranceSepQAPairCut(){
  /* no-op */
}
//__________________
bool AliFemtoShareQualityTPCEntranceSepQAPairCut::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC entrance separation and
  // quality and sharity
  bool pass = true;
  
  double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
  double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
  double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
  double dist = sqrt(distx*distx + disty*disty + distz*distz);

  if (fDTPCQASwitch) {
    pass = ((dist > fDTPCMin) && (dist < fDTPCQAExclusionZone[0])) ||
           ((dist > fDTPCQAExclusionZone[1]) && (dist < fDTPCMax));
  }
  else {
    pass = (dist > fDTPCMin) && (dist < fDTPCMax);
  }

  if (pass) {
    pass = AliFemtoShareQualityQAPairCut::Pass(pair);
  }
  else
    fNPairsFailed++;

  return pass;
}
//__________________
AliFemtoString AliFemtoShareQualityTPCEntranceSepQAPairCut::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoShareQualityTPCEntranceSep Pair Cut - remove shared and split pairs and pairs with small separation at the entrance to the TPC\n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with TPC entrance separation more that %f",fDTPCMin);
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoShareQualityTPCEntranceSepQAPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoShareQualityQAPairCut::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoShareQualityTPCEntranceSepQAPairCut.tpcentsepmin=%f", fDTPCMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoShareQualityTPCEntranceSepQAPairCut::SetTPCEntranceSepMinimum(double dtpc)
{
  fDTPCMin = dtpc;
}

void AliFemtoShareQualityTPCEntranceSepQAPairCut::SetTPCEntranceSepMaximum(double dtpc)
{
  fDTPCMax = dtpc;
}

void AliFemtoShareQualityTPCEntranceSepQAPairCut::SetTPCEntranceSepQASwitch(bool Switch)
{
  fDTPCQASwitch = Switch;
}

void AliFemtoShareQualityTPCEntranceSepQAPairCut::SetTPCEntranceSepQAExclusionZone(double lo, double hi)
{
  fDTPCQAExclusionZone[0] = lo;
  fDTPCQAExclusionZone[1] = hi;
}

