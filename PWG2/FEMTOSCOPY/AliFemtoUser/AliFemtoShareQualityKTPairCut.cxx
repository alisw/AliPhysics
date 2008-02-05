/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityKTPairCut - a pair cut which checks for some pair   //
// qualities that attempt to identify slit/doubly reconstructed tracks     //
// and selects pairs based on their transverse momentum kT                 //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityKTPairCut.cxx,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
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

#include "AliFemtoShareQualityKTPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoShareQualityKTPairCut)
#endif

//__________________
AliFemtoShareQualityKTPairCut::AliFemtoShareQualityKTPairCut():
  AliFemtoShareQualityPairCut(),
  fKTMin(0),
  fKTMax(1.0e6)
{
}
//__________________
AliFemtoShareQualityKTPairCut::AliFemtoShareQualityKTPairCut(const AliFemtoShareQualityKTPairCut& c) : 
  AliFemtoShareQualityPairCut(c),
  fKTMin(0),
  fKTMax(1.0e6)
{ 
  fKTMin = c.fKTMin;
  fKTMax = c.fKTMax;
}

//__________________
AliFemtoShareQualityKTPairCut::~AliFemtoShareQualityKTPairCut(){
  /* no-op */
}
//__________________
bool AliFemtoShareQualityKTPairCut::Pass(const AliFemtoPair* pair){
  // Accept a pair base on its Kt and sharity and quality
  bool temp = true;
  
  if (pair->KT() < fKTMin)
    temp = false;

  if (pair->KT() > fKTMax)
    temp = false;

  if (temp) {
    temp = AliFemtoShareQualityPairCut::Pass(pair);
  }
  else
    fNPairsFailed++;

  return temp;
}
//__________________
AliFemtoString AliFemtoShareQualityKTPairCut::Report(){
  // Prepare a report from execution
  string stemp = "AliFemtoShareQuality Pair Cut - remove shared and split pairs\n";  char ctemp[100];
  sprintf(ctemp,"Accept pair with kT in range %lf , %lf",fKTMin,fKTMax);
  stemp += ctemp;
  sprintf(ctemp,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoShareQualityKTPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoShareQualityPairCut::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoShareQualityKTPairCut.ktmax=%lf", fKTMax);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoShareQualityKTPairCut.ktmin=%lf", fKTMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoShareQualityKTPairCut::SetKTRange(double ktmin, double ktmax)
{
  // Set the accepted kT range
  fKTMin = ktmin;
  fKTMax = ktmax;
}
