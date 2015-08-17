/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoDeltaPtPairCut - a pair cut which selects pairs based on their       //
// transverse momentum kT                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoDeltaPtPairCut.cxx,v 1.1.2.2 2007/11/09 11:20:35 akisiel Exp $
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
#include "AliFemtoDeltaPtPairCut.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoDeltaPtPairCut)
#endif

//__________________
AliFemtoDeltaPtPairCut::AliFemtoDeltaPtPairCut():
  AliFemtoPairCut(),
  fDeltaPtMin(0),
  fDeltaPtMax(0)
{

}
//__________________
AliFemtoDeltaPtPairCut::AliFemtoDeltaPtPairCut(double lo, double hi) :
  AliFemtoPairCut(),
  fDeltaPtMin(lo),
  fDeltaPtMax(hi)
{
}
//__________________
AliFemtoDeltaPtPairCut::AliFemtoDeltaPtPairCut(const AliFemtoDeltaPtPairCut& c) : 
  AliFemtoPairCut(c),
  fDeltaPtMin(0),
  fDeltaPtMax(0)
{ 
  fDeltaPtMin = c.fDeltaPtMin;
  fDeltaPtMax = c.fDeltaPtMax;
}

//__________________
AliFemtoDeltaPtPairCut::~AliFemtoDeltaPtPairCut(){
  /* no-op */
}
AliFemtoDeltaPtPairCut& AliFemtoDeltaPtPairCut::operator=(const AliFemtoDeltaPtPairCut& c)
{
  if (this != &c) {
    fDeltaPtMin = c.fDeltaPtMin;
    fDeltaPtMax = c.fDeltaPtMax;
  }

  return *this;
}
//__________________
/*bool AliFemtoDeltaPtPairCut::Pass(const AliFemtoPair* pair){
  bool temp = true;
  
  if (pair->KT() < fKTMin)
    temp = false;

  if (pair->KT() > fKTMax)
    temp = false;

  return temp;
}*/
//__________________
AliFemtoString AliFemtoDeltaPtPairCut::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoKT Pair Cut \n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with DeltaPt in range %f , %f",fDeltaPtMin,fDeltaPtMax);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoDeltaPtPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoDeltaPtPairCut.ktmax=%f", fDeltaPtMax);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoDeltaPtPairCut.ktmin=%f", fDeltaPtMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoDeltaPtPairCut::SetDeltaPtRange(double ptmin, double ptmax)
{
  fDeltaPtMin = ptmin;
  fDeltaPtMax = ptmax;
}


//______________________________________________________
bool AliFemtoDeltaPtPairCut::Pass(const AliFemtoPair* pair)
{
  bool temp = true;

  double pT1 = pair->Track1()->Track()->Pt();
  double pT2 = pair->Track2()->Track()->Pt();
  double DeltaPt = TMath::Abs(pT1-pT2);

  if (DeltaPt >= fDeltaPtMin && DeltaPt <= fDeltaPtMax)
    temp = false;

  return temp;
}
