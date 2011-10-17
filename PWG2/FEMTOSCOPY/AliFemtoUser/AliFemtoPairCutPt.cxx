/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoPairCutPt - a pair cut which checks if the sum of the transverse        //
// momenta of two particles fit within given range 		                   //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//  	       				                                           //
/////////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoPairCutPt.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutPt)
#endif

//__________________
AliFemtoPairCutPt::AliFemtoPairCutPt():
  AliFemtoPairCut(),
  fSumPtMin(0),
  fSumPtMax(10000),
  fNPairsFailed(0),
  fNPairsPassed(0)
{

}
//__________________
AliFemtoPairCutPt::AliFemtoPairCutPt(double lo, double hi):
  AliFemtoPairCut(),
  fSumPtMin(lo),
  fSumPtMax(hi),
  fNPairsFailed(0),
  fNPairsPassed(0)
{
  fSumPtMin=lo;
  fSumPtMax=hi;
}
//__________________
AliFemtoPairCutPt::AliFemtoPairCutPt(const AliFemtoPairCutPt& c) : 
  AliFemtoPairCut(c),
  fSumPtMin(0),
  fSumPtMax(0),
  fNPairsFailed(0),
  fNPairsPassed(0)
{ 
  fSumPtMin = c.fSumPtMin;
  fSumPtMax = c.fSumPtMax;
}

//__________________
AliFemtoPairCutPt::~AliFemtoPairCutPt(){
  /* no-op */
}
//__________________
bool AliFemtoPairCutPt::Pass(const AliFemtoPair* pair){

  bool temp = true;

  double pt1 = pair->Track1()->Track()->Pt();
  double pt2 = pair->Track2()->Track()->Pt();

  double pt_sum = pt1 + pt2;

  if(pt_sum >= fSumPtMin && pt_sum <= fSumPtMax)
    temp = true;
  else
    temp = false;

  if(temp) 
    fNPairsPassed++;
  else fNPairsFailed++;


  return temp;
  
}
//__________________
AliFemtoString AliFemtoPairCutPt::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoPairCutPt Pair Cut\n";  
  char ctemp[100];
  stemp += ctemp;
  snprintf(ctemp,100,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",(long int) fNPairsPassed,(long int) fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutPt::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutPt.sumptmin=%f", fSumPtMin);
  snprintf(buf, 200, "AliFemtoPairCutPr.sumptmax=%f", fSumPtMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutPt::SetMinSumPt(Double_t sumptmin)
{
  fSumPtMin = sumptmin;
}

 
void AliFemtoPairCutPt::SetMaxSumPt(Double_t sumptmax)
{
  fSumPtMax = sumptmax;
}

