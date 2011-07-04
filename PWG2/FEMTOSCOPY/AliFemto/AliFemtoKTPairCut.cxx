/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoKTPairCut - a pair cut which selects pairs based on their       //
// transverse momentum kT                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoKTPairCut.cxx,v 1.1.2.2 2007/11/09 11:20:35 akisiel Exp $
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
#include "AliFemtoKTPairCut.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoKTPairCut)
#endif

//__________________
AliFemtoKTPairCut::AliFemtoKTPairCut():
  AliFemtoPairCut(),
  fKTMin(0),
  fKTMax(1.0e6),
  fPhiMin(0),
  fPhiMax(360.0),
  fPtMin(0.0),
  fPtMax(1000.0)
{
  fKTMin = 0;
   fKTMax = 1.0e6;
}
//__________________
AliFemtoKTPairCut::AliFemtoKTPairCut(double lo, double hi) :
  AliFemtoPairCut(),
  fKTMin(lo),
  fKTMax(hi),
  fPhiMin(0),
  fPhiMax(360),
  fPtMin(0.0),
  fPtMax(1000.0)
{
}
//__________________
AliFemtoKTPairCut::AliFemtoKTPairCut(const AliFemtoKTPairCut& c) : 
  AliFemtoPairCut(c),
  fKTMin(0),
  fKTMax(1.0e6),
  fPhiMin(0),
  fPhiMax(360),
  fPtMin(0.0),
  fPtMax(1000.0)
{ 
  fKTMin = c.fKTMin;
  fKTMax = c.fKTMax;
  fPhiMin = c.fPhiMin;
  fPhiMax = c.fPhiMax;
  fPtMin = c.fPtMin;
  fPtMax = c.fPtMax;
}

//__________________
AliFemtoKTPairCut::~AliFemtoKTPairCut(){
  /* no-op */
}
//__________________
/*bool AliFemtoKTPairCut::Pass(const AliFemtoPair* pair){
  bool temp = true;
  
  if (pair->KT() < fKTMin)
    temp = false;

  if (pair->KT() > fKTMax)
    temp = false;

  return temp;
}*/
//__________________
AliFemtoString AliFemtoKTPairCut::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoKT Pair Cut \n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with kT in range %f , %f",fKTMin,fKTMax);
  snprintf(ctemp , 100, "Accept pair with angle in range %f , %f",fPhiMin,fPhiMax);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoKTPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoKTPairCut.ktmax=%f", fKTMax);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCut.ktmin=%f", fKTMin);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCut.phimax=%f", fPhiMax);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCut.phimin=%f", fPhiMin);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCut.ptmin=%f", fPtMin);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCut.ptmax=%f", fPtMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoKTPairCut::SetKTRange(double ktmin, double ktmax)
{
  fKTMin = ktmin;
  fKTMax = ktmax;
}

void AliFemtoKTPairCut::SetPhiRange(double phimin, double phimax)
{
  fPhiMin = phimin;
  fPhiMax = phimax;
}

void AliFemtoKTPairCut::SetPTMin(double ptmin, double ptmax)
{
  fPtMin = ptmin;
  fPtMax = ptmax;
}

//______________________________________________________
bool AliFemtoKTPairCut::Pass(const AliFemtoPair* pair)
{
  bool temp = true;

//Taking care of the Kt cut
  if (pair->KT() < fKTMin)
    temp = false;

  if (pair->KT() > fKTMax)
    temp = false;

  if (!temp) return temp;

  if ((fPtMin > 0.0) || (fPtMax<1000.0)) {
//     double px1 = pair->Track1()->Track()->P().x();
//     double py1 = pair->Track1()->Track()->P().y();

//     double px2 = pair->Track2()->Track()->P().x();
//     double py2 = pair->Track2()->Track()->P().y();
    
//     double pt1 = TMath::Hypot(px1, py1);
//     double pt2 = TMath::Hypot(px2, py2);
    
//     if ((pt1<fPtMin) || (pt1>fPtMax)) return false;
//     if ((pt2<fPtMin) || (pt2>fPtMax)) return false;
    if ((pair->Track1()->Track()->Pt()<fPtMin) || (pair->Track1()->Track()->Pt()>fPtMax)) return false;
    if ((pair->Track2()->Track()->Pt()<fPtMin) || (pair->Track2()->Track()->Pt()>fPtMax)) return false;
  }

//Taking care of the Phi cut
//   double rpangle = (pair->GetPairAngleEP())*180/TMath::Pi();
  double rpangle = pair->GetPairAngleEP();

  if (rpangle > 180.0) rpangle -= 180.0;
  if (rpangle < 0.0) rpangle += 180.0;
  
  if (fPhiMin < 0) {
    if ((rpangle > fPhiMax) && (rpangle < 180.0+fPhiMin)) 
      temp = false;
  }
  else {
    if ((rpangle < fPhiMin) || (rpangle > fPhiMax))
      temp = false;
  }
  return temp;
}

//_____________________________________
bool AliFemtoKTPairCut::Pass(const AliFemtoPair* pair, double aRPAngle)
{
//The same as above, but it is defined with RP Angle as input in all the Correlatin function classes.

  bool temp = (aRPAngle > 0.);
  aRPAngle = true;
   
  if (!Pass(pair))
	temp = false;

  return temp;
}
