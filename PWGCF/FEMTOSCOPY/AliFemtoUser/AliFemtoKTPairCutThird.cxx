/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoKTPairCutThird - a pair cut which selects pairs based on their       //
// transverse momentum kT                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoKTPairCutThird.cxx,v 1.1.2.2 2007/11/09 11:20:35 akisiel Exp $
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
#include "AliFemtoKTPairCutThird.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoKTPairCutThird)
#endif

//__________________
AliFemtoKTPairCutThird::AliFemtoKTPairCutThird():
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
AliFemtoKTPairCutThird::AliFemtoKTPairCutThird(double lo, double hi) :
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
AliFemtoKTPairCutThird::AliFemtoKTPairCutThird(const AliFemtoKTPairCutThird& c) : 
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
AliFemtoKTPairCutThird::~AliFemtoKTPairCutThird(){
  /* no-op */
}
AliFemtoKTPairCutThird& AliFemtoKTPairCutThird::operator=(const AliFemtoKTPairCutThird& c)
{
  if (this != &c) {
    fKTMin = c.fKTMin;
    fKTMax = c.fKTMax;
    fPhiMin = c.fPhiMin;
    fPhiMax = c.fPhiMax;
    fPtMin = c.fPtMin;
    fPtMax = c.fPtMax;
  }

  return *this;
}
//__________________
/*bool AliFemtoKTPairCutThird::Pass(const AliFemtoPair* pair){
  bool temp = true;
  
  if (pair->KT() < fKTMin)
    temp = false;

  if (pair->KT() > fKTMax)
    temp = false;

  return temp;
}*/
//__________________
AliFemtoString AliFemtoKTPairCutThird::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoKT Pair Cut \n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with kT in range %f , %f",fKTMin,fKTMax);
  snprintf(ctemp , 100, "Accept pair with angle in range %f , %f",fPhiMin,fPhiMax);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoKTPairCutThird::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoKTPairCutThird.ktmax=%f", fKTMax);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCutThird.ktmin=%f", fKTMin);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCutThird.phimax=%f", fPhiMax);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCutThird.phimin=%f", fPhiMin);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCutThird.ptmin=%f", fPtMin);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKTPairCutThird.ptmax=%f", fPtMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoKTPairCutThird::SetKTRange(double ktmin, double ktmax)
{
  fKTMin = ktmin;
  fKTMax = ktmax;
}

void AliFemtoKTPairCutThird::SetPhiRange(double phimin, double phimax)
{
  fPhiMin = phimin;
  fPhiMax = phimax;
}

void AliFemtoKTPairCutThird::SetPTMin(double ptmin, double ptmax)
{
  fPtMin = ptmin;
  fPtMax = ptmax;
}

//______________________________________________________
bool AliFemtoKTPairCutThird::Pass(const AliFemtoPair* pair)
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

  if (rpangle > 120.0) rpangle -= 120.0;
  if (rpangle < 0.0) rpangle += 120.0;
  
  if (fPhiMin < 0) {
    if ((rpangle > fPhiMax) && (rpangle < 120.0+fPhiMin))
      temp = false;
  }
  else {
    if ((rpangle < fPhiMin) || (rpangle > fPhiMax))
      temp = false;
  }
  return temp;
}

//_____________________________________
bool AliFemtoKTPairCutThird::Pass(const AliFemtoPair* pair, double aRPAngle)
{
//The same as above, but it is defined with RP Angle as input in all the Correlatin function classes.

  bool temp = (aRPAngle > 0.);
  aRPAngle = true;
   
  if (!Pass(pair))
	temp = false;

  return temp;
}
