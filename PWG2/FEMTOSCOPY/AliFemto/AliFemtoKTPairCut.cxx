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
  fPhiMax(360.0)
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
  fPhiMax(360)
{
}
//__________________
AliFemtoKTPairCut::AliFemtoKTPairCut(const AliFemtoKTPairCut& c) : 
  AliFemtoPairCut(c),
  fKTMin(0),
  fKTMax(1.0e6),
  fPhiMin(0),
  fPhiMax(360)
{ 
  fKTMin = c.fKTMin;
  fKTMax = c.fKTMax;
  fPhiMin = c.fPhiMin;
  fPhiMax = c.fPhiMax;
}

//__________________
AliFemtoKTPairCut::~AliFemtoKTPairCut(){
  /* no-op */
}
//__________________
bool AliFemtoKTPairCut::Pass(const AliFemtoPair* pair){
  bool temp = true;
  
  if (pair->KT() < fKTMin)
    temp = false;

  if (pair->KT() > fKTMax)
    temp = false;

  return temp;
}
//__________________
AliFemtoString AliFemtoKTPairCut::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoKT Pair Cut \n";  char ctemp[100];
  sprintf(ctemp,"Accept pair with kT in range %f , %f",fKTMin,fKTMax);
  sprintf(ctemp,"Accept pair with angle in range %f , %f",fPhiMin,fPhiMax);
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

bool AliFemtoKTPairCut::Pass(const AliFemtoPair* pair, double aRPAngle)
{
  if (!(Pass(pair))) return false;

  //  cout << "Got pair angle RP " << pair->EmissionAngle() << "   " << aRPAngle << endl;

  bool temp = true;
  double rpangle = pair->EmissionAngle();
  if (rpangle > 180.0) rpangle -= 180.0;
  rpangle -= aRPAngle*180/TMath::Pi();
  if (rpangle > 180.0) rpangle -= 180.0;
  if (rpangle < 0.0) rpangle += 180.0;

  //  cout << "Got difference " << rpangle << endl;

  if (fPhiMin < 0) {
    if ((rpangle > fPhiMax) && (rpangle < 180.0+fPhiMin)) 
      temp = false;
  }
  else {
    if ((rpangle < fPhiMin) || (rpangle > fPhiMax))
      temp = false;
  }
      
  //  if (temp) cout << "Accepted !" << endl;

  return temp;
}
