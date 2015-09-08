/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoBetaTPairCut - a pair cut which selects pairs based on their       //
// transverse momentum kT                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliFemtoBetaTPairCut.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoBetaTPairCut)
#endif

//__________________
AliFemtoBetaTPairCut::AliFemtoBetaTPairCut():
  AliFemtoPairCut(),
  fBetaTMin(0.0),
  fBetaTMax(1.0e6)
{
  fBetaTMin = 0;
  fBetaTMax = 1.0e6;
}

//__________________
AliFemtoBetaTPairCut::AliFemtoBetaTPairCut(double minbetat, double maxbetat) :
  AliFemtoPairCut(),
  fBetaTMin(minbetat),
  fBetaTMax(maxbetat)
{
}

//__________________
AliFemtoBetaTPairCut::AliFemtoBetaTPairCut(const AliFemtoBetaTPairCut& c) : 
  AliFemtoPairCut(c),
  fBetaTMin(0.0),
  fBetaTMax(1.0e6)
{ 
  fBetaTMin = c.fBetaTMin;
  fBetaTMax = c.fBetaTMax;
}

//__________________
AliFemtoBetaTPairCut::~AliFemtoBetaTPairCut() {
  /* no-op */
}

//__________________
AliFemtoBetaTPairCut& AliFemtoBetaTPairCut::operator=(const AliFemtoBetaTPairCut& c) {
  if (this != &c) {
    fBetaTMin = c.fBetaTMin;
    fBetaTMax = c.fBetaTMax;
  }
  return *this;
}

//__________________
AliFemtoString AliFemtoBetaTPairCut::Report() {
  // Prepare a report from the execution
  string stemp = "AliFemtoBetaT Pair Cut \n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with betaT in range %f , %f",fBetaTMin,fBetaTMax);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;
}

//__________________
TList *AliFemtoBetaTPairCut::ListSettings() {
  // return a list of settings in a writable form
  TList *tListSetttings =  new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoBetaTPairCut.betatmax=%f", fBetaTMax);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBetaTPairCut.betatmin=%f", fBetaTMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

//__________________
void AliFemtoBetaTPairCut::SetBetaTRange(double minbetat, double maxbetat) {
  fBetaTMin = minbetat;
  fBetaTMax = maxbetat;
}

//______________________________________________________
bool AliFemtoBetaTPairCut::Pass(const AliFemtoPair* pair) {

  bool pairpass = true;

  // Prepare variables:
  double pt1 = pair->Track1()->Track()->Pt();
  double pt2 = pair->Track2()->Track()->Pt();
  double pTpair = pt1 + pt2;
  double px1 = pair->Track1()->Track()->P().x();
  double px2 = pair->Track2()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  double py2 = pair->Track2()->Track()->P().y();
  double pz1 = pair->Track1()->Track()->P().z();
  double pz2 = pair->Track2()->Track()->P().z();
  double pxpair = px1 + px2;
  double pypair = py1 + py2;
  double pzpair = pz1 + pz2;
  double p1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);
  double p2 = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2);
  double ppair = TMath::Sqrt(pxpair*pxpair + pypair*pypair + pzpair*pzpair);
  double m1 = pair->Track1()->Track()->GetMass();
  double m2 = pair->Track2()->Track()->GetMass();
  double e1 = TMath::Sqrt(p1*p1 + m1*m1);
  double e2 = TMath::Sqrt(p2*p2 + m2*m2);
  double epair = e1 + e2;
  double mTpair = TMath::Sqrt(epair*epair - pzpair*pzpair);
  double betaT = TMath::Abs(pTpair) / mTpair;

  // Check if betaT is in range:
  if(betaT < fBetaTMin)
    pairpass = false;
  if(betaT > fBetaTMax)
    pairpass = false;

  return pairpass;
}
