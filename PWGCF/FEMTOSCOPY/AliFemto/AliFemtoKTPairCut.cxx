///
/// \class AliFemtoKTPairCut.cxx
///

#include "AliFemtoKTPairCut.h"

#include <TMath.h>
#include <string>


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoKTPairCut);
  /// \endcond
#endif

AliFemtoKTPairCut::AliFemtoKTPairCut():
  AliFemtoKTPairCut(0, 1.0e6)
{
}

AliFemtoKTPairCut::AliFemtoKTPairCut(double lo, double hi):
  AliFemtoPairCut()
  , fKTMin(lo)
  , fKTMax(hi)
  , fPhiMin(0)
  , fPhiMax(360)
  , fPtMin(0.0)
  , fPtMax(1000.0)
{
}

AliFemtoKTPairCut::AliFemtoKTPairCut(const AliFemtoKTPairCut& c):
  AliFemtoPairCut(c)
  , fKTMin(c.fKTMin)
  , fKTMax(c.fKTMax)
  , fPhiMin(c.fPhiMin)
  , fPhiMax(c.fPhiMax)
  , fPtMin(c.fPtMin)
  , fPtMax(c.fPtMax)
{
}

AliFemtoKTPairCut::~AliFemtoKTPairCut()
{
}

AliFemtoKTPairCut& AliFemtoKTPairCut::operator=(const AliFemtoKTPairCut& c)
{
  if (this != &c) {
    AliFemtoPairCut::operator=(c);
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
/*bool AliFemtoKTPairCut::Pass(const AliFemtoPair* pair){
  bool temp = true;

  if (pair->KT() < fKTMin)
    temp = false;

  if (pair->KT() > fKTMax)
    temp = false;

  return temp;
}*/

AliFemtoString AliFemtoKTPairCut::Report()
{
  // Prepare a report from the execution
  AliFemtoString report = "AliFemtoKT Pair Cut \n";
  report += Form("Accept pair with kT in range %f , %f",fKTMin,fKTMax);
  report += Form("Accept pair with pT in range %f , %f", fPtMin,fPtMax);
  report += Form("Accept pair with angle in range %f , %f",fPhiMin,fPhiMax);

  return report;
}

TList *AliFemtoKTPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *settings = new TList();

  settings->AddVector(
    new TObjString(TString::Format("AliFemtoKTPairCut.ktmax=%f", fKTMax)),
    new TObjString(TString::Format("AliFemtoKTPairCut.ktmin=%f", fKTMin)),
    new TObjString(TString::Format("AliFemtoKTPairCut.phimax=%f", fPhiMax)),
    new TObjString(TString::Format("AliFemtoKTPairCut.phimin=%f", fPhiMin)),
    new TObjString(TString::Format("AliFemtoKTPairCut.ptmin=%f", fPtMin)),
    new TObjString(TString::Format("AliFemtoKTPairCut.ptmax=%f", fPtMax)),
    nullptr
  );

  return settings;
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

bool AliFemtoKTPairCut::Pass(const AliFemtoPair* pair)
{
  // Taking care of the kT cut
  if (pair->KT() < fKTMin || fKTMax <= pair->KT()) {
    return false;
  }

  // Taking care of the pT cut
  if ((fPtMin > 0.0) || (fPtMax < 1000.0)) {
    const double pt1 = pair->Track1()->Track()->Pt(),
                 pt2 = pair->Track2()->Track()->Pt();

    if (pt1 < fPtMin || fPtMax <= pt1) {
      return false;
    }

    if (pt2 < fPtMin || fPtMax <= pt2) {
      return false;
    }
  }

  // Taking care of the Phi cut
  double rpangle = pair->GetPairAngleEP();

  if (rpangle > 180.0) rpangle -= 180.0;
  if (rpangle < 0.0) rpangle += 180.0;

  // note: handle "wrap around" if the minimum phi is negative
  if (fPhiMin < 0.0) {
    return (rpangle < fPhiMax) || (180.0+fPhiMin <= rpangle);
  }

  // return whether angle is within phi-range
  return (fPhiMin <= rpangle) && (rpangle < fPhiMax);
}

bool AliFemtoKTPairCut::Pass(const AliFemtoPair* pair, double aRPAngle)
{
  // The same as above, but it is defined with RP Angle as input in
  // all the Correlation function classes.

  bool passes = (aRPAngle > 0.) && Pass(pair);

  return passes;
}
