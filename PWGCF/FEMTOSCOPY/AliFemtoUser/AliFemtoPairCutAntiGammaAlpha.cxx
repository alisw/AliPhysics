/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoPairCutAntiGammaAlpha - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoPairCutAntiGammaAlpha.cxx,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
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

#include "AliFemtoPairCutAntiGammaAlpha.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutAntiGammaAlpha)
#endif

//__________________
AliFemtoPairCutAntiGammaAlpha::AliFemtoPairCutAntiGammaAlpha():
  AliFemtoShareQualityPairCut(),
  fMaxEEMinv(0.0),
  fMaxDAlpha(0.0),
  fDTPCMin(0),
  fDataType(kESD)
{
}
//__________________
AliFemtoPairCutAntiGammaAlpha::AliFemtoPairCutAntiGammaAlpha(const AliFemtoPairCutAntiGammaAlpha& c) :
  AliFemtoShareQualityPairCut(c),
  fMaxEEMinv(c.fMaxEEMinv),
  fMaxDAlpha(c.fMaxDAlpha),
  fDTPCMin(c.fDTPCMin),
  fDataType(c.fDataType)
{
}

AliFemtoPairCutAntiGammaAlpha& AliFemtoPairCutAntiGammaAlpha::operator=(const AliFemtoPairCutAntiGammaAlpha& c)
{
  if (this != &c) {
    AliFemtoShareQualityPairCut::operator=(c);
    fMaxEEMinv = c.fMaxEEMinv;
    fMaxDAlpha = c.fMaxDAlpha;
    fDTPCMin = c.fDTPCMin;
    fDataType = c.fDataType;
  }

  return *this;

}
//__________________
AliFemtoPairCutAntiGammaAlpha::~AliFemtoPairCutAntiGammaAlpha(){
  /* no-op */
}
//__________________
bool AliFemtoPairCutAntiGammaAlpha::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC entrance separation and
  // quality and sharity
  bool temp = true;

  if(fDataType==kKine)
    return true;

  double me = 0.000511;

  AliFemtoTrack *track1 = pair->Track1()->Track(),
                *track2 = pair->Track2()->Track();

  AliFemtoThreeVector p1 = track1->P(),
                      p2 = track2->P();

  if ((track1->Charge() * track2->Charge()) < 0.0) {

    double dalpha = TMath::Abs(p1.Dot(p2) / TMath::Sqrt(p1.Mag2() * p2.Mag2()));

    double e1 = TMath::Sqrt(me*me + p1.Mag2());
    double e2 = TMath::Sqrt(me*me + p2.Mag2());

    double minv = 2*me*me + 2*(e1*e2 - p1.Dot(p2));

    double sminv = TMath::Sqrt(minv);

    if ((sminv < fMaxEEMinv) && (dalpha > fMaxDAlpha)) {
      temp = false;
    }
  }

  bool tempTPCEntrance = true;

  if(fDataType==kESD || fDataType==kAOD)
    {
      auto &epoint1 = track1->NominalTpcEntrancePoint(),
           &epoint2 = track2->NominalTpcEntrancePoint();

      tempTPCEntrance = (epoint1 - epoint2).Mag2() > fDTPCMin * fDTPCMin;
    }


  if (temp && tempTPCEntrance) {

    temp = AliFemtoShareQualityPairCut::Pass(pair);
    if (temp) {
      fNPairsPassed++;
    }
    else fNPairsFailed++;
    return temp;
  }
  else {
    fNPairsFailed++;
    return false;
  }
}
//__________________
AliFemtoString AliFemtoPairCutAntiGammaAlpha::Report()
{
  // Prepare a report from the execution
  AliFemtoString report = "AliFemtoPairCutAntiGammaAlpha Pair Cut - remove pairs possibly coming from Gamma conversions\n";
  report += Form("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);
  return report;
}

TList *AliFemtoPairCutAntiGammaAlpha::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoShareQualityPairCut::ListSettings();

  tListSetttings->AddVector(
      new TObjString(Form("AliFemtoPairCutAntiGammaAlpha.maxeeminv=%f", fMaxEEMinv)),
      new TObjString(Form("AliFemtoPairCutAntiGammaAlpha.maxdalpha=%f", fMaxDAlpha)),
      new TObjString(Form("AliFemtoPairCutAntiGammaAlpha.mintpcentrancesep=%f", fDTPCMin)),
      nullptr
  );

  return tListSetttings;
}

void AliFemtoPairCutAntiGammaAlpha::SetMaxEEMinv(Double_t maxeeminv)
{
  fMaxEEMinv = maxeeminv;
}

void AliFemtoPairCutAntiGammaAlpha::SetMaxAlphaDiff(Double_t maxdalpha)
{
  fMaxDAlpha = maxdalpha;
}

void AliFemtoPairCutAntiGammaAlpha::SetTPCEntranceSepMinimum(double dtpc)
{
  fDTPCMin = dtpc;
}

void AliFemtoPairCutAntiGammaAlpha::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}
