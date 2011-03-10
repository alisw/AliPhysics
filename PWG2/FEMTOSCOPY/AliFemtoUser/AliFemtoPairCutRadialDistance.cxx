/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoPairCutRadialDistance - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoPairCutRadialDistance.cxx,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
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

#include "AliFemtoPairCutRadialDistance.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutRadialDistance)
#endif

//__________________
AliFemtoPairCutRadialDistance::AliFemtoPairCutRadialDistance():
  AliFemtoShareQualityPairCut(),
  fDRadMin(0), 
  fRadius(100),
  fEtaMin(0)
{
}
//__________________
AliFemtoPairCutRadialDistance::AliFemtoPairCutRadialDistance(const AliFemtoPairCutRadialDistance& c) : 
  AliFemtoShareQualityPairCut(c),
  fDRadMin(0), 
  fRadius(100),
  fEtaMin(0)
{ 
  fDRadMin = c.fDRadMin;
  fRadius = c.fRadius;
  fEtaMin = c.fEtaMin;
}

//__________________
AliFemtoPairCutRadialDistance::~AliFemtoPairCutRadialDistance(){
  /* no-op */
}
//__________________
bool AliFemtoPairCutRadialDistance::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC entrance separation and
  // quality and sharity
  bool temp = true;
  
//   double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
//   double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
//   double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
//   double dist = sqrt(distx*distx + disty*disty + distz*distz);

//   temp = dist > fDRadMin;

  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double chg1 = pair->Track1()->Track()->Charge();
  double chg2 = pair->Track2()->Track()->Charge();
  double ptv1 = pair->Track1()->Track()->Pt();
  double ptv2 = pair->Track2()->Track()->Pt();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

  double dist = phi2 - phi1 + TMath::ASin(-0.3 * 0.5 * chg2 * fRadius/(2*ptv2)) - TMath::ASin(-0.3 * 0.5 * chg1 * fRadius/(2*ptv1));
  double etad = eta2 - eta1;

  temp = ((TMath::Abs(dist) > fDRadMin) || (TMath::Abs(etad) > fEtaMin));
  
  if (temp) {
    temp = AliFemtoShareQualityPairCut::Pass(pair);
  }
  else
    fNPairsFailed++;

  return temp;
}
//__________________
AliFemtoString AliFemtoPairCutRadialDistance::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoRadialDistance Pair Cut - remove shared and split pairs and pairs with small separation at the specified radius\n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with separation more that %f",fDRadMin);
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutRadialDistance::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoShareQualityPairCut::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutRadialDistance.radialsepmin=%f", fDRadMin);
  snprintf(buf, 200, "AliFemtoPairCutRadialDistance.radius=%f", fRadius);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutRadialDistance::SetRadialDistanceMinimum(double radius, double dtpc)
{
  fDRadMin = dtpc;
  fRadius = radius;
}

void AliFemtoPairCutRadialDistance::SetEtaDifferenceMinimum(double etpc) 
{
  fEtaMin = etpc;
}
