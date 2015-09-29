///
/// \file AliFemtoPairCutGeneralizedRadialDistance.cxx
///

#include "AliFemtoPairCutGeneralizedRadialDistance.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutGeneralizedRadialDistance)
#endif

//__________________
AliFemtoPairCutGeneralizedRadialDistance::AliFemtoPairCutGeneralizedRadialDistance():
AliFemtoPairCutAntiGamma(),
  fDPhiStarMin(0.0),
  fDEtaStarMin(0.0),
  fRadius(1.2)
{
}
//__________________
AliFemtoPairCutGeneralizedRadialDistance::AliFemtoPairCutGeneralizedRadialDistance(const AliFemtoPairCutGeneralizedRadialDistance& c) :
  AliFemtoPairCutAntiGamma(c),
  fDPhiStarMin(0.0),
  fDEtaStarMin(0.0),
  fRadius(1.2)
{
  fDPhiStarMin = c.fDPhiStarMin;
  fDEtaStarMin = c.fDEtaStarMin;
  fRadius = c.fRadius;
}

//__________________
AliFemtoPairCutGeneralizedRadialDistance::~AliFemtoPairCutGeneralizedRadialDistance(){
  /* no-op */
}
AliFemtoPairCutGeneralizedRadialDistance& AliFemtoPairCutGeneralizedRadialDistance::operator=(const AliFemtoPairCutGeneralizedRadialDistance& c)
{
  if (this != &c) {
    fDPhiStarMin = c.fDPhiStarMin;
    fDEtaStarMin = c.fDEtaStarMin;
    fRadius = c.fRadius;
  }

  return *this;
}
//__________________
bool AliFemtoPairCutGeneralizedRadialDistance::Pass(const AliFemtoPair* pair){

  const AliFemtoTrack *track1 = pair->Track1()->Track();
  const AliFemtoTrack *track2 = pair->Track2()->Track();

  double thetas1 = TMath::Pi()/2. - TMath::ATan(track1->NominalTpcPointShifted().z()/(fRadius*1e2));
  double thetas2 = TMath::Pi()/2. - TMath::ATan(track2->NominalTpcPointShifted().z()/(fRadius*1e2));
  double etas1 = -TMath::Log( TMath::Tan(thetas1/2.) );
  double etas2 = -TMath::Log( TMath::Tan(thetas2/2.) );
  double detas = TMath::Abs(etas1 - etas2);
  double distSft= TMath::Sqrt(TMath::Power(track1->NominalTpcPointShifted().x() -
					   track2->NominalTpcPointShifted().x(),2) +
			      TMath::Power(track1->NominalTpcPointShifted().y() -
					   track2->NominalTpcPointShifted().y(),2));
  double dPhiS = 2.0 * TMath::ATan(distSft/2./((fRadius*1e2)));

  if (TMath::Abs(detas) < fDEtaStarMin && TMath::Abs(dPhiS) < fDPhiStarMin)
    return false;
  return true;
}
//__________________
AliFemtoString AliFemtoPairCutGeneralizedRadialDistance::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoRadialDistance Pair Cut - remove shared and split pairs and pairs with small separation at the specified radius\n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with separation more that %f",fDPhiStarMin);
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutGeneralizedRadialDistance::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoPairCutAntiGamma::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutGeneralizedRadialDistance.phistarsepmin=%f", fDPhiStarMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutGeneralizedRadialDistance::SetPhiStarDifferenceMinimum(Double_t dps)
{
  fDPhiStarMin = dps;
}

void AliFemtoPairCutGeneralizedRadialDistance::SetEtaStarDifferenceMinimum(Double_t des)
{
  fDEtaStarMin = des;
}


void AliFemtoPairCutGeneralizedRadialDistance::SetRadius(Double_t rad)
{
  fRadius = rad;
}
