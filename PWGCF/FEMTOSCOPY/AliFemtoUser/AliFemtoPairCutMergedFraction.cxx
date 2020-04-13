/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutMergedFraction - a pair cut which checks                     //
// for some pair qualities that attempt to identify split/doubly               //
// reconstructed tracks, calculating the fraction of "merged" points (number   //
// of points in selected range of radius from primary vertex where tracks are  //
// closer than selected distance divided by total number of calculated points) //
//                                                                             //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch            //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoPairCutMergedFraction.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutMergedFraction)
#endif

//__________________
AliFemtoPairCutMergedFraction::AliFemtoPairCutMergedFraction(Double_t aDistanceMax, Double_t aMergedFractionLimit, Double_t aDEtaMax, Double_t aRadiusMin, Double_t aRadiusMax) :
AliFemtoPairCutAntiGamma(),
  fDistanceMax(0.03),
  fMergedFractionLimit(0.01),
  fDEtaMax(0.01),
  fRadiusMin(0.8),
  fRadiusMax(2.5),
  fMagSign(1),
  fMagFieldVal(0.5),
  fMergedFractionDataType(kESD)
{
  fDistanceMax = aDistanceMax;
  fMergedFractionLimit = aMergedFractionLimit;
  fDEtaMax = aDEtaMax;
  fRadiusMin = aRadiusMin;
  fRadiusMax = aRadiusMax;
}

//__________________
AliFemtoPairCutMergedFraction::AliFemtoPairCutMergedFraction(const AliFemtoPairCutMergedFraction& cPairCut) :
  AliFemtoPairCutAntiGamma(cPairCut),
  fDistanceMax(cPairCut.fDistanceMax),
  fMergedFractionLimit(cPairCut.fMergedFractionLimit),
  fDEtaMax(cPairCut.fDEtaMax),
  fRadiusMin(cPairCut.fRadiusMin),
  fRadiusMax(cPairCut.fRadiusMax),
  fMagSign(cPairCut.fMagSign),
  fMagFieldVal(cPairCut.fMagFieldVal),
  fMergedFractionDataType(cPairCut.fMergedFractionDataType)
{
}

//__________________
AliFemtoPairCutMergedFraction::~AliFemtoPairCutMergedFraction() {

}

AliFemtoPairCutMergedFraction& AliFemtoPairCutMergedFraction::operator=(const AliFemtoPairCutMergedFraction& cPairCut) {
  if(this != &cPairCut) {
    AliFemtoPairCutAntiGamma::operator=(cPairCut);
    fDistanceMax = cPairCut.fDistanceMax;
    fMergedFractionLimit = cPairCut.fMergedFractionLimit;
    fDEtaMax = cPairCut.fDEtaMax;
    fRadiusMin = cPairCut.fRadiusMin;
    fRadiusMax = cPairCut.fRadiusMax;
    fMagSign = cPairCut.fMagSign;
    fMagFieldVal = cPairCut.fMagFieldVal;
    fMergedFractionDataType = cPairCut.fMergedFractionDataType;
  }
  return *this;
}

//__________________
bool AliFemtoPairCutMergedFraction::Pass(const AliFemtoPair* pair) {

  if(fMergedFractionDataType == kKine)
    return true;

  // Prepare variables:
  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double chg1 = pair->Track1()->Track()->Charge();
  double chg2 = pair->Track2()->Track()->Charge();
  double pt1 = pair->Track1()->Track()->Pt();
  double pt2 = pair->Track2()->Track()->Pt();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();
  double magval = fMagFieldVal;

  // Check magnetic field sign:
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  Double_t magsign = 0.0;
  if (!aodH) {
    return false;
  }
  else {
    AliAODEvent *fAOD;
    fAOD = aodH->GetEvent();
    magsign = fAOD->GetMagneticField();
  }
  if (magsign > 1)
    fMagSign = 1;
  else if ( magsign < 1)
    fMagSign = -1;
  else
    fMagSign = magsign;

  // Calculate dEta:
  double deta = eta2 - eta1;

  // Check the correctness of parameters:
  if(fRadiusMin != fRadiusMax) {
    if(fRadiusMin > fRadiusMax) {
      Double_t radiustmp = fRadiusMin;
      fRadiusMin = fRadiusMax;
      fRadiusMax = radiustmp;
    }
  }
  if(fMergedFractionLimit < 0.0) {
    fMergedFractionLimit = - fMergedFractionLimit;
  }
  if(fMergedFractionLimit > 1.0) {
    return true;
  }
  if(fDistanceMax < 0.0) {
    fDistanceMax = - fDistanceMax;
  }

  Bool_t pairpass = kTRUE;
  Double_t badpoints = 0.0;
  Double_t allpoints = 0.0;

  // If dEta is not in range:
  if(TMath::Abs(deta) > TMath::Abs(fDEtaMax))
    return kTRUE;
  // Iterate through all radii in range (fRadiusMin, fRadiusMax):
  for(double irad = fRadiusMin; irad < fRadiusMax; irad += 0.01) {

    // Calculate radius:
    Double_t rad = irad;

    // Calculate dPhiStar:
    //double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
    //double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
    double afsi0b = -0.15*magval*chg1*fMagSign*rad/pt1;
    double afsi1b = -0.15*magval*chg2*fMagSign*rad/pt2; 
    Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
    dphistar = TVector2::Phi_mpi_pi(dphistar); // returns phi angle in the interval [-PI,PI)

    // Calculate distance:
    //double distance = TMath::Sqrt(rad * rad * (2 - 2 * TMath::Cos(dphistar)));
    double distance = 2 * TMath::Sin(TMath::Abs(dphistar) * 0.5) * rad;

    // Check if pair parameters meet the requirements:
    if(distance < fDistanceMax) {
      badpoints += 1.0;
    }
    allpoints += 1.0;
  }

  if(allpoints != 0.0) {
    // Calculate fraction:
    Double_t fraction = badpoints / allpoints;

    // Remove pair if the fraction is above limit:
    if(fraction > fMergedFractionLimit) {
      pairpass = kFALSE;
    }
  }
  else {
    pairpass = kTRUE;
  }

  // Check the antigamma cut:
  if (pairpass) {
    pairpass = AliFemtoPairCutAntiGamma::Pass(pair);
  }
  else {
    fNPairsFailed++;
  }

  return pairpass;
}

//__________________
AliFemtoString AliFemtoPairCutMergedFraction::Report()
{
  // Prepare a report from the execution
  AliFemtoString report = "AliFemtoPairCutMergedFraction Pair Cut";
  report += "- remove shared and split pairs and pairs with small separation at the specified radius\n";
  report += Form("Accept pair with separation more than %f fraction in %f m distance", fMergedFractionLimit, fDistanceMax);
  return report;
}

//__________________
TList *AliFemtoPairCutMergedFraction::ListSettings() {
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoPairCut::ListSettings();
  tListSetttings->AddLast(new TObjString(Form("AliFemtoPairCutMergedFraction.radiusrange=(%f,%f)", fRadiusMin, fRadiusMax)));

  return tListSetttings;
}

//__________________
void AliFemtoPairCutMergedFraction::SetDistanceMax(double maxdistance) {
  fDistanceMax = maxdistance;
}

//__________________
void AliFemtoPairCutMergedFraction::SetMergedFractionLimit(double fractionlimit) {
  fMergedFractionLimit = fractionlimit;
}

//__________________
void AliFemtoPairCutMergedFraction::SetDEtaMax(double maxeta) {
  fDEtaMax = maxeta;
}

//__________________
void AliFemtoPairCutMergedFraction::SetMagneticFieldSign(int magsign) {
  if(magsign>1)
    fMagSign = 1;
  else if(magsign<1)
    fMagSign = -1;
  else
    fMagSign = magsign;
}

void AliFemtoPairCutMergedFraction::SetMagneticFieldValue(double magval) {
  fMagFieldVal = magval;
}

//__________________
void AliFemtoPairCutMergedFraction::SetRadiusMin(double radmin) {
  fRadiusMin = radmin;
}

//__________________
void AliFemtoPairCutMergedFraction::SetRadiusMax(double radmax) {
  fRadiusMax = radmax;
}

//__________________
void AliFemtoPairCutMergedFraction::SetMergedFractionDataType(AliFemtoDataType datatype) {
  fMergedFractionDataType = datatype;
}
