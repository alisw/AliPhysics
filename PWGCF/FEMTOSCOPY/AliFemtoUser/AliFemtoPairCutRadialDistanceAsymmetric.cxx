/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistance - a pair cut which checks                     //
// for some pair qualities that attempt to identify split/doubly               //
// reconstructed tracks and also selects pairs based on their separation       //
// at the entrance to the TPC, adjusted to perform asymmetric cut for pair     //
//                                                                             //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch            //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoPairCutRadialDistanceAsymmetric.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutRadialDistanceAsymmetric)
#endif

//__________________
AliFemtoPairCutRadialDistanceAsymmetric::AliFemtoPairCutRadialDistanceAsymmetric(Double_t aRadiusMin, Double_t aRadiusMax, Bool_t aCalculateRadiusRange, Double_t aDEtaRangeLow, Double_t aDEtaRangeUp, Double_t aDPhiStarRangeLow, Double_t aDPhiStarRangeUp) :
AliFemtoPairCutAntiGamma(),
  fRadiusMin(0.8),
  fRadiusMax(2.5),
  fCalculateRadiusRange(true),
  fDEtaRangeLow(0.0),
  fDEtaRangeUp(0.0),
  fDPhiStarRangeLow(0.0),
  fDPhiStarRangeUp(0.0),
  fMagSign(1)
  
{
  fRadiusMin = aRadiusMin;
  fRadiusMax = aRadiusMax;
  fCalculateRadiusRange = aCalculateRadiusRange;
  fDEtaRangeLow = aDEtaRangeLow;
  fDEtaRangeUp = aDEtaRangeUp;
  fDPhiStarRangeLow = aDPhiStarRangeLow;
  fDPhiStarRangeUp = aDPhiStarRangeUp;
}

//__________________
AliFemtoPairCutRadialDistanceAsymmetric::AliFemtoPairCutRadialDistanceAsymmetric(const AliFemtoPairCutRadialDistanceAsymmetric& cPairCut) :
  AliFemtoPairCutAntiGamma(cPairCut),
  fRadiusMin(0.8),
  fRadiusMax(2.5),
  fCalculateRadiusRange(true),
  fDEtaRangeLow(0.0),
  fDEtaRangeUp(0.0),
  fDPhiStarRangeLow(0.0),
  fDPhiStarRangeUp(0.0),
  fMagSign(1)
{
  fRadiusMin = cPairCut.fRadiusMin;
  fRadiusMax = cPairCut.fRadiusMax;
  fCalculateRadiusRange = cPairCut.fCalculateRadiusRange;
  fDEtaRangeLow = cPairCut.fDEtaRangeLow;
  fDEtaRangeUp = cPairCut.fDEtaRangeUp;
  fDPhiStarRangeLow = cPairCut.fDPhiStarRangeLow;
  fDPhiStarRangeUp = cPairCut.fDPhiStarRangeUp;
  fMagSign = cPairCut.fMagSign;
}

//__________________
  AliFemtoPairCutRadialDistanceAsymmetric::~AliFemtoPairCutRadialDistanceAsymmetric(){
  
}

AliFemtoPairCutRadialDistanceAsymmetric& AliFemtoPairCutRadialDistanceAsymmetric::operator=(const AliFemtoPairCutRadialDistanceAsymmetric& cPairCut) {
  if(this != &cPairCut) {
    fRadiusMin = cPairCut.fRadiusMin;
    fRadiusMax = cPairCut.fRadiusMax;
    fCalculateRadiusRange = cPairCut.fCalculateRadiusRange;
    fDEtaRangeLow = cPairCut.fDEtaRangeLow;
    fDEtaRangeUp = cPairCut.fDEtaRangeUp;
    fDPhiStarRangeLow = cPairCut.fDPhiStarRangeLow;
    fDPhiStarRangeUp = cPairCut.fDPhiStarRangeUp;
    fMagSign = cPairCut.fMagSign;
  }
  return *this;
}

//__________________
bool AliFemtoPairCutRadialDistanceAsymmetric::Pass(const AliFemtoPair* pair) {

  if(AliFemtoPairCutAntiGamma::fDataType == kKine)
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
  if(fDPhiStarRangeLow != fDPhiStarRangeUp) {
    if(fDPhiStarRangeLow > fDPhiStarRangeUp) {
      Double_t dphistartmp = fDPhiStarRangeLow;
      fDPhiStarRangeLow = fDPhiStarRangeUp;
      fDPhiStarRangeUp = dphistartmp;
    }
  }
  if(fDEtaRangeLow != fDEtaRangeUp) {
    if(fDEtaRangeLow > fDEtaRangeUp) {
      Double_t detatmp = fDEtaRangeLow;
      fDEtaRangeLow = fDEtaRangeUp;
      fDEtaRangeUp = detatmp;
    }
  }

  Bool_t pairpass = kTRUE;
  
  // Iterate through all radii in range (fRadiusMin, fRadiusMax):
  if(fCalculateRadiusRange == true) {
    for(double irad = fRadiusMin; irad < fRadiusMax; irad = irad + 0.01) {
    
      // Calculate radius:
      Double_t rad = irad;

      // Calculate dPhiStar:
      double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
      double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
      Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
      dphistar = TVector2::Phi_mpi_pi(dphistar); // returns phi angle in the interval [-PI,PI)

      // Check if pair parameters meet the requirements:
      if(deta > fDEtaRangeLow && deta < fDEtaRangeUp) {
	if(dphistar > fDPhiStarRangeLow && dphistar < fDPhiStarRangeUp) {
	  pairpass = kFALSE;
	  break;
	}
	else {
	  pairpass = kTRUE;
	}
      }
      else {
	pairpass = kTRUE;
      }
    }
  }
  // Perform cut only for fRadiusMin:
  else {
    // Calculate radius:
    Double_t rad = fRadiusMin;

    // Calculate dPhiStar:
    double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
    double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
    Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
    dphistar = TVector2::Phi_mpi_pi(dphistar); // returns phi angle in the interval [-PI,PI)

    // Check if pair parameters meet the requirements:
    if(deta > fDEtaRangeLow && deta < fDEtaRangeUp) {
      if(dphistar > fDPhiStarRangeLow && dphistar < fDPhiStarRangeUp) {
	pairpass = kFALSE;
      }
      else {
	pairpass = kTRUE;
      }
    }
    else {
      pairpass = kTRUE;
    }
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
AliFemtoString AliFemtoPairCutRadialDistanceAsymmetric::Report() {
  // Prepare a report from the execution
  string stemp = "AliFemtoRadialDistanceAsymmetric Pair Cut - remove shared and split pairs and pairs with small separation at the specified radius\n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with separation more that %f and less than %f", fDPhiStarRangeLow, fDPhiStarRangeUp);
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;
}

//__________________
TList *AliFemtoPairCutRadialDistanceAsymmetric::ListSettings() {
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoPairCutAntiGamma::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutRadialDistanceAsymmetric.phistarrange=(%f,%f)", fDPhiStarRangeLow, fDPhiStarRangeUp);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

//__________________
void AliFemtoPairCutRadialDistanceAsymmetric::SetDPhiStarRangeLow(double minphistar) {
  fDPhiStarRangeLow = minphistar;
}

//__________________
void AliFemtoPairCutRadialDistanceAsymmetric::SetDPhiStarRangeUp(double maxphistar) {
  fDPhiStarRangeUp = maxphistar;
}

//__________________
void AliFemtoPairCutRadialDistanceAsymmetric::SetDEtaRangeLow(double mineta) {
  fDEtaRangeLow = mineta;
}

//__________________
void AliFemtoPairCutRadialDistanceAsymmetric::SetDEtaRangeUp(double maxeta) {
  fDEtaRangeUp = maxeta;
}

//__________________
void AliFemtoPairCutRadialDistanceAsymmetric::SetMagneticFieldSign(int magsign) {
  if(magsign>1) 
    fMagSign = 1;
  else if(magsign<1) 
    fMagSign = -1;
  else 
    fMagSign = magsign;
}

//__________________
void AliFemtoPairCutRadialDistanceAsymmetric::SetRadiusMin(double radmin) {
  fRadiusMin = radmin;
}

//__________________
void AliFemtoPairCutRadialDistanceAsymmetric::SetRadiusMax(double radmax) {
  fRadiusMax = radmax;
}

//__________________
void AliFemtoPairCutRadialDistanceAsymmetric::SetCalculateRadiusRange(bool calculaterange) {
  fCalculateRadiusRange = calculaterange;
}
