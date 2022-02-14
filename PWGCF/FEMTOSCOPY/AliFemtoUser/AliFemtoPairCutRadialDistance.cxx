/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistance - a pair cut which checks                     //
// for some pair qualities that attempt to identify slit/doubly                //
// reconstructed tracks and also selects pairs based on their separation       //
// at the entrance to the TPC                                                  //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////
/********************************************************************************
 *
 * Author: Johanna Gramling, University of Heidelberg, jgramlin@cern.ch
 *         Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch
 *         Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch
 *
 ********************************************************************************/

#include "AliFemtoPairCutRadialDistance.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutRadialDistance)
#endif

//__________________
AliFemtoPairCutRadialDistance::AliFemtoPairCutRadialDistance():
AliFemtoPairCutAntiGamma(),
  fDPhiStarMin(0),
  fEtaMin(0),
  fMinRad(0.8),
  fMaxRad(2.5),
  fMagSign(1),
  fMagFieldVal(0.5),
  fPhistarmin(kTRUE)
{
}
//__________________
AliFemtoPairCutRadialDistance::AliFemtoPairCutRadialDistance(const AliFemtoPairCutRadialDistance& c) :
  AliFemtoPairCutAntiGamma(c),
  fDPhiStarMin(0),
  fEtaMin(0),
  fMinRad(0.8),
  fMaxRad(2.5),
  fMagSign(1),
  fMagFieldVal(0.5),
  fPhistarmin(kTRUE)
{
  fDPhiStarMin = c.fDPhiStarMin;
  fEtaMin = c.fEtaMin;
  fMinRad = c.fMinRad;
  fMaxRad = c.fMaxRad;
  fMagSign = c.fMagSign;
  fMagFieldVal = c.fMagFieldVal;
  fPhistarmin = c.fPhistarmin;
}

//__________________
AliFemtoPairCutRadialDistance::~AliFemtoPairCutRadialDistance(){
  /* no-op */
}
AliFemtoPairCutRadialDistance& AliFemtoPairCutRadialDistance::operator=(const AliFemtoPairCutRadialDistance& c)
{
  if (this != &c) {
    fDPhiStarMin = c.fDPhiStarMin;
    fEtaMin = c.fEtaMin;
    fMinRad = c.fMinRad;
    fMaxRad = c.fMaxRad;
    fMagSign = c.fMagSign;
    fMagFieldVal = c.fMagFieldVal;
    fPhistarmin = c.fPhistarmin;

  }

  return *this;
}
//__________________
bool AliFemtoPairCutRadialDistance::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC entrance separation and
  // quality and sharity
  //  bool temp = true;

//    double pih = 3.14159265358979312;
//    double pit = 6.28318530717958623;

  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double chg1 = pair->Track1()->Track()->Charge();
  double chg2 = pair->Track2()->Track()->Charge();
  double ptv1 = pair->Track1()->Track()->Pt();
  double ptv2 = pair->Track2()->Track()->Pt();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();


  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  Double_t magsign = 0.0;

  if (!aodH) {
    // AliWarning("Could not get AODInputHandler");
    return false;
  }
  else {
    AliAODEvent *fAOD;
    fAOD = aodH->GetEvent();
    //cout<<fAOD<<endl;
    magsign = fAOD->GetMagneticField();
  }

  if (magsign > 1)
    fMagSign = 1;
  else if ( magsign < 1)
    fMagSign = -1;
  else
    fMagSign = magsign;

  //cout << "mag sign = " << magsign << endl;

  Double_t rad;
  Bool_t pass5 = kTRUE;

  rad = fMinRad;

  if (fPhistarmin) {
    for (rad = fMinRad; rad < fMaxRad; rad += 0.01) {
      Double_t dps = (phi2-phi1+(TMath::ASin(-0.15*fMagFieldVal*chg2*fMagSign*rad/ptv2))-(TMath::ASin(-0.15*fMagFieldVal*chg1*fMagSign*rad/ptv1)));
      dps = TVector2::Phi_mpi_pi(dps);
      Double_t etad = eta2 - eta1;
      if (fabs(etad)<fEtaMin && fabs(dps)<fDPhiStarMin) {
        // cout << "5% cut is not passed - returning" << endl;
        pass5 = kFALSE;
        break;
      }
    }
  }
  else {

    double afsi0b = 0.15*fMagFieldVal*chg1*fMagSign*rad/ptv1;
    double afsi1b = 0.15*fMagFieldVal*chg2*fMagSign*rad/ptv2;

    if (fabs(afsi0b) >=1.) return kTRUE;
    if (fabs(afsi1b) >=1.) return kTRUE;

    Double_t dps =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
    dps = TVector2::Phi_mpi_pi(dps);

    Double_t etad = eta2 - eta1;
    if (fabs(etad)<fEtaMin && fabs(dps)<fDPhiStarMin) {
      pass5 = kFALSE;
    }

  }

  if (pass5) {
    pass5 = AliFemtoPairCutAntiGamma::Pass(pair);
  }
  else {
    fNPairsFailed++;
  }

  return pass5;
}
//__________________
AliFemtoString AliFemtoPairCutRadialDistance::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoRadialDistance Pair Cut - remove shared and split pairs and pairs with small separation at the specified radius\n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with separation more that %f",fDPhiStarMin);
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutRadialDistance::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoPairCutAntiGamma::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutRadialDistance.phistarsepmin=%f", fDPhiStarMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutRadialDistance::SetPhiStarDifferenceMinimum(double dtpc)
{
  fDPhiStarMin = dtpc;
}

void AliFemtoPairCutRadialDistance::SetEtaDifferenceMinimum(double etpc)
{
  fEtaMin = etpc;
}


void AliFemtoPairCutRadialDistance::SetMinimumRadius(double minrad)
{
  fMinRad = minrad;
}

void AliFemtoPairCutRadialDistance::SetMaximumRadius(double maxrad)
{
  fMaxRad = maxrad;
}

void AliFemtoPairCutRadialDistance::SetMagneticFieldSign(int magsign)
{
  if(magsign>1) fMagSign = 1;
  else if(magsign<1) fMagSign = -1;
  else fMagSign = magsign;
}

void AliFemtoPairCutRadialDistance::SetMagneticFieldValue(double magval)
{
  fMagFieldVal = magval;
}

void AliFemtoPairCutRadialDistance::SetPhiStarMin(Bool_t phistarmin)
{
  fPhistarmin = phistarmin;
}
