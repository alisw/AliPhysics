/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistanceLM - a pair cut which checks                   //
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
 *         Jorge Mercado, University of Heidelberg, jmercado@cern.ch
 *
 ********************************************************************************/

#include "AliFemtoPairCutRadialDistanceLM.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutRadialDistanceLM)
#endif

//__________________
AliFemtoPairCutRadialDistanceLM::AliFemtoPairCutRadialDistanceLM():
  AliFemtoPairCutAntiGamma(),
  fDPhiStarMin(0),
  fEtaMin(0),
  fMinRad(0.8),
  fMagSign(1)
{
}
//__________________
AliFemtoPairCutRadialDistanceLM::AliFemtoPairCutRadialDistanceLM(const AliFemtoPairCutRadialDistanceLM& c) : 
  AliFemtoPairCutAntiGamma(c),
  fDPhiStarMin(0), 
  fEtaMin(0),
  fMinRad(0.8),
  fMagSign(1)
{ 
  fDPhiStarMin = c.fDPhiStarMin;
  fEtaMin = c.fEtaMin;
  fMinRad = c.fMinRad;
  fMagSign = c.fMagSign;
}

//__________________
AliFemtoPairCutRadialDistanceLM::~AliFemtoPairCutRadialDistanceLM(){
  /* no-op */
}
//__________________
bool AliFemtoPairCutRadialDistanceLM::Pass(const AliFemtoPair* pair){
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


  Double_t rad;
  Bool_t pass5 = kTRUE;
  
//   cout << "min radius: " << fMinRad << endl;
//   cout << "min deta: " << fEtaMin << endl;
//   cout << "min dphi: " << fDPhiStarMin << endl;
  
    double etad = eta2 - eta1;
    
    if (fabs(etad)<fEtaMin)
    {
      rad = fMinRad;
      for (Double_t iter=fMinRad*100; iter<251; iter+=1.0) {
	Double_t dps = (phi1-phi2+(TMath::ASin(-0.075*chg1*fMagSign*rad/ptv1))-(TMath::ASin(-0.075*chg2*fMagSign*rad/ptv2)));
	if (fabs(dps)<fDPhiStarMin) {
	  //       cout << "5% cut is not passed - returning" << endl;
	  pass5 = kFALSE;
	  break;
	}
	rad+=0.01;
      }
    }
  

  if (pass5) {
    pass5 = AliFemtoPairCutAntiGamma::Pass(pair);
  }
  else
    fNPairsFailed++;

  return pass5;
}
//__________________
AliFemtoString AliFemtoPairCutRadialDistanceLM::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoRadialDistance Pair Cut - remove shared and split pairs and pairs with small separation at the specified radius\n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with separation more that %f",fDPhiStarMin);
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutRadialDistanceLM::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoPairCutAntiGamma::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutRadialDistanceLM.phistarsepmin=%f", fDPhiStarMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutRadialDistanceLM::SetPhiStarDifferenceMinimum(double dtpc)
{
  fDPhiStarMin = dtpc;
}

void AliFemtoPairCutRadialDistanceLM::SetEtaDifferenceMinimum(double etpc) 
{
  fEtaMin = etpc;
}


void AliFemtoPairCutRadialDistanceLM::SetMinimumRadius(double minrad) 
{
  fMinRad = minrad;
}

void AliFemtoPairCutRadialDistanceLM::SetMagneticFieldSign(int magsign)
{
  if(magsign>1) fMagSign = 1;
  else if(magsign<1) fMagSign = -1;
  else fMagSign = magsign;
}
