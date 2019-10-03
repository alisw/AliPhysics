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

#include "AliFemtoPairCutRadialDistanceKK.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutRadialDistanceKK)
#endif

//__________________
AliFemtoPairCutRadialDistanceKK::AliFemtoPairCutRadialDistanceKK():
  AliFemtoPairCutAntiGamma(),
  fDPhiStarMin(0),
  fEtaMin(0),
  fMinRad(0.8),
  fMagSign(1)
{
}
//__________________
AliFemtoPairCutRadialDistanceKK::AliFemtoPairCutRadialDistanceKK(const AliFemtoPairCutRadialDistanceKK& c) : 
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
AliFemtoPairCutRadialDistanceKK::~AliFemtoPairCutRadialDistanceKK(){
  /* no-op */
}
AliFemtoPairCutRadialDistanceKK& AliFemtoPairCutRadialDistanceKK::operator=(const AliFemtoPairCutRadialDistanceKK& c)
{
  if (this != &c) {
    fDPhiStarMin = c.fDPhiStarMin;
    fEtaMin = c.fEtaMin;
    fMinRad = c.fMinRad;
    fMagSign = c.fMagSign;
  }

  return *this;
}
//__________________
bool AliFemtoPairCutRadialDistanceKK::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC entrance separation and
  // quality and sharity
  //  bool temp = true;
  
    double PI = 3.14159265358979312;
//    double pit = 6.28318530717958623;


 Bool_t pass5 = kTRUE;
  
  double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double chg1 = pair->Track1()->Track()->Charge();
  double chg2 = pair->Track2()->Track()->Charge();
  double ptv1 = pair->Track1()->Track()->Pt();
  double ptv2 = pair->Track2()->Track()->Pt();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();


    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());


    AliAODEvent *fAOD;

    if (!aodH) {
        // AliWarning("Could not get AODInputHandler");
        return false;
    }
    else {

        fAOD = aodH->GetEvent();
    }
    
     Double_t Bfield = fAOD->GetMagneticField();
  

    if(fabs(eta1-eta2) < fEtaMin) {
  
     // propagate through B field to r=1m
    Double_t phistar1 = phi1 - asin(chg1*(0.1*Bfield)*0.15/ptv1);// 0.15 for D=1m
    if(phistar1 > 2*PI) phistar1 -= 2*PI;
    if(phistar1 < 0) phistar1 += 2*PI;
    Double_t phistar2 = phi2 - asin(chg2*(0.1*Bfield)*0.15/ptv2);// 0.15 for D=1m 
    if(phistar2 > 2*PI) phistar2 -= 2*PI;
    if(phistar2 < 0) phistar2 += 2*PI;

    Double_t deltaphi = phistar1 - phistar2;
    if(deltaphi > PI) deltaphi -= 2*PI;
    if(deltaphi < -PI) deltaphi += 2*PI;
    deltaphi = fabs(deltaphi);

   
     if(deltaphi < fDPhiStarMin){ 
    //  cout<<"---1--- phi1 = "<<phi1<<" phi2 = "<<phi2<<endl;
    //  cout<<" eta1 = "<<eta1<<" eta2 = "<<eta2<<endl;
    //  cout<<" ptv1 = "<<ptv1<<" ptv2 = "<<ptv2<<endl;
    //   cout<<"----1-- "<< deltaphi<<" deltaeta"<< fabs(eta1-eta2)<<endl;
       pass5 = kFALSE;
//	break;
     }
     else
     {
    // propagate through B field to r=1.6m
    phistar1 = phi1 - asin(chg1*(0.1*Bfield)*0.24/ptv1);// mine. 0.24 for D=1.6m
    if(phistar1 > 2*PI) phistar1 -= 2*PI;
    if(phistar1 < 0) phistar1 += 2*PI;
    phistar2 = phi2 - asin(chg2*(0.1*Bfield)*0.24/ptv2);// mine. 0.24 for D=1.6m 
    if(phistar2 > 2*PI) phistar2 -= 2*PI;
    if(phistar2 < 0) phistar2 += 2*PI;

    deltaphi = phistar1 - phistar2;
    if(deltaphi > PI) deltaphi -= 2*PI;
    if(deltaphi < -PI) deltaphi += 2*PI;
    deltaphi = fabs(deltaphi);

    if(deltaphi < fDPhiStarMin){ 	
//      cout<<"---2--- phi1 = "<<phi1<<" phi2 = "<<phi2<<endl;
//      cout<<" eta1 = "<<eta1<<" eta2 = "<<eta2<<endl;
//      cout<<" ptv1 = "<<ptv1<<" ptv2 = "<<ptv2<<endl;
//      cout<<"----2-- "<< deltaphi<<" deltaeta"<< fabs(eta1-eta2)<<endl;
        pass5 = kFALSE;
//	break;
      }
      else
      {
       pass5 = kTRUE;
      }
    
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
AliFemtoString AliFemtoPairCutRadialDistanceKK::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoRadialDistance Pair Cut - remove shared and split pairs and pairs with small separation at the specified radius\n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with separation more that %f",fDPhiStarMin);
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutRadialDistanceKK::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoPairCutAntiGamma::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutRadialDistanceKK.phistarsepmin=%f", fDPhiStarMin);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutRadialDistanceKK::SetPhiStarDifferenceMinimum(double dtpc)
{
  fDPhiStarMin = dtpc;
}

void AliFemtoPairCutRadialDistanceKK::SetEtaDifferenceMinimum(double etpc) 
{
  fEtaMin = etpc;
}


void AliFemtoPairCutRadialDistanceKK::SetMinimumRadius(double minrad) 
{
  fMinRad = minrad;
}

void AliFemtoPairCutRadialDistanceKK::SetMagneticFieldSign(int magsign)
{
  if(magsign>1) fMagSign = 1;
  else if(magsign<1) fMagSign = -1;
  else fMagSign = magsign;
}
