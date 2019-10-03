/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoPairCutResonances - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliFemtoPairCutResonances.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutResonances)
#endif

//__________________
AliFemtoPairCutResonances::AliFemtoPairCutResonances():
   AliFemtoShareQualityPairCut(),
   fMaxEEMinv(0.0),
   fMaxDTheta(0.0),
   fDataType(kAOD),
   fSwitchPassFail(0)
{
}
//__________________
AliFemtoPairCutResonances::AliFemtoPairCutResonances(const AliFemtoPairCutResonances& c) : 
  AliFemtoShareQualityPairCut(c),
  fMaxEEMinv(0.0),
  fMaxDTheta(0.0),
  fDataType(kAOD),
  fSwitchPassFail(0)
{
  fMaxEEMinv = c.fMaxEEMinv;
  fMaxDTheta = c.fMaxDTheta;
  fDataType = c.fDataType;
  fSwitchPassFail=c.fSwitchPassFail;
}

AliFemtoPairCutResonances& AliFemtoPairCutResonances::operator=(const AliFemtoPairCutResonances& c)
{
  if (this != &c) {
    fMaxEEMinv = c.fMaxEEMinv;
    fMaxDTheta = c.fMaxDTheta;
    fDataType = c.fDataType;
    fSwitchPassFail=c.fSwitchPassFail;
  }

  return *this;

}
//__________________
AliFemtoPairCutResonances::~AliFemtoPairCutResonances(){
}
//__________________
bool AliFemtoPairCutResonances::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC entrance separation and
  // quality and sharity
  bool temp = true;

  if(fDataType==kKine)
    return true;

  double me = 0.000511;
  double mPi = 0.13957018;
  double mp = 0.938272046;

  double mK0min = 0.00049;
  double mK0max = 0.00051;
  //double mK0 = 0.000497614;
  double mRhomin = 0.000765;
  double mRhomax = 0.000785;
  //double mRho = 0.00077526;
  double mLmin = 1.095;
  double mLmax = 1.135;
  //double mL = 1.115683;

  if ((pair->Track1()->Track()->Charge() * pair->Track2()->Track()->Charge()) < 0.0) {
    // double theta1 = pair->Track1()->Track()->P().Theta();
    // double theta2 = pair->Track2()->Track()->P().Theta();
    // double dtheta = TMath::Abs(theta1 - theta2);

    // check on ee pairs (gamma)
    double e1 = TMath::Sqrt(me*me + pair->Track1()->Track()->P().Mag2());
    double e2 = TMath::Sqrt(me*me + pair->Track2()->Track()->P().Mag2());
    double minvGamma = 2*me*me + 2*(e1*e2 -
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			       pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());
    if ( minvGamma < fMaxEEMinv )
       temp = false;
    //check on resonances
    double pi1 =  TMath::Sqrt(mPi*mPi + pair->Track1()->Track()->P().Mag2());
    double pi2 =  TMath::Sqrt(mPi*mPi + pair->Track2()->Track()->P().Mag2());
    double p1 =  TMath::Sqrt(mp*mp + pair->Track1()->Track()->P().Mag2());
    double p2 =  TMath::Sqrt(mp*mp + pair->Track2()->Track()->P().Mag2());
    //check on K0 and Rho
    double minv2pi = 2*mPi*mPi + 2*(pi1*pi2 -
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			       pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());
    if ( ((minv2pi>mK0min && minv2pi<mK0max) || (minv2pi>mRhomin && minv2pi<mRhomax)) )
       temp = false;
    //check on L0
    double minvpPi = 2*mp*mPi + 2*(p1*pi2 -
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			       pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());
    double minvPip = 2*mPi*mp + 2*(pi1*p2 -
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			       pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());
    if( ((minvpPi>mLmin) && (minvpPi<mLmax)) || ((minvPip>mLmin) && (minvPip<mLmax)) )
       temp = false;
  }
  if (fSwitchPassFail) // choose only resonances
  {
     if (!temp) {
        temp = AliFemtoShareQualityPairCut::Pass(pair);
        if (temp) {fNPairsPassed++;}
        else fNPairsFailed++;
        return temp;
     }
     else
     {
        fNPairsFailed++;
        return false;
     }
  }
  else // cut resonances
  {
     if (temp) {
        temp = AliFemtoShareQualityPairCut::Pass(pair);
        if (temp) {fNPairsPassed++;}
        else fNPairsFailed++;
        return temp;
     }
     else
     {
        fNPairsFailed++;
        return false;
     }
  }
}
//__________________
AliFemtoString AliFemtoPairCutResonances::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoPairCutResonances Pair Cut - remove pairs possibly coming from Gamma conversions\n";  
  char ctemp[100];
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutResonances::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoShareQualityPairCut::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutResonances.maxeeminv=%f", fMaxEEMinv);
  snprintf(buf, 200, "AliFemtoPairCutResonances.maxdtheta=%f", fMaxDTheta);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutResonances::SetMaxEEMinv(Double_t maxeeminv)
{
  fMaxEEMinv = maxeeminv;
}

void AliFemtoPairCutResonances::SetMaxThetaDiff(Double_t maxdtheta)
{
  fMaxDTheta = maxdtheta;
}

void AliFemtoPairCutResonances::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}

void AliFemtoPairCutResonances::SetChooseResonances(bool onlyResonances)
{
  fSwitchPassFail = onlyResonances;
}
