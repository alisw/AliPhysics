/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoPairCutMInv - a pair cut which checks if the sum of the                 //
//                       invariant mass of two particles fit within given range 	  //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//          Piotr Modzelewski, Warsaw University of Technology, pmodzele@cern.ch   //
//  	       				                                                           //
/////////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoPairCutMInv.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutMInv)
#endif

//__________________
AliFemtoPairCutMInv::AliFemtoPairCutMInv():
  AliFemtoPairCut(),
  fNPairsFailed(0),
  fNPairsPassed(0),
  fMInvMin(0),
  fMInvMax(0),
  fM1(0),
  fM2(0)
{

}
//__________________
AliFemtoPairCutMInv::AliFemtoPairCutMInv(double m1, double m2, double minvmin, double minvmax):
  AliFemtoPairCut(),
  fNPairsFailed(0),
  fNPairsPassed(0),
  fMInvMin(minvmin),
  fMInvMax(minvmax),
  fM1(m1),
  fM2(m2)
{
}
//__________________
AliFemtoPairCutMInv::AliFemtoPairCutMInv(const AliFemtoPairCutMInv& c) : 
  AliFemtoPairCut(c),
  fNPairsFailed(0),
  fNPairsPassed(0),
  fMInvMin(0),
  fMInvMax(0),
  fM1(0),
  fM2(0)
{ 
  fMInvMin = c.fMInvMin;
  fMInvMax = c.fMInvMax;
  fM1 = c.fM1;
  fM2 = c.fM2;
}
AliFemtoPairCutMInv& AliFemtoPairCutMInv::operator=(const AliFemtoPairCutMInv& c)
{
  if (this != &c) {
    fMInvMin = c.fMInvMin;
    fMInvMax = c.fMInvMax;
    fM1 = c.fM1;
    fM2 = c.fM2;
  }

  return *this;

}
// AliFemtoPairCutMInv::~AliFemtoPairCutMInv()
// {
//   cout<<Form("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",(long int) fNPairsPassed,(long int) fNPairsFailed)<<endl;
// }

//__________________
AliFemtoPairCutMInv::~AliFemtoPairCutMInv(){
  /* no-op */
}
//__________________
bool AliFemtoPairCutMInv::Pass(const AliFemtoPair* pair){

  bool temp = true;
  if (temp){
  double phi1 = pair->Track1()->FourMomentum().Phi();
  double phi2 = pair->Track2()->FourMomentum().Phi();
  double eta1 = pair->Track1()->FourMomentum().PseudoRapidity();
  double eta2 = pair->Track2()->FourMomentum().PseudoRapidity();

  double px1 = pair->Track1()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  //double pz1 = pair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->Track()->P().x();
  double py2 = pair->Track2()->Track()->P().y();
  //double pz2 = pair->Track2()->Track()->P().z();

  double pt11 = TMath::Hypot(px1, py1);
  double pt22 = TMath::Hypot(px2, py2);

  //Invariant Mass
  double Invpx1=pt11*cos(phi1);
  double Invpy1=pt11*sin(phi1);
  double Invpz1=pt11*TMath::SinH(eta1);

  double Invpx2=pt22*cos(phi2);
  double Invpy2=pt22*sin(phi2);
  double Invpz2=pt22*TMath::SinH(eta2);

  double p21 = Invpx1*Invpx1+Invpy1*Invpy1+Invpz1*Invpz1;
  double p22 = Invpx2*Invpx2+Invpy2*Invpy2+Invpz2*Invpz2;

  //double KaonMass = 0.493677;

  double e1 = TMath::Sqrt(fM1*fM1 + p21);
  double e2 = TMath::Sqrt(fM2*fM2 + p22);
   
  double minv = TMath::Sqrt(fM1*fM1 + fM2*fM2 + 2*(e1*e2 - Invpx1*Invpx2 - Invpy1*Invpy2 - Invpz1*Invpz2));
   // cout<<"Cut: "<<minv<<" masses "<<fM1<<" "<<fM2<<endl;
  if (minv > fMInvMin && minv < fMInvMax) 
    {
      temp = false; 
      // cout<<"Invariant Mass cut -> Pair rejected!"<<endl;
    }
  }

  if(temp) 
    fNPairsPassed++;
  else fNPairsFailed++;

  return temp;
  
}
//__________________
AliFemtoString AliFemtoPairCutMInv::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoPairCutMInv Pair Cut\n";  
  char ctemp[100];
  stemp += ctemp;
  snprintf(ctemp,100,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",(long int) fNPairsPassed,(long int) fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutMInv::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutMInv.fMInvMin=%f", fMInvMin);
  snprintf(buf, 200, "AliFemtoPairCutMInv.fMInvMax=%f", fMInvMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

