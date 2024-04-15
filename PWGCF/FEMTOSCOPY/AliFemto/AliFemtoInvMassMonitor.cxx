////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoInvMassMonitor - the cut monitor for invariant mass of pairs       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoInvMassMonitor.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoInvMassMonitor::AliFemtoInvMassMonitor():
  fInvMass(0),
  fM1(0),
  fM2(0)
{
  // Default constructor
  fInvMass = new TH1D("fInvMass", "Invariant mass for kaons",80, 0.99, 1.07);
}

AliFemtoInvMassMonitor::AliFemtoInvMassMonitor(const char *aName,double m1, double m2):
 fInvMass(0),
 fM1(m1),
 fM2(m2)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "fInvMass%s", aName);
  fInvMass = new TH1D(name, "Invariant mass for kaons", 80, 0.99, 1.07);
}

AliFemtoInvMassMonitor::AliFemtoInvMassMonitor(const AliFemtoInvMassMonitor &aCut):
 fInvMass(0),
 fM1(0),
 fM2(0)
{
  // copy constructor
  fInvMass = new TH1D(*aCut.fInvMass);
  fM1 = aCut.fM1;
  fM2 = aCut.fM2;
}

AliFemtoInvMassMonitor::~AliFemtoInvMassMonitor()
{
  // Destructor
  delete fInvMass;
  
}

AliFemtoInvMassMonitor& AliFemtoInvMassMonitor::operator=(const AliFemtoInvMassMonitor& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fInvMass) delete fInvMass;
  fInvMass = new TH1D(*aCut.fInvMass);
  fM1 = aCut.fM1;
  fM2 = aCut.fM2;
 
  return *this;
}

AliFemtoString AliFemtoInvMassMonitor::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoInvMassMonitor report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoInvMassMonitor::Fill(const AliFemtoPair* pair)
{

  double phi1 = pair->Track1()->FourMomentum().Phi();
  double phi2 = pair->Track2()->FourMomentum().Phi();
  double eta1 = pair->Track1()->FourMomentum().PseudoRapidity();
  double eta2 = pair->Track2()->FourMomentum().PseudoRapidity();

  double px1 = pair->Track1()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  //double pz1 = aPair->Track1()->Track()->P().z();

  double px2 = pair->Track2()->Track()->P().x();
  double py2 = pair->Track2()->Track()->P().y();
  //double pz2 = aPair->Track2()->Track()->P().z();

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

  fInvMass->Fill(minv);

}

void AliFemtoInvMassMonitor::Write()
{
  // Write out the relevant histograms
  fInvMass->Write();
}

TList *AliFemtoInvMassMonitor::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fInvMass);
  

  return tOutputList;
}
