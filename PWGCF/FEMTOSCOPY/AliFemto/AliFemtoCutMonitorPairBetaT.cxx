//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// AliFemtoCutMonitorPairBetaT - the cut monitor for particles to study     //
// the betaT bins                                                           //
//                                                                          //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorPairBetaT.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoPair.h"
#include <TH1D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorPairBetaT::AliFemtoCutMonitorPairBetaT():
  fHistBetaT(0),
  fBinsBetaT(100),
  fMinBetaT(0.0),
  fMaxBetaT(1.0)
{
  // Default constructor
  fMinBetaT = 0.0;
  fMaxBetaT = 1.0;
  fBinsBetaT = 100;
  fHistBetaT = new TH1D("BetaT", "BetaT distribution", fBinsBetaT, 0.0, 1.0);
}

AliFemtoCutMonitorPairBetaT::AliFemtoCutMonitorPairBetaT(const char *aName, const int aBinsBetaT, double aMinBetaT, double aMaxBetaT):
  fHistBetaT(0),
  fBinsBetaT(100),
  fMinBetaT(0.0),
  fMaxBetaT(1.0)
{
  // Normal constructor
  fBinsBetaT = aBinsBetaT;
  fMinBetaT = aMinBetaT;
  fMaxBetaT = aMaxBetaT;
  char name[200];
  snprintf(name, 200, "BetaT%s", aName);
  fHistBetaT = new TH1D(name, "BetaT distribution", 100, 0.0, 1.0);
}

AliFemtoCutMonitorPairBetaT::AliFemtoCutMonitorPairBetaT(const AliFemtoCutMonitorPairBetaT &c):
  fHistBetaT(0),
  fBinsBetaT(100),
  fMinBetaT(0.0),
  fMaxBetaT(1.0)
{
  // Copy constructor
  fBinsBetaT = c.fBinsBetaT;
  fMinBetaT = c.fMinBetaT;
  fMaxBetaT = c.fMaxBetaT;
  fHistBetaT = new TH1D(*c.fHistBetaT);
}

AliFemtoCutMonitorPairBetaT::~AliFemtoCutMonitorPairBetaT() {
  // Destructor
  delete fHistBetaT;
}

AliFemtoCutMonitorPairBetaT& AliFemtoCutMonitorPairBetaT::operator=(const AliFemtoCutMonitorPairBetaT& c) {
  // Assignment operator
  if (this == &c) 
    return *this;

  if (fHistBetaT) delete fHistBetaT;
  fHistBetaT = new TH1D(*c.fHistBetaT);
  fBinsBetaT = c.fBinsBetaT;
  fMinBetaT = c.fMinBetaT;
  fMaxBetaT = c.fMaxBetaT;
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorPairBetaT::Report() { 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorPairBetaT report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorPairBetaT::Fill(const AliFemtoPair* aPair) {
  // Calculate transverse momentum of the pair:
  double px1 = aPair->Track1()->Track()->P().x();
  double px2 = aPair->Track2()->Track()->P().x();
  double py1 = aPair->Track1()->Track()->P().y();
  double py2 = aPair->Track2()->Track()->P().y();
  double pxpair = px1 + px2;
  double pypair = py1 + py2;
  double pTpair = TMath::Sqrt(pxpair*pxpair + pypair*pypair);
  // Calculate energy of the pair:
  double pz1 = aPair->Track1()->Track()->P().z();
  double pz2 = aPair->Track2()->Track()->P().z();
  double pzpair = pz1 + pz2;
  double p1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);
  double p2 = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2);
  double m1 = aPair->Track1()->Track()->GetMass();
  double m2 = aPair->Track2()->Track()->GetMass();
  double e1 = TMath::Sqrt(p1*p1 + m1*m1);
  double e2 = TMath::Sqrt(p2*p2 + m2*m2);
  double epair = e1 + e2;
  // Calculate transverse mass of the pair:
  double mTpair = TMath::Sqrt(epair*epair - pzpair*pzpair);
  // Calculate betaT:
  double betaT = pTpair / mTpair;
  
  // Fill in the monitor histograms with the values from the current pair
  fHistBetaT->Fill(betaT);
}

void AliFemtoCutMonitorPairBetaT::Write() {
  // Write out the relevant histograms
  fHistBetaT->Write();
}

TList *AliFemtoCutMonitorPairBetaT::GetOutputList() {
  TList *tOutputList = new TList();
  tOutputList->Add(fHistBetaT);

  return tOutputList;
}
