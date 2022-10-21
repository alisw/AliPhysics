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
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorPairBetaT::AliFemtoCutMonitorPairBetaT():
  fHistBetaT(0),
  fHistBetaTpT1(0),
  fHistBetaTpT2(0),
  fBinsBetaT(100),
  fMinBetaT(0.0),
  fMaxBetaT(1.0),
  fMassPart1(0.13957018),
  fMassPart2(0.13957018)
{
  // Default constructor
  fMinBetaT = 0.0;
  fMaxBetaT = 1.0;
  fBinsBetaT = 100;
  fHistBetaT = new TH1D("BetaT", "BetaT distribution", fBinsBetaT, 0.0, 1.0);
  fHistBetaTpT1 = new TH2D("histBetaT1pT","BetaT vs pT for part 1", 20, 0.0, 1.0, 245, 0.0, 5.0);
  fHistBetaTpT2 = new TH2D("histBetaT2pT","BetaT vs pT for part 2", 20, 0.0, 1.0, 245, 0.0, 5.0);
  fMassPart1 = 0.13957018;
  fMassPart2 = 0.13957018;
}

AliFemtoCutMonitorPairBetaT::AliFemtoCutMonitorPairBetaT(const char *aName, const int aBinsBetaT, double aMinBetaT, double aMaxBetaT, double aMassPart1, double aMassPart2):
  fHistBetaT(0),
  fHistBetaTpT1(0),
  fHistBetaTpT2(0),
  fBinsBetaT(100),
  fMinBetaT(0.0),
  fMaxBetaT(1.0),
  fMassPart1(0.13957018),
  fMassPart2(0.13957018)
{
  // Normal constructor
  fBinsBetaT = aBinsBetaT;
  fMinBetaT = aMinBetaT;
  fMaxBetaT = aMaxBetaT;
  fMassPart1 = aMassPart1;
  fMassPart2 = aMassPart2;
  char name[200];
  snprintf(name, 200, "BetaT%s", aName);
  char name1[200];
  snprintf(name1, 200, "BetaT1pT%s", aName);
  char name2[200];
  snprintf(name2, 200, "BetaT2pT%s", aName);
  fHistBetaT = new TH1D(name, "BetaT distribution", fBinsBetaT, 0.0, 1.0);
  fHistBetaTpT1 = new TH2D(name1, "BetaT vs pT for part 1", 20, 0.0, 1.0, 245, 0.0, 5.0);
  fHistBetaTpT2 = new TH2D(name2, "BetaT vs pT for part 2", 20, 0.0, 1.0, 245, 0.0, 5.0);
}

AliFemtoCutMonitorPairBetaT::AliFemtoCutMonitorPairBetaT(const AliFemtoCutMonitorPairBetaT &c):
  fHistBetaT(0),
  fHistBetaTpT1(0),
  fHistBetaTpT2(0),  
  fBinsBetaT(100),
  fMinBetaT(0.0),
  fMaxBetaT(1.0),
  fMassPart1(0.13957018),
  fMassPart2(0.13957018)
{
  // Copy constructor
  fBinsBetaT = c.fBinsBetaT;
  fMinBetaT = c.fMinBetaT;
  fMaxBetaT = c.fMaxBetaT;
  fMassPart1 = c.fMassPart1;
  fMassPart2 = c.fMassPart2;
  fHistBetaT = new TH1D(*c.fHistBetaT);
  fHistBetaTpT1 = new TH2D(*c.fHistBetaTpT1);
  fHistBetaTpT2 = new TH2D(*c.fHistBetaTpT2);
}

AliFemtoCutMonitorPairBetaT::~AliFemtoCutMonitorPairBetaT() {
  // Destructor
  delete fHistBetaT;
  delete fHistBetaTpT1;
  delete fHistBetaTpT2;
}

AliFemtoCutMonitorPairBetaT& AliFemtoCutMonitorPairBetaT::operator=(const AliFemtoCutMonitorPairBetaT& c) {
  // Assignment operator
  if (this == &c) 
    return *this;

  if (fHistBetaT) delete fHistBetaT;
  if (fHistBetaTpT1) delete fHistBetaTpT1;
  if (fHistBetaTpT2) delete fHistBetaTpT2;
  fHistBetaT = new TH1D(*c.fHistBetaT);
  fHistBetaTpT1 = new TH2D(*c.fHistBetaTpT1);
  fHistBetaTpT2 = new TH2D(*c.fHistBetaTpT2);
  fBinsBetaT = c.fBinsBetaT;
  fMinBetaT = c.fMinBetaT;
  fMaxBetaT = c.fMaxBetaT;
  fMassPart1 = c.fMassPart1;
  fMassPart2 = c.fMassPart2;
  
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
  // Calculate energies of particles:
  double pz1 = aPair->Track1()->Track()->P().z();
  double pz2 = aPair->Track2()->Track()->P().z();
  double pzpair = pz1 + pz2;
  double p1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);
  double p2 = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2);
  double m1 = fMassPart1;
  double m2 = fMassPart2;
  double e1 = TMath::Sqrt(p1*p1 + m1*m1);
  double e2 = TMath::Sqrt(p2*p2 + m2*m2);
  // Calculate transverse mass of the pair:
  double mInvpair_2 = m1*m1 + m2*m2 + 2*(e1*e2 - px1*px2 - py1*py2 - pz1*pz2);
  double mTpair = TMath::Sqrt(mInvpair_2 + pTpair*pTpair);
  // Calculate betaT:
  double betaT = pTpair / mTpair;
  
  // Fill in the monitor histograms with the values from the current pair
  fHistBetaT->Fill(betaT);
  fHistBetaTpT1->Fill(betaT,TMath::Sqrt(px1*px1 + py1*py1));
  fHistBetaTpT2->Fill(betaT,TMath::Sqrt(px2*px2 + py2*py2));
}

void AliFemtoCutMonitorPairBetaT::Write() {
  // Write out the relevant histograms
  fHistBetaT->Write();
  fHistBetaTpT1->Write();
  fHistBetaTpT2->Write();
}

TList *AliFemtoCutMonitorPairBetaT::GetOutputList() {
  TList *tOutputList = new TList();
  tOutputList->Add(fHistBetaT);
  tOutputList->Add(fHistBetaTpT1);
  tOutputList->Add(fHistBetaTpT2);

  return tOutputList;
}
