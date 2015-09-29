//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// AliFemtoCutMonitorPairBetaT - the cut monitor for particles to study     //
// the betaT bins                                                           //
//                                                                          //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOCUTMONITORPAIRBETAT_H
#define ALIFEMTOCUTMONITORPAIRBETAT_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair;
class TH1D;
class TList;
//#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorPairBetaT : public AliFemtoCutMonitor {

 public:
  AliFemtoCutMonitorPairBetaT();
  AliFemtoCutMonitorPairBetaT(const char *aName, const int aBinsBetaT, double aMinBetaT, double aMaxBetaT, double aMassPart1, double aMassParrt2);
  AliFemtoCutMonitorPairBetaT(const AliFemtoCutMonitorPairBetaT& c);
  virtual ~AliFemtoCutMonitorPairBetaT();

  AliFemtoCutMonitorPairBetaT& operator=(const AliFemtoCutMonitorPairBetaT& c);

  virtual AliFemtoString Report();
  void Write();
  virtual TList *GetOutputList();
  virtual void Fill(const AliFemtoEvent* aEvent) { AliFemtoCutMonitor::Fill(aEvent); }
  virtual void Fill(const AliFemtoTrack* aTrack) { AliFemtoCutMonitor::Fill(aTrack); }
  virtual void Fill(const AliFemtoV0* aV0) { AliFemtoCutMonitor::Fill(aV0); }
  virtual void Fill(const AliFemtoXi* aXi) {AliFemtoCutMonitor::Fill(aXi);}
  virtual void Fill(const AliFemtoKink* aKink) { AliFemtoCutMonitor::Fill(aKink); }
  virtual void Fill(const AliFemtoPair* aPair);
  virtual void Fill(const AliFemtoParticleCollection* aCollection) { AliFemtoCutMonitor::Fill(aCollection); }
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection) { AliFemtoCutMonitor::Fill(aEvent, aCollection); }
  virtual void Fill(const AliFemtoParticleCollection* aCollection1, const AliFemtoParticleCollection* aCollection2) { AliFemtoCutMonitor::Fill(aCollection1, aCollection2); }

 private:
  TH1D *fHistBetaT;      // BetaT plot
  int fBinsBetaT;        // Number of bins in betaT plot
  double fMinBetaT;      // Minimum betaT
  double fMaxBetaT;      // Maximum betaT
  double fMassPart1;     // Mass of the first particle in pair [GeV]
  double fMassPart2;     // Mass of the second particle in pair [GeV]
};

#endif
