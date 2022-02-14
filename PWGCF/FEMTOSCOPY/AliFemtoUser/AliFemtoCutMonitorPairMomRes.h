#ifndef AliFemtoCutMonitorPairMomRes_hh
#define AliFemtoCutMonitorPairMomRes_hh

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair;
class TH1D;
class TH2D;
class TList;

#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorPairMomRes : public AliFemtoCutMonitor {
public:
  AliFemtoCutMonitorPairMomRes();
  AliFemtoCutMonitorPairMomRes(const char *aName, double massPart1, double massPart2, double qmin, double qmax, int nbins);
  AliFemtoCutMonitorPairMomRes(const AliFemtoCutMonitorPairMomRes &aCut);
  virtual ~AliFemtoCutMonitorPairMomRes();

  AliFemtoCutMonitorPairMomRes& operator=(const AliFemtoCutMonitorPairMomRes& aCut);

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
  TH2D *fMomRes;
  TH2D *fMomResTrueMass;
  TH2D *fMomRes_KPpairOnly;
  TH2D *fMomResTrueMass_KPpairOnly;
  double fMassPart1;
  double fMassPart2;
  
};

#endif
