////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticlePID - the cut monitor for particles to study   ///
/// various aspects of the PID determination                                 ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticlePID_hh
#define AliFemtoCutMonitorParticlePID_hh

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
class TH1D;
class TH2D;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorParticlePID : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticlePID();
  AliFemtoCutMonitorParticlePID(const char *aName);
  AliFemtoCutMonitorParticlePID(const AliFemtoCutMonitorParticlePID &aCut);
  virtual ~AliFemtoCutMonitorParticlePID();

  AliFemtoCutMonitorParticlePID& operator=(const AliFemtoCutMonitorParticlePID& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack);
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}


  void Write();

  virtual TList *GetOutputList();

private:
  TH2D *fTPCdEdx;   // TPC dEdx information

};

#endif
