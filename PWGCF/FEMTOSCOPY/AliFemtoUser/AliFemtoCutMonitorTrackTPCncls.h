////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorTrackTPCncls - the cut monitor for tracks to study     ///
/// the number of TPC Clusters distribution.                                 ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorTrackTPCncls_hh
#define AliFemtoCutMonitorTrackTPCncls_hh

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

class AliFemtoCutMonitorTrackTPCncls : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorTrackTPCncls();
  AliFemtoCutMonitorTrackTPCncls(const char *aName);
  AliFemtoCutMonitorTrackTPCncls(const AliFemtoCutMonitorTrackTPCncls &aCut);
  virtual ~AliFemtoCutMonitorTrackTPCncls();

  AliFemtoCutMonitorTrackTPCncls& operator=(const AliFemtoCutMonitorTrackTPCncls& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack); 
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2) {AliFemtoCutMonitor::Fill(aCollection1, aCollection2);}

  void Write();

  virtual TList *GetOutputList();

private:
  TH1D *fTrTPCncls;    // TPC track TPC clusters distribution
};

#endif
