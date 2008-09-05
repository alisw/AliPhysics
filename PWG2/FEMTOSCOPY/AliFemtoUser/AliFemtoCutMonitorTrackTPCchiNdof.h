////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorTrackTPCchiNdof - the cut monitor for tracks to study  ///
/// the number of TPC Clusters distribution.                                 ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorTrackTPCchiNdof_hh
#define AliFemtoCutMonitorTrackTPCchiNdof_hh

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

class AliFemtoCutMonitorTrackTPCchiNdof : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorTrackTPCchiNdof();
  AliFemtoCutMonitorTrackTPCchiNdof(const char *aName);
  AliFemtoCutMonitorTrackTPCchiNdof(const AliFemtoCutMonitorTrackTPCchiNdof &aCut);
  virtual ~AliFemtoCutMonitorTrackTPCchiNdof();

  AliFemtoCutMonitorTrackTPCchiNdof& operator=(const AliFemtoCutMonitorTrackTPCchiNdof& aCut);

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
  TH1D *fTrTPCchiNdof;    // TPC track TPC clusters distribution
};

#endif
