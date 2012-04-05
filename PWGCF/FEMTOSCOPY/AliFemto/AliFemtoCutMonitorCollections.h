////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorCollections - the cut monitor for events to study        ///
/// the multiplicity distribution of events                                  ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorCollections_hh
#define AliFemtoCutMonitorCollections_hh

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; 
class TH1D;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorCollections : public AliFemtoCutMonitor{
  
 public:
  AliFemtoCutMonitorCollections();
  AliFemtoCutMonitorCollections(const char *aName);
  AliFemtoCutMonitorCollections(const AliFemtoCutMonitorCollections &aCut);
  virtual ~AliFemtoCutMonitorCollections();

  AliFemtoCutMonitorCollections& operator=(const AliFemtoCutMonitorCollections& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack) {AliFemtoCutMonitor::Fill(aTrack);}
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection){AliFemtoCutMonitor::Fill(aEvent, aCollection);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2);

  void Write();

  virtual TList *GetOutputList();

 private:
  TH1D *fCollection1Mult;     // Collection 1 multiplicity distribution
  TH1D *fCollection2Mult; // Collection 2 multiplicity distribution


};

#endif
