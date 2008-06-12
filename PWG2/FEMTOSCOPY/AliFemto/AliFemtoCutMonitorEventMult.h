////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorEventMult - the cut monitor for events to study        ///
/// the multiplicity distribution of events                                  ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorEventMult_hh
#define AliFemtoCutMonitorEventMult_hh

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

class AliFemtoCutMonitorEventMult : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorEventMult();
  AliFemtoCutMonitorEventMult(const char *aName);
  AliFemtoCutMonitorEventMult(const AliFemtoCutMonitorEventMult &aCut);
  virtual ~AliFemtoCutMonitorEventMult();

  AliFemtoCutMonitorEventMult& operator=(const AliFemtoCutMonitorEventMult& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent);
  virtual void Fill(const AliFemtoTrack* aTrack) {AliFemtoCutMonitor::Fill(aTrack);}
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}

  void Write();

  virtual TList *GetOutputList();

private:
  TH1D *fEvMult;    // Multiplicity distribution
};

#endif
