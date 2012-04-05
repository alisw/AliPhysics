////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorEventVertex - the cut monitor for events to study      ///
/// the distribution and error of the primary vertex                         ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorEventVertex_hh
#define AliFemtoCutMonitorEventVertex_hh

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

class AliFemtoCutMonitorEventVertex : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorEventVertex();
  AliFemtoCutMonitorEventVertex(const char *aName);
  AliFemtoCutMonitorEventVertex(const AliFemtoCutMonitorEventVertex &aCut);
  virtual ~AliFemtoCutMonitorEventVertex();

  AliFemtoCutMonitorEventVertex& operator=(const AliFemtoCutMonitorEventVertex& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent);
  virtual void Fill(const AliFemtoTrack* aTrack) {AliFemtoCutMonitor::Fill(aTrack);}
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
  TH1D *fEvVertRad;     // Vertex position in radial direction
  TH2D *fEvVertXY;      // Vertex position in XY plane
  TH1D *fEvVertSigXY;   // Sigma in XY plane
  TH1D *fEvVertZ;       // Vertex position in Z direction
  TH1D *fEvVertSigZ;    // Sigma in Z direction
};

#endif
