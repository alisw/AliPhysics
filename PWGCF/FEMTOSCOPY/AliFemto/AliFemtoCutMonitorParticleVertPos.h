////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticleVertPos - the cut monitor for particles to study  ///
/// the difference between reconstructed and true momentum    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticleVertPos_hh
#define AliFemtoCutMonitorParticleVertPos_hh

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

class AliFemtoCutMonitorParticleVertPos : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticleVertPos();
  AliFemtoCutMonitorParticleVertPos(const char *aName);
  AliFemtoCutMonitorParticleVertPos(const AliFemtoCutMonitorParticleVertPos &aCut);
  virtual ~AliFemtoCutMonitorParticleVertPos();

  AliFemtoCutMonitorParticleVertPos& operator=(const AliFemtoCutMonitorParticleVertPos& aCut);

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
  TH2D *fVertPos;    // Vertex position x vs. y monitor
  TH2D *fEtaZ;       // Vertex z position vs. eta monitor
  TH1D *fRadPos;     // Radial position close to vertex
  TH1D *fEmPointX;   // Emission point - x
  TH1D *fEmPointY;   // Emission point - y
  TH1D *fEmPointZ;   // Emission point - z
  TH1D *fEmPointT;   // Emission point - t
  
};

#endif
