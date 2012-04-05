////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleMomRes - the cut monitor for particles to study  //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticleMomRes_H
#define AliFemtoCutMonitorParticleMomRes_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
class TH1D;
class TH2D;
class TH3D;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorParticleMomRes : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticleMomRes();
  AliFemtoCutMonitorParticleMomRes(const char *aName);
  AliFemtoCutMonitorParticleMomRes(const AliFemtoCutMonitorParticleMomRes &aCut);
  virtual ~AliFemtoCutMonitorParticleMomRes();

  AliFemtoCutMonitorParticleMomRes& operator=(const AliFemtoCutMonitorParticleMomRes& aCut);

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
  TH3D *fMomRes3D;    // 3D momentum resolution
  TH2D *fMomResXvsP;  // X resolution vs momentum
  TH2D *fMomResYvsP;  // Y resolution vs momentum
  TH2D *fMomResZvsP;  // Z resolution vs momentum
  TH2D *fImpactXY;    // XY impact parameter
  TH2D *fImpactZ;     // Z impact parameter
  TH2D *fSigma;       // Sigma to vertex vs momentum
};

#endif
