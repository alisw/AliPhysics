///
/// \file AliFemtoCutMonitorEventSphericity.h
///

#ifndef AliFemtoCutMonitorEventSphericity_hh
#define AliFemtoCutMonitorEventSphericity_hh

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

/// \class AliFemtoCutMonitorEventSphericity
/// \brief The cut monitor to study the multiplicity distribution of events.
///
class AliFemtoCutMonitorEventSphericity : public AliFemtoCutMonitor {
public:
  AliFemtoCutMonitorEventSphericity();
  AliFemtoCutMonitorEventSphericity(const char *aName, int nBins=5000, double multMax=5000.5);
  AliFemtoCutMonitorEventSphericity(const AliFemtoCutMonitorEventSphericity &aCut);
  virtual ~AliFemtoCutMonitorEventSphericity();

  AliFemtoCutMonitorEventSphericity& operator=(const AliFemtoCutMonitorEventSphericity& aCut);

  virtual AliFemtoString Report();

  virtual void Fill(const AliFemtoEvent* aEvent);
  virtual void Fill(const AliFemtoTrack* aTrack) {AliFemtoCutMonitor::Fill(aTrack);}
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoXi* aXi) {AliFemtoCutMonitor::Fill(aXi);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aEvent, aCollection);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2) {AliFemtoCutMonitor::Fill(aCollection1, aCollection2);}

  //  void SetReadMC(Bool_t mc);
  // void AdditionalMultHistsOn(Bool_t addhists);
  void Write();

  virtual TList *GetOutputList();

private:
  TH1D *fEvSpher;      ///< Sphericity distribution
  TH2D *fEvSpherMult;  ///< Sphericity versus event multiplicity distribution
  TH2D *fMultSumPt;   ///< Event total pT vs. multiplicity

};

#endif
