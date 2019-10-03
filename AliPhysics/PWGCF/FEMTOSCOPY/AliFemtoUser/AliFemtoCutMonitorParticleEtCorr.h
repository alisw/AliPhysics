////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleEtCorr - the cut monitor for particles           //
// which saves particles' et histogram and makes the bin-by-bin correlation   //
//                                                                            //
// Author: Adam.Kisiel@cern.ch                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticleEtCorr_hh
#define AliFemtoCutMonitorParticleEtCorr_hh

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

class AliFemtoCutMonitorParticleEtCorr : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticleEtCorr();
  AliFemtoCutMonitorParticleEtCorr(const char *aName, int aPhiBins);
  AliFemtoCutMonitorParticleEtCorr(const AliFemtoCutMonitorParticleEtCorr &aCut);
  virtual ~AliFemtoCutMonitorParticleEtCorr();

  AliFemtoCutMonitorParticleEtCorr& operator=(const AliFemtoCutMonitorParticleEtCorr& aCut);

  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);

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
  Double_t fPtSumEvent[200];      // A table where Pt sum per event is stored
  Double_t fMultSumEvent[200];    // A table where mult sum per event is stored
  Int_t    fPhiBins;              // Number of Phi bins
  TH1D    *fPtPerPhi;             // Histogram storing per-bin sum pt
  TH2D    *fPtCovPerPhi;          // Histogram storing per-bin covariance
  TH2D    *fPtMultPerPhi;         // Histogram storing per-bin multiplicity
  Int_t    fNEventsProcessed;     // Count processed events
};

#endif
