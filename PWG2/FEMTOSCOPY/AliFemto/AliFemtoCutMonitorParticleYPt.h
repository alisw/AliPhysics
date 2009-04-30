////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticleYPt - the cut monitor for particles to study  ///
/// the difference between reconstructed and true momentum    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticleYPt_hh
#define AliFemtoCutMonitorParticleYPt_hh

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

class AliFemtoCutMonitorParticleYPt : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticleYPt();
  AliFemtoCutMonitorParticleYPt(const char *aName, float aMass);
  AliFemtoCutMonitorParticleYPt(const AliFemtoCutMonitorParticleYPt &aCut);
  virtual ~AliFemtoCutMonitorParticleYPt();

  AliFemtoCutMonitorParticleYPt& operator=(const AliFemtoCutMonitorParticleYPt& aCut);

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
  TH2D *fYPt;    // Rapidity vs. Pt monitor
  TH2D *fYPhi;   // Rapidity cs. Phi monitor
  TH2D *fPtPhi;  // Pt vs. Phi monitor
  TH2D *fEtaPhi; // Pseudorapidity vs. Phi monitor
  TH2D *fEtaPt;  // Pseudorapidity vs. Pt monitor
  TH2D *fEtaPhiW;// Pseudorapidity vs. Phi monitor chi2 weighted
  TH2D *fEtaPtW; // Pseudorapidity vs. Pt monitor chi2 weighted
  float fMass;   // Mass hypothesis
};

#endif
