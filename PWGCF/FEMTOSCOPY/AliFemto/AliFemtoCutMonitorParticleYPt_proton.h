///
/// \file AliFemtoCutMonitorParticleYPt_proton.h
///
/// \class AliFemtoCutMonitorParticleYPt_proton
/// \brief A cut monitor for particles storing the difference between
///        reconstructed and true momentum
///

#ifndef AliFemtoCutMonitorParticleYPt_proton_hh
#define AliFemtoCutMonitorParticleYPt_proton_hh

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

class AliFemtoCutMonitorParticleYPt_proton : public AliFemtoCutMonitor {
public:
  AliFemtoCutMonitorParticleYPt_proton();
  AliFemtoCutMonitorParticleYPt_proton(const char *aName, float aMass);
  AliFemtoCutMonitorParticleYPt_proton(const AliFemtoCutMonitorParticleYPt_proton &aCut);
  virtual ~AliFemtoCutMonitorParticleYPt_proton();

  AliFemtoCutMonitorParticleYPt_proton& operator=(const AliFemtoCutMonitorParticleYPt_proton& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack);
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoXi* aXi) {AliFemtoCutMonitor::Fill(aXi);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2) {AliFemtoCutMonitor::Fill(aCollection1, aCollection2);}

  void Write();

  virtual TList *GetOutputList();

private:
  TH2D *fYPt;     // Rapidity vs. Pt monitor
  TH2D *fYPhi;    // Rapidity cs. Phi monitor
  TH2D *fPtPhi;   // Pt vs. Phi monitor
  TH2D *fEtaPhi;  // Pseudorapidity vs. Phi monitor
  TH2D *fEtaPt;   // Pseudorapidity vs. Pt monitor
  TH2D *fEtaPhiW; // Pseudorapidity vs. Phi monitor chi2 weighted
  TH2D *fEtaPtW;  // Pseudorapidity vs. Pt monitor chi2 weighted
  TH2D *fDCARPt;  // Pt vs. DCA XY
  TH2D *fDCAZPt;  // Pt vs. DCA Z
  float fMass;    // Mass hypothesis
};

#endif
