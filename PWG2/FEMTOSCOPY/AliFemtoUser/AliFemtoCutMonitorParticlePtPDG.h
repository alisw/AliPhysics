////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticlePtPDG - the cut monitor for particles to study  ///
/// the difference between reconstructed and true momentum    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticlePtPDG_hh
#define AliFemtoCutMonitorParticlePtPDG_hh

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

class AliFemtoCutMonitorParticlePtPDG : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticlePtPDG();
  AliFemtoCutMonitorParticlePtPDG(const char *aName, float aMass);
  AliFemtoCutMonitorParticlePtPDG(const AliFemtoCutMonitorParticlePtPDG &aCut);
  virtual ~AliFemtoCutMonitorParticlePtPDG();

  AliFemtoCutMonitorParticlePtPDG& operator=(const AliFemtoCutMonitorParticlePtPDG& aCut);

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
  TH2D *fPtPDG;    // Rapidity vs. Pt monitor
  TH2D *ftpcHist;
  TH1D *fPtGoodPi;
  TH1D *fPtFakePi;
  TH1D *fPtGoodK;
  TH1D *fPtFakeK;
  TH1D *fPtGoodP;
  TH1D *fPtFakeP;
  TH1D *fPtRPi;
  TH1D *fPtRK;
  TH1D *fPtRP;  
  
  TH1D *fPtContP;
  TH1D *fPtContPi;
  TH1D *fPtContMup;
  TH1D *fPtContElp;
  
  
  //TH2D *fYPhi;   // Rapidity cs. Phi monitor
  //TH2D *fPtPhi;  // Pt vs. Phi monitor
  //TH2D *fEtaPhi; // Pseudorapidity vs. Phi monitor
  //TH2D *fEtaPt;  // Pseudorapidity vs. Pt monitor
  float fMass;   // Mass hypothesis
};

#endif

