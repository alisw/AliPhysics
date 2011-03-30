////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticlePtPDG - the cut monitor for particles to study  ///
/// the difference between reconstructed and true momentum    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOCUTMONITORPARTICLEPTPDG_H
#define ALIFEMTOCUTMONITORPARTICLEPTPDG_H

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
  TH2D *ftpcHist;  // TPC histogram
  TH1D *fPtGoodPi; // Good pions
  TH1D *fPtFakePi; // Fake pions
  TH1D *fPtGoodK;  // Good kaons
  TH1D *fPtFakeK;  // Fake kaons
  TH1D *fPtGoodP;  // Good protons
  TH1D *fPtFakeP;  // Fake protons
  TH1D *fPtRPi;    // Pions pt 
  TH1D *fPtRK;     // Kaons pt
  TH1D *fPtRP;     // Protons pt
  
  TH1D *fPtContP;  // Contamination protons
  TH1D *fPtContPi; // Contamination pions
  TH1D *fPtContMup;// Contamination muons
  TH1D *fPtContElp;// Contamination electrons
  
  
  //TH2D *fYPhi;   // Rapidity cs. Phi monitor
  //TH2D *fPtPhi;  // Pt vs. Phi monitor
  //TH2D *fEtaPhi; // Pseudorapidity vs. Phi monitor
  //TH2D *fEtaPt;  // Pseudorapidity vs. Pt monitor
  float fMass;   // Mass hypothesis
};

#endif

