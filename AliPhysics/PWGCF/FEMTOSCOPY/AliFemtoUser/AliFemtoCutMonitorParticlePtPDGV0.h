////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticlePtPDG - the cut monitor for particles to study  ///
/// the difference between reconstructed and true momentum    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOCUTMONITORPARTICLEPTPDGV0_H
#define ALIFEMTOCUTMONITORPARTICLEPTPDGV0_H

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

class AliFemtoCutMonitorParticlePtPDGV0 : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticlePtPDGV0();
  AliFemtoCutMonitorParticlePtPDGV0(const char *aName, float aMass);
  AliFemtoCutMonitorParticlePtPDGV0(const AliFemtoCutMonitorParticlePtPDGV0 &aCut);
  virtual ~AliFemtoCutMonitorParticlePtPDGV0();

  AliFemtoCutMonitorParticlePtPDGV0& operator=(const AliFemtoCutMonitorParticlePtPDGV0& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack)  {AliFemtoCutMonitor::Fill(aTrack);};
  virtual void Fill(const AliFemtoV0* aV0);
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2) {AliFemtoCutMonitor::Fill(aCollection1, aCollection2);}


  void Write();

  virtual TList *GetOutputList();

private:
  TH2D *fPtPDG;    // Rapidity vs. Pt monitor
  TH2D *ftpcHist;  // TPC histogram
  TH1D *fPtMostProbable; // V0s we have
  TH1D *fPtFakeLambdas; // Fake lambdas
  TH1D *fFakeProtonDaughters;  // Fake proton daughters
  TH1D *fFakeAntiProtonDaughters;  // Fake antiproton daughters
  TH1D *fFakePionPlusDaughters;  // Fake pion daughters
  TH1D *fFakePionMinusDaughters;  // Fake pion minus daughters


  TH1D *fPtV0;  // Pt spectra of real lambdas
  TH1D *fPtPosProton;// Pt spectra of proton positive daughters
  TH1D *fPtNegProton;// Pt spectra of proton negative daughters
  TH1D *fPtPosPion; //  Pt spectra of pion positive daughters
  TH1D *fPtNegPion;//  Pt spectra of pion negative daughters

  
  //TH2D *fYPhi;   // Rapidity cs. Phi monitor
  //TH2D *fPtPhi;  // Pt vs. Phi monitor
  //TH2D *fEtaPhi; // Pseudorapidity vs. Phi monitor
  //TH2D *fEtaPt;  // Pseudorapidity vs. Pt monitor
  float fMass;   // Mass hypothesis
};

#endif

