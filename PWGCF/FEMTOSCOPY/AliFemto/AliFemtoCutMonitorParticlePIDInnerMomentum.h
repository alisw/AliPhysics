////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticlePID - the cut monitor for particles to study   ///
/// various aspects of the PID determination                                 ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticlePIDInnerMomentum_hh
#define AliFemtoCutMonitorParticlePIDInnerMomentum_hh

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

class AliFemtoCutMonitorParticlePIDInnerMomentum: public AliFemtoCutMonitor{

public:
  AliFemtoCutMonitorParticlePIDInnerMomentum();
  AliFemtoCutMonitorParticlePIDInnerMomentum(const char *aName, Int_t aTOFParticle, Double_t yTOFTimeMin= -4000.0, Double_t yTOFTimeMax=4000.0);
  AliFemtoCutMonitorParticlePIDInnerMomentum(const AliFemtoCutMonitorParticlePIDInnerMomentum &aCut);
  virtual ~AliFemtoCutMonitorParticlePIDInnerMomentum();

  AliFemtoCutMonitorParticlePIDInnerMomentum& operator=(const AliFemtoCutMonitorParticlePIDInnerMomentum& aCut);

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
 /* void SetTOFParticle(Int_t ipart); */

  void Write();
  void SetUsePt(Bool_t usept){fIfUsePt=usept;}
  virtual TList *GetOutputList();

protected:
  Int_t fTOFParticle; ///< Select TOF time hypothesis, 0-pion, 1-kaon, 2-proton
  Bool_t fIfUsePt;    ///< Plot pT instead of p in all momentum histograms

  TH2D *fTPCdEdx;        ///< TPC dEdx information
  TH2D *fTOFTime;        ///< TOF time
  TH2D *fTOFNSigma;      ///< TOF NSigma values vs mom
  TH2D *fTPCNSigma;      ///< TPC NSigma values vs mom
  TH2D *fTPCTOFNSigma;   ///< TPC^2+ TOF^2 NSigma values vs mom
  TH2D *fTPCvsTOFNSigma; ///< TPC vs TOF

  TH1D *fParticleOrigin; ///< particle origin from MC
  TH1D *fParticleId;     ///< true particle identification from MC
};

#endif
