////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticleNumber - the cut monitor for particles to study  ///
/// the difference between reconstructed and true momentum    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOCUTMONITORPARTICLENUMBER_H
#define ALIFEMTOCUTMONITORPARTICLENUMBER_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
class TH3F;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorParticleNumber : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticleNumber();
  AliFemtoCutMonitorParticleNumber(const char *title, const int &pT1Bins, const double& pT1min, const double& pT1max, const int &zvtxBins, const double& zvtxmin, const double& zvtxmax, const int &multBins, const int& multmin, const int& multmax);
  AliFemtoCutMonitorParticleNumber(const AliFemtoCutMonitorParticleNumber &aCut);
  virtual ~AliFemtoCutMonitorParticleNumber();

  AliFemtoCutMonitorParticleNumber& operator=(const AliFemtoCutMonitorParticleNumber& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack);
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
  TH3F *fNtrig;        //Number of triggers
  
  
  //TH2D *fYPhi;   // Rapidity cs. Phi monitor
  //TH2D *fPtPhi;  // Pt vs. Phi monitor
  //TH2D *fEtaPhi; // Pseudorapidity vs. Phi monitor
  //TH2D *fEtaPt;  // Pseudorapidity vs. Pt monitor
};

#endif

