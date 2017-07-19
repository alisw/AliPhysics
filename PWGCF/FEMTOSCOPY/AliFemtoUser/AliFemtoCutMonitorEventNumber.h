////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorEventNumber - the cut monitor for particles to study  ///
/// the difference between reconstructed and true momentum    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOCUTMONITOREVENTNUMBER_H
#define ALIFEMTOCUTMONITOREVENTENUMBER_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
class TH2F;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorEventNumber : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorEventNumber();
  AliFemtoCutMonitorEventNumber(const char *title, const int &zvtxBins, const double& zvtxmin, const double& zvtxmax, const int &multBins, const int& multmin, const int& multmax);
  AliFemtoCutMonitorEventNumber(const AliFemtoCutMonitorEventNumber &aCut);
  virtual ~AliFemtoCutMonitorEventNumber();

  AliFemtoCutMonitorEventNumber& operator=(const AliFemtoCutMonitorEventNumber& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) ;
  virtual void Fill(const AliFemtoTrack* aTrack) {AliFemtoCutMonitor::Fill(aTrack);}
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
  TH2F *fNevent;        //Number of events
  
  
  //TH2D *fYPhi;   // Rapidity cs. Phi monitor
  //TH2D *fPtPhi;  // Pt vs. Phi monitor
  //TH2D *fEtaPhi; // Pseudorapidity vs. Phi monitor
  //TH2D *fEtaPt;  // Pseudorapidity vs. Pt monitor
};

#endif

