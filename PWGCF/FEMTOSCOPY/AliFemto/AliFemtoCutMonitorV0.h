////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorV0 -                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorV0_H
#define AliFemtoCutMonitorV0_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
class TH1F;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorV0 : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorV0();
  AliFemtoCutMonitorV0(const char *aName);
  AliFemtoCutMonitorV0(const AliFemtoCutMonitorV0 &aCut);
  virtual ~AliFemtoCutMonitorV0();

  AliFemtoCutMonitorV0& operator=(const AliFemtoCutMonitorV0& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack){AliFemtoCutMonitor::Fill(aTrack);}
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
  TH1F *fLambdaMass;     // Mass assuming lambda hypothesis	
  TH1F *fAntiLambdaMass; // Mass assuming antilambda hypothesis
  TH1F *fK0ShortMass;    // Mass assuming k-short hypothesis	
  TH1F *fDcaDaughters;   // DCA of v0 daughters at Decay vertex
  TH1F *fDcaV0ToPrimVertex;// DCA of v0 to primary vertex
  TH1F *fCosPointingAngle; 
  TH1F *fEtaV0;
  TH1F *fPtV0;
  TH1F *fPtPosDaughter;
  TH1F *fPtNegDaughter;
};

#endif
