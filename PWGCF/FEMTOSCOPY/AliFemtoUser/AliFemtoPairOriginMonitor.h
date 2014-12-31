////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Author: Dominik Arominski, WUT                                            //
//  dominik.arominski@cern.ch                                                 //
//                                                                            //
//  AliFemtoPairOriginMonitor - allows to extract origin and indentification  //
//  of particles after both one-particle and pair-wise cuts.                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoPairOriginMonitor_H
#define AliFemtoPairOriginMonitor_H

#include "AliFemtoPair.h"
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
class TH1D;
class TH2D;
class TList;

class AliFemtoPairOriginMonitor : public AliFemtoCutMonitor{

public:
  AliFemtoPairOriginMonitor();
  AliFemtoPairOriginMonitor(const char *aName);
  AliFemtoPairOriginMonitor(const AliFemtoPairOriginMonitor &aCut);
  virtual ~AliFemtoPairOriginMonitor();

  AliFemtoPairOriginMonitor& operator=(const AliFemtoPairOriginMonitor& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack){AliFemtoCutMonitor::Fill(aTrack);}
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair);
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2) {AliFemtoCutMonitor::Fill(aCollection1, aCollection2);}
  void Write();

  virtual TList *GetOutputList();

private:
  TH1D *fParticle1Origin; //first particle in pair origin from MC
  TH1D *fParticle2Origin; //second particle in pair origin from MC

  TH1D *fParticle1Id;     //true particle identification from MC
  TH1D *fParticle2Id;     //true particle identification from MC
};

#endif
