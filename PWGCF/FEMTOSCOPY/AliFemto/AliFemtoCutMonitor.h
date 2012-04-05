////////////////////////////////////////////////////////////////////////////////
/// AliFemtoCutMonitor - the  base class for cut monitor                     ///
/// A cut monitor saves the entities that passed and failed the given cut    ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOCUTMONITOR_H
#define ALIFEMTOCUTMONITOR_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h" // Gael 19/06/02
#include <TList.h>

class AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitor(){/* no-op */};
  virtual ~AliFemtoCutMonitor(){/* no-op */};
  virtual AliFemtoString Report(){ 
    string Stemp = "*** no user defined Fill(const AliFemtoEvent*), take from base class"; 
    AliFemtoString returnThis = Stemp;
    return returnThis; 
  }
  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);
  virtual TList *GetOutputList() { TList *tOutputList = new TList(); return tOutputList; };
  virtual void Fill(const AliFemtoEvent* aEvent);
  virtual void Fill(const AliFemtoTrack* aTrack);
  virtual void Fill(const AliFemtoV0* aV0);
  virtual void Fill(const AliFemtoKink* aKink);
  virtual void Fill(const AliFemtoPair* aPair);
  virtual void Fill(const AliFemtoParticleCollection* aCollection);
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection);
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2);
  virtual void Finish() { 
#ifdef STHBTDEBUG
    cout << " *** no user defined Finish(), take from base class" << endl;
#endif
  }
  virtual void Init() { 
#ifdef STHBTDEBUG
    cout << " *** no user defined Init(), take from base class" << endl;
#endif
  }
};

#endif
