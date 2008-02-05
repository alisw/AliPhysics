////////////////////////////////////////////////////////////////////////////////
/// AliFemtoCutMonitor - the  base class for cut monitor                     ///
/// A cut monitor saves the entities that passed and failed the given cut    ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitor_hh
#define AliFemtoCutMonitor_hh

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h" // Gael 19/06/02
#include <TList.h>

class AliFemtoCutMonitor{
  
private:
  
public:
  AliFemtoCutMonitor(){/* no-op */};
  virtual ~AliFemtoCutMonitor(){/* no-op */};
  virtual AliFemtoString Report(){ 
    string Stemp = "*** no user defined Fill(const AliFemtoEvent*), take from base class"; 
    AliFemtoString returnThis = Stemp;
    return returnThis; 
  }
  virtual void EventBegin(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void EventEnd(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual TList *GetOutputList() { TList *tOutputList = new TList(); return tOutputList; };
  virtual void Fill(const AliFemtoEvent* aEvent) { 

#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoEvent*), take from base class" << endl;
#endif
  }
  virtual void Fill(const AliFemtoTrack* aTrack) { 
#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoTrack*), take from base class" << endl;
#endif
  }
  virtual void Fill(const AliFemtoV0* aV0) { 
#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoV0Track*), take from base class" << endl;
#endif
  }
  virtual void Fill(const AliFemtoKink* aKink) { 
#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoKink*), take from base class" << endl;
#endif
  }
  //-----------------------------------Gael 12/04/02------------------------------------
  virtual void Fill(const AliFemtoPair* aPair) { 
#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoPair*), take from base class" << endl;
#endif
  }
  //-----------------------------------Gael 19/06/02------------------------------------
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {
#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoParticleCollection*), take from base class" << endl;
#endif
  }
  //-----------------------------------Gael 19/06/02------------------------------------
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection) {
#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoEvent*,const AliFemtoParticleCollection*), take from base class" << endl;
#endif
  }
  // -------------------------------------------------------------------------------------
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
