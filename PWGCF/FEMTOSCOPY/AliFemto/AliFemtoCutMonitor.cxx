////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitor - the  base class for cut monitor                       //
// A cut monitor saves the entities that passed and failed the given cut      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitor.h"

inline void AliFemtoCutMonitor::EventBegin(const AliFemtoEvent* /* aEvent */ ) 
{ /* no-op */ }

inline void AliFemtoCutMonitor::EventEnd(const AliFemtoEvent* /* aEvent */ ) 
{ /* no-op */ }

inline void AliFemtoCutMonitor::Fill(const AliFemtoEvent* /* aEvent */) { 
  // cut event
#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoEvent*), take from base class" << endl;
#endif
  }

inline void AliFemtoCutMonitor::Fill(const AliFemtoTrack* /* aTrack */) { 
  // cut track
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoTrack*), take from base class" << endl;
#endif
}

inline void AliFemtoCutMonitor::Fill(const AliFemtoV0* /* aV0 */ ) { 
  // cut V0
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoV0Track*), take from base class" << endl;
#endif
}
inline void AliFemtoCutMonitor::Fill(const AliFemtoKink* /* aKink */) { 
  // cut Kink
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoKink*), take from base class" << endl;
#endif
}
 
//-----------------------------------Gael 12/04/02------------------------------------
inline void AliFemtoCutMonitor::Fill(const AliFemtoPair* /* aPair */) { 
  // cut pair
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoPair*), take from base class" << endl;
#endif
}
//-----------------------------------Gael 19/06/02------------------------------------
inline void AliFemtoCutMonitor::Fill(const AliFemtoParticleCollection* /* aCollection */) {
  // cut particle collection
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoParticleCollection*), take from base class" << endl;
#endif
}
//-----------------------------------Gael 19/06/02------------------------------------
inline void AliFemtoCutMonitor::Fill(const AliFemtoEvent* /* aEvent */,const AliFemtoParticleCollection* /* aCollection */) {
  // cut event and particle collection
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoEvent*,const AliFemtoParticleCollection*), take from base class" << endl;
#endif
}
inline void AliFemtoCutMonitor::Fill(const AliFemtoParticleCollection* /* aCollection */,const AliFemtoParticleCollection* /* aCollection */) {
  // cut event and particle collection
  #ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoParticleCollection*,const AliFemtoParticleCollection*), take from base class" << endl;
  #endif
}
