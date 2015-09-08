///
/// \file AliFemtoCutMonitor.cxx
///

#include "AliFemtoCutMonitor.h"

void AliFemtoCutMonitor::EventBegin(const AliFemtoEvent* /* aEvent */ ) 
{ /* no-op */ }

void AliFemtoCutMonitor::EventEnd(const AliFemtoEvent* /* aEvent */ ) 
{ /* no-op */ }

void AliFemtoCutMonitor::Fill(const AliFemtoEvent* /* aEvent */) { 
  // cut event
#ifdef STHBTDEBUG
    cout << " *** no user defined Fill(const AliFemtoEvent*), take from base class" << endl;
#endif
  }

void AliFemtoCutMonitor::Fill(const AliFemtoTrack* /* aTrack */) { 
  // cut track
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoTrack*), take from base class" << endl;
#endif
}

void AliFemtoCutMonitor::Fill(const AliFemtoV0* /* aV0 */ ) { 
  // cut V0
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoV0Track*), take from base class" << endl;
#endif
}
void AliFemtoCutMonitor::Fill(const AliFemtoXi* /* aXi */ ) { 
  // cut Xi
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoV0Track*), take from base class" << endl;
#endif
}
void AliFemtoCutMonitor::Fill(const AliFemtoKink* /* aKink */) { 
  // cut Kink
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoKink*), take from base class" << endl;
#endif
}
 
//-----------------------------------Gael 12/04/02------------------------------------
void AliFemtoCutMonitor::Fill(const AliFemtoPair* /* aPair */) { 
  // cut pair
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoPair*), take from base class" << endl;
#endif
}
//-----------------------------------Gael 19/06/02------------------------------------
void AliFemtoCutMonitor::Fill(const AliFemtoParticleCollection* /* aCollection */) {
  // cut particle collection
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoParticleCollection*), take from base class" << endl;
#endif
}
//-----------------------------------Gael 19/06/02------------------------------------
void AliFemtoCutMonitor::Fill(const AliFemtoEvent* /* aEvent */,const AliFemtoParticleCollection* /* aCollection */) {
  // cut event and particle collection
#ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoEvent*,const AliFemtoParticleCollection*), take from base class" << endl;
#endif
}
void AliFemtoCutMonitor::Fill(const AliFemtoParticleCollection* /* aCollection */,const AliFemtoParticleCollection* /* aCollection */) {
  // cut event and particle collection
  #ifdef STHBTDEBUG
  cout << " *** no user defined Fill(const AliFemtoParticleCollection*,const AliFemtoParticleCollection*), take from base class" << endl;
  #endif
}
