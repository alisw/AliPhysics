///
/// \file AliFemtoCutMonitor.cxx
///

#include "AliFemtoCutMonitor.h"

void AliFemtoCutMonitor::EventBegin(const AliFemtoEvent* /* aEvent */ )
{ /* no-op */
}

void AliFemtoCutMonitor::EventEnd(const AliFemtoEvent* /* aEvent */)
{ /* no-op */
}

void AliFemtoCutMonitor::Fill(const AliFemtoEvent* /* aEvent */)
{
  // cut event
}

void AliFemtoCutMonitor::Fill(const AliFemtoTrack* /* aTrack */)
{
  // cut track
}

void AliFemtoCutMonitor::Fill(const AliFemtoV0* /* aV0 */)
{
  // cut V0
}

void AliFemtoCutMonitor::Fill(const AliFemtoXi* /* aXi */)
{
  // cut Xi
}
void AliFemtoCutMonitor::Fill(const AliFemtoKink* /* aKink */)
{
  // cut Kink
}

void AliFemtoCutMonitor::Fill(const AliFemtoPair* /* aPair */)
{
  // cut pair
}

void AliFemtoCutMonitor::Fill(const AliFemtoParticleCollection* /* aCollection */)
{
  // cut particle collection
}

void AliFemtoCutMonitor::Fill(const AliFemtoEvent* /* aEvent */,
                              const AliFemtoParticleCollection* /* aCollection */)
{
  // cut event and particle collection
}

void AliFemtoCutMonitor::Fill(const AliFemtoParticleCollection* /* aCollection */,
                              const AliFemtoParticleCollection* /* aCollection */)
{
  // cut event and particle collection
}
