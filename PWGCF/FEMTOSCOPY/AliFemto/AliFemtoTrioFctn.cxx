///
/// \file AliFemtoTrioFctn.cxx
///

#include "AliFemtoTrioFctn.h"

AliFemtoTrioFctn::AliFemtoTrioFctn():
  fyAnalysis(nullptr),
  fTrioCut(nullptr)
{ /* no-op */
}

AliFemtoTrioFctn::AliFemtoTrioFctn(const AliFemtoTrioFctn& c):
  fyAnalysis(c.fyAnalysis),
  fTrioCut(c.fTrioCut)
{
}

AliFemtoTrioFctn& AliFemtoTrioFctn::operator=(const AliFemtoTrioFctn& aTrioFctn)
{
  if (this != &aTrioFctn) {
    fyAnalysis = aTrioFctn.fyAnalysis;
    fTrioCut = aTrioFctn.fTrioCut; 
  }
  return *this;
}

void AliFemtoTrioFctn::SetAnalysis(AliFemtoTrioAnalysis* aAnalysis)
{
  fyAnalysis = aAnalysis;
}

void AliFemtoTrioFctn::SetTrioCut(AliFemtoTrioCut* aCut)
{
  fTrioCut = aCut;
}

void AliFemtoTrioFctn::EventBegin(const AliFemtoEvent* /* event */)
{ // no-op
}

void AliFemtoTrioFctn::EventEnd(const AliFemtoEvent* /* event */)
{  // no-op
}

void AliFemtoTrioFctn::AddRealTrio(AliFemtoTrio*)
{
  cout << "AliFemtoTrioFctn::AddRealPair -- Not implemented\n";
}
void AliFemtoTrioFctn::AddMixedTrio(AliFemtoTrio*)
{
  cout << "AliFemtoTrioFctn::AddMixedPair -- Not implemented\n";
}

void AliFemtoTrioFctn::AddFirstParticle(AliFemtoParticle*, bool)
{
  cout << "AliFemtoTrioFctn::AddFirstParticle -- Not implemented\n";
}
void AliFemtoTrioFctn::AddSecondParticle(AliFemtoParticle*)
{
  cout << "AliFemtoTrioFctn::AddSecondParticle -- Not implemented\n";
}
void AliFemtoTrioFctn::AddThirdParticle(AliFemtoParticle*)
{
  cout << "AliFemtoTrioFctn::AddSecondParticle -- Not implemented\n";
}
void AliFemtoTrioFctn::CalculateAnglesForEvent()
{
  cout << "AliFemtoTrioFctn::CalculateAnglesForEvent -- Not implemented\n";
}
