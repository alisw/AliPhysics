///
/// \file AliFemtoCorrFctn.cxx
///

#include "AliFemtoCorrFctn.h"

AliFemtoCorrFctn::AliFemtoCorrFctn():
  fyAnalysis(nullptr),
  fPairCut(nullptr)
{ /* no-op */
}

AliFemtoCorrFctn::AliFemtoCorrFctn(const AliFemtoCorrFctn& c):
  fyAnalysis(c.fyAnalysis),
  fPairCut(c.fPairCut)
{
}

AliFemtoCorrFctn& AliFemtoCorrFctn::operator=(const AliFemtoCorrFctn& aCorrFctn)
{
  if (this != &aCorrFctn) {
    fyAnalysis = aCorrFctn.fyAnalysis;
    fPairCut = aCorrFctn.fPairCut;
  }
  return *this;
}

void AliFemtoCorrFctn::AddRealPair(AliFemtoPair*)
{
  cout << "Not implemented" << endl;
}
void AliFemtoCorrFctn::AddMixedPair(AliFemtoPair*)
{
  cout << "Not implemented" << endl;
}

void AliFemtoCorrFctn::AddFirstParticle(AliFemtoParticle*, bool)
{
  cout << "Not implemented" << endl;
}
void AliFemtoCorrFctn::AddSecondParticle(AliFemtoParticle*)
{
  cout << "Not implemented" << endl;
}
void AliFemtoCorrFctn::CalculateAnglesForEvent()
{
  cout << "Not implemented" << endl;
}
