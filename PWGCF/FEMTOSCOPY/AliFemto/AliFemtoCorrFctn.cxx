////////////////////////////////////////////////////////////////////////////////
/// AliFemtoCorrFctn - the pure virtual base class for correlation function  ///
/// All correlation function classes must inherit from this one              ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCorrFctn.h"

void AliFemtoCorrFctn::AddRealPair(AliFemtoPair*) { cout << "Not implemented" << endl; }
void AliFemtoCorrFctn::AddMixedPair(AliFemtoPair*) { cout << "Not implemented" << endl; }

void AliFemtoCorrFctn::AddFirstParticle(AliFemtoParticle*,bool) { cout << "Not implemented" << endl; }
void AliFemtoCorrFctn::AddSecondParticle(AliFemtoParticle*) { cout << "Not implemented" << endl; }
void AliFemtoCorrFctn::CalculateAnglesForEvent(){ cout << "Not implemented" << endl; }

AliFemtoCorrFctn::AliFemtoCorrFctn(const AliFemtoCorrFctn& /* c */):fyAnalysis(0),fPairCut(0x0) {}
AliFemtoCorrFctn::AliFemtoCorrFctn(): fyAnalysis(0),fPairCut(0x0) {/* no-op */}
void AliFemtoCorrFctn::SetAnalysis(AliFemtoAnalysis* analysis) { fyAnalysis = analysis; }
AliFemtoCorrFctn& AliFemtoCorrFctn::operator=(const AliFemtoCorrFctn& aCorrFctn) { if (this == &aCorrFctn) return *this; fyAnalysis = aCorrFctn.fyAnalysis; fPairCut = aCorrFctn.fPairCut; return *this; }

void AliFemtoCorrFctn::EventBegin(const AliFemtoEvent* /* aEvent */) { /* no-op */ }
void AliFemtoCorrFctn::EventEnd(const AliFemtoEvent* /* aEvent */) { /* no-op */ }
void AliFemtoCorrFctn::SetPairSelectionCut(AliFemtoPairCut* aCut) { fPairCut =  aCut; }

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctn);
  /// \endcond
#endif
