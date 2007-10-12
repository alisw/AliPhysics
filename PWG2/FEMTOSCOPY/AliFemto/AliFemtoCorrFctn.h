////////////////////////////////////////////////////////////////////////////////
/// AliFemtoCorrFctn - the pure virtual base class for correlation function  ///
/// All correlation function classes must inherit from this one              ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCorrFctn_hh
#define AliFemtoCorrFctn_hh

#include "AliFemtoAnalysis.h"
#include "AliFemtoEvent.h"
#include "AliFemtoPair.h"

class AliFemtoCorrFctn{

  friend class AliFemtoAnalysis;

public:
  AliFemtoCorrFctn();
  AliFemtoCorrFctn(const AliFemtoCorrFctn& aCorrFctn);
  virtual ~AliFemtoCorrFctn(){/* no-op */};
  AliFemtoCorrFctn& operator=(const AliFemtoCorrFctn& aCorrFctn);

  virtual AliFemtoString Report() = 0;

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPir);

  virtual void EventBegin(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void EventEnd(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void Finish() = 0;

  virtual AliFemtoCorrFctn* Clone() { return 0;}

  AliFemtoAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoAnalysis* aAnalysis);

protected:
  AliFemtoAnalysis* fyAnalysis;

  private:

};

inline void AliFemtoCorrFctn::AddRealPair(AliFemtoPair*) { cout << "Not implemented" << endl; }
inline void AliFemtoCorrFctn::AddMixedPair(AliFemtoPair*) { cout << "Not implemented" << endl; }

inline AliFemtoCorrFctn::AliFemtoCorrFctn(const AliFemtoCorrFctn& c):fyAnalysis(0) {}
inline AliFemtoCorrFctn::AliFemtoCorrFctn(): fyAnalysis(0) {/* no-op */}
inline void AliFemtoCorrFctn::SetAnalysis(AliFemtoAnalysis* analysis) { fyAnalysis = analysis; }
inline AliFemtoCorrFctn& AliFemtoCorrFctn::operator=(const AliFemtoCorrFctn& aCorrFctn) { if (this == &aCorrFctn) return *this; fyAnalysis = aCorrFctn.fyAnalysis; return *this; }

#endif
