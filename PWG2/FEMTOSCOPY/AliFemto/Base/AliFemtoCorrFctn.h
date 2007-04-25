////////////////////////////////////////////////////////////////////////////////
/// AliFemtoCorrFctn - the pure virtual base class for correlation function  ///
/// All correlation function classes must inherit from this one              ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCorrFctn_hh
#define AliFemtoCorrFctn_hh

#include "Base/AliFemtoBaseAnalysis.h"
#include "Infrastructure/AliFemtoEvent.h"
#include "Infrastructure/AliFemtoPair.h"

class AliFemtoCorrFctn{

  friend class AliFemtoBaseAnalysis;

public:
  AliFemtoCorrFctn(){/* no-op */};
  AliFemtoCorrFctn(const AliFemtoCorrFctn& aCorrFctn);
  virtual ~AliFemtoCorrFctn(){/* no-op */};

  virtual AliFemtoString Report() = 0;

  virtual void AddRealPair(const AliFemtoPair* aPair);
  virtual void AddMixedPair(const AliFemtoPair* aPir);

  virtual void EventBegin(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void EventEnd(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void Finish() = 0;

  virtual AliFemtoCorrFctn* Clone() { return 0;}

  AliFemtoBaseAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoBaseAnalysis* aAnalysis);

protected:
  AliFemtoBaseAnalysis* fyAnalysis;

private:

};

inline void AliFemtoCorrFctn::AddRealPair(const AliFemtoPair*) { cout << "Not implemented" << endl; }
inline void AliFemtoCorrFctn::AddMixedPair(const AliFemtoPair*) { cout << "Not implemented" << endl; }

inline AliFemtoCorrFctn::AliFemtoCorrFctn(const AliFemtoCorrFctn& c) { fyAnalysis =0; }
inline void AliFemtoCorrFctn::SetAnalysis(AliFemtoBaseAnalysis* analysis) { fyAnalysis = analysis; }

#endif
