////////////////////////////////////////////////////////////////////////////////
/// AliFemtoPairCut - the pure virtual base class for the pair cut           ///
/// All pair cuts must inherit from this one                                 ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoPairCut_hh
#define AliFemtoPairCut_hh

#include <string>

class AliFemtoBaseAnalysis;

#include "Infrastructure/AliFemtoString.h"
#include "Infrastructure/AliFemtoEvent.h"
#include "Infrastructure/AliFemtoPair.h"
#include "Infrastructure/AliFemtoCutMonitorHandler.h"

class AliFemtoPairCut : public AliFemtoCutMonitorHandler {

  friend class AliFemtoBaseAnalysis;

public:

  AliFemtoPairCut();   // default constructor. - Users should write their own
  AliFemtoPairCut(const AliFemtoPairCut& c); // copy constructor
  virtual ~AliFemtoPairCut(){/* no-op */};  // destructor
  AliFemtoPairCut& operator=(const AliFemtoPairCut &aCut);

  virtual bool Pass(const AliFemtoPair* pair) =0;  // true if passes, false if not

  virtual AliFemtoString Report() =0;    // user-written method to return string describing cuts
  virtual void EventBegin(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void EventEnd(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual AliFemtoPairCut* Clone() { return 0;}

  // the following allows "back-pointing" from the CorrFctn to the "parent" Analysis
  AliFemtoBaseAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoBaseAnalysis* aAnalysis);    // Set Back pointer to Analysis

protected:
  AliFemtoBaseAnalysis* fyAnalysis; // Link to the base analysis class

#ifdef __ROOT__
  ClassDef(AliFemtoPairCut, 0)
#endif
};


inline AliFemtoPairCut::AliFemtoPairCut(const AliFemtoPairCut& c) :  AliFemtoCutMonitorHandler(), fyAnalysis(0) {  }
inline void AliFemtoPairCut::SetAnalysis(AliFemtoBaseAnalysis* analysis) { fyAnalysis = analysis; }
inline AliFemtoPairCut::AliFemtoPairCut(): AliFemtoCutMonitorHandler(), fyAnalysis(0) {};   // default constructor. - Users should write their own
inline AliFemtoPairCut& AliFemtoPairCut::operator=(const AliFemtoPairCut &aCut) { if (this == &aCut) return *this; fyAnalysis = aCut.fyAnalysis; return *this; }

#endif
