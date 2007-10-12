////////////////////////////////////////////////////////////////////////////////
/// AliFemtoPairCut - the pure virtual base class for the pair cut           ///
/// All pair cuts must inherit from this one                                 ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoPairCut_hh
#define AliFemtoPairCut_hh

#include <string>

class AliFemtoAnalysis;

#include "AliFemtoString.h"
#include "AliFemtoEvent.h"
#include "AliFemtoPair.h"
#include "AliFemtoCutMonitorHandler.h"
#include <TList.h>
#include <TObjString.h>

class AliFemtoPairCut : public AliFemtoCutMonitorHandler {

  friend class AliFemtoAnalysis;

public:

  AliFemtoPairCut();   // default constructor. - Users should write their own
  AliFemtoPairCut(const AliFemtoPairCut& c); // copy constructor
  virtual ~AliFemtoPairCut(){/* no-op */};  // destructor
  AliFemtoPairCut& operator=(const AliFemtoPairCut &aCut);

  virtual bool Pass(const AliFemtoPair* pair) =0;  // true if passes, false if not

  virtual AliFemtoString Report() =0;    // user-written method to return string describing cuts
  virtual TList *ListSettings() =0;
  virtual void EventBegin(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void EventEnd(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual AliFemtoPairCut* Clone() { return 0;}

  // the following allows "back-pointing" from the CorrFctn to the "parent" Analysis
  AliFemtoAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoAnalysis* aAnalysis);    // Set Back pointer to Analysis

protected:
  AliFemtoAnalysis* fyAnalysis; // Link to the base analysis class

#ifdef __ROOT__
  ClassDef(AliFemtoPairCut, 0)
#endif
};


inline AliFemtoPairCut::AliFemtoPairCut(const AliFemtoPairCut& c) :  AliFemtoCutMonitorHandler(), fyAnalysis(0) {  }
inline void AliFemtoPairCut::SetAnalysis(AliFemtoAnalysis* analysis) { fyAnalysis = analysis; }
inline AliFemtoPairCut::AliFemtoPairCut(): AliFemtoCutMonitorHandler(), fyAnalysis(0) {}   // default constructor. - Users should write their own
inline AliFemtoPairCut& AliFemtoPairCut::operator=(const AliFemtoPairCut &aCut) { if (this == &aCut) return *this; fyAnalysis = aCut.fyAnalysis; return *this; }

#endif
