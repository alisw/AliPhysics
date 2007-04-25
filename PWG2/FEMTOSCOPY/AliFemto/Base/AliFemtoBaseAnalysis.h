////////////////////////////////////////////////////////////////////////////////
/// AliFemtoBaseAnalysis - the pure virtual base class for femto analysis    ///
/// All analysis classes must inherit from this one                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoBaseAnalysis_hh
#define AliFemtoBaseAnalysis_hh

#include "Infrastructure/AliFemtoTypes.h"

class AliFemtoEvent;

class AliFemtoBaseAnalysis{

public:

  AliFemtoBaseAnalysis() { /* noop */ };
  virtual ~AliFemtoBaseAnalysis() { /* noop */ };

#ifdef __ROOT__
  ClassDef(AliFemtoBaseAnalysis, 0)
#endif 
  
  virtual AliFemtoString Report() = 0;       //! returns reports of all cuts applied and correlation functions being done

  virtual void ProcessEvent(const AliFemtoEvent* aEventToAnalyze) = 0;

  virtual void Finish() = 0;

};

#endif
