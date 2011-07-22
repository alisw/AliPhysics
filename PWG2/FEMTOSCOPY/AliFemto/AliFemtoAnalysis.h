////////////////////////////////////////////////////////////////////////////////
/// AliFemtoAnalysis - the pure virtual base class for femto analysis        ///
/// All analysis classes must inherit from this one                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoAnalysis_hh
#define AliFemtoAnalysis_hh

#include "AliFemtoTypes.h"
#include <TList.h>
#include <TObjString.h>

class AliFemtoEvent;

class AliFemtoAnalysis{

public:

  AliFemtoAnalysis() { /* noop */ };
  virtual ~AliFemtoAnalysis() { /* noop */ };

  
  virtual AliFemtoString Report() = 0;       //! returns reports of all cuts applied and correlation functions being done
  virtual TList* ListSettings() = 0;         // return list of cut settings for the analysis

  virtual TList* GetOutputList() = 0;        // Return a TList of objects to be written as output
  virtual void ProcessEvent(const AliFemtoEvent* aEventToAnalyze) = 0;

  virtual void Finish() = 0;

};

#endif
