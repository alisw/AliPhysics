/***************************************************************************
 * Collection and analysis for vertex dependent event mixing
 * Frank Laue, Ohio State, 2000
 *
 **************************************************************************/

#ifndef AliFemtoVertexAnalysis_hh
#define AliFemtoVertexAnalysis_hh

#include "AliFemtoSimpleAnalysis.h"        // base analysis class

class AliFemtoVertexAnalysis : public AliFemtoSimpleAnalysis {

public:

  AliFemtoVertexAnalysis(unsigned int =10, double =-100., double=+100.);
  AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& OriginalAnalysis);  // copy constructor
  virtual void ProcessEvent(const AliFemtoEvent* ProcessThisEvent);
  virtual ~AliFemtoVertexAnalysis();
  virtual AliFemtoString Report();       //! returns reports of all cuts applied and correlation functions being done
  virtual unsigned int Overflow() { return fOverFlow;}
  virtual unsigned int Underflow() { return fUnderFlow;}
protected:
  double fVertexZ[2];            /* min/max z-vertex position allowed to be processed */
  unsigned int fVertexBins;      /* number of mixing bins in z-vertex in EventMixing Buffer */
  unsigned int fOverFlow;        /* number of events encountered which had too large z-vertex */
  unsigned int fUnderFlow;       /* number of events encountered which had too small z-vertex */
  
#ifdef __ROOT__
  ClassDef(AliFemtoVertexAnalysis, 0)
#endif
    
};

#endif
