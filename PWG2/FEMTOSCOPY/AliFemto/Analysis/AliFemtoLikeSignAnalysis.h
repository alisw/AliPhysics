/***************************************************************************
 *      This is an analysis which calculated the background from like sign
 *      pairs in the same event
 *      Frank Laue, Ohio State, 2000
 ***************************************************************************/


#ifndef AliFemtoLikeSignAnalysis_hh
#define AliFemtoLikeSignAnalysis_hh
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "Base/AliFemtoBaseAnalysis.h"        // base analysis class
#include "Infrastructure/AliFemtoTypes.h"
#include "Base/AliFemtoEventCut.h"             // base class 
#include "Base/AliFemtoParticleCut.h"          // base class
#include "Base/AliFemtoPairCut.h"              // base class
#include "Base/AliFemtoLikeSignCorrFctn.h"    // base class
#include "Analysis/AliFemtoAnalysis.h"
#include "Infrastructure/AliFemtoCorrFctnCollection.h"


class AliFemtoLikeSignAnalysis : public AliFemtoAnalysis {

public: 

  AliFemtoLikeSignAnalysis(unsigned int bins=20, double min=-100., double max=100.);
  AliFemtoLikeSignAnalysis(const AliFemtoLikeSignAnalysis& OriginalAnalysis);  // copy constructor
  virtual ~AliFemtoLikeSignAnalysis();

  virtual void ProcessEvent(const AliFemtoEvent* TheEventToBeProcessed);
  virtual AliFemtoString Report();
  virtual unsigned int Overflow() { return fOverFlow;}
  virtual unsigned int Underflow() { return fUnderFlow;}

protected:
  double fVertexZ[2];           /* min/max z-vertex position allowed to be processed */
  unsigned int fVertexBins;     /* number of mixing bins in z-vertex in EventMixing Buffer */
  unsigned int fOverFlow;       /* number of events encountered which had too large z-vertex */
  unsigned int fUnderFlow;      /* number of events encountered which had too small z-vertex */

#ifdef __ROOT__
  ClassDef(AliFemtoLikeSignAnalysis, 0)
#endif

};


#endif
