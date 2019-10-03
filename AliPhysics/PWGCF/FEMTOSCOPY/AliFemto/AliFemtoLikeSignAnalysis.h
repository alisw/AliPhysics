/// \class AliFemtoLikeSignAnalysis
/// \brief This is an analysis which calculated the background from like sign
///
/// pairs in the same event

#ifndef ALIFEMTOLIKESIGNANALYSIS_H
#define ALIFEMTOLIKESIGNANALYSIS_H
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoAnalysis.h"        // base analysis class
#include "AliFemtoTypes.h"
//#include "AliFemtoEventCut.h"             // base class
//#include "AliFemtoParticleCut.h"          // base class
//#include "AliFemtoPairCut.h"              // base class
#include "AliFemtoLikeSignCorrFctn.h"    // base class
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoCorrFctnCollection.h"

class AliFemtoEventCut;
class AliFemtoParticleCut;
class AliFemtoPairCut;

class AliFemtoLikeSignAnalysis : public AliFemtoSimpleAnalysis {

public:

  AliFemtoLikeSignAnalysis(unsigned int bins=20, double min=-100., double max=100.);
  AliFemtoLikeSignAnalysis(const AliFemtoLikeSignAnalysis& OriginalAnalysis);  // copy constructor
  AliFemtoLikeSignAnalysis& operator=(const AliFemtoLikeSignAnalysis& OriginalAnalysis);
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
  /// \cond CLASSIMP
  ClassDef(AliFemtoLikeSignAnalysis, 0);
  /// \endcond
#endif

};


#endif
