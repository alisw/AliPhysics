///
/// \file AliFemtoVertexAnalysis.h
///

#ifndef ALIFEMTOVERTEXANALYSIS_H
#define ALIFEMTOVERTEXANALYSIS_H

#include "AliFemtoSimpleAnalysis.h"        // base analysis class

/// \class AliFemtoVertexAnalysis
/// \brief Femtoscopic analysis which mixes events with respect to the z
/// position of the primary vertex
///
/// \author Frank Laue, Ohio State, 2000
///
class AliFemtoVertexAnalysis : public AliFemtoSimpleAnalysis {
public:
  AliFemtoVertexAnalysis(unsigned int bins=10, double min=-100., double max=+100.);
  AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& OriginalAnalysis);  /// copy constructor
  AliFemtoVertexAnalysis& operator=(const AliFemtoVertexAnalysis& OriginalAnalysis);
  virtual void ProcessEvent(const AliFemtoEvent* ProcessThisEvent);
  virtual ~AliFemtoVertexAnalysis();

  virtual AliFemtoString Report();     ///< returns reports of all cuts applied and correlation functions being done
  virtual unsigned int Overflow();
  virtual unsigned int Underflow();

protected:
  double fVertexZ[2];       ///< min/max z-vertex position allowed to be processed
  unsigned int fVertexBins; ///< number of mixing bins in z-vertex in EventMixing Buffer
  unsigned int fOverFlow;   ///< number of events encountered which had too large z-vertex
  unsigned int fUnderFlow;  ///< number of events encountered which had too small z-vertex

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoVertexAnalysis, 0);
  /// \endcond
#endif
};

inline unsigned int AliFemtoVertexAnalysis::Overflow()
{
  return fOverFlow;
}
inline unsigned int AliFemtoVertexAnalysis::Underflow()
{
  return fOverFlow;
}


#endif
