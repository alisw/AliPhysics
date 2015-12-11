/// \class AliFemtoAnalysisReactionPlane
/// \brief AliFemtoAnalysisReactionPlane - Femtoscopic analysis which mixes event
///
/// with respect to the z position of the primary vertex and event total
/// multiplicity and uses only events in certain reaction plane angle bin

#ifndef ALIFEMTOANALYSISREACTIONPLANE_H
#define ALIFEMTOANALYSISREACTIONPLANE_H

#include "AliFemtoSimpleAnalysis.h"        // base analysis class

class AliFemtoAnalysisReactionPlane : public AliFemtoSimpleAnalysis {

public:

  AliFemtoAnalysisReactionPlane(unsigned int binsVertex=10, double minVertex=-100., double maxVertex=+100., unsigned int binsMult=10, double minMult=-1.e9, double maxMult=+1.e9, unsigned short binsRP=10);
  AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& TheOriginalAnalysis);  // copy constructor
  AliFemtoAnalysisReactionPlane& operator=(const AliFemtoAnalysisReactionPlane& TheOriginalAnalysis);
  virtual void ProcessEvent(const AliFemtoEvent* ProcessThisEvent);
  virtual ~AliFemtoAnalysisReactionPlane();
  virtual AliFemtoString Report();       //! returns reports of all cuts applied and correlation functions being done
  virtual unsigned int OverflowVertexZ() const { return fOverFlowVertexZ;}
  virtual unsigned int UnderflowVertexZ() const { return fUnderFlowVertexZ;}
  virtual unsigned int OverflowMult() const { return fOverFlowMult;}
  virtual unsigned int UnderflowMult() const { return fUnderFlowMult;}
  double GetCurrentReactionPlane();

protected:
  double fVertexZ[2];                 /* min/max z-vertex position allowed to be processed */
  unsigned int fVertexZBins;          /* number of VERTEX mixing bins in z-vertex in EventMixing Buffer */
  unsigned int fOverFlowVertexZ;      /* number of events encountered which had too large z-vertex */
  unsigned int fUnderFlowVertexZ;     /* number of events encountered which had too small z-vertex */
  double fMult[2];                    /* min/max multiplicity allowed for event to be processed */
  unsigned int fMultBins;             /* number of MULTIPLICITY mixing bins in z-vertex in EventMixing Buffer */
  unsigned int fOverFlowMult;         /* number of events encountered which had too large multiplicity */
  unsigned int fUnderFlowMult;        /* number of events encountered which had too small multiplicity */
  unsigned short fRPBins;             // Number of reaction plane angle orientation bins
  double fCurrentRP;                  // Reaction plane angle of the current event

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoAnalysisReactionPlane, 0);
  /// \endcond
#endif

};

#endif
