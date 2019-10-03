///
/// \file AliFemtoVertexMultAnalysis.h
///

#ifndef ALIFEMTOVERTEXMULTANALYSIS_H
#define ALIFEMTOVERTEXMULTANALYSIS_H

#include "AliFemtoSimpleAnalysis.h"

/// \class AliFemtoVertexMultAnalysis
/// \brief Femtoscopic analysis which mixes events binned by vertices'
///        z-position and multiplicity
///
///
/// Femtoscopic analysis which mixes events with respect to the logitudinal (z)
/// position of the primary vertex and event total multiplicity. These
/// parameters are only set via the constructor.
///
/// Binning is done by using the fPicoEventCollectionVectorHideAway member
/// provided by the superclass AliFemtoSimpleAnalysis. It should be noted that
/// this member is not created or deleted by the superclass, so this class
/// deletes the member in its destructor.
///
class AliFemtoVertexMultAnalysis : public AliFemtoSimpleAnalysis {
public:

  /// Construct with parameters for event-mixing bins
  ///
  /// If the min/max values are backwards or equal, this class prints a warning
  /// but continues onward, potentially ignoring all events.
  ///
  AliFemtoVertexMultAnalysis(UInt_t binsVertex=10,
                             Double_t minVertex=-100.0,
                             Double_t maxVertex=+100.0,
                             UInt_t binsMult=10,
                             Double_t minMult=-1.0e9,
                             Double_t maxMult=+1.0e9);

  /// Copies the event-mixing binning parameters, creates a new event
  /// collection and resets the overflow & underflow members to 0.
  AliFemtoVertexMultAnalysis(const AliFemtoVertexMultAnalysis& TheOriginalAnalysis);
  AliFemtoVertexMultAnalysis& operator=(const AliFemtoVertexMultAnalysis& TheOriginalAnalysis);

  /// Deletes the event collection hideaway, as superclass does not do this
  virtual ~AliFemtoVertexMultAnalysis();

  /// Passes the event to AliFemtoSimpleAnalysis::ProcessEvent after
  /// determining which (if any) mixing buffer to use.
  virtual void ProcessEvent(const AliFemtoEvent* ProcessThisEvent);
  virtual AliFemtoString Report();       ///< returns reports of all cuts applied and correlation functions being done

  /// Return a TList of analysis settings.
  ///
  /// The TList comprises TObjStrings describing the settings provided by the
  /// AliFemtoSimpleAnalysis::ListSettings class followed by all event-mixing
  /// binning parameters.
  virtual TList* ListSettings();

  virtual UInt_t OverflowVertexZ() const;   ///< Number of events above vertex-z range
  virtual UInt_t UnderflowVertexZ() const;  ///< Number of events below vertex-z range
  virtual UInt_t OverflowMult() const;      ///< Number of events above multiplicity range
  virtual UInt_t UnderflowMult() const;     ///< Number of events below multiplicity range

protected:

  Double_t fVertexZ[2];     ///< min/max z-vertex position allowed to be processed
  UInt_t fVertexZBins;      ///< number of VERTEX mixing bins in z-vertex in EventMixing Buffer
  UInt_t fOverFlowVertexZ;  ///< number of events encountered which had too large z-vertex
  UInt_t fUnderFlowVertexZ; ///< number of events encountered which had too small z-vertex

  Double_t fMult[2];        ///< min/max multiplicity allowed for event to be processed
  UInt_t fMultBins;         ///< number of MULTIPLICITY mixing bins in z-vertex in EventMixing Buffer
  UInt_t fOverFlowMult;     ///< number of events encountered which had too large multiplicity
  UInt_t fUnderFlowMult;    ///< number of events encountered which had too small multiplicity

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoVertexMultAnalysis, 1);
  /// \endcond
#endif

};

inline UInt_t AliFemtoVertexMultAnalysis::OverflowVertexZ() const
{
  return fOverFlowVertexZ;
}

inline UInt_t AliFemtoVertexMultAnalysis::UnderflowVertexZ() const
{
  return fUnderFlowVertexZ;
}

inline UInt_t AliFemtoVertexMultAnalysis::OverflowMult() const
{
  return fOverFlowMult;
}

inline UInt_t AliFemtoVertexMultAnalysis::UnderflowMult() const
{
  return fUnderFlowMult;
}

#endif
