///
/// \file AliFemtoVertexMultAnalysis.cxx
/// \author Frank Laue, Ohio State, laue@mps.ohio-state.edu
///

#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoXiCut.h"
#include "AliFemtoKinkCut.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include "AliFemtoPicoEventCollectionVectorHideAway.h"


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoVertexMultAnalysis);
  /// \endcond
#endif


//____________________________
AliFemtoVertexMultAnalysis::AliFemtoVertexMultAnalysis(UInt_t binsVertex,
                                                       Double_t minVertex,
                                                       Double_t maxVertex,
                                                       UInt_t binsMult,
                                                       Double_t minMult,
                                                       Double_t maxMult):
  AliFemtoSimpleAnalysis(),
  fVertexZBins(binsVertex),
  fOverFlowVertexZ(0),
  fUnderFlowVertexZ(0),
  fMultBins(binsMult),
  fOverFlowMult(0),
  fUnderFlowMult(0)
{
  fVertexZ[0] = minVertex;
  fVertexZ[1] = maxVertex;
  fMult[0] = minMult;
  fMult[1] = maxMult;

  // We COULD swap these automatically, but we will just print out a warning
  // so users may potentially fix bugs of a larger scale.
  if (minVertex >= maxVertex) {
    cout << "W-AliFemtoVertexMultAnalysis: Provided minVertex >= maxVertex "
            "(" << minVertex << " >= " << maxVertex << "). "
            "No events are expected to pass.";
  }

  if (minMult >= maxMult) {
    cout << "W-AliFemtoVertexMultAnalysis: Provided minMult >= maxMult "
            "(" << minMult << " >= " << maxMult << "). "
            "No events are expected to pass.";
  }

  // create if the correlation function collection was not set in previous
  // constructor (though it SHOULD be)
  if (fCorrFctnCollection == NULL) {
    fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  }

  // always keep fMixingBuffer NULL unless in ProcessEvent()
  if (fMixingBuffer) {
    delete fMixingBuffer;
    fMixingBuffer = NULL;
  }

  // if the event collection was already create (it should NOT be) delete
  // before we allocate a new one
  if (fPicoEventCollectionVectorHideAway) {
    delete fPicoEventCollectionVectorHideAway;
  }

  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(
    fVertexZBins,
    fVertexZ[0],
    fVertexZ[1],
    fMultBins,
    fMult[0],
    fMult[1]
  );

}
//____________________________

AliFemtoVertexMultAnalysis::AliFemtoVertexMultAnalysis(const AliFemtoVertexMultAnalysis& orig):
  AliFemtoSimpleAnalysis(orig),
  fVertexZBins(orig.fVertexZBins),
  fOverFlowVertexZ(0),
  fUnderFlowVertexZ(0),
  fMultBins(orig.fMultBins),
  fOverFlowMult(0),
  fUnderFlowMult(0)
{
  fVertexZ[0] = orig.fVertexZ[0];
  fVertexZ[1] = orig.fVertexZ[1];
  fMult[0] = orig.fMult[0];
  fMult[1] = orig.fMult[1];

  if (fMixingBuffer) {
    delete fMixingBuffer;
    fMixingBuffer = NULL;
  }

  // This *should* be NULL from AliFemtoSimpleAnalysis constructor - but delete just in case
  delete fPicoEventCollectionVectorHideAway;

  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(
    fVertexZBins,
    fVertexZ[0],
    fVertexZ[1],
    fMultBins,
    fMult[0],
    fMult[1]
  );

  if (fVerbose) {
    cout << " AliFemtoVertexMultAnalysis::AliFemtoVertexMultAnalysis(const AliFemtoVertexMultAnalysis&) - analysis copied " << endl;
  }
}

AliFemtoVertexMultAnalysis& AliFemtoVertexMultAnalysis::operator=(const AliFemtoVertexMultAnalysis& rhs)
{
  if (this == &rhs) {
    return *this;
  }

  // Allow parent class to copy the cuts and correlation functions
  AliFemtoSimpleAnalysis::operator=(rhs);

  fVertexZBins = rhs.fVertexZBins;
  fMultBins = rhs.fMultBins;

  fVertexZ[0] = rhs.fVertexZ[0];
  fVertexZ[1] = rhs.fVertexZ[1];
  fUnderFlowVertexZ = 0;
  fOverFlowVertexZ = 0;

  fMult[0] = rhs.fMult[0];
  fMult[1] = rhs.fMult[1];
  fUnderFlowMult = 0;
  fOverFlowMult = 0;

  if (fMixingBuffer) {
    delete fMixingBuffer;
    fMixingBuffer = NULL;
  }

  delete fPicoEventCollectionVectorHideAway;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(
    fVertexZBins,
    fVertexZ[0],
    fVertexZ[1],
    fMultBins,
    fMult[0],
    fMult[1]
  );

  return *this;
}
//____________________________
AliFemtoVertexMultAnalysis::~AliFemtoVertexMultAnalysis()
{
  // Superclass deletes all pointer memebers except these:
  delete fPicoEventCollectionVectorHideAway;
}

//____________________________
AliFemtoString AliFemtoVertexMultAnalysis::Report()
{
  /// Prepare a report of the execution
  if (fVerbose) {
    cout << "AliFemtoVertexMultAnalysis - constructing Report...\n";
  }

  TString report("-----------\nHbt AliFemtoVertexMultAnalysis Report:\n");

  report += TString::Format("Events are mixed in %d VertexZ bins in the range %E cm to %E cm.\n", fVertexZBins, fVertexZ[0], fVertexZ[1])
          + TString::Format("Events underflowing: %d\n", fUnderFlowVertexZ)
          + TString::Format("Events overflowing: %d\n", fOverFlowVertexZ)
          + TString::Format("Events are mixed in %d Mult bins in the range %E cm to %E cm.\n", fMultBins, fMult[0], fMult[1])
          + TString::Format("Events underflowing: %d\n", fUnderFlowMult)
          + TString::Format("Events overflowing: %d\n", fOverFlowMult)
          + TString::Format("Now adding AliFemtoSimpleAnalysis(base) Report\n")
          + AliFemtoSimpleAnalysis::Report();

  return AliFemtoString((const char *)report);
}

TList* AliFemtoVertexMultAnalysis::ListSettings()
{
  TList *settings = AliFemtoSimpleAnalysis::ListSettings();

  settings->AddVector(

    new TObjString(
      TString::Format("AliFemtoVertexMultAnalysis.vertex_z.bins=%d", fVertexZBins)
    ),

    new TObjString(
      TString::Format("AliFemtoVertexMultAnalysis.vertex_z.min=%f", fVertexZ[0])
    ),

    new TObjString(
      TString::Format("AliFemtoVertexMultAnalysis.vertex_z.max=%f", fVertexZ[1])
    ),

    new TObjString(
      TString::Format("AliFemtoVertexMultAnalysis.multiplicity.bins=%d", fMultBins)
    ),

    new TObjString(
      TString::Format("AliFemtoVertexMultAnalysis.multiplicity.min=%f", fMult[0])
    ),

    new TObjString(
      TString::Format("AliFemtoVertexMultAnalysis.multiplicity.max=%f", fMult[1])
    ),

  NULL);

  return settings;
}

//_________________________
void AliFemtoVertexMultAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent)
{
  // Perform event processing in bins of z vertex and multiplicity

  // find the correct mixing buffer
  const Double_t vertexZ = hbtEvent->PrimVertPos().z(),
                    mult = hbtEvent->UncorrectedNumberOfPrimaries();

  fMixingBuffer = fPicoEventCollectionVectorHideAway->PicoEventCollection(vertexZ, mult);

  if (!fMixingBuffer) {

    if (vertexZ < fVertexZ[0]) {
      fUnderFlowVertexZ++;
    }
    else if (vertexZ > fVertexZ[1]) {
      fOverFlowVertexZ++;
    }

    if (mult < fMult[0]) {
      fUnderFlowMult++;
    }
    else if (mult > fMult[1]) {
      fOverFlowMult++;
    }

    return;
  }

  // now that fMixingBuffer has been set - call superclass ProcessEvent()
  AliFemtoSimpleAnalysis::ProcessEvent(hbtEvent);

  // NULL out the mixing buffer after event processed
  fMixingBuffer = NULL;
}
