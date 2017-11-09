///
/// \file AliFemtoVertexAnalysis.cxx
///

#include "AliFemtoVertexAnalysis.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoKinkCut.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include "AliFemtoPicoEventCollectionVectorHideAway.h"


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoVertexAnalysis);
  /// \endcond
#endif

//____________________________
AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(unsigned int bins, double min, double max):
  AliFemtoSimpleAnalysis(),
  fVertexBins(bins),
  fOverFlow(0),
  fUnderFlow(0)
{
  fVertexZ[0] = min;
  fVertexZ[1] = max;

  if (fMixingBuffer) {
    delete fMixingBuffer;
    fMixingBuffer = NULL;
  }
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexBins, fVertexZ[0], fVertexZ[1]);
}
//____________________________

AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a):
  AliFemtoSimpleAnalysis(a),
  fVertexBins(a.fVertexBins),
  fOverFlow(0),
  fUnderFlow(0)
{ // copy constructor
  fVertexZ[0] = a.fVertexZ[0];
  fVertexZ[1] = a.fVertexZ[1];

  if (fMixingBuffer) {
    delete fMixingBuffer;
    fMixingBuffer = NULL;
  }
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexBins, fVertexZ[0], fVertexZ[1]);

  if (fVerbose) {
    cout << " AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) - analysis copied\n";
  }
}

AliFemtoVertexAnalysis& AliFemtoVertexAnalysis::operator=(const AliFemtoVertexAnalysis& OriginalAnalysis)
{
  if (this == &OriginalAnalysis) {
    return *this;
  }

  AliFemtoSimpleAnalysis::operator=(OriginalAnalysis);

  fVertexBins = OriginalAnalysis.fVertexBins;
  fVertexZ[0] = OriginalAnalysis.fVertexZ[0];
  fVertexZ[1] = OriginalAnalysis.fVertexZ[1];
  fUnderFlow = 0;
  fOverFlow = 0;

  if (fMixingBuffer) {
    delete fMixingBuffer;
    fMixingBuffer = NULL;
  }
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexBins, fVertexZ[0], fVertexZ[1]);

  return *this;
}

//____________________________
AliFemtoVertexAnalysis::~AliFemtoVertexAnalysis()
{
  /// now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself
  delete fPicoEventCollectionVectorHideAway;
}

//____________________________
AliFemtoString AliFemtoVertexAnalysis::Report()
{
  /// prepare report fromt the execution
  if (fVerbose) {
    cout << "AliFemtoVertexAnalysis - constructing Report...\n";
  }
  TString report("-----------\nHbt AliFemtoVertexAnalysis Report:\n");

  report += TString::Format("Events are mixed in %d bins in the range %E cm to %E cm.\n", fVertexBins, fVertexZ[0], fVertexZ[1])
          + TString::Format("Events underflowing: %d\n", fUnderFlow)
          + TString::Format("Events overflowing: %d\n",fOverFlow)
          + TString::Format("Now adding AliFemtoSimpleAnalysis(base) Report\n");

  report += AliFemtoSimpleAnalysis::Report();

  return AliFemtoString((const char *)report);
}
//_________________________
void AliFemtoVertexAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent)
{
  /// Put the event though the analysis
  if (fVerbose) {
    cout << " AliFemtoVertexAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent)\n";
  }
  // get right mixing buffer
  const double vertexZ = hbtEvent->PrimVertPos().z();
  fMixingBuffer = fPicoEventCollectionVectorHideAway->PicoEventCollection(vertexZ);

  if (fMixingBuffer == NULL) {
    if (vertexZ < fVertexZ[0]) fUnderFlow++;
    if (vertexZ > fVertexZ[1]) fOverFlow++;
    return;
  }

  // call ProcessEvent() from AliFemtoSimpleAnalysis-base
  AliFemtoSimpleAnalysis::ProcessEvent(hbtEvent);

  // NULL out the mixing buffer after event processed
  fMixingBuffer = NULL;
}
