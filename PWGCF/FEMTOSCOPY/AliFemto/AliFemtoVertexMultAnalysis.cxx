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
ClassImp(AliFemtoVertexMultAnalysis)
/// \endcond
#endif


//____________________________
AliFemtoVertexMultAnalysis::AliFemtoVertexMultAnalysis(unsigned int binsVertex, double minVertex, double maxVertex,
                                                       unsigned int binsMult, double minMult, double maxMult):
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

  // create if the correlation function collection was not set in previous constructor (though it SHOULD be)
  if (fCorrFctnCollection == NULL) {
    fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  }

  if (fMixingBuffer) {
    delete fMixingBuffer;
    fMixingBuffer = NULL;
  }

  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexZBins, fVertexZ[0], fVertexZ[1],
										        fMultBins,    fMult[0],    fMult[1]);
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

  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexZBins, fVertexZ[0], fVertexZ[1],
                                                                                        fMultBins,    fMult[0],    fMult[1]);

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
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexZBins, fVertexZ[0], fVertexZ[1],
                                                                                        fMultBins,    fMult[0],    fMult[1]);
  return *this;
}
//____________________________
AliFemtoVertexMultAnalysis::~AliFemtoVertexMultAnalysis()
{
  /// now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself

  delete fPicoEventCollectionVectorHideAway;
}

//____________________________
AliFemtoString AliFemtoVertexMultAnalysis::Report()
{
  /// Prepare a report of the execution

  // cout << "AliFemtoVertexMultAnalysis - constructing Report..."<<endl;

  TString report("-----------\nHbt AliFemtoVertexMultAnalysis Report:\n");

  report += TString::Format("Events are mixed in %d VertexZ bins in the range %E cm to %E cm.\n", fVertexZBins, fVertexZ[0], fVertexZ[1]);
  report += TString::Format("Events underflowing: %d\n", fUnderFlowVertexZ);
  report += TString::Format("Events overflowing: %d\n", fOverFlowVertexZ);
  report += TString::Format("Events are mixed in %d Mult bins in the range %E cm to %E cm.\n", fMultBins, fMult[0], fMult[1]);
  report += TString::Format("Events underflowing: %d\n", fUnderFlowMult);
  report += TString::Format("Events overflowing: %d\n", fOverFlowMult);
  report += TString::Format("Now adding AliFemtoSimpleAnalysis(base) Report\n");
  report += AliFemtoSimpleAnalysis::Report();

  return AliFemtoString(report);
}
//_________________________
void AliFemtoVertexMultAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent)
{
  /// Perform event processing in bins of z vertex and multiplicity

  // find the correct mixing buffer
  const double vertexZ = hbtEvent->PrimVertPos().z(),
                  mult = hbtEvent->UncorrectedNumberOfPrimaries();

  fMixingBuffer = fPicoEventCollectionVectorHideAway->PicoEventCollection(vertexZ, mult);

  if (!fMixingBuffer) {
    if ( vertexZ < fVertexZ[0] ) fUnderFlowVertexZ++;
    if ( vertexZ > fVertexZ[1] ) fOverFlowVertexZ++;
    if ( mult < fMult[0] ) fUnderFlowMult++;
    if ( mult > fMult[1] ) fOverFlowMult++;
    return;
  }

  // now that fMixingBuffer has been set - call superclass ProcessEvent()
  AliFemtoSimpleAnalysis::ProcessEvent(hbtEvent);

  // NULL out the mixing buffer after event processed
  fMixingBuffer = NULL;
}
