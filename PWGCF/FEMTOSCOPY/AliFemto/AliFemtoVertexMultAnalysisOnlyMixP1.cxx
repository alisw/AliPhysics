//===============================================================================//
// dongfang.wang@cern.ch
// AliFemtoVertexMultAnalysisOnlyMixP1: Making pool which keep update/save
// as long as the first particle has been found in this event!
// Only difference with AliFemtoVertexMultAnalysis.cxx in line 217.
//===============================================================================//

#include "AliFemtoVertexMultAnalysisOnlyMixP1.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoXiCut.h"
#include "AliFemtoKinkCut.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include "AliFemtoPicoEventCollectionVectorHideAway.h"


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoVertexMultAnalysisOnlyMixP1);
  /// \endcond
#endif

AliFemtoVertexMultAnalysisOnlyMixP1::AliFemtoVertexMultAnalysisOnlyMixP1(   UInt_t binsVertex,
                                                                            Double_t minVertex,
                                                                            Double_t maxVertex,
                                                                            UInt_t binsMult,
                                                                            Double_t minMult,
                                                                            Double_t maxMult):
    AliFemtoSimpleAnalysisOnlyMixP1(),
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
AliFemtoVertexMultAnalysisOnlyMixP1::AliFemtoVertexMultAnalysisOnlyMixP1(const AliFemtoVertexMultAnalysisOnlyMixP1 &OriAnalysis):
    AliFemtoSimpleAnalysisOnlyMixP1(OriAnalysis),
    fVertexZBins(OriAnalysis.fVertexZBins),
    fOverFlowVertexZ(0),
    fUnderFlowVertexZ(0),
    fMultBins(OriAnalysis.fMultBins),
    fOverFlowMult(0),
    fUnderFlowMult(0)
{
    //copy constructor 
    fVertexZ[0] = OriAnalysis.fVertexZ[0];
    fVertexZ[1] = OriAnalysis.fVertexZ[1];
    fMult[0] = OriAnalysis.fMult[0];
    fMult[1] = OriAnalysis.fMult[1];

    if (fMixingBuffer) {
        delete fMixingBuffer;
        fMixingBuffer = NULL;
    }

    // This *should* be NULL from AliFemtoVertexMultAnalysisOnlyMixP1 constructor - but delete just in case
    delete fPicoEventCollectionVectorHideAway;

    fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(
        fVertexZBins,
        fVertexZ[0],
        fVertexZ[1],
        fMultBins,
        fMult[0],
        fMult[1]
    );

    // if (fVerbose) {
    //     cout << " AliFemtoVertexMultAnalysisOnlyMixP1::AliFemtoVertexMultAnalysisOnlyMixP1(const AliFemtoVertexMultAnalysisOnlyMixP1&) - analysis copied " << endl;
    // }
}


AliFemtoVertexMultAnalysisOnlyMixP1& AliFemtoVertexMultAnalysisOnlyMixP1::operator=(const AliFemtoVertexMultAnalysisOnlyMixP1& OriAnalysis)
{
    // assignment operator
    if (this == &OriAnalysis) {
        return *this;
    }

    // Allow parent class to copy the cuts and correlation functions
    AliFemtoVertexMultAnalysisOnlyMixP1::operator=(OriAnalysis);

    fVertexZBins = OriAnalysis.fVertexZBins;
    fMultBins = OriAnalysis.fMultBins;

    fVertexZ[0] = OriAnalysis.fVertexZ[0];
    fVertexZ[1] = OriAnalysis.fVertexZ[1];
    fUnderFlowVertexZ = 0;
    fOverFlowVertexZ = 0;

    fMult[0] = OriAnalysis.fMult[0];
    fMult[1] = OriAnalysis.fMult[1];
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

AliFemtoVertexMultAnalysisOnlyMixP1::~AliFemtoVertexMultAnalysisOnlyMixP1()
{
    // Destructor
    delete fPicoEventCollectionVectorHideAway;
}
AliFemtoString AliFemtoVertexMultAnalysisOnlyMixP1::Report()
{
    TString report("-----------\nHbt AliFemtoVertexMultAnalysisOnlyMixP1 Report:\n");

    report += TString::Format("Events are mixed in %d VertexZ bins in the range %E cm to %E cm.\n", fVertexZBins, fVertexZ[0], fVertexZ[1])
            + TString::Format("Events underflowing: %d\n", fUnderFlowVertexZ)
            + TString::Format("Events overflowing: %d\n", fOverFlowVertexZ)
            + TString::Format("Events are mixed in %d Mult bins in the range %E cm to %E cm.\n", fMultBins, fMult[0], fMult[1])
            + TString::Format("Events underflowing: %d\n", fUnderFlowMult)
            + TString::Format("Events overflowing: %d\n", fOverFlowMult)
            + TString::Format("Now adding AliFemtoVertexMultAnalysisOnlyMixP1 Report\n");
            //+ AliFemtoSimpleAnalysis::Report();

    return AliFemtoString((const char *)report);

}
TList* AliFemtoVertexMultAnalysisOnlyMixP1::ListSettings()
{
    TList *settings = new TList;
    return settings;
}

void AliFemtoVertexMultAnalysisOnlyMixP1::ProcessEvent(const AliFemtoEvent* HbtEventToProcess){
    const Double_t  vertexZ = HbtEventToProcess->PrimVertPos().z(),
                    mult = HbtEventToProcess->UncorrectedNumberOfPrimaries();

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
    AliFemtoSimpleAnalysisOnlyMixP1::ProcessEvent(HbtEventToProcess);

    // NULL out the mixing buffer after event processed
    fMixingBuffer = NULL;

}
