////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFemtoAnalysisReactionPlane - Femtoscopic analysis which mixes event //
// with respect to the z position of the primary vertex and event total   //
// multiplicity and uses only events in certain reaction plane angle bin  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include "AliFemtoAnalysisReactionPlane.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoKinkCut.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include "AliFemtoPicoEventCollectionVectorHideAway.h"

#ifdef __ROOT__ 
ClassImp(AliFemtoAnalysisReactionPlane)
#endif

extern void FillHbtParticleCollection(AliFemtoParticleCut*         partCut,
				      AliFemtoEvent*               hbtEvent,
				      AliFemtoParticleCollection*  partCollection);


//____________________________
AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(unsigned int binsVertex, double minVertex, double maxVertex,
						       unsigned int binsMult, double minMult, double maxMult, unsigned short binsRP) 
  : 
  fVertexZBins(binsVertex), 
  fOverFlowVertexZ(0), 
  fUnderFlowVertexZ(0),
  fMultBins(binsMult) ,
  fOverFlowMult(0),    
  fUnderFlowMult(0),
  fRPBins(binsRP),
  fCurrentRP(0)
{
  //  mControlSwitch     = 0;
  fEventCut          = 0;
  fFirstParticleCut  = 0;
  fSecondParticleCut = 0;
  fPairCut           = 0;
  fCorrFctnCollection= 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fVertexZ[0] = minVertex;
  fVertexZ[1] = maxVertex;
  fUnderFlowVertexZ = 0; 
  fOverFlowVertexZ = 0; 
  fMult[0] = minMult;
  fMult[1] = maxMult;
  fUnderFlowMult = 0; 
  fOverFlowMult = 0; 
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexZBins,fVertexZ[0],fVertexZ[1],
										     fMultBins,fMult[0],fMult[1],
										     fRPBins,0.0,TMath::Pi());
}
//____________________________

AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& a) : 
  AliFemtoSimpleAnalysis(),
  fVertexZBins(a.fVertexZBins), 
  fOverFlowVertexZ(0), 
  fUnderFlowVertexZ(0),
  fMultBins(a.fMultBins) ,
  fOverFlowMult(0),    
  fUnderFlowMult(0),
  fRPBins(a.fRPBins),
  fCurrentRP(0)
{
  //AliFemtoAnalysisReactionPlane();
  fEventCut          = 0;
  fFirstParticleCut  = 0;
  fSecondParticleCut = 0;
  fPairCut           = 0;
  fCorrFctnCollection= 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fVertexZ[0] = a.fVertexZ[0]; 
  fVertexZ[1] = a.fVertexZ[1];
  fUnderFlowVertexZ = 0; 
  fOverFlowVertexZ = 0; 
  fMult[0] = a.fMult[0]; 
  fMult[1] = a.fMult[1];
  fUnderFlowMult = 0; 
  fOverFlowMult = 0; 
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexZBins,fVertexZ[0],fVertexZ[1],
										     fMultBins,fMult[0],fMult[1],
										     fRPBins,0.0,TMath::Pi());

  // find the right event cut
  fEventCut = a.fEventCut->Clone();
  // find the right first particle cut
  fFirstParticleCut = a.fFirstParticleCut->Clone();
  // find the right second particle cut
  if (a.fFirstParticleCut==a.fSecondParticleCut) 
    SetSecondParticleCut(fFirstParticleCut); // identical particle hbt
  else
  fSecondParticleCut = a.fSecondParticleCut->Clone();

  fPairCut = a.fPairCut->Clone();
  
  if ( fEventCut ) {
      SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& a) - event cut set " << endl;
  }
  if ( fFirstParticleCut ) {
      SetFirstParticleCut(fFirstParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& a) - first particle cut set " << endl;
  }
  if ( fSecondParticleCut ) {
      SetSecondParticleCut(fSecondParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& a) - second particle cut set " << endl;
  }  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& a) - pair cut set " << endl;
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=a.fCorrFctnCollection->begin(); iter!=a.fCorrFctnCollection->end();iter++){
    cout << " AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& a) - looking for correlation functions " << endl;
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
    else cout << " AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& a) - correlation function not found " << endl;
  }

  fNumEventsToMix = a.fNumEventsToMix;

  cout << " AliFemtoAnalysisReactionPlane::AliFemtoAnalysisReactionPlane(const AliFemtoAnalysisReactionPlane& a) - analysis copied " << endl;

}
//____________________________
AliFemtoAnalysisReactionPlane::~AliFemtoAnalysisReactionPlane(){
  // now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself
  delete fPicoEventCollectionVectorHideAway;
}

//____________________________
AliFemtoString AliFemtoAnalysisReactionPlane::Report()
{
  // Prepare a report of the execution
  cout << "AliFemtoAnalysisReactionPlane - constructing Report..."<<endl;
  char ctemp[200];
  AliFemtoString temp = "-----------\nHbt AliFemtoAnalysisReactionPlane Report:\n";
  snprintf(ctemp , 200, "Events are mixed in %d VertexZ bins in the range %E cm to %E cm.\n",fVertexZBins,fVertexZ[0],fVertexZ[1]);
  temp += ctemp;
  snprintf(ctemp , 200, "Events underflowing: %d\n",fUnderFlowVertexZ);
  temp += ctemp;
  snprintf(ctemp , 200, "Events overflowing: %d\n",fOverFlowVertexZ);
  temp += ctemp;
  snprintf(ctemp , 200, "Events are mixed in %d Mult bins in the range %E cm to %E cm.\n",fMultBins,fMult[0],fMult[1]);
  temp += ctemp;
  snprintf(ctemp , 200, "Events underflowing: %d\n",fUnderFlowMult);
  temp += ctemp;
  snprintf(ctemp , 200, "Events overflowing: %d\n",fOverFlowMult);
  temp += ctemp;
  snprintf(ctemp , 200, "Now adding AliFemtoSimpleAnalysis(base) Report\n");
  temp += ctemp;
  temp += AliFemtoSimpleAnalysis::Report();
  AliFemtoString returnThis=temp;
  return returnThis;
}
//_________________________
void AliFemtoAnalysisReactionPlane::ProcessEvent(const AliFemtoEvent* hbtEvent) {
  // Perform event processing
  // in bins of z vertex and multiplicity

  // cout << " AliFemtoAnalysisReactionPlane::ProcessEvent(const AliFemtoEvent* hbtEvent) " << endl;
  // get right mixing buffer
  double vertexZ = hbtEvent->PrimVertPos().z();
  double mult = hbtEvent->UncorrectedNumberOfPrimaries();
  fCurrentRP = hbtEvent->ReactionPlaneAngle();
  
  fMixingBuffer = fPicoEventCollectionVectorHideAway->PicoEventCollection(vertexZ,mult,fCurrentRP); 
  if (!fMixingBuffer) {
    if ( vertexZ < fVertexZ[0] ) fUnderFlowVertexZ++;
    if ( vertexZ > fVertexZ[1] ) fOverFlowVertexZ++;
    if ( mult < fMult[0] ) fUnderFlowMult++;
    if ( mult > fMult[1] ) fOverFlowMult++;
    return;
  }
  // call ProcessEvent() from AliFemtoSimpleAnalysis-base
  AliFemtoSimpleAnalysis::ProcessEvent(hbtEvent);
}

double AliFemtoAnalysisReactionPlane::GetCurrentReactionPlane()
{
  return fCurrentRP;
}
