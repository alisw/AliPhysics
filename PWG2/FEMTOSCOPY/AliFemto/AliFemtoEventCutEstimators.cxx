////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventCutEstimators - the basic cut for events.                          //
// Only cuts on event multiplicity and z-vertex position                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoEventCutEstimators.h"
//#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoEventCutEstimators)
#endif

AliFemtoEventCutEstimators::AliFemtoEventCutEstimators() :
  AliFemtoEventCut(),
  fEventMultEst1(),
  fEventMultEst2(),
  fEventMultEst3(),
  fUseMultEst1(0), 
  fUseMultEst2(0), 
  fUseMultEst3(0), 
  fEventCentEst1(),
  fEventCentEst2(),
  fEventCentEst3(),
  fEventCentEst4(),
  fUseCentEst1(0), 
  fUseCentEst2(0), 
  fUseCentEst3(0), 
  fUseCentEst4(0),
  fNEventsPassed(0), 
  fNEventsFailed(0) 
{
  // Default constructor
  fEventMultEst1[0] = 0;  fEventMultEst1[1] = 10000;
  fEventMultEst2[0] = 0;  fEventMultEst2[1] = 10000;
  fEventMultEst3[0] = 0;  fEventMultEst3[1] = 10000;
  fEventCentEst1[0] = 0;  fEventCentEst1[1] = 1000.0;
  fEventCentEst2[0] = 0;  fEventCentEst2[1] = 1000.0;
  fEventCentEst3[0] = 0;  fEventCentEst3[1] = 1000.0;
  fEventCentEst4[0] = 0;  fEventCentEst4[1] = 1000.0;
  fVertZPos[0] = -100.0;
  fVertZPos[1] = 100.0;
} 
//------------------------------
AliFemtoEventCutEstimators::~AliFemtoEventCutEstimators(){
  // Default destructor
}
//------------------------------
bool AliFemtoEventCutEstimators::Pass(const AliFemtoEvent* event){
  // Pass events if they fall within the multiplicity and z-vertex
  // position range. Fail otherwise
  //  int mult =  event->NumberOfTracks();
  
  bool goodEvent = true;

  printf("Cutting event with %i %i %i - %i %i %i %i\n", fUseMultEst1, fUseMultEst2, fUseMultEst3, fUseCentEst1, fUseCentEst2, fUseCentEst3, fUseCentEst4);
  printf("  On %i %i %i - %f %f %f %f\n", event->MultiplicityEstimateTracklets(), event->MultiplicityEstimateITSTPC(), event->MultiplicityEstimateITSPure(),
	 event->CentralityV0(), event->CentralityFMD(), event->CentralitySPD1(), event->CentralityTrk());

  if (fUseMultEst1) { goodEvent &= ((event->MultiplicityEstimateTracklets() >= fEventMultEst1[0]) &&
				    (event->MultiplicityEstimateTracklets() <= fEventMultEst1[1])); }
  if (fUseMultEst2) { goodEvent &= ((event->MultiplicityEstimateITSTPC() >= fEventMultEst2[0]) &&
				    (event->MultiplicityEstimateITSTPC() <= fEventMultEst2[1])); }
  if (fUseMultEst3) { goodEvent &= ((event->MultiplicityEstimateITSPure() >= fEventMultEst3[0]) &&
				    (event->MultiplicityEstimateITSPure() <= fEventMultEst3[1])); }

  if (fUseCentEst1) { goodEvent &= ((event->CentralityV0() > fEventCentEst1[0]) &&
				    (event->CentralityV0() < fEventCentEst1[1])); }
  if (fUseCentEst2) { goodEvent &= ((event->CentralityFMD() > fEventCentEst2[0]) &&
				    (event->CentralityFMD() < fEventCentEst2[1])); }
  if (fUseCentEst3) { goodEvent &= ((event->CentralitySPD1() > fEventCentEst3[0]) &&
				    (event->CentralitySPD1() < fEventCentEst3[1])); }
  if (fUseCentEst4) { goodEvent &= ((event->CentralityTrk() > fEventCentEst4[0]) &&
				    (event->CentralityTrk() < fEventCentEst4[1])); }
  double vertexZPos = event->PrimVertPos().z();
  //   cout << "AliFemtoEventCutEstimators:: mult:       " << fEventMult[0] << " < " << mult << " < " << fEventMult[1] << endl;
  //   cout << "AliFemtoEventCutEstimators:: VertexZPos: " << fVertZPos[0] << " < " << vertexZPos << " < " << fVertZPos[1] << endl;
  //   cout << "AliFemtoEventCutEstimators:: VertexZErr: " << event->PrimVertCov()[4] << endl;
  goodEvent &=
    ((vertexZPos > fVertZPos[0]) &&
     (vertexZPos < fVertZPos[1]));
  goodEvent ? fNEventsPassed++ : fNEventsFailed++ ;
  //   cout << "AliFemtoEventCutEstimators:: return : " << goodEvent << endl;
  //     (fAcceptBadVertex || (event->PrimVertCov()[4] > -1000.0)) &&
  return (goodEvent);
}
  //------------------------------
AliFemtoString AliFemtoEventCutEstimators::Report(){
  // Prepare report
  string stemp;
  char ctemp[100];
  snprintf(ctemp , 100, "\nMultiplicity:\t %d-%d",fEventMultEst2[0],fEventMultEst2[1]);
  stemp = ctemp;
  snprintf(ctemp , 100, "\nVertex Z-position:\t %E-%E",fVertZPos[0],fVertZPos[1]);
  stemp += ctemp;
  snprintf(ctemp , 100, "\nNumber of events which passed:\t%ld  Number which failed:\t%ld",fNEventsPassed,fNEventsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;
}
