////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoBasicEventCut - the basic cut for events.                          //
// Only cuts on event multiplicity and z-vertex position                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoBasicEventCut.h"
//#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoBasicEventCut)
#endif

AliFemtoBasicEventCut::AliFemtoBasicEventCut() :
  AliFemtoEventCut(),
  fEventMult(),
  fVertZPos(),
  fAcceptBadVertex(false), 
  fNEventsPassed(0), 
  fNEventsFailed(0),
  fAcceptOnlyPhysics(true),
  fSelectTrigger(0)
{
  // Default constructor
  fEventMult[0] = 0;
  fEventMult[1] = 100000;
  fVertZPos[0] = -100.0;
  fVertZPos[1] = 100.0;
} 
//------------------------------
AliFemtoBasicEventCut::~AliFemtoBasicEventCut(){
  // Default destructor
}
//------------------------------
bool AliFemtoBasicEventCut::Pass(const AliFemtoEvent* event){
  // Pass events if they fall within the multiplicity and z-vertex
  // position range. Fail otherwise
  //  int mult =  event->NumberOfTracks();
  int mult = (int) event->UncorrectedNumberOfPrimaries();
  double vertexZPos = event->PrimVertPos().z();
//   cout << "AliFemtoBasicEventCut:: mult:       " << fEventMult[0] << " < " << mult << " < " << fEventMult[1] << endl;
//   cout << "AliFemtoBasicEventCut:: VertexZPos: " << fVertZPos[0] << " < " << vertexZPos << " < " << fVertZPos[1] << endl;
//   cout << "AliFemtoBasicEventCut:: VertexZErr: " << event->PrimVertCov()[4] << endl;
  bool goodEvent =
    ((mult >= fEventMult[0]) && 
     (mult <= fEventMult[1]) && 
     (vertexZPos > fVertZPos[0]) &&
     (vertexZPos < fVertZPos[1]) &&
     ((!fAcceptBadVertex) || (event->ZDCParticipants() > 1.0)) &&
     ((!fAcceptOnlyPhysics) || (event->IsCollisionCandidate())) &&
     ((!fSelectTrigger) || (event->TriggerCluster() == fSelectTrigger)));
  goodEvent ? fNEventsPassed++ : fNEventsFailed++ ;
//   cout << "AliFemtoBasicEventCut:: return : " << goodEvent << endl;
//     (fAcceptBadVertex || (event->PrimVertCov()[4] > -1000.0)) &&
  return (goodEvent);
}
//------------------------------
AliFemtoString AliFemtoBasicEventCut::Report(){
  // Prepare report
  string stemp;
  char ctemp[100];
  snprintf(ctemp , 100, "\nMultiplicity:\t %d-%d",fEventMult[0],fEventMult[1]);
  stemp = ctemp;
  snprintf(ctemp , 100, "\nVertex Z-position:\t %E-%E",fVertZPos[0],fVertZPos[1]);
  stemp += ctemp;
  snprintf(ctemp , 100, "\nNumber of events which passed:\t%ld  Number which failed:\t%ld",fNEventsPassed,fNEventsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;
}
void AliFemtoBasicEventCut::SetAcceptBadVertex(bool b)
{
  fAcceptBadVertex = b;
}
bool AliFemtoBasicEventCut::GetAcceptBadVertex()
{
  return fAcceptBadVertex;
}
void AliFemtoBasicEventCut::SetAcceptOnlyPhysics(bool b)
{
  fAcceptOnlyPhysics = b;
}
bool AliFemtoBasicEventCut::GetAcceptOnlyPhysics()
{
  return fAcceptOnlyPhysics;
}
