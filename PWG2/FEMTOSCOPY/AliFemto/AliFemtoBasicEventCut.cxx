////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoBasicEventCut - the basic cut for events.                          //
// Only cuts on event multiplicity and z-vertex position                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoBasicEventCut.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoBasicEventCut)
#endif

AliFemtoBasicEventCut::AliFemtoBasicEventCut() :
  fNEventsPassed(0), fNEventsFailed(0)
{
  /* no-op */
} 
//------------------------------
//AliFemtoBasicEventCut::~AliFemtoBasicEventCut(){
//  /* noop */
//}
//------------------------------
bool AliFemtoBasicEventCut::Pass(const AliFemtoEvent* event){
  int mult =  event->NumberOfTracks();
  double VertexZPos = event->PrimVertPos().z();
  cout << "AliFemtoBasicEventCut:: mult:       " << fEventMult[0] << " < " << mult << " < " << fEventMult[1] << endl;
  cout << "AliFemtoBasicEventCut:: VertexZPos: " << fVertZPos[0] << " < " << VertexZPos << " < " << fVertZPos[1] << endl;
  bool goodEvent =
    ((mult > fEventMult[0]) && 
     (mult < fEventMult[1]) && 
     (VertexZPos > fVertZPos[0]) &&
     (VertexZPos < fVertZPos[1]));
  goodEvent ? fNEventsPassed++ : fNEventsFailed++ ;
  cout << "AliFemtoBasicEventCut:: return : " << goodEvent << endl;
  return (goodEvent);
}
//------------------------------
AliFemtoString AliFemtoBasicEventCut::Report(){
  string Stemp;
  char Ctemp[100];
  sprintf(Ctemp,"\nMultiplicity:\t %d-%d",fEventMult[0],fEventMult[1]);
  Stemp = Ctemp;
  sprintf(Ctemp,"\nVertex Z-position:\t %E-%E",fVertZPos[0],fVertZPos[1]);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nNumber of events which passed:\t%ld  Number which failed:\t%ld",fNEventsPassed,fNEventsFailed);
  Stemp += Ctemp;
  AliFemtoString returnThis = Stemp;
  return returnThis;
}
