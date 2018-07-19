///
/// \file AliFemtoUser/AliFemtoQAEventCut.cxx
///

#include "AliFemtoQAEventCut.h"
#include <TString.h>

#ifdef __ROOT__
  ClassImp(AliFemtoQAEventCut)
#endif

  AliFemtoQAEventCut::AliFemtoQAEventCut() :
    AliFemtoEventCut(),
    fEventMult(),
    fVertZPos(),
    fAcceptBadVertex(false),
    fNEventsPassed(0),
    fNEventsFailed(0),
    fHighOrLowSwitch(0),
    fEventMultQASwitch(kFALSE),
    fEventZPosQASwitch(kFALSE)
  {
    // Default constructor
    fEventMult[0] = 0;
    fEventMult[1] = 100000;
    fVertZPos[0] = -100.0;
    fVertZPos[1] = 100.0;

    fHighOrLowSwitch = 1;
    fEventMultQASwitch = false;
    fEventZPosQASwitch = false;
    fEventMultQAExclusionZone[0] = 0;
    fEventMultQAExclusionZone[1] = 100000;
    fEventZPosQAExclusionZone[0] = -100.0;
    fEventZPosQAExclusionZone[1] = 100.0;

  }
  //------------------------------
  AliFemtoQAEventCut::~AliFemtoQAEventCut(){
    // Default destructor
  }
  //------------------------------
  AliFemtoQAEventCut& AliFemtoQAEventCut::operator=(const AliFemtoQAEventCut& c)
  {
    if (this != &c) {
      AliFemtoEventCut::operator=(c);

      fEventMult[0] = c.fEventMult[0];
      fEventMult[1] = c.fEventMult[1];
      fVertZPos[0] = c.fVertZPos[0];
      fVertZPos[1] = c.fVertZPos[1];

      fHighOrLowSwitch = c.fHighOrLowSwitch;
      fEventMultQASwitch = c.fEventMultQASwitch;
      fEventZPosQASwitch = c.fEventZPosQASwitch;
      fEventMultQAExclusionZone[0] = c.fEventMultQAExclusionZone[0];
      fEventMultQAExclusionZone[1] = c.fEventMultQAExclusionZone[1];
      fEventZPosQAExclusionZone[0] = c.fEventZPosQAExclusionZone[0];
      fEventZPosQAExclusionZone[1] = c.fEventZPosQAExclusionZone[1];
    }

    return *this;
  }
  //------------------------------
  bool AliFemtoQAEventCut::Pass(const AliFemtoEvent* event)
  {
    // Pass events if they fall within the multiplicity and z-vertex
    // position range. If QA cutting on quantity, pass if outside
    // exclusion zone between low and high cut values. Fail otherwise.
    int mult =  event->NumberOfTracks();
  double vertexZPos = event->PrimVertPos().z();
  cout << "AliFemtoQAEventCut:: mult:       " << fEventMult[0] << " < " << mult << " < " << fEventMult[1] << endl;
  cout << "AliFemtoQAEventCut:: VertexZPos: " << fVertZPos[0] << " < " << vertexZPos << " < " << fVertZPos[1] << endl;

  bool goodEvent;

  if (fEventMultQASwitch) {
    goodEvent =
      ( (((mult < fEventMultQAExclusionZone[0]) && (fHighOrLowSwitch > 0))  ||
	 ((mult > fEventMultQAExclusionZone[1]) && (fHighOrLowSwitch < 0))) &&
      (mult > fEventMult[0]) &&
      (mult < fEventMult[1]) &&
      (vertexZPos > fVertZPos[0]) &&
      (vertexZPos < fVertZPos[1]) &&
      (fAcceptBadVertex || (event->PrimVertCov()[4] > -1000.0)));
  }
  else if (fEventZPosQASwitch) {
    goodEvent =
      ((((vertexZPos < fEventZPosQAExclusionZone[0]) && (fHighOrLowSwitch > 0))  ||
	((vertexZPos > fEventZPosQAExclusionZone[1]) && (fHighOrLowSwitch < 0))) &&
      (mult > fEventMult[0]) &&
      (mult < fEventMult[1]) &&
      (vertexZPos > fVertZPos[0]) &&
      (vertexZPos < fVertZPos[1]) &&
      (fAcceptBadVertex || (event->PrimVertCov()[4] > -1000.0)));
  }
  else {
  goodEvent =
    ((mult > fEventMult[0]) &&
     (mult < fEventMult[1]) &&
     (vertexZPos > fVertZPos[0]) &&
     (vertexZPos < fVertZPos[1]) &&
     (fAcceptBadVertex || (event->PrimVertCov()[4] > -1000.0)));
  }

  if (goodEvent) fHighOrLowSwitch *= -1;
  goodEvent ? fNEventsPassed++ : fNEventsFailed++ ;
  //cout << "AliFemtoQAEventCut:: return : " << goodEvent << endl;
  return (goodEvent);
}
//------------------------------
AliFemtoString AliFemtoQAEventCut::Report()
{
  // Prepare report
  TString report;
  report += TString::Format("\nMultiplicity:\t %d-%d",fEventMult[0],fEventMult[1])
          + TString::Format("\nVertex Z-position:\t %E-%E",fVertZPos[0],fVertZPos[1])
          + TString::Format("\nNumber of events which passed:\t%ld  Number which failed:\t%ld\n",fNEventsPassed,fNEventsFailed);
  return AliFemtoString(report.Data());
}

void AliFemtoQAEventCut::SetAcceptBadVertex(bool b)
{
  fAcceptBadVertex = b;
}
bool AliFemtoQAEventCut::GetAcceptBadVertex()
{
  return fAcceptBadVertex;
}
