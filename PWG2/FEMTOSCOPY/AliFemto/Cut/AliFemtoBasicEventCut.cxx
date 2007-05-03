/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   A simple event-wise cut that selects on multiplicity and z-position
 *   of primary vertex           
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.7  2000/02/18 21:27:10  laue
 * franksTrackCut changed. If mCharge is set to '0' there will be no cut
 * on charge. This is important for front-loaded cuts.
 *
 * Revision 1.6  2000/01/25 17:35:02  laue
 * I. In order to run the stand alone version of the AliFemtoMaker the following
 * changes have been done:
 * a) all ClassDefs and ClassImps have been put into #ifdef __ROOT__ statements
 * b) unnecessary includes of StMaker.h have been removed
 * c) the subdirectory AliFemtoMaker/doc/Make has been created including everything
 * needed for the stand alone version
 *
 * II. To reduce the amount of compiler warning
 * a) some variables have been type casted
 * b) some destructors have been declared as virtual
 *
 * Revision 1.5  1999/10/15 01:57:03  lisa
 * Important enhancement of AliFemtoMaker - implement Franks CutMonitors
 * ----------------------------------------------------------
 * This means 3 new files in Infrastructure area (CutMonitor),
 * several specific CutMonitor classes in the Cut area
 * and a new base class in the Base area (AliFemtoCutMonitor).
 * This means also changing all Cut Base class header files from .h to .h
 * so we have access to CutMonitor methods from Cint command line.
 * This last means
 * 1) files which include these header files are slightly modified
 * 2) a side benefit: the TrackCuts and V0Cuts no longer need
 * a SetMass() implementation in each Cut class, which was stupid.
 * Also:
 * -----
 * Include Franks AliFemtoAssociationReader
 * ** None of these changes should affect any user **
 *
 * Revision 1.4  1999/07/24 16:24:20  lisa
 * adapt AliFemtoMaker to dev version of library - solaris still gives problems with strings
 *
 * Revision 1.3  1999/07/19 14:24:04  hardtke
 * modifications to implement uDST
 *
 * Revision 1.2  1999/07/06 22:33:21  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:56  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "Cut/AliFemtoBasicEventCut.h"
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
