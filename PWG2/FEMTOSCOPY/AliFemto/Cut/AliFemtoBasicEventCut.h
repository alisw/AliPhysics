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
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.5  2000/03/23 22:57:28  laue
 * Clone() function implemented
 *
 * Revision 1.4  2000/01/25 17:35:02  laue
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
 * Revision 1.3  1999/10/15 01:57:04  lisa
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
 * Revision 1.2  1999/07/06 22:33:21  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:56  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#ifndef AliFemtoBasicEventCut_hh
#define AliFemtoBasicEventCut_hh

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "Base/AliFemtoEventCut.h"

class AliFemtoBasicEventCut : public AliFemtoEventCut {

public:

  AliFemtoBasicEventCut();
  AliFemtoBasicEventCut(AliFemtoBasicEventCut&);
  //~AliFemtoBasicEventCut();

  void SetEventMult(const int& lo,const int& hi);
  void SetVertZPos(const float& lo, const float& hi);
  int NEventsPassed();
  int NEventsFailed();

  virtual AliFemtoString Report();
  virtual bool Pass(const AliFemtoEvent*);

  AliFemtoBasicEventCut* Clone();

private:   // here are the quantities I want to cut on...

  int fEventMult[2];      // range of multiplicity
  float fVertZPos[2];     // range of z-position of vertex

  long fNEventsPassed;
  long fNEventsFailed;

#ifdef __ROOT__
  ClassDef(AliFemtoBasicEventCut, 1)
#endif

};

inline void AliFemtoBasicEventCut::SetEventMult(const int& lo, const int& hi){fEventMult[0]=lo; fEventMult[1]=hi;}
inline void AliFemtoBasicEventCut::SetVertZPos(const float& lo, const float& hi){fVertZPos[0]=lo; fVertZPos[1]=hi;}
inline int  AliFemtoBasicEventCut::NEventsPassed() {return fNEventsPassed;}
inline int  AliFemtoBasicEventCut::NEventsFailed() {return fNEventsFailed;}
inline AliFemtoBasicEventCut* AliFemtoBasicEventCut::Clone() { AliFemtoBasicEventCut* c = new AliFemtoBasicEventCut(*this); return c;}
inline AliFemtoBasicEventCut::AliFemtoBasicEventCut(AliFemtoBasicEventCut& c) : AliFemtoEventCut(c) {
  fEventMult[0] = c.fEventMult[0];
  fEventMult[1] = c.fEventMult[1];
  fVertZPos[0] = c.fVertZPos[0];
  fVertZPos[1] = c.fVertZPos[1];
  fNEventsPassed = 0;
  fNEventsFailed = 0;
}


#endif
