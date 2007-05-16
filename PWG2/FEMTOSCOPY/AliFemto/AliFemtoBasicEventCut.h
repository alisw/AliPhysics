////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoBasicEventCut - the basic cut for events.                          //
// Only cuts on event multiplicity and z-vertex position                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoBasicEventCut_hh
#define AliFemtoBasicEventCut_hh

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoEventCut.h"

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
inline AliFemtoBasicEventCut::AliFemtoBasicEventCut(AliFemtoBasicEventCut& c) : AliFemtoEventCut(c), fNEventsPassed(0), fNEventsFailed(0) {
  fEventMult[0] = c.fEventMult[0];
  fEventMult[1] = c.fEventMult[1];
  fVertZPos[0] = c.fVertZPos[0];
  fVertZPos[1] = c.fVertZPos[1];
}


#endif
