////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoBasicEventCut - the basic cut for events.                          //
// Only cuts on event multiplicity and z-vertex position                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOBASICEVENTCUT_H
#define ALIFEMTOBASICEVENTCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoEventCut.h"

class AliFemtoBasicEventCut : public AliFemtoEventCut {

public:

  AliFemtoBasicEventCut();
  AliFemtoBasicEventCut(AliFemtoBasicEventCut& c);
  virtual ~AliFemtoBasicEventCut();
  AliFemtoBasicEventCut& operator=(AliFemtoBasicEventCut& c);

  void SetEventMult(const int& lo,const int& hi);
  void SetVertZPos(const float& lo, const float& hi);
  void SetAcceptBadVertex(bool b);
  int NEventsPassed() const;
  int NEventsFailed() const;
  bool GetAcceptBadVertex();
  bool GetAcceptOnlyPhysics();
  void SetTriggerSelection(int trig);

  void SetEPVZERO(const float& lo, const float& hi);

  virtual AliFemtoString Report();
  virtual bool Pass(const AliFemtoEvent* event);

  AliFemtoBasicEventCut* Clone();

private:   // here are the quantities I want to cut on...

  int fEventMult[2];      // range of multiplicity
  float fVertZPos[2];     // range of z-position of vertex
  float fPsiEP[2];     // range of vzero ep angle
  bool fAcceptBadVertex;  // Set to true to accept events with bad vertex
  long fNEventsPassed;    // Number of events checked by this cut that passed
  long fNEventsFailed;    // Number of events checked by this cut that failed
  bool fAcceptOnlyPhysics;// Accept only physics events
  int  fSelectTrigger;    // If set, only given trigger will be selected

#ifdef __ROOT__
  ClassDef(AliFemtoBasicEventCut, 1)
#endif

};

inline void AliFemtoBasicEventCut::SetEventMult(const int& lo, const int& hi){fEventMult[0]=lo; fEventMult[1]=hi;}
inline void AliFemtoBasicEventCut::SetVertZPos(const float& lo, const float& hi){fVertZPos[0]=lo; fVertZPos[1]=hi;}
inline void AliFemtoBasicEventCut::SetEPVZERO(const float& lo, const float& hi){fPsiEP[0]=lo; fPsiEP[1]=hi;}
inline int  AliFemtoBasicEventCut::NEventsPassed() const {return fNEventsPassed;}
inline int  AliFemtoBasicEventCut::NEventsFailed() const {return fNEventsFailed;}
inline void AliFemtoBasicEventCut::SetTriggerSelection(int trig) { fSelectTrigger = trig; }
inline AliFemtoBasicEventCut* AliFemtoBasicEventCut::Clone() { AliFemtoBasicEventCut* c = new AliFemtoBasicEventCut(*this); return c;}
inline AliFemtoBasicEventCut::AliFemtoBasicEventCut(AliFemtoBasicEventCut& c) : AliFemtoEventCut(c), fAcceptBadVertex(false), fNEventsPassed(0), fNEventsFailed(0), fAcceptOnlyPhysics(false), fSelectTrigger(0) {
  fEventMult[0] = c.fEventMult[0];
  fEventMult[1] = c.fEventMult[1];
  fVertZPos[0] = c.fVertZPos[0];
  fVertZPos[1] = c.fVertZPos[1];
  fPsiEP[0] = c.fPsiEP[0];
  fPsiEP[1] = c.fPsiEP[1];
}

inline AliFemtoBasicEventCut& AliFemtoBasicEventCut::operator=(AliFemtoBasicEventCut& c) {   
  if (this != &c) {
    AliFemtoEventCut::operator=(c);
    fEventMult[0] = c.fEventMult[0];
    fEventMult[1] = c.fEventMult[1];
    fVertZPos[0] = c.fVertZPos[0];
    fVertZPos[1] = c.fVertZPos[1];
    fPsiEP[0] = c.fPsiEP[0];
    fPsiEP[1] = c.fPsiEP[1];
  }

  return *this;
}


#endif
