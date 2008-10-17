////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoQAEventCut - the basic cut to check QA for event cuts.             //
// Only cuts on event multiplicity and z-vertex position                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOQAEVENTCUT_H
#define ALIFEMTOQAEVENTCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoEventCut.h"

class AliFemtoQAEventCut : public AliFemtoEventCut {

public:

  AliFemtoQAEventCut();
  AliFemtoQAEventCut(AliFemtoQAEventCut& c);
  virtual ~AliFemtoQAEventCut();

  void SetEventMult(const int& lo,const int& hi);
  void SetEventMultQASwitch(const bool Switch);
  void SetEventMultQAExclusionZone(const int& lo, const int& hi);
  void SetVertZPos(const float& lo, const float& hi);
  void SetEventZPosQASwitch(const bool Switch);
  void SetEventZPosQAExclusionZone(const float& lo, const float& hi);
  void SetAcceptBadVertex(bool b);
  int NEventsPassed() const;
  int NEventsFailed() const;
  bool GetAcceptBadVertex();

  virtual AliFemtoString Report();
  virtual bool Pass(const AliFemtoEvent* event);

  AliFemtoQAEventCut* Clone();

private:   // here are the quantities I want to cut on...

  int   fEventMult[2];     // range of multiplicity
  float fVertZPos[2];      // range of z-position of vertex
  bool  fAcceptBadVertex;  // Set to true to accept events with bad vertex
  long  fNEventsPassed;    // Number of events checked by this cut that passed
  long  fNEventsFailed;    // Number of events checked by this cut that failed
  
  int   fHighOrLowSwitch;             // if 1, then previous hbtEvent was high; if -1, then previous event was low.
  bool  fEventMultQASwitch;           // Turn on multiplicity exclusion zone (true=on)
  int   fEventMultQAExclusionZone[2]; // Set limits of the multiplicity exclusion zone
  bool  fEventZPosQASwitch;           // Turn on Zpos exclusion zone (true=on)
  float fEventZPosQAExclusionZone[2]; // Set limits of the Zpos exclusion zone

#ifdef __ROOT__
  ClassDef(AliFemtoQAEventCut, 1)
#endif

};

inline void AliFemtoQAEventCut::SetEventMult(const int& lo, const int& hi){fEventMult[0]=lo; fEventMult[1]=hi;}
inline void AliFemtoQAEventCut::SetEventMultQASwitch(const bool Switch) { fEventMultQASwitch = Switch; }
inline void AliFemtoQAEventCut::SetEventMultQAExclusionZone(const int& lo, const int& hi) { fEventMultQAExclusionZone[0]=lo; fEventMultQAExclusionZone[1]=hi; }
inline void AliFemtoQAEventCut::SetVertZPos(const float& lo, const float& hi){fVertZPos[0]=lo; fVertZPos[1]=hi;}
inline void AliFemtoQAEventCut::SetEventZPosQASwitch(const bool Switch) { fEventZPosQASwitch = Switch; }
inline void AliFemtoQAEventCut::SetEventZPosQAExclusionZone(const float& lo, const float& hi)  { fEventZPosQAExclusionZone[0]=lo; fEventZPosQAExclusionZone[1]=hi; }
inline int  AliFemtoQAEventCut::NEventsPassed() const {return fNEventsPassed;}
inline int  AliFemtoQAEventCut::NEventsFailed() const {return fNEventsFailed;}
inline AliFemtoQAEventCut* AliFemtoQAEventCut::Clone() { AliFemtoQAEventCut* c = new AliFemtoQAEventCut(*this); return c;}
inline AliFemtoQAEventCut::AliFemtoQAEventCut(AliFemtoQAEventCut& c) : AliFemtoEventCut(c), fAcceptBadVertex(kFALSE), fNEventsPassed(0), fNEventsFailed(0), fHighOrLowSwitch(0), fEventMultQASwitch(kFALSE), fEventZPosQASwitch(kFALSE) {
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


#endif
