////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoSphericityEventCutKK1 - the basic cut for events.                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoSphericityEventCutKK1_H
#define AliFemtoSphericityEventCutKK1_H

#include "AliFemtoEventCut.h"

class AliFemtoSphericityEventCutKK1 : public AliFemtoEventCut {

public:

  AliFemtoSphericityEventCutKK1();
  AliFemtoSphericityEventCutKK1(const AliFemtoSphericityEventCutKK1& c);
  virtual ~AliFemtoSphericityEventCutKK1();
  AliFemtoSphericityEventCutKK1& operator=(const AliFemtoSphericityEventCutKK1& c);

  void SetEventMult(const int lo,const int hi);
  void SetVertZPos(const float lo, const float hi);
  void SetAcceptBadVertex(bool b);
  int NEventsPassed() const;
  int NEventsFailed() const;
  bool GetAcceptBadVertex();
  bool GetAcceptOnlyPhysics() {return kFALSE;} // Not implemented
  void SetStMin(double stMin );
  void SetStMax(double stMax );
  void SetTriggerSelection(int trig);

  void SetEPVZERO(const float lo, const float hi);

  virtual AliFemtoString Report();
  virtual bool Pass(const AliFemtoEvent* event);

  virtual AliFemtoEventCut* Clone() const;

private:   // here are the quantities I want to cut on...

  int fEventMult[2];      // range of multiplicity
  float fVertZPos[2];     // range of z-position of vertex
  float fPsiEP[2];        // range of vzero ep angle
  bool fAcceptBadVertex;  // Set to true to accept events with bad vertex
  long fNEventsPassed;    // Number of events checked by this cut that passed
  long fNEventsFailed;    // Number of events checked by this cut that failed
  bool fAcceptOnlyPhysics;// Accept only physics events
  double fStCutMin;       // transverse sphericity minimum
  double fStCutMax;       // transverse sphericity maximum
  int  fSelectTrigger;    // If set, only given trigger will be selected

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoSphericityEventCutKK1, 1);
  /// \endcond
#endif

};

inline void AliFemtoSphericityEventCutKK1::SetEventMult(const int lo, const int hi){fEventMult[0]=lo; fEventMult[1]=hi;}
inline void AliFemtoSphericityEventCutKK1::SetVertZPos(const float lo, const float hi){fVertZPos[0]=lo; fVertZPos[1]=hi;}
inline void AliFemtoSphericityEventCutKK1::SetEPVZERO(const float lo, const float hi){fPsiEP[0]=lo; fPsiEP[1]=hi;}
inline int  AliFemtoSphericityEventCutKK1::NEventsPassed() const {return fNEventsPassed;}
inline int  AliFemtoSphericityEventCutKK1::NEventsFailed() const {return fNEventsFailed;}
inline void  AliFemtoSphericityEventCutKK1::SetStMin(double stMin) {fStCutMin=stMin;}
inline void  AliFemtoSphericityEventCutKK1::SetStMax(double stMax) {fStCutMax=stMax;}
inline void AliFemtoSphericityEventCutKK1::SetTriggerSelection(int trig) { fSelectTrigger = trig; }
inline AliFemtoEventCut* AliFemtoSphericityEventCutKK1::Clone() const { AliFemtoSphericityEventCutKK1* c = new AliFemtoSphericityEventCutKK1(*this); return c;}
inline AliFemtoSphericityEventCutKK1::AliFemtoSphericityEventCutKK1(const AliFemtoSphericityEventCutKK1& c):
  AliFemtoEventCut(c),
  fAcceptBadVertex(c.fAcceptBadVertex),
  fNEventsPassed(0),
  fNEventsFailed(0),
  fAcceptOnlyPhysics(c.fAcceptOnlyPhysics),
  fStCutMin(c.fStCutMin),
  fStCutMax(c.fStCutMax),
  fSelectTrigger(c.fSelectTrigger)
{
  fEventMult[0] = c.fEventMult[0];
  fEventMult[1] = c.fEventMult[1];
  fVertZPos[0] = c.fVertZPos[0];
  fVertZPos[1] = c.fVertZPos[1];
  fPsiEP[0] = c.fPsiEP[0];
  fPsiEP[1] = c.fPsiEP[1];
}

inline AliFemtoSphericityEventCutKK1& AliFemtoSphericityEventCutKK1::operator=(const AliFemtoSphericityEventCutKK1& c)
{
  if (this != &c) {
    AliFemtoEventCut::operator=(c);
    fEventMult[0] = c.fEventMult[0];
    fEventMult[1] = c.fEventMult[1];
    fVertZPos[0] = c.fVertZPos[0];
    fVertZPos[1] = c.fVertZPos[1];
    fPsiEP[0] = c.fPsiEP[0];
    fPsiEP[1] = c.fPsiEP[1];
    fAcceptBadVertex = c.fAcceptBadVertex;
    // fNEventsPassed = c.fNEventsPassed;
    // fNEventsFailed = c.fNEventsFailed;
    fAcceptOnlyPhysics = c.fAcceptOnlyPhysics;
    fStCutMin = c.fStCutMin;
    fStCutMax = c.fStCutMax;
    fSelectTrigger = c.fSelectTrigger;
  }

  return *this;
}


#endif
