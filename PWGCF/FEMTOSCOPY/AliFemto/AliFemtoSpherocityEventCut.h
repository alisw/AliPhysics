////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoSpherocityEventCut - the basic cut for events.                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoSpherocityEventCUT_H
#define AliFemtoSpherocityEventCUT_H

#include "AliFemtoEventCut.h"

class AliFemtoSpherocityEventCut : public AliFemtoEventCut {

public:

  AliFemtoSpherocityEventCut();
  AliFemtoSpherocityEventCut(const AliFemtoSpherocityEventCut& c);
  virtual ~AliFemtoSpherocityEventCut();
  AliFemtoSpherocityEventCut& operator=(const AliFemtoSpherocityEventCut& c);

  void SetEventMult(const int& lo,const int& hi);
  void SetVertZPos(const float& lo, const float& hi);
  void SetAcceptBadVertex(bool b);
  int NEventsPassed() const;
  int NEventsFailed() const;
  bool GetAcceptBadVertex();
  bool GetAcceptOnlyPhysics() {return kFALSE;} // Not implemented
  void SetSoMin(double soMin );
  void SetSoMax(double soMax );
  void SetTriggerSelection(int trig);

  void SetEPVZERO(const float& lo, const float& hi);

  virtual AliFemtoString Report();
  virtual bool Pass(const AliFemtoEvent* event);

  virtual AliFemtoEventCut* Clone() const;

private:   // here are the quantities I want to cut on...

  int fEventMult[2];       ///< range of multiplicity
  float fVertZPos[2];      ///< range of z-position of vertex
  float fPsiEP[2];         ///< range of vzero ep angle
  bool fAcceptBadVertex;   ///< Set to true to accept events with bad vertex
  long fNEventsPassed;     ///< Number of events checked by this cut that passed
  long fNEventsFailed;     ///< Number of events checked by this cut that failed
  bool fAcceptOnlyPhysics; ///< Accept only physics events
  double fSoCutMin;        ///< transverse sphericity minimum
  double fSoCutMax;        ///< transverse sphericity maximum
  int  fSelectTrigger;     ///< If set, only given trigger will be selected

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoSpherocityEventCut, 1);
  /// \endcond
#endif

};

inline void AliFemtoSpherocityEventCut::SetEventMult(const int& lo, const int& hi){fEventMult[0]=lo; fEventMult[1]=hi;}
inline void AliFemtoSpherocityEventCut::SetVertZPos(const float& lo, const float& hi){fVertZPos[0]=lo; fVertZPos[1]=hi;}
inline void AliFemtoSpherocityEventCut::SetEPVZERO(const float& lo, const float& hi){fPsiEP[0]=lo; fPsiEP[1]=hi;}
inline int  AliFemtoSpherocityEventCut::NEventsPassed() const {return fNEventsPassed;}
inline int  AliFemtoSpherocityEventCut::NEventsFailed() const {return fNEventsFailed;}
inline void AliFemtoSpherocityEventCut::SetSoMin(double soMin) {fSoCutMin=soMin;}
inline void AliFemtoSpherocityEventCut::SetSoMax(double soMax) {fSoCutMax=soMax;}
inline void AliFemtoSpherocityEventCut::SetTriggerSelection(int trig) { fSelectTrigger = trig; }
inline AliFemtoEventCut* AliFemtoSpherocityEventCut::Clone() const { AliFemtoSpherocityEventCut* c = new AliFemtoSpherocityEventCut(*this); return c;}
inline AliFemtoSpherocityEventCut::AliFemtoSpherocityEventCut(const AliFemtoSpherocityEventCut& c):
  AliFemtoEventCut(c),
  fAcceptBadVertex(c.fAcceptBadVertex),
  fNEventsPassed(0),
  fNEventsFailed(0),
  fAcceptOnlyPhysics(c.fAcceptOnlyPhysics),
  fSoCutMin(c.fSoCutMin),
  fSoCutMax(c.fSoCutMax),
  fSelectTrigger(c.fSelectTrigger)
{
  fEventMult[0] = c.fEventMult[0];
  fEventMult[1] = c.fEventMult[1];
  fVertZPos[0] = c.fVertZPos[0];
  fVertZPos[1] = c.fVertZPos[1];
  fPsiEP[0] = c.fPsiEP[0];
  fPsiEP[1] = c.fPsiEP[1];
}

inline AliFemtoSpherocityEventCut& AliFemtoSpherocityEventCut::operator=(const AliFemtoSpherocityEventCut& c) {
  if (this != &c) {
    AliFemtoEventCut::operator=(c);
    fEventMult[0] = c.fEventMult[0];
    fEventMult[1] = c.fEventMult[1];
    fVertZPos[0] = c.fVertZPos[0];
    fVertZPos[1] = c.fVertZPos[1];
    fPsiEP[0] = c.fPsiEP[0];
    fPsiEP[1] = c.fPsiEP[1];

    fAcceptBadVertex = c.fAcceptBadVertex;
    fNEventsPassed = 0;
    fNEventsFailed = 0;
    fAcceptOnlyPhysics = c.fAcceptOnlyPhysics;
    fSoCutMin = c.fSoCutMin;
    fSoCutMax = c.fSoCutMax;
    fSelectTrigger = c.fSelectTrigger;
  }

  return *this;
}


#endif
