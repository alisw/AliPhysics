////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventCutEstimators - the basic cut for events.                          //
// Only cuts on event multiplicity and z-vertex position                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOEVENTCUTESTIMATORS_H
#define ALIFEMTOEVENTCUTESTIMATORS_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoEventCut.h"
#include "AliFemtoEventReaderESDChain.h"

class AliFemtoEventCutEstimators : public AliFemtoEventCut {

public:

  AliFemtoEventCutEstimators();
  AliFemtoEventCutEstimators(AliFemtoEventCutEstimators& c);
  virtual ~AliFemtoEventCutEstimators();

  void SetMultEst1Range(const unsigned short &lo, const unsigned short &hi);
  void SetMultEst2Range(const unsigned short &lo, const unsigned short &hi);
  void SetMultEst3Range(const unsigned short &lo, const unsigned short &hi);
  
  void SetCentEst1Range(const float &lo, const float &hi);
  void SetCentEst2Range(const float &lo, const float &hi);
  void SetCentEst3Range(const float &lo, const float &hi);
  void SetCentEst4Range(const float &lo, const float &hi);

  void SetVertZPos(const float& lo, const float& hi);
  int NEventsPassed() const;
  int NEventsFailed() const;

  virtual AliFemtoString Report();
  virtual bool Pass(const AliFemtoEvent* event);

  AliFemtoEventCutEstimators* Clone();

private:   // here are the quantities I want to cut on...

  unsigned short fEventMultEst1[2];      // range of multiplicity
  unsigned short fEventMultEst2[2];      // range of multiplicity
  unsigned short fEventMultEst3[2];      // range of multiplicity
  unsigned char fUseMultEst1;  // if 1 cut on Mult Est 1
  unsigned char fUseMultEst2;  // if 1 cut on Mult Est 2
  unsigned char fUseMultEst3;  // if 1 cut on Mult Est 3
  
  float fEventCentEst1[2];      // range of multiplicity
  float fEventCentEst2[2];      // range of multiplicity
  float fEventCentEst3[2];      // range of multiplicity
  float fEventCentEst4[2];      // range of multiplicity
  unsigned char fUseCentEst1;  // if 1 cut on Mult Est 1
  unsigned char fUseCentEst2;  // if 1 cut on Mult Est 2
  unsigned char fUseCentEst3;  // if 1 cut on Mult Est 3
  unsigned char fUseCentEst4;  // if 1 cut on Mult Est 4

  float fVertZPos[2];     // range of z-position of vertex
  long fNEventsPassed;    // Number of events checked by this cut that passed
  long fNEventsFailed;    // Number of events checked by this cut that failed

#ifdef __ROOT__
  ClassDef(AliFemtoEventCutEstimators, 1)
#endif

};

inline void AliFemtoEventCutEstimators::SetMultEst1Range(const unsigned short& lo, const unsigned short& hi){fEventMultEst1[0]=lo; fEventMultEst1[1]=hi; fUseMultEst1=1;}
inline void AliFemtoEventCutEstimators::SetMultEst2Range(const unsigned short& lo, const unsigned short& hi){fEventMultEst2[0]=lo; fEventMultEst2[1]=hi; fUseMultEst2=1;}
inline void AliFemtoEventCutEstimators::SetMultEst3Range(const unsigned short& lo, const unsigned short& hi){fEventMultEst3[0]=lo; fEventMultEst3[1]=hi; fUseMultEst3=1;}
inline void AliFemtoEventCutEstimators::SetCentEst1Range(const float& lo, const float& hi){fEventCentEst1[0]=lo; fEventCentEst1[1]=hi; fUseCentEst1=1;}
inline void AliFemtoEventCutEstimators::SetCentEst2Range(const float& lo, const float& hi){fEventCentEst2[0]=lo; fEventCentEst2[1]=hi; fUseCentEst2=1;}
inline void AliFemtoEventCutEstimators::SetCentEst3Range(const float& lo, const float& hi){fEventCentEst3[0]=lo; fEventCentEst3[1]=hi; fUseCentEst3=1;}
inline void AliFemtoEventCutEstimators::SetCentEst4Range(const float& lo, const float& hi){fEventCentEst4[0]=lo; fEventCentEst4[1]=hi; fUseCentEst4=1;}
inline void AliFemtoEventCutEstimators::SetVertZPos(const float& lo, const float& hi){fVertZPos[0]=lo; fVertZPos[1]=hi;}
inline int  AliFemtoEventCutEstimators::NEventsPassed() const {return fNEventsPassed;}
inline int  AliFemtoEventCutEstimators::NEventsFailed() const {return fNEventsFailed;}
inline AliFemtoEventCutEstimators* AliFemtoEventCutEstimators::Clone() { AliFemtoEventCutEstimators* c = new AliFemtoEventCutEstimators(*this); return c;}
inline AliFemtoEventCutEstimators::AliFemtoEventCutEstimators(AliFemtoEventCutEstimators& c) : 
  AliFemtoEventCut(c), 
  fUseMultEst1(0), fUseMultEst2(0), fUseMultEst3(0), 
  fUseCentEst1(0), fUseCentEst2(0), fUseCentEst3(0), fUseCentEst4(0),
  fNEventsPassed(0), fNEventsFailed(0) {
  fEventMultEst1[0] = c.fEventMultEst1[0];  fEventMultEst1[1] = c.fEventMultEst1[1];
  fEventMultEst2[0] = c.fEventMultEst2[0];  fEventMultEst2[1] = c.fEventMultEst2[1];
  fEventMultEst3[0] = c.fEventMultEst3[0];  fEventMultEst3[1] = c.fEventMultEst3[1];
  fEventCentEst1[0] = c.fEventCentEst1[0];  fEventCentEst1[1] = c.fEventCentEst1[1];
  fEventCentEst2[0] = c.fEventCentEst2[0];  fEventCentEst2[1] = c.fEventCentEst2[1];
  fEventCentEst3[0] = c.fEventCentEst3[0];  fEventCentEst3[1] = c.fEventCentEst3[1];
  fEventCentEst4[0] = c.fEventCentEst4[0];  fEventCentEst4[1] = c.fEventCentEst4[1];
  fVertZPos[0] = c.fVertZPos[0];
  fVertZPos[1] = c.fVertZPos[1];
}


#endif
