/***************************************************************************
 *
 * $Id:
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *    a simple particle cut that selects on phasespace, #hits, DCA, and PID          
 *
 ***************************************************************************
 *
 * $Log:
 **************************************************************************/

#ifndef AliFemtoBasicTrackCut_hh
#define AliFemtoBasicTrackCut_hh

//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "Base/AliFemtoTrackCut.h"

class AliFemtoBasicTrackCut : public AliFemtoTrackCut {

public:

  AliFemtoBasicTrackCut();
  //~mikesTrackCut();

  virtual bool Pass(const AliFemtoTrack*);

  virtual AliFemtoString Report();


  void SetNSigmaPion(const float& lo, const float& hi);
  void SetNSigmaKaon(const float& lo, const float& hi);
  void SetNSigmaProton(const float& lo, const float& hi);

  void SetNHits(const int& lo, const int& hi);
  void SetPt(const float& lo, const float& hi);
  void SetRapidity(const float& lo, const float& hi);
  void SetDCA(const float& lo, const float& hi);
  void SetCharge(const int&);


private:   // here are the quantities I want to cut on...

  int               fCharge;
  float             fNSigmaPion[2];
  float             fNSigmaKaon[2];
  float             fNSigmaProton[2];
  int               fNHits[2];
  float             fPt[2];
  float             fRapidity[2];
  float             fDCA[2];

  long              fNTracksPassed;
  long              fNTracksFailed;

#ifdef __ROOT__ 
  ClassDef(AliFemtoBasicTrackCut, 1)
#endif
};


inline void AliFemtoBasicTrackCut::SetNSigmaPion(const float& lo, const float& hi){fNSigmaPion[0]=lo; fNSigmaPion[1]=hi;}
inline void AliFemtoBasicTrackCut::SetNSigmaKaon(const float& lo, const float& hi){fNSigmaKaon[0]=lo; fNSigmaKaon[1]=hi;}
inline void AliFemtoBasicTrackCut::SetNSigmaProton(const float& lo, const float& hi){fNSigmaProton[0]=lo; fNSigmaProton[1]=hi;}

inline void AliFemtoBasicTrackCut::SetNHits(const int& lo, const int& hi){fNHits[0]=lo;fNHits[1]=hi;}
inline void AliFemtoBasicTrackCut::SetPt(const float& lo, const float& hi){fPt[0]=lo; fPt[1]=hi;}
inline void AliFemtoBasicTrackCut::SetRapidity(const float& lo,const float& hi){fRapidity[0]=lo; fRapidity[1]=hi;}
inline void AliFemtoBasicTrackCut::SetDCA(const float& lo,const float& hi){fDCA[0]=lo; fDCA[1]=hi;}
inline void AliFemtoBasicTrackCut::SetCharge(const int& ch){fCharge = ch;}

#endif
