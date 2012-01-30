////////////////////////////////////////////////////////////////////////////////
/// AliFemtoTrackCut - the pure virtual base class for the track cut         ///
/// All track cuts must inherit from this one                                ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoTrackCut_hh
#define AliFemtoTrackCut_hh

#include "AliFemtoTypes.h"
#include "AliFemtoTrack.h"
#include "AliFemtoParticleCut.h"

class AliFemtoTrackCut : public AliFemtoParticleCut {

public:

  AliFemtoTrackCut(){/* no-op */};                       // default constructor. - Users should write their own
  AliFemtoTrackCut(const AliFemtoTrackCut&);                // copy constructor
  virtual ~AliFemtoTrackCut(){/* no-op */};              // destructor
  AliFemtoTrackCut& operator=(const AliFemtoTrackCut&);                // copy constructor

  virtual bool Pass(const AliFemtoTrack* track)=0;       // true if passes, false if not
  virtual AliFemtoParticleType Type(){return hbtTrack;}
  virtual AliFemtoTrackCut* Clone() { return 0;}

#ifdef __ROOT__
  ClassDef(AliFemtoTrackCut, 0)
#endif
};

inline AliFemtoTrackCut::AliFemtoTrackCut(const AliFemtoTrackCut& c) : AliFemtoParticleCut(c) {
#ifdef STHBTDEBUG
  cout << " AliFemtoTrackCut::AliFemtoTrackCut(const AliFemtoTrackCut& c) : AliFemtoParticleCut(c) " << endl;
#endif
}
inline AliFemtoTrackCut& AliFemtoTrackCut::operator=(const AliFemtoTrackCut& c) { if (this != &c) { AliFemtoParticleCut::operator=(c); } return *this; } 
#endif
