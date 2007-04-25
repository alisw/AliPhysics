////////////////////////////////////////////////////////////////////////////////
/// AliFemtoV0Cut - the pure virtual base class for the V0 cut           ///
/// All V0 cuts must inherit from this one                                 ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoV0Cut_hh
#define AliFemtoV0Cut_hh

#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoV0.h"
#include "Base/AliFemtoParticleCut.h"

class AliFemtoV0Cut : public AliFemtoParticleCut {

public:

  AliFemtoV0Cut(){/* no-op */};                             // default constructor. - Users should write their own
  AliFemtoV0Cut(const AliFemtoV0Cut& aCut);                 // copy constructor
  virtual ~AliFemtoV0Cut(){/* no-op */};                    // destructor

  virtual bool Pass(const AliFemtoV0* aV0)=0;               // true if passes, false if not

  virtual AliFemtoParticleType Type(){return hbtV0;}
  virtual AliFemtoV0Cut* Clone() { return 0;}

#ifdef __ROOT__
  ClassDef(AliFemtoV0Cut, 0)
#endif
};

inline AliFemtoV0Cut::AliFemtoV0Cut(const AliFemtoV0Cut& c) : AliFemtoParticleCut(c) { /* no-op */ } 

#endif
