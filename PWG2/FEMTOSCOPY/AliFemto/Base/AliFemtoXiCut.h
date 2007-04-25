////////////////////////////////////////////////////////////////////////////////
/// AliFemtoXiCut - the pure virtual base class for the Xi cut               ///
/// All Xi cuts must inherit from this one                                   ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoXiCut_hh
#define AliFemtoXiCut_hh

#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoXi.h"
#include "Base/AliFemtoParticleCut.h"

class AliFemtoXiCut : public AliFemtoParticleCut {

public:

  AliFemtoXiCut(){/* no-op */};                          // default constructor. - Users should write their own
  AliFemtoXiCut(const AliFemtoXiCut& aCut);              // copy constructor
  virtual ~AliFemtoXiCut(){/* no-op */};                 // destructor

  virtual bool Pass(const AliFemtoXi* aCut)=0;               // true if passes, false if not

  virtual AliFemtoParticleType Type(){return hbtXi;}
  virtual AliFemtoXiCut* Clone() { return 0;}

#ifdef __ROOT__
  ClassDef(AliFemtoXiCut, 0)
#endif
};

inline AliFemtoXiCut::AliFemtoXiCut(const AliFemtoXiCut& c) : AliFemtoParticleCut(c) { /* no-op */ } 

#endif
