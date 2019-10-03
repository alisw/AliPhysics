////////////////////////////////////////////////////////////////////////////////
/// AliFemtoKinkCut - the pure virtual base class for the kink cut           ///
/// All kink cuts must inherit from this one                                 ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoKinkCut_hh
#define AliFemtoKinkCut_hh

#include "AliFemtoTypes.h"
#include "AliFemtoKink.h"
#include "AliFemtoParticleCut.h"

class AliFemtoKinkCut : public AliFemtoParticleCut {

public:

  AliFemtoKinkCut(){/* no-op */};                       // default constructor. - Users should write their own
  AliFemtoKinkCut(const AliFemtoKinkCut&);                         // copy constructor
  virtual ~AliFemtoKinkCut(){/* no-op */};              // destructor
  AliFemtoKinkCut& operator=(const AliFemtoKinkCut&);                         // copy constructor

  virtual bool Pass(const AliFemtoKink* aKink)=0;               // true if passes, false if not

  virtual AliFemtoParticleType Type(){return hbtKink;}
  virtual AliFemtoKinkCut* Clone() { return 0;}

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoKinkCut, 0);
  /// \endcond
#endif
};
//_____________________________
inline AliFemtoKinkCut::AliFemtoKinkCut(const AliFemtoKinkCut& c) : AliFemtoParticleCut(c) { /* no-op */ } 
inline AliFemtoKinkCut& AliFemtoKinkCut::operator=(const AliFemtoKinkCut& aCorrFctn) {   if (this != &aCorrFctn) { AliFemtoParticleCut::operator=(aCorrFctn); } return *this; }

#endif
