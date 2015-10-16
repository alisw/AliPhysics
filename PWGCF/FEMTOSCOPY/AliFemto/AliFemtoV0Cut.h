///
/// \file AliFemtoV0Cut.h
///

#ifndef AliFemtoV0Cut_hh
#define AliFemtoV0Cut_hh

#include "AliFemtoTypes.h"
#include "AliFemtoV0.h"
#include "AliFemtoXi.h"
#include "AliFemtoParticleCut.h"

/// \class AliFemtoV0Cut
/// \brief The pure virtual base class for the V0 cut
///
/// All V0 cuts must inherit from this one.
///
class AliFemtoV0Cut : public AliFemtoParticleCut {
public:

  AliFemtoV0Cut();                          ///< default constructor. - Users should write their own
  AliFemtoV0Cut(const AliFemtoV0Cut& aCut); ///< copy constructor
  virtual ~AliFemtoV0Cut();                 ///< destructor
  AliFemtoV0Cut& operator=(const AliFemtoV0Cut& aCut); ///< copy constructor

  virtual bool Pass(const AliFemtoV0* aV0) = 0; ///< true if V0 passes, false if not
  virtual bool Pass(const AliFemtoXi* aXi) = 0; ///< true if Xi passes, false if not

  virtual AliFemtoParticleType Type() { return hbtV0; };
  virtual AliFemtoV0Cut* Clone() { return NULL; }; ///< WARNING - default implementation returns NULL

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0Cut, 0);
  /// \endcond
#endif
};

inline AliFemtoV0Cut::AliFemtoV0Cut()
{ // no-op
}

inline AliFemtoV0Cut::AliFemtoV0Cut(const AliFemtoV0Cut& c):
  AliFemtoParticleCut(c)
{ // no-op
}

inline AliFemtoV0Cut::~AliFemtoV0Cut()
{ // no-op
}

inline AliFemtoV0Cut& AliFemtoV0Cut::operator=(const AliFemtoV0Cut& c)
{
  if (this != &c) {
    AliFemtoParticleCut::operator=(c);
  }
  return *this;
}

#endif
