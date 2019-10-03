///
/// \file AliFemtoTrackCut.h
///

#ifndef AliFemtoTrackCut_hh
#define AliFemtoTrackCut_hh

#pragma once

#include "AliFemtoTypes.h"
#include "AliFemtoTrack.h"
#include "AliFemtoParticleCut.h"


/// \class AliFemtoTrackCut
/// \brief The pure virtual base class for track cuts. All track cuts must
///        inherit from this one.
///
class AliFemtoTrackCut : public AliFemtoParticleCut {
public:

  AliFemtoTrackCut();                                   ///< default constructor - Users should write their own
  AliFemtoTrackCut(const AliFemtoTrackCut&);            ///< copy constructor
  virtual ~AliFemtoTrackCut();                          ///< destructor
  AliFemtoTrackCut& operator=(const AliFemtoTrackCut&); ///< Assignment operator

  virtual bool Pass(const AliFemtoTrack* track) = 0;    ///< Returns true if passes, false if not
  virtual AliFemtoParticleType Type();                  ///< Always returns hbtTrack
  virtual AliFemtoTrackCut* Clone();                    ///< Returns NULL - subclasses must overload this method to use

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoTrackCut, 0);
  /// \endcond
#endif
};


inline AliFemtoTrackCut::AliFemtoTrackCut():
  AliFemtoParticleCut()
{ // no-op
}

inline AliFemtoTrackCut::AliFemtoTrackCut(const AliFemtoTrackCut& c):
  AliFemtoParticleCut(c)
{ // no-op
}

inline AliFemtoTrackCut::~AliFemtoTrackCut()
{ // no-op
}

inline AliFemtoTrackCut& AliFemtoTrackCut::operator=(const AliFemtoTrackCut& c)
{
  if (this == &c) {
    return *this;
  }
  AliFemtoParticleCut::operator=(c);
  return *this;
}

inline AliFemtoParticleType AliFemtoTrackCut::Type()
{
  return hbtTrack;
}

inline AliFemtoTrackCut* AliFemtoTrackCut::Clone()
{
  return NULL;
}

#endif
