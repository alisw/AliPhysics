///
/// \file AliFemtoParticleCut.h
/// \class AliFemtoParticleCut
/// \brief The pure virtual base class for the particle cut.
///
/// All particle cuts must inherit from this one
///


#ifndef ALIFEMTOPARTICLECUT_H
#define ALIFEMTOPARTICLECUT_H

#include "AliFemtoTypes.h"
#include "AliFemtoCutMonitorHandler.h"
#include <TObjString.h>
#include <TList.h>

class AliFemtoAnalysis;

class AliFemtoParticleCut : public AliFemtoCutMonitorHandler {

  friend class AliFemtoAnalysis;

public:
  AliFemtoParticleCut();                            ///< Default constructor. - Users should write their own
  AliFemtoParticleCut(const AliFemtoParticleCut &); ///< Copy constructor
  virtual ~AliFemtoParticleCut() {/* no-op */ };    ///< Destructor
  AliFemtoParticleCut &operator=(const AliFemtoParticleCut &aCut);

  virtual AliFemtoString Report() = 0;    ///< User-written method to return string describing cuts
  virtual TList *ListSettings();          ///< User-written list of settings which is stored in the result file

  double Mass() { return fMass; };        ///< Mass of the particle being selected
  virtual void SetMass(const double &mass) { fMass = mass; };

  virtual AliFemtoParticleCut* Clone() { return NULL; }

  virtual AliFemtoParticleType Type() = 0;  ///< Pure virtual function which returns the particle type

  ///< The following allows "back-pointing" from the CorrFctn to the "parent" Analysis
  AliFemtoAnalysis* HbtAnalysis() { return fyAnalysis; };
  void SetAnalysis(AliFemtoAnalysis *anAnalysis) { fyAnalysis = anAnalysis; };

protected:
  AliFemtoAnalysis *fyAnalysis; ///< Link to the base analysis class
  double fMass;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoParticleCut, 0);
  /// \endcond
#endif
};

inline AliFemtoParticleCut::AliFemtoParticleCut():
  AliFemtoCutMonitorHandler(),
  fyAnalysis(NULL),
  fMass(0.0)
{
    /* no-op */
}

inline AliFemtoParticleCut::AliFemtoParticleCut(const AliFemtoParticleCut &c):
  AliFemtoCutMonitorHandler(c),
  fyAnalysis(c.fyAnalysis),
  fMass(c.fMass)
{
#ifdef STHBTDEBUG
  cout << " AliFemtoParticleCut::AliFemtoParticleCut(const AliFemtoParticleCut& c) - fMass: " << fMass << endl;
#endif
}

inline AliFemtoParticleCut &AliFemtoParticleCut::operator=(const AliFemtoParticleCut &aCut)
{
  if (this == &aCut) return *this;
  fyAnalysis = aCut.fyAnalysis;
  fMass = aCut.fMass;
  return *this;
}

inline TList *AliFemtoParticleCut::ListSettings()
{
  TList *tListSetttings = new TList();
  char buf[100];
  snprintf(buf, 100, "AliFemtoParticleCut.mass=%f", fMass);
  TObjString *str = new TObjString(buf);
  tListSetttings->Add(str);
  return tListSetttings;
}

#endif
