////////////////////////////////////////////////////////////////////////////////
/// AliFemtoParticleCut - the pure virtual base class for the particle cut   ///
/// All particle cuts must inherit from this one                             ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoParticleCut_hh
#define AliFemtoParticleCut_hh

#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoCutMonitorHandler.h"

class AliFemtoBaseAnalysis;

class AliFemtoParticleCut : public AliFemtoCutMonitorHandler {

  friend class AliFemtoBaseAnalysis;

public:
  AliFemtoParticleCut(){/* no-op */};   // default constructor. - Users should write their own
  AliFemtoParticleCut(const AliFemtoParticleCut&); // copy constructor
  virtual ~AliFemtoParticleCut(){/* no-op */};  // destructor

  virtual AliFemtoString Report() =0;    // user-written method to return string describing cuts

  double Mass(){return fMass;};       // mass of the particle being selected
  virtual void SetMass(const double& mass) {fMass = mass;};

  virtual void EventBegin(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void EventEnd(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual AliFemtoParticleCut* Clone() { return 0;}

  virtual AliFemtoParticleType Type()=0;

  // the following allows "back-pointing" from the CorrFctn to the "parent" Analysis
  AliFemtoBaseAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoBaseAnalysis*);

protected:
  double fMass;
  AliFemtoBaseAnalysis* fyAnalysis; // Link to the base analysis class

#ifdef __ROOT__
  ClassDef(AliFemtoParticleCut, 0)
#endif
};

inline AliFemtoParticleCut::AliFemtoParticleCut(const AliFemtoParticleCut& c) : AliFemtoCutMonitorHandler() { 
  fMass = c.fMass; fyAnalysis = 0; 
#ifdef STHBTDEBUG
  cout << " AliFemtoParticleCut::AliFemtoParticleCut(const AliFemtoParticleCut& c) - fMass: " << fMass << endl;
#endif
}
inline void AliFemtoParticleCut::SetAnalysis(AliFemtoBaseAnalysis* analysis) { fyAnalysis = analysis; }
#endif
