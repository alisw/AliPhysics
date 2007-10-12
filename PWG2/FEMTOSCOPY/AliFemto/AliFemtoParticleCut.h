////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoParticleCut - the pure virtual base class for the particle cut     //
// All particle cuts must inherit from this one                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

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
  AliFemtoParticleCut();   // default constructor. - Users should write their own
  AliFemtoParticleCut(const AliFemtoParticleCut&); // copy constructor
  virtual ~AliFemtoParticleCut(){/* no-op */};  // destructor
  AliFemtoParticleCut& operator=(const AliFemtoParticleCut& aCut);

  virtual AliFemtoString Report() =0;    // user-written method to return string describing cuts
  virtual TList* ListSettings();      // user-written list of settings which is stored in the result file

  double Mass(){return fMass;};       // mass of the particle being selected
  virtual void SetMass(const double& mass) {fMass = mass;};

  virtual void EventBegin(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual void EventEnd(const AliFemtoEvent* aEvent) { /* no-op */ }
  virtual AliFemtoParticleCut* Clone() { return 0;}

  virtual AliFemtoParticleType Type()=0;

  // the following allows "back-pointing" from the CorrFctn to the "parent" Analysis
  AliFemtoAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoAnalysis*);

protected:
  AliFemtoAnalysis* fyAnalysis; // Link to the base analysis class
  double fMass;

#ifdef __ROOT__
  ClassDef(AliFemtoParticleCut, 0)
#endif
};

inline AliFemtoParticleCut::AliFemtoParticleCut(): AliFemtoCutMonitorHandler(), fyAnalysis(0), fMass(0){}   // default constructor. - Users should write their own
inline AliFemtoParticleCut::AliFemtoParticleCut(const AliFemtoParticleCut& c): AliFemtoCutMonitorHandler(), fyAnalysis(0), fMass(0) { 
  fMass = c.fMass; fyAnalysis = 0; 
#ifdef STHBTDEBUG
  cout << " AliFemtoParticleCut::AliFemtoParticleCut(const AliFemtoParticleCut& c) - fMass: " << fMass << endl;
#endif
}
inline void AliFemtoParticleCut::SetAnalysis(AliFemtoAnalysis* analysis) { fyAnalysis = analysis; }
inline AliFemtoParticleCut& AliFemtoParticleCut::operator=(const AliFemtoParticleCut& aCut) { if (this == &aCut) return *this; fyAnalysis = aCut.fyAnalysis; fMass=aCut.fMass; return *this; }
  inline TList *AliFemtoParticleCut::ListSettings() { 
    TList *tListSetttings = new TList(); 
    char buf[100];
    snprintf(buf, 100, "AliFemtoParticleCut.mass=%lf", fMass);
    TObjString *str = new TObjString(buf);
    tListSetttings->Add(str);
    return tListSetttings;
  }
  
#endif
