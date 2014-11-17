//-----------------------------------------------------------------------------------------------
//
// Module: EvtHybNonLepton.hh
//
// Desription: Routine to implement Hyperon(s=1/2) -> Baryon(s=1/2) + Scalar decays accroding to
//             Review Of Particle Physics 2004, Phys.Lett.B, Vol.592, p.864
//
// Modification history:
//
//  20/02/2005  PR   Module created according to PHSP and Lb2Lll model
//
//-----------------------------------------------------------------------------------------------

#ifndef EVTHYBNONLEPTON_HH
#define EVTHYBNONLEPTON_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtComplex.hh"

class EvtHypNonLepton:public EvtDecayAmp {

public:

  EvtHypNonLepton() {}
  virtual ~EvtHypNonLepton();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void init();
  void initProbMax();

  void calcAmp(EvtAmp *amp,EvtParticle *parent);

private:

  double     m_alpha;
  double     m_phi;
  EvtComplex m_B_to_A;
  long       m_noTries;

};

#endif
