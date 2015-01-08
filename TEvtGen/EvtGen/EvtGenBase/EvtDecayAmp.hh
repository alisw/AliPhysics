//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtDecayAmp.hh
//
// Description: Baseclass for models that calculates amplitudes
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTDECAYAMP_HH
#define EVTDECAYAMP_HH

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtAmp.hh"

class EvtDecayAmp : public EvtDecayBase{

public:

  void makeDecay(EvtParticle* p, bool recursive=true);
  inline void setWeight(double weight) {_weight=weight;}

  /**
  * sets the amplitudes calculated in the decay objects
  */
  void vertex(const EvtComplex& amp){_amp2.vertex(amp);}

  /**
  * sets the amplitudes calculated in the decay objects
  */
  void vertex(int i1, const EvtComplex& amp){_amp2.vertex(i1,amp);}

  /**
  * sets the amplitudes calculated in the decay objects
  */
  void vertex(int i1, int i2, const EvtComplex& amp)
  {_amp2.vertex(i1,i2,amp);}


  /**
  * sets the amplitudes calculated in the decay objects
  */
  void vertex(int i1, int i2, int i3, const EvtComplex& amp)
  {_amp2.vertex(i1,i2,i3,amp);}

  /**
  * sets the amplitudes calculated in the decay objects
  */
  void vertex(int *i1, const EvtComplex& amp)
  { _amp2.vertex(i1,amp);}

  /**
   *  Provide access to the amplitude
   */
  const EvtAmp & amplitude() const 
  { return _amp2;}




  virtual ~EvtDecayAmp() {}

protected:
  EvtAmp _amp2;

private:
  double _weight;


};



#endif
