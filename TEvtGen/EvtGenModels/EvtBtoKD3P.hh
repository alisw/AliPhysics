//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtBtoKD3P.hh,v 1.1 2009-03-16 16:49:00 robbep Exp $
// 
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2003, Colorado State University
//
// Module creator:
//      Abi soffer, CSU, 2003
//-----------------------------------------------------------------------

// Decay model that does the decay B+ -> K+ D , D -> 3 psudoscalars.
//
// The B- daughters specified in the decay file should be K-, D0, D0,
// where the first D0 is produced via b->c decay and the second via b->u.
// In reality, only one D daughter exists, so the first two
// daughters must be defined to decay to the same final state using
// the EvtPto3P model, but with CP-conjugate amplitudes. 
//
// For a given point in the Pto3P Dalitz plot,
// the total amplitude is \propto [A1 + A2 r exp(i(phase))], where
//
// A1 & A2 are the amplitudes of the D0 and D0bar to decay into that
// Dalitz plot point, 
//
// r is the (positive) ratio between the A(B->B0bar K) and A(B->D0 K)
// B decay amplitudes,
// 
// phase is the total phase difference (weak phase + strong phase) between
// A(B->D0bar K) and A(B->B0 K).
//
// Note that this model knows nothing about your convention for the
// sign of the phase, so when specifying the decay of a B- you need to
// change the order of D0 and D0bar and change the total phase so that
// the sign of the weak phase flips with respect to the parameters of B+.
// 

#ifndef EVT_BTOKD3P
#define EVT_BTOKD3P

class EvtParticle;
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtComplex.hh"


class EvtBtoKD3P : public  EvtDecayAmp {

public:
  EvtBtoKD3P();
  EvtBtoKD3P(const EvtBtoKD3P & other);
  virtual ~EvtBtoKD3P();
  EvtDecayBase* clone();
  
  // Initialize model
  virtual void init();      
  virtual void initProbMax();
  virtual void decay(EvtParticle *p);
  
  // we really have two daughters, although three are listed in the .dec file:
  virtual int nRealDaughters() { return 2;}

  std::string getName();
  
protected:
  // parameters:
  double _r;
  EvtComplex _exp; 

  // other:
  const EvtDecayBase * _model1;
  const EvtDecayBase * _model2;
  bool _decayedOnce;

};


#endif





