//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtSVPCP.cc
//
// Description: Routine to decay scalar -> vectors+photon
//              including CP violation effects
//
// Modification history:
//
//    Maurizio pierini   Nov 11, 2003       Module created
//
//------------------------------------------------------------------------
// 

#ifndef EVTSVPCP_HH
#define EVTSVPCP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

//Class to handle decays of the form SCALAR ->VECTOR PHOTON
//where the helicity amplitudes must be specified.  The
//first and third arguements are the magnetudes of the H+
//and H- helicity amplitudes respectively.  The second and
//fourth arguements are the phases.
//Calls EvtSVPHel.

class EvtSVPCP:public  EvtDecayAmp  {

public:

  EvtSVPCP() {}
  virtual ~EvtSVPCP();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();
  void decay(EvtParticle *p); 
  static void SVPHel(EvtParticle *parent,EvtAmp& amp,EvtId n_v1,EvtId n_v2,
	      const EvtComplex& hp,const EvtComplex& hm);

  std::string getParamName(int i);
  std::string getParamDefault(int i);
};

#endif
