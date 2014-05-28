// $Id: EvtSSD_DirectCP.hh,v 1.2 2009-03-16 16:31:53 robbep Exp $
// Generation of direct CP violation in hadronic environment
// Patrick Robbe, LHCb,  08 Nov 2006
// 

#ifndef EVTSSD_DirectCP_HH
#define EVTSSD_DirectCP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtSSD_DirectCP : public  EvtDecayAmp  {

public:

  EvtSSD_DirectCP() {}
  virtual ~EvtSSD_DirectCP();
  
  virtual std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();
  void decay(EvtParticle *p); 

  std::string getParamName(int i);

private:
  bool isB0Mixed( EvtParticle * p ) ;
  bool isBsMixed( EvtParticle * p ) ;

  //Arguments

  double _acp;
};

#endif
