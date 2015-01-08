#ifndef EvtBcVNpi_HH
#define EvtBcVNpi_HH

#include <iostream>
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenModels/EvtBCVFF.hh"

#include "EvtGenModels/EvtWnPi.hh"


using std::endl;
using std::fstream;
using std::ifstream;
using std::cout;
using std::string;

class EvtBcVNpi:public  EvtDecayAmp  {
public:
    EvtBcVNpi() {  };
    virtual ~EvtBcVNpi();
    std::string getName();
    EvtDecayBase* clone();
    void initProbMax();
    void init();
    void decay(EvtParticle *p); 
protected:
  int nCall;
  int whichfit, idVector;
  EvtBCVFF *ffmodel;
  EvtWnPi *wcurr;

  EvtComplex Fpi( EvtVector4R q1, EvtVector4R q2);
  EvtComplex BWa( EvtVector4R q);
  EvtComplex BWf( EvtVector4R q);
  EvtComplex BWr( EvtVector4R q);
  EvtVector4C JB(EvtVector4R q1, EvtVector4R q2, EvtVector4R q3, EvtVector4R q4, EvtVector4R q5); 

};
#endif

