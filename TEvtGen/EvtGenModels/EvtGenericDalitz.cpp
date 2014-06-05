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
// Module: EvtGen/EvtGenericDalitz.hh
//
// Description: Model to describe a generic dalitz decay
//
// Modification history:
//
//    DCC     16 December, 2011         Module created
//
//------------------------------------------------------------------------

#include "EvtGenModels/EvtGenericDalitz.hh"
#include "EvtGenModels/EvtDalitzTable.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtPDL.hh"

std::string EvtGenericDalitz::getName() {
  return "GENERIC_DALITZ";
}

EvtDecayBase* EvtGenericDalitz::clone() {
  return new EvtGenericDalitz();
}

void EvtGenericDalitz::init() {

  checkNArg(1);

  EvtId parnum=getParentId();
  EvtId d1=getDaug(0);
  EvtId d2=getDaug(1);
  EvtId d3=getDaug(2);

  std::vector<EvtDalitzDecayInfo> decays = EvtDalitzTable::getInstance(getArgStr(0))->getDalitzTable(parnum);

  std::vector<EvtDalitzDecayInfo>::iterator i = decays.begin();
  for( ; i != decays.end(); i++) {
    EvtId daughter1 = (*i).daughter1();
    EvtId daughter2 = (*i).daughter2();
    EvtId daughter3 = (*i).daughter3();

    if(d1 == daughter1 && d2 == daughter2 && d3 == daughter3) {
      _d1 = 0;
      _d2 = 1;
      _d3 = 2;
    } else if(d1 == daughter1 && d2 == daughter3 && d3 == daughter2) {
      _d1 = 0;
      _d2 = 2;
      _d3 = 1;
    } else if(d1 == daughter2 && d2 == daughter1 && d3 == daughter3) {
      _d1 = 1;
      _d2 = 0;
      _d3 = 2;
    } else if(d1 == daughter2 && d2 == daughter3 && d3 == daughter1) {
      _d1 = 1;
      _d2 = 2;
      _d3 = 0;
    } else if(d1 == daughter3 && d2 == daughter1 && d3 == daughter2) {
      _d1 = 2;
      _d2 = 0;
      _d3 = 1;
    } else if(d1 == daughter3 && d2 == daughter2 && d3 == daughter1) {
      _d1 = 2;
      _d2 = 1;
      _d3 = 0;
    } else {
      continue;
    }

    _resonances = (*i).getResonances();
    setProbMax((*i).getProbMax());
    return;
  }
}

void EvtGenericDalitz::decay(EvtParticle *p) {

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtVector4R p4_d1 = p->getDaug(_d1)->getP4();
  EvtVector4R p4_d2 = p->getDaug(_d2)->getP4();
  EvtVector4R p4_d3 = p->getDaug(_d3)->getP4();

  double mA = p->getDaug(_d1)->mass();
  double mB = p->getDaug(_d2)->mass();
  double mC = p->getDaug(_d3)->mass();

  double m2AB = ( p4_d1 + p4_d2 ).mass2();
  double m2CA = ( p4_d1 + p4_d3 ).mass2();
  double m2BC = ( p4_d2 + p4_d3 ).mass2();

  EvtDalitzPoint point( mA, mB, mC, m2AB, m2BC, m2CA );

  EvtComplex amp(0,0);
  std::vector<std::pair<EvtComplex,EvtDalitzReso> >::iterator i = _resonances.begin();
  for( ; i!= _resonances.end(); i++) {
    std::pair<EvtComplex,EvtDalitzReso> res = (*i);
    amp += res.first * res.second.evaluate( point );
  }

  vertex(amp);
  return;
}

std::string EvtGenericDalitz::getParamName(int i) {
  switch(i) {
  case 0:
    return "xmlFile";
  default:
    return "";
  }
}
