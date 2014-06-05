//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtPto3P.cpp,v 1.4 2009-03-16 15:46:31 robbep Exp $
// 
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998 Caltech, UCSB
//
// Module creator:
//      Alexei Dvoretskii, Caltech, 2001-2002.
//-----------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtPto3P.hh"
#include "EvtGenBase/EvtPto3PAmpFactory.hh"
using namespace EvtCyclic3;

EvtDalitzPlot EvtPto3P::dp()
{
  // There must be 3 daughters. All particles must be pseudoscalars.
  // Charge must be conserved. Number of arguments must be non-zero.

  EvtId parent = getParentId();
  assert(getNDaug() == 3);
  EvtId dau0 = getDaug(0);
  EvtId dau1 = getDaug(1);
  EvtId dau2 = getDaug(2);
  
  assert(EvtPDL::getSpinType(parent) == EvtSpinType::SCALAR);
  assert(EvtPDL::getSpinType(dau0) == EvtSpinType::SCALAR);
  assert(EvtPDL::getSpinType(dau1) == EvtSpinType::SCALAR);
  assert(EvtPDL::getSpinType(dau2) == EvtSpinType::SCALAR);  
  assert(EvtPDL::chg3(parent) == EvtPDL::chg3(dau0) + EvtPDL::chg3(dau1) + EvtPDL::chg3(dau2));
  assert(getNArg() > 0);

  return EvtDalitzPlot(EvtPDL::getMass(dau0),EvtPDL::getMass(dau1),EvtPDL::getMass(dau2),EvtPDL::getMass(parent));
}

EvtAmpFactory<EvtDalitzPoint>* EvtPto3P::createFactory(const EvtMultiChannelParser& parser)
{
  // Compute the interval size

  EvtDalitzPlot plot = dp();
  EvtAmpFactory<EvtDalitzPoint>* fact = new EvtPto3PAmpFactory(plot);
  fact->build(parser,10000);
  return fact;
}


std::vector<EvtVector4R> EvtPto3P::initDaughters(const EvtDalitzPoint& x) const
{
  std::vector<EvtVector4R> v;
  assert(x.isValid());
  
  // Calculate in the r.f. of AB
                                                                              
  double eA = x.e(A,AB);
  double eB = x.e(B,AB);
  double eC = x.e(C,AB);
  double pA = x.p(A,AB);
  double pC = x.p(C,AB);
  double cos = x.cosTh(CA,AB);
  double sin = sqrt(1.0-cos*cos);
                                                                                  
  EvtVector4R vA(eA,0,0,pA);                                                   
  EvtVector4R vB(eB,0,0,-pA); 
  EvtVector4R vC(eC,0,pC*sin,pC*cos);

  // Boost from rest frame of AB to rest-frame of decaying particle
  // vboost is the 4-momentum of frame being boosted from in the frame
  // being boosted into.
 
  EvtVector4R vboost = vA + vB + vC;
  vboost.set(1,-vboost.get(1));
  vboost.set(2,-vboost.get(2));
  vboost.set(3,-vboost.get(3));
  vA.applyBoostTo(vboost);
  vB.applyBoostTo(vboost);
  vC.applyBoostTo(vboost);
  
  // Rotate
                                                                                
  double alpha = EvtRandom::Flat( EvtConst::twoPi );
  double beta = acos(EvtRandom::Flat( -1.0, 1.0 ));
  double gamma = EvtRandom::Flat( EvtConst::twoPi );
                                                                                
  vA.applyRotateEuler( alpha, beta, gamma );
  vB.applyRotateEuler( alpha, beta, gamma );
  vC.applyRotateEuler( alpha, beta, gamma );
                                                                                 
   // Fill vector
                                                                                 
  assert(v.size() == 0);  
  v.push_back(vA);
  v.push_back(vB);
  v.push_back(vC);

  return v;
}                                                                             











