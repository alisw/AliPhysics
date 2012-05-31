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
// Module: EvtSemiLeptonicBaryonAmp.cc
//
// Description: Routine to implement semileptonic decays to vector
//              mesons. 
//
// Modification history:
//
//    Lange    Oct 20, 2004   Module created.
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtSemiLeptonicBaryonAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
using std::endl;

void EvtSemiLeptonicBaryonAmp::CalcAmp( EvtParticle *parent,
					EvtAmp& amp,
					EvtSemiLeptonicFF *FormFactors ) {

  static EvtId EM=EvtPDL::getId("e-");
  static EvtId MUM=EvtPDL::getId("mu-");
  static EvtId TAUM=EvtPDL::getId("tau-");
  static EvtId EP=EvtPDL::getId("e+");
  static EvtId MUP=EvtPDL::getId("mu+");
  static EvtId TAUP=EvtPDL::getId("tau+");

 
  //Add the lepton and neutrino 4 momenta to find q2

  EvtVector4R q = parent->getDaug(1)->getP4() 
                    + parent->getDaug(2)->getP4(); 
  double q2 = (q.mass2());

  double f1v,f1a,f2v,f2a;
  double m_meson = parent->getDaug(0)->mass();

  FormFactors->getbaryonff(parent->getId(),
			   parent->getDaug(0)->getId(),
                           q2,
                           m_meson,
                           &f1v, 
                           &f1a, 
                           &f2v, 
                           &f2a);

  EvtVector4R p4b;
  p4b.set(parent->mass(),0.0,0.0,0.0);
  
  EvtVector4C temp_00_term1;
  EvtVector4C temp_00_term2;
  
  EvtVector4C temp_01_term1;
  EvtVector4C temp_01_term2;
  
  EvtVector4C temp_10_term1;
  EvtVector4C temp_10_term2;
  
  EvtVector4C temp_11_term1;
  EvtVector4C temp_11_term2;
  
  EvtDiracSpinor p0=parent->sp(0);
  EvtDiracSpinor p1=parent->sp(1);
  
  EvtDiracSpinor d0=parent->getDaug(0)->spParent(0);
  EvtDiracSpinor d1=parent->getDaug(0)->spParent(1);
  
  temp_00_term1.set(0,f1v*(d0*(EvtGammaMatrix::g0()*p0)));
  temp_00_term2.set(0,f1a*(d0*((EvtGammaMatrix::g0()*EvtGammaMatrix::g5())*p0)));
  temp_01_term1.set(0,f1v*(d0*(EvtGammaMatrix::g0()*p1)));
  temp_01_term2.set(0,f1a*(d0*((EvtGammaMatrix::g0()*EvtGammaMatrix::g5())*p1)));
  temp_10_term1.set(0,f1v*(d1*(EvtGammaMatrix::g0()*p0)));
  temp_10_term2.set(0,f1a*(d1*((EvtGammaMatrix::g0()*EvtGammaMatrix::g5())*p0)));
  temp_11_term1.set(0,f1v*(d1*(EvtGammaMatrix::g0()*p1)));
  temp_11_term2.set(0,f1a*(d1*((EvtGammaMatrix::g0()*EvtGammaMatrix::g5())*p1)));
  
  temp_00_term1.set(1,f1v*(d0*(EvtGammaMatrix::g1()*p0)));
  temp_00_term2.set(1,f1a*(d0*((EvtGammaMatrix::g1()*EvtGammaMatrix::g5())*p0)));
  temp_01_term1.set(1,f1v*(d0*(EvtGammaMatrix::g1()*p1)));
  temp_01_term2.set(1,f1a*(d0*((EvtGammaMatrix::g1()*EvtGammaMatrix::g5())*p1)));
  temp_10_term1.set(1,f1v*(d1*(EvtGammaMatrix::g1()*p0)));
  temp_10_term2.set(1,f1a*(d1*((EvtGammaMatrix::g1()*EvtGammaMatrix::g5())*p0)));
  temp_11_term1.set(1,f1v*(d1*(EvtGammaMatrix::g1()*p1)));
  temp_11_term2.set(1,f1a*(d1*((EvtGammaMatrix::g1()*EvtGammaMatrix::g5())*p1)));
  
  temp_00_term1.set(2,f1v*(d0*(EvtGammaMatrix::g2()*p0)));
  temp_00_term2.set(2,f1a*(d0*((EvtGammaMatrix::g2()*EvtGammaMatrix::g5())*p0)));
  temp_01_term1.set(2,f1v*(d0*(EvtGammaMatrix::g2()*p1)));
  temp_01_term2.set(2,f1a*(d0*((EvtGammaMatrix::g2()*EvtGammaMatrix::g5())*p1)));
  temp_10_term1.set(2,f1v*(d1*(EvtGammaMatrix::g2()*p0)));
  temp_10_term2.set(2,f1a*(d1*((EvtGammaMatrix::g2()*EvtGammaMatrix::g5())*p0)));
  temp_11_term1.set(2,f1v*(d1*(EvtGammaMatrix::g2()*p1)));
  temp_11_term2.set(2,f1a*(d1*((EvtGammaMatrix::g2()*EvtGammaMatrix::g5())*p1)));
  
  temp_00_term1.set(3,f1v*(d0*(EvtGammaMatrix::g3()*p0)));
  temp_00_term2.set(3,f1a*(d0*((EvtGammaMatrix::g3()*EvtGammaMatrix::g5())*p0)));
  temp_01_term1.set(3,f1v*(d0*(EvtGammaMatrix::g3()*p1)));
  temp_01_term2.set(3,f1a*(d0*((EvtGammaMatrix::g3()*EvtGammaMatrix::g5())*p1)));
  temp_10_term1.set(3,f1v*(d1*(EvtGammaMatrix::g3()*p0)));
  temp_10_term2.set(3,f1a*(d1*((EvtGammaMatrix::g3()*EvtGammaMatrix::g5())*p0)));
  temp_11_term1.set(3,f1v*(d1*(EvtGammaMatrix::g3()*p1)));
  temp_11_term2.set(3,f1a*(d1*((EvtGammaMatrix::g3()*EvtGammaMatrix::g5())*p1)));
  


  EvtVector4C l1,l2;

  EvtId l_num = parent->getDaug(1)->getId();
  if (l_num==EM||l_num==MUM||l_num==TAUM){

    l1=EvtLeptonVACurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParentNeutrino());
    l2=EvtLeptonVACurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParentNeutrino());
  }
  else{
    if (l_num==EP||l_num==MUP||l_num==TAUP){
    l1=EvtLeptonVACurrent(parent->getDaug(2)->spParentNeutrino(),
			    parent->getDaug(1)->spParent(0));
    l2=EvtLeptonVACurrent(parent->getDaug(2)->spParentNeutrino(),
			    parent->getDaug(1)->spParent(1));
    }
    else{
      report(ERROR,"EvtGen") << "Wrong lepton number"<<endl;
    }
  }

  amp.vertex(0,0,0,l1.cont(temp_00_term1+temp_00_term2));
  amp.vertex(0,0,1,l2.cont(temp_00_term1+temp_00_term2));

  amp.vertex(0,1,0,l1.cont(temp_01_term1+temp_01_term2));
  amp.vertex(0,1,1,l2.cont(temp_01_term1+temp_01_term2));

  amp.vertex(1,0,0,l1.cont(temp_10_term1+temp_10_term2));
  amp.vertex(1,0,1,l2.cont(temp_10_term1+temp_10_term2));

  amp.vertex(1,1,0,l1.cont(temp_11_term1+temp_11_term2));
  amp.vertex(1,1,1,l2.cont(temp_11_term1+temp_11_term2));

  return;
}

