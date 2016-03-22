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
// Module: EvtSemiLeptonicTensorAmp.cc
//
// Description: Routine to implement semileptonic decays to pseudo-scalar
//              mesons. 
//
// Modification history:
//
//    DJL       April 17,1998       Module created
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
#include "EvtGenBase/EvtSemiLeptonicTensorAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"

void EvtSemiLeptonicTensorAmp::CalcAmp( EvtParticle *parent,
					EvtAmp& amp,
					EvtSemiLeptonicFF *FormFactors ) {
  static EvtId EM=EvtPDL::getId("e-");
  static EvtId MUM=EvtPDL::getId("mu-");
  static EvtId TAUM=EvtPDL::getId("tau-");
  static EvtId EP=EvtPDL::getId("e+");
  static EvtId MUP=EvtPDL::getId("mu+");
  static EvtId TAUP=EvtPDL::getId("tau+");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");
  static EvtId DSM=EvtPDL::getId("D_s-");
  static EvtId DSP=EvtPDL::getId("D_s+");

  //Add the lepton and neutrino 4 momenta to find q2

  EvtVector4R q = parent->getDaug(1)->getP4() 
                    + parent->getDaug(2)->getP4(); 
  double q2 = (q.mass2());

  double hf,kf,bpf,bmf;

  FormFactors->gettensorff(parent->getId(),
                           parent->getDaug(0)->getId(),
                           q2,
			   parent->getDaug(0)->mass(),
                           &hf, 
                           &kf, 
                           &bpf, 
                           &bmf);


  double costhl_flag = 1.0;

  if(parent->getId()==D0||parent->getId()==D0B||
     parent->getId()==DP||parent->getId()==DM) {
    costhl_flag = -1.0;
  }
  if(parent->getId()==DSP||parent->getId()==DSM) {
    costhl_flag = -1.0;
  }
  hf = hf * costhl_flag;

  EvtVector4R p4b;
  p4b.set(parent->mass(),0.0,0.0,0.0);
 
  EvtVector4R p4meson = parent->getDaug(0)->getP4();
 
  EvtVector4C l1,l2;

  EvtId l_num = parent->getDaug(1)->getId();

  EvtVector4C ep_meson_b[5];

  ep_meson_b[0] = ((parent->getDaug(0)->epsTensorParent(0)).cont2(p4b)).conj();
  ep_meson_b[1] = ((parent->getDaug(0)->epsTensorParent(1)).cont2(p4b)).conj();
  ep_meson_b[2] = ((parent->getDaug(0)->epsTensorParent(2)).cont2(p4b)).conj();
  ep_meson_b[3] = ((parent->getDaug(0)->epsTensorParent(3)).cont2(p4b)).conj();
  ep_meson_b[4] = ((parent->getDaug(0)->epsTensorParent(4)).cont2(p4b)).conj();

  EvtVector4R pp,pm;

  pp=p4b+p4meson;
  pm=p4b-p4meson;

  //lange - October 31,2002 - try to lessen the mass dependence of probmax
  double q2max = p4b.mass2() + p4meson.mass2() - 2.0*p4b.mass()*p4meson.mass();
  double q2maxin=1.0/q2max;

  EvtComplex ep_meson_bb[5];

  ep_meson_bb[0]=ep_meson_b[0]*(p4b);
  ep_meson_bb[1]=ep_meson_b[1]*(p4b);
  ep_meson_bb[2]=ep_meson_b[2]*(p4b);
  ep_meson_bb[3]=ep_meson_b[3]*(p4b);
  ep_meson_bb[4]=ep_meson_b[4]*(p4b);


  EvtVector4C tds0,tds1,tds2,tds3,tds4;

  EvtTensor4C tds;
  if (l_num==EM||l_num==MUM||l_num==TAUM){
    EvtTensor4C tdual=EvtComplex(0.0,hf)*dual(EvtGenFunctions::directProd(pp,pm));
    tds0=tdual.cont2(ep_meson_b[0])
      -kf*ep_meson_b[0]
      -bpf*ep_meson_bb[0]*pp-bmf*ep_meson_bb[0]*pm;
    tds0*=q2maxin;

    tds1=tdual.cont2(ep_meson_b[1])
      -kf*ep_meson_b[1]
      -bpf*ep_meson_bb[1]*pp-bmf*ep_meson_bb[1]*pm;
    tds1*=q2maxin;

    tds2=tdual.cont2(ep_meson_b[2])
      -kf*ep_meson_b[2]
      -bpf*ep_meson_bb[2]*pp-bmf*ep_meson_bb[2]*pm;
    tds2*=q2maxin;

    tds3=tdual.cont2(ep_meson_b[3])
      -kf*ep_meson_b[3]
      -bpf*ep_meson_bb[3]*pp-bmf*ep_meson_bb[3]*pm;
    tds3*=q2maxin;

    tds4=tdual.cont2(ep_meson_b[4])
      -kf*ep_meson_b[4]
      -bpf*ep_meson_bb[4]*pp-bmf*ep_meson_bb[4]*pm;
    tds4*=q2maxin;


    l1=EvtLeptonVACurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParentNeutrino());
    l2=EvtLeptonVACurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParentNeutrino());
  }
  else{
    if (l_num==EP||l_num==MUP||l_num==TAUP){
     EvtTensor4C tdual=EvtComplex(0.0,-hf)*dual(EvtGenFunctions::directProd(pp,pm));
      tds0=tdual.cont2(ep_meson_b[0])
        -kf*ep_meson_b[0]
        -bpf*ep_meson_bb[0]*pp-bmf*ep_meson_bb[0]*pm;
      tds0*=q2maxin;

      tds1=tdual.cont2(ep_meson_b[1])
        -kf*ep_meson_b[1]
        -bpf*ep_meson_bb[1]*pp-bmf*ep_meson_bb[1]*pm;
      tds1*=q2maxin;

      tds2=tdual.cont2(ep_meson_b[2])
        -kf*ep_meson_b[2]
        -bpf*ep_meson_bb[2]*pp-bmf*ep_meson_bb[2]*pm;
      tds2*=q2maxin;

      tds3=tdual.cont2(ep_meson_b[3])
        -kf*ep_meson_b[3]
        -bpf*ep_meson_bb[3]*pp-bmf*ep_meson_bb[3]*pm;
      tds3*=q2maxin;

      tds4=tdual.cont2(ep_meson_b[4])
        -kf*ep_meson_b[4]
        -bpf*ep_meson_bb[4]*pp-bmf*ep_meson_bb[4]*pm;
      tds4*=q2maxin;

      l1=EvtLeptonVACurrent(parent->getDaug(2)->spParentNeutrino(),
			    parent->getDaug(1)->spParent(0));
      l2=EvtLeptonVACurrent(parent->getDaug(2)->spParentNeutrino(),
			    parent->getDaug(1)->spParent(1));
    }
    else{
      report(Severity::Error,"EvtGen") << "dfnb89agngri wrong lepton number\n";
    }
  }
 
  amp.vertex(0,0,l1*tds0);
  amp.vertex(0,1,l2*tds0);

  amp.vertex(1,0,l1*tds1);
  amp.vertex(1,1,l2*tds1);

  amp.vertex(2,0,l1*tds2);
  amp.vertex(2,1,l2*tds2);

  amp.vertex(3,0,l1*tds3);
  amp.vertex(3,1,l2*tds3);

  amp.vertex(4,0,l1*tds4);
  amp.vertex(4,1,l2*tds4);

  return;
 
}

