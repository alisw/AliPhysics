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
// Module: EvtISGW2.cc
//
// Description: Routine to implement semileptonic decays according
//              to the model ISGW2
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtISGW2.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include <string>
#include "EvtGenModels/EvtISGW2FF.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicTensorAmp.hh"

EvtISGW2::EvtISGW2():
   isgw2ffmodel(0)
  ,calcamp(0)
{}


EvtISGW2::~EvtISGW2() {
  delete isgw2ffmodel;
  isgw2ffmodel=0;
  delete calcamp;
  calcamp=0;
}

std::string EvtISGW2::getName(){

  return "ISGW2";     

}



EvtDecayBase* EvtISGW2::clone(){

  return new EvtISGW2;

}

void EvtISGW2::decay( EvtParticle *p ){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  calcamp->CalcAmp(p,_amp2,isgw2ffmodel);

}

void EvtISGW2::initProbMax() {

  //added by Lange Jan4,2000
  static EvtId EM=EvtPDL::getId("e-");
  static EvtId EP=EvtPDL::getId("e+");
  static EvtId MUM=EvtPDL::getId("mu-");
  static EvtId MUP=EvtPDL::getId("mu+");
  static EvtId TAUM=EvtPDL::getId("tau-");
  static EvtId TAUP=EvtPDL::getId("tau+");

  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");
  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BCP=EvtPDL::getId("B_c+");
  static EvtId BCM=EvtPDL::getId("B_c-");

  static EvtId DST0=EvtPDL::getId("D*0");
  static EvtId DSTB=EvtPDL::getId("anti-D*0");
  static EvtId DSTP=EvtPDL::getId("D*+");
  static EvtId DSTM=EvtPDL::getId("D*-");
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId D1P1P=EvtPDL::getId("D_1+");
  static EvtId D1P1N=EvtPDL::getId("D_1-");
  static EvtId D1P10=EvtPDL::getId("D_10");
  static EvtId D1P1B=EvtPDL::getId("anti-D_10");

  static EvtId D3P2P=EvtPDL::getId("D_2*+");
  static EvtId D3P2N=EvtPDL::getId("D_2*-");
  static EvtId D3P20=EvtPDL::getId("D_2*0");
  static EvtId D3P2B=EvtPDL::getId("anti-D_2*0");

  static EvtId D3P1P=EvtPDL::getId("D'_1+");
  static EvtId D3P1N=EvtPDL::getId("D'_1-");
  static EvtId D3P10=EvtPDL::getId("D'_10");
  static EvtId D3P1B=EvtPDL::getId("anti-D'_10");

  static EvtId D3P0P=EvtPDL::getId("D_0*+");
  static EvtId D3P0N=EvtPDL::getId("D_0*-");
  static EvtId D3P00=EvtPDL::getId("D_0*0");
  static EvtId D3P0B=EvtPDL::getId("anti-D_0*0");

  static EvtId D21S0P=EvtPDL::getId("D(2S)+");
  static EvtId D21S0N=EvtPDL::getId("D(2S)-");
  static EvtId D21S00=EvtPDL::getId("D(2S)0");
  static EvtId D21S0B=EvtPDL::getId("anti-D(2S)0");

  static EvtId D23S1P=EvtPDL::getId("D*(2S)+");
  static EvtId D23S1N=EvtPDL::getId("D*(2S)-");
  static EvtId D23S10=EvtPDL::getId("D*(2S)0");
  static EvtId D23S1B=EvtPDL::getId("anti-D*(2S)0");

  static EvtId RHO2S0=EvtPDL::getId("rho(2S)0");
  static EvtId RHO2SP=EvtPDL::getId("rho(2S)+");
  static EvtId RHO2SM=EvtPDL::getId("rho(2S)-");
  static EvtId OMEG2S=EvtPDL::getId("omega(2S)");
  static EvtId ETA2S=EvtPDL::getId("eta(2S)");

  static EvtId PI2S0=EvtPDL::getId("pi(2S)0");
  static EvtId PI2SP=EvtPDL::getId("pi(2S)+");
  static EvtId PI2SM=EvtPDL::getId("pi(2S)-");

  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PIM=EvtPDL::getId("pi-");
  static EvtId PI0=EvtPDL::getId("pi0");

  static EvtId RHOP=EvtPDL::getId("rho+");
  static EvtId RHOM=EvtPDL::getId("rho-");
  static EvtId RHO0=EvtPDL::getId("rho0");

  static EvtId A2P=EvtPDL::getId("a_2+");
  static EvtId A2M=EvtPDL::getId("a_2-");
  static EvtId A20=EvtPDL::getId("a_20");

  static EvtId A1P=EvtPDL::getId("a_1+");
  static EvtId A1M=EvtPDL::getId("a_1-");
  static EvtId A10=EvtPDL::getId("a_10");

  static EvtId A0P=EvtPDL::getId("a_0+");
  static EvtId A0M=EvtPDL::getId("a_0-");
  static EvtId A00=EvtPDL::getId("a_00");

  static EvtId B1P=EvtPDL::getId("b_1+");
  static EvtId B1M=EvtPDL::getId("b_1-");
  static EvtId B10=EvtPDL::getId("b_10");

  static EvtId H1=EvtPDL::getId("h_1");
  static EvtId H1PR=EvtPDL::getId("h'_1");

  static EvtId F1=EvtPDL::getId("f_1");
  static EvtId F1PR=EvtPDL::getId("f'_1");
  static EvtId F0=EvtPDL::getId("f_0");
  static EvtId F0PR=EvtPDL::getId("f'_0");
  static EvtId F2=EvtPDL::getId("f_2");
  static EvtId F2PR=EvtPDL::getId("f'_2");

  static EvtId ETA=EvtPDL::getId("eta");
  static EvtId ETAPR=EvtPDL::getId("eta'");
  static EvtId OMEG=EvtPDL::getId("omega");

  static EvtId KP=EvtPDL::getId("K+");
  static EvtId KM=EvtPDL::getId("K-");
  static EvtId K0=EvtPDL::getId("K0");
  static EvtId KB=EvtPDL::getId("anti-K0");
  static EvtId K0S=EvtPDL::getId("K_S0");
  static EvtId K0L=EvtPDL::getId("K_L0");

  static EvtId KSTP=EvtPDL::getId("K*+");
  static EvtId KSTM=EvtPDL::getId("K*-");
  static EvtId KST0=EvtPDL::getId("K*0");
  static EvtId KSTB=EvtPDL::getId("anti-K*0");

  static EvtId K1P=EvtPDL::getId("K_1+");
  static EvtId K1M=EvtPDL::getId("K_1-");
  static EvtId K10=EvtPDL::getId("K_10");
  static EvtId K1B=EvtPDL::getId("anti-K_10");

  static EvtId K1STP=EvtPDL::getId("K'_1+");
  static EvtId K1STM=EvtPDL::getId("K'_1-");
  static EvtId K1ST0=EvtPDL::getId("K'_10");
  static EvtId K1STB=EvtPDL::getId("anti-K'_10");

  static EvtId K2STP=EvtPDL::getId("K_2*+");
  static EvtId K2STM=EvtPDL::getId("K_2*-");
  static EvtId K2ST0=EvtPDL::getId("K_2*0");
  static EvtId K2STB=EvtPDL::getId("anti-K_2*0");

  static EvtId PHI=EvtPDL::getId("phi");
  static EvtId DSP=EvtPDL::getId("D_s+");
  static EvtId DSM=EvtPDL::getId("D_s-");

  static EvtId DSSTP=EvtPDL::getId("D_s*+");
  static EvtId DSSTM=EvtPDL::getId("D_s*-");
  static EvtId DS1P=EvtPDL::getId("D_s1+");
  static EvtId DS1M=EvtPDL::getId("D_s1-");
  static EvtId DS0STP=EvtPDL::getId("D_s0*+");
  static EvtId DS0STM=EvtPDL::getId("D_s0*-");
  static EvtId DPS1P=EvtPDL::getId("D'_s1+");
  static EvtId DPS1M=EvtPDL::getId("D'_s1-");
  static EvtId DS2STP=EvtPDL::getId("D_s2*+");
  static EvtId DS2STM=EvtPDL::getId("D_s2*-");


EvtId parnum,mesnum,lnum;

parnum = getParentId();
mesnum = getDaug(0);
lnum = getDaug(1);


if ( parnum==BP||parnum==BM||parnum==B0||parnum==B0B||parnum==BS0||parnum==BSB ) {

  if ( mesnum==DST0||mesnum==DSTP||mesnum==DSTB||mesnum==DSTM||mesnum==DSSTP||mesnum==DSSTM) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(10000.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(7000.0);
       return;
    }
  }


  if ( mesnum==D0||mesnum==DP||mesnum==D0B||mesnum==DM||mesnum==DSP||mesnum==DSM) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(4000.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(3500.0);
       return;
    }
  }


  if ( mesnum==D1P1P||mesnum==D1P1N||mesnum==D1P10||mesnum==D1P1B||mesnum==DS1P||mesnum==DS1M) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(1300.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(480.0);
       return;
    }
  }

  if ( mesnum==D3P1P||mesnum==D3P1N||mesnum==D3P10||mesnum==D3P1B||mesnum==DS0STP||mesnum==DS0STM) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(450.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
      setProbMax(73.0);//???
       return;
    }
  }

  if ( mesnum==D3P0P||mesnum==D3P0N||mesnum==D3P00||mesnum==D3P0B||mesnum==DPS1P||mesnum==DPS1M) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(200.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(90.0);
       return;
    }
  }
  if ( mesnum==D3P2P||mesnum==D3P2N||mesnum==D3P20||mesnum==D3P2B||mesnum==DS2STP||mesnum==DS2STM) {

    if ( mesnum==DS2STP|| mesnum==DS2STM) {
      setProbMax(550.0);
      return;
    }
    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(400.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(220.0);
       return;
    }
  }

  if ( mesnum==D21S0P||mesnum==D21S0N||mesnum==D21S00||mesnum==D21S0B) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(16.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(3.0);
       return;
    }
  }

  if ( mesnum==D23S1P||mesnum==D23S1N||mesnum==D23S10||mesnum==D23S1B) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(500.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(250.0);
       return;
    }
  }

  if ( mesnum==RHOP||mesnum==RHOM||mesnum==RHO0) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(6500.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(6000.0);
       return;
    }
  }

  if ( mesnum==OMEG) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(6800.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(6000.0);
       return;
    }
  }

  if ( mesnum==PIP||mesnum==PIM||mesnum==PI0) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(1200.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(1150.0);
       return;
    }
  }

  if ( mesnum==ETA) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(1800.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(1900.0);
       return;
    }
  }

  if ( mesnum==ETAPR) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(3000.0);
       return;
    } 
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(3000.0);
       return;
    }
  }


  if ( mesnum==B1P||mesnum==B1M||mesnum==B10) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(2500.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(1700.0);
       return;
    }
  }

  if ( mesnum==A0P||mesnum==A0M||mesnum==A00) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(80.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(62.0);
       return;
    }
  }

  if ( mesnum==A1P||mesnum==A1M||mesnum==A10) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(4500.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(3500.0);
       return;
    }
  }

  if ( mesnum==A2P||mesnum==A2M||mesnum==A20) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(1200.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(1000.0);
       return;
    }
  }

  if ( mesnum==H1) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(2600.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(2900.0);
       return;
    }
  }

  if ( mesnum==H1PR) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(1400.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(1500.0);
       return;
    }
  }

  if ( mesnum==F2) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(1100.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(1100.0);
       return;
    }
  }

  if ( mesnum==F2PR) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(804.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(600.0);
       return;
    }
  }

  if ( mesnum==F1) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(2500.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(2000.0) ;
       return;
    }
  }

  if ( mesnum==F1PR) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(2400.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(1700.0);
       return;
    }
  }

  if ( mesnum==F0) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 80.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(63.0);
       return;
    }
  }

  if ( mesnum==F0PR) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(120.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(120.0);
       return;
    }
  }


  if ( mesnum==RHO2SP||mesnum==RHO2SM||mesnum==RHO2S0) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 2400.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(2000.0);
       return;
    }
  }

  if ( mesnum==OMEG2S) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(1600.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(1400.0) ;
       return;
    }
  }

  if ( mesnum==PI2SP||mesnum==PI2SM||mesnum==PI2S0) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 500.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(300.0);
       return;
    }
  }

  if ( mesnum==ETA2S) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(344.0);
       return;
    }
    if ( lnum==TAUP||lnum==TAUM ) {
       setProbMax(300.0);
       return;
    }
  }

   if ( mesnum==KP||mesnum==KM||
        mesnum==K1P||mesnum==K1M||mesnum==K1STP||mesnum==K1STM) {

     if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
        setProbMax(2000.0);
        return;
     }
     if ( lnum==TAUP||lnum==TAUM ) {
        setProbMax(1000.0);
        return;
     }
   }

   if ( mesnum==KSTP||mesnum==KSTM ) {

     if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
        setProbMax(10000.0);
        return;
     }
     if ( lnum==TAUP||lnum==TAUM ) {
        setProbMax(7000.0);
        return;
     }
   }


}

if ( parnum==D0||parnum==DP||parnum==DM||parnum==D0B ) {


  if ( mesnum==RHOP||mesnum==RHOM||mesnum==RHO0) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(110.0);
       return;
    }
  }

  if ( mesnum==OMEG) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(75.0);
       return;
    }
  }

  if ( mesnum==PIP||mesnum==PIM||mesnum==PI0) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(40.0);
       return;
    }
  }

  if ( mesnum==ETA) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 65.0);
       return;
    }
  }

  if ( mesnum==ETAPR) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 60.0);
       return;
    }
  }

  if ( mesnum==KP||mesnum==KM||mesnum==K0||
       mesnum==K0S||mesnum==K0L||mesnum==KB) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 70.0);
       return;
    }
  }

  if ( mesnum==K1STP||mesnum==K1STM||mesnum==K1ST0||mesnum==K1STB) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 3.3);
       return;
    }
  }

  if ( mesnum==K1P||mesnum==K1M||mesnum==K10||mesnum==K1B) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 100.0);
       return;
    }
  }

  if ( mesnum==KSTP||mesnum==KSTM||mesnum==KST0||mesnum==KSTB) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 135.0);
       return;
    }
  }

  if ( mesnum==K2STP||mesnum==K2STM||mesnum==K2ST0||mesnum==K2STB) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
      //Lange - Oct 26,2001 - increasing from 0.75 to 
      //accomodate
      setProbMax( 9.0);
      // setProbMax( 0.75);
       return;
    }
  }

  if ( mesnum==F0) {
    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax(1.0);
       return;
    }
  }


}

if ( parnum==DSP||parnum==DSM ) {


  if ( mesnum==PHI ) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 90.0 );
       return;
    }
  }

  if ( mesnum==ETA ) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 75.0 );
       return;
    }
  }

  if ( mesnum==ETAPR ) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 80.0) ;
       return;
    }
  }

  if ( mesnum==KST0||mesnum==KSTB ) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
      setProbMax( 100.0) ;
       return;	
    }
  }


  if ( mesnum==K0 || mesnum==KB || mesnum==K0S || mesnum==K0L ) {

    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
       setProbMax( 45.0 );
       return;
    }
  }

  if ( mesnum==F0) {
    if ( lnum==EP||lnum==EM||lnum==MUP||lnum==MUM ) {
      setProbMax(1.0);
      return;
    }
  }


}

if ( parnum==BCP||parnum==BCM ) {
   setProbMax(1000.0 );
   return;
}



//This is a real cludge.. (ryd)
 setProbMax(0.0);

}

void EvtISGW2::init(){

  checkNArg(0);
  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton neutrino

  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

  EvtSpinType::spintype mesontype=EvtPDL::getSpinType(getDaug(0));

  isgw2ffmodel = new EvtISGW2FF;
  
  if ( mesontype==EvtSpinType::SCALAR ) { 
    calcamp = new EvtSemiLeptonicScalarAmp; 
  }
  if ( mesontype==EvtSpinType::VECTOR ) { 
    calcamp = new EvtSemiLeptonicVectorAmp; 
  }
  if ( mesontype==EvtSpinType::TENSOR ) { 
    calcamp = new EvtSemiLeptonicTensorAmp; 
  }
  
}







