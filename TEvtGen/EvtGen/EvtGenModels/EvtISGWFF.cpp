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
// Module: EvtISGWFF.cc
//
// Description: Routine to implement semileptonic form factors
//              according to the model ISGW
//
// Modification history:
//
//    DJL     April 17, 1998        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtISGWFF.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include <math.h>
#include <stdlib.h>
using std::endl;

void EvtISGWFF::getscalarff(EvtId parent,EvtId daught,
                       double t, double mass, double *fpf,
			    double *f0f ) {

  //added by Lange Jan4,2000
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId D3P0P=EvtPDL::getId("D_0*+");
  static EvtId D3P0N=EvtPDL::getId("D_0*-");
  static EvtId D3P00=EvtPDL::getId("D_0*0");
  static EvtId D3P0B=EvtPDL::getId("anti-D_0*0");

  static EvtId D21S0P=EvtPDL::getId("hi");
  static EvtId D21S0N=EvtPDL::getId("hi");
  static EvtId D21S00=EvtPDL::getId("hi");
  static EvtId D21S0B=EvtPDL::getId("hi");

  static EvtId ETA2S=EvtPDL::getId("eta(2S)");

  static EvtId PI2S0=EvtPDL::getId("pi(2S)0");
  static EvtId PI2SP=EvtPDL::getId("pi(2S)+");
  static EvtId PI2SM=EvtPDL::getId("pi(2S)-");

  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PIM=EvtPDL::getId("pi-");
  static EvtId PI0=EvtPDL::getId("pi0");

  static EvtId A0P=EvtPDL::getId("a_0+");
  static EvtId A0M=EvtPDL::getId("a_0-");
  static EvtId A00=EvtPDL::getId("a_00"); 

  static EvtId F0=EvtPDL::getId("f_0");
  static EvtId F0PR=EvtPDL::getId("f'_0");

  static EvtId ETA=EvtPDL::getId("eta");
  static EvtId ETAPR=EvtPDL::getId("eta'");

  static EvtId KP=EvtPDL::getId("K+");
  static EvtId KM=EvtPDL::getId("K-");
  static EvtId K0=EvtPDL::getId("K0");
  static EvtId KB=EvtPDL::getId("anti-K0");
  static EvtId K0S=EvtPDL::getId("K_S0");
  static EvtId K0L=EvtPDL::getId("K_L0");

  static EvtId K0STP=EvtPDL::getId("K_0*+");
  static EvtId K0STM=EvtPDL::getId("K_0*-");
  static EvtId K0ST0=EvtPDL::getId("K_0*0");
  static EvtId K0STB=EvtPDL::getId("anti-K_0*0");

  static EvtId DSP=EvtPDL::getId("D_s+");
  static EvtId DSM=EvtPDL::getId("D_s-");

  static EvtId D3P0SP=EvtPDL::getId("D_s0*+");
  static EvtId D3P0SN=EvtPDL::getId("D_s0*-");

  double fmf;
  double mb=EvtPDL::getMeanMass(parent);

  if (daught==PI0||daught==PIP||daught==PIM||daught==ETA||
      daught==ETAPR||daught==D0||daught==D0B||daught==DP||
      daught==DM||daught==KP||daught==KM||daught==K0||daught==K0L||
      daught==KB||daught==DSP||daught==DSM||daught==K0S) {

      EvtISGW1FF1S0(parent,daught,t,mass,fpf,&fmf);
  }
      
  if (daught==PI2S0||daught==PI2SP||daught==PI2SM||daught==ETA2S||
      daught==D21S0P||daught==D21S0B||daught==D21S0N||daught==D21S00){
      EvtISGW1FF21S0(parent,daught,t,mass,fpf,&fmf);
  }
  
  if (daught==A00||daught==A0P||daught==A0M||daught==F0||
      daught==F0PR||daught==D3P0P||daught==D3P00||daught==D3P0B||
      daught==D3P0N||daught==K0STM||daught==K0STB||daught==K0STP||
      daught==D3P0SP||daught==D3P0SN||
      daught==K0ST0) {
    EvtISGW1FF3P0(parent,daught,t,mass,fpf,&fmf);
  }

  *f0f = (fmf/((mb*mb-mass*mass)/t))+(*fpf);

  return ;
}

 void EvtISGWFF::gettensorff(EvtId parent,EvtId daught,
                       double t, double mass, double *hf,
			     double *kf, double *bpf, double *bmf ){

  //added by Lange Jan4,2000
   EvtISGW1FF3P2(parent,daught,t,mass,hf,kf,bpf,bmf);
   
   return;

 }


 void EvtISGWFF::getvectorff(EvtId parent,EvtId daught,
                       double t, double mass, double *a1f,
			     double *a2f, double *vf, double *a0f ){
 
  //added by Lange Jan4,2000
  static EvtId DST0=EvtPDL::getId("D*0");
  static EvtId DSTB=EvtPDL::getId("anti-D*0");
  static EvtId DSTP=EvtPDL::getId("D*+");
  static EvtId DSTM=EvtPDL::getId("D*-");

  static EvtId D1P1P=EvtPDL::getId("D_1+");
  static EvtId D1P1N=EvtPDL::getId("D_1-");
  static EvtId D1P10=EvtPDL::getId("D_10");
  static EvtId D1P1B=EvtPDL::getId("anti-D_10");

  static EvtId D3P1P=EvtPDL::getId("D'_1+");
  static EvtId D3P1N=EvtPDL::getId("D'_1-");
  static EvtId D3P10=EvtPDL::getId("D'_10");
  static EvtId D3P1B=EvtPDL::getId("anti-D'_10");

  static EvtId D23S1P=EvtPDL::getId("hi");
  static EvtId D23S1N=EvtPDL::getId("hi");
  static EvtId D23S10=EvtPDL::getId("hi");
  static EvtId D23S1B=EvtPDL::getId("hi");

  static EvtId RHO2S0=EvtPDL::getId("rho(2S)0");
  static EvtId RHO2SP=EvtPDL::getId("rho(2S)+");
  static EvtId RHO2SM=EvtPDL::getId("rho(2S)-");
  static EvtId OMEG2S=EvtPDL::getId("omega(2S)");

  static EvtId RHOP=EvtPDL::getId("rho+");
  static EvtId RHOM=EvtPDL::getId("rho-");
  static EvtId RHO0=EvtPDL::getId("rho0");

  static EvtId A1P=EvtPDL::getId("a_1+");
  static EvtId A1M=EvtPDL::getId("a_1-");
  static EvtId A10=EvtPDL::getId("a_10");

  static EvtId B1P=EvtPDL::getId("b_1+");
  static EvtId B1M=EvtPDL::getId("b_1-");
  static EvtId B10=EvtPDL::getId("b_10");

  static EvtId H1=EvtPDL::getId("h_1");
  static EvtId H1PR=EvtPDL::getId("h'_1");

  static EvtId F1=EvtPDL::getId("f_1");
  static EvtId F1PR=EvtPDL::getId("f'_1");

  static EvtId OMEG=EvtPDL::getId("omega");

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

  static EvtId PHI=EvtPDL::getId("phi");

  static EvtId D1P1SP=EvtPDL::getId("D_s1+");
  static EvtId D1P1SN=EvtPDL::getId("D_s1-");

  static EvtId D3P1SP=EvtPDL::getId("D'_s1*+");
  static EvtId D3P1SN=EvtPDL::getId("D'_s1*-");

  static EvtId DSSTP=EvtPDL::getId("D_s*+");
  static EvtId DSSTM=EvtPDL::getId("D_s*-");

   double ff,gf,apf,amf;

   if (daught==DST0||daught==DSTP||daught==DSTM||daught==DSTB||
       daught==OMEG||daught==RHO0||daught==RHOM||daught==RHOP||
       daught==KSTP||daught==KSTM||daught==KST0||daught==KSTB||
       daught==PHI||daught==DSSTP||daught==DSSTM) {
     EvtISGW1FF3S1(parent,daught,t,mass,&ff,&gf,&apf,&amf);
   }
   if (daught==B10||daught==B1P||daught==B1M||daught==H1||
       daught==H1PR||daught==D1P1P||daught==D1P10||daught==D1P1B||
       daught==D1P1SP||daught==D1P1SN||
       daught==D1P1N||daught==K10||daught==K1B||daught==K1P||
       daught==K1M) {
     EvtISGW1FF1P1(parent,daught,t,mass,&ff,&gf,&apf,&amf);
   }
   if (daught==RHO2S0||daught==RHO2SP||daught==RHO2SM||daught==OMEG2S||
       daught==D23S1P||daught==D23S1B||daught==D23S1N||daught==D23S10){
     EvtISGW1FF23S1(parent,daught,t,mass,&ff,&gf,&apf,&amf);
   }
   if (daught==A10||daught==A1P||daught==A1M||daught==F1||
                  daught==F1PR||daught==D3P1P||daught==D3P10||daught==D3P1B||
       daught==D3P1N||daught==K1STM||daught==K1STB||daught==K1STP||
       daught==D3P1SP||daught==D3P1SN||
       daught==K1ST0) {
     EvtISGW1FF3P1(parent,daught,t,mass,&ff,&gf,&apf,&amf);
   }

   // Need to stuff in some factors to make these the ffs that
   // is used elsewhere...

   double mb=EvtPDL::getMeanMass(parent);
  

   *vf = (gf)*(mb+mass);
   *a1f = (ff)/(mb+mass);
   *a2f = -1.0*(apf)*(mb+mass);
   double a3f = ((mb+mass)/(2.0*mass))*(*a1f) -
        ((mb-mass)/(2.0*mass))*(*a2f);
   
   *a0f = a3f - ( (t*amf)/(2.0*mass));

   return;
 }

void EvtISGWFF::EvtISGW1FF3P2 (EvtId parent,EvtId daugt,
       double t, double mass,
       double *hf, double *kf, double *bpf, double *bmf ) {

  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D3P2P=EvtPDL::getId("D_2*+");
  static EvtId D3P2N=EvtPDL::getId("D_2*-");
  static EvtId D3P20=EvtPDL::getId("D_2*0");
  static EvtId D3P2B=EvtPDL::getId("anti-D_2*0");

  static EvtId A2P=EvtPDL::getId("a_2+");
  static EvtId A2M=EvtPDL::getId("a_2-");
  static EvtId A20=EvtPDL::getId("a_20");

  static EvtId F2=EvtPDL::getId("f_2");
  static EvtId F2PR=EvtPDL::getId("f'_2");

  double mtb;
  double msd(0.0), mx(0.0), mb(0.0); 
  double msq(0.0), bx2(0.0),mtx,f5;
  double mum,mup,tm,bb2(0.0),bbx2;
  double msb(0.0), kap;

  if (parent==BM||parent==BP||parent==B0||parent==B0B) { 
    msb=5.2;
    msd=0.33;
    bb2=0.41*0.41;
    if (daugt==A20||daugt==A2P||daugt==A2M||daugt==F2||daugt==F2PR) {
      msq=0.33;
      bx2=0.27*0.27;
    }
    else{
      if (daugt==D3P2P||daugt==D3P2N||daugt==D3P2B||daugt==D3P20) {
	msq=1.82;
	bx2=0.34*0.34;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_3P1.\n";
      }
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_3P1.\n";
  }
  
  mtb = msb + msd;
  mtx = msq + msd;

  mb = EvtPDL::getMeanMass( parent );
  mx = mass;

  mup=1.0/(1.0/msq+1.0/msb);
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);

  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;
  kap = 0.7*0.7;

  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0)*
       exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));

  *hf = f5*(msd/(sqrt(8.0*bb2)*mtb))*((1.0/msq)-(msd*bb2/(2.0*mum*
        mtx*bbx2)));
  
  *kf = f5*msd*sqrt(2.0/bb2);

  *bpf = (-1.0*f5*msd/(sqrt(8.0*bb2)*msb*mtx))*(1.0-(msd*msb*bx2/(
         2.0*mup*mtb*bbx2))+(msd*msb*bx2*(1.0-(msd*bx2/(2.0*mtb*bbx2)))/
         (4.0*mtb*mum*bbx2)));
  *bmf = 0.0;
  return;
} //get_ff_isgw_1p1

void EvtISGWFF::EvtISGW1FF1S0 ( EvtId parent, EvtId daugt,
       double t, double mass, double *fpf, double *fmf ) {
  
  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PIM=EvtPDL::getId("pi-");
  static EvtId PI0=EvtPDL::getId("pi0");

  static EvtId ETA=EvtPDL::getId("eta");
  static EvtId ETAPR=EvtPDL::getId("eta'");

  double mtb;
  double msd(0.0), mx(0.0), mb(0.0); 
  double msq(0.0), bx2(0.0),mtx;
  double f3,kap; 
  double msb(0.0),bb2(0.0),mup,mum,bbx2,tm;

  if (parent==BM||parent==BP||parent==B0||parent==B0B) {
    msb=5.2;
    msd=0.33;
    bb2=0.41*0.41;
    if (daugt==PIP||daugt==PIM||daugt==PI0||daugt==ETA||daugt==ETAPR) {
      msq=0.33;
      bx2=0.31*0.31;
    }
    else{
      if (daugt==D0||daugt==DP||daugt==DM||daugt==D0B) {      
	msq=1.82;
	bx2=0.39*0.39;
      }      
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_1S0.\n";
      }
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_1S0.\n";
    report(Severity::Error,"EvtGen") << "Parent:"<<parent.getId()<<endl;
  }
  
  mtb = msb + msd;
  mtx = msq + msd;
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;
  mup=1.0/(1.0/msq+1.0/msb);
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if ( t>tm ) t=0.99*tm;
  
  kap = 0.7*0.7;
  f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,3.0/2.0)*
    exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));

  *fpf = f3*(1+(msb/(2.0*mum))-(msb*msq*msd*bb2/(4.0*mup*mum*mtx*bbx2)));
  *fmf = f3*(1.0-(mtb+mtx)*(0.5/msq-(msd*bb2/(4.0*mup*mtx*bbx2))));

  return;
} //get_ff_isgw_1s0



void  EvtISGWFF::EvtISGW1FF3S1(EvtId parent,EvtId daugt,double t,
      double mass, double *f,double *g,double *ap,double *am){
 
  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId DST0=EvtPDL::getId("D*0");
  static EvtId DSTB=EvtPDL::getId("anti-D*0");
  static EvtId DSTP=EvtPDL::getId("D*+");
  static EvtId DSTM=EvtPDL::getId("D*-");

  static EvtId RHOP=EvtPDL::getId("rho+");
  static EvtId RHOM=EvtPDL::getId("rho-");
  static EvtId RHO0=EvtPDL::getId("rho0");

  static EvtId OMEG=EvtPDL::getId("omega");
 
  double msd(0.0),msq(0.0),bb2(0.0),mum,mtx,bbx2;
  double bx2(0.0),msb(0.0),tm;
  double mb,mx,f3, kap;

  if (parent==BM||parent==BP||parent==B0||parent==B0B) { 
    msb=5.2;
    msd=0.33;
    bb2=0.41*0.41;
    if (daugt==DSTP||daugt==DSTM||daugt==DSTB||daugt==DST0) {
      msq=1.82;
      bx2=0.39*0.39;
    }
    else{
      if (daugt==RHOP||daugt==RHOM||daugt==RHO0||daugt==OMEG) {
	msq=0.33;
	bx2=0.31*0.31;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_3S1.\n";
      }
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_3S1.\n";
  }
  
  double mtb;

  mtb=msb+msd;
  mtx=msq+msd;
  
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  mb=EvtPDL::getMeanMass(parent);
  mx=mass;
  tm=(mb-mx)*(mb-mx);
  if ( t > tm ) t = 0.99*tm;

  kap = 0.7*0.7;
  f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,3.0/2.0)*
       exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));
  
  *f=2.0*mtb*f3;
  *g=0.5*f3*((1/msq)-(msd*bb2/(2.0*mum*mtx*bbx2)));
  *ap=(-1.0*f3/(2.0*mtx))*(1.0+(msd*(bb2-bx2)/(msb
      *(bb2+bx2)))-(msd*msd*bx2*bx2/(4.0*mum*mtb*bbx2*bbx2)));
  *am=0.0;

}

void EvtISGWFF::EvtISGW1FF23S1 (EvtId parent,EvtId daugt,
       double t, double mass, double *fpf, double *gpf, 
       double *appf, double *apmf ) {
  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D23S1P=EvtPDL::getId("hi");
  static EvtId D23S1N=EvtPDL::getId("hi");
  static EvtId D23S10=EvtPDL::getId("hi");
  static EvtId D23S1B=EvtPDL::getId("hi");

  static EvtId RHO2S0=EvtPDL::getId("rho(2S)0");
  static EvtId RHO2SP=EvtPDL::getId("rho(2S)+");
  static EvtId RHO2SM=EvtPDL::getId("rho(2S)-");
  static EvtId OMEG2S=EvtPDL::getId("omega(2S)");

  double mtb;
  double msd(0.0), mx(0.0), mb(0.0); 
  double msq(0.0), bx2(0.0),mtx;
  double f3,f5,tt;
  double mum,bb2(0.0),bbx2,tm,msb(0.0);
    
  if (parent==BM||parent==BP||parent==B0||parent==B0B) {
    msb=5.2;
    msd=0.33;
    bb2=0.41*0.41;
    if (daugt==RHO2SP||daugt==RHO2SM||daugt==RHO2S0||daugt==OMEG2S) {
      msq=0.33;
      bx2=0.31*0.31;
    }
    else{
      if (daugt==D23S1N||daugt==D23S10||daugt==D23S1P||daugt==D23S1B) {
      msq=1.82;
      bx2=0.39*0.39;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_23P1.\n";
      }
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_23P1.\n";
  }
  
  mtb = msb + msd;
  mtx = msq + msd;
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;
  
  double kap = 0.7*0.7;
  f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,3.0/2.0)*
    exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));
  
  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0)*
    exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));
  
  *fpf = sqrt(6.0)*f3*mtb*( ((bb2-bx2)/(bb2+bx2)) + ((msd*msd*bx2*
         (tm-t))/(6.0*mtx*mtb*bbx2*kap*bbx2)));
       
  *gpf = sqrt(3.0/8.0)*f3*( ((((bb2-bx2)/(bb2+bx2)) + ((msd*msd*bx2*
         (tm-t))/(6.0*mtx*mtb*bbx2*kap*bbx2)))*
         ((1.0/msq)-((msd*bb2)/(2.0*mum*mtx*bbx2)))) +
         ((msd*bb2*bx2)/(3.0*mum*mtx*bbx2*bbx2)));

  tt = (msd*msd*bx2*(tm-t))/(mtx*mtb*bb2*kap*bbx2);

  *appf = (f5/(sqrt(6.0)*mtx))* ( ((3.0*mtb*bbx2/(2.0*msb*sqrt(bb2*bx2)))*
          (1.0 - ( (msd*msd*msb*bx2*bx2)/(4.0*mtb*mtb*mum*bbx2*bbx2)))) -
          ( (3.0*msd*sqrt(bx2/bb2))/(2.0*msb)) + 
          ( (5.0*msd*sqrt(bx2*bb2)*(1.0 + 0.1*tt))/(2.0*msb*bbx2)) -
          ( (3.0*mtb*sqrt(bb2/bx2)*(1.0 + (tt/6.0)))/(2.0*msb)) +
          ( (7.0*msd*msd*sqrt(bb2/bx2)*bx2*bx2*(1.0 + (tt/14.0))) /
            (8.0*mtb*mum*bbx2*bbx2))); 

  *apmf = 0.0;
  return;
} //get_ff_isgw_23s1


void EvtISGWFF::EvtISGW1FF3P1 (EvtId parent,EvtId daugt,
       double t, double mass,
       double *lf, double *qf, double *cpf, double *cmf ) {
  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D3P1P=EvtPDL::getId("D'_1+");
  static EvtId D3P1N=EvtPDL::getId("D'_1-");
  static EvtId D3P10=EvtPDL::getId("D'_10");
  static EvtId D3P1B=EvtPDL::getId("anti-D'_10");

  static EvtId A1P=EvtPDL::getId("a_1+");
  static EvtId A1M=EvtPDL::getId("a_1-");
  static EvtId A10=EvtPDL::getId("a_10");

  static EvtId F1=EvtPDL::getId("f_1");
  static EvtId F1PR=EvtPDL::getId("f'_1");

  double mtb;
  double msd(0.0), mx(0.0), mb(0.0); 
  double msq(0.0), bx2(0.0),mtx,f5;
  double msb(0.0),bb2(0.0),mum,bbx2,tm;
  double kap;

  if (parent==BM||parent==BP||parent==B0||parent==B0B) {  
    msb=5.2;
    msd=0.33;
    bb2=0.41*0.41;
    if (daugt==A10||daugt==A1P||daugt==A1M||daugt==F1||daugt==F1PR) {
      msq=0.33;
      bx2=0.27*0.27;
    }
    else{
      if (daugt==D3P1P||daugt==D3P1N||daugt==D3P1B||daugt==D3P10) {
      msq=1.82;
      bx2=0.34*0.34;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_3P1.\n";
      }
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_3P1.\n";
  }
  
  mtb = msb + msd;
  mtx = msq + msd;
  
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;
  
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;
  
  kap = 0.7*0.7;
  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0)*
       exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));

  *qf = (f5*msd)/(2.0*mtx*sqrt(bb2));
  
  *lf = -1.0*mtb*sqrt(bb2)*f5*(1/mum+(msd*(tm-t)/(2.0*mtb*
        kap*bb2))*((1.0/msq)-(1.0*msd*bb2/(2.0*mum*mtx*bbx2)))); 

  *cpf = (f5*msd*msb/(4.0*mtb*sqrt(bb2)*mum))*(1.0-(msd*msq*bb2/(
         2.0*mtx*mum*bbx2)));
  *cmf = 0.0;
  return;
} //get_ff_isgw_3p1



void EvtISGWFF::EvtISGW1FF3P0 (EvtId parent,EvtId daugt,
       double t, double mass, double *upf, double *umf ) {
  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D3P0P=EvtPDL::getId("D_0*+");
  static EvtId D3P0N=EvtPDL::getId("D_0*-");
  static EvtId D3P00=EvtPDL::getId("D_0*0");
  static EvtId D3P0B=EvtPDL::getId("anti-D_0*0");

  static EvtId A0P=EvtPDL::getId("a_0+");
  static EvtId A0M=EvtPDL::getId("a_0-");
  static EvtId A00=EvtPDL::getId("a_00");

  static EvtId F0=EvtPDL::getId("f_0");
  static EvtId F0PR=EvtPDL::getId("f'_0");

  double mtb;
  double msd(0.0), mx(0.0), mb(0.0); 
  double msq(0.0), bx2(0.0),mtx;
  double f5;
  double mum,bb2(0.0),bbx2,msb(0.0),tm;

  if (parent==BM||parent==BP||parent==B0||parent==B0B) {    
    msb=5.2;
    msd=0.33;
    bb2=0.41*0.41;
    if (daugt==A00||daugt==A0P||daugt==A0M||daugt==F0||daugt==F0PR) {
      msq=0.33;
      bx2=0.27*0.27;
    }
    else{
      if (daugt==D3P0P||daugt==D3P0N||daugt==D3P0B||daugt==D3P00) {
	msq=1.82;
	bx2=0.34*0.34;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_3P0.\n";
      }
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_3P0.\n";
  }
  
  mtb = msb + msd;
  mtx = msq + msd;
  
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;
  
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;
  
  double kap = 0.7*0.7;
  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0)*
       exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));

  *upf = f5*msd*msq*msb/(sqrt(6.0*bb2)*mtx*mum);
  *umf = 0.0;
  return;
} //get_ff_isgw_3p0



void EvtISGWFF::EvtISGW1FF1P1 (EvtId parent,EvtId daugt,
       double t, double mass,
       double *vf, double *rf, double *spf, double *smf ) {
    //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D1P1P=EvtPDL::getId("D_1+");
  static EvtId D1P1N=EvtPDL::getId("D_1-");
  static EvtId D1P10=EvtPDL::getId("D_10");
  static EvtId D1P1B=EvtPDL::getId("anti-D_10");

  static EvtId B1P=EvtPDL::getId("b_1+");
  static EvtId B1M=EvtPDL::getId("b_1-");
  static EvtId B10=EvtPDL::getId("b_10");

  static EvtId H1=EvtPDL::getId("h_1");
  static EvtId H1PR=EvtPDL::getId("h'_1");

  double mtb;
  double msd(0.0), mx(0.0), mb(0.0); 
  double msq(0.0), bx2(0.0),mtx,f5;
  double mup,mum,kap;
  double msb(0.0),bb2(0.0),bbx2,tm;

  if (parent==BM||parent==BP||parent==B0||parent==B0B) {

    msb=5.2;
    msd=0.33;
    bb2=0.41*0.41;
    if (daugt==H1||daugt==H1PR||daugt==B10||daugt==B1P||daugt==B1M) {
      msq=0.33;
      bx2=0.27*0.27;
    }
    else{
      if (daugt==D1P1P||daugt==D1P1N||daugt==D1P10||daugt==D1P1B) {
	msq=1.82;
	bx2=0.34*0.34;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_3P1.\n";
      }
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_3P1.\n";
  }
  
  mtb = msb + msd;
  mtx = msq + msd;

  mb = EvtPDL::getMeanMass( parent );
  mx = mass;

  mup=1.0/(1.0/msq+1.0/msb);
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;

  kap = 0.7*0.7;
  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0)*
       exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));

  *vf = f5*(((mtb*sqrt(bb2))/(4.0*sqrt(2.0)*msb*msq*mtx))); 
  *rf = f5*mtb*sqrt(bb2/2)*((1.0/mup));

  *spf = (f5*msd/(sqrt(2.0*bb2)*mtb))*(1.0+(msb/(2.0*mum))-
         (msb*msq*msd*bb2/(4.0*mup*mum*mtx*bbx2)));
  *smf = 0.0;

  return;
//get_ff_isgw_1p1

}

void EvtISGWFF::EvtISGW1FF21S0 (EvtId parent,EvtId daugt,
       double t, double mass, double *fppf, double *fpmf ) {
  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D21S0P=EvtPDL::getId("hi");
  static EvtId D21S0N=EvtPDL::getId("hi");
  static EvtId D21S00=EvtPDL::getId("hi");
  static EvtId D21S0B=EvtPDL::getId("hi");

  static EvtId ETA2S=EvtPDL::getId("eta(2S)");

  static EvtId PI2S0=EvtPDL::getId("pi(2S)0");
  static EvtId PI2SP=EvtPDL::getId("pi(2S)+");
  static EvtId PI2SM=EvtPDL::getId("pi(2S)-");

  double mtb;
  double msd(0.0), mx(0.0), mb(0.0); 
  double msq(0.0), bx2(0.0),mtx;
  double f3;
  double msb(0.0);
  double mum,mup,tm,bb2(0.0),bbx2;
  
  if (parent==BM||parent==BP||parent==B0||parent==B0B) {
    msb=5.2;
    msd=0.33;
    bb2=0.41*0.41;
    if (daugt==PI2S0||daugt==PI2SP||daugt==PI2SM||daugt==ETA2S){
      msq=0.33;
      bx2=0.31*0.31;
    }
    else{
      if (daugt==D21S00||daugt==D21S0P||daugt==D21S0N||daugt==D21S0B) {
	msq=1.82;
	bx2=0.39*0.39;
      }
      else{  
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw1_ff_21S0.\n";
      }
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw1_ff_21S0.\n";
  }

  mtb = msb + msd;
  mtx = msq + msd;
  
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;
  
  mup=1.0/(1.0/msq+1.0/msb);
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;

  double kap = 0.7*0.7;
  f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,3.0/2.0)*
       exp(-1.0*((msd*msd*(tm-t)/(4.0*mtb*mtx*kap*bbx2))));
  
  *fppf = f3*sqrt(3.0/8.0)*(msb/mup)*( ((bb2-bx2)/(bb2+bx2)) +
          (((msq*msd*bb2)/(3.0*mum*mtx*bbx2))*((7.0*bx2-3.0*bb2)/
          (4.0*bbx2))) + 
          (((msd*msd*bx2*(tm-t))/(6.0*mtx*mtb*bbx2*kap*bbx2))*
          (1.0 - ((msq*msd*bb2)/(2.0*mum*mtx*bbx2)))));

  *fpmf = 0.0;
  return;
} //get_ff_isgw_21s0


void EvtISGWFF::getbaryonff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtISGWFF.\n";  
  ::abort();

}

void EvtISGWFF::getdiracff(EvtId, EvtId, double, double, double*, double*,
			   double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtISGWFF.\n";
  ::abort();

}

void EvtISGWFF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
			    double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtISGWFF.\n";
  ::abort();

}
