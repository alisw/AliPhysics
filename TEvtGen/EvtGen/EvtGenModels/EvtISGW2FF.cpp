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
// Module: EvtISGW2FF.cc
//
// Description: Routine to implement semileptonic form factors
//              according to the model ISGW2
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
#include "EvtGenModels/EvtISGW2FF.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtConst.hh"
#include <string>
#include <math.h>
#include <stdlib.h>
using std::endl;

void EvtISGW2FF::getscalarff(EvtId parent,EvtId daught,
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

  static EvtId D21S0P=EvtPDL::getId("D(2S)+");
  static EvtId D21S0N=EvtPDL::getId("D(2S)-");
  static EvtId D21S00=EvtPDL::getId("D(2S)0");
  static EvtId D21S0B=EvtPDL::getId("anti-D(2S)0");

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
      daught==DM||daught==KP||daught==KM||daught==K0||daught==K0S||
      daught==K0L||daught==KB||daught==DSP||daught==DSM) {

      EvtISGW2FF1S0(parent,daught,t,mass,fpf,&fmf);
  }
      
  if (daught==PI2S0||daught==PI2SP||daught==PI2SM||daught==ETA2S||
      daught==D21S0P||daught==D21S0B||daught==D21S0N||daught==D21S00){
      EvtISGW2FF21S0(parent,daught,t,mass,fpf,&fmf);
  }
  
  if (daught==A00||daught==A0P||daught==A0M||daught==F0||
      daught==F0PR||daught==D3P0P||daught==D3P00||daught==D3P0B||
      daught==D3P0N||daught==K0STM||daught==K0STB||daught==K0STP||
      daught==D3P0SP||daught==D3P0SN||
      daught==K0ST0) {
    EvtISGW2FF3P0(parent,daught,t,mass,fpf,&fmf);
  }

  *f0f = (fmf/((mb*mb-mass*mass)/t))+(*fpf);

  return ;
}

 void EvtISGW2FF::gettensorff(EvtId parent,EvtId daught,
                       double t, double mass, double *hf,
			     double *kf, double *bpf, double *bmf ){

  //added by Lange Jan4,2000
    EvtISGW2FF3P2(parent,daught,t,mass,hf,kf,bpf,bmf);

    return;

 }


 void EvtISGW2FF::getvectorff(EvtId parent,EvtId daught,
                       double t, double mass, double *a1f,
			     double *a2f, double *vf, double *a0f ){
   double ff,gf,apf,amf;

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

  static EvtId D23S1P=EvtPDL::getId("D*(2S)+");
  static EvtId D23S1N=EvtPDL::getId("D*(2S)-");
  static EvtId D23S10=EvtPDL::getId("D*(2S)0");
  static EvtId D23S1B=EvtPDL::getId("anti-D*(2S)0");

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

  static EvtId D3P1SP=EvtPDL::getId("D'_s1+");
  static EvtId D3P1SN=EvtPDL::getId("D'_s1-");

  static EvtId DSSTP=EvtPDL::getId("D_s*+");
  static EvtId DSSTM=EvtPDL::getId("D_s*-");

   if (daught==DST0||daught==DSTP||daught==DSTM||daught==DSTB||
       daught==OMEG||daught==RHO0||daught==RHOM||daught==RHOP||
       daught==KSTP||daught==KSTM||daught==KST0||daught==KSTB||
       daught==PHI||daught==DSSTP||daught==DSSTM) {
     EvtISGW2FF3S1(parent,daught,t,mass,&ff,&gf,&apf,&amf);
   }
   if (daught==B10||daught==B1P||daught==B1M||daught==H1||
       daught==H1PR||daught==D1P1P||daught==D1P10||daught==D1P1B||
       daught==D1P1SP||daught==D1P1SN||
       daught==D1P1N||daught==K10||daught==K1B||daught==K1P||
       daught==K1M) {
     EvtISGW2FF1P1(parent,daught,t,mass,&ff,&gf,&apf,&amf);
   }
   if (daught==RHO2S0||daught==RHO2SP||daught==RHO2SM||daught==OMEG2S||
       daught==D23S1P||daught==D23S1B||daught==D23S1N||daught==D23S10){
     EvtISGW2FF23S1(parent,daught,t,mass,&ff,&gf,&apf,&amf);
   }
   if (daught==A10||daught==A1P||daught==A1M||daught==F1||
                  daught==F1PR||daught==D3P1P||daught==D3P10||daught==D3P1B||
       daught==D3P1N||daught==K1STM||daught==K1STB||daught==K1STP||
       daught==D3P1SP||daught==D3P1SN||
       daught==K1ST0) {
     EvtISGW2FF3P1(parent,daught,t,mass,&ff,&gf,&apf,&amf);
   }

   // Need to stuff in some factors to make these the ffs that
   // is used elsewhere...

   double mb=EvtPDL::getMeanMass(parent);
  

   *vf = (gf)*(mb+mass);
   *a1f = (ff)/(mb+mass);
   *a2f = -1.0*(apf)*(mb+mass);

   double a3f = ((mb+mass)/(2.0*mass))*(*a1f) -
        ((mb-mass)/(2.0*mass))*(*a2f);

   *a0f = a3f + ( (t*amf)/(2.0*mass));

   return;
 }



void EvtISGW2FF::EvtISGW2FF1S0 (EvtId parent,EvtId daugt,
			  double t, double mass, double *fpf, double *fmf ) {

  double mtb, mbb(0.0);
  double msd(0.0), mx,mb,nf(0.0),nfp(0.0); 
  double msq(0.0),bx2(0.0),mbx(0.0),mtx;
  double zji,cji,gammaji,chiji,betaji_fppfm;
  double rfppfm,rfpmfm,f3fppfm,f3fpmfm,fppfm,fpmfm,ai,f3; 
  double mqm,msb(0.0),bb2(0.0),mup,bbx2,tm,r2,betaji_fpmfm;

  EvtId prnt=parent;
  EvtId dgt=daugt;

  //added by Lange Jan4,2000
  static EvtIdSet theB("B+","B-","B0","anti-B0");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PIM=EvtPDL::getId("pi-");
  static EvtId PI0=EvtPDL::getId("pi0");

  static EvtId ETA=EvtPDL::getId("eta");
  static EvtId ETAPR=EvtPDL::getId("eta'");

  static EvtId KP=EvtPDL::getId("K+");
  static EvtId KM=EvtPDL::getId("K-");
  static EvtId K0=EvtPDL::getId("K0");
  static EvtId KB=EvtPDL::getId("anti-K0");
  static EvtId K0S=EvtPDL::getId("K_S0");
  static EvtId K0L=EvtPDL::getId("K_L0");

  static EvtId DSP=EvtPDL::getId("D_s+");
  static EvtId DSM=EvtPDL::getId("D_s-");

  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BS0=EvtPDL::getId("B_s0");

  if (theB.contains(prnt)) {
    
    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=5.31;
    nf = 4.0;
   
    if (dgt==PI0||dgt==PIP||dgt==PIM||dgt==ETA||dgt==ETAPR) {

      msq=0.33;
      bx2=0.406*0.406;
      mbx=0.75*0.770+0.25*0.14;
      nfp = 0.0;
    }
    else{
      if (dgt==D0||dgt==D0B||dgt==DP||dgt==DM) {
	msq=1.82;
	bx2=0.45*0.45;
	mbx=0.75*2.01+0.25*1.87;
	nfp = 3.0;
      }
      else{
      report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_1S0.\n";
      }
    }
  }
  else{
    if (prnt==D0||prnt==D0B||prnt==DP||prnt==DM) {
      msb=1.82;
      msd=0.33;
      bb2=0.45*0.45;
      mbb=1.963;
      nf = 3.0;

      if (dgt==PIP||dgt==PIM||dgt==PI0||dgt==ETA||dgt==ETAPR) {
	msq=0.33;
	bx2=0.406*0.406;
	mbx=0.75*0.770+0.25*0.14;
	nfp = 0.0;
      }
      else{
	if (dgt==K0||dgt==K0S||dgt==K0L||dgt==KB||dgt==KP||dgt==KM) {
	  msq=0.55;
	  bx2=0.44*0.44;
	  mbx=0.75*0.892+0.25*0.49767;
	  nfp = 2.0;
	}
	else{
      report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_1S0.\n";
	}
      }
    }
    else{
      if (prnt==DSP||prnt==DSM){
	msb=1.82;
	msd=0.55;
	bb2=0.56*0.56;
	mbb=1.968;
	nf = 3.0;
    
	if  (dgt==K0||dgt==K0S||dgt==K0L||dgt==KB) {
      
	  msq=0.33;
	  bx2=0.44*0.44;
	  mbx=0.75*0.770+0.25*0.14;
	  nfp = 0.0;
	}
	else{
	  if  (dgt==PI0||dgt==ETA||dgt==ETAPR) {
	    msq=0.33;
	    bx2=0.53*0.53;
	    mbx=0.75*0.892+0.25*0.49767;
	    nfp = 0.0;
	  }
	  else{

	    report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_1S0.\n";
	  }
	}
      }
    else{
      //BS -> cs constants added by djl on Jan. 21,1998
      if (prnt==BS0||prnt==BSB){

	msb=5.2;
	msd=0.55;
	bb2=0.54*0.54;
	mbb=5.38;
	nf = 4.0;
    
	if  (dgt==DSP||dgt==DSM) {
      
	  msq=1.82;
	  bx2=0.56*0.56;
	  mbx=0.75*2.11+0.25*1.97;
	  nfp = 3.0;
	}
	else if  (dgt==KP||dgt==KM) {
      
	  msq=0.55;
	  bx2=0.44*0.44;
	  mbx=0.75*0.892+0.25*0.49767;
	  nfp = 2.0;
	}
	else{
	    report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_1S0.\n";
	  }
	}
    
      else{
	report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_1S0.\n";
	report(Severity::Error,"EvtGen") << "Parent:"<<parent.getId()<<endl;
      }
    }
    }
  }

  mtb = msb + msd;
  mtx = msq + msd;
  mb=EvtPDL::getMeanMass(parent);
  mx=mass;
  
  mup=1.0/(1.0/msq+1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if ( t>tm ) t=0.99*tm;
  
  mqm = 0.1;
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm,mqm)/EvtGetas(msq,msq));
  
  f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,1.5) /
    ((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0));
  
//  for w use wt def with physical masses.
//  report(Severity::Error,"EvtGen") << "before w\n";
  
  ai = -1.0* ( 6.0/( 33.0 - 2.0*nf));  
  cji = pow(( EvtGetas( msb,msb ) / EvtGetas( msq,msq ) ),ai);
  
  zji = msq / msb;
  
  gammaji = EvtGetGammaji( zji );
  chiji = -1.0 - ( gammaji / ( 1- zji ));
  betaji_fppfm = gammaji - (2.0/3.0)*chiji;
  betaji_fpmfm = gammaji + (2.0/3.0)*chiji;
  rfppfm = cji *(1.0 + betaji_fppfm*EvtGetas( msq,sqrt(msb*msq) )/EvtConst::pi);
  rfpmfm = cji *(1.0 + betaji_fpmfm*EvtGetas( msq,sqrt(msb*msq) )/EvtConst::pi);
  f3fppfm = f3*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),0.5);
  f3fpmfm = f3*pow(( mbb / mtb ),0.5)*pow((mbx/mtx),-0.5);
  fppfm = f3fppfm* rfppfm * ( 2.0 - ( ( mtx/msq)*(1- ( (msd*msq*bb2)
						       /(2.0*mup*mtx*bbx2)))));
  fpmfm = f3fpmfm* rfpmfm * ( mtb/msq) * ( 1 - ( ( msd*msq*bb2)/
						 ( 2.0*mup*mtx*bbx2)));
  
  *fpf = (fppfm + fpmfm)/2.0;
  *fmf = (fppfm - fpmfm)/2.0;
  
  return;  
} //get_ff_isgw_1s0





void  EvtISGW2FF::EvtISGW2FF3S1(EvtId parent,EvtId daugt,double t,double mass,
      double *f,double *g,double *ap,double *am){

  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId DST0=EvtPDL::getId("D*0");
  static EvtId DSTB=EvtPDL::getId("anti-D*0");
  static EvtId DSTP=EvtPDL::getId("D*+");
  static EvtId DSTM=EvtPDL::getId("D*-");
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId RHOP=EvtPDL::getId("rho+");
  static EvtId RHOM=EvtPDL::getId("rho-");
  static EvtId RHO0=EvtPDL::getId("rho0");
  static EvtId OMEG=EvtPDL::getId("omega");

  static EvtId KSTP=EvtPDL::getId("K*+");
  static EvtId KSTM=EvtPDL::getId("K*-");
  static EvtId KST0=EvtPDL::getId("K*0");
  static EvtId KSTB=EvtPDL::getId("anti-K*0");

  static EvtId PHI=EvtPDL::getId("phi");
  static EvtId DSP=EvtPDL::getId("D_s+");
  static EvtId DSM=EvtPDL::getId("D_s-");

  static EvtId DSSTP=EvtPDL::getId("D_s*+");
  static EvtId DSSTM=EvtPDL::getId("D_s*-");

  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BS0=EvtPDL::getId("B_s0");

  double cf(0.0),mtb,wt,msd(0.0),mup,f3f,msq(0.0),bb2(0.0),mum,mtx,bbx2,f3g;
  double cji,bx2(0.0),f3appam,msb(0.0),tm,mbb(0.0),mbx(0.0);
  double f3apmam,appam,apmam,mb,mx,f3;
  double r_f,r_g,r_apmam, betaji_f,betaji_g;
  double betaji_appam, betaji_apmam;
  double mqm,r2,chiji,zji,ai,nf(0.0),nfp(0.0),gammaji;

  EvtId prnt=parent;
  EvtId dgt=daugt;

  if (parent==B0||parent==B0B||parent==BP||parent==BM) {

    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=5.31;
    nf = 4.0;
    
    if (dgt==DST0||dgt==DSTP||dgt==DSTM||dgt==DSTB) {

      cf=0.989;
      msq=1.82;
      bx2=0.38*0.38;
      mbx=0.75*2.01+0.25*1.87;
      nfp = 3.0;
    }  
    else{
      if (dgt==OMEG||dgt==RHO0||dgt==RHOM||dgt==RHOP) {
	
	cf=0.905;
	msq=0.33;
	bx2=0.299*0.299;
	mbx=0.75*0.770+0.25*0.14;
	nfp = 0.0;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_3S1.\n";
      }
    }
  }
  else{
    if (prnt==D0||prnt==D0B||prnt==DP||prnt==DM) {
      
      msb=1.82;
      msd=0.33;
      bb2=0.45*0.45;
      mbb=1.963;
      nf = 3.0;

      if (dgt==KSTP||dgt==KSTM||dgt==KST0||dgt==KSTB) {
	
	cf=0.928;
	msq=0.55;
	bx2=0.33*0.33;
	mbx=0.75*0.892+0.25*0.494;
	nfp = 2.0;
      }
      else{
	if (dgt==RHO0||dgt==OMEG||dgt==RHOM||dgt==RHOP) {
	  cf=0.889;
	  msq=0.33;
	  bx2=0.299*0.299;
	  mbx=0.75*0.770+0.25*0.14;
	  nfp = 0.0;
	}
	else{
	  report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_3S1.\n";
	}
      }
    }
    else{
      if (prnt==DSP||prnt==DSM){
    
	msb=1.82;
	msd=0.55;
	bb2=0.56*0.56;
	mbb=1.968;
	nf = 3.0;
	
	if  (dgt==KSTB||dgt==KST0) {

	  cf=0.8731;
	  msq=0.55;
	  bx2=0.33*0.33;
	  mbx=0.87;
	  nfp = 2.0;
	}
	else{
	  if(dgt==PHI){
	    cf=0.911;
	    msq=0.55;
	    bx2=0.37*0.37;
	    mbx=0.97;
	    nfp = 2.0;
	  }
	  else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_3S1.\n";
	  }
	}
      }
    else{
      //BS -> cs constants added by djl on Jan. 21,1998
      if (prnt==BS0||prnt==BSB){

	msb=5.2;
	msd=0.55;
	bb2=0.54*0.54;
	mbb=5.38;
	nf = 4.0;
    
	if  (dgt==DSSTP||dgt==DSSTM) {
      
          cf=0.984;
	  msq=1.82;
	  bx2=0.49*0.49;
	  mbx=0.75*2.11+0.25*1.97;
	  nfp = 3.0;
	}
	else if (dgt==KSTP||dgt==KSTM||dgt==KST0||dgt==KSTB) {
	
	  cf=0.928;
	  msq=0.55;
	  bx2=0.33*0.33;
	  mbx=0.75*0.892+0.25*0.494;
	  nfp = 2.0;
	}
	else{
	    report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_1S0.\n";
	  }
	}
      
      else{
	report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw2_ff_3S1.\n";
      }
    }
    }
  }


  mtb=msb+msd;
  mtx=msq+msd;
  
  mup=1.0/(1.0/msq+1.0/msb);
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  mb=EvtPDL::getMeanMass(parent);
  mx=mass;
  tm=(mb-mx)*(mb-mx);
  if ( t > tm ) t = 0.99*tm;

  wt=1.0+(tm-t)/(2.0*mbb*mbx);
  mqm = 0.1;
  
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm,mqm)/EvtGetas(msq,msq));
 
  ai = -1.0* ( 6.0/( 33.0 - 2.0*nf));  
  
  cji = pow(( EvtGetas( msb,msb ) / EvtGetas( msq,msq ) ),ai);
  zji = msq / msb;

  gammaji = EvtGetGammaji( zji );

  chiji = -1.0 - ( gammaji / ( 1- zji ));
  
  betaji_g = (2.0/3.0)+gammaji;
  betaji_f = (-2.0/3.0)+gammaji;
  betaji_appam = -1.0-chiji+(4.0/(3.0*(1.0-zji)))+
                 (2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)));
  
  betaji_apmam = (1.0/3.0)-chiji-(4.0/(3.0*(1.0-zji)))-
                 (2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)))+
                 gammaji;

  r_g = cji*(1+(betaji_g*EvtGetas( msq,sqrt(mb*msq) )/(EvtConst::pi)));
  r_f = cji*(1+(betaji_f*EvtGetas( msq,sqrt(mb*msq) )/(EvtConst::pi)));
  r_apmam = cji*(1+(betaji_apmam*EvtGetas( msq,sqrt(mb*msq) )/(EvtConst::pi)));

  
  f3=sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,1.5)/
    ((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0));
  
  f3f=sqrt(mbx*mbb/(mtx*mtb))*f3;
  f3g=sqrt(mtx*mtb/(mbx*mbb))*f3;
  f3appam=sqrt(mtb*mtb*mtb*mbx/(mbb*mbb*mbb*mtx))*f3;
  f3apmam=sqrt(mtx*mtb/(mbx*mbb))*f3;
  *f=cf*mtb*(1+wt+msd*(wt-1)/(2*mup))*f3f*r_f;
  *g=0.5*(1/msq-msd*bb2/(2*mum*mtx*bbx2))*f3g*r_g;
  
  appam=cji*(msd*bx2*(1-msd*bx2/(2*mtb*bbx2))/ 
	     ((1+wt)*msq*msb*bbx2)-
	     betaji_appam*EvtGetas( msq,sqrt(msq*mb) )/
	     (mtb*EvtConst::pi))*f3appam;
  
  apmam=-1.0*(mtb/msb-msd*bx2/(2*mup*bbx2)+wt*msd*mtb*bx2*
	      (1-msd*bx2/(2*mtb*bbx2))/((wt+1)*msq*msb*bbx2))*
    f3apmam*r_apmam/mtx;
  
  *ap=0.5*(appam+apmam);
  *am=0.5*(appam-apmam);
  return;
}


void EvtISGW2FF::EvtISGW2FF21S0 (EvtId parent,EvtId daugt,
       double t, double mass, double *fppf, double *fpmf ) {

  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId D21S0P=EvtPDL::getId("D(2S)+");
  static EvtId D21S0N=EvtPDL::getId("D(2S)-");
  static EvtId D21S00=EvtPDL::getId("D(2S)0");
  static EvtId D21S0B=EvtPDL::getId("anti-D(2S)0");
  
  static EvtId ETA2S=EvtPDL::getId("eta(2S)");

  static EvtId PI2S0=EvtPDL::getId("pi(2S)0");
  static EvtId PI2SP=EvtPDL::getId("pi(2S)+");
  static EvtId PI2SM=EvtPDL::getId("pi(2S)-");

  double mtb, mbb(0.0);
  double msd(0.0), mx,mb,nfp(0.0); 
  double msq(0.0),bx2(0.0),mbx(0.0),mtx;
  double f3fppfm,f3fpmfm,fppfm,fpmfm,f3;
  double mqm,msb(0.0);
  double r2,wt,tm,bb2(0.0),bbx2;
  double tau,udef,vdef;

  EvtId prnt=parent;
  EvtId dgt=daugt;

  if (prnt==B0||prnt==B0B||prnt==BP||prnt==BM) {
    
    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=0.75*5.325+0.25*5.279;

    if (dgt==PI2S0||dgt==PI2SP||dgt==PI2SM||dgt==ETA2S) {

      msq=0.33;
      bx2=0.406*0.406;
      mbx=0.75*1.45+0.25*1.300;
      nfp = 0.0;
    }
    else{
      if (dgt==D21S0P||dgt==D21S0B||dgt==D21S0N||dgt==D21S00) {
	msq=1.82;
	bx2=0.45*0.45;
	mbx=0.75*2.64+0.25*2.58;
	nfp=3.0;
      }
      else{

      report(Severity::Error,"EvtGen") << "Not implemented daugt in get_EvtISGW2_ff_21S0.\n";
      }
    }
  }
  else{
    if (prnt==D0||prnt==D0B||prnt==DP||prnt==DM) {
      msb=1.82;
      msd=0.33;
      bb2=0.45*0.45;
      mbb=1.963;
      if (dgt==PI2SP||dgt==PI2SM||dgt==PI2S0||dgt==ETA2S) {
	msq=0.33;
	bx2=0.406*0.406;
	mbx=0.75*1.45+0.25*1.300;
	nfp = 0.0;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_EvtISGW2_ff_21S0.\n";
      }
    }
    else{
      report(Severity::Error,"EvtGen") << "Not implemented parent in get_EvtISGW2_ff_21S0.\n";
    }
  }
  
  mtb = msb + msd;
  mtx = msq + msd;
  
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;

  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;
  wt=1.0+(tm-t)/(2.0*mbb*mbx);
  
  mqm = 0.1;
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm)/EvtGetas(msq));
  
  f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,3.0/2.0) /
    (pow((1.0+r2*(tm-t)/24.0),4.0));
  
  f3fppfm = f3*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),0.5);
  f3fpmfm = f3*pow(( mbb / mtb ),0.5)*pow((mbx/mtx),-0.5);
  
  tau = msd*msd*bx2*(wt-1)/(bb2*bbx2);
  udef = (( bb2-bx2)/(2.0*bbx2)) + ((bb2*tau)/(3.0*bbx2));
  vdef = (bb2*(1.0+(msq/msb))/(6.0*bbx2))*(7.0 - ((bb2/bbx2)*(5+tau)));

  fppfm = f3fppfm*sqrt(1.5)*((1.0-(msd/msq))*udef-(msd*vdef/msq));
  fpmfm = f3fpmfm*sqrt(1.5)*(mtb/msq)*(udef+(msd*vdef/mtx));

  *fppf = (fppfm + fpmfm) /2.0;
  *fpmf = (fppfm - fpmfm) /2.0;
  return;

} //get_ff_isgw_21s0


void EvtISGW2FF::EvtISGW2FF23S1 (EvtId parent,EvtId daugt,
       double t, double mass, double *fpf, double *gpf, 
       double *appf, double *apmf ) {

  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId D23S1P=EvtPDL::getId("D*(2S)+");
  static EvtId D23S1N=EvtPDL::getId("D*(2S)-");
  static EvtId D23S10=EvtPDL::getId("D*(2S)0");
  static EvtId D23S1B=EvtPDL::getId("anti-D*(2S)0");

  static EvtId RHO2S0=EvtPDL::getId("rho(2S)0");
  static EvtId RHO2SP=EvtPDL::getId("rho(2S)+");
  static EvtId RHO2SM=EvtPDL::getId("rho(2S)-");
  static EvtId OMEG2S=EvtPDL::getId("omega(2S)");

  double mtb,  mbb(0.0);
  double msd(0.0), mx,mb,nfp(0.0); 
  double msq(0.0),bx2(0.0),mbx(0.0),mtx;
  double f3appam,f3apmam,f3,appam,apmam,f3fp,f3gp;
  double udef,tau,mum,bb2(0.0),bbx2,tm,wt,mqm,r2,msb(0.0);
  double cfp(0.0);

  EvtId prnt=parent;
  EvtId dgt=daugt;

  if (prnt==B0||prnt==B0B||prnt==BP||prnt==BM) {

    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=0.75*5.325+0.25*5.279;

    if (dgt==RHO2S0||dgt==RHO2SP||dgt==RHO2SM||dgt==OMEG2S) {
      
      cfp=0.776;
      msq=0.33;
      bx2=0.299*0.299;
      mbx=0.75*1.45+0.25*1.300;
      nfp = 0.0;

    }
    else{
      if (dgt==D23S1N||dgt==D23S1P||dgt==D23S1B||dgt==D23S10) {
	cfp=0.929;
	msq=1.82;
	bx2=0.38*0.38;
	mbx=0.75*2.64+0.25*2.58;
	nfp=3.0;
      }
      else{
      report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_23P1.\n";
      }
    }
  }
  else{
    if (prnt==D0||prnt==D0B||prnt==DP||prnt==DM) {
      msb=1.82;
      msd=0.33;
      bb2=0.45*0.45;
      mbb=1.963;

      if (dgt==RHO2S0||dgt==RHO2SP||dgt==RHO2SM||dgt==OMEG2S) {
	cfp=0.74;
	msq=0.33;
	bx2=0.299*0.299;
	mbx=0.75*1.45+0.25*1.300;
	nfp = 0.0;
      }
      else{      
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_23P1.\n";
      }
    }
    else{
      report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_23P1.\n";
    }
  }

  mtb = msb + msd;
  mtx = msq + msd;
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;
  
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);

  if (t>tm) t = 0.99*tm;
  wt=1.0+(tm-t)/(2.0*mbb*mbx);
  
  mqm = 0.1;
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
     (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
     log(EvtGetas(mqm)/EvtGetas(msq));

  f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,3.0/2.0) /
       (pow((1.0+r2*(tm-t)/24.0),4.0));
  
  f3fp = f3*pow(( mbb / mtb ),0.5)*pow((mbx/mtx),0.5);
  f3gp = f3*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),-0.5);
  f3appam = f3*pow(( mbb / mtb ),-1.5)*pow((mbx/mtx),0.5);
  f3apmam = f3*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),-0.5);

  tau = msd*msd*bx2*(wt-1.0)/(bb2*bbx2);
  udef = (( bb2-bx2)/(2.0*bbx2));
  udef = udef + ((bb2*tau)/(3.0*bbx2));

  *fpf = cfp*sqrt(1.5)*mtb*(1.0+wt)*udef*f3fp;

  *gpf = sqrt(3.0/8.0)*f3gp*(((1.0/msq)-((msd*bb2)/(2.0*mum*mtx*bbx2)))*
        udef + ( (msd*bb2*bx2)/(3.0*mum*mtx*bbx2*bbx2)));
       
  appam = f3appam*sqrt(2.0/3.0)*(bb2/(msq*msb*bbx2))*((-7.0*msd*msd*bx2*
          bx2*(1.0+(tau/7.0))/(8.0*mtb*bbx2*bbx2))+(5.0*msd*bx2*(1.0+
          (tau/5.0))/(4.0*bbx2))+(3.0*msd*msd*bx2*bx2/(8.0*mtb*bb2*bbx2))-
          (3.0*msd*bx2/(4.0*bb2)));
         
  apmam = f3apmam*sqrt(3.0/2.0)*(mtb/(msb*mtx))*(1.0-(bb2*(1.0+(tau/7.0))/
          bbx2)-(msd*bx2*(1.0-(5.0*bb2*(1.0+(tau/5.0))/(3.0*bbx2)))
          /(2.0*mtb*bbx2))-(7.0*msd*msd*bb2*bx2/(12.0*msq*mtb*bbx2*bbx2))*
          (1.0-(bx2/bbx2)+(bb2*tau/(7.0*bbx2))));

  *appf = (appam + apmam) /2.0;
  *apmf = (appam - apmam) /2.0;
  return;
} //get_ff_isgw_23s1

void EvtISGW2FF::EvtISGW2FF1P1 (EvtId parent,EvtId daugt,
       double t, double mass, double *rf, double *vf, 
       double *spf, double *smf ) {
  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId D1P1P=EvtPDL::getId("D_1+");
  static EvtId D1P1N=EvtPDL::getId("D_1-");
  static EvtId D1P10=EvtPDL::getId("D_10");
  static EvtId D1P1B=EvtPDL::getId("anti-D_10");

  static EvtId B1P=EvtPDL::getId("b_1+");
  static EvtId B1M=EvtPDL::getId("b_1-");
  static EvtId B10=EvtPDL::getId("b_10");

  static EvtId H1=EvtPDL::getId("h_1");
  static EvtId H1PR=EvtPDL::getId("h'_1");

  static EvtId K1P=EvtPDL::getId("K_1+");
  static EvtId K1M=EvtPDL::getId("K_1-");
  static EvtId K10=EvtPDL::getId("K_10");
  static EvtId K1B=EvtPDL::getId("anti-K_10");

  static EvtId D1P1SP=EvtPDL::getId("D_s1+");
  static EvtId D1P1SN=EvtPDL::getId("D_s1-");

  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BS0=EvtPDL::getId("B_s0");

  double mtb, mbb(0.0);
  double msd(0.0), mx,mb,nfp(0.0); 
  double msq(0.0),bx2(0.0),mbx(0.0),mtx,f5;
  double f5sppsm,f5spmsm;
  double f5v,f5r,mup,mum,vv,rr,spmsm,sppsm;
  double mqm,msb(0.0),bb2(0.0),bbx2,tm,wt,r2;
  EvtId prnt=parent;
  EvtId dgt=daugt;
  if (prnt==B0||prnt==B0B||prnt==BP||prnt==BM) {
    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=5.31;
    if (dgt==B10||dgt==B1P||dgt==B1M||dgt==H1||dgt==H1PR) {
      msq=0.33;
      bx2=0.275*0.275;
      mbx=(3.0*1.123+0.98+5.0*1.32+3.0*1.26)/12.0;
      nfp = 0.0;
    }
    else{
      if (dgt==D1P1P||dgt==D1P10||dgt==D1P1B||dgt==D1P1N) {
	msq=1.82;
	bx2=0.33*0.33;
	mbx=(5.0*2.46+3.0*2.42)/8.0;
	nfp = 3.0;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_1P1.\n";
      }
    }
  }
  else{
    if (prnt==DM||prnt==DP||prnt==D0B||prnt==D0) {
      msb=1.82;
      msd=0.33;
      bb2=0.45*0.45;
      mbb=1.963;
      if (dgt==B10||dgt==B1P||dgt==B1M||dgt==H1||dgt==H1PR) {
	msq=0.33;
	bx2=0.275*0.275;
	mbx=(3.0*1.123+0.98+5.0*1.32+3.0*1.26)/12.0;
	nfp = 0.0;
      }
      else{
	if (dgt==K10||dgt==K1B||dgt==K1P||dgt==K1M) {
	  msq=0.55;
	  bx2=0.30*0.30;
	  mbx=(3.0*1.27+1.43+5.0*1.43+3.0*1.4)/12.0;
	  nfp = 2.0;
	}
	else{
	  report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_1P1.\n";
	}
      }
    }
    else{
      //BS -> cs constants added by djl on Jan. 21,1998
      if (prnt==BS0||prnt==BSB){

	msb=5.2;
	msd=0.55;
	bb2=0.54*0.54;
	mbb=5.38;
    
	if  (dgt==D1P1SP||dgt==D1P1SN) {
      
	  msq=1.82;
	  bx2=0.41*0.41;
	  mbx=(5.0*2.61+3.0*2.54)/8.0;
	  nfp = 3.0;
	}
	else if (dgt==K10||dgt==K1B||dgt==K1P||dgt==K1M) {
	  msq=0.55;
	  bx2=0.30*0.30;
	  mbx=(3.0*1.27+1.43+5.0*1.43+3.0*1.4)/12.0;
	  nfp = 2.0;
	}
	else{
	    report(Severity::Error,"EvtGen") << "Not implemented daugt:"
                         <<daugt.getId()<<" in get_isgw_ff_1S0.\n";
	}
      }
    
    else{
      report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_1P1.\n";
    }
    }
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
  wt=1.0+(tm-t)/(2.0*mbb*mbx);
  
  mqm = 0.1;
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2)+
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm,mqm)/EvtGetas(msq,msq));
  
  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0) /
    (pow((1.0+r2*(tm-t)/18.0),3.0));
  
  f5v = f5*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),-0.5);
  f5r = f5*pow(( mbb / mtb ),0.5)*pow((mbx/mtx),0.5);
  f5sppsm = f5*pow(( mbb / mtb ),-1.5)*pow((mbx/mtx),0.5);
  f5spmsm = f5*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),-0.5);
  
  if (msq == msd) { 
    vv = f5v*(((mtb*sqrt(bb2))/(4.0*sqrt(2.0)*msb*msq*mtx)) +
              (((wt-1)*msd)/(6.0*sqrt(2.0*bb2)*mtx)));
    
    rr = f5r*mtb*sqrt(bb2/2)*((1.0/mup)+((msd*mtx*(wt-1)*(wt-1))/
				 (3.0*msq*bb2)));
    
    sppsm = msd*f5sppsm/(sqrt(2.0*bb2)*mtb)*(1.0-(msd/msq)+((msd*bb2)/
				    (2.0*mup*bbx2)));
    
    spmsm = msd*f5spmsm/(sqrt(2.0*bb2)*msq)*(((4-wt)/3.0)- ( (msd*msq*bb2)/
				     (2.0*mtx*mup*bbx2)));
    
  } else {
    vv = -1.0*msd*f5v/(2.0*sqrt(3.0*bb2)*mtx)*
      ((wt+1)/2.0+bb2*mtb/(2.0*msd*msq*msb));
    
    rr = -2.0*mtb*sqrt(bb2/3.0)*f5r*(1.0/msq + mtx*msd*(wt-1)/(2.0*bb2)*
		     ((wt+1)/(2.0*msq)-msd*bb2/(2.0*mum*mtx*bbx2)));
    
    sppsm = -1.0*sqrt(3.0)*msd*f5sppsm/(2.0*sqrt(bb2)*mtb)*(1 - msd/(3.0*msq) -
			    msd*bb2/(3.0*bbx2)*(1.0/(2.0*mum)-1.0/mup));
    
    spmsm = -1.0*msd*f5spmsm/(2.0*sqrt(3.0*bb2)*mtx)*((2-wt)*mtx/msq +
			   msd*bb2/bbx2*(1.0/(2.0*mum)-1.0/mup));        

  }

  //smooth out the mass(meson) dependence a little
  double parMass=EvtPDL::getMeanMass(prnt);
  double q2max = parMass*parMass + mass*mass - 2.0*parMass*mass;
  double massNom= EvtPDL::getMeanMass(dgt);
  double q2maxNom = parMass*parMass + massNom*massNom - 2.0*parMass*massNom;
  double q2maxin=sqrt(q2maxNom/q2max);
  if ( q2maxin > 1000. ) q2maxin=1000.;

  vv*=q2maxin;
  rr*=q2maxin;
  sppsm*=q2maxin;
  spmsm*=q2maxin;

  *vf = vv;
  *rf = rr;
  *spf = (sppsm + spmsm)/2.0;
  *smf = (sppsm - spmsm)/2.0;
  return;
} //get_ff_isgw_1p1


void EvtISGW2FF::EvtISGW2FF3P1 (EvtId parent,EvtId daugt, 
       double t, double mass, double *lf, double *qf, 
       double *cpf, double *cmf ) {

  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId D3P1P=EvtPDL::getId("D'_1+");
  static EvtId D3P1N=EvtPDL::getId("D'_1-");
  static EvtId D3P10=EvtPDL::getId("D'_10");
  static EvtId D3P1B=EvtPDL::getId("anti-D'_10");

  static EvtId A1P=EvtPDL::getId("a_1+");
  static EvtId A1M=EvtPDL::getId("a_1-");
  static EvtId A10=EvtPDL::getId("a_10");

  static EvtId F1=EvtPDL::getId("f_1");
  static EvtId F1PR=EvtPDL::getId("f'_1");

  static EvtId K1STP=EvtPDL::getId("K'_1+");
  static EvtId K1STM=EvtPDL::getId("K'_1-");
  static EvtId K1ST0=EvtPDL::getId("K'_10");
  static EvtId K1STB=EvtPDL::getId("anti-K'_10");

  static EvtId D3P1SP=EvtPDL::getId("D'_s1+");
  static EvtId D3P1SN=EvtPDL::getId("D'_s1-");

  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BS0=EvtPDL::getId("B_s0");

  double mtb,  mbb(0.0);
  double msd(0.0), mx,mb,nfp(0.0); 
  double msq(0.0),bx2(0.0),mbx(0.0),mtx;
  double f5cppcm,f5cpmcm,f5,ql,ll,cppcm,cpmcm,f5q,f5l;
  double mqm,msb(0.0),bb2(0.0),mum,bbx2,tm,wt,r2;
  EvtId prnt=parent;
  EvtId dgt=daugt;

  if (prnt==B0||prnt==B0B||prnt==BP||prnt==BM) {

    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=5.31;

    if (dgt==A10||dgt==A1P||dgt==A1M||dgt==F1||dgt==F1PR) {

      msq=0.33;
      bx2=0.275*0.275;
      mbx=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
      nfp = 0.0;
    }
    else{
      if (dgt==D3P1P||dgt==D3P1N||dgt==D3P10||dgt==D3P1B) {
	msq=1.82;
	bx2=0.33*0.33;
	mbx=(3.0*2.49+2.40)/4.0;
	nfp = 3.0;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_3P1.\n";
      }
    }
  }
  else{
    if (prnt==D0||prnt==D0B||prnt==DP||prnt==DM) {
    
      msb=1.82;
      msd=0.33;
      bb2=0.45*0.45;
      mbb=1.963;

      if (dgt==F1||dgt==F1PR||dgt==A10||dgt==A1P||dgt==A1M) {

	msq=0.33;
	bx2=0.275*0.275;
	mbx=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
	nfp = 0.0;
      }
      else{
	if (dgt==K1STM||dgt==K1STB||dgt==K1STP||dgt==K1ST0) {
	  msq=0.55;
	  bx2=0.30*0.30;
	  mbx=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
	  nfp = 2.0;
	}
	else{
	  report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_3P1.\n";
	}
      }
    }
    else{
      //BS -> cs constants added by djl on Jan. 21,1998
      if (prnt==BS0||prnt==BSB){

	msb=5.2;
	msd=0.55;
	bb2=0.54*0.54;
	mbb=5.38;
    
	if  (dgt==D3P1SP||dgt==D3P1SN) {
      
	  msq=1.82;
	  bx2=0.41*0.41;
	  mbx=(3.0*2.54+2.46)/4.0;
	  nfp = 3.0;
	}
	else if (dgt==K1STM||dgt==K1STB||dgt==K1STP||dgt==K1ST0) {
	  msq=0.55;
	  bx2=0.30*0.30;
	  mbx=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
	  nfp = 2.0;
	}
	else{
	    report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_1S0.\n";
	  }
	}
    
    else{
      report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_3P1.\n";
    }
    }
  }


  
  mtb = msb + msd;
  mtx = msq + msd;
  
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;
  
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;
  wt=1.0+(tm-t)/(2.0*mbb*mbx);
  
  mqm = 0.1;
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm)/EvtGetas(msq));
  
  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0) /
    (pow((1.0+r2*(tm-t)/18.0),3.0));

  f5q = f5*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),-0.5);
  f5l = f5*pow(( mbb / mtb ),0.5)*pow((mbx/mtx),0.5);
  f5cppcm = f5*pow(( mbb / mtb ),-1.5)*pow((mbx/mtx),0.5);
  f5cpmcm = f5*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),-0.5);
  
  if (msq == msd) { 
    
    ql = -1.0*(msd*(5.0+wt)*f5q/(2.0*mtx*sqrt(bb2)*6.0));
  
    ll = -1.0*mtb*sqrt(bb2)*f5l*(1/mum+ ( (msd*mtx*(wt-1)/bb2)*
         ( (5.0+wt)/(6.0*msq)-(msd*bb2)/(2.0*mum*mtx*bbx2))));
    
    cppcm = (-1.0*(msd*mtx*f5cppcm/(2.0*msq*mtb*sqrt(bb2)))*
            (1-(msd*msq*bb2)/(2.0*mtx*mum*bbx2)));
    
    cpmcm = 1.0*(msd*mtx*f5cpmcm/(2.0*msq*mtb*sqrt(bb2)))*
            (((wt+2.0)/3.0)-(msd*msq*bb2)/(2.0*mtx*mum*bbx2))
            *(mtb/mtx);
  } else {

    ql = f5q*sqrt(1.0/6.0)*msd/(sqrt(bb2)*mtx)*
        (1.0-bb2*mtb/(4.0*msd*msq*msb));
    ll = f5l*sqrt(2.0/3.0)*mtb*sqrt(bb2)*(1.0/(2.0*msq) - 3.0/(2.0*msb) +
         msd*mtx*(wt-1)/bb2*(1.0/msq-msd*bb2/(2.0*mum*mtx*bbx2)));  
    cppcm = msd*msd*bx2*f5cppcm/(sqrt(6.0)*mtb*msq*sqrt(bb2)*bbx2);
    cpmcm = -sqrt(2.0/3.0)*msd*f5cpmcm/(sqrt(bb2)*mtx)*
      (1+msd*bx2/(2.0*msq*bbx2));
  }

  //smooth out the mass(meson) dependence a little
  double parMass=EvtPDL::getMeanMass(prnt);
  double q2max = parMass*parMass + mass*mass - 2.0*parMass*mass;
  double massNom= EvtPDL::getMeanMass(dgt);
  double q2maxNom = parMass*parMass + massNom*massNom - 2.0*parMass*massNom;
  double q2maxin=sqrt(q2maxNom/q2max);
  if ( q2maxin > 1000. ) q2maxin=1000.;
  ql*=q2maxin;
  ll*=q2maxin;
  cppcm*=q2maxin;
  cpmcm*=q2maxin;

  *qf = ql;
  *lf = ll;
  *cpf = (cppcm + cpmcm)/2.0;
  *cmf = (cppcm - cpmcm)/2.0;
  return;
} //get_ff_isgw_3p1


void EvtISGW2FF::EvtISGW2FF3P0 (EvtId parent,EvtId daugt,
       double t, double mass, double *upf, double *umf ) {

  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId D3P0P=EvtPDL::getId("D_0*+");
  static EvtId D3P0N=EvtPDL::getId("D_0*-");
  static EvtId D3P00=EvtPDL::getId("D_0*0");
  static EvtId D3P0B=EvtPDL::getId("anti-D_0*0");

  static EvtId D3P0SP=EvtPDL::getId("D_s0*+");
  static EvtId D3P0SN=EvtPDL::getId("D_s0*-");

  static EvtId A0P=EvtPDL::getId("a_0+");
  static EvtId A0M=EvtPDL::getId("a_0-");
  static EvtId A00=EvtPDL::getId("a_00");

  static EvtId F0=EvtPDL::getId("f_0");
  static EvtId F0PR=EvtPDL::getId("f'_0");

  static EvtId K0STP=EvtPDL::getId("K_0*+");
  static EvtId K0STM=EvtPDL::getId("K_0*-");
  static EvtId K0ST0=EvtPDL::getId("K_0*0");
  static EvtId K0STB=EvtPDL::getId("anti-K_0*0");

  static EvtId DSP=EvtPDL::getId("D_s+");
  static EvtId DSM=EvtPDL::getId("D_s-");

  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BS0=EvtPDL::getId("B_s0");

  double mtb, mbb(0.0);
  double msd(0.0), mx,mb,nfp(0.0); 
  double msq(0.0),bx2(0.0),mbx(0.0),mtx;
  double f5uppum,f5upmum,uppum,upmum,f5;
  double mqm,r2,bb2(0.0),bbx2,msb(0.0),tm;

  EvtId prnt=parent;
  EvtId dgt=daugt;

  if (prnt==B0||prnt==B0B||prnt==BP||prnt==BM) {
      
    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=5.31;
    if (dgt==A00||dgt==A0P||dgt==A0M||dgt==F0||dgt==F0PR) {

      msq=0.33;
      bx2=0.275*0.275;
      mbx=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
      nfp = 0.0;
    }
    else{
      if (dgt==D3P0P||dgt==D3P0N||dgt==D3P00||dgt==D3P0B) {
	msq=1.82;
	bx2=0.33*0.33;
	mbx=(3.0*2.49+2.40)/4.0;
	nfp = 3.0;
      }
      else{
	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_EvtISGW2_ff_3P0.\n";
      }
    }
  }
  else{
    if (prnt==D0||prnt==D0B||prnt==DP||prnt==DM) {

      msb=1.82;
      msd=0.33;
      bb2=0.45*0.45;
      mbb=1.963;
      if (dgt==F0||dgt==F0PR||dgt==A00||dgt==A0P||dgt==A0M) {
	msq=0.33;
	bx2=0.275*0.275;
	mbx=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
	nfp = 0.0;
      }
      else{
	if (dgt==K0STM||dgt==K0STB||dgt==K0STP||dgt==K0ST0) {
	  msq=0.55;
	  bx2=0.30*0.30;
	  mbx=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
	  nfp = 2.0;
	}
	else{
	  report(Severity::Error,"EvtGen") << "Not implemented daugt in get_EvtISGW2_ff_3P0.\n";
	}
      }
    }
    else{
      if (prnt==DSP||prnt==DSM){
	msb=1.82;
	msd=0.55;
	bb2=0.56*0.56;
	mbb=1.968;
	
	if (dgt==F0||dgt==F0PR||dgt==A00||dgt==A0P||dgt==A0M) {
	  msq=0.55;
	  bx2=0.33*0.33;
	  mbx=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
	  nfp = 2.0;
	}
	else{
	  if (dgt==K0STM||dgt==K0STB||dgt==K0STP||dgt==K0ST0) {
	    msq=0.33;
	    bx2=0.30*0.30;
	    mbx=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
	    nfp = 0.0;
	  }
	  else{
	    report(Severity::Error,"EvtGen") << "Not implemented daugt in get_EvtISGW2_ff_3P0.\n";
	  }
	}
      }
      else{
	//BS -> cs constants added by djl on Jan. 21,1998
	if (prnt==BS0||prnt==BSB){

	  msb=5.2;
	  msd=0.55;
	  bb2=0.54*0.54;
	  mbb=5.38;
	  
	  if  (dgt==D3P0SP||dgt==D3P0SN) {
	    
	    msq=1.82;
	    bx2=0.41*0.41;
	    mbx=(3.0*2.54+2.46)/4.0;
	    nfp = 3.0;
	  }
	  else if (dgt==K0STM||dgt==K0STB||dgt==K0STP||dgt==K0ST0) {
	    msq=0.55;
	    bx2=0.30*0.30;
	    mbx=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
	    nfp = 2.0;
	  }
	  else{
	    report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_1S0.\n";
	  }
	}
	else{
	  report(Severity::Error,"EvtGen") << "Not implemented parent in get_EvtISGW2_ff_3P0.\n";
	}
      }
    }
  }
  

  mtb = msb + msd;
  mtx = msq + msd;
  
  mb = EvtPDL::getMeanMass( parent );
  mx = mass;

  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if (t>tm) t = 0.99*tm;
  
  mqm = 0.1;
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
     (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
     log(EvtGetas(mqm)/EvtGetas(msq));

  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0) /
       (pow((1.0+r2*(tm-t)/18.0),3.0));

  f5uppum = f5*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),0.5);
  f5upmum = f5*pow(( mbb / mtb ),0.5)*pow((mbx/mtx),-0.5);

  uppum = -1.0*f5uppum*sqrt(2.0/(3.0*bb2))*msd;
  upmum = 1.0*f5upmum*sqrt(2.0/(3.0*bb2))*msd*mtb/mtx;

  *upf = (uppum + upmum)/2.0;
  *umf = (uppum - upmum)/2.0;

  return;

}


void EvtISGW2FF::EvtISGW2FF3P2 (EvtId parent,EvtId daugt,
       double t, double mass, double *hf, double *kf, 
       double *bpf, double *bmf ) {
  
  //added by Lange Jan4,2000
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");

  static EvtId D3P2P=EvtPDL::getId("D_2*+");
  static EvtId D3P2N=EvtPDL::getId("D_2*-");
  static EvtId D3P20=EvtPDL::getId("D_2*0");
  static EvtId D3P2B=EvtPDL::getId("anti-D_2*0");

  static EvtId A2P=EvtPDL::getId("a_2+");
  static EvtId A2M=EvtPDL::getId("a_2-");
  static EvtId A20=EvtPDL::getId("a_20");

  static EvtId F2=EvtPDL::getId("f_2");
  static EvtId F2PR=EvtPDL::getId("f'_2");

  static EvtId K2STP=EvtPDL::getId("K_2*+");
  static EvtId K2STM=EvtPDL::getId("K_2*-");
  static EvtId K2ST0=EvtPDL::getId("K_2*0");
  static EvtId K2STB=EvtPDL::getId("anti-K_2*0");

  static EvtId D3P2SP=EvtPDL::getId("D_s2*+");
  static EvtId D3P2SN=EvtPDL::getId("D_s2*-");

  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BS0=EvtPDL::getId("B_s0");


  double mtb, mbb(0.0);
  double msd(0.0), mx,mb,nfp(0.0); 
  double msq(0.0),bx2(0.0),mbx(0.0),mtx,f5;
  double f5h,f5k,f5bppbm,f5bpmbm,bppbm,bpmbm;
  double mqm,mum,mup,tm,wt,r2,bb2(0.0),bbx2;
  double msb(0.0);
  EvtId prnt=parent;
  EvtId dgt=daugt;

  if (prnt==B0||prnt==B0B||prnt==BP||prnt==BM) {
    
    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=5.31;

    if (dgt==A20||dgt==A2P||dgt==A2M||dgt==F2||dgt==F2PR) {

      msq=0.33;
      bx2=0.275*0.275;
      mbx=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
      nfp = 0.0;

    }

    else{
      if (dgt==D3P2P||dgt==D3P2N||dgt==D3P20||dgt==D3P2B) {

	msq=1.82;
	bx2=0.33*0.33;
	mbx=(5.0*2.46+3.0*2.42)/8.0;
	nfp = 3.0;
      }
      else{

	report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_3P2.\n";
      }
    }
  }
  else{
    if (prnt==D0||prnt==D0B||prnt==DP||prnt==DM) {
      
      msb=1.82;
      msd=0.33;
      bb2=0.45*0.45;
      mbb=1.963;
      if (dgt==F2||dgt==F2PR||dgt==A20||dgt==A2P||dgt==A2M) {
	msq=0.33;
	bx2=0.275*0.275;
	mbx=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
	nfp = 0.0;
      }
      else{
	if (dgt==K2STM||dgt==K2STB||dgt==K2STP||dgt==K2ST0) {
	  msq=0.55;
	  bx2=0.30*0.30;
	  mbx=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
	  nfp = 2.0;
	}
	else{
	  report(Severity::Error,"EvtGen") << "Not implemented daugt in get_isgw_ff_3P2.\n";
	}
      }
    }
    else{
      //BS -> cs constants added by djl on Jan. 21,1998
      if (prnt==BS0||prnt==BSB){

	msb=5.2;
	msd=0.55;
	bb2=0.54*0.54;
	mbb=5.38;
    
	if  (dgt==D3P2SP||dgt==D3P2SN) {
      
	  msq=1.82;
	  bx2=0.41*0.41;
	  mbx=(5.0*2.61+3.0*2.54)/8.0;
	  nfp = 3.0;
	}
	else if (dgt==K2STM||dgt==K2STB||dgt==K2STP||dgt==K2ST0) {
	  msq=0.55;
	  bx2=0.30*0.30;
	  mbx=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
	  nfp = 2.0;
	}
	else{
	    report(Severity::Error,"EvtGen") << "Not implemented daugt:"<<daugt.getId()<<" in get_isgw_ff_1S0.\n";
	  }
      }
    
    else{
      report(Severity::Error,"EvtGen") << "Not implemented parent in get_isgw_ff_3P2.\n";
    }
    }
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
  wt=1.0+(tm-t)/(2.0*mbb*mbx);
  
  mqm = 0.1;
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2)+
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm)/EvtGetas(msq));

  f5 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,5.0/2.0) /
       (pow((1.0+r2*(tm-t)/18.0),3.0));
  
  f5h = f5*pow(( mbb / mtb ),-1.5)*pow((mbx/mtx),-0.5);
  f5k = f5*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),0.5);
  f5bppbm = f5*pow(( mbb / mtb ),-2.5)*pow((mbx/mtx),0.5);
  f5bpmbm = f5*pow(( mbb / mtb ),-1.5)*pow((mbx/mtx),-0.5);
  
  *hf = f5h*(msd/(sqrt(8.0*bb2)*mtb))*((1.0/msq)-(msd*bb2/(2.0*mum*
        mtx*bbx2)));
  
  *kf = f5k*(msd/(sqrt(2.0*bb2)))*(1.0+wt);
  
  bppbm = ((msd*msd*f5bppbm*bx2)/(sqrt(32.0*bb2)*msq*msb*mtb*bbx2))*
          (1.0-(msd*bx2/(2.0*mtb*bbx2)));

  bpmbm = -1.0*(msd*f5bpmbm/(sqrt(2.0*bb2)*msb*mtx))*(1.0-
          ((msd*msb*bx2)/(2.0*mup*mtb*bbx2))+((msd*bx2*(1.0-
          ((msd*bx2)/(2.0*mtb*bbx2))))/(4.0*msq*bbx2)));

  *bpf = (bppbm + bpmbm)/2.0;
  *bmf = (bppbm - bpmbm)/2.0;
  return;
} //get_ff_isgw_1p1


double EvtISGW2FF::EvtGetGammaji ( double z )

{
double temp;

   temp = 2+((2.0*z)/(1-z))*log(z);
   temp = -1.0*temp;

   return temp;

} //EvtGetGammaji



double EvtISGW2FF::EvtGetas ( double massq, double massx )
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  
  if ( massx > 0.6 ) {
    if ( massq < 1.85 ) {
      nflav = 3.0;}
    
    temp = 12.0*EvtConst::pi / ( 33.0 - 2.0*nflav) /
      log( massx*massx/lqcd2);
  }
  return temp;
  
} //EvtGetas

double EvtISGW2FF::EvtGetas ( double mass )
     
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  
  if ( mass > 0.6 ) {
    if ( mass < 1.85 ) {
      nflav = 3.0;}
    
    temp = 12.0*EvtConst::pi / ( 33.0 - 2.0*nflav) /
      log( mass*mass/lqcd2);
  }
  return temp;
  
} //EvtGetas


void EvtISGW2FF::getbaryonff(EvtId, EvtId, double, double, double*, 
			     double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtISGW2FF.\n";  

  ::abort();

}

void EvtISGW2FF::getdiracff(EvtId, EvtId, double, double, double*, double*,
			    double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtISGW2FF.\n";
  ::abort();

}

void EvtISGW2FF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
			     double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtISGW2FF.\n";
  ::abort();

}
