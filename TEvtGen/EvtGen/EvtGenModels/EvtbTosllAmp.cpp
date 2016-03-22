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
// Module: EvtbTosllAmp.cc
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
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenModels/EvtbTosllAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtDiLog.hh"

double EvtbTosllAmp::CalcMaxProb( EvtId parent, EvtId meson, 
				  EvtId lepton1, EvtId lepton2,
				  EvtbTosllFF *FormFactors,
				  double& poleSize) {

  //This routine takes the arguements parent, meson, and lepton
  //number, and a form factor model, and returns a maximum
  //probability for this semileptonic form factor model.  A
  //brute force method is used.  The 2D cos theta lepton and
  //q2 phase space is probed.

  //Start by declaring a particle at rest.

  //It only makes sense to have a scalar parent.  For now. 
  //This should be generalized later.

  EvtScalarParticle *scalar_part;
  EvtParticle *root_part;

  scalar_part=new EvtScalarParticle;

  //cludge to avoid generating random numbers!
  scalar_part->noLifeTime();

  EvtVector4R p_init;
  
  p_init.set(EvtPDL::getMass(parent),0.0,0.0,0.0);
  scalar_part->init(parent,p_init);
  root_part=(EvtParticle *)scalar_part;
  root_part->setDiagonalSpinDensity();      

  EvtParticle *daughter, *lep1, *lep2;
  
  EvtAmp amp;

  EvtId listdaug[3];
  listdaug[0] = meson;
  listdaug[1] = lepton1;
  listdaug[2] = lepton2;

  amp.init(parent,3,listdaug);

  root_part->makeDaughters(3,listdaug);
  daughter=root_part->getDaug(0);
  lep1=root_part->getDaug(1);
  lep2=root_part->getDaug(2);

  //cludge to avoid generating random numbers!
  daughter->noLifeTime();
  lep1->noLifeTime();
  lep2->noLifeTime();


  //Initial particle is unpolarized, well it is a scalar so it is 
  //trivial
  EvtSpinDensity rho;
  rho.setDiag(root_part->getSpinStates());
  
  double mass[3];
  
  double m = root_part->mass();
  
  EvtVector4R p4meson, p4lepton1, p4lepton2, p4w;
  double q2max;

  double q2, elepton, plepton;
  int i,j;
  double erho,prho,costl;

  double maxfoundprob = 0.0;
  double prob = -10.0;
  int massiter;

  double maxpole=0;

  for (massiter=0;massiter<3;massiter++){

    mass[0] = EvtPDL::getMeanMass(meson);
    mass[1] = EvtPDL::getMeanMass(lepton1);
    mass[2] = EvtPDL::getMeanMass(lepton2);
    if ( massiter==1 ) {
      mass[0] = EvtPDL::getMinMass(meson);
    }
    if ( massiter==2 ) {
      mass[0] = EvtPDL::getMaxMass(meson);
      if ( (mass[0]+mass[1]+mass[2])>m) mass[0]=m-mass[1]-mass[2]-0.00001; 
    }

    q2max = (m-mass[0])*(m-mass[0]);
    
    //loop over q2
    //cout << "m " << m << "mass[0] " << mass[0] << " q2max "<< q2max << endl;
    for (i=0;i<25;i++) {
      //want to avoid picking up the tail of the photon propagator
      q2 = ((i+1.5)*q2max)/26.0;

      if (i==0) q2=4*(mass[1]*mass[1]);

      erho = ( m*m + mass[0]*mass[0] - q2 )/(2.0*m);
      
      prho = sqrt(erho*erho-mass[0]*mass[0]);
      
      p4meson.set(erho,0.0,0.0,-1.0*prho);
      p4w.set(m-erho,0.0,0.0,prho);
      
      //This is in the W rest frame
      elepton = (q2+mass[1]*mass[1])/(2.0*sqrt(q2));
      plepton = sqrt(elepton*elepton-mass[1]*mass[1]);
      
      double probctl[3];

      for (j=0;j<3;j++) {
	
        costl = 0.99*(j - 1.0);
	
	//These are in the W rest frame. Need to boost out into
	//the B frame.
        p4lepton1.set(elepton,0.0,
		  plepton*sqrt(1.0-costl*costl),plepton*costl);
        p4lepton2.set(elepton,0.0,
		 -1.0*plepton*sqrt(1.0-costl*costl),-1.0*plepton*costl);

	EvtVector4R boost((m-erho),0.0,0.0,1.0*prho);
        p4lepton1=boostTo(p4lepton1,boost);
        p4lepton2=boostTo(p4lepton2,boost);

	//Now initialize the daughters...

        daughter->init(meson,p4meson);
        lep1->init(lepton1,p4lepton1);
        lep2->init(lepton2,p4lepton2);

        CalcAmp(root_part,amp,FormFactors);

	//Now find the probability at this q2 and cos theta lepton point
        //and compare to maxfoundprob.

	//Do a little magic to get the probability!!

	//cout <<"amp:"<<amp.getSpinDensity()<<endl;

	prob = rho.normalizedProb(amp.getSpinDensity());

	//cout << "prob:"<<q2<<" "<<costl<<" "<<prob<<endl;

	probctl[j]=prob;
      }

      //probclt contains prob at ctl=-1,0,1.
      //prob=a+b*ctl+c*ctl^2

      double a=probctl[1];
      double b=0.5*(probctl[2]-probctl[0]);
      double c=0.5*(probctl[2]+probctl[0])-probctl[1];

      prob=probctl[0];
      if (probctl[1]>prob) prob=probctl[1];
      if (probctl[2]>prob) prob=probctl[2];

      if (fabs(c)>1e-20){
	double ctlx=-0.5*b/c;
	if (fabs(ctlx)<1.0){
	  double probtmp=a+b*ctlx+c*ctlx*ctlx;
	  if (probtmp>prob) prob=probtmp;
	} 

      }

      //report(Severity::Debug,"EvtGen") << "prob,probctl:"<<prob<<" "
      //			    << probctl[0]<<" "
      //			    << probctl[1]<<" "
      //			    << probctl[2]<<endl;

      if (i==0) {
	maxpole=prob;
	continue;
      }

      if ( prob > maxfoundprob ) {
	maxfoundprob = prob; 
      }

      //cout << "q2,maxfoundprob:"<<q2<<" "<<maxfoundprob<<endl;

    }
    if ( EvtPDL::getWidth(meson) <= 0.0 ) {
      //if the particle is narrow dont bother with changing the mass.
      massiter = 4;
    }

  }

  root_part->deleteTree();  

  poleSize=0.04*(maxpole/maxfoundprob)*4*(mass[1]*mass[1]);

  //poleSize=0.002;

  //cout <<"maxfoundprob,maxpole,poleSize:"<<maxfoundprob<<" "
  //     <<maxpole<<" "<<poleSize<<endl;

  maxfoundprob *=1.15;

  return maxfoundprob;
  
}


EvtComplex EvtbTosllAmp::GetC7Eff(double q2, bool nnlo) 
{

  if (!nnlo) return -0.313;
  double mbeff = 4.8;
  double shat = q2/mbeff/mbeff;
  double logshat;
  logshat = log(shat);
  
  double muscale;
  muscale = 2.5;
  double alphas;
  alphas = 0.267;
  double A7;
  A7 = -0.353 + 0.023;
  double A8;
  A8 = -0.164;
  double C1;
  C1 = -0.697;
  double C2;
  C2 = 1.046;
  
  double Lmu;
  Lmu = log(muscale/mbeff);

  EvtComplex uniti(0.0,1.0);

  EvtComplex c7eff;
  if (shat > 0.25)
  { 
   c7eff = A7;
   return c7eff;
  }

  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  muscale = 5.0;
  alphas = 0.215;
  A7 = -0.312 + 0.008;
  A8 = -0.148;
  C1 = -0.487;
  C2 = 1.024;
  Lmu = log(muscale/mbeff);

  EvtComplex F71;
  EvtComplex f71;
  EvtComplex k7100(-0.68192,-0.074998);
  EvtComplex k7101(0.0,0.0);
  EvtComplex k7110(-0.23935,-0.12289);
  EvtComplex k7111(0.0027424,0.019676);
  EvtComplex k7120(-0.0018555,-0.175);
  EvtComplex k7121(0.022864,0.011456);
  EvtComplex k7130(0.28248,-0.12783);
  EvtComplex k7131(0.029027,-0.0082265);
  f71 = k7100 + k7101*logshat + shat*(k7110 + k7111*logshat) +
        shat*shat*(k7120 + k7121*logshat) + 
        shat*shat*shat*(k7130 + k7131*logshat); 
  F71 = (-208.0/243.0)*Lmu + f71;

  EvtComplex F72;
  EvtComplex f72;
  EvtComplex k7200(4.0915,0.44999);
  EvtComplex k7201(0.0,0.0);
  EvtComplex k7210(1.4361,0.73732);
  EvtComplex k7211(-0.016454,-0.11806);
  EvtComplex k7220(0.011133,1.05);
  EvtComplex k7221(-0.13718,-0.068733);
  EvtComplex k7230(-1.6949,0.76698);
  EvtComplex k7231(-0.17416,0.049359);
  f72 = k7200 + k7201*logshat + shat*(k7210 + k7211*logshat) +
        shat*shat*(k7220 + k7221*logshat) + 
        shat*shat*shat*(k7230 + k7231*logshat); 
  F72 = (416.0/81.0)*Lmu + f72;
  
  EvtComplex F78;
  F78 = (-32.0/9.0)*Lmu + 8.0*EvtConst::pi*EvtConst::pi/27.0 + (-44.0/9.0) 
        + (-8.0*EvtConst::pi/9.0)*uniti +
        (4.0/3.0*EvtConst::pi*EvtConst::pi - 40.0/3.0)*shat +
        (32.0*EvtConst::pi*EvtConst::pi/9.0 - 316.0/9.0)*shat*shat +
        (200.0*EvtConst::pi*EvtConst::pi/27.0 - 658.0/9.0)*shat*shat*shat +
    (-8.0*logshat/9.0)*(shat + shat*shat + shat*shat*shat);
        
  c7eff = A7 - alphas/(4.0*EvtConst::pi)*(C1*F71 + C2*F72 + A8*F78);

  return c7eff;
}


EvtComplex EvtbTosllAmp::GetC9Eff(double q2, bool nnlo, bool btod) 
{

  if (!nnlo) return 4.344;
  double mbeff = 4.8;
  double shat = q2/mbeff/mbeff;
  double logshat;
  logshat = log(shat);
  double mchat = 0.29;

  
  double muscale;
  muscale = 2.5;
  double alphas;
  alphas = 0.267;
  double A8;
  A8 = -0.164;
  double A9;
  A9 = 4.287 + (-0.218);
  double C1;
  C1 = -0.697;
  double C2;
  C2 = 1.046;
  double T9;
  T9 = 0.114 + 0.280;
  double U9;
  U9 = 0.045 + 0.023;
  double W9;
  W9 = 0.044 + 0.016;
  
  double Lmu;
  Lmu = log(muscale/mbeff);


  EvtComplex uniti(0.0,1.0);

  EvtComplex hc;
  double xarg;
  xarg = 4.0*mchat/shat;
  hc = -4.0/9.0*log(mchat*mchat) + 8.0/27.0 + 4.0*xarg/9.0;

if (xarg < 1.0)
  {
    hc = hc - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
      (log(fabs((sqrt(1.0 - xarg)+1.0)/(sqrt(1.0 - xarg) - 1.0))) -
       uniti*EvtConst::pi);
  }
  else
  {
    hc = hc - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
      2.0*atan(1.0/sqrt(xarg - 1.0));
  }
                                                                                                                                                             
  EvtComplex h1;
  xarg = 4.0/shat;
  h1 = 8.0/27.0 + 4.0*xarg/9.0;
  if (xarg < 1.0)
  {
    h1 = h1 - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
      (log(fabs((sqrt(1.0 - xarg)+1.0)/(sqrt(1.0 - xarg) - 1.0))) -
       uniti*EvtConst::pi);
  }
  else
  {
    h1 = h1 - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
      2.0*atan(1.0/sqrt(xarg - 1.0));
  }


  EvtComplex h0;
  h0 = 8.0/27.0 - 4.0*log(2.0)/9.0 + 4.0*uniti*EvtConst::pi/9.0;


  // X=V_{ud}^* V_ub / V_{td}^* V_tb * (4/3 C_1 +C_2) * (h(\hat m_c^2, hat s)-
  // h(\hat m_u^2, hat s))
  EvtComplex Vudstar(1.0 - 0.2279*0.2279/2.0, 0.0);
  EvtComplex Vub((0.118+0.273)/2.0, -1.0*(0.305+0.393)/2.0);
  EvtComplex Vtdstar(1.0 - (0.118+0.273)/2.0,(0.305+0.393)/2.0);
  EvtComplex Vtb(1.0,0.0);

  EvtComplex Xd;
  Xd = (Vudstar * Vub / Vtdstar * Vtb) * (4.0/3.0*C1 + C2) * (hc - h0);


  EvtComplex c9eff=4.344;
  if (shat > 0.25)
  { 
   c9eff =  A9 + T9*hc + U9*h1 + W9*h0;
   if (btod)
   {
    c9eff += Xd; 
   }

   return c9eff;
  }

  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  muscale = 5.0;
  alphas = 0.215;
  A9 = 4.174 + (-0.035);
  C1 = -0.487;
  C2 = 1.024;
  A8 = -0.148;
  T9 = 0.374 + 0.252;
  U9 = 0.033 + 0.015;
  W9 = 0.032 + 0.012;
  Lmu = log(muscale/mbeff);

  EvtComplex F91;
  EvtComplex f91;
  EvtComplex k9100(-11.973,0.16371);
  EvtComplex k9101(-0.081271,-0.059691);
  EvtComplex k9110(-28.432,-0.25044);
  EvtComplex k9111(-0.040243,0.016442);
  EvtComplex k9120(-57.114,-0.86486);
  EvtComplex k9121(-0.035191,0.027909);
  EvtComplex k9130(-128.8,-2.5243);
  EvtComplex k9131(-0.017587,0.050639);
  f91 = k9100 + k9101*logshat + shat*(k9110 + k9111*logshat) +
        shat*shat*(k9120 + k9121*logshat) + 
        shat*shat*shat*(k9130 + k9131*logshat); 
  F91 = (-1424.0/729.0 + 16.0*uniti*EvtConst::pi/243.0 
         + 64.0/27.0*log(mchat))*Lmu - 16.0*Lmu*logshat/243.0 +
        (16.0/1215.0 - 32.0/135.0/mchat/mchat)*Lmu*shat +
        (4.0/2835.0 - 8.0/315.0/mchat/mchat/mchat/mchat)*Lmu*shat*shat +
    (16.0/76545.0 - 32.0/8505.0/mchat/mchat/mchat/mchat/mchat/mchat)*
    Lmu*shat*shat*shat -256.0*Lmu*Lmu/243.0 + f91;

  EvtComplex F92;
  EvtComplex f92;
  EvtComplex k9200(6.6338,-0.98225);
  EvtComplex k9201(0.48763,0.35815);
  EvtComplex k9210(3.3585,1.5026);
  EvtComplex k9211(0.24146,-0.098649);
  EvtComplex k9220(-1.1906,5.1892);
  EvtComplex k9221(0.21115,-0.16745);
  EvtComplex k9230(-17.12,15.146);
  EvtComplex k9231(0.10552,-0.30383);
  f92 = k9200 + k9201*logshat + shat*(k9210 + k9211*logshat) +
        shat*shat*(k9220 + k9221*logshat) + 
        shat*shat*shat*(k9230 + k9231*logshat); 
  F92 = (256.0/243.0 - 32.0*uniti*EvtConst::pi/81.0 
         - 128.0/9.0*log(mchat))*Lmu + 32.0*Lmu*logshat/81.0 +
        (-32.0/405.0 + 64.0/45.0/mchat/mchat)*Lmu*shat +
        (-8.0/945.0 + 16.0/105.0/mchat/mchat/mchat/mchat)*Lmu*shat*shat +
    (-32.0/25515.0 + 64.0/2835.0/mchat/mchat/mchat/mchat/mchat/mchat)*
    Lmu*shat*shat*shat + 512.0*Lmu*Lmu/81.0 + f92;
  
  EvtComplex F98;
  F98 = 104.0/9.0 - 32.0*EvtConst::pi*EvtConst::pi/27.0 + 
        (1184.0/27.0 - 40.0*EvtConst::pi*EvtConst::pi/9.0)*shat +
        (14212.0/135.0 - 32.0*EvtConst::pi*EvtConst::pi/3.0)*shat*shat +
    (193444.0/945.0 - 560.0*EvtConst::pi*EvtConst::pi/27.0)*shat*shat*shat +
        16.0*logshat/9.0*(1.0 + shat + shat*shat + shat*shat*shat);

  Xd = (Vudstar * Vub / Vtdstar * Vtb) * (4.0/3.0*C1 + C2) * (hc - h0);

  c9eff = A9 + T9*hc + U9*h1 + W9*h0 -             
    alphas/(4.0*EvtConst::pi)*(C1*F91 + C2*F92 + A8*F98);
  if (btod)
  {
   c9eff += Xd; 
  }

  return c9eff;
}

EvtComplex EvtbTosllAmp::GetC10Eff(double /*q2*/, bool nnlo) 
{

  if (!nnlo) return -4.669;
  double A10;
  A10 = -4.592 + 0.379;

  EvtComplex c10eff;
  c10eff = A10;

  return c10eff;
}

double EvtbTosllAmp::dGdsProb(double mb, double ms, double ml,
                                double s)
{
  // Compute the decay probability density function given a value of s
  // according to Ali's paper


  double delta, lambda, prob;
  double f1, f2, f3, f4;
  double msh, mlh, sh;

  mlh = ml / mb;
  msh = ms / mb;
  sh  = s  / (mb*mb);

  EvtComplex c9eff = EvtbTosllAmp::GetC9Eff(sh*mb);
  EvtComplex c7eff = EvtbTosllAmp::GetC7Eff(sh*mb);
  EvtComplex c10eff = EvtbTosllAmp::GetC10Eff(sh*mb);

  double alphas = 0.119/
     (1 + 0.119*log(pow(4.8,2)/pow(91.1867,2))*23.0/12.0/EvtConst::pi);
  double omega9 = -2.0/9.0*EvtConst::pi*EvtConst::pi - 4.0/3.0*EvtDiLog::DiLog(sh)
                 - 2.0/3.0*log(sh)*log(1.0-sh)
                 - (5.0+4.0*sh)/(3.0*(1.0+2.0*sh)) * log(1.0-sh)
                 - 2.0*sh*(1.0+sh)*(1.0-2.0*sh)
                 /(3.0*pow(1.0-sh,2)*(1.0+2.0*sh)) * log(sh)
                 + (5.0+9.0*sh-6.0*sh*sh)/(6.0*(1.0-sh)*(1.0+2.0*sh));
  double eta9 = 1.0 + alphas*omega9/EvtConst::pi;
  double omega7 = -8.0/3.0*log(4.8/mb)
                  -4.0/3.0*EvtDiLog::DiLog(sh) 
                  -2.0/9.0*EvtConst::pi*EvtConst::pi
                  -2.0/3.0*log(sh)*log(1.0-sh)
                  -log(1-sh)*(8.0+sh)/(2.0+sh)/3.0 
    -2.0/3.0*sh*(2.0 - 2.0*sh - sh*sh)*log(sh)/pow((1.0 - sh),2)/(2.0 + sh)
    -(16.0 - 11.0*sh - 17.0*sh*sh)/18.0/(2.0 + sh)/(1.0 - sh);
  double eta7 = 1.0 + alphas*omega7/EvtConst::pi;

  double omega79 = -4.0/3.0*log(4.8/mb)
                   -4.0/3.0*EvtDiLog::DiLog(sh) 
                   -2.0/9.0*EvtConst::pi*EvtConst::pi
                   -2.0/3.0*log(sh)*log(1.0-sh)
                   -1.0/9.0*(2.0+7.0*sh)*log(1.0 - sh)/sh
                   -2.0/9.0*sh*(3.0 - 2.0*sh)*log(sh)/pow((1.0 - sh),2) 
                   +1.0/18.0*(5.0 - 9.0*sh)/(1.0 - sh);
  double eta79 = 1.0 + alphas*omega79/EvtConst::pi;

  double c7c9 = abs(c7eff)*real(c9eff);
  c7c9 *= pow(eta79,2); 
  double c7c7 = pow(abs(c7eff),2);
  c7c7 *= pow(eta7,2); 

  double c9c9plusc10c10 = pow(abs(c9eff),2) + pow(abs(c10eff),2);
  c9c9plusc10c10 *= pow(eta9,2);
  double c9c9minusc10c10 = pow(abs(c9eff),2) - pow(abs(c10eff),2);
  c9c9minusc10c10 *= pow(eta9,2);
 
  lambda = 1.0 + sh*sh + pow(msh,4) - 2.0*(sh + sh*msh*msh + msh*msh);

  f1 = pow(1.0-msh*msh,2) - sh*(1.0 + msh*msh);
  f2 = 2.0*(1.0 + msh*msh) * pow(1.0-msh*msh,2)
       - sh*(1.0 + 14.0*msh*msh + pow(msh,4)) - sh*sh*(1.0 + msh*msh);
  f3 = pow(1.0-msh*msh,2) + sh*(1.0 + msh*msh) - 2.0*sh*sh
       + lambda*2.0*mlh*mlh/sh;
  f4 = 1.0 - sh + msh*msh;

  delta = (  12.0*c7c9*f1 + 4.0*c7c7*f2/sh ) * (1.0 + 2.0*mlh*mlh/sh)
            + c9c9plusc10c10*f3 
            + 6.0*mlh*mlh*c9c9minusc10c10*f4;

  prob =  sqrt(lambda*(1.0 - 4.0*mlh*mlh/sh)) * delta;

  return prob;
}

double EvtbTosllAmp::dGdsdupProb(double mb, double ms, double ml,
                                   double s,  double u)
{
  // Compute the decay probability density function given a value of s and u
  // according to Ali's paper

  double prob;
  double f1sp, f2sp, f3sp;

  double sh = s / (mb*mb);

  EvtComplex c9eff = EvtbTosllAmp::GetC9Eff(sh*mb);
  EvtComplex c7eff = EvtbTosllAmp::GetC7Eff(sh*mb);
  EvtComplex c10eff = EvtbTosllAmp::GetC10Eff(sh*mb);

  double alphas = 0.119/
     (1 + 0.119*log(pow(4.8,2)/pow(91.1867,2))*23.0/12.0/EvtConst::pi);
  double omega9 = - 2.0/9.0*EvtConst::pi*EvtConst::pi - 4.0/3.0*EvtDiLog::DiLog(sh)
                 - 2.0/3.0*log(sh)*log(1.0-sh)
                 - (5.0+4.0*sh)/(3.0*(1.0+2.0*sh)) * log(1.0-sh)
                 - 2.0*sh*(1.0+sh)*(1.0-2.0*sh)
                 /(3.0*pow(1.0-sh,2)*(1.0+2.0*sh)) * log(sh)
                 + (5.0+9.0*sh-6.0*sh*sh)/(6.0*(1.0-sh)*(1.0+2.0*sh));
  double eta9 = 1.0 + alphas*omega9/EvtConst::pi;
  double omega7 = -8.0/3.0*log(4.8/mb)
                  -4.0/3.0*EvtDiLog::DiLog(sh) 
                  -2.0/9.0*EvtConst::pi*EvtConst::pi
                  -2.0/3.0*log(sh)*log(1.0-sh)
                  -log(1-sh)*(8.0+sh)/(2.0+sh)/3.0 
    -2.0/3.0*sh*(2.0 - 2.0*sh - sh*sh)*log(sh)/pow((1.0 - sh),2)/(2.0 + sh)
    -(16.0 - 11.0*sh - 17.0*sh*sh)/18.0/(2.0 + sh)/(1.0 - sh);
  double eta7 = 1.0 + alphas*omega7/EvtConst::pi;

  double omega79 = -4.0/3.0*log(4.8/mb)
                   -4.0/3.0*EvtDiLog::DiLog(sh) 
                   -2.0/9.0*EvtConst::pi*EvtConst::pi
                   -2.0/3.0*log(sh)*log(1.0-sh)
                   -1.0/9.0*(2.0+7.0*sh)*log(1.0 - sh)/sh
                   -2.0/9.0*sh*(3.0 - 2.0*sh)*log(sh)/pow((1.0 - sh),2) 
                   +1.0/18.0*(5.0 - 9.0*sh)/(1.0 - sh);
  double eta79 = 1.0 + alphas*omega79/EvtConst::pi;

  double c7c9 = abs(c7eff)*real(c9eff);
  c7c9 *= pow(eta79,2); 
  double c7c7 = pow(abs(c7eff),2);
  c7c7 *= pow(eta7,2); 

  double c9c9plusc10c10 = pow(abs(c9eff),2) + pow(abs(c10eff),2);
  c9c9plusc10c10 *= pow(eta9,2);
  double c9c9minusc10c10 = pow(abs(c9eff),2) - pow(abs(c10eff),2);
  c9c9minusc10c10 *= pow(eta9,2);
  double c7c10 = abs(c7eff)*real(c10eff);
  c7c10 *= eta7; c7c10 *= eta9;
  double c9c10 = real(c9eff)*real(c10eff);
  c9c10 *= pow(eta9,2); 

  f1sp = ( pow(mb*mb-ms*ms,2) - s*s) * c9c9plusc10c10 
         + 4.0*( pow(mb,4) - ms*ms*mb*mb - pow(ms,4)*(1.0 - ms*ms/(mb*mb))
         - 8.0*s*ms*ms - s*s*(1.0 + ms*ms/(mb*mb) ))*mb*mb*c7c7/s
    // kludged mass term
         *(1.0 + 2.0*ml*ml/s)
         - 8.0*(s*(mb*mb + ms*ms) - pow(mb*mb-ms*ms,2)) * c7c9
    // kludged mass term
         *(1.0 + 2.0*ml*ml/s);

  f2sp = 4.0*s*c9c10 + 8.0*(mb*mb + ms*ms)*c7c10;
  f3sp = - (c9c9plusc10c10)
         + 4.0*(1.0 + pow(ms/mb,4)) * mb*mb*c7c7/s
    // kludged mass term
         *(1.0 + 2.0*ml*ml/s);

  prob = (f1sp + f2sp*u + f3sp*u*u)/ pow(mb,3);

  return prob;
}













