//----------------------------------------------------------------------------------
//
// Module: EvtLb2Lll.cpp
//
// Desription: Routine to implement Lambda_b0 -> Lambda_0 l+ l- decays accroding to
//             several models: Chen. Geng.
//                             Aliev. Ozpineci. Savci.
//
// Modification history:
//
//  10/07/2012  MK   Fix calculation of N1, N2; based on hep-ph/021144
//  09/02/2009  PR   Commented check for (anti-)Lambda0 names
//  15/09/2004  PR   Module created according to PHSP model
//  20/02/2005  PR   Added parameters, created matrix element (without polarization)
//  04/03/2005  PR   LD contrib., corrected WC eff. according to Chen. Geng.
//
// Todo list:
//
//   - Properly handle antiparticles, needs change of u, ubar to v, vbar in
//     hadronic current, or other way of putting that in
//----------------------------------------------------------------------------------

#ifdef WIN32
#pragma warning( disable : 4786 )
// Disable anoying warning about symbol size
#endif

#include "EvtGenModels/EvtLb2Lll.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include <stdio.h>
#include <string.h>

EvtLb2Lll::~EvtLb2Lll() {}

EvtDecayBase* EvtLb2Lll::clone(){
  return new EvtLb2Lll;
}

std::string EvtLb2Lll::getName(){
  return "Lb2Lll";
}

void EvtLb2Lll::init(){

  if(getNArg()>8){ // Decay parameters
    report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll generator expected max. 8 arguments but found: " << getNArg() << std::endl;
    report(Severity::Info ,"EvtGen") << "  1. Lambda_b0 polarization - zero is default" << std::endl;
    report(Severity::Info ,"EvtGen") << "  2. Model type - \"SM\" for Standard Model is default" << std::endl;
    report(Severity::Info ,"EvtGen") << "  3. Form-Factors - \"HQET\" is used by default" << std::endl;
    report(Severity::Info ,"EvtGen") << "  4. How to set polarization - \"ModifiedSpinors\" is default" << std::endl;
    report(Severity::Info ,"EvtGen") << "  5. Include long distance (LD) effects - \"SD\" (no) is default" << std::endl;
    report(Severity::Info ,"EvtGen") << "  6. NonFactorizable contribution (omega) to b->sg decay at q2=0 " << std::endl;
    report(Severity::Info ,"EvtGen") << "  7. Note on every x-th decay" << std::endl;
    report(Severity::Info ,"EvtGen") << "  8. Maximum probability - automatic by default" << std::endl;
    ::abort();
  }

  if(getNDaug()!=3){ // Check that there are 3 daughters only
    report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll generator expected 3 daughters but found: " << getNDaug() << std::endl;
    ::abort();
  }

  // TODO: better check based on spin and falvour is needed to allow usage of aliases !
  if(EvtPDL::name(getParentId())=="Lambda_b0"){ // Check daughters of Lambda_b0
    report(Severity::Info,"EvtGen") << " EvtLb2Lll generator found Lambda_b0" << std::endl;
    //if(EvtPDL::name(getDaug(0))!="Lambda0"){
    //  report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll generator expected Lambda0 daughter but found: " << EvtPDL::name(getDaug(0)) << std::endl;
    //  ::abort();
    //}
    if(EvtPDL::name(getDaug(1))=="e-" && EvtPDL::name(getDaug(2))=="e+"){
      m_decayName="Lambda_b0 -> Lambda0 e- e+";
      report(Severity::Info,"EvtGen") << " EvtLb2Lll generator found decay:  Lambda_b0 -> Lambda0 e- e+" << std::endl;
    }else if(EvtPDL::name(getDaug(1))=="mu-" && EvtPDL::name(getDaug(2))=="mu+"){
      m_decayName="Lambda_b0 -> Lambda0 mu- mu+";
      report(Severity::Info,"EvtGen") << " EvtLb2Lll generator found decay:  Lambda_b0 -> Lambda0 mu- mu+" << std::endl;
    }else if(EvtPDL::name(getDaug(1))=="tau-" && EvtPDL::name(getDaug(2))=="tau+"){
      m_decayName="Lambda_b0 -> Lambda0 tau- tau+";
      report(Severity::Info,"EvtGen") << " EvtLb2Lll generator found decay:  Lambda_b0 -> Lambda0 tau- tau+" << std::endl;
    }else{
      report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll generator expected lepton pair daughters but found:  " << EvtPDL::name(getDaug(1)) << " " << EvtPDL::name(getDaug(2)) << std::endl;
      ::abort();
    }
  //TODO: The model is known not to work correctly for anti-Lambda_b0 (A_FB does not change its sign)
  }else if(EvtPDL::name(getParentId())=="anti-Lambda_b0"){ // Check daughters of anti-Lambda_b0
    report(Severity::Info,"EvtGen") << " EvtLb2Lll generator found anti-Lambda_b0" << std::endl;
    //if(EvtPDL::name(getDaug(0))!="anti-Lambda0"){
    //  report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll generator expected anti-Lambda0 daughter but found: " << EvtPDL::name(getDaug(0)) << std::endl;
    //  ::abort();
    //}
    if(EvtPDL::name(getDaug(1))=="e+" && EvtPDL::name(getDaug(2))=="e-"){
      m_decayName="anti-Lambda_b0 -> anti-Lambda0 e+ e-";
      report(Severity::Info,"EvtGen") << " EvtLb2Lll generator found decay:  anti-Lambda_b0 -> anti-Lambda0 e+ e-" << std::endl;
    }else if(EvtPDL::name(getDaug(1))=="mu+" && EvtPDL::name(getDaug(2))=="mu-"){
      m_decayName="anti-Lambda_b0 -> anti-Lambda0 mu+ mu-";
      report(Severity::Info,"EvtGen") << " EvtLb2Lll generator found decay:  anti-Lambda_b0 -> anti-Lambda0 mu+ mu-" << std::endl;
    }else if(EvtPDL::name(getDaug(1))=="tau-" && EvtPDL::name(getDaug(2))=="tau+"){
      m_decayName="anti-Lambda_b0 -> anti-Lambda0 tau+ tau-";
      report(Severity::Info,"EvtGen") << " EvtLb2Lll generator found decay:  anti-Lambda_b0 -> anti-Lambda0 tau+ tau-" << std::endl;
    }else{
      report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll generator expected lepton pair daughters but found:  " << EvtPDL::name(getDaug(1)) << " " << EvtPDL::name(getDaug(2)) << std::endl;
      ::abort();
    }
  }else{ // This model is not intended for decay of anything else than (anti-)Lambda_b0
    report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll generator expected (anti-)Lambda_b0 parent but found:  " << EvtPDL::name(getParentId()) << std::endl;
    ::abort();
  }

  // Read and check all parameters
  if(getNArg()>0){
    if(getArg(0)>1. || getArg(0)<-1.){
      report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll expects polarization to be in interval <-1,1>, not " << getArg(0) << std::endl;
      ::abort();
    }
    m_polarizationLambdab0 = getArg(0);
  }else{
    m_polarizationLambdab0 = 0;
  }
  report(Severity::Info,"EvtGen") << " EvtLb2Lll set Lambda_b0 polarization to " << m_polarizationLambdab0 << std::endl;

  if(getNArg()>1){
    if(getArgStr(1).substr(1,getArgStr(1).size()-2)!="SM" &&
       getArgStr(1).substr(1,getArgStr(1).size()-2)!="-C7_SM" &&
       getArgStr(1).substr(1,getArgStr(1).size()-2)!="SUSY-ChenGeng"){
      report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll doesn't know this physics model: " << getArgStr(1) << std::endl;
      ::abort();
    }
    m_HEPmodel = getArgStr(1).substr(1,getArgStr(1).size()-2);
  }else{
    m_HEPmodel = "SM";
  }
  report(Severity::Info,"EvtGen") << " EvtLb2Lll will use this physics model: " << m_HEPmodel << std::endl;

  if(getNArg()>2){
    if(getArgStr(2).substr(1,getArgStr(2).size()-2)!="HQET" &&
       getArgStr(2).substr(1,getArgStr(2).size()-2)!="HQET-noF2" &&
       getArgStr(2).substr(1,11)                   !="HQET-delta="){
      report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll doesn't know this Form-Factors model: " << getArgStr(2) << std::endl;
      ::abort();
    }
    m_FFtype = getArgStr(2).substr(1,getArgStr(2).size()-2);
  }else{
    m_FFtype = "HQET";
  }
  report(Severity::Info,"EvtGen") << " EvtLb2Lll will use this Form-Factors model: " << m_FFtype << std::endl;

  if(getNArg()>3){
    if(getArgStr(3).substr(1,getArgStr(3).size()-2)!="Unpolarized"){
      report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll doesn't know kind of introducing polarization: " << getArgStr(3) << std::endl;
      ::abort();
    }
    m_polarizationIntroduction = getArgStr(3).substr(1,getArgStr(3).size()-2);
  }else{
    m_polarizationIntroduction = "Unpolarized";
  }
  report(Severity::Info,"EvtGen") << " EvtLb2Lll will use this kind of introducing polarization: " << m_polarizationIntroduction << std::endl;

  if(getNArg()>4){
    if(getArgStr(4).substr(1,getArgStr(4).size()-2)!="SD" && getArgStr(4).substr(1,getArgStr(4).size()-2)!="LD"){
      report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll didn't find SD or LD parameter: " << getArgStr(4) << std::endl;
      ::abort();
    }
    m_effectContribution = getArgStr(5).substr(1,getArgStr(4).size()-2);
  }else{
    m_effectContribution = "SD";
  }
  report(Severity::Info,"EvtGen") << " EvtLb2Lll will include contribution from these effects: " << m_effectContribution << std::endl;

  if(getNArg()>5){
    if(fabs(getArg(5))>0.15){
      report(Severity::Warning,"EvtGen") << " WARNING: EvtLb2Lll found very high contribution to b->sg decay at q2=0: " << getArg(5) << std::endl;
    }
    m_omega = getArg(5);
  }else{
    m_omega = 0;
  }
  report(Severity::Info,"EvtGen") << " EvtLb2Lll will use this contribution to b->sg decay at q2=0: " << m_omega << std::endl;

  if(getNArg()>6) m_noTries = (long)(getArg(6));
  else            m_noTries = 0;

  if(getNArg()>7){
    if(getArg(7)<0.){
      report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll expects positive maximum probability not : " << getArg(7) << std::endl;
      ::abort();
    }
    m_maxProbability = getArg(7);
  }else{
    m_maxProbability = 0.;
  }
  report(Severity::Info,"EvtGen") << " EvtLb2Lll maximum probability was set to " << m_maxProbability << std::endl;
  m_poleSize=0;

  // Initialize Wilson coefficients by Buras and Munz
  // TODO: should have common W.C. source for all decays in EvtGen
  m_WC.CalculateAllCoefficients();

}

void EvtLb2Lll::initProbMax(){

  report(Severity::Info,"EvtGen") << " EvtLb2Lll is finding maximum probability ... " << std::endl;

  if(m_maxProbability==0){

    EvtDiracParticle *parent = new EvtDiracParticle;
    parent->noLifeTime();
    parent->init(getParentId(),EvtVector4R(EvtPDL::getMass(getParentId()),0,0,0));
    parent->setDiagonalSpinDensity();

    EvtAmp amp;
    EvtId daughters[3] = {getDaug(0),getDaug(1),getDaug(2)};
    amp.init(getParentId(),3,daughters);
    parent->makeDaughters(3,daughters);
    EvtParticle *lambda = parent->getDaug(0);
    EvtParticle *lep1   = parent->getDaug(1);
    EvtParticle *lep2   = parent->getDaug(2);
    lambda -> noLifeTime();
    lep1   -> noLifeTime();
    lep2   -> noLifeTime();

    EvtSpinDensity rho;
    rho.setDiag(parent->getSpinStates());

    double M0 = EvtPDL::getMass(getParentId());
    double mL = EvtPDL::getMass(getDaug(0));
    double m1 = EvtPDL::getMass(getDaug(1));
    double m2 = EvtPDL::getMass(getDaug(2));

    double q2,pstar,elambda,theta;
    double q2min = (m1+m2)*(m1+m2);
    double q2max = (M0-mL)*(M0-mL);

    EvtVector4R p4lambda,p4lep1,p4lep2,boost;

    report(Severity::Info,"EvtGen") << " EvtLb2Lll is probing whole phase space ..." << std::endl;

    int i,j;
    double prob=0;
    for(i=0;i<=100;i++){
      q2 = q2min+i*(q2max-q2min)/100.;
      elambda = (M0*M0+mL*mL-q2)/2/M0;
      if(i==0) pstar = 0;
      else     pstar = sqrt(q2-(m1+m2)*(m1+m2))*sqrt(q2-(m1-m2)*(m1-m2))/2/sqrt(q2);
      boost.set(M0-elambda,0,0,+sqrt(elambda*elambda-mL*mL));
      p4lambda.set(elambda,0,0,-sqrt(elambda*elambda-mL*mL));
      for(j=0;j<=45;j++){
        theta = j*EvtConst::pi/45;
        p4lep1.set(sqrt(pstar*pstar+m1*m1),0,+pstar*sin(theta),+pstar*cos(theta));
        p4lep2.set(sqrt(pstar*pstar+m2*m2),0,-pstar*sin(theta),-pstar*cos(theta));
	//std::cout << "p1: " << p4lep1 << " p2: " << p4lep2 << " pstar: " << pstar << std::endl;
	p4lep1 = boostTo(p4lep1,boost);
	p4lep2 = boostTo(p4lep1,boost);
	lambda -> init(getDaug(0),p4lambda);
	lep1   -> init(getDaug(1),p4lep1  );
	lep2   -> init(getDaug(2),p4lep2  );
	calcAmp(&amp,parent);
	prob = rho.normalizedProb(amp.getSpinDensity());
	//std::cout << "q2:  " << q2 << " \t theta:  " << theta << " \t prob:  " << prob << std::endl;
	//std::cout << "p1: " << p4lep1 << " p2: " << p4lep2 << " q2-q2min: " << q2-(m1+m2)*(m1+m2) << std::endl;
	if(prob>m_maxProbability){
	  report(Severity::Info,"EvtGen") << "  - probability " << prob << " found at q2 = " << q2 << " (" << 100*(q2-q2min)/(q2max-q2min) << " %) and theta = " << theta*180/EvtConst::pi << std::endl;
	  m_maxProbability=prob;
	}
      }
      //::abort();
    }

    //m_poleSize = 0.04*q2min;
    m_maxProbability *= 1.2;
    delete parent;
  }

  setProbMax(m_maxProbability);
  report(Severity::Info,"EvtGen") << " EvtLb2Lll set up maximum probability to " << m_maxProbability << std::endl;

}

void EvtLb2Lll::decay(EvtParticle* parent){

  //setWeight(parent->initializePhaseSpace(getNDaug(),getDaugs(),m_poleSize,1,2));
  parent->initializePhaseSpace(getNDaug(),getDaugs());
  calcAmp(&_amp2,parent);
  
}

void EvtLb2Lll::calcAmp(EvtAmp *amp,EvtParticle *parent){

  static long noTries=0;
  static double delta=0;

  EvtComplex Matrix[2][2][2][2];

  EvtComplex i1(0,1);

  int i,j,spins[4];
  char ch;

  double r,M_L,M_Lb,M_s,M_c,M_b,q2,alpha,M_W,M_t;
  double M_psi[2]={0,0},Gamma_psi[2]={0,0},k_psi[2]={0,0};
  double F0_1,F0_2,a_F1,a_F2,b_F1,b_F2,F1,F2;
  double f_1,f_2,f_3,g_1,g_2,g_3,f_1T,f_2T,f_3T,g_1T,g_2T,g_3T,f_TV,f_TS,g_TV(0.0),g_TS,f_T,g_T;
  EvtComplex A1,A2,A3,B1,B2,B3,D1,D2,D3,E1,E2,E3,N1,N2,H1,H2;
  EvtComplex C_SL,C_BR,C_LLtot,C_LRtot,C_LL,C_LR,C_RL,C_RR,C_LRLR,C_RLLR,C_LRRL,C_RLRL,C_T,C_TE;
  EvtComplex Yld,C_7eff,C_9eff;
  EvtComplex V_ts,V_tb;

  EvtVector4C lbar_Gmu_l[2][2],lbar_GmuG5_l[2][2],hbar_GmuPlusG5_h[2][2],hbar_GmuMinusG5_h[2][2],hbar_Gmu_h[2][2];
  EvtComplex  lbar_l[2][2],lbar_G5_l[2][2],hbar_1PlusG5_h[2][2],hbar_1MinusG5_h[2][2],hbar_G5_h[2][2],hbar_h[2][2];
  EvtTensor4C lbar_Smunu_l[2][2],lbar_ESmunu_l[2][2],hbar_SmunuPlusG5_h[2][2],hbar_SmunuMinusG5_h[2][2],hbar_Smunu_h[2][2];
  EvtVector4R q_mu,P_mu;

  EvtDiracSpinor parent__spParent[2];

  M_Lb    = parent->mass();
  M_L     = parent->getDaug(0)->mass();
  M_s     = 0.13;
  M_c     = 1.35;
  M_b     = 4.8;
  alpha   = 1./137.036;
  M_W     = 80.425;
  M_t     = 174.3;
  M_psi[0] = 3.096916;
  M_psi[1] = 3.686093;
  if(m_decayName=="Lambda_b0 -> Lambda0 e- e+" || m_decayName=="anti-Lambda_b0 -> anti-Lambda0 e+ e-"){
    Gamma_psi[0] = 5.40;
    Gamma_psi[1] = 2.12;
  }
  if(m_decayName=="Lambda_b0 -> Lambda0 mu- mu+" || m_decayName=="anti-Lambda_b0 -> anti-Lambda0 mu+ mu-"){
    Gamma_psi[0] = 5.35;
    Gamma_psi[1] = 2.05;
  }
  if(m_decayName=="Lambda_b0 -> Lambda0 tau- tau+" || m_decayName=="anti-Lambda_b0 -> anti-Lambda0 tau+ tau-"){
    Gamma_psi[0] = 0.00;
    Gamma_psi[1] = 0.79;
  }
  if(m_effectContribution=="LD"){
    k_psi[0] = 1.65;
    k_psi[1] = 1.65;
  }
  //G_F   = 1.16637e-5;
  //V_tb  = sqrt(1-pow(0.0413,2))*sqrt(1-pow(0.0037,2));
  //V_ts  = -sqrt(1-pow(0.2243,2))*0.0413-0.2243*sqrt(1-pow(0.0413,2))*0.0037*(cos(60*EvtConst::pi/180)+i1*sin(60*EvtConst::pi/180));

  P_mu = parent->getP4Restframe()+parent->getDaug(0)->getP4();
  q_mu = parent->getP4Restframe()-parent->getDaug(0)->getP4();
  q2   = q_mu.mass2();

  if(m_noTries>0) if(!((++noTries)%m_noTries)) report(Severity::Debug,"EvtGen") << " EvtLb2Lll already finished " << noTries << " matrix element calculations" << std::endl;

  if(m_FFtype=="HQET"){
    r = M_L*M_L/M_Lb/M_Lb;
    F0_1 = +0.462;
    F0_2 = -0.077;
    a_F1 = -0.0182;
    a_F2 = -0.0685;
    b_F1 = -0.000176;
    b_F2 = +0.001460;
    F1 = F0_1/(1.0-(q2/M_Lb/M_Lb)*(a_F1-b_F1*(q2/M_Lb/M_Lb)));
    F2 = F0_2/(1.0-(q2/M_Lb/M_Lb)*(a_F2-b_F2*(q2/M_Lb/M_Lb)));
    g_1 = f_1 = f_2T = g_2T = F1+sqrt(r)*F2;
    //std::cout << " F1: " << F1 << "  F2: " << F2 << "  r: " << r << "  M_L: " << M_L << "  M_Lb: " << M_Lb << std::endl;
    //std::cout << " sqrt(q2): " << sqrt(q2) << "  q2: " << q2 << "  M_Lb^2" << M_Lb*M_Lb << std::endl;
    g_2 = f_2 = g_3 = f_3 = g_TV = f_TV = F2/M_Lb;
    g_TS = f_TS = 0;
    g_1T = f_1T = F2/M_Lb*q2;
    g_3T = +F2/M_Lb*(M_Lb+M_L);
    f_3T = -F2/M_Lb*(M_Lb-M_L);
    f_T = f_2T-f_TS*q2;
    g_T = g_2T-g_TS*q2;
  }else if(strstr(m_FFtype.c_str(),"HQET-delta=")==m_FFtype.c_str()){
    //report(Severity::Warning,"EvtGen") << " WARNING: HQET-delta FF model should be checked for correctness" << std::endl;
    if(delta==0) sscanf(m_FFtype.c_str(),"%c%c%c%c%c%c%c%c%c%c%c%lf",&ch,&ch,&ch,&ch,&ch,&ch,&ch,&ch,&ch,&ch,&ch,&delta);
    r = M_L*M_L/M_Lb/M_Lb;
    F0_1 = +0.462;
    F0_2 = -0.077;
    a_F1 = -0.0182;
    a_F2 = -0.0685;
    b_F1 = -0.000176;
    b_F2 = +0.001460;
    F1 = F0_1/(1.0-(q2/M_Lb/M_Lb)*(a_F1-b_F1*(q2/M_Lb/M_Lb)));
    F2 = F0_2/(1.0-(q2/M_Lb/M_Lb)*(a_F2-b_F2*(q2/M_Lb/M_Lb)));
    g_1 = f_1 = f_2T = g_2T = F1+sqrt(r)*F2;
    g_1 += delta*g_1;
    f_1 -= delta*f_1;
    g_2 = f_2 = g_3 = f_3 = g_TV = f_TV = F2/M_Lb;
    g_TS = f_TS = 0;
    g_1T = f_1T = F2/M_Lb*q2;
    g_3T = +F2/M_Lb*(M_Lb+M_L);
    f_3T = -F2/M_Lb*(M_Lb-M_L);
    f_T = f_2T-f_TS*q2;
    g_T = g_2T-g_TS*q2;
  }else if(m_FFtype=="HQET-noF2"){
    //report(Severity::Warning,"EvtGen") << " WARNING: HQET-noF2 FF model should be checked for correctness" << std::endl;
    r = M_L*M_L/M_Lb/M_Lb;
    F0_1 = +0.462;
    a_F1 = -0.0182;
    b_F1 = -0.000176;
    F1 = F0_1/(1.0-(q2/M_Lb/M_Lb)*(a_F1-b_F1*(q2/M_Lb/M_Lb)));
    g_1 = f_1 = f_2T = g_2T = F1;
    g_2 = f_2 = g_3 = f_3 = g_TV = f_TV = 0;
    g_TS = f_TS = 0;
    g_1T = f_1T = 0;
    g_3T = 0;
    f_3T = 0;
    f_T = f_2T-f_TS*q2;
    g_T = g_2T-g_TS*q2;
  }else{ // general relations for Form-Factors
    f_1 = f_2 = f_3 = g_1 = g_2 = g_3 = f_3T = g_3T = f_TS = g_TS = f_T = g_T = f_TV = 0;
    f_2T = f_T+f_TS*q2;
    f_1T = (f_TV+f_TS*(M_L+M_Lb))*q2;
    f_1T = -q2/(M_Lb-M_L)*f_3T;
    g_2T = g_T+g_TS*q2;
    g_1T = (g_TV-g_TS*(M_L-M_Lb))*q2;
    g_1T = +q2/(M_Lb+M_L)*g_3T;
    report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll - unknown Form-Factors model: " << m_FFtype << " - this should never happen !" << std::endl;
    ::abort();
  }

  if(m_HEPmodel=="SM"){
    C_LL=C_LR=C_RL=C_RR=C_LRLR=C_RLLR=C_LRRL=C_RLRL=C_T=C_TE=EvtComplex(0,0);
    Yld = m_WC.Yld(q2,k_psi,Gamma_psi,M_psi,2,m_WC.GetC1(),m_WC.GetC2(),m_WC.GetC3(),m_WC.GetC4(),m_WC.GetC5(),m_WC.GetC6(),1./alpha);
    C_7eff = m_WC.GetC7eff0() + m_WC.C7b2sg(m_WC.GetStrongCouplingConst(),m_WC.GetEta(),m_WC.GetC2(),M_t,M_W) + m_omega*(m_WC.hzs(M_c/M_b,q2/M_b/M_b,M_b,M_b) + Yld);
    C_9eff = Yld + m_WC.C9efftilda(M_c/M_b,q2/M_b/M_b,m_WC.GetStrongCouplingConst(),m_WC.GetC1(),m_WC.GetC2(),m_WC.GetC3(),m_WC.GetC4(),m_WC.GetC5(),m_WC.GetC6(),m_WC.GetC9tilda(),m_WC.GetRenormSchemePar());
    C_SL = -2*M_s*C_7eff;
    C_BR = -2*M_b*C_7eff;
    C_LLtot = C_9eff-m_WC.GetC10tilda()+C_LL;
    C_LRtot = C_9eff+m_WC.GetC10tilda()+C_LR;
    //std::cout << "Yld: " << Yld << "  C7eff: " << C_7eff << "  C_9eff: " << C_9eff << "  Diff7: " << C_7eff-m_WC.GetC7eff0() << "  Diff9: " << C_9eff-m_WC.GetC9tilda() << std::endl;
  }else if(m_HEPmodel=="-C7_SM"){
    C_LL=C_LR=C_RL=C_RR=C_LRLR=C_RLLR=C_LRRL=C_RLRL=C_T=C_TE=EvtComplex(0,0);
    Yld = m_WC.Yld(q2,k_psi,Gamma_psi,M_psi,2,m_WC.GetC1(),m_WC.GetC2(),m_WC.GetC3(),m_WC.GetC4(),m_WC.GetC5(),m_WC.GetC6(),1./alpha);
    C_7eff = m_WC.GetC7eff0() + m_WC.C7b2sg(m_WC.GetStrongCouplingConst(),m_WC.GetEta(),m_WC.GetC2(),M_t,M_W) + m_omega*(m_WC.hzs(M_c/M_b,q2/M_b/M_b,M_b,M_b) + Yld);
    C_9eff = Yld + m_WC.C9efftilda(M_c/M_b,q2/M_b/M_b,m_WC.GetStrongCouplingConst(),m_WC.GetC1(),m_WC.GetC2(),m_WC.GetC3(),m_WC.GetC4(),m_WC.GetC5(),m_WC.GetC6(),m_WC.GetC9tilda(),m_WC.GetRenormSchemePar());
    C_SL = +2*M_s*C_7eff;
    C_BR = +2*M_b*C_7eff;
    C_LLtot = C_9eff-m_WC.GetC10tilda()+C_LL;
    C_LRtot = C_9eff+m_WC.GetC10tilda()+C_LR;
    //std::cout << "Yld: " << Yld << "  C7eff: " << C_7eff << "  C_9eff: " << C_9eff << "  Diff7: " << C_7eff-m_WC.GetC7eff0() << "  Diff9: " << C_9eff-m_WC.GetC9tilda() << std::endl;
  }else if(m_HEPmodel=="SUSY-ChenGeng"){
    //report(Severity::Warning,"EvtGen") << " WARNING: SUSY-ChenGeng model should be checked for correctness" << std::endl;
    C_LL=C_LR=C_RL=C_RR=C_LRLR=C_RLLR=C_LRRL=C_RLRL=C_T=C_TE=EvtComplex(0,0);
    EvtComplex d_u23LL = 0.1;
    EvtComplex d_u33RL = 0.65;
    EvtComplex d_d23LR = 0.03*exp(i1*EvtConst::pi*2/5);
    EvtComplex d_u23LR = -0.8*exp(i1*EvtConst::pi/4);
    EvtComplex C_7susy = -1.75*d_u23LL - 0.25*d_u23LR - 10.3*d_d23LR;
    EvtComplex C_9susy = 0.82*d_u23LR;
    EvtComplex C_10susy = -9.37*d_u23LR + 1.4*d_u23LR*d_u33RL + 2.7*d_u23LL;
    Yld = m_WC.Yld(q2,k_psi,Gamma_psi,M_psi,2,m_WC.GetC1(),m_WC.GetC2(),m_WC.GetC3(),m_WC.GetC4(),m_WC.GetC5(),m_WC.GetC6(),1./alpha);
    C_7eff = m_WC.GetC7eff0() + C_7susy*pow(m_WC.GetEta(),16./23.) + m_WC.C7b2sg(m_WC.GetStrongCouplingConst(),m_WC.GetEta(),m_WC.GetC2(),M_t,M_W) + m_omega*(m_WC.hzs(M_c/M_b,q2/M_b/M_b,M_b,M_b) + Yld);
    C_9eff = Yld + m_WC.C9efftilda(M_c/M_b,q2/M_b/M_b,m_WC.GetStrongCouplingConst(),m_WC.GetC1(),m_WC.GetC2(),m_WC.GetC3(),m_WC.GetC4(),m_WC.GetC5(),m_WC.GetC6(),m_WC.GetC9tilda()+C_9susy,m_WC.GetRenormSchemePar());
    C_SL = -2*M_s*C_7eff;
    C_BR = -2*M_b*C_7eff;
    C_LLtot = C_9eff-m_WC.GetC10tilda()-C_10susy+C_LL;
    C_LRtot = C_9eff+m_WC.GetC10tilda()+C_10susy+C_LR;
    //std::cout << "Yld: " << Yld << "  C7eff: " << C_7eff << "  C_9eff: " << C_9eff << "  Diff7: " << C_7eff-m_WC.GetC7eff0() << "  Diff9: " << C_9eff-m_WC.GetC9tilda() << std::endl;
  }else{
    report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll - unknown physics model: " << m_HEPmodel << " - this should never happen !" << std::endl;
    ::abort();
  }

  A1 = (f_1T-g_1T)*C_SL/q2 + (f_1T+g_1T)*C_BR/q2 + 0.5*(f_1-g_1)*(C_LLtot+C_LRtot) + 0.5*(f_1+g_1)*(C_RL+C_RR);
  //std::cout << "f_1T: " << f_1T << "  g_1T: " << g_1T << "  C_SL: " << C_SL << "  C_BR: " << C_BR << "  f_1: " << f_1 << "  g_1: " << g_1 << "  C_LLtot: " << C_LLtot << "  C_LRtot: " << C_LRtot << "  C_RL: " << C_RL << "  C_RR: " << C_RR << std::endl;
  A2 = (f_2T-g_2T)*C_SL/q2 + (f_2T+g_2T)*C_BR/q2 + 0.5*(f_2-g_2)*(C_LLtot+C_LRtot) + 0.5*(f_2+g_2)*(C_RL+C_RR);
  A3 = (f_3T-g_3T)*C_SL/q2 + (f_3T+g_3T)*C_BR/q2 + 0.5*(f_3-g_3)*(C_LLtot+C_LRtot) + 0.5*(f_3+g_3)*(C_RL+C_RR);

  B1 = (f_1T+g_1T)*C_SL/q2 + (f_1T-g_1T)*C_BR/q2 + 0.5*(f_1+g_1)*(C_LLtot+C_LRtot) + 0.5*(f_1-g_1)*(C_RL+C_RR);
  B2 = (f_2T+g_2T)*C_SL/q2 + (f_2T-g_2T)*C_BR/q2 + 0.5*(f_2+g_2)*(C_LLtot+C_LRtot) + 0.5*(f_2-g_2)*(C_RL+C_RR);
  B3 = (f_3T+g_3T)*C_SL/q2 + (f_3T-g_3T)*C_BR/q2 + 0.5*(f_3+g_3)*(C_LLtot+C_LRtot) + 0.5*(f_3-g_3)*(C_RL+C_RR);

  D1 = 0.5*(C_RR-C_RL)*(f_1+g_1) + 0.5*(C_LRtot-C_LLtot)*(f_1-g_1);
  D2 = 0.5*(C_RR-C_RL)*(f_2+g_2) + 0.5*(C_LRtot-C_LLtot)*(f_2-g_2);
  D3 = 0.5*(C_RR-C_RL)*(f_3+g_3) + 0.5*(C_LRtot-C_LLtot)*(f_3-g_3);

  E1 = 0.5*(C_RR-C_RL)*(f_1-g_1) + 0.5*(C_LRtot-C_LLtot)*(f_1+g_1);
  E2 = 0.5*(C_RR-C_RL)*(f_2-g_2) + 0.5*(C_LRtot-C_LLtot)*(f_2+g_2);
  E3 = 0.5*(C_RR-C_RL)*(f_3-g_3) + 0.5*(C_LRtot-C_LLtot)*(f_3+g_3);

  N1 = (f_1*(M_Lb-M_L)+f_3*q2)/M_b*(C_LRLR+C_RLLR+C_LRRL+C_RLRL);  // Should be mLb - mL
  N2 = (f_1*(M_Lb-M_L)+f_3*q2)/M_b*(C_LRLR+C_RLLR-C_LRRL-C_RLRL);

  H1 = (g_1*(M_Lb+M_L)-g_3*q2)/M_b*(C_LRLR-C_RLLR+C_LRRL-C_RLRL);
  H2 = (g_1*(M_Lb+M_L)-g_3*q2)/M_b*(C_LRLR-C_RLLR-C_LRRL+C_RLRL);


  for(i=0;i<4;i++){
    lbar_Gmu_l   [i/2][i%2] = EvtLeptonVCurrent(parent->getDaug(1)->spParent(i/2),parent->getDaug(2)->spParent(i%2));
    lbar_GmuG5_l [i/2][i%2] = EvtLeptonACurrent(parent->getDaug(1)->spParent(i/2),parent->getDaug(2)->spParent(i%2));
    lbar_l       [i/2][i%2] = EvtLeptonSCurrent(parent->getDaug(1)->spParent(i/2),parent->getDaug(2)->spParent(i%2));
    lbar_G5_l    [i/2][i%2] = EvtLeptonPCurrent(parent->getDaug(1)->spParent(i/2),parent->getDaug(2)->spParent(i%2));
    lbar_Smunu_l [i/2][i%2] = EvtLeptonTCurrent(parent->getDaug(1)->spParent(i/2),parent->getDaug(2)->spParent(i%2));
    lbar_ESmunu_l[i/2][i%2] = dual(EvtLeptonTCurrent(parent->getDaug(1)->spParent(i/2),parent->getDaug(2)->spParent(i%2)));
  }

  // TODO: polarization not yet introduced
  if(m_polarizationIntroduction=="SpinDensityMatrix"){
    //parent->setSpinDensityForward();
    parent__spParent[0]=parent->sp(0);
    parent__spParent[1]=parent->sp(1);
  }else if(m_polarizationIntroduction=="ModifiedSpinors"){
    parent__spParent[0]=parent->sp(0);
    parent__spParent[1]=parent->sp(1);
  }else if(m_polarizationIntroduction=="Unpolarized"){
    parent__spParent[0]=parent->sp(0);
    parent__spParent[1]=parent->sp(1);
  }else{
    report(Severity::Error,"EvtGen") << " ERROR: EvtLb2Lll - unknown polarization: " << m_polarizationIntroduction << " - this should never happen !" << std::endl;
    ::abort();
  }

  for(i=0;i<4;i++){
    hbar_GmuPlusG5_h   [i/2][i%2] = EvtLeptonVCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]) + EvtLeptonACurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_GmuMinusG5_h  [i/2][i%2] = EvtLeptonVACurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_SmunuPlusG5_h [i/2][i%2] = EvtLeptonTCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]) + EvtLeptonTG5Current(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_SmunuMinusG5_h[i/2][i%2] = EvtLeptonTCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]) - EvtLeptonTG5Current(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_1PlusG5_h     [i/2][i%2] = EvtLeptonSCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]) + EvtLeptonPCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_1MinusG5_h    [i/2][i%2] = EvtLeptonSCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]) - EvtLeptonPCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_G5_h          [i/2][i%2] = EvtLeptonPCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_h             [i/2][i%2] = EvtLeptonSCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_Smunu_h       [i/2][i%2] = EvtLeptonTCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
    hbar_Gmu_h         [i/2][i%2] = EvtLeptonVCurrent(parent->getDaug(0)->spParent(i/2),parent__spParent[i%2]);
  }

  for(i=0;i<4;i++) for(j=0;j<4;j++) {
    //std::cout << "--------------------------------------------------" << std::endl;
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += lbar_Gmu_l[i/2][i%2]    * (A1*hbar_GmuPlusG5_h[j/2][j%2]+B1*hbar_GmuMinusG5_h[j/2][j%2]);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    //std::cout << "A1: " << A1 << " B1: " << B1 << " lbar_Gmu_l: " << lbar_Gmu_l[i/2][i%2] << 
    //             " hbar_GmuPlusG5_h: " << hbar_GmuPlusG5_h[j/2][j%2] << " hbar_GmuMinusG5_h: " << hbar_GmuMinusG5_h[j/2][j%2] <<
    //	           " sp1: " << parent->getDaug(1)->spParent(i/2) << " sp2: " << parent->getDaug(1)->spParent(i%2) << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += lbar_Gmu_l[i/2][i%2]    * (i1*A2*(hbar_SmunuPlusG5_h[j/2][j%2].cont2(q_mu))+B2*(hbar_SmunuMinusG5_h[j/2][j%2].cont2(q_mu)));
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += lbar_Gmu_l[i/2][i%2]    * ((A3*hbar_1PlusG5_h[j/2][j%2]+B3*hbar_1MinusG5_h[j/2][j%2])*q_mu);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += lbar_GmuG5_l[i/2][i%2]  * (D1*hbar_GmuPlusG5_h[j/2][j%2]+E1*hbar_GmuMinusG5_h[j/2][j%2]);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += lbar_GmuG5_l[i/2][i%2]  * (i1*D2*(hbar_SmunuPlusG5_h[j/2][j%2].cont2(q_mu))+E2*(hbar_SmunuMinusG5_h[j/2][j%2].cont2(q_mu)));
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += lbar_GmuG5_l[i/2][i%2]  * ((D3*hbar_1PlusG5_h[j/2][j%2]+E3*hbar_1MinusG5_h[j/2][j%2])*q_mu);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += lbar_l[i/2][i%2]        * (N1*hbar_h[j/2][j%2]+H1*hbar_G5_h[j/2][j%2]);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += lbar_G5_l[i/2][i%2]     * (N2*hbar_h[j/2][j%2]+H2*hbar_G5_h[j/2][j%2]);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += cont(lbar_Smunu_l[i/2][i%2]  , 4*C_T*f_T*hbar_Smunu_h[j/2][j%2]);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += 
      cont(lbar_Smunu_l[i/2][i%2]  , 
           -4*C_T*f_TV*i1*(EvtGenFunctions::directProd(q_mu,hbar_Gmu_h[j/2][j%2])-
                           EvtGenFunctions::directProd(hbar_Gmu_h[j/2][j%2],q_mu)));
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += cont(lbar_Smunu_l[i/2][i%2]  , -4*C_T*f_TS*i1*(EvtGenFunctions::directProd(P_mu,q_mu)-EvtGenFunctions::directProd(q_mu,P_mu))*hbar_h[j/2][j%2]);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += cont(lbar_ESmunu_l[i/2][i%2] , 4*C_TE*f_T*i1*hbar_Smunu_h[j/2][j%2]);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += cont(lbar_ESmunu_l[i/2][i%2] , 4*C_TE*f_TV*(EvtGenFunctions::directProd(q_mu,hbar_Gmu_h[j/2][j%2])-EvtGenFunctions::directProd(hbar_Gmu_h[j/2][j%2],q_mu)));
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    Matrix[j/2][j%2][i/2][i%2] += cont(lbar_ESmunu_l[i/2][i%2] , 4*C_TE*f_TS*(EvtGenFunctions::directProd(P_mu,q_mu)-EvtGenFunctions::directProd(q_mu,P_mu))*hbar_h[j/2][j%2]);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    //Matrix[j/2][j%2][i/2][i%2] *= G_F*alpha/4/sqrt(2)/EvtConst::pi*V_tb*conj(V_ts);
    //std::cout << "Matrix = " << Matrix[j/2][j%2][i/2][i%2] << std::endl;
    //std::cout << "--------------------------------------------------" << std::endl;
    spins[0]=j%2;
    spins[1]=j/2;
    spins[2]=i/2;
    spins[3]=i%2;
    amp->vertex(spins,Matrix[j/2][j%2][i/2][i%2]);
  }

  //std::cout << "==================================================" << std::endl;
  //std::cout << "Lambda_b0:  " << parent->getP4Restframe() << std::endl;
  //std::cout << "Lambda0:  " << parent->getDaug(0)->getP4() << std::endl;
  //std::cout << "mu-:  " << parent->getDaug(1)->getP4() << std::endl;
  //std::cout << "mu+:  " << parent->getDaug(2)->getP4() << std::endl;
  //std::cout << "P_mu:  " << P_mu << std::endl;
  //std::cout << "q_mu:  " << q_mu << std::endl;
  //std::cout << "q2:  " << q2 << std::endl;
  //std::cout << "==================================================" << std::endl;

  return;
}

EvtTensor4C EvtLb2Lll::EvtLeptonTG5Current(const EvtDiracSpinor &d,const EvtDiracSpinor &dp){
// <u|sigma^munu*gamma^5|v>

  EvtTensor4C temp;
  temp.zero();
  EvtComplex i2(0,0.5);

  static EvtGammaMatrix mat01=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g0()*EvtGammaMatrix::g1()-
     EvtGammaMatrix::g1()*EvtGammaMatrix::g0())*EvtGammaMatrix::g5();
  static EvtGammaMatrix mat02=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g0()*EvtGammaMatrix::g2()-
     EvtGammaMatrix::g2()*EvtGammaMatrix::g0())*EvtGammaMatrix::g5();
  static EvtGammaMatrix mat03=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g0()*EvtGammaMatrix::g3()-
     EvtGammaMatrix::g3()*EvtGammaMatrix::g0())*EvtGammaMatrix::g5();
  static EvtGammaMatrix mat12=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g1()*EvtGammaMatrix::g2()-
     EvtGammaMatrix::g2()*EvtGammaMatrix::g1())*EvtGammaMatrix::g5();
  static EvtGammaMatrix mat13=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g1()*EvtGammaMatrix::g3()-
     EvtGammaMatrix::g3()*EvtGammaMatrix::g1())*EvtGammaMatrix::g5();
  static EvtGammaMatrix mat23=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g2()*EvtGammaMatrix::g3()-
     EvtGammaMatrix::g3()*EvtGammaMatrix::g2())*EvtGammaMatrix::g5();

  temp.set(0,1,i2*(d*(mat01*dp)));
  temp.set(1,0,-temp.get(0,1));

  temp.set(0,2,i2*(d*(mat02*dp)));
  temp.set(2,0,-temp.get(0,2));

  temp.set(0,3,i2*(d*(mat03*dp)));
  temp.set(3,0,-temp.get(0,3));

  temp.set(1,2,i2*(d*(mat12*dp)));
  temp.set(2,1,-temp.get(1,2));

  temp.set(1,3,i2*(d*(mat13*dp)));
  temp.set(3,1,-temp.get(1,3));

  temp.set(2,3,i2*(d*(mat23*dp)));
  temp.set(3,2,-temp.get(2,3));

  return temp;
}
