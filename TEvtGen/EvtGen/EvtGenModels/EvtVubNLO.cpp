//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtVubNLO.cc
//
// Description: Routine to decay B->Xulnu according to Bosch, Lange, Neubert, and Paz hep-ph/0402094
//              Equation numbers refer to this paper
//
// Modification history:
//
//    Riccardo Faccini       Feb. 11, 2004       
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtVubNLO.hh"
#include <string>
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenModels/EvtItgSimpsonIntegrator.hh"
#include "EvtGenModels/EvtBtoXsgammaFermiUtil.hh"
#include "EvtGenModels/EvtItgPtrFunction.hh"
#include "EvtGenModels/EvtPFermi.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtDiLog.hh"

using std::cout;
using std::endl;

EvtVubNLO::~EvtVubNLO() {
  delete [] _masses;
  delete [] _weights;
  cout <<" max pdf : "<<_gmax<<endl;
  cout <<" efficiency : "<<(float)_ngood/(float)_ntot<<endl;

}


std::string EvtVubNLO::getName(){

  return "VUB_NLO";     

}

EvtDecayBase* EvtVubNLO::clone(){

  return new EvtVubNLO;

}


void EvtVubNLO::init(){

  // max pdf
  _gmax=0;
  _ntot=0;
  _ngood=0;
  _lbar=-1000;
  _mupi2=-1000;

  // check that there are at least 6 arguments
  int npar = 8;
  if (getNArg()<npar) {

    report(Severity::Error,"EvtGen") << "EvtVubNLO generator expected "
                           << " at least npar arguments  but found: "
			   <<getNArg()<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();

  }
  // this is the shape function parameter
  _mb      	= getArg(0);
  _b       	= getArg(1);
  _lambdaSF  	= getArg(2);// shape function lambda is different from lambda
  _mui  	= 1.5;// GeV (scale)
  _kpar         = getArg(3);// 0
  _idSF         = abs((int)getArg(4));// type of shape function 1: exponential (from Neubert)
  _nbins        = abs((int)getArg(5));
  _masses       = new double[_nbins];
  _weights      = new double[_nbins];

  // Shape function normalization
  _mB=5.28;// temporary B meson mass for normalization

  std::vector<double> sCoeffs(11);
  sCoeffs[3] = _b;
  sCoeffs[4] = _mb;
  sCoeffs[5] = _mB;
  sCoeffs[6] = _idSF;
  sCoeffs[7] =  lambda_SF();
  sCoeffs[8] =  mu_h();
  sCoeffs[9] =  mu_i();
  sCoeffs[10] = 1.;
  _SFNorm = SFNorm(sCoeffs) ; // SF normalization;


  cout << " pdf 0.66, 1.32 , 4.32 "<<tripleDiff(0.66, 1.32 , 4.32)<<endl;
  cout << " pdf 0.23,0.37,3.76 "<<tripleDiff(0.23,0.37,3.76)<<endl;
  cout << " pdf 0.97,4.32,4.42 "<<tripleDiff(0.97,4.32,4.42)<<endl;
  cout << " pdf 0.52,1.02,2.01 "<<tripleDiff(0.52,1.02,2.01)<<endl;
  cout << " pdf 1.35,1.39,2.73 "<<tripleDiff(1.35,1.39,2.73)<<endl;

  
  if (getNArg()-npar+2 != 2*_nbins) {
    report(Severity::Error,"EvtGen") << "EvtVubNLO generator expected " 
                           << _nbins << " masses and weights but found: "
			   <<(getNArg()-npar)/2 <<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }
  int i,j = npar-2;
  double maxw = 0.;
  for (i=0;i<_nbins;i++) {
    _masses[i] = getArg(j++);
    if (i>0 && _masses[i] <= _masses[i-1]) {
      report(Severity::Error,"EvtGen") << "EvtVubNLO generator expected " 
			     << " mass bins in ascending order!"
			     << "Will terminate execution!"<<endl;
      ::abort();
    }
    _weights[i] = getArg(j++);
    if (_weights[i] < 0) {
      report(Severity::Error,"EvtGen") << "EvtVubNLO generator expected " 
			     << " weights >= 0, but found: " 
			     <<_weights[i] <<endl;
      report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
      ::abort();
    }
    if ( _weights[i] > maxw ) maxw = _weights[i];
  }
  if (maxw == 0) {
    report(Severity::Error,"EvtGen") << "EvtVubNLO generator expected at least one " 
			   << " weight > 0, but found none! " 
			   << "Will terminate execution!"<<endl;
    ::abort();
  }
  for (i=0;i<_nbins;i++) _weights[i]/=maxw;

  // the maximum dGamma*p2 value depends on alpha_s only:


  //  _dGMax = 0.05;
   _dGMax = 150.;

  // for the Fermi Motion we need a B-Meso\n mass - but it's not critical
  // to get an exact value; in order to stay in the phase space for
  // B+- and B0 use the smaller mass

  
  // check that there are 3 daughters
  checkNDaug(3);
}

void EvtVubNLO::initProbMax(){
  
  noProbMax();
  
}

void EvtVubNLO::decay( EvtParticle *p ){
  
  int j;
  // B+ -> u-bar specflav l+ nu
  
  EvtParticle *xuhad, *lepton, *neutrino;
  EvtVector4R p4;
  
  double pp,pm,pl,ml,El(0.0),Eh(0.0),sh(0.0);
  
  
  
  p->initializePhaseSpace(getNDaug(),getDaugs());
  
  xuhad=p->getDaug(0);
  lepton=p->getDaug(1);
  neutrino=p->getDaug(2);
  
  _mB = p->mass();
  ml = lepton->mass();
  
  bool tryit = true;
  
  while (tryit) {
    // pm=(E_H+P_H)
    pm= EvtRandom::Flat(0.,1);
    pm= pow(pm,1./3.)*_mB;
    // pl=mB-2*El
    pl = EvtRandom::Flat(0.,1);
    pl=sqrt(pl)*pm;
    // pp=(E_H-P_H)    
    pp = EvtRandom::Flat(0.,pl);

    _ntot++;
    
    El = (_mB-pl)/2.;      
    Eh = (pp+pm)/2;
    sh = pp*pm;
     
    double pdf(0.);
    if (pp<pl && El>ml&& sh > _masses[0]*_masses[0]&& _mB*_mB + sh - 2*_mB*Eh > ml*ml) {
      double xran = EvtRandom::Flat(0,_dGMax);
      pdf = tripleDiff(pp,pl,pm); // triple differential distribution
      //      cout <<" P+,P-,Pl,Pdf= "<<pp <<" "<<pm<<" "<<pl<<" "<<pdf<<endl;
      if(pdf>_dGMax){
	report(Severity::Error,"EvtGen") << "EvtVubNLO pdf above maximum: " <<pdf
			       <<" P+,P-,Pl,Pdf= "<<pp <<" "<<pm<<" "<<pl<<" "<<pdf<<endl;
	//::abort();
	
      }
      if ( pdf >= xran ) tryit = false;
	
      if(pdf>_gmax)_gmax=pdf;
    } else {
      //      cout <<" EvtVubNLO incorrect kinematics  sh= "<<sh<<"EH "<<Eh<<endl;
    }	
      
    
    // reweight the Mx distribution
    if(!tryit && _nbins>0){
      _ngood++;
      double xran1 = EvtRandom::Flat();
      double m = sqrt(sh);j=0;
      while ( j < _nbins && m > _masses[j] ) j++; 
      double w = _weights[j-1]; 
      if ( w < xran1 ) tryit = true;// through away this candidate
    }
  }

  //  cout <<" max prob "<<gmax<<" " << pp<<" "<<y<<" "<<x<<endl;
  
  // o.k. we have the three kineamtic variables 
  // now calculate a flat cos Theta_H [-1,1] distribution of the 
  // hadron flight direction w.r.t the B flight direction 
  // because the B is a scalar and should decay isotropic.
  // Then chose a flat Phi_H [0,2Pi] w.r.t the B flight direction 
  // and and a flat Phi_L [0,2Pi] in the W restframe w.r.t the 
  // W flight direction.
  
  double ctH = EvtRandom::Flat(-1,1);
  double phH = EvtRandom::Flat(0,2*M_PI);
  double phL = EvtRandom::Flat(0,2*M_PI);

  // now compute the four vectors in the B Meson restframe
    
  double ptmp,sttmp;
  // calculate the hadron 4 vector in the B Meson restframe
  
  sttmp = sqrt(1-ctH*ctH);
  ptmp = sqrt(Eh*Eh-sh);
  double pHB[4] = {Eh,ptmp*sttmp*cos(phH),ptmp*sttmp*sin(phH),ptmp*ctH};
  p4.set(pHB[0],pHB[1],pHB[2],pHB[3]);
  xuhad->init( getDaug(0), p4);


  // calculate the W 4 vector in the B Meson restrframe

  double apWB = ptmp;
  double pWB[4] = {_mB-Eh,-pHB[1],-pHB[2],-pHB[3]};

  // first go in the W restframe and calculate the lepton and
  // the neutrino in the W frame

  double mW2   = _mB*_mB + sh - 2*_mB*Eh;
  //  if(mW2<0.1){
  //  cout <<" low Q2! "<<pp<<" "<<epp<<" "<<x<<" "<<y<<endl;
  //}
  double beta  = ptmp/pWB[0];
  double gamma = pWB[0]/sqrt(mW2);

  double pLW[4];
    
  ptmp = (mW2-ml*ml)/2/sqrt(mW2);
  pLW[0] = sqrt(ml*ml + ptmp*ptmp);

  double ctL = (El - gamma*pLW[0])/beta/gamma/ptmp;
  if ( ctL < -1 ) ctL = -1;
  if ( ctL >  1 ) ctL =  1;
  sttmp = sqrt(1-ctL*ctL);

  // eX' = eZ x eW
  double xW[3] = {-pWB[2],pWB[1],0};
  // eZ' = eW
  double zW[3] = {pWB[1]/apWB,pWB[2]/apWB,pWB[3]/apWB};
  
  double lx = sqrt(xW[0]*xW[0]+xW[1]*xW[1]);
  for (j=0;j<2;j++) 
    xW[j] /= lx;

  // eY' = eZ' x eX'
  double yW[3] = {-pWB[1]*pWB[3],-pWB[2]*pWB[3],pWB[1]*pWB[1]+pWB[2]*pWB[2]};
  double ly = sqrt(yW[0]*yW[0]+yW[1]*yW[1]+yW[2]*yW[2]);
  for (j=0;j<3;j++) 
    yW[j] /= ly;

  // p_lep = |p_lep| * (  sin(Theta) * cos(Phi) * eX'
  //                    + sin(Theta) * sin(Phi) * eY'
  //                    + cos(Theta) *            eZ')
  for (j=0;j<3;j++)
    pLW[j+1] = sttmp*cos(phL)*ptmp*xW[j] 
      +        sttmp*sin(phL)*ptmp*yW[j]
      +          ctL         *ptmp*zW[j];

  double apLW = ptmp;
    
  // boost them back in the B Meson restframe
  
  double appLB = beta*gamma*pLW[0] + gamma*ctL*apLW;
 
  ptmp = sqrt(El*El-ml*ml);
  double ctLL = appLB/ptmp;

  if ( ctLL >  1 ) ctLL =  1;
  if ( ctLL < -1 ) ctLL = -1;
    
  double pLB[4] = {El,0,0,0};
  double pNB[8] = {pWB[0]-El,0,0,0};

  for (j=1;j<4;j++) {
    pLB[j] = pLW[j] + (ctLL*ptmp - ctL*apLW)/apWB*pWB[j];
    pNB[j] = pWB[j] - pLB[j];
  }

  p4.set(pLB[0],pLB[1],pLB[2],pLB[3]);
  lepton->init( getDaug(1), p4);

  p4.set(pNB[0],pNB[1],pNB[2],pNB[3]);
  neutrino->init( getDaug(2), p4);

  return ;
}

double
EvtVubNLO::tripleDiff (  double pp, double pl, double pm){

  std::vector<double> sCoeffs(11);
  sCoeffs[0] = pp;
  sCoeffs[1] = pl;
  sCoeffs[2] = pm;
  sCoeffs[3] = _b;
  sCoeffs[4] =  _mb;
  sCoeffs[5] = _mB;
  sCoeffs[6] = _idSF;
  sCoeffs[7] =  lambda_SF();
  sCoeffs[8] =  mu_h();
  sCoeffs[9] =  mu_i();
  sCoeffs[10] =  _SFNorm; // SF normalization;
  

  double c1=(_mB+pl-pp-pm)*(pm-pl);
  double c2=2*(pl-pp)*(pm-pl);
  double c3=(_mB-pm)*(pm-pp);
  double aF1=F10(sCoeffs);
  double aF2=F20(sCoeffs);
  double aF3=F30(sCoeffs);
  double td0=c1*aF1+c2*aF2+c3*aF3;


  EvtItgPtrFunction *func = new EvtItgPtrFunction(&integrand, 0., _mB, sCoeffs);
  EvtItgAbsIntegrator *jetSF = new EvtItgSimpsonIntegrator(*func,0.01,25);
  double smallfrac=0.000001;// stop a bit before the end to avoid problems with numerical integration
  double tdInt = jetSF->evaluate(0,pp*(1-smallfrac));
  delete jetSF;
  
  double SU=U1lo(mu_h(),mu_i())*pow((pm-pp)/(_mB-pp),alo(mu_h(),mu_i()));
  double TD=(_mB-pp)*SU*(td0+tdInt);

  return TD;

}

double
EvtVubNLO::integrand(double omega, const std::vector<double> &coeffs){
  //double pp=coeffs[0];
  double c1=(coeffs[5]+coeffs[1]-coeffs[0]-coeffs[2])*(coeffs[2]-coeffs[1]);
  double c2=2*(coeffs[1]-coeffs[0])*(coeffs[2]-coeffs[1]);
  double c3=(coeffs[5]-coeffs[2])*(coeffs[2]-coeffs[0]);

  return  c1*F1Int(omega,coeffs)+c2*F2Int(omega,coeffs)+c3*F3Int(omega,coeffs);  
}

double 
EvtVubNLO::F10(const std::vector<double> &coeffs){
  double pp=coeffs[0];
  double y=(coeffs[2]-coeffs[0])/(coeffs[5]-coeffs[0]);
  double mui=coeffs[9];
  double muh=coeffs[8];
  double z=1-y;
  double result= U1nlo(muh,mui)/ U1lo(muh,mui);

  result += anlo(muh,mui)*log(y);

  result += C_F(muh)*(-4*pow(log(y*coeffs[4]/muh),2)+10*log(y*coeffs[4]/muh)-4*log(y)-2*log(y)/(1-y)-4.0*EvtDiLog::DiLog(z)-pow(EvtConst::pi,2)/6.-12  );

  result += C_F(mui)*(2*pow(log(y*coeffs[4]*pp/pow(mui,2)),2)-3*log(y*coeffs[4]*pp/pow(mui,2))+7-pow(EvtConst::pi,2) );
  result *=shapeFunction(pp,coeffs);
  // changes due to SSF 
  result +=  (-subS(coeffs)+2*subT(coeffs)+(subU(coeffs)-subV(coeffs))*(1/y-1.))/(coeffs[5]-pp);
  result += shapeFunction(pp,coeffs)/pow((coeffs[5]-coeffs[0]),2)*(-5*(lambda1()+3*lambda2())/6+2*(2*lambda1()/3-lambda2())/pow(y,2));
  //  result +=  (subS(coeffs)+subT(coeffs)+(subU(coeffs)-subV(coeffs))/y)/(coeffs[5]-pp);
  // this part has been added after Feb '05
  
  //result += shapeFunction(pp,coeffs)/pow((coeffs[5]-coeffs[0]),2)*((lambda1()+3*lambda2())/6+2*(2*lambda1()/3-lambda2())/pow(y,2));
  return result;
}

double 
EvtVubNLO::F1Int(double omega,const std::vector<double> &coeffs){
  double pp=coeffs[0];
  double y=(coeffs[2]-coeffs[0])/(coeffs[5]-coeffs[0]);
  // mubar == mui
  return C_F(coeffs[9])*(
			 (shapeFunction(omega,coeffs)-shapeFunction(pp,coeffs))*(4*log(y*coeffs[4]*(pp-omega)/pow(coeffs[9],2))-3)/(pp-omega)+
			 (g1(y,(pp-omega)/(coeffs[5]-coeffs[0]))/(coeffs[5]-pp)*shapeFunction(omega,coeffs))
			 );  
}

double 
EvtVubNLO::F20(const std::vector<double> &coeffs){
  double pp=coeffs[0];
  double y=(coeffs[2]-coeffs[0])/(coeffs[5]-coeffs[0]);
  double result= C_F(coeffs[8])*log(y)/(1-y)*shapeFunction(pp,coeffs)-
    1/y*(subS(coeffs)+2*subT(coeffs)-(subT(coeffs)+subV(coeffs))/y)/(coeffs[5]-pp);
  // added after Feb '05
  result += shapeFunction(pp,coeffs)/pow((coeffs[5]-coeffs[0])*y,2)*(2*lambda1()/3+4*lambda2()-y*(7/6*lambda1()+3*lambda2()));
  return result;
}

double 
EvtVubNLO::F2Int(double omega,const std::vector<double> &coeffs){
  double pp=coeffs[0];
  double y=(coeffs[2]-coeffs[0])/(coeffs[5]-coeffs[0]);
  return C_F(coeffs[9])*g3(y,(pp-omega)/(coeffs[5]-coeffs[0]))*shapeFunction(omega,coeffs)/(coeffs[5]-pp);
}

double 
EvtVubNLO::F30(const std::vector<double> &coeffs){
  double y=(coeffs[2]-coeffs[0])/(coeffs[5]-coeffs[0]);
  return shapeFunction(coeffs[0],coeffs)/pow((coeffs[5]-coeffs[0])*y,2)*(-2*lambda1()/3+lambda2());
}

double 
EvtVubNLO::F3Int(double omega,const std::vector<double> &coeffs){
  double pp=coeffs[0];
  double y=(coeffs[2]-coeffs[0])/(coeffs[5]-coeffs[0]);
  return C_F(coeffs[9])*g3(y,(pp-omega)/(coeffs[5]-coeffs[0]))/2*shapeFunction(omega,coeffs)/(coeffs[2]-coeffs[0]);
}

double 
EvtVubNLO::g1(double y, double x){
  double result=(y*(-9+10*y)+x*x*(-12+13*y)+2*x*(-8+6*y+3*y*y))/y/pow(1+x,2)/(x+y);
  result -= 4*log((1+1/x)*y)/x;
  result -=2*log(1+y/x)*(3*pow(x,4)*(-2+y)-2*pow(y,3)-4*pow(x,3)*(2+y)-2*x*y*y*(4+y)-x*x*y*(12+4*y+y*y))/x/pow((1+x)*y,2)/(x+y);
  return result;
}

double 
EvtVubNLO::g2(double y, double x){
  double result=y*(10*pow(x,4)+y*y+3*x*x*y*(10+y)+pow(x,3)*(12+19*y)+x*y*(8+4*y+y*y));
  result -= 2*x*log(1+y/x)*(5*pow(x,4)+2*y*y+6*pow(x,3)*(1+2*y)+4*y*x*(1+2*y)+x*x*y*(18+5*y));
  result *= 2/(pow(y*(1+x),2)*y*(x+y));
  return result;
}

double 
EvtVubNLO::g3(double y, double x){
  double result=(2*pow(y,3)*(-11+2*y)-10*pow(x,4)*(6-6*y+y*y)+x*y*y*(-94+29*y+2*y*y)+2*x*x*y*(-72+18*y+13*y*y)-pow(x,3)*(72+42*y-70*y*y+3*pow(y,3)))/(pow(y*(1+x),2)*y*(x+y));
  result += 2*log(1+y/x)*(-6*x*pow(y,3)*(-5+y)+4*pow(y,4)+5*pow(x,5)*(6-6*y+y*y)-4*pow(x*y,2)*(-20+6*y+y*y)+pow(x,3)*y*(90-10*y-28*y*y+pow(y,3))+pow(x,4)*(36+36*y-50*y*y+4*pow(y,3)))/(pow((1+x)*y*y,2)*(x+y));
  return result;
}

/* old version (before Feb 05 notebook from NNeubert

double 
EvtVubNLO::F1Int(double omega,const std::vector<double> &coeffs){
  double pp=coeffs[0];
  double y=(coeffs[2]-coeffs[0])/(coeffs[5]-coeffs[0]);
  // mubar == mui
  return C_F(coeffs[9])*(
			 (shapeFunction(omega,coeffs)-shapeFunction(pp,coeffs))*(4*log(y*coeffs[4]*(pp-omega)/pow(coeffs[9],2))-3)/(pp-omega)-
			 (1./y/(coeffs[5]-pp)*shapeFunction(omega,coeffs)*(5-6*y+4*(3-y)*log((pp-omega)/y/coeffs[4])))
			 );  
}


double 
EvtVubNLO::F2Int(double omega,const std::vector<double> &coeffs){
  double pp=coeffs[0];
  double y=(coeffs[2]-coeffs[0])/(coeffs[5]-coeffs[0]);
  return C_F(coeffs[9])*shapeFunction(omega,coeffs)*(2-11/y-4/y*log((pp-omega)/y/coeffs[4]))/(coeffs[5]-pp);
}

double 
EvtVubNLO::F3(const std::vector<double> &coeffs){
  return C_F(coeffs[9])*shapeFunction(omega,coeffs)/(coeffs[2]-coeffs[0]);
}
*/

double EvtVubNLO::SFNorm(  const std::vector<double> &/*coeffs*/){
  
  double omega0=1.68;//normalization scale (mB-2*1.8)
  if(_idSF==1){ // exponential SF
    double omega0=1.68;//normalization scale (mB-2*1.8)
    return M0(mu_i(),omega0)*pow(_b,_b)/lambda_SF()/ (Gamma(_b)-Gamma(_b,_b*omega0/lambda_SF()));
  } else if(_idSF==2){ // Gaussian SF
    double c=cGaus(_b);
    return M0(mu_i(),omega0)*2/lambda_SF()/pow(c,-(1+_b)/2.)/
      (Gamma((1+_b)/2)-Gamma((1+_b)/2,pow(omega0/lambda_SF(),2)*c));
  } else {
    report(Severity::Error,"EvtGen")  << "unknown SF "<<_idSF<<endl;
    return -1;
  }
}

double
EvtVubNLO::shapeFunction ( double omega, const std::vector<double> &sCoeffs){
  if( sCoeffs[6]==1){
    return sCoeffs[10]*expShapeFunction(omega, sCoeffs);
  } else if( sCoeffs[6]==2) {
    return sCoeffs[10]*gausShapeFunction(omega, sCoeffs);
  } else {
 report(Severity::Error,"EvtGen") << "EvtVubNLO : unknown shape function # "
			   <<sCoeffs[6]<<endl;
  }
  return -1.;
}


// SSF
double 
EvtVubNLO::subS(const std::vector<double> &c){ return (lambda_bar(1.68)-c[0])*shapeFunction(c[0],c);}
double 
EvtVubNLO::subT(const std::vector<double> &c){  return -3*lambda2()*subS(c)/mu_pi2(1.68);}
double 
EvtVubNLO::subU(const std::vector<double> &c){ return -2*subS(c);}
double 
EvtVubNLO::subV(const std::vector<double> &c){ return -subT(c);}


double 
EvtVubNLO::lambda_bar(double omega0){
  if(_lbar<0){
    if(_idSF==1){ // exponential SF
      double rat=omega0*_b/lambda_SF();
      _lbar=lambda_SF()/_b*(Gamma(1+_b)-Gamma(1+_b,rat))/(Gamma(_b)-Gamma(_b,rat));
    } else if(_idSF==2){ // Gaussian SF
      double c=cGaus(_b);
      _lbar=lambda_SF()*(Gamma(1+_b/2)-Gamma(1+_b/2,pow(omega0/lambda_SF(),2)*c))/(Gamma((1+_b)/2)-Gamma((1+_b)/2,pow(omega0/lambda_SF(),2)*c))/sqrt(c);
    }
  }
  return _lbar;
}


double 
EvtVubNLO::mu_pi2(double omega0){
  if(_mupi2<0){
    if(_idSF==1){ // exponential SF
      double rat=omega0*_b/lambda_SF();
      _mupi2= 3*(pow(lambda_SF()/_b,2)*(Gamma(2+_b)-Gamma(2+_b,rat))/(Gamma(_b)-Gamma(_b,rat))-pow(lambda_bar(omega0),2));
    } else if(_idSF==2){ // Gaussian SF
      double c=cGaus(_b);
      double m1=Gamma((3+_b)/2)-Gamma((3+_b)/2,pow(omega0/lambda_SF(),2)*c);
      double m2=Gamma(1+_b/2)-Gamma(1+_b/2,pow(omega0/lambda_SF(),2)*c);
      double m3=Gamma((1+_b)/2)-Gamma((1+_b)/2,pow(omega0/lambda_SF(),2)*c);
      _mupi2= 3*pow(lambda_SF(),2)*(m1/m3-pow(m2/m3,2))/c;
    }
  }
  return _mupi2;
}

double 
EvtVubNLO::M0(double mui,double omega0){
  double mf=omega0-lambda_bar(omega0);
  return 1+4*C_F(mui)*(-pow(log(mf/mui),2)-log(mf/mui)-pow(EvtConst::pi/2,2)/6.+mu_pi2(omega0)/3/pow(mf,2)*(log(mf/mui)-0.5));
}

double 
EvtVubNLO::alphas(double mu){
  double Lambda4=0.302932;
  double lg=2*log(mu/Lambda4);
  return 4*EvtConst::pi/lg/beta0()*(1-beta1()*log(lg)/pow(beta0(),2)/lg+pow(beta1()/lg,2)/pow(beta0(),4)*(pow(log(lg)-0.5,2)-1.25+beta2()*beta0()/pow(beta1(),2)));
}

double
EvtVubNLO::gausShapeFunction ( double omega, const std::vector<double> &sCoeffs){
  double b=sCoeffs[3];
  double l=sCoeffs[7];
  double wL=omega/l;

  return pow(wL,b)*exp(-cGaus(b)*wL*wL);
}

double
EvtVubNLO::expShapeFunction ( double omega, const std::vector<double> &sCoeffs){
  double b=sCoeffs[3];
  double l=sCoeffs[7];
  double wL=omega/l;

  return pow(wL,b-1)*exp(-b*wL);
}

double 
EvtVubNLO::Gamma(double z) {

  std::vector<double> gammaCoeffs(6);
  gammaCoeffs[0]=76.18009172947146;
  gammaCoeffs[1]=-86.50532032941677;
  gammaCoeffs[2]=24.01409824083091;
  gammaCoeffs[3]=-1.231739572450155;
  gammaCoeffs[4]=0.1208650973866179e-2;
  gammaCoeffs[5]=-0.5395239384953e-5;
  
  //Lifted from Numerical Recipies in C
  double x, y, tmp, ser;
 
  int j;
  y = z;
  x = z;
  
  tmp = x + 5.5;
  tmp = tmp - (x+0.5)*log(tmp);
  ser=1.000000000190015;

  for (j=0;j<6;j++) {
    y = y +1.0;
    ser = ser + gammaCoeffs[j]/y;
  }

  return exp(-tmp+log(2.5066282746310005*ser/x));

}



double 
EvtVubNLO::Gamma(double z, double tmin) {
  std::vector<double> c(1);
  c[0]=z;
  EvtItgPtrFunction *func = new EvtItgPtrFunction(&dgamma, tmin, 100., c);
  EvtItgAbsIntegrator *jetSF = new EvtItgSimpsonIntegrator(*func,0.001);
  return jetSF->evaluate(tmin,100.);
}
