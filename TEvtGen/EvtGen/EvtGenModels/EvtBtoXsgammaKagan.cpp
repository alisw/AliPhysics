//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtBtoXsgammaKagan.cc
//
// Description:
//       Routine to perform two-body non-resonant B->Xs,gamma decays.
//       The X_s mass spectrum generated is based on the Kagan-Neubert model. 
//       See hep-ph/9805303 for the model details and input parameters.
//
//       The input parameters are 1:fermi_model, 2:mB, 3:mb, 4:mu, 5:lam1, 
//       6:delta, 7:z, 8:nIntervalS, 9:nIntervalmH. Choosing fermi_model=1 
//       uses an exponential shape function, fermi_model=2 uses a gaussian 
//       shape function and fermi_model=3 a roman shape function. The complete mass
//       spectrum for a given set of input parameters is calculated from 
//       scratch in bins of nIntervalmH. The s22, s27 and s28 coefficients are calculated
//       in bins of nIntervalS. As the program includes lots of integration, the
//       theoretical hadronic mass spectra is computed for the first time
//       the init method is called. Then, all the other times (eg if we want to decay a B0 
//       as well as an anti-B0) the vector mass info stored the first time is used again.
//
// Modification history:
//
//      Jane Tinslay, Francesca Di Lodovico  March 21, 2001       Module created
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"

#include <stdlib.h>
#include "EvtGenModels/EvtBtoXsgamma.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenModels/EvtBtoXsgammaKagan.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtItgSimpsonIntegrator.hh"
#include "EvtGenModels/EvtItgFunction.hh"
#include "EvtGenModels/EvtItgPtrFunction.hh"
#include "EvtGenModels/EvtItgTwoCoeffFcn.hh"
#include "EvtGenModels/EvtItgThreeCoeffFcn.hh"
#include "EvtGenModels/EvtItgFourCoeffFcn.hh"
#include "EvtGenModels/EvtItgAbsIntegrator.hh"
#include "EvtGenModels/EvtBtoXsgammaFermiUtil.hh"

#include <fstream>
using std::endl;
using std::fstream;

bool EvtBtoXsgammaKagan::bbprod = false;
double EvtBtoXsgammaKagan::intervalMH = 0;

EvtBtoXsgammaKagan::~EvtBtoXsgammaKagan(){
  delete [] massHad;
  delete [] brHad;
}

void EvtBtoXsgammaKagan::init(int nArg, double* args){

  if ((nArg) > 12 || (nArg > 1 && nArg <10) || nArg == 11){
  
  report(Severity::Error,"EvtGen") << "EvtBtoXsgamma generator model "
			 << "EvtBtoXsgammaKagan expected " 
			 << "either 1(default config) or " 
			 << "10 (default mass range) or " 
			 << "12 (user range) arguments but found: "
			 <<nArg<<endl;
  report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();  
  }
  
  if(nArg == 1){
    bbprod = true;
    getDefaultHadronicMass();
  }else{
    bbprod = false;
    computeHadronicMass(nArg, args);
  }

  double mHminLimit=0.6373;
  double mHmaxLimit=4.5;

  if (nArg>10){
    _mHmin = args[10];
    _mHmax = args[11]; 
    if (_mHmin > _mHmax){
      report(Severity::Error,"EvtGen") << "Minimum hadronic mass exceeds maximum " 
			     << endl;
      report(Severity::Error,"EvtGen") << "Will terminate execution!" << endl;
      ::abort();
    }
    if (_mHmin < mHminLimit){
      report(Severity::Error,"EvtGen") << "Minimum hadronic mass below K pi threshold" 
			     << endl;
      report(Severity::Error,"EvtGen") << "Resetting to K pi threshold" << endl;
      _mHmin = mHminLimit;
    }     
    if (_mHmax > mHmaxLimit){
      report(Severity::Error,"EvtGen") << "Maximum hadronic mass above 4.5 GeV/c^2" 
			     << endl;
      report(Severity::Error,"EvtGen") << "Resetting to 4.5 GeV/c^2" << endl;
      _mHmax = mHmaxLimit;
    }     
  }else{
    _mHmin=mHminLimit; //  usually just above K pi threshold for Xsd/u
    _mHmax=mHmaxLimit;    
  }  
  
}

void EvtBtoXsgammaKagan::getDefaultHadronicMass(){

    massHad = new double[81];
    brHad = new double[81];
  
    double mass[81] = { 0, 0.0625995, 0.125199, 0.187798, 0.250398, 0.312997, 0.375597, 0.438196, 0.500796, 0.563395, 0.625995, 0.688594, 0.751194, 0.813793, 0.876392, 0.938992, 1.00159, 1.06419, 1.12679, 1.18939, 1.25199, 1.31459, 1.37719, 1.43979, 1.50239, 1.56499, 1.62759, 1.69019, 1.75278, 1.81538, 1.87798, 1.94058, 2.00318, 2.06578, 2.12838, 2.19098, 2.25358, 2.31618, 2.37878, 2.44138, 2.50398, 2.56658, 2.62918, 2.69178, 2.75438, 2.81698, 2.87958, 2.94217, 3.00477, 3.06737, 3.12997, 3.19257, 3.25517, 3.31777, 3.38037, 3.44297, 3.50557, 3.56817, 3.63077, 3.69337, 3.75597, 3.81857, 3.88117, 3.94377, 4.00637, 4.06896, 4.13156, 4.19416, 4.25676, 4.31936, 4.38196, 4.44456, 4.50716, 4.56976, 4.63236, 4.69496, 4.75756, 4.82016, 4.88276, 4.94536, 5.00796};

    double br[81] = { 0, 1.03244e-09, 3.0239e-08, 1.99815e-07, 7.29392e-07, 1.93129e-06, 4.17806e-06, 7.86021e-06, 1.33421e-05, 2.09196e-05, 3.07815e-05, 4.29854e-05, 5.74406e-05, 7.3906e-05, 9.2003e-05, 0.000111223, 0.000130977, 0.000150618, 0.000169483, 0.000186934, 0.000202392, 0.000215366, 0.000225491, 0.000232496, 0.000236274, 0.000236835, 0.000234313, 0.000228942, 0.000221042, 0.000210994, 0.000199215, 0.000186137, 0.000172194, 0.000157775, 0.000143255, 0.000128952, 0.000115133, 0.000102012, 8.97451e-05, 7.84384e-05, 6.81519e-05, 5.89048e-05, 5.06851e-05, 4.34515e-05, 3.71506e-05, 3.1702e-05, 2.70124e-05, 2.30588e-05, 1.96951e-05, 1.68596e-05, 1.44909e-05, 1.25102e-05, 1.08596e-05, 9.48476e-06, 8.34013e-06, 7.38477e-06, 6.58627e-06, 5.91541e-06, 5.35022e-06, 4.87047e-06, 4.46249e-06, 4.11032e-06, 3.80543e-06, 3.54051e-06, 3.30967e-06, 3.10848e-06, 2.93254e-06, 2.78369e-06, 2.65823e-06, 2.55747e-06, 2.51068e-06, 2.57179e-06, 2.74684e-06, 3.02719e-06, 3.41182e-06, 3.91387e-06, 4.56248e-06, 5.40862e-06, 6.53915e-06, 8.10867e-06, 1.04167e-05 };

  for(int i=0; i<81; i++){
    massHad[i] = mass[i];
    brHad[i] = br[i];
  }
  intervalMH=80;
}

void EvtBtoXsgammaKagan::computeHadronicMass(int /*nArg*/, double* args){

  //Input parameters
  int fermiFunction = (int)args[1];
  _mB = args[2];
  _mb = args[3];
  _mu = args[4];
  _lam1 = args[5];
  _delta = args[6];
  _z = args[7];
  _nIntervalS = args[8];
  _nIntervalmH = args[9];
  std::vector<double> mHVect(int(_nIntervalmH+1.0));
  massHad = new double[int(_nIntervalmH+1.0)];
  brHad = new double[int(_nIntervalmH+1.0)];
  intervalMH=_nIntervalmH;

  //Going to have to add a new entry into the data file - takes ages...
  report(Severity::Warning,"EvtGen") << "EvtBtoXsgammaKagan: calculating new hadronic mass spectra. This takes a while..." << endl;
  
  //Now need to compute the mHVect vector for
  //the current parameters
  
  //A few more parameters
  double _mubar = _mu;
  _mW = 80.33;
  _mt = 175.0;
  _alpha = 1./137.036;
  _lambdabar = _mB - _mb;
  _kappabar = 3.382 - 4.14*(sqrt(_z) - 0.29);
  _fz=Fz(_z);
  _rer8 = (44./9.) - (8./27.)*pow(EvtConst::pi,2.);
  _r7 = (-10./3.) - (8./9.)*pow(EvtConst::pi,2.);
  _rer2 = -4.092 + 12.78*(sqrt(_z) -.29);
  _gam77 = 32./3.;
  _gam27 = 416./81.;
  _gam87 = -32./9.;
  _lam2 = .12;
  _beta0 = 23./3.;
  _beta1 = 116./3.;
  _alphasmZ = .118;
  _mZ = 91.187;
  _ms = _mb/50.;
  
  double eGammaMin = 0.5*_mB*(1. - _delta);
  double eGammaMax = 0.5*_mB;
  double yMin = 2.*eGammaMin/_mB;
  double yMax = 2.*eGammaMax/_mB;
  double _CKMrat= 0.976;
  double Nsl = 1.0;
  
  //Calculate alpha the various scales
  _alphasmW = CalcAlphaS(_mW);
  _alphasmt = CalcAlphaS(_mt);
  _alphasmu = CalcAlphaS(_mu);
  _alphasmubar = CalcAlphaS(_mubar);
  
  //Calculate the Wilson Coefficients and Delta 
  _etamu = _alphasmW/_alphasmu;
  _kSLemmu = (12./23.)*((1./_etamu) -1.);
  CalcWilsonCoeffs();
  CalcDelta();
  
  //Build s22 and s27 vector - saves time because double
  //integration is required otherwise
  std::vector<double> s22Coeffs(int(_nIntervalS+1.0));
  std::vector<double> s27Coeffs(int(_nIntervalS+1.0));
  std::vector<double> s28Coeffs(int(_nIntervalS+1.0));
  
  double dy = (yMax - yMin)/_nIntervalS;
  double yp = yMin;
  
  std::vector<double> sCoeffs(1);
  sCoeffs[0] = _z;
  
  //Define s22 and s27 functions
  EvtItgPtrFunction *mys22Func = new EvtItgPtrFunction(&s22Func, 0., yMax+0.1, sCoeffs);
  EvtItgPtrFunction *mys27Func = new EvtItgPtrFunction(&s27Func, 0., yMax+0.1, sCoeffs);
  
  //Use a simpson integrator
  EvtItgAbsIntegrator *mys22Simp = new EvtItgSimpsonIntegrator(*mys22Func, 1.0e-4, 20);
  EvtItgAbsIntegrator *mys27Simp = new EvtItgSimpsonIntegrator(*mys27Func, 1.0e-4, 50);

  int i;

  for (i=0;i<int(_nIntervalS+1.0);i++) {
    
    s22Coeffs[i] = (16./27.)*mys22Simp->evaluate(1.0e-20,yp);
    s27Coeffs[i] = (-8./9.)*_z*mys27Simp->evaluate(1.0e-20,yp);
    s28Coeffs[i] = -s27Coeffs[i]/3.;
    yp = yp + dy;
    
  }

  delete mys22Func;
  delete mys27Func;
  delete mys22Simp;
  delete mys27Simp;
  
  //Define functions and vectors used to calculate mHVect. Each function takes a set
  //of vectors which are used as the function coefficients
  std::vector<double> FermiCoeffs(6);
  std::vector<double> varCoeffs(3);
  std::vector<double> DeltaCoeffs(1);
  std::vector<double> s88Coeffs(2);
  std::vector<double> sInitCoeffs(3);
  
  varCoeffs[0] = _mB;
  varCoeffs[1] = _mb;
  varCoeffs[2] = 0.;
  
  DeltaCoeffs[0] = _alphasmu;
  
  s88Coeffs[0] = _mb;
  s88Coeffs[1] = _ms;
  
  sInitCoeffs[0] = _nIntervalS;
  sInitCoeffs[1] = yMin;
  sInitCoeffs[2] = yMax;
  
  FermiCoeffs[0]=fermiFunction;
  FermiCoeffs[1]=0.0;
  FermiCoeffs[2]=0.0;
  FermiCoeffs[3]=0.0;
  FermiCoeffs[4]=0.0;
  FermiCoeffs[5]=0.0;
  
  //Coefficients for gamma function
  std::vector<double> gammaCoeffs(6);
  gammaCoeffs[0]=76.18009172947146;
  gammaCoeffs[1]=-86.50532032941677;
  gammaCoeffs[2]=24.01409824083091;
  gammaCoeffs[3]=-1.231739572450155;
  gammaCoeffs[4]=0.1208650973866179e-2;
  gammaCoeffs[5]=-0.5395239384953e-5;
  
  //Calculate quantities for the fermi function to be used
  //Distinguish among the different shape functions
  if (fermiFunction == 1) {
    
    FermiCoeffs[1]=_lambdabar;
    FermiCoeffs[2]=(-3.*pow(_lambdabar,2.)/_lam1) - 1.;
    FermiCoeffs[3]=_lam1;
    FermiCoeffs[4]=1.0;
    
    EvtItgPtrFunction *myNormFunc = new EvtItgPtrFunction(&EvtBtoXsgammaFermiUtil::FermiExpFunc, -_mb, _mB-_mb, FermiCoeffs);
    EvtItgAbsIntegrator *myNormSimp = new EvtItgSimpsonIntegrator(*myNormFunc, 1.0e-4, 40);
    FermiCoeffs[4]=myNormSimp->normalisation();
    delete myNormFunc; myNormFunc=0;
    delete myNormSimp; myNormSimp=0;
    
  } else if (fermiFunction == 2) {
    
    double a = EvtBtoXsgammaFermiUtil::FermiGaussFuncRoot(_lambdabar, _lam1, _mb, gammaCoeffs);
    FermiCoeffs[1]=_lambdabar;
    FermiCoeffs[2]=a;
    FermiCoeffs[3]= EvtBtoXsgammaFermiUtil::Gamma((2.0 + a)/2., gammaCoeffs)/
      EvtBtoXsgammaFermiUtil::Gamma((1.0 + a)/2., gammaCoeffs);
    FermiCoeffs[4]=1.0;
    
    EvtItgPtrFunction *myNormFunc = new EvtItgPtrFunction(&EvtBtoXsgammaFermiUtil::FermiGaussFunc, -_mb, _mB-_mb, FermiCoeffs);
    EvtItgAbsIntegrator *myNormSimp = new EvtItgSimpsonIntegrator(*myNormFunc, 1.0e-4, 40);
    FermiCoeffs[4]=myNormSimp->normalisation();
    delete myNormFunc; myNormFunc=0;
    delete myNormSimp; myNormSimp=0;
    
  }
  else if (fermiFunction == 3) {
    
    double rho = EvtBtoXsgammaFermiUtil::FermiRomanFuncRoot(_lambdabar, _lam1);
    FermiCoeffs[1]=_mB;
    FermiCoeffs[2]=_mb;
    FermiCoeffs[3]= rho;
    FermiCoeffs[4]=_lambdabar;
    FermiCoeffs[5]=1.0;
    
    EvtItgPtrFunction *myNormFunc = new EvtItgPtrFunction(&EvtBtoXsgammaFermiUtil::FermiRomanFunc, -_mb, _mB-_mb, FermiCoeffs);
    EvtItgAbsIntegrator *myNormSimp = new EvtItgSimpsonIntegrator(*myNormFunc, 1.0e-4, 40);
    FermiCoeffs[5]=myNormSimp->normalisation();
    delete myNormFunc; myNormFunc=0;
    delete myNormSimp; myNormSimp=0;
    
  }
  
  //Define functions
  EvtItgThreeCoeffFcn* myDeltaFermiFunc = new EvtItgThreeCoeffFcn(&DeltaFermiFunc, -_mb, _mB-_mb, FermiCoeffs, varCoeffs, DeltaCoeffs);
  EvtItgThreeCoeffFcn* mys88FermiFunc = new EvtItgThreeCoeffFcn(&s88FermiFunc, -_mb, _mB-_mb, FermiCoeffs, varCoeffs, s88Coeffs);
  EvtItgTwoCoeffFcn* mys77FermiFunc = new EvtItgTwoCoeffFcn(&s77FermiFunc, -_mb, _mB-_mb, FermiCoeffs, varCoeffs);
  EvtItgTwoCoeffFcn* mys78FermiFunc = new EvtItgTwoCoeffFcn(&s78FermiFunc, -_mb, _mB-_mb, FermiCoeffs, varCoeffs);
  EvtItgFourCoeffFcn* mys22FermiFunc = new EvtItgFourCoeffFcn(&sFermiFunc, -_mb, _mB-_mb, FermiCoeffs, varCoeffs, sInitCoeffs, s22Coeffs);
  EvtItgFourCoeffFcn* mys27FermiFunc = new EvtItgFourCoeffFcn(&sFermiFunc, -_mb, _mB-_mb, FermiCoeffs, varCoeffs, sInitCoeffs, s27Coeffs);
  EvtItgFourCoeffFcn* mys28FermiFunc = new EvtItgFourCoeffFcn(&sFermiFunc, -_mb, _mB-_mb, FermiCoeffs, varCoeffs, sInitCoeffs, s28Coeffs);
  
  //Define integrators
  EvtItgSimpsonIntegrator* myDeltaFermiSimp = 
    new EvtItgSimpsonIntegrator(*myDeltaFermiFunc, 1.0e-4, 40);
  EvtItgSimpsonIntegrator* mys77FermiSimp = 
    new EvtItgSimpsonIntegrator(*mys77FermiFunc, 1.0e-4, 40);
  EvtItgSimpsonIntegrator* mys88FermiSimp = 
    new EvtItgSimpsonIntegrator(*mys88FermiFunc, 1.0e-4, 40);
  EvtItgSimpsonIntegrator* mys78FermiSimp = 
    new EvtItgSimpsonIntegrator(*mys78FermiFunc, 1.0e-4, 40);
  EvtItgSimpsonIntegrator* mys22FermiSimp = 
    new EvtItgSimpsonIntegrator(*mys22FermiFunc, 1.0e-4, 40);
  EvtItgSimpsonIntegrator* mys27FermiSimp = 
    new EvtItgSimpsonIntegrator(*mys27FermiFunc, 1.0e-4, 40);
  EvtItgSimpsonIntegrator* mys28FermiSimp = 
    new EvtItgSimpsonIntegrator(*mys28FermiFunc, 1.0e-4, 40);
  
  //Finally calculate mHVect for the range of hadronic masses
  double mHmin = sqrt(_mB*_mB - 2.*_mB*eGammaMax);
  double mHmax = sqrt(_mB*_mB - 2.*_mB*eGammaMin);
  double dmH = (mHmax - mHmin)/_nIntervalmH;
  
  double mH=mHmin;

  //Calculating the Branching Fractions
  for (i=0;i<int(_nIntervalmH+1.0);i++) {
    
    double ymH = 1. - ((mH*mH)/(_mB*_mB));
    
    //Need to set ymH as one of the input parameters
    myDeltaFermiFunc->setCoeff(2, 2, ymH);
    mys77FermiFunc->setCoeff(2, 2, ymH);
    mys88FermiFunc->setCoeff(2, 2, ymH);
    mys78FermiFunc->setCoeff(2, 2, ymH);
    mys22FermiFunc->setCoeff(2, 2, ymH);
    mys27FermiFunc->setCoeff(2, 2, ymH);
    mys28FermiFunc->setCoeff(2, 2, ymH);
    
    //Integrate
    
    double deltaResult = myDeltaFermiSimp->evaluate((_mB*ymH-_mb),_mB-_mb);
    double s77Result = mys77FermiSimp->evaluate((_mB*ymH-_mb),_mB-_mb);
    double s88Result = mys88FermiSimp->evaluate((_mB*ymH-_mb),_mB-_mb);
    double s78Result = mys78FermiSimp->evaluate((_mB*ymH-_mb),_mB-_mb);
    double s22Result = mys22FermiSimp->evaluate((_mB*ymH-_mb),_mB-_mb);
    double s27Result = mys27FermiSimp->evaluate((_mB*ymH-_mb),_mB-_mb);
    mys28FermiSimp->evaluate((_mB*ymH-_mb),_mB-_mb);
    
    double py = (pow(_CKMrat,2.)*(6./_fz)*(_alpha/EvtConst::pi)*(deltaResult*_cDeltatot  + (_alphasmu/EvtConst::pi)*(s77Result*pow(_c70mu,2.) + s27Result*_c2mu*(_c70mu  - _c80mu/3.) + s78Result*_c70mu*_c80mu + s22Result*_c2mu*_c2mu  + s88Result*_c80mu*_c80mu )  ) );
    
    mHVect[i] = 2.*(mH/(_mB*_mB))*0.105*Nsl*py;

    massHad[i] = mH;
    brHad[i] =  2.*(mH/(_mB*_mB))*0.105*Nsl*py;

    mH = mH+dmH;
    
  }

  //Clean up
  delete  myDeltaFermiFunc; myDeltaFermiFunc=0;
  delete mys88FermiFunc; mys88FermiFunc=0;
  delete mys77FermiFunc; mys77FermiFunc=0;
  delete mys78FermiFunc; mys78FermiFunc=0;
  delete mys22FermiFunc; mys22FermiFunc=0;
  delete mys27FermiFunc; mys27FermiFunc=0;
  delete mys28FermiFunc; mys28FermiFunc=0;
  
  delete myDeltaFermiSimp; myDeltaFermiSimp=0;
  delete mys77FermiSimp; mys77FermiSimp=0;
  delete mys88FermiSimp; mys88FermiSimp=0;
  delete mys78FermiSimp; mys78FermiSimp=0;
  delete mys22FermiSimp; mys22FermiSimp=0;
  delete mys27FermiSimp; mys27FermiSimp=0;
  delete mys28FermiSimp; mys28FermiSimp=0;
  
}

double EvtBtoXsgammaKagan::GetMass( int /*Xscode*/ ){
 
//  Get hadronic mass for the event according to the hadronic mass spectra computed in computeHadronicMass
  double mass=0.0;
  double min=_mHmin;
  if(bbprod)min=1.1;
  //  double max=4.5;
  double max=_mHmax;
  double xbox(0), ybox(0);
  double boxheight(0);
  double trueHeight(0);
  double boxwidth=max-min;
  double wgt(0.);

  for (int i=0;i<int(intervalMH+1.0);i++) {
    if(brHad[i]>boxheight)boxheight=brHad[i];
  }
  while ((mass > max) || (mass < min)){
    xbox = EvtRandom::Flat(boxwidth)+min;
    ybox=EvtRandom::Flat(boxheight);
    trueHeight=0.0;
    // Correction by Peter Richardson
    for( int i = 1 ; i < int( intervalMH + 1.0 ) ; ++i ) {
      if ( ( massHad[i] >= xbox ) && ( 0.0 == trueHeight ) ) {
	wgt=(xbox-massHad[i-1])/(massHad[i]-massHad[i-1]);
	trueHeight=brHad[i-1]+wgt*(brHad[i]-brHad[i-1]);
      }
    }

    if (ybox>trueHeight) {
      mass=0.0;
    } else {
      mass=xbox;
    }
  }
 
  return mass;
}

double EvtBtoXsgammaKagan::CalcAlphaS(double scale) {

  double v = 1. -_beta0*(_alphasmZ/(2.*EvtConst::pi))*(log(_mZ/scale));
  return (_alphasmZ/v)*(1. - ((_beta1/_beta0)*(_alphasmZ/(4.*EvtConst::pi))*(log(v)/v)));

}

void EvtBtoXsgammaKagan::CalcWilsonCoeffs( ){
  
   double mtatmw=_mt*pow((_alphasmW/_alphasmt),(12./23.))*(1 + (12./23.)*((253./18.) - (116./23.))*((_alphasmW - _alphasmt)/(4.0*EvtConst::pi)) - (4./3.)*(_alphasmt/EvtConst::pi));
  double xt=pow(mtatmw,2.)/pow(_mW,2.);
 

  
  /////LO
  _c2mu = .5*pow(_etamu,(-12./23.)) + .5*pow(_etamu,(6./23.));
  
  double c7mWsm = ((3.*pow(xt,3.) - 2.*pow(xt,2.))/(4.*pow((xt - 1.),4.)))*log(xt)
    + ((-8.*pow(xt,3.) - 5.*pow(xt,2.) + 7.*xt)/(24.*pow((xt - 1.),3.) )) ;
  
  double c8mWsm =  ((-3.*pow(xt,2.))/(4.*pow((xt - 1.),4.)))*log(xt)
    + ((- pow(xt,3.) + 5.*pow(xt,2.) + 2.*xt)/(8.*pow((xt - 1.),3.)));
  
  double c7constmu = (626126./272277.)*pow(_etamu,(14./23.))
    - (56281./51730.)*pow(_etamu,(16./23.)) - (3./7.)*pow(_etamu,(6./23.)) 
    - (1./14.)*pow(_etamu,(-12./23.)) - .6494*pow(_etamu,.4086) - .038*pow(_etamu,-.423) 
    - .0186*pow(_etamu,-.8994) - .0057*pow(_etamu,.1456);
  
  _c70mu = c7mWsm*pow(_etamu,(16./23.)) + (8./3.)*(pow(_etamu,(14./23.))
    -pow(_etamu,(16./23.)))*c8mWsm + c7constmu;
  
  double c8constmu =  (313063./363036.)*pow(_etamu,(14./23.))
    -.9135*pow(_etamu,.4086) + .0873*pow(_etamu,-.423) - .0571*pow(_etamu,-.8994)
    + .0209*pow(_etamu,.1456);

  _c80mu = c8mWsm*pow(_etamu,(14./23.)) + c8constmu;

 //Compute the dilogarithm (PolyLog(2,x)) with the Simpson integrator
 //The dilogarithm is defined as: Li_2(x)=Int_0^x(-log(1.-z)/z)
 //however, Mathematica implements it as  Sum[z^k/k^2,{k,1,Infinity}], so, althought the two
 //results are similar and both implemented in the program, we prefer to use the
 //one closer to the Mathematica implementation as identical to what used by the theorists.
  
 // EvtItgFunction *myDiLogFunc = new EvtItgFunction(&diLogFunc, 0., 1.-1./xt);
 //EvtItgAbsIntegrator *myDiLogSimp = new EvtItgSimpsonIntegrator(*myDiLogFunc, 1.0e-4, 50);
 //double li2 = myDiLogSimp->evaluate(1.0e-20,1.-1./xt);

 double li2=diLogMathematica(1.-1./xt);

double c7mWsm1 = ( (-16. *pow(xt,4.) -122. *pow(xt,3.) + 80. *pow(xt,2.) -8. *xt)/
(9. *pow((xt -1.),4.)) * li2 +
(6. *pow(xt,4.) + 46. *pow(xt,3.) -28. *pow(xt,2.))/(3. *pow((xt-1.),5.)) *pow(log(xt),2.)
+ (-102. *pow(xt,5.) -588. *pow(xt,4.) -2262. *pow(xt,3.) + 3244. *pow(xt,2.) -1364. *xt
+ 208.)/(81. *pow((xt-1),5.)) *log(xt)
+ (1646. *pow(xt,4.) + 12205. *pow(xt,3.) -10740. *pow(xt,2.) + 2509. *xt -436.)/
(486. *pow((xt-1),4.)) );

double c8mWsm1 = ((-4. *pow(xt,4.) + 40. *pow(xt,3.) + 41. *pow(xt,2.) + xt)/
(6. *pow((xt-1.),4.))  * li2
+ (-17. *pow(xt,3.) -31. *pow(xt,2.))/(2. *pow((xt-1.),5.) ) *pow(log(xt),2.)
+ (-210. *pow(xt,5.) + 1086. *pow(xt,4.) + 4893. *pow(xt,3.) + 2857. *pow(xt,2.)
-1994. *xt + 280.)/(216. *pow((xt-1),5.)) *log(xt)
+ (737. *pow(xt,4.) -14102. *pow(xt,3.) -28209. *pow(xt,2.) + 610. *xt -508.)/
(1296. *pow((xt-1),4.)) );

double E1 = (xt *(18. -11. *xt -pow(xt,2.))/(12.*pow( (1. -xt),3.))
+ pow(xt,2.)* (15. -16. *xt + 4. *pow(xt,2.))/(6. *pow((1. -xt),4.)) *log(xt)
-2./3. *log(xt) );

double e1 = 4661194./816831.;
double e2 = -8516./2217. ;
double e3 = 0.;
double e4 = 0.;
double e5 = -1.9043;
double e6 = -.1008;
double e7 = .1216;
double e8 = .0183;

double f1 = -17.3023;
double f2 = 8.5027;
double f3 = 4.5508;
double f4 = .7519;
double f5 =  2.004;
double f6 = .7476;
double f7 = -.5385;
double f8 = .0914;

double g1 = 14.8088;
double g2 = -10.809;
double g3 = -.874;
double g4 = .4218;
double g5 = -2.9347;
double g6 = .3971;
double g7 = .1600;
double g8 = .0225;


double c71constmu  = ((e1 *_etamu *E1 + f1 + g1 *_etamu) *pow(_etamu,(14./23.))
+ (e2 *_etamu *E1 + f2 + g2 *_etamu) *pow(_etamu,(16./23.))
+ (e3 *_etamu *E1 + f3 + g3 *_etamu) *pow(_etamu,(6./23.))
+ (e4 *_etamu *E1 + f4 + g4 *_etamu) *pow(_etamu,(-12./23.))
+ (e5 *_etamu *E1 + f5 + g5 *_etamu) *pow(_etamu,.4086)
+ (e6 *_etamu *E1 + f6 + g6 *_etamu) *pow(_etamu,(-.423))
+ (e7 *_etamu *E1 + f7 + g7 *_etamu) *pow(_etamu,(-.8994))
+ (e8 *_etamu *E1 + f8 + g8 *_etamu) *pow(_etamu,.1456 ));

double c71pmu = ( ((297664./14283. *pow(_etamu,(16./23.))
-7164416./357075. *pow(_etamu,(14./23.))
+ 256868./14283. *pow(_etamu,(37./23.)) - 6698884./357075. *pow(_etamu,(39./23.)))
*(c8mWsm))
+ 37208./4761. *(pow(_etamu,(39./23.)) - pow(_etamu,(16./23.))) *(c7mWsm)
+ c71constmu );

_c71mu = (_alphasmW/_alphasmu *(pow(_etamu,(16./23.))* c7mWsm1 + 8./3. *(pow(_etamu,(14./23.))
- pow(_etamu,(16./23.)) ) *c8mWsm1 ) + c71pmu);

_c7emmu = ((32./75. *pow(_etamu,(-9./23.)) - 40./69. *pow(_etamu,(-7./23.)) +
               88./575. *pow(_etamu,(16./23.))) *c7mWsm + (-32./575. *pow(_etamu,(-9./23.)) +
               32./1449. *pow(_etamu,(-7./23.)) + 640./1449.*pow(_etamu,(14./23.)) -
               704./1725.*pow(_etamu,(16./23.)) ) *c8mWsm
         - 190./8073.*pow(_etamu,(-35./23.))  - 359./3105. *pow(_etamu,(-17./23.)) +
         4276./121095. *pow(_etamu,(-12./23.)) + 350531./1009125.*pow(_etamu,(-9./23.))
         + 2./4347. *pow(_etamu,(-7./23.)) - 5956./15525. *pow(_etamu,(6./23.)) +
         38380./169533. *pow(_etamu,(14./23.))   - 748./8625. *pow(_etamu,(16./23.)));

// Wilson coefficients values as according to Kagan's program
// _c2mu=1.10566;
//_c70mu=-0.314292;
// _c80mu=-0.148954; 
// _c71mu=0.480964;
// _c7emmu=0.0323219;
 
}

void EvtBtoXsgammaKagan::CalcDelta() {

  double cDelta77 = (1. + (_alphasmu/(2.*EvtConst::pi)) *(_r7 - (16./3.) + _gam77*log(_mb/_mu)) + ( (pow((1. - _z),4.)/_fz) - 1.)*(6.*_lam2/pow(_mb,2.)) + (_alphasmubar/(2.*EvtConst::pi))*_kappabar )*pow(_c70mu,2.);
  
  double cDelta27 = ((_alphasmu/(2.*EvtConst::pi))*(_rer2 + _gam27*log(_mb/_mu)) - (_lam2/(9.*_z*pow(_mb,2.))))*_c2mu*_c70mu;
  
  double cDelta78 = (_alphasmu/(2.*EvtConst::pi))*(_rer8 + _gam87*log(_mb/_mu))*_c70mu*_c80mu;
  
  _cDeltatot = cDelta77  + cDelta27 + cDelta78 + (_alphasmu/(2.*EvtConst::pi))*_c71mu*_c70mu + (_alpha/_alphasmu)*(2.*_c7emmu*_c70mu - _kSLemmu*pow(_c70mu,2.));
  
}

double EvtBtoXsgammaKagan::Delta(double y, double alphasMu) {
  
  //Fix for singularity at endpoint
  if (y >= 1.0) y = 0.9999999999;
  
  return ( - 4.*(alphasMu/(3.*EvtConst::pi*(1. - y)))*(log(1. - y) + 7./4.)*
	   exp(-2.*(alphasMu/(3.*EvtConst::pi))*(pow(log(1. - y),2) + (7./2.)*log(1. - y))));
  
}

double EvtBtoXsgammaKagan::s77(double y) {
  
  //Fix for singularity at endpoint
  if (y >= 1.0) y = 0.9999999999;
  
  return ((1./3.)*(7. + y - 2.*pow(y,2) - 2.*(1. + y)*log(1. - y)));
}

double EvtBtoXsgammaKagan::s88(double y, double mb, double ms) {
  
  //Fix for singularity at endpoint
  if (y >= 1.0) y = 0.9999999999;
  
  return ((1./27.)*((2.*(2. - 2.*y + pow(y,2))/y)*(log(1. - y) + 2.*log(mb/ms))
		    - 2.*pow(y,2) - y - 8.*((1. - y)/y)));
}

double EvtBtoXsgammaKagan::s78(double y) {
  
  //Fix for singularity at endpoint
  if (y >= 1.0) y = 0.9999999999;
  
  return ((8./9.)*(((1. - y)/y)*log(1. - y) + 1. + (pow(y,2)/4.)));
}

double EvtBtoXsgammaKagan::ReG(double y) {
  
  if (y < 4.) return -2.*pow(atan(sqrt(y/(4. - y))),2.);
  else {
    return 2.*(pow(log((sqrt(y) + sqrt(y - 4.))/2.),2.)) - (1./2.)*pow(EvtConst::pi,2.);
  }
  
}

double EvtBtoXsgammaKagan::ImG(double y) {
  
  if (y < 4.) return 0.0;
  else {
    return (-2.*EvtConst::pi*log((sqrt(y) + sqrt(y - 4.))/2.));
  }
}

double EvtBtoXsgammaKagan::s22Func(double y, const std::vector<double> &coeffs) {

  //coeffs[0]=z
  return (1. - y)*((pow(coeffs[0],2.)/pow(y,2.))*(pow(ReG(y/coeffs[0]),2.) + pow(ImG(y/coeffs[0]),2.)) + (coeffs[0]/y)*ReG(y/coeffs[0]) + (1./4.));
  
}

double EvtBtoXsgammaKagan::s27Func(double y, const std::vector<double> &coeffs) {
 
  //coeffs[0] = z
  return (ReG(y/coeffs[0]) + y/(2.*coeffs[0]));

}

double EvtBtoXsgammaKagan::DeltaFermiFunc(double y, const std::vector<double> &coeffs1, 
					const std::vector<double> &coeffs2, const std::vector<double> &coeffs3) {
 
  //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb, 
  //coeffs2[2]=ymH, coeffs3[0]=DeltaCoeff (alphasmu)
 
  return FermiFunc(y,coeffs1)*(coeffs2[0]/(coeffs2[1]+y))*
    Delta((coeffs2[0]*coeffs2[2])/(coeffs2[1]+y),coeffs3[0]);

}

double EvtBtoXsgammaKagan::s77FermiFunc(double y, const std::vector<double> &coeffs1, 
				      const std::vector<double> &coeffs2) {

  //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb, 
  //coeffs2[2]=ymH
  return FermiFunc(y,coeffs1)*(coeffs2[0]/(coeffs2[1]+y))*
    s77((coeffs2[0]*coeffs2[2])/(coeffs2[1]+y));

}

double EvtBtoXsgammaKagan::s88FermiFunc(double y, const std::vector<double> &coeffs1,  
				      const std::vector<double> &coeffs2, const std::vector<double> &coeffs3) {

  //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb, 
  //coeffs2[2]=ymH, coeffs3=s88 coeffs
  return FermiFunc(y,coeffs1)*(coeffs2[0]/(coeffs2[1]+y))*
   s88((coeffs2[0]*coeffs2[2])/(coeffs2[1]+y),coeffs3[0], coeffs3[1]);

}

double EvtBtoXsgammaKagan::s78FermiFunc(double y, const std::vector<double> &coeffs1, 
				      const std::vector<double> &coeffs2) {

  //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb, 
  //coeffs2[2]=ymH
  return FermiFunc(y,coeffs1)*(coeffs2[0]/(coeffs2[1]+y))*
    s78((coeffs2[0]*coeffs2[2])/(coeffs2[1]+y));

}

double EvtBtoXsgammaKagan::sFermiFunc(double y, const std::vector<double> &coeffs1, 
				      const std::vector<double> &coeffs2, const std::vector<double> &coeffs3, 
				      const std::vector<double> &coeffs4) {

  //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb, 
  //coeffs2[2]=ymH, coeffs3[0]=nIntervals in s22 or s27 array, coeffs3[1]=yMin,
  //coeffs3[2]=yMax, coeffs4=s22 or s27 array
  return FermiFunc(y,coeffs1)*(coeffs2[0]/(coeffs2[1]+y))*
    GetArrayVal(coeffs2[0]*coeffs2[2]/(coeffs2[1]+y), coeffs3[0], coeffs3[1], coeffs3[2], coeffs4);

}

double EvtBtoXsgammaKagan::Fz(double z) {

  return (1. -8.*z + 8.*pow(z,3.) - pow(z,4.) - 12.*pow(z,2.)*log(z));
}

double EvtBtoXsgammaKagan::GetArrayVal(double xp, double nInterval, double xMin, double xMax, std::vector<double> array) {
 
  double dx = (xMax - xMin)/nInterval;
  int bin1 = int(((xp-xMin)/(xMax - xMin))*nInterval);

  double x1 = double(bin1)*dx + xMin;

  if (xp == x1) return array[bin1];

  int bin2(0);
  if (xp > x1) {
    bin2 = bin1 + 1;
  }
  else if (xp < x1) {
    bin2 = bin1 - 1;
  }

  if (bin1 <= 0) {
    bin1=0;
    bin2 = 1;
  }

  //If xp is in the last bin, always interpolate between the last two bins
  if (bin1 == (int)nInterval){
    bin2 = (int)nInterval;
    bin1 = (int)nInterval - 1;
    x1 = double(bin1)*dx + xMin;
  }
 
  double x2 = double(bin2)*dx + xMin;
  double y1 = array[bin1];

  double y2 = array[bin2];
  double m = (y2 - y1)/(x2 - x1);
  double c =  y1 - m*x1;
  double result = m*xp + c;

  return result;
  
}

double EvtBtoXsgammaKagan::FermiFunc(double y, const std::vector<double> &coeffs) {

  //Fermi shape functions :1=exponential, 2=gaussian, 3=roman
  if (int(coeffs[0]) == 1) return EvtBtoXsgammaFermiUtil::FermiExpFunc(y, coeffs);
  if (int(coeffs[0]) == 2) return EvtBtoXsgammaFermiUtil::FermiGaussFunc(y, coeffs);
  if (int(coeffs[0]) == 3) return EvtBtoXsgammaFermiUtil::FermiRomanFunc(y, coeffs);
  return 1.;

}

double EvtBtoXsgammaKagan::diLogFunc(double y) {

  return -log(fabs(1. - y))/y;

}


double EvtBtoXsgammaKagan::diLogMathematica(double y) {

  double li2(0);
  for(int i=1; i<1000; i++){ //the value 1000 should actually be Infinite...
    li2+=pow(y,i)/(i*i);
  }
  return li2;
}
