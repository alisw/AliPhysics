//----------------------------------------------------------------------
//                              AliBWFunc
//
// This class implements several function useful to fit pt spectra,
// including but not limited to blast wave models.
//
// It can return the same functional for as a function of different
// variables: dNdpt vs pt, 1/pt dNdpt vs pt, 1/mt dNdmt vs mt. 
//
// Before getting the function you need, you have to chose the
// variable you want to use calling AliBWFunc::SetVarType with one of
// the elements of the VarType_t enum.
//
// Warning: not all variables are implemented for all the functions.
//
// Author: M. Floris, CERN
//----------------------------------------------------------------------

#include "AliBWFunc.h"
#include "TMath.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1.h"
#include "TSpline.h"
#include "AliLog.h"

ClassImp(AliBWFunc)

AliBWFunc::AliBWFunc () : fLastFunc(0),  fLineWidth(1), fVarType(kdNdpt) {

  // ctor
  fLineWidth = 1;
}
AliBWFunc::~AliBWFunc(){
  
  // dtor
  if (fLastFunc) delete fLastFunc;

}


TF1 * AliBWFunc::GetHistoFunc(TH1 * h, const char * name) {

  // Regardless of the variable type, this returns a function made
  // from the histo * a multiplicative normalization.
  // This uses a bad hack...

  fLastFunc = new TF1 (name, StaticHistoFunc, 0.0, 10, 2);
  fLastFunc->SetParameter(0,1);
  fLastFunc->FixParameter(1,Double_t(Long64_t(h)));
  fLastFunc->SetParNames("norm", "pointer to histo");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
  


}
TF1 * AliBWFunc::GetGraphFunc(TGraph * g, const char * name) {

  // Regardless of the variable type, this returns a function made
  // from the graph * a multiplicative normalization.
  // This uses a bad hack...

  fLastFunc = new TF1 (name, StaticHistoFunc, 0.0, 10, 2);
  fLastFunc->SetParameter(0,1);
  fLastFunc->FixParameter(1,Double_t(Long64_t(g)));
  fLastFunc->SetParNames("norm", "pointer to histo");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
  

}


TF1 * AliBWFunc::GetBGBW(Double_t mass, Double_t beta, Double_t T,
			 Double_t n, Double_t norm, const char * name){

  // Boltzmann-Gibbs blast wave

  switch (fVarType) {
  case kdNdpt:
    return GetBGBWdNdptTimesPt(mass,beta,T,n,norm,name);
    break;
  case kOneOverPtdNdpt:
    return GetBGBWdNdpt(mass,beta,T,n,norm,name);
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not implemented");
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;

}
  

TF1 * AliBWFunc::GetBoltzmann(Double_t mass, Double_t T, Double_t norm, const char * name){
  // Boltzmann
  switch (fVarType) {
  case kdNdpt:
    return GetBoltzmanndNdptTimesPt(mass, T, norm, name);
  case kOneOverPtdNdpt:
    AliFatal("Not implemented");
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not implemented");
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;

}


TF1 * AliBWFunc::GetTsallisBW(Double_t mass, Double_t beta, Double_t T, Double_t q,
			      Double_t norm, Double_t ymax, const char * name){

  // Tsallis blast wave
  switch (fVarType) {
  case kdNdpt:
    return GetTsallisBWdNdptTimesPt(mass,beta,T,q,norm,ymax,name);
    break;
  case kOneOverPtdNdpt:
    return GetTsallisBWdNdpt(mass,beta,T,q,norm,ymax,name);
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not implemented");
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;
  
}


TF1 * AliBWFunc::GetMTExp(Double_t mass, Double_t T, Double_t norm, const char * name){

  // Simple exponential in 1/mt*MT
  switch (fVarType) {
  case kdNdpt:
    return GetMTExpdNdptTimesPt(mass,T,norm,name);
    break;
  case kOneOverPtdNdpt:
    return GetMTExpdNdpt(mass,T,norm,name);
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not implemented");
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;


}


TF1 * AliBWFunc::GetBoseEinstein(Double_t mass, Double_t T, Double_t norm, const char * name){

  // Bose einstein
  switch (fVarType) {
  case kdNdpt:
    return GetBoseEinsteindNdptTimesPt(mass,T,norm,name);
    break;
  case kOneOverPtdNdpt:
    return GetBoseEinsteindNdpt(mass,T,norm,name);
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not implemented");
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;


}

TF1 * AliBWFunc::GetFermiDirac(Double_t mass, Double_t T, Double_t norm, const char * name){

  // Simple exponential in 1/mt*MT
  switch (fVarType) {
  case kdNdpt:
    return GetFermiDiracdNdptTimesPt(mass,T,norm,name);
    break;
  case kOneOverPtdNdpt:
    return GetFermiDiracdNdpt(mass,T,norm,name);
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not implemented");
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;


}


TF1 * AliBWFunc::GetPTExp(Double_t T, Double_t norm, const char * name){

  // Simple exponential in 1/mt*MT
  switch (fVarType) {
  case kdNdpt:
    return GetPTExpdNdptTimesPt(T,norm,name);
    break;
  case kOneOverPtdNdpt:
    AliFatal("Not implemented");
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not implemented");
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;


}


TF1 * AliBWFunc::GetLevi(Double_t mass, Double_t T, Double_t n, Double_t norm, const char * name){
  // Levi function (aka Tsallis)
  switch (fVarType) {
  case kdNdpt:
    return GetLevidNdptTimesPt(mass,T,n,norm,name);
    break;
  case kOneOverPtdNdpt:
    return GetLevidNdpt(mass,T,n,norm,name);
    break;
  case kOneOverMtdNdmt:
    return GetLevidNdmt(mass,T,n,norm,name,kOneOverMtdNdmt);
    break;
  case kdNdmt:
    return GetLevidNdmt(mass,T,n,norm,name,kdNdmt);
    break;
  case kOneOverMtdNdmtMinusM:
    return GetLevidNdmt(mass,T,n,norm,name,kOneOverMtdNdmtMinusM);
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;


}

TF1 * AliBWFunc::GetPowerLaw(Double_t pt0, Double_t n, Double_t norm, const char * name){
  // power law Nuclear Physics B, Vol. 335, No. 2. (7 May 1990), pp. 261-287.
  // This is sometimes also called Hagedorn or modified Hagedorn

  switch (fVarType) {
  case kdNdpt:
    return GetPowerLawdNdptTimesPt(pt0,n,norm,name);
    break;
  case kOneOverPtdNdpt:
    return GetPowerLawdNdpt(pt0,n,norm,name);
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not Implemented");
    //    return GetUA1dNdmt(mass,T,n,norm,name);
    break;
  default:
    AliFatal("Not implemented");
  }
  
  return 0;


}

TF1 * AliBWFunc::GetUA1(Double_t mass, Double_t p0star, Double_t pt0, Double_t n, Double_t T, Double_t norm, const char * name) {
  // UA1 parametrization Nuclear Physics B, Vol. 335, No. 2. (7 May 1990), pp. 261-287.

  switch (fVarType) {
  case kdNdpt:

    fLastFunc = new TF1 (name, StaticUA1Func, 0.0, 10, 6);
    fLastFunc->FixParameter(0,mass);
    fLastFunc->SetParameter(1,p0star);
    fLastFunc->SetParameter(2,pt0);
    fLastFunc->SetParameter(3,n);
    fLastFunc->SetParameter(4,T);
    fLastFunc->SetParameter(5,norm);
    fLastFunc->SetParLimits(1,0.01,1);
    fLastFunc->SetParLimits(2,0.01,100);
    fLastFunc->SetParLimits(3,0.01,100);
    fLastFunc->SetParLimits(4,0.01,100);
    fLastFunc->SetParNames("mass","p0star","pt0","n","T","norm");
    fLastFunc->SetNpx(5000);
    fLastFunc->SetLineWidth(fLineWidth);
    return fLastFunc;

    break;
  case kOneOverPtdNdpt:
    AliFatal("Not Implemented");
    break;
  case kOneOverMtdNdmt:
    AliFatal("Not Implemented");
    //    return GetUA1dNdmt(mass,T,n,norm,name);
    break;
  default:
    AliFatal("Not implemented");
  }

  return 0;
}




// ________________________________________________________________________

// Backend (private functions and support functions for numerical integration)

Double_t AliBWFunc::StaticHistoFunc(const double * x, const double* p){

  // provides a function interpolating a histo with a spline; 
  // using double to store a pointer... This is a bad hack. To be replaced

  double norm = p[0];
  
  TObject * h     = (TObject*) Long64_t(p[1]);

//    Int_t bin = h->FindBin(x[0]);
//    double value = h->GetBinContent(bin);


  // static TH1 * oldptr = 0;
  // static TSpline3 * spl = 0;
  // if (h!=oldptr) {
  // FIXME: recheck static pointers
  TSpline3 * spl  = 0;
  if(h->InheritsFrom("TH1")) {
    if ( ((TH1*)h)->FindBin(x[0]) > ((TH1*)h)->GetNbinsX()) return 0;
    spl= new TSpline3((TH1*)h);
  }
  else if(h->InheritsFrom("TGraph")) spl= new TSpline3("fGraph",(TGraph*)h);
  else {
    Printf("AliBWFunc::StaticHistoFunc: Unsupported type");
    return 0;
  }
    //  }
  double value = spl->Eval(x[0]);
  delete spl;

  return value*norm;
  
}

Double_t AliBWFunc::StaticUA1Func(const double * x, const double* p) {
  

  // "mass","p0star","pt0","n","T","norm"
  Double_t mass   = p[0];
  Double_t p0star = p[1];
  Double_t pt0    = p[2];
  Double_t n      = p[3];
  Double_t temp   = p[4];
  Double_t norm   = p[5];
  
  Double_t xx = x[0];

  static AliBWFunc * self = new AliBWFunc;
  static TF1 * fPLaw   = self->GetPowerLawdNdptTimesPt(pt0, n, norm, "fLocalPLawUA1");
  static TF1 * fPMTExp = self->GetMTExpdNdptTimesPt   (mass, temp, norm, "fLocalMTexpUA1");

  fPLaw->SetParameters(norm,pt0,n);
  fPMTExp->SetParameters(1,temp);
  

  Double_t normMT =fPMTExp->Eval(p0star) >0 ? fPLaw->Eval(p0star) / fPMTExp->Eval(p0star) *  fPMTExp->GetParameter(0) : 1;
  fPMTExp->SetParameter(0,normMT);
  
  
  if (TMath::Abs(fPMTExp->Eval(p0star) - fPLaw->Eval(p0star)) > 0.0001 ) {
    Printf("AliBWFunc::StaticUA1Func - Wrong norm") ; 
    Printf(" p0* %f  NMT: %f  N: %f  PL: %f  MT: %f", p0star, normMT, norm, fPLaw->Eval(p0star), fPMTExp->Eval(p0star));
  }

  if (xx > p0star)  return fPLaw->Eval(xx);
  return fPMTExp->Eval(xx);    
  
  

}


Double_t AliBWFunc::IntegrandBG(const double * x, const double* p){
  // integrand for boltzman-gibbs blast wave
     // x[0] -> r (radius)
     // p[0] -> mass
     // p[1] -> pT (transverse momentum)
     // p[2] -> beta_max (surface velocity)
     // p[3] -> T (freezout temperature)
     // p[4] -> n (velocity profile)


  double x0 = x[0]; 
  
  double mass     = p[0];
  double pT       = p[1];
  double beta_max = p[2];
  double temp     = p[3];
  Double_t n      = p[4];

  // Keep beta within reasonable limits
  Double_t beta = beta_max * TMath::Power(x0, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;

  double mT      = TMath::Sqrt(mass*mass+pT*pT);

  double rho0   = TMath::ATanH(beta);  
  double arg00 = pT*TMath::SinH(rho0)/temp;
  if (arg00 > 700.) arg00 = 700.; // avoid FPE
  double arg01 = mT*TMath::CosH(rho0)/temp;
  double f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);

  //  printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", x0, pT, beta_max, temp, n, mT, beta, rho0, arg00, arg01);

  return f0;
}



Double_t AliBWFunc::StaticBGdNdPt(const double * x, const double* p) {

  // implementation of BGBW (1/pt dNdpt)

  double pT = x[0];;
  

  double mass    = p[0];
  double beta    = p[1];
  double temp    = p[2];
  double n       = p[3];
  double norm    = p[4];

  static TF1 * fIntBG = 0;
  if(!fIntBG)
    fIntBG = new TF1 ("fIntBG", IntegrandBG, 0, 1, 5);

  fIntBG->SetParameters(mass, pT, beta, temp,n);
  double result = fIntBG->Integral(0,1);
  //  printf ("[%4.4f], Int :%f\n", pT, result);
  return result*norm;//*1e30;;

}

Double_t AliBWFunc::StaticBGdNdPtTimesPt(const double * x, const double* p) {
  // BGBW dNdpt implementation
  return x[0]*StaticBGdNdPt(x,p);
}


TF1 * AliBWFunc::GetBGBWdNdpt(Double_t mass, Double_t beta, Double_t temp,
			      Double_t n, Double_t norm, const char * name){
  
  // BGBW 1/pt dNdpt

  fLastFunc = new TF1 (name, StaticBGdNdPt, 0.0, 10, 5);
  fLastFunc->SetParameters(mass,beta,temp,n,norm);    
  fLastFunc->FixParameter(0,mass);
  fLastFunc->SetParNames("mass", "#beta", "T", "n", "norm");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
  
}


//_____________________________________________________________________
// Tsallis

Double_t AliBWFunc::IntegrandTsallis(const double * x, const double* p){

  // integrand for numerical integration (tsallis)

  Double_t r   = x[0]; 
  Double_t phi = x[1];
  Double_t y   = x[2];

  Double_t mass = p[0];
  Double_t pt   = p[1];
  Double_t beta = p[2];
  Double_t temp    = p[3];
  Double_t q    = p[4];
  
  Double_t mt      = TMath::Sqrt(mass*mass+pt*pt);

  Double_t rho    = TMath::ATanH(beta*r); // TODO: implement different velocity profiles  

  Double_t res = mt*
    r*TMath::CosH(y) *TMath::Power( (
				     1+(q-1)/temp * (
						  mt*TMath::CosH(y)*TMath::CosH(rho) -
						  pt*TMath::SinH(rho)*TMath::Cos(phi)
						  )
				     ),
				       -1/(q-1)
				    );			


  return res;
}



Double_t AliBWFunc::StaticTsallisdNdPt(const double * x, const double* p) {

  // tsallis BW implementation 1/pt dNdpt

  double pT = x[0];;
  

  double mass = p[0];
  double beta = p[1];
  double temp    = p[2];
  double q    = p[3];

  Double_t ymax = p[5];


  static TF3 * fInt = 0;
  if(!fInt){
    fInt = new TF3 ("fIntTsa", IntegrandTsallis, 0, 1, -TMath::Pi(), TMath::Pi(), -ymax, ymax, 5);
//     fInt->SetNpx(10000);
//     fInt->SetNpy(10000);
//     fInt->SetNpz(10000);
  }
  
  fInt->SetParameters(mass, pT, beta, temp, q);
  double result = fInt->Integral(0,1, -TMath::Pi(), TMath::Pi(), -ymax, ymax);
  //  double result = fInt->Integral(0,1, -2, 2, -ymax, ymax);
  
  return result*p[4];//*1e30;;

}

Double_t AliBWFunc::StaticTsallisdNdPtTimesPt(const double * x, const double* p) {

  // tsallis BW , implementatio of dNdpt
  return x[0]*StaticTsallisdNdPt(x,p);

}

TF1 * AliBWFunc::GetTsallisBWdNdpt(Double_t mass, Double_t beta, Double_t temp, Double_t q,
				   Double_t norm, Double_t ymax,const char * name){
  

  // tsallis BW, 1/pt dNdpt

  fLastFunc = new TF1 (name, StaticTsallisdNdPt, 0.0, 10, 6);
  fLastFunc->SetParameters(mass,beta,temp,q,norm,ymax);
  fLastFunc->SetParLimits(1,0.0,0.99);
  fLastFunc->SetParLimits(2,0.01,0.99);
  fLastFunc->SetParLimits(3,1.0001,1.9);
  fLastFunc->SetParNames("mass", "#beta", "temp", "q", "norm", "ymax");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
  
}

// Times Pt funcs
// Boltzmann-Gibbs Blast Wave
TF1 * AliBWFunc::GetBGBWdNdptTimesPt(Double_t mass, Double_t beta, Double_t temp, Double_t n,
				     Double_t norm, const char * name){

  // BGBW, dNdpt

  fLastFunc = new TF1 (name, StaticBGdNdPtTimesPt, 0.0, 10, 5);
  fLastFunc->SetParameters(mass,beta,temp,n,norm);    
  fLastFunc->FixParameter(0,mass);
  fLastFunc->SetParNames("mass", "#beta", "temp", "n", "norm");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}



TF1 * AliBWFunc::GetTsallisBWdNdptTimesPt(Double_t mass, Double_t beta, Double_t temp, Double_t q,
					  Double_t norm, Double_t ymax, const char * name){

// Tsallis blast wave, dNdpt

  fLastFunc = new TF1 (name, StaticTsallisdNdPtTimesPt, 0.0, 10, 6);
  fLastFunc->SetParameters(mass,beta,temp,q,norm,ymax);    
  fLastFunc->SetParNames("mass", "#beta", "temp", "q", "norm", "ymax");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
 

}



TF1 * AliBWFunc::GetMTExpdNdptTimesPt(Double_t mass, Double_t temp, Double_t norm, const char * name){

  // Simple exponential in 1/mt*MT, as a function of dNdpt
  char formula[500];
  snprintf(formula,500,"[0]*x*exp(-sqrt(x**2+%f**2)/[1])", mass);
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, temp);
  fLastFunc->SetParLimits(1, 0.01, 10);
  fLastFunc->SetParNames("norm", "T");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}

TF1 * AliBWFunc::GetBoseEinsteindNdptTimesPt(Double_t mass, Double_t temp, Double_t norm, const char * name){

  // Bose einstein distribution as a function of dNdpt
  char formula[500];
  snprintf(formula,500,"[0]*x*1./(exp(sqrt(x**2+%f**2)/[1])-1)", mass);
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, temp);
  fLastFunc->SetParLimits(1, 0.01, 10);
  fLastFunc->SetParNames("norm", "T");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}

TF1 * AliBWFunc::GetFermiDiracdNdptTimesPt(Double_t mass, Double_t temp, Double_t norm, const char * name){

  // Bose einstein distribution as a function of dNdpt
  char formula[500];
  snprintf(formula,500,"[0]*x*1./(exp(sqrt(x**2+%f**2)/[1])+1)", mass);
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, temp);
  fLastFunc->SetParLimits(1, 0.01, 10);
  fLastFunc->SetParNames("norm", "T");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}



TF1 * AliBWFunc::GetPTExpdNdptTimesPt(Double_t temp, Double_t norm, const char * name){

  // Simple exponential in 1/pt*dNdpT, as a function of dNdpt
  char formula[500];
  snprintf(formula,500,"[0]*x*exp(-x/[1])");
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, temp);
  fLastFunc->SetParLimits(1, 0.01, 10);
  fLastFunc->SetParNames("norm", "T");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}


TF1 * AliBWFunc::GetBoltzmanndNdptTimesPt(Double_t mass, Double_t temp, Double_t norm, const char * name){
  // Boltzmann (exp in 1/mt*dNdmT times mt) as a function of dNdpt
 char formula[500];
 snprintf(formula,500,"[0]*x*sqrt(x**2+%f**2)*exp(-sqrt(x**2+%f**2)/[1])", mass,mass);
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, temp);
  fLastFunc->SetParLimits(1, 0.01, 10);
  fLastFunc->SetParNames("norm", "T");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}


// Tsallis (no BW, a la CMS)
// TF1 * AliBWFunc::GetTsallisdNdptTimesPt(Double_t mass, Double_t T, Double_t q, Double_t norm, const char * name){

//   char formula[500];
//   //  sprintf(formula,"[0]*x*pow((1+(([2]-1)/[1])*(sqrt(x**2+%f**2)-%f)),(-1/([2]-1)))", mass, mass); //CMS
//   sprintf(formula,"[0]*x*pow((1+(([2]-1)/[1])*(sqrt(x**2+%f**2))),(-1/([2]-1)))", mass);  // STAR
//   //sprintf(formula,"[0]*x*sqrt(x**2+%f**2)*pow((1+(([2]-1)/[1])*(sqrt(x**2+%f**2))),(-1/([2]-1)))", mass,mass);  // STAR * mt
//   fLastFunc=new TF1(name,formula,0,10);
//   fLastFunc->SetParameters(norm, T, q);
//   fLastFunc->SetParLimits(1, 0.001, 10);
//   fLastFunc->SetParNames("norm", "T", "q");
//   fLastFunc->SetLineWidth(fLineWidth);
//   return fLastFunc;


// }


TF1 * AliBWFunc::GetLevidNdptTimesPt(Double_t mass, Double_t temp, Double_t n, Double_t norm, const char * name){

  // Levi function, dNdpt
  char formula[500];

  snprintf(formula,500,"( x*[0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + (sqrt([3]*[3]+x*x) -[3])/([1]*[2])  )^(-[1])");
  //  sprintf(formula,"( x*[0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + (sqrt([3]*[3]+x*x))/([1]*[2])  )^(-[1])");
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, n, temp,mass);
  fLastFunc->SetParLimits(2, 0.01, 10);
  fLastFunc->SetParNames("norm (dN/dy)", "n", "T", "mass");
  fLastFunc->FixParameter(3,mass);
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}

TF1 * AliBWFunc::GetPowerLawdNdptTimesPt(Double_t pt0, Double_t n, Double_t norm, const char * name){

  // PowerLaw function, dNdpt
  char formula[500];

  snprintf(formula,500,"x*[0]*( 1 + x/[1] )^(-[2])");
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, pt0, n);
  fLastFunc->SetParLimits(1, 0.01, 10);
  //fLastFunc->SetParLimits(2, 0.01, 50);
  fLastFunc->SetParNames("norm", "pt0", "n");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}

TF1 * AliBWFunc::GetPowerLawdNdpt(Double_t pt0, Double_t n, Double_t norm, const char * name){

  // PowerLaw function, 1/pt dNdpt
  char formula[500];

  snprintf(formula,500," [0]*( 1 + x/[1] )^(-[2])");
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, pt0, n);
  //  fLastFunc->SetParLimits(2, 0.01, 10);
  fLastFunc->SetParNames("norm", "pt0", "n");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}


TF1 * AliBWFunc::GetLevidNdpt(Double_t mass, Double_t temp, Double_t n, Double_t norm, const char * name){

  // Levi function, dNdpt
  char formula[500];

  snprintf(formula,500,"( [0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + (sqrt([3]*[3]+x*x) -[3])/([1]*[2])  )^(-[1])");
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, n, temp,mass);
  fLastFunc->SetParLimits(2, 0.01, 10);
  fLastFunc->SetParNames("norm (dN/dy)", "n", "T", "mass");
  fLastFunc->FixParameter(3,mass);
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}

TF1 * AliBWFunc::GetLevidNdmt(Double_t mass, Double_t temp, Double_t n, Double_t norm, const char * name, VarType_t var){

  // Levi function, 1/mt dNdmt
  char formula[500];
  if (var == kOneOverMtdNdmt)
    snprintf(formula,500,"( [0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + (x -[3])/([1]*[2])  )^(-[1])");
  else if (var == kdNdmt) 
    snprintf(formula,500,"( x*[0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + (x-[3])/([1]*[2])  )^(-[1])");
  if (var == kOneOverMtdNdmtMinusM)
    snprintf(formula,500,"( [0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + (x)/([1]*[2])  )^(-[1])");

  //sprintf(formula,"( [0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + x/([1]*[2])  )^(-[1])");
  //  sprintf(formula,"[0] * ( 1 + x/([1]*[2])  )^(-[1])");
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, n, temp,mass);
  fLastFunc->SetParLimits(2, 0.01, 10);
  fLastFunc->SetParNames("norm", "n", "T", "mass");
  fLastFunc->FixParameter(3,mass);
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;


}




// Test Function
Double_t AliBWFunc::IntegrandTest(const double * x, const double* p){

  // test function

  Double_t y = x[0];

  Double_t mass = p[0];
  Double_t pt   = p[1];
  Double_t temp    = p[2];

  Double_t mt      = TMath::Sqrt(mass*mass+pt*pt);    
  
  return mt*TMath::CosH(y)*TMath::Exp(-mt*TMath::CosH(y)/temp);

}

Double_t AliBWFunc::StaticTest(const double * x, const double* p) {

  // test function

  double pT = x[0];;
  

  double mass = p[0];
  double temp    = p[1];
  Double_t ymax = p[3];


  static TF3 * fIntTest = 0;
  if(!fIntTest){
    fIntTest = new TF3 ("fIntTest", IntegrandTest, 0, 1, -TMath::Pi(), TMath::Pi(), -ymax, ymax, 5);
    //    fInt->SetNpx(10000);
  }
  
  fIntTest->SetParameters(mass, pT, temp);
  double result = fIntTest->Integral(-ymax, ymax);
  
  return result*p[2];//*1e30;;

}

TF1 * AliBWFunc::GetTestFunc(Double_t mass, Double_t temp, Double_t norm, Double_t ymax, const char * name){
  
  // test function
  
  fLastFunc = new TF1 (name, StaticTest, 0.0, 10, 4);
  fLastFunc->SetParameters(mass,temp,norm,ymax);    
  fLastFunc->SetParNames("mass", "#beta", "T", "q", "norm", "ymax");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
  
}


//___________________________________________________________


TF1 * AliBWFunc::GetMTExpdNdpt(Double_t mass, Double_t temp, Double_t norm, const char * name){
  // Simple exp in 1/mt dNdmt, as a function of dNdpt
  // mt scaling
  char formula[500];
  snprintf(formula,500,"[0]*exp(-sqrt(x**2+%f**2)/[1])", mass);
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, temp);
  fLastFunc->SetParLimits(1, 0.01, 10);
  fLastFunc->SetParNames("norm", "T");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
}

TF1 * AliBWFunc::GetBoseEinsteindNdpt(Double_t mass, Double_t temp, Double_t norm, const char * name){
  // bose einstein
  char formula[500];
  snprintf(formula,500,"[0]*1./(exp(sqrt(x**2+%f**2)/[1])-1)", mass);
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, temp);
  fLastFunc->SetParLimits(1, 0.01, 10);
  fLastFunc->SetParNames("norm", "T");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
}

TF1 * AliBWFunc::GetFermiDiracdNdpt(Double_t mass, Double_t temp, Double_t norm, const char * name){
  // bose einstein
  char formula[500];
  snprintf(formula,500,"[0]*1./(exp(sqrt(x**2+%f**2)/[1])+1)", mass);
  fLastFunc=new TF1(name,formula,0,10);
  fLastFunc->SetParameters(norm, temp);
  fLastFunc->SetParLimits(1, 0.01, 10);
  fLastFunc->SetParNames("norm", "T");
  fLastFunc->SetLineWidth(fLineWidth);
  return fLastFunc;
}


// // Simple tsallis (a la CMS)
// TF1 * AliBWFunc::GetTsallisdNdpt(Double_t mass, Double_t temp, Double_t q, Double_t norm, const char * name){
  
//   char formula[500];
//   sprintf(formula,"[0]*sqrt(x**2+%f**2)*pow((1+(([2]-1)/[1])*(sqrt(x**2+%f**2))),(-1/([2]-1)))", mass,mass); 
//   fLastFunc=new TF1(name,formula,0,10);
//   fLastFunc->SetParameters(norm, temp, q);
//   fLastFunc->SetParLimits(1, 0.01, 10);
//   fLastFunc->SetParNames("norm", "T", "q");
//   fLastFunc->SetLineWidth(fLineWidth);
//   return fLastFunc;
// }
