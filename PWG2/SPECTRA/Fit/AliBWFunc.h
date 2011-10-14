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
// Author: M. Floris, CERN
//----------------------------------------------------------------------

#ifndef ALIBWFUNC_H
#define ALIBWFUNC_H

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TObject.h"

class TF1;
class TH1;
class TGraph;
#endif



class AliBWFunc : public TObject {


public:
  // define the variables used for the function
  typedef enum {kdNdpt,kOneOverPtdNdpt,kOneOverMtdNdmt,kdNdmt,kOneOverMtdNdmtMinusM} VarType_t;

  AliBWFunc();
  ~AliBWFunc();

  // Boltzmann-Gibbs Blast Wave
  TF1 * GetBGBW(Double_t mass, Double_t beta, Double_t T, Double_t n,
		Double_t norm, const char * name = "fBGBW");

  
  // Boltzmann
  TF1 * GetBoltzmann(Double_t mass, Double_t T, Double_t norm, const char * name ="fBoltzmann");

  // Tsallis blast wave
  TF1 * GetTsallisBW(Double_t mass, Double_t beta, Double_t T, Double_t q,
		     Double_t norm, Double_t ymax = 0.5, const char * name = "fTsallisBW");

  // Simple exponential in 1/mt*dNdmt
  TF1 * GetMTExp(Double_t mass, Double_t T, Double_t norm, const char * name ="fMtExp");

  // Simple exponential in 1/pt*dNdpt
  TF1 * GetPTExp(Double_t T, Double_t norm, const char * name ="fPtExp");

  // Tsallis (no BW, a la CMS)
  TF1 * GetTsallis(Double_t mass, Double_t T, Double_t q, Double_t norm, const char * name="fTsallis") 
  {return GetLevi (mass,T,1/(q-1),norm,name);}
  
  // Levi function (aka Tsallis)
  TF1 * GetLevi(Double_t mass, Double_t T, Double_t n, Double_t norm, const char * name="fLevi");

  // UA1 function 
  TF1 * GetUA1(Double_t mass, Double_t p0star, Double_t pt0, Double_t n, Double_t T, Double_t norm, const char * name="fUA1");

  // Function derived from a histo
  TF1 * GetHistoFunc(TH1 * h, const char * name = "fHisto");

  // Function derived from a graph
  TF1 * GetGraphFunc(TGraph * h, const char * name = "fHisto");

  // Power law
  TF1 * GetPowerLaw(Double_t pt0, Double_t n, Double_t norm, const char * name="fPowerLaw");
  

  // Bose-Einstein
  TF1 * GetBoseEinstein(Double_t mass, Double_t T, Double_t norm, const char * name="fBoseEinstein");

  // Fermi-Dirac
  TF1 * GetFermiDirac(Double_t mass, Double_t T, Double_t norm, const char * name="fFermiDirac");


  void SetVarType(VarType_t tp) {fVarType=tp;}

protected:

  // dNdpt here means 1/pt dN/dpt
  

  // Boltzmann-Gibbs blast wave
  TF1 * GetBGBWdNdpt(Double_t mass, Double_t beta, Double_t temp,
		     Double_t n, Double_t norm, const char * name= "fBGBW");

  // Tsallis blast wave
  TF1 * GetTsallisBWdNdpt(Double_t mass, Double_t beta, Double_t T, Double_t q,
			  Double_t norm, Double_t ymax = 0.5, const char * name = "fTsallisBW");

  // Simple exponential in 1/mt*MT
  TF1 * GetMTExpdNdpt(Double_t mass, Double_t T, Double_t norm, const char * name ="fExp");

  // Bose-Einstein
  TF1 * GetBoseEinsteindNdpt(Double_t mass, Double_t T, Double_t norm, const char * name="fBoseEinstein");

  // Fermi-Dirac
  TF1 * GetFermiDiracdNdpt(Double_t mass, Double_t T, Double_t norm, const char * name="fFermiDirac");


  // Tsallis (no BW, a la CMS)
  TF1 * GetTsallisdNdpt(Double_t mass, Double_t T, Double_t q, Double_t norm, const char * name="fTsallis");

  // Levi function
  TF1 * GetLevidNdpt(Double_t mass, Double_t T, Double_t n, Double_t norm, const char * name="fLevi");

  // Power Law function
  TF1 * GetPowerLawdNdpt(Double_t pt0, Double_t n, Double_t norm, const char * name="fLevi");

  // UA1 function
  TF1 * GetUA1dNdpt(Double_t pt0, Double_t n, Double_t norm, const char * name="fLevi");

  // TimesPt means dNdpt

  // Boltzmann-Gibbs Blast Wave
  TF1 * GetBGBWdNdptTimesPt(Double_t mass, Double_t beta, Double_t T, Double_t n,
			    Double_t norm, const char * name = "fBGBWTimesPt");

  // Tsallis blast wave
  TF1 * GetTsallisBWdNdptTimesPt(Double_t mass, Double_t beta, Double_t T, Double_t q,
				 Double_t norm, Double_t ymax = 0.5, const char * name = "fTsallisBWTimesPt");
  // Levi function
  TF1 * GetLevidNdptTimesPt(Double_t mass, Double_t T, Double_t n, Double_t norm, const char * name="fLevi");

  // UA1 function
  TF1 * GetUA1dNdptTimesPt(Double_t pt0, Double_t n, Double_t norm, const char * name="fLevi");

  // PowerLaw function
  TF1 * GetPowerLawdNdptTimesPt(Double_t pt0, Double_t n, Double_t norm, const char * name="fLevi");

  // Simple exponential in 1/mt*dNdmT
  TF1 * GetMTExpdNdptTimesPt(Double_t mass, Double_t T, Double_t norm, const char * name ="fMtExpTimesPt");

  // Bose-Einstein
  TF1 * GetBoseEinsteindNdptTimesPt(Double_t mass, Double_t T, Double_t norm, const char * name="fBoseEinstein");

  // Fermi-Dirac
  TF1 * GetFermiDiracdNdptTimesPt(Double_t mass, Double_t T, Double_t norm, const char * name="fFermiDirac");


  // Simple exponential in 1/mp*dNdpT
  TF1 * GetPTExpdNdptTimesPt(Double_t T, Double_t norm, const char * name ="fPtExpTimesPt");

  // Boltzmann (exp in 1/mt*dNdmT times mt)
  TF1 * GetBoltzmanndNdptTimesPt(Double_t mass, Double_t T, Double_t norm, const char * name ="fBoltzmannTimesPt");

  // Tsallis (no BW, a la CMS)
  TF1 * GetTsallisdNdptTimesPt(Double_t mass, Double_t T, Double_t q, Double_t norm, const char * name="fTsallisTimesPt");

  // 1/mt dNdmt
  
  // Levi
  TF1 * GetLevidNdmt(Double_t mass, Double_t T, Double_t n, Double_t norm, const char * name = "fLeviMt", VarType_t var = kOneOverMtdNdmt);

  // gereral setters
  void SetLineWidth(Width_t width) { fLineWidth = width;}

  TF1 * GetTestFunc(Double_t mass, Double_t T, Double_t norm, Double_t ymax, const char * name ="fTest") ;

  // static functions for TF1
  // Boltzmann-Gibbs Blast Wave
  static Double_t StaticBGdNdPt(const double * x, const double* p);
  static Double_t StaticBGdNdPtTimesPt(const double * x, const double* p);
  // Tsallis blast wave
  static Double_t StaticTsallisdNdPt(const double * x, const double* p);
  static Double_t StaticTsallisdNdPtTimesPt(const double * x, const double* p);
  // Helper funcs for numeric integration
  static Double_t IntegrandBG(const double * x, const double* p);
  static Double_t IntegrandTsallis(const double * x, const double* p);

  // Test func
  static Double_t StaticTest(const double * x, const double* p);
  static Double_t IntegrandTest(const double * x, const double* p);

  // histo func
  static Double_t StaticHistoFunc(const double * x, const double* p);

  // UA1 parametrization
  static Double_t StaticUA1Func(const double * x, const double* p);
  

private:


  TF1 * fLastFunc;     // Last function returned
  Width_t fLineWidth;  // Line width
  VarType_t fVarType;  // Variable types (e.g. dNdpt vs pt, 1/mt dNdmt vs mt...) 
  
  AliBWFunc(const AliBWFunc&);            // not implemented
  AliBWFunc& operator=(const AliBWFunc&); // not implemented


  ClassDef(AliBWFunc, 1)


};

#endif
