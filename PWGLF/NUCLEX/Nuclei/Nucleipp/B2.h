#ifndef B2_H
#define B2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// some common functions
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TH1D.h>
#include <TF1.h>
#include <cstdlib>
#endif

template <class T>
inline T* FindObj(TFile* f, const TString& name)
{
//
// check if the object exists
//
	T* obj = dynamic_cast<T*>(f->Get(name.Data()));
	if(obj == 0)
	{
		f->Error("Get","%s not found", name.Data());
		exit(1);
	}
	
	return obj;
}

template <class T>
inline T* FindObj(TFile* f, const TString& dir, const TString& name)
{
//
// check if the object exists
//
	if(dir=="") return FindObj<T>(f,name);
	return FindObj<T>(f, dir + "/" + name + ";1");
}

template <class T>
inline T* FindObj(const TList* l, const TString& name)
{
//
// check if the object exists
//
	T* obj = dynamic_cast<T*>(l->FindObject(name.Data()));
	
	if(obj == 0)
	{
		l->Error("FindObject","%s not found", name.Data());
		exit(1);
	}
	
	return obj;
}

inline TH1D* Divide(const TH1* hX, const TH1* hY, const TString& name)
{
//
// clone and divide
//
	TH1D* q = dynamic_cast<TH1D*>(hX->Clone(name.Data()));
	
	if(q == 0)
	{
		hX->Warning("Clone", "could not clone %s", hX->GetName());
		return q;
	}
	
	q->Sumw2();
	q->Divide(hY);
	
	return q;
}

inline Double_t GetMass(const TString& name)
{
//
// return particle mass
//
	TString particle=name;
	particle.ToLower();
	
	if(particle.Contains("electron")) return 0.000511;
	if(particle.Contains("muon")) return 0.10566;
	if(particle.Contains("pion")) return 0.13957;
	if(particle.Contains("kaon")) return 0.49368;
	if(particle.Contains("proton")) return 0.93827;
	if(particle.Contains("deuteron")) return 1.87561;
	if(particle.Contains("triton")) return 2.80925;
	if(particle.Contains("he3")) return 2.80923;
	if(particle.Contains("alpha")) return 3.727417;
	
	//cerr << "Warning unknown particle " << particle << endl;
	
	return 0;
}

inline TF1* Tsallis(Double_t m0, const TString& name, Double_t xmin=0, Double_t xmax=10)
{
//
// Tsallis distribution
// Journal of Statistical Physics, Vol. 52, Nos. 1/2, 1988
// J. Phys. G: Nucl. Part. Phys. 39 (2012)
//
	TF1* fnc = new TF1(name.Data(), Form("[0]*sqrt(x*x+%f)*pow(1+([1]-1)*sqrt(x*x+%f)/[2],[1]/(1-[1]))/pow(2*TMath::Pi(),3)",m0*m0,m0*m0),xmin,xmax);
	
	fnc->SetParNames("gV","q","T");
	fnc->SetParameters(100., 1.1, 0.07);
	
	fnc->SetParLimits(0, 0., 1.e+7);
	fnc->SetParLimits(1, 1.0001, 3.);
	fnc->SetParLimits(2, 0.001, 0.3);
	
	return fnc;
}

inline TF1* TsallisDYield(Double_t m0, const TString& name, Double_t xmin=0, Double_t xmax=10)
{
//
// Tsallis distribution for differential yield
//
	TF1* fnc = new TF1(name.Data(), Form("[0]*x*sqrt(x*x+%f)*pow(1+([1]-1)*sqrt(x*x+%f)/[2],[1]/(1-[1]))/pow(2.*TMath::Pi(),2)",m0*m0,m0*m0),xmin,xmax);
	
	fnc->SetParNames("gV","q","T");
	fnc->SetParameters(100., 1.1, 0.07);
	
	fnc->SetParLimits(0, 0., 1.e+7);
	fnc->SetParLimits(1, 1.0001, 3.);
	fnc->SetParLimits(2, 0.001, 0.3);
	
	return fnc;
}

inline TF1* PtTsallisDYield(Double_t m0, const TString& name, Double_t xmin=0, Double_t xmax=10)
{
//
// Tsallis distribution for mean pt
//
	TF1* fnc = new TF1(name.Data(), Form("[0]*x*x*sqrt(x*x+%f)*pow(1+([1]-1)*sqrt(x*x+%f)/[2],[1]/(1-[1]))/pow(2.*TMath::Pi(),2)",m0*m0,m0*m0),xmin,xmax);
	
	fnc->SetParNames("gV","q","T");
	fnc->SetParameters(100., 1.1, 0.07);
	
	fnc->SetParLimits(0, 0., 1.e+7);
	fnc->SetParLimits(1, 1.0001, 3.);
	fnc->SetParLimits(2, 0.001, 0.3);
	
	return fnc;
}

inline TF1* QTsallis(Double_t m0, const TString& name, Double_t xmin=0, Double_t xmax=10)
{
//
// q-Tsallis distribution
// Phys. Rev. C 83, 064903 (2011)
// Phys. Rev. C 75, 064901 (2007)
//
	TF1* fnc = new TF1(name.Data(), Form("[0]*([1]-1)*([1]-2)*pow(1+(sqrt(x*x+%f)-%f)/([1]*[2]),-[1])/(2*TMath::Pi()*[1]*[2]*([1]*[2]+%f*([1]-2)))", m0*m0, m0, m0), xmin, xmax);
	
	fnc->SetParNames("dN/dy","n","C");
	fnc->SetParameters(0.05, 7, 0.3);
	
	fnc->SetParLimits(0, 0, 10);
	fnc->SetParLimits(1, 4, 50);
	fnc->SetParLimits(2, 0.01, 10);
	
	return fnc;
}

inline TF1* QTsallisDYield(Double_t m0, const TString& name, Double_t xmin=0, Double_t xmax=10)
{
//
// q-Tsallis distribution for differential yield
//
	TF1* fnc = new TF1(name.Data(), Form("x*[0]*([1]-1)*([1]-2)*pow(1+(sqrt(x*x+%f)-%f)/([1]*[2]),-[1])/([1]*[2]*([1]*[2]+%f*([1]-2)))", m0*m0, m0, m0), xmin, xmax);
	
	fnc->SetParNames("dN/dy","n","C");
	fnc->SetParameters(0.05, 7, 0.3);
	
	fnc->SetParLimits(0, 0, 1);
	fnc->SetParLimits(1, 4, 50);
	fnc->SetParLimits(2, 0.01, 10);
	
	return fnc;
}

inline TF1* PtQTsallisDYield(Double_t m0, const TString& name, Double_t xmin=0, Double_t xmax=10)
{
//
// q-Tsallis distribution for mean pt
//
	TF1* fnc = new TF1(name.Data(), Form("x*x*[0]*([1]-1)*([1]-2)*pow(1+(sqrt(x*x+%f)-%f)/([1]*[2]),-[1])/([1]*[2]*([1]*[2]+%f*([1]-2)))", m0*m0, m0, m0), xmin, xmax);
	
	fnc->SetParNames("dN/dy","n","C");
	fnc->SetParameters(0.05, 7, 0.3);
	
	fnc->SetParLimits(0, 0, 1);
	fnc->SetParLimits(1, 4, 50);
	fnc->SetParLimits(2, 0.01, 10);
	
	return fnc;
}

#endif // B2_H
