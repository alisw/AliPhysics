#ifndef B2_H
#define B2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// some common functions
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TH1D.h>
#include <cstdlib>

TObject* FindObj(TFile* f, const TString& name)
{
//
// check if the object exists
//
	TObject* obj = f->Get(name.Data());
	if(obj == 0)
	{
		f->Error("Get","%s not found",name.Data());
		exit(1);
	}
	
	return obj;
}

TObject* FindObj(TFile* f, const TString& dir, const TString& name)
{
//
// check if the object exists
//
	if(dir=="") return FindObj(f,name);
	
	TObject* obj;
	f->GetObject(Form("%s/%s;1",dir.Data(),name.Data()),obj);
	if(obj == 0)
	{
		f->Error("GetObject","%s/%s not found",dir.Data(),name.Data());
		exit(1);
	}
	
	return obj;
}

TObject* FindObj(const TList* l, const TString& name)
{
//
// check if the object exists
//
	TObject* obj = l->FindObject(name.Data());
	if(obj == 0)
	{
		l->Error("FindObject","%s not found",name.Data());
		exit(1);
	}
	
	return obj;
}

TH1D* Divide(const TH1* hX, const TH1* hY, const TString& name)
{
//
// clone and divide
//
	TH1D* q = (TH1D*)hX->Clone(name.Data());
	q->Sumw2();
	q->Divide(hY);
	return q;
}

Double_t GetMass(const TString& name)
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

#endif // B2_H
