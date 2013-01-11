/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// class for handling histograms
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TString.h>
#include <TMap.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObjString.h>
#include "AliLnHistoMap.h"

ClassImp(AliLnHistoMap)

AliLnHistoMap::AliLnHistoMap()
: TObject()
, fHistoMap(0)
{
//
// Default constructor
//
	fHistoMap = new TMap();
	fHistoMap->SetOwner(kTRUE);
}

AliLnHistoMap::~AliLnHistoMap()
{
//
// Default destructor
//
	delete fHistoMap;
	
}

Int_t AliLnHistoMap::Write(const char *name, Int_t option, Int_t bsize) const
{
//
// write only the histograms
//
	Int_t nbytes = 0;
	
	TObjString* key;
	TIter hIter(fHistoMap);
	while( (key = (TObjString*)hIter.Next()) )
	{
		if(fHistoMap->GetValue(key)->InheritsFrom("TH1"))
		{
			nbytes += ((TH1*)fHistoMap->GetValue(key))->Write(name,option,bsize);
		}
	}
	
	return nbytes;
}

Int_t AliLnHistoMap::Write(const char *name, Int_t option, Int_t bsize)
{
// write only the histograms
//
	return ((const AliLnHistoMap*)this)->Write(name,option,bsize);
}

TObject* AliLnHistoMap::Add(const TString& keyname, TObject* value)
{
//
// Add a TObject to the histogram map
//
	if(fHistoMap->Contains(keyname.Data()))
	{
		fHistoMap->Warning("Add", "object %s already exists", keyname.Data());
		return 0;
	}
	
	TObjString* key = new TObjString(keyname.Data());
	
	fHistoMap->Add((TObject*)key, value);
	return value;
}

TH1D* AliLnHistoMap::Add(const TString& name, Int_t nbins, Double_t xmin, Double_t xmax, const TString& title, const TString& xlabel, const TString& ylabel)
{
//
// Add a TH1D histogram to the histogram map
//
	TH1D* value = 0;
	if(fHistoMap->Contains(name.Data()))
	{
		fHistoMap->Warning("Add", "histogram %s already exists", name.Data());
	}
	else
	{
		TObjString* key = new TObjString(name.Data());
		value = new TH1D(name.Data(),title.Data(),nbins,xmin,xmax);
		value->SetXTitle(xlabel.Data());
		value->SetYTitle(ylabel.Data());
		
		fHistoMap->Add((TObject*)key, (TObject*)value);
	}
	
	return value;
}

TH2D* AliLnHistoMap::Add(const TString& name, Int_t xbins, Double_t xmin, Double_t xmax, Int_t ybins, Double_t ymin, Double_t ymax, const TString& title, const TString& xlabel, const TString& ylabel)
{
//
// Add a TH2D histogram to the output map
//
	TH2D* value = 0;
	if(fHistoMap->Contains(name.Data()))
	{
		fHistoMap->Warning("Add", "histogram %s already exists", name.Data());
	}
	else
	{
		TObjString* key = new TObjString(name.Data());
		
		value = new TH2D(name.Data(),title.Data(),xbins,xmin,xmax,ybins,ymin,ymax);
		value->SetXTitle(xlabel.Data());
		value->SetYTitle(ylabel.Data());
		
		fHistoMap->Add((TObject*)key, (TObject*)value);
	}
	
	return value;
}
