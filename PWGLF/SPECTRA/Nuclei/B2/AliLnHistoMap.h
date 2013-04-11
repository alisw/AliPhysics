#ifndef ALILNHISTOMAP_H
#define ALILNHISTOMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// class for handling histograms
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TMap.h>

class TString;
class TH1;
class TH1D;
class TH2D;
class TAxis;
class TObjString;

class AliLnHistoMap: public TObject
{
  public:
	AliLnHistoMap();
	
	virtual ~AliLnHistoMap();
	
	//virtual Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
	//virtual Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;
	
	TObject* Get(const TString& keyname) const { return fHistoMap->GetValue(keyname.Data()); }
	TObject* Get(const char* keyname) const { return fHistoMap->GetValue(keyname); }
	TObject* Get(const TObjString* key) const { return fHistoMap->GetValue((TObject*)key); }
	
	TMap* GetMap() const { return fHistoMap; }
	
	TObject* Add(const TString& keyname, TObject* value);
	
	TH1D* Add(const TString& name, Int_t nbins, Double_t xmin, Double_t xmax, const TString& title="", const TString& xlabel="", const TString& ylabel="", Bool_t logx=0);
	
	TH2D* Add(const TString& name, Int_t xbins, Double_t xmin, Double_t xmax, Int_t ybins, Double_t ymin, Double_t ymax, const TString& title="", const TString& xlabel="", const TString& ylabel="", Bool_t logx=0, Bool_t logy=0);
	
	Bool_t SetLogXaxis(TH1* h);
	Bool_t SetLogYaxis(TH1* h);
	
  private:

	AliLnHistoMap(const AliLnHistoMap& other);
	AliLnHistoMap& operator=(const AliLnHistoMap& other);
	
	Bool_t SetLogBins(TAxis* axis);

  private:
 
	TMap* fHistoMap; //-> histogram map
	
	ClassDef(AliLnHistoMap, 1)
};

#endif // ALILNHISTOMAP_H
