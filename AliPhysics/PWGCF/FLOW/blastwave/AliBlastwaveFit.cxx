#include<stdio.h>

#include"AliBlastwaveFit.h"

#include"TMath.h"

ClassImp(AliBlastwaveFit);

AliBlastwaveFit::AliBlastwaveFit(const char *name,Double_t mass) :
  TNamed(name,name),
  fMass(mass),
  fFunctionYield(NULL),
  fFunctionV2(NULL),
  fSpectraObj(NULL),
  fV2Obj(NULL),
  fSpectraObjCopy(NULL),
  fXmin(0.2),
  fXmax(3.0)
{
}
//------------------------------------------------------------------------------
AliBlastwaveFit::AliBlastwaveFit() :
  TNamed("BlastwaveFit","BlastwaveFit"),
  fMass(0.0),
  fFunctionYield(NULL),
  fFunctionV2(NULL),
  fSpectraObj(NULL),
  fV2Obj(NULL),
  fSpectraObjCopy(NULL),
  fXmin(0.2),
  fXmax(3.0)
{  
}
//------------------------------------------------------------------------------
AliBlastwaveFit::~AliBlastwaveFit(){
}
//------------------------------------------------------------------------------
void AliBlastwaveFit::SetSpectrumObj(TObject *obj){
  fSpectraObj = obj;
  if(fSpectraObj && fSpectraObj->InheritsFrom("TH1")) 
    fSpectraObjCopy=(TH1 *)fSpectraObj;
  else if(fSpectraObj && fSpectraObj->InheritsFrom("TGraphErrors")){
    TGraphErrors *g = (TGraphErrors *)fSpectraObj;
    Int_t np = g->GetN();
    Float_t xbin[1000];
    for(Int_t i=0;i<np;i++){
      Float_t x = g->GetX()[i];
      Float_t binwidth = g->GetEX()[i];
      if(i==0 && np > 1){
	Float_t binwidth2 = (g->GetX()[i+1] - x)/2;
	if(binwidth2 > binwidth) binwidth = binwidth2;
      }
      else if(i != 0){
	Float_t binwidth2 = (x - g->GetX()[i-1])/2;
	if(binwidth2 > binwidth) binwidth = binwidth2;
      }
      xbin[i] = x - binwidth;
      xbin[i+1] = x + binwidth;
    }
    fSpectraObjCopy= new TH1D(Form("%s%s",g->GetName(),"Copy"),g->GetTitle(),np,xbin);
    for(Int_t i=0;i<np;i++){
      fSpectraObjCopy->SetBinContent(i+1,g->GetY()[i]);
      fSpectraObjCopy->SetBinError(i+1,g->GetEY()[i]);
    }   
  }
}
