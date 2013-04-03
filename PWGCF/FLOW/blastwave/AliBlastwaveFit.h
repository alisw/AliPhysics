#ifndef ALIBLASTWAVEFIT_H
#define ALIBLASTWAVEFIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliBlastwaveFit.h 49869 2012-05-17 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//         Blastwave base Fit Class            //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TMinuit.h"

class AliBlastwaveFit : public TNamed
{
public:
    AliBlastwaveFit(const char *name,Double_t mass);
    AliBlastwaveFit();
    ~AliBlastwaveFit();
    
// virtual method to be filled
    virtual void Initialize() = 0; // should be filled with the initialization of the TF1
    virtual Int_t SetParameter(Int_t ipar,Double_t val) = 0;
    virtual Int_t SetNormalization() = 0;
    virtual Int_t GetNpar() = 0;
    virtual const Float_t GetParStart(Int_t ipar) = 0;
    virtual const Float_t GetParStep(Int_t ipar) = 0;
    virtual const Float_t GetParMin(Int_t ipar) = 0;
    virtual const Float_t GetParMax(Int_t ipar) = 0;
    
    virtual void SetMass(Double_t mass) {fMass=mass;};
    virtual const char *GetParName(Int_t i) {if(i>=0) return "parNoName"; else return "";};
    virtual void SwitchOffFlow(TMinuit *) const {};
    

    virtual const Float_t GetMeanBeta() = 0;
    virtual const Float_t GetMeanBeta(Double_t par[]) = 0;

    Double_t GetMass() const {return fMass;};
    
    TF1 *GetSpectraFit() const {return fFunctionYield;};
    TF1 *GetV2Fit() const {return fFunctionV2;};
    
    Double_t EvalYield(Double_t pt){if(fFunctionYield) return fFunctionYield->Eval(pt); else return 0;};
    Double_t EvalV2(Double_t pt){if(fFunctionV2) return fFunctionV2->Eval(pt); else return 0;};
    
    void SetSpectrumObj(TObject *obj){fSpectraObj = obj;if(fSpectraObj && fSpectraObj->InheritsFrom("TH1")) fSpectraObjCopy=(TH1 *)fSpectraObj;else if(fSpectraObj && fSpectraObj->InheritsFrom("TGraphErrors")) fSpectraObjCopy=((TGraphErrors *)fSpectraObj)->GetHistogram();};
    void SetV2Obj(TObject *obj){fV2Obj = obj;};
   
    TObject *GetSpectrumObj(){return fSpectraObj;};
    TObject *GetV2Obj(){return fV2Obj;};
    TH1 *GetSpectrumObjCopy(){return fSpectraObjCopy;};

    void SetMinPt(Float_t pt){fXmin = pt;};
    void SetMaxPt(Float_t pt){fXmax = pt;};

    Float_t GetMinPt() const {return fXmin;};
    Float_t GetMaxPt() const {return fXmax;};

    void Terminate(){if(fFunctionYield) fFunctionYield->SetRange(fXmin,fXmax); if(fFunctionV2) fFunctionV2->SetRange(fXmin,fXmax);};


private:
    AliBlastwaveFit(const AliBlastwaveFit & old);
    AliBlastwaveFit& operator=(const AliBlastwaveFit &/*source*/); // ass. op.

    Double_t fMass; // mass value
    
    TF1 *fFunctionYield;  //! function to describe the yield behaviuor for a fixed fMass
    TF1 *fFunctionV2;     //! function to describe the v2 behaviuor for a fixed fMass
    
    TObject *fSpectraObj; // object to be fitted for yield
    TObject *fV2Obj; // object to be fitted for v2

    TH1 *fSpectraObjCopy; //! object to be fitted for yield as histo

    Float_t fXmin; // mininimum pt for the fit
    Float_t fXmax; // maximum pt for the fit 
    
    ClassDef(AliBlastwaveFit,1)  // blastwave base class
};
#endif


