#ifndef ALIITSRESPONSESSD_H
#define ALIITSRESPONSESSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSresponse.h"


/////////////////////////////////////////////////// 
// Response class for SSD                        //
//                                               //
///////////////////////////////////////////////////

class AliITSresponseSSD : public AliITSresponse {

 public:
    AliITSresponseSSD();
    virtual ~AliITSresponseSSD() {;}
    virtual void    SetParamOptions(const char *opt1, const char *opt2) {
	// parameters: "SetInvalid" to simulate the invalid strips
	fOption1=opt1; fOption2=opt2;
    }
    virtual void    ParamOptions(char *opt1,char *opt2) const {
	// options
	strcpy(opt1,fOption1.Data());  strcpy(opt2,fOption2.Data());
    }
    void SetADCpereV(Double_t a=200./30000.0){fADCpereV = a;}
    Double_t DEvToADC(Double_t eV) const {return eV*fADCpereV;}
    Int_t IEvToADC(Double_t eV) const { // Converts electron-hole pairs to
      return ((Int_t) DEvToADC(eV)); }
      
    Double_t  GetCouplingPR() const {// couplings
      return fCouplingPR;}
    Double_t  GetCouplingPL() const {// couplings
      return fCouplingPL;}
    Double_t  GetCouplingNR() const {// couplings
      return fCouplingNR;}
    Double_t  GetCouplingNL() const {// couplings
      return fCouplingNL;}
    virtual void SetCouplings(Double_t pr, Double_t pl, Double_t nr, Double_t nl) {
      fCouplingPR=pr; fCouplingPL=pl; fCouplingNR=nr; fCouplingNL=nl; }

    Int_t GetZSThreshold() const { // ZS threshold
      return fZSThreshold; }
    virtual void SetZSThreshold(Int_t zsth) { fZSThreshold = zsth; }

protected:
    static const Float_t fgkDiffCoeffDefault; //default for fDiffCoeff
    static const TString fgkOption1Default; // default for fOption1
    static const TString fgkOption2Default; // default for fOption2

    static const Double_t fgkfCouplingPR;
    static const Double_t fgkfCouplingPL;
    static const Double_t fgkfCouplingNR;
    static const Double_t fgkfCouplingNL;

    Double_t fADCpereV;        // Constant to convert eV to ADC.

    Double_t  fCouplingPR;  // couplings
    Double_t  fCouplingPL;  // couplings
    Double_t  fCouplingNR;  // couplings
    Double_t  fCouplingNL;  // couplings   

    static const Int_t fgkZSThreshold; // threshold for the zero suppresion
    Int_t fZSThreshold; 

     TString fOption1;         // Simulate invalid strips option
    TString fOption2;         // Not used for the moment

 private:
    AliITSresponseSSD(const AliITSresponseSSD &source); // copy constructor
    AliITSresponseSSD& operator=(const AliITSresponseSSD &source); // ass. op.

    ClassDef(AliITSresponseSSD,4) //Response class for SSD
};
#endif
