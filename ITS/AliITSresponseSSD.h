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
    void SetADCpereV(Double_t a=50./30000.0){fADCpereV = a;}
    Double_t DEvToADC(Double_t eV) const {return eV*fADCpereV;}
    Int_t IEvToADC(Double_t eV) const { // Converts electron-hole pairs to
      return ((Int_t) DEvToADC(eV)); }
      
 protected:
    static const Float_t fgkDiffCoeffDefault; //default for fDiffCoeff
    static const TString fgkOption1Default; // default for fOption1
    static const TString fgkOption2Default; // default for fOption2

    Double_t fADCpereV;        // Constant to convert eV to ADC.

    TString fOption1;         // Simulate invalid strips option
    TString fOption2;         // Not used for the moment

 private:
    AliITSresponseSSD(const AliITSresponseSSD &source); // copy constructor
    AliITSresponseSSD& operator=(const AliITSresponseSSD &source); // ass. op.

    ClassDef(AliITSresponseSSD,4) //Response class for SSD
};
#endif
