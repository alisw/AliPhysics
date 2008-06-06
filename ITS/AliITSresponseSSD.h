#ifndef ALIITSRESPONSESSD_H
#define ALIITSRESPONSESSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSresponse.h"
#include <Riostream.h>

/////////////////////////////////////////////////// 
// Response class for SSD                        //
//                                               //
///////////////////////////////////////////////////

class AliITSresponseSSD : public AliITSresponse {

 public:
    AliITSresponseSSD();
    virtual ~AliITSresponseSSD() {;}

    void SetADCpereV(Double_t a=120./24888.9){fADCpereV = a;}
    Double_t DEvToADC(Double_t eV) const {return eV*fADCpereV;}
    Int_t IEvToADC(Double_t eV) const { // Converts electron-hole pairs to
      return ((Int_t) DEvToADC(eV)); }
      
    void SetKeVperADC(Double_t a=86.4/120.){fKeVperADC = a;}
    Double_t ADCToKeV(Double_t adc) const {return adc*fKeVperADC;}

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

    static const Double_t fgkfCouplingPR;  // default value for couplings
    static const Double_t fgkfCouplingPL;  // default value for couplings
    static const Double_t fgkfCouplingNR;  // default value for couplings
    static const Double_t fgkfCouplingNL;  // default value for couplings

    Double_t fADCpereV;        // Constant to convert eV to ADC.
    Double_t fKeVperADC;       // Constant to convert ADC to keV

    Double_t  fCouplingPR;  // couplings
    Double_t  fCouplingPL;  // couplings
    Double_t  fCouplingNR;  // couplings
    Double_t  fCouplingNL;  // couplings   

    static const Int_t fgkZSThreshold; // threshold for the zero suppresion
    Int_t fZSThreshold; 


 private:
    AliITSresponseSSD(const AliITSresponseSSD &source); // copy constructor
    AliITSresponseSSD& operator=(const AliITSresponseSSD &source); // ass. op.

    ClassDef(AliITSresponseSSD,5) //Response class for SSD
};
#endif
