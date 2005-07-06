#ifndef ALIGENGSILIB_H
#define ALIGENGSILIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Implementation of AliGenLib for GSI simulations.                        //
// It is an extension of AliMUONLib providing the option for different     //
// parametrisations of pt, y for every particle type                       //
//                                                                         //
// Responsible: Andres.Sandoval@cern.ch                                    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliGenLib.h"
class TRandom;

class AliGenGSIlib :public AliGenLib {
 public:
    GenFunc   GetPt(Int_t param, const char * tname=0) const;
    GenFunc   GetY(Int_t param, const char * tname=0) const;
    GenFuncIp GetIp(Int_t param, const char * tname=0) const;    

    enum constants{kUpsilon, kJPsi, kCharm, kBeauty, kEta, kEtaprime, kOmega, kRho, kKaon, kPion, kPhi, kLambda, kBaryons};

 private:

    static Double_t PtScal(Double_t pt, Int_t np);

// Upsilon
    static Int_t    IpUpsilon(TRandom *ran);
// Upsilon RITMAN   
    static Double_t PtUpsilonRitman( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonRitman(Double_t *py, Double_t *dummy);
// Upsilon FLAT   
    static Double_t PtUpsilonFlat( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonFlat(Double_t *py, Double_t *dummy);
// Upsilon Karel
    static Double_t PtUpsilonKarel( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonKarel(Double_t *py, Double_t *dummy);
// Upsilon MUONlib
    static Double_t PtUpsilonMUON( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonMUON(Double_t *py, Double_t *dummy);


// JPsi 
    static Int_t    IpJpsi(TRandom *ran);
// JPsi FLAT   
    static Double_t PtJpsiFlat( Double_t *px, Double_t *dummy );
    static Double_t YJpsiFlat(Double_t *py, Double_t *dummy);
// JPsi from MUONlib
    static Double_t PtJpsiMUON( Double_t *px, Double_t *dummy );
    static Double_t YJpsiMUON(Double_t *py, Double_t *dummy);
// JPsi from Ritman
    static Double_t PtJpsiRitman( Double_t *px, Double_t *dummy );

    // JPsi from Sergei
    //    static Double_t PtJpsi( Double_t *px, Double_t *dummy );
    //    static Double_t YJpsi(Double_t *py, Double_t *dummy);
    //    static Int_t    IpJpsi(TRandom *ran);


// Charm 
    static Int_t IpCharm(TRandom *ran);
    static Double_t PtCharmFlat( Double_t *px, Double_t *dummy );
    static Double_t PtCharmMUON( Double_t *px, Double_t *dummy );
    static Double_t PtCharmGSI( Double_t *px, Double_t *dummy );
    static Double_t YCharm(Double_t *py, Double_t *dummy);


// Beauty
    static Int_t IpBeauty(TRandom *ran);
    static Double_t PtBeautyFlat( Double_t *px, Double_t *dummy );
    static Double_t PtBeautyMUON( Double_t *px, Double_t *dummy );
    static Double_t PtBeautyGSI( Double_t *px, Double_t *dummy );
    static Double_t YBeauty(Double_t *py, Double_t *dummy);


// Eta
    static Int_t IpEta(TRandom *ran);
    static Double_t PtEtaPHOS( Double_t *px, Double_t *dummy );
    static Double_t YEtaPHOS(Double_t *py, Double_t *dummy);


// Etaprime
    static Int_t IpEtaprime(TRandom *ran);
    static Double_t PtEtaprimePHOS( Double_t *px, Double_t *dummy );
    static Double_t YEtaprimePHOS(Double_t *py, Double_t *dummy);


// Omega
    static Int_t IpOmega(TRandom *ran);
    static Double_t PtOmega( Double_t *px, Double_t *dummy );
    static Double_t YOmega(Double_t *py, Double_t *dummy);


// Rho
   static Int_t IpRho(TRandom *ran);
   static Double_t PtRho( Double_t *px, Double_t *dummy );
   static Double_t YRho(Double_t *py, Double_t *dummy);



// Kaon
    static Int_t IpKaonPHOS(TRandom *ran);
    static Double_t PtKaonPHOS( Double_t *px, Double_t *dummy );
    static Double_t YKaonPHOS(Double_t *py, Double_t *dummy);


// Pion
    static Int_t IpPionPHOS(TRandom *ran);
    static Double_t PtPion( Double_t *px, Double_t *dummy );
    static Double_t YPion(Double_t *py, Double_t *dummy);


// Phi
    static Int_t IpPhi(TRandom *ran);
    static Double_t PtPhiPHOS( Double_t *px, Double_t *dummy );
    static Double_t YPhiPHOS(Double_t *py, Double_t *dummy);


// Lambda
    //    static Double_t PtLambda( Double_t *px, Double_t *dummy );
    //    static Double_t YLambda(Double_t *py, Double_t *dummy);
    //    static Int_t IpLambda(TRandom *ran);


// Baryons
    static Int_t IpBaryons(TRandom *ran);
    static Double_t PtBaryons( Double_t *px, Double_t *dummy );
    static Double_t YBaryons(Double_t *py, Double_t *dummy);



  ClassDef(AliGenGSIlib,0)
};

#endif







