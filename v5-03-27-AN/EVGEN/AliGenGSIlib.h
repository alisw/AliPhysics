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
    static Double_t PtUpsilonRitman( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonRitman(const Double_t *py, const Double_t *dummy);
// Upsilon FLAT   
    static Double_t PtUpsilonFlat( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonFlat(const Double_t *py, const Double_t *dummy);
// Upsilon Karel
    static Double_t PtUpsilonKarel( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonKarel(const Double_t *py, const Double_t *dummy);
// Upsilon MUONlib
    static Double_t PtUpsilonMUON( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonMUON(const Double_t *py, const Double_t *dummy);


// JPsi 
    static Int_t    IpJpsi(TRandom *ran);
// JPsi FLAT   
    static Double_t PtJpsiFlat( const Double_t *px, const Double_t *dummy );
    static Double_t YJpsiFlat(const Double_t *py, const Double_t *dummy);
// JPsi from MUONlib
    static Double_t PtJpsiMUON( const Double_t *px, const Double_t *dummy );
    static Double_t YJpsiMUON(const Double_t *py, const Double_t *dummy);
// JPsi from Ritman
    static Double_t PtJpsiRitman( const Double_t *px, const Double_t *dummy );

    // JPsi from Sergei
    //    static Double_t PtJpsi( Double_t *px, Double_t *dummy );
    //    static Double_t YJpsi(Double_t *py, Double_t *dummy);
    //    static Int_t    IpJpsi(TRandom *ran);


// Charm 
    static Int_t IpCharm(TRandom *ran);
    static Double_t PtCharmFlat( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmMUON( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmGSI( const Double_t *px, const Double_t *dummy );
    static Double_t YCharm(const Double_t *py, const Double_t *dummy);


// Beauty
    static Int_t IpBeauty(TRandom *ran);
    static Double_t PtBeautyFlat( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyMUON( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyGSI( const Double_t *px, const Double_t *dummy );
    static Double_t YBeauty(const Double_t *py, const Double_t *dummy);


// Eta
    static Int_t IpEta(TRandom *ran);
    static Double_t PtEtaPHOS( const Double_t *px, const Double_t *dummy );
    static Double_t YEtaPHOS(const Double_t *py, const Double_t *dummy);


// Etaprime
    static Int_t IpEtaprime(TRandom *ran);
    static Double_t PtEtaprimePHOS( const Double_t *px, const Double_t *dummy );
    static Double_t YEtaprimePHOS(const Double_t *py, const Double_t *dummy);


// Omega
    static Int_t IpOmega(TRandom *ran);
    static Double_t PtOmega( const Double_t *px, const Double_t *dummy );
    static Double_t YOmega(const Double_t *py, const Double_t *dummy);


// Rho
   static Int_t IpRho(TRandom *ran);
   static Double_t PtRho( const Double_t *px, const Double_t *dummy );
   static Double_t YRho(const Double_t *py, const Double_t *dummy);



// Kaon
    static Int_t IpKaonPHOS(TRandom *ran);
    static Double_t PtKaonPHOS( const Double_t *px, const Double_t *dummy );
    static Double_t YKaonPHOS(const Double_t *py, const Double_t *dummy);


// Pion
    static Int_t IpPionPHOS(TRandom *ran);
    static Double_t PtPion( const Double_t *px, const Double_t *dummy );
    static Double_t YPion(const Double_t *py, const Double_t *dummy);


// Phi
    static Int_t IpPhi(TRandom *ran);
    static Double_t PtPhiPHOS( const Double_t *px, const Double_t *dummy );
    static Double_t YPhiPHOS(const Double_t *py, const Double_t *dummy);


// Lambda
    //    static Double_t PtLambda( Double_t *px, Double_t *dummy );
    //    static Double_t YLambda(Double_t *py, Double_t *dummy);
    //    static Int_t IpLambda(TRandom *ran);


// Baryons
    static Int_t IpBaryons(TRandom *ran);
    static Double_t PtBaryons( const Double_t *px, const Double_t *dummy );
    static Double_t YBaryons(const Double_t *py, const Double_t *dummy);



  ClassDef(AliGenGSIlib,0)
};

#endif







