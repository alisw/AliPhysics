#ifndef _AliGenMUONlib_H
#define _AliGenMUONlib_H
#include <TROOT.h>
#include "GenTypeDefs.h"
class AliGenMUONlib :
public TObject
{
 public:
// pions
    static Double_t PtPion(Double_t *px, Double_t *);
    static Double_t PtScal(Double_t pt, Int_t np);
    static Double_t YPion( Double_t *py, Double_t *);
    static Int_t    IpPion();
// kaons
    static Double_t PtKaon(Double_t *px, Double_t *);
    static Double_t YKaon( Double_t *py, Double_t *);
    static Int_t    IpKaon();
// Phi
    static Double_t PtPhi( Double_t *px, Double_t *);
    static Double_t YPhi( Double_t *px, Double_t *);
    static Int_t    IpPhi();
// J/Psi     
    static Double_t PtJpsi( Double_t *px, Double_t *);
    static Double_t YJpsi(Double_t *py, Double_t *);
    static Int_t    IpJpsi();
// Upsilon    
    static Double_t PtUpsilon( Double_t *px, Double_t * );
    static Double_t YUpsilon(Double_t *py, Double_t *);
    static Int_t    IpUpsilon();
//
// Charm    
    static Double_t PtCharm( Double_t *px, Double_t * );
    static Double_t YCharm(Double_t *py, Double_t *);
    static Int_t    IpCharm();
//
// Beauty
    static Double_t PtBeauty( Double_t *px, Double_t * );
    static Double_t YBeauty(Double_t *py, Double_t *);
    static Int_t    IpBeauty();
//
    typedef Double_t (*GenFunc)  (Double_t *, Double_t *);
    typedef Int_t    (*GenFuncIp)();    
    static GenFunc   GetPt(Param_t param);
    static GenFunc   GetY(Param_t param);
    static GenFuncIp GetIp(Param_t param);    
    ClassDef(AliGenMUONlib,1)
};
#endif







