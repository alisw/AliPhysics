#ifndef ALIHBTFITS_H
#define ALIHBTFITS_H
//_________________________________________________
///////////////////////////////////////////////////////////////////////////////////
//
// class AliHBTFits
//
// Sets of methods for fittig correlation functions
//
//
//
// 
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class TF1;
class TString;
class AliHBTFits: public TObject
{
  public:
    AliHBTFits(){;}
    virtual ~AliHBTFits();
    static void FitQOutCylSurf (const TString& hname, Option_t* fopt = "R", Float_t max = 0.0);
    static void FitQSideCylSurf(const TString& hname, Option_t* fopt = "R", Float_t max = 0.0);

    static void FitQOutQSideQLongCylSurf(const TString& hname,Option_t* fopt = "R",
                 Float_t xmax = 0,Float_t ymax = 0, Float_t zmax = 0);

    static void FitQOutQSideCylSurf(const TString& hname,Option_t* fopt = "R",
                 Float_t xmax = 0,Float_t ymax = 0);
    
    static Double_t QOutCylSurf(Double_t *x, Double_t *par);
    static Double_t QSideCylSurf(Double_t *x, Double_t *par);
    
    static Double_t QOutQSideQLongCylSurf(Double_t *x, Double_t *par);
    static Double_t QOutQSideCylSurf(Double_t *x, Double_t *par);

    
//    static double Qromo(double (*fun)(double), double a, double b,double(*choose)())
//    static double Midexp(double (*fun)(double), double aa, double bb, int n)
    
  protected:
  private:
    static TF1*     fgF1; //functions
    static TF1*     fgF2; //functions
    ClassDef(AliHBTFits,1)
};
#endif

