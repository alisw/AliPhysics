#ifndef AliRICHConst_h
#define AliRICHConst_h
#include <TMath.h>
static const Double_t d2r=TMath::Pi()/180;
static const Double_t r2d=57.2957795130823229;
static const Double_t deg=TMath::Pi()/180;
static const Double_t rad=1;
static const Double_t mm=0.1;
static const Double_t cm=1;
static const Double_t m=100;

static const Float_t adc_satm  = 4096; // dynamic range (10 bits)
static const Int_t kMaxNeighbours = 24; // max number of neighbours
#endif
