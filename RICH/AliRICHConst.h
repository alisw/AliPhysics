#ifndef AliRICHConst_h
#define AliRICHConst_h
#include <TMath.h>
const Double_t d2r=TMath::Pi()/180;
const Double_t r2d=57.2957795130823229;
const Double_t deg=TMath::Pi()/180;
const Double_t rad=1;
const Double_t mm=0.1;
const Double_t cm=1;
const Double_t m=100;

//const Float_t zend = 511.+0.15-2*0.001; // z-out position of first chamber???????
//const Float_t zero_supm = 6.; // zero suppression?????????????
//const Float_t sig_noise = 500.; // electronics noise (no. of electrons)?????????

const Float_t adc_satm  = 4096; // dynamic range (10 bits)
const Int_t kMaxNeighbours = 24; // max number of neighbours
#endif
