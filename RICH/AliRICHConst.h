#ifndef AliRICHConst_h
#define AliRICHConst_h
#include <TMath.h>
const Double_t d2r=TMath::Pi()/180;
const Double_t r2d=57.2957795130823229;
const Double_t deg=TMath::Pi()/180;
const Double_t rad=1;
const Double_t cm=1;
const Double_t mm=0.1;
const Double_t m=100;
const Double_t nm=1e-7;
const Double_t fm=1e-13;
const Double_t fermi=1e-13;
const Int_t    kBad=-101;        //useful const to mark initial (uninitalised) values


const int kNCH=7;                //number of RICH chambers 
const Float_t adc_satm  = 4096;  //dynamic range (10 bits)
const Int_t kMaxNeighbours = 24; //max number of neighbours
const Int_t kCerenkov=50000050;  //??? go to something more general like TPDGCode
const Int_t kFeedback=50000051;  //??? go to something more general like TPDGCode
#endif
