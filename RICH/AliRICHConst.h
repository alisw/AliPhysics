#ifndef AliRICHConst_h
#define AliRICHConst_h
#include <TMath.h>

static const Int_t kNpadsX = 144;
static const Int_t kNpadsY = 160;
static const Double_t kD2r=TMath::Pi()/180;
static const Double_t kR2d=57.2957795130823229;
static const Double_t kDeg=TMath::Pi()/180;
static const Double_t kRad=1;
static const Double_t kcm=1;
static const Double_t kmm=0.1;
static const Double_t km=100;
static const Double_t knm=1e-7;
static const Double_t kfm=1e-13;
static const Double_t kfermi=1e-13;
static const Int_t    kBad=-101;        //useful static const to mark initial (uninitalised) values


static const int kNCH=7;                //number of RICH chambers 
static const Float_t kadc_satm  = 4096;  //dynamic range (10 bits)
static const Int_t kMaxNeighbours = 24; //max number of neighbours
static const Int_t kCerenkov=50000050;  //??? go to something more general like TPDGCode
static const Int_t kFeedback=50000051;  //??? go to something more general like TPDGCode
#endif
