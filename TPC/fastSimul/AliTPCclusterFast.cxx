/*
  gSystem->Load("libSTAT.so");
  .x ~/NimStyle.C
 
 .L $ALICE_ROOT/TPC/fastSimul/AliTPCclusterFast.cxx+
 //
 AliTPCclusterFast::fPRF = new TF1("fprf","gausn");
 AliTPCclusterFast::fTRF = new TF1("ftrf","gausn");
 AliTPCclusterFast::fPRF->SetParameters(1,0,0.5);
 AliTPCclusterFast::fTRF->SetParameters(1,0,0.5);
 //

 AliTPCclusterFast::Simul("aaa.root",50000); 
 gSystem->Load("libSTAT.so");

 TFile f("aaa.root");
 TTree * tree = (TTree*)f.Get("simul");

*/

#include "TObject.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TTreeStream.h"

class AliTPCclusterFast: public TObject {
public:
  AliTPCclusterFast();
  virtual ~AliTPCclusterFast();
  void SetParam(Float_t mnprim, Float_t diff, Float_t y, Float_t z, Float_t ky, Float_t kz);
  void GenerElectrons();
  void Digitize();
  Double_t GetNsec();
  static void Simul(const char* simul, Int_t npoints);
public:
  Float_t fMNprim;     // mean number of primary electrons
  //                   //electrons part input
  Int_t   fNprim;      // mean number of primary electrons
  Int_t   fNtot;       // total number of primary electrons 
  Float_t fDiff;       // diffusion sigma
  Float_t fY;          // y position 
  Float_t fZ;          // z postion 
  Float_t fAngleY;     // y angle - tan(y)
  Float_t fAngleZ;     // z angle - tan z
  //
  //
  //                   // electron part simul
  TVectorD fSec;       // number of secondary electrons
  TVectorD fPosY;      //! position y for each electron
  TVectorD fPosZ;      //! position z for each electron
  TVectorD fGain;      //! gg for each electron
  //
  TVectorD fStatY;     // stat Y  
  TVectorD fStatZ;     // stat Y
  //
  // digitization part
  //
  TMatrixD fDigits;    // response matrix
  static TF1* fPRF;    // Pad response
  static TF1* fTRF;    // Time response function 
  ClassDef(AliTPCclusterFast,1)  // container for
};
ClassImp(AliTPCclusterFast)

TF1 *AliTPCclusterFast::fPRF=0;
TF1 *AliTPCclusterFast::fTRF=0;


AliTPCclusterFast::AliTPCclusterFast(){
  //
  //
  fDigits.ResizeTo(5,7);
}

AliTPCclusterFast::~AliTPCclusterFast(){
}


void AliTPCclusterFast::SetParam(Float_t mnprim, Float_t diff, Float_t y, Float_t z, Float_t ky, Float_t kz){
  //
  //
  fMNprim = mnprim; fDiff = diff;
  fY=y; fZ=z; 
  fAngleY=ky; fAngleZ=kz;
}
Double_t AliTPCclusterFast::GetNsec(){
  //
  // Generate number of secondary electrons
  // copy of procedure implemented in geant
  //
  const Double_t FPOT=20.77E-9, EEND=10E-6, EEXPO=2.2, EEND1=1E-6;
  const Double_t XEXPO=-EEXPO+1, YEXPO=1/XEXPO;
  const Double_t W=20.77E-9;
  Float_t RAN = gRandom->Rndm();
  return TMath::Nint(TMath::Power((TMath::Power(FPOT,XEXPO)*(1-RAN)+TMath::Power(EEND,XEXPO)*RAN),YEXPO)/W);
}

void AliTPCclusterFast::GenerElectrons(){
  //
  //
  //
  //
  const Int_t knMax=1000;
  if (fPosY.GetNrows()<knMax){
    fPosY.ResizeTo(knMax);
    fPosZ.ResizeTo(knMax);
    fGain.ResizeTo(knMax);
    fSec.ResizeTo(knMax);
    fStatY.ResizeTo(3);
    fStatZ.ResizeTo(3);
  }
  fNprim = gRandom->Poisson(fMNprim);  //number of primary electrons
  fNtot=0;
  //
  Double_t sumQ=0;
  Double_t sumYQ=0;
  Double_t sumZQ=0;
  Double_t sumY2Q=0;
  Double_t sumZ2Q=0;
  for (Int_t i=0;i<knMax;i++){ 
    fSec[i]=0;
  }
  for (Int_t iprim=0; iprim<fNprim;iprim++){
    Float_t dN   =  GetNsec();
    fSec[iprim]=dN;
    Double_t yc = fY+(gRandom->Rndm()-0.5)*fAngleY;
    Double_t zc = fZ+(gRandom->Rndm()-0.5)*fAngleZ;
    for (Int_t isec=0;isec<=dN;isec++){
      //
      //
      Double_t y = gRandom->Gaus(0,fDiff)+yc;
      Double_t z = gRandom->Gaus(0,fDiff)+zc;
      Double_t gg = -TMath::Log(gRandom->Rndm());
      fPosY[fNtot]=y;
      fPosZ[fNtot]=z;
      fGain[fNtot]=gg;
      fNtot++;
      sumQ+=gg;
      sumYQ+=gg*y;
      sumY2Q+=gg*y*y;
      sumZQ+=gg*z;
      sumZ2Q+=gg*z*z;
      if (fNtot>=knMax) break;
    }
    if (fNtot>=knMax) break;
  }
  if (sumQ>0){
    fStatY[0]=sumQ;
    fStatY[1]=sumYQ/sumQ;
    fStatY[2]=sumY2Q/sumQ-fStatY[1]*fStatY[1];
    fStatZ[0]=sumQ;
    fStatZ[1]=sumZQ/sumQ;
    fStatZ[2]=sumZ2Q/sumQ-fStatZ[1]*fStatZ[1];
  }
}

void AliTPCclusterFast::Digitize(){
  //
  //
  //
  // 1. Clear digits
  for (Int_t i=0; i<5;i++)
    for (Int_t j=0; j<7;j++){
      fDigits(i,j)=0;
    }
  //
  // Fill digits
  for (Int_t iel = 0; iel<fNtot; iel++){
    for (Int_t di=-2; di<=2;di++)
      for (Int_t dj=-2; dj<=2;dj++){
	Float_t fac = fPRF->Eval(di-fPosY[iel])*fTRF->Eval(dj-fPosZ[iel]);
	fac*=fGain[iel];
	fDigits(2+di,3+dj)+=fac;
      }
  }
}



void AliTPCclusterFast::Simul(const char* fname, Int_t npoints){
  AliTPCclusterFast fast;
  TTreeSRedirector cstream(fname);
  for (Int_t icl=0; icl<npoints; icl++){
    Float_t nprim=(10+40*gRandom->Rndm());
    Float_t diff =0.01 +0.3*gRandom->Rndm();
    Float_t posY = gRandom->Rndm()-0.5;
    Float_t posZ = gRandom->Rndm()-0.5;
    Float_t ky   = 1.*(gRandom->Rndm()-0.5);
    Float_t kz   = 1.*(gRandom->Rndm()-0.5);
    fast.SetParam(nprim,diff,posY,posZ,ky,kz);
    fast.GenerElectrons();
    fast.Digitize();
    cstream<<"simul"<<
      "s.="<<&fast<<
      "\n";
  }
}


/*
  TH2F *hisL = new TH2F("hisL","hisL",10,10,50,100,0,10)

*/
