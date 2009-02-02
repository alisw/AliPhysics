#include "AliRieman.h"
#include "TTreeStream.h"
#include "TRandom.h"
#include "AliExternalTrackParam.h"

/*
  //
  // Test Progaram for AliRieman class
  //
  // How to use it:
  // 1. Load compiled macros
  // 2. Run function TestRieman()
  // 3. Open file and tree with results
  // 4. Check results - residuals - pulls
  // See example bellow:
  
  .L AliGenInfo.C+
  .L TestRieman.C+
  TestRieman();

  TFile  f("TestRieman.root");
  TTree * tree = (TTree*)f.Get("Test");
  tree->Draw("Rieman.fZ-RiemanRef.fZ:Rieman.fX")

  AliComparisonDraw comp;
  comp->fTree = tree;
  comp->DrawXY("Par0.fP[2]","(Par0.fP[4]-ParR.fP[4])/sqrt(ParR.fC[14])","1","1",5,-0.1,0.1,-5.1,5.1);
  comp->fRes->Draw();

*/



void GetProlongation(Double_t xk, Double_t &x, Double_t *param, Double_t &y, Double_t &z){
  Double_t dx=xk-x;  
  Double_t f1=param[2], f2=f1 + param[4]*dx;
  Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
  y  = param[0] + dx*(f1+f2)/(r1+r2);
  //z  = param[1] + dx*(r2 + f2*(f1+f2)/(r1+r2))*param[3];
  //z    = param[1] + dx*(f1+f2)/(f1*r2 + f2*r1)*param[3];
  Double_t dy   = y-param[0];
  Double_t dfi  = 2*TMath::ASin(0.5*TMath::Sqrt(dx*dx+dy*dy)*TMath::Abs(param[4]));
  Double_t sign = (dx>0) ?  1:-1;
  z  = param[1] + sign*param[3]*dfi/param[4];  
}


void TestRieman(){
  const Double_t kB2C=0.299792458e-3;
  TTreeSRedirector cstream("TestRieman.root");
  AliRieman rieman(1000);
  AliRieman riemanR(1000);
  AliRieman riemanRef(1000);
  Int_t npoints =150;
  Double_t param[5], paramIdeal[5], paramR[5];

  for (Int_t i=0; i<10000; i++){
    //
    // random samples
    //
    Double_t r     = 600*(1+0.5*(gRandom->Rndm()-0.5));  
    Double_t kz    = 2*(gRandom->Rndm()-0.5);
    Double_t snp   = (gRandom->Rndm()-0.5)*0.3;
    Double_t sy    =0.1, sz =0.1;
    Double_t sign  = (gRandom->Rndm()>0.5? -1:1);

    Double_t x0    = 0;
    Double_t y0    = snp*x0+10*(gRandom->Rndm()-0.5);
    Double_t z0    = kz*x0+20*(gRandom->Rndm()-0.5);
    param[0] = y0;
    param[1] = z0;
    param[2] = snp;
    param[3] = kz;
    param[4] = sign/r;
    Double_t covar[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Double_t covarI[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Double_t covarR[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    AliExternalTrackParam trackParam(x0, 0, param,  covar);
    //    
    //
    rieman.Reset();
    riemanR.Reset();
    riemanRef.Reset();
    //    for (Float_t dx =0; dx<90; dx+=5){trackParam.PropagateTo(dx,1./kB2C);}
    for (Int_t ipoint =0; ipoint<npoints; ipoint++){
      Double_t x   =  90 + ipoint;
      Double_t y,z; 
      GetProlongation(x, x0, param,y,z);  
      //       trackParam.PropagateTo(x,1./kB2C);     
      //       y = trackParam.GetY();
      //       z = trackParam.GetZ();
      //      for (Float_t dx =0; dx<1; dx+=1){trackParam.PropagateTo(x+dx,1./kB2C);}
      //      Double_t xyz[3];
      // trackParam.GetXYZAt(x,1./kB2C,xyz);
      //       y = xyz[1];
      //       z = xyz[2];
      //
      rieman.AddPoint(x,y,z,sy,sz);
      riemanR.AddPoint(x, gRandom->Gaus(y,sy), gRandom->Gaus(z,sz), sy,sz);
    }
    rieman.Update();
    riemanR.Update();
    //
    // track extrapolation errrs
    //
    for (Int_t ipoint =0; ipoint<npoints; ipoint++){
      Double_t x =  rieman.GetX()[ipoint];
      Double_t y = riemanR.GetYat(x);
      Double_t z = riemanR.GetZat(x);
      Double_t sty = riemanR.GetErrY(x);
      Double_t stz = riemanR.GetErrZ(x);
      riemanRef.AddPoint(x, y,z,sty,stz);
    }
    riemanRef.Update();
    rieman.GetExternalParameters(x0, paramIdeal,covarI);
    riemanR.GetExternalParameters(x0, paramR,covarR);

    AliExternalTrackParam trackParamI(x0, 0, paramIdeal,  covarI);
    AliExternalTrackParam trackParamR(x0, 0, paramR,  covarR);
    AliRieman * res = rieman.MakeResiduals(); // residuals
    AliRieman * resR = riemanR.MakeResiduals(); // residuals
    cstream<<"Test"<<
      "x0="<<x0<<
      "Par0.="<<&trackParam<<
      "ParI.="<<&trackParamI<<
      "ParR.="<<&trackParamR<<
      "Rieman.="<<&rieman<<  
      "RiemanR.="<<&riemanR<<  
      "RiemanRef.="<<&riemanRef<<  
      "Res.="<<res<<
      "ResR.="<<resR<<
      "\n";
    delete res;
    delete resR;
  }  
}

