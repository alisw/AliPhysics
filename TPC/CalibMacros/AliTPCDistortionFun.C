/// \file AliTPCDistortionFun.C
///
/// Simple compiled macro for declaration of static distortion function
/// on top of the AliTPCDistortion class.
/// Why:
/// 1. Use static function in the fitting procedure
/// 2. Usage in TFormual, TF1, Tf2 ... for visualization.
/// 3. Usage in  tree->Draw() for visualization
/// 4. Simple visualization of fit residuals in multidemension - using tree Draw functionality
///
/// Usage:
/// ~~~
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
/// .L $ALICE_ROOT/TPC/CalibMacros/AliTPCDistortions.cxx+
/// .L $ALICE_ROOT/TPC/CalibMacros/AliTPCDistortionFun.C+
/// ~~~
///
/// Example:
///
/// ~~~{.cpp}
/// // Draw integrated distortion in local x
///
/// TF2 fdistIFCXvZX("fdistIFCXvZX","GetIFCDistortion(y,0,x,0)*sign(x)",-250,250,80,250);
/// fdistIFCXvZX.SetNpx(200);
/// fdistIFCXvZX.SetNpy(200);
/// fdistIFCXvZX.GetXaxis()->SetTitle("Z (cm)");
/// fdistIFCXvZX.GetYaxis()->SetTitle("local X (cm)");
/// fdistIFCXvZX->Draw("colz");
///
/// // Draw local distortion angle dx/dz  in mrad
/// TF2 fangleIFCXvZX("fangleIFCXvZX","1000*(GetIFCDistortion(y,0,x,0)-GetIFCDistortion(y,0,x-1,0))",-250,250,85,245);
/// fangleIFCXvZX.SetNpx(200);
/// fangleIFCXvZX.SetNpy(200);
/// fangleIFCXvZX.GetXaxis()->SetTitle("Z (cm)");
/// fangleIFCXvZX.GetYaxis()->SetTitle("local X (cm)");
/// fangleIFCXvZX->Draw("colz");
///
/// TF2 fangleGGXvZX("fangleGGXvZX","1000*(GetGGDistortion(y,0,x,0,1*400,1*400)-GetGGDistortion(y,0,x-1,0,1*400,1*400))*sign(x)",-250,250,85,245);
/// fangleGGXvZX.SetNpx(200);
/// fangleGGXvZX.SetNpy(200);
/// fangleGGXvZX.GetXaxis()->SetTitle("Z (cm)");
/// fangleGGXvZX.GetYaxis()->SetTitle("local X (cm)");
/// fangleGGXvZX->Draw("colz");
/// ~~~

#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "AliTPCDistortions.h"

TObjArray *arrayPic=new TObjArray;

Double_t GetIFCDistortion(Double_t lx, Double_t ly, Double_t lz, Int_t icoord, Double_t shift=1.){
  ///

  static AliTPCDistortions transform;
  static Bool_t doInit=kTRUE;
  if (doInit){
    transform.SetIFCShift(1);
    transform.InitIFCShiftDistortion();
    doInit=kFALSE;
  }

  Double_t xyzIn[3]={lx, ly, lz};
  Double_t xyzOut[3]={lx, ly, lz};
  Int_t dummyROC=0;
  if (lz<0) { dummyROC=36;}
  transform.UndoIFCShiftDistortion(xyzIn,xyzOut,dummyROC);
  Double_t result=0;
  if (icoord<3) result=xyzOut[icoord]-xyzIn[icoord];
  return result*shift;
}



Double_t GetGGDistortion(Double_t lx, Double_t ly, Double_t lz, Int_t icoord, Double_t deltaVGGA=1., Double_t deltaVGGC=1.){
  /// GG distortion induced distortions

  static AliTPCDistortions transform;
  static Bool_t doInit=kTRUE;
  if (doInit){
    transform.SetDeltaVGGA(1.);
    transform.SetDeltaVGGC(1.);
    transform.InitGGVoltErrorDistortion();
    doInit=kFALSE;
  }
  Double_t xyzIn[3]={lx, ly, lz};
  Double_t xyzOut[3]={lx, ly, lz};
  Int_t dummyROC=0;
  if (lz<0) { dummyROC=36;}
  transform.UndoGGVoltErrorDistortion(xyzIn,xyzOut,dummyROC);
  Double_t result=0;
  if (icoord<3) result=xyzOut[icoord]-xyzIn[icoord];
  if (lz<0) result*=deltaVGGA;
  if (lz>0) result*=deltaVGGC;
  return result;
}


void MakePicIFCDX(){
  ///

  TCanvas *canvasIFC1D= new TCanvas("IFC radial shift", "IFC radial shift");
  TLegend *legend = new TLegend(0.1,0.60,0.5,0.9, "Radial distortion due IFC shift 1 mm ");

  for (Int_t iradius=0; iradius<=10; iradius++){
    Double_t radius=85+iradius*(245-85)/10.;
    TF1 *f1= new TF1(Form("fIFC_ZX%f",radius),Form("0.1*GetIFCDistortion(%f,0,x,0)*sign(x)",radius),-250,250);
    f1->SetMaximum( 0.6);
    f1->SetMinimum(-0.6);
    f1->SetNpx(200);
    f1->SetLineColor(1+((20+iradius)%20));
    f1->SetLineWidth(1);
    if (iradius==0) f1->Draw("p");
    f1->Draw("samep");
    legend->AddEntry(f1,Form("R=%f",radius));    
  }
  legend->Draw();  
}

void MakePicIFCDXangle(){
  ///

  TCanvas *canvasIFC1D= new TCanvas("IFC radial shift - angle", "IFC radial shift angle");
  TLegend *legend = new TLegend(0.1,0.60,0.9,0.9, "Radial distortion due IFC shift 1 mm ");

  for (Int_t iradius=0; iradius<=10; iradius++){
    Double_t radius=85+iradius*(245-85)/10.;
    TF1 *f1= new TF1(Form("fIFC_ZX%f",radius),Form("0.1*1000*(GetIFCDistortion(%f,0,x,0)-GetIFCDistortion(%f,0,x+1,0))",radius,radius),-250,250);
    f1->SetMaximum( 10);
    f1->SetMinimum(-10);
    f1->SetNpx(200);
    f1->SetLineColor(1+(iradius%10));
    f1->SetLineWidth(1);
    //f1->SetLineStyle(2+(iradius%6));
    if (iradius==0) f1->Draw("");
    f1->Draw("same");
    legend->AddEntry(f1,Form("R=%f",radius));    
  }
  legend->Draw();  
}

