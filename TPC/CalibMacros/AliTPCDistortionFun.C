/*
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
  .L $ALICE_ROOT/TPC/AliTPCDistortions.cxx+
  .L $ALICE_ROOT/TPC/CalibMacros/AliTPCDistortionFun.C+

  //
  // integrated distortion
  //
  TF2 fdistIFCXvZX("fdistIFCXvZX","GetIFCDistortion(y,0,x,0)",-250,250,80,250);
  fdistIFCXvZX.GetXaxis()->SetTitle("Z (cm)");
  fdistIFCXvZX.GetYaxis()->SetTitle("local X (cm)");  
  fdistIFCXvZX->Draw("colz");
  //
  // local distortion angle mrad
  //
  TF2 fangleIFCXvZX("fangleIFCXvZX","1000*(GetIFCDistortion(y,0,x,0)-GetIFCDistortion(y,0,x-1,0))*sign(x)",-250,250,85,245);
  fangleIFCXvZX.GetXaxis()->SetTitle("Z (cm)");
  fangleIFCXvZX.GetYaxis()->SetTitle("local X (cm)");  
  fangleIFCXvZX->Draw("colz");

*/


#include "AliTPCDistortions.h"


Double_t GetIFCDistortion(Double_t lx, Double_t ly, Double_t lz, Int_t icoord, Double_t shift=1.){
  //
  static AliTPCDistortions transform;
  Double_t xyzIn[3]={lx, ly, lz};
  Double_t xyzOut[3]={lx, ly, lz};
  transform.SetIFCShift(1);
  Int_t dummyROC=0;
  if (lz<0) { dummyROC=36;}
  transform.UndoIFCShiftDistortion(xyzIn,xyzOut,dummyROC);
  Double_t result=0;
  if (icoord<3) result=xyzOut[icoord]-xyzIn[icoord];
  return result*shift;
}
