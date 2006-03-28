#include "AliRICHHelix.h" //class header
#include <TPolyLine3D.h>  //Draw()

ClassImp(AliRICHHelix)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliRICHHelix::AliRICHHelix(Double_t p,Double_t theta,Double_t phi,Double_t bz):TObject() 
{
  fX0.SetXYZ(0,0,0);                                                       fX=fX0;
  fP0.SetMagThetaPhi(p,theta*TMath::DegToRad(),phi*TMath::DegToRad());     fP=fP0;
  fLen=0; fQ=1;fBz=bz;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void  AliRICHHelix::Print(Option_t *opt) const
{
// Debug printout
  Printf("%s helix for Q=%i in B=(0,0,%.2f) tesla",opt,fQ,fBz);
  Printf("Helix parametrised  @ x0=(%6.2f,%6.2f,%6.2f) cm p0=(%6.2f,%6.2f,%6.2f) GeV  P=%.2f GeV Theta=%.2f Phi=%.2f degrees", 
                                fX0.X(),fX0.Y(),fX0.Z(),    fP0.Px(),fP0.Py(),fP0.Pz(), fP0.Mag(),fP0.Theta()*TMath::RadToDeg(),fP0.Phi()*TMath::RadToDeg());
  Printf("@ %7.2f cm           x=(%6.2f,%6.2f,%6.2f) cm  p=(%6.2f,%6.2f,%6.2f) GeV  P=%.2f GeV Theta=%.2f Phi=%.2f degrees",
           fLen,              fX.X(), fX.Y(), fX.Z(),       fP.Px(),fP.Py(),fP.Pz(), fP.Mag(), fP.Theta()*TMath::RadToDeg(),fP.Phi()*TMath::RadToDeg());
  Printf("             in LORS   rad=(%5.2f,%5.2f)  pc=(%5.2f,%5.2f)",fPosRad.X(), fPosRad.Y(),fPosPc.X(),fPosPc.Y());
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHHelix::Draw(const Option_t *)
{
// Draw helix by a set of points seperated by 1 cm distance 
  
  const Int_t kNpoints=500;
  TPolyLine3D *pHelDraw = new TPolyLine3D(kNpoints); pHelDraw->SetLineColor(kGreen);
  for(Int_t i=0;i<kNpoints;i++) {
    Propagate(i);
    pHelDraw->SetPoint(i,fX.X(),fX.Y(),fX.Z());
  }  
  pHelDraw->Draw();
}
