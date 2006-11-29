#include "AliHMPIDHelix.h" //class header
#include <TPolyLine3D.h>  //Draw()

ClassImp(AliHMPIDHelix)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void  AliHMPIDHelix::Print(Option_t *opt) const
{
// Debug printout
  Printf("%s helix for Q=%i in B=(0,0,%.2f) tesla",opt,fQ,fBz);
  Printf("Helix parametrised  @ x0=(%6.2f,%6.2f,%6.2f) cm p0=(%6.2f,%6.2f,%6.2f) GeV  P=%.2f GeV Theta=%.2f Phi=%.2f degrees", 
                                fX0.X(),fX0.Y(),fX0.Z(),    fP0.Px(),fP0.Py(),fP0.Pz(), fP0.Mag(),fP0.Theta()*TMath::RadToDeg(),fP0.Phi()*TMath::RadToDeg());
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDHelix::Draw(const Option_t *)
{
// Draw helix by a set of points seperated by 1 cm distance 
  TVector3 pos,mom;
  const Int_t kNpoints=500;
  TPolyLine3D *pHelDraw = new TPolyLine3D(kNpoints); pHelDraw->SetLineColor(kGreen);
  for(Int_t i=0;i<kNpoints;i++) {
    Propagate(i,pos,mom);
    pHelDraw->SetPoint(i,pos.X(),pos.Y(),pos.Z());
  }  
  pHelDraw->Draw();
}
