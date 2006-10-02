////////////////////////////////////////////////
//  Reconstructed point class for set:ITS     //
////////////////////////////////////////////////
// This class is used to write the tracks information into a file
// The authors  thank Mariana Bondila to have help them to resolve some problems

#include "AliITSIOTrack.h"

ClassImp(AliITSIOTrack)

AliITSIOTrack::AliITSIOTrack():
fLab(-3),
fTPCLab(-3),
fX(0.),
fY(0.),
fZ(0.),
fPx(0.),
fPy(0.),
fPz(0.),
fStateVPhi(0.),
fStateVZ(0.),
fStateVD(0.),
fStateVTgl(0.),
fStateVC(0.),
fRadius(0.),
fPid(0),
fCharge(0),
fMass(0.),
fDz(0.),
fdEdx(0.),
fC00(0.),
fC10(0.),
fC11(0.),
fC20(0.),
fC21(0.),
fC22(0.),
fC30(0.),
fC31(0.),
fC32(0.),
fC33(0.),
fC40(0.),
fC41(0.),
fC42(0.),
fC43(0.),
fC44(0.) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// default constructor

  for (Int_t i=0;i<6;i++) {fIdModules[i]=fIdPoints[i]=-1; fIdPoints[i]=-1;}
}
   

void AliITSIOTrack::SetCovMatrix(Double_t C00, Double_t C10, Double_t C11, Double_t C20, Double_t C21, 
Double_t C22, Double_t C30, Double_t C31, Double_t C32, Double_t C33, Double_t C40, 
Double_t C41, Double_t C42, Double_t C43, Double_t C44){
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// sets the elements of the covariance matrix

  fC00=C00; fC10=C10; fC11=C11; fC20=C20; fC21=C21; fC22=C22; fC30=C30; fC31=C31;
  fC32=C32; fC33=C33; fC40=C40; fC41=C41; fC42=C42; fC43=C43; fC44=C44; 
}

void AliITSIOTrack::GetCovMatrix(Double_t &C00, Double_t &C10, Double_t &C11,
Double_t &C20, Double_t &C21, Double_t &C22, Double_t &C30, Double_t &C31, 
Double_t &C32, Double_t &C33, Double_t &C40, Double_t &C41, Double_t &C42, 
Double_t &C43, Double_t &C44) const {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// gets the elements of the covariance matrix

  C00=fC00; C10=fC10; C11=fC11; C20=fC20; C21=fC21; 
  C22=fC22; C30=fC30; C31=fC31; C32=fC32; C33=fC33;
  C40=fC40; C41=fC41; C42=fC42; C43=fC43; C44=fC44;  
}

