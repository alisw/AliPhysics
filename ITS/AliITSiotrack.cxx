////////////////////////////////////////////////
//  Reconstructed point class for set:ITS     //
////////////////////////////////////////////////

#include <TMatrix.h>

#include "AliITSiotrack.h"

ClassImp(AliITSiotrack)

AliITSiotrack::AliITSiotrack() {

  // default creator
  fLab=-3; 
  fX=fZ=fY=0.; 
  fPx=fPy=fPz=0.;
  for (Int_t i=0;i<6;i++) {fIdModules[i]=fIdPoints[i]=-1; fIdPoints[i]=-1;}
  fStateVPhi=0.; fStateVZ=0.; fStateVD=0.; fStateVTgl=0.; fStateVC=0.;
  fRadius=0.; fCharge=0; 
  fC00=fC10=fC11=fC20=fC21=fC22=fC30=fC31=fC32=fC33=fC40=fC41=fC42=fC43=fC44=0.;
      
}
   

void AliITSiotrack::SetCovMatrix(TMatrix *cov) {

  fC00=(*cov)(0,0);
  fC10=(*cov)(1,0);
  fC11=(*cov)(1,1);
  fC20=(*cov)(2,0);
  fC21=(*cov)(2,1);
  fC22=(*cov)(2,2);
  fC30=(*cov)(3,0);
  fC31=(*cov)(3,1);
  fC32=(*cov)(3,2);
  fC33=(*cov)(3,3);
  fC40=(*cov)(4,0);
  fC41=(*cov)(4,1);
  fC42=(*cov)(4,2);
  fC43=(*cov)(4,3);
  fC44=(*cov)(4,4);

}


Double_t * AliITSiotrack::GetCovMatrix() {
// Return Covariance Matrix in array format.

  fcovar[0]=fC00;
  fcovar[1]=fC10;
  fcovar[2]=fC20;
  fcovar[3]=fC30;
  fcovar[4]=fC40;
  fcovar[5]=fC10;
  fcovar[6]=fC11;
  fcovar[7]=fC21;
  fcovar[8]=fC31;
  fcovar[9]=fC41;
  fcovar[10]=fC20;
  fcovar[11]=fC21;
  fcovar[12]=fC22;
  fcovar[13]=fC32;
  fcovar[14]=fC42;
  fcovar[15]=fC30;
  fcovar[16]=fC31;
  fcovar[17]=fC32;
  fcovar[18]=fC33;
  fcovar[19]=fC43;
  fcovar[20]=fC40;
  fcovar[21]=fC41;
  fcovar[22]=fC42;
  fcovar[23]=fC43;
  fcovar[24]=fC44;

  return fcovar;
  
}
/*
Double_t * AliITSiotrack::GetStateVector() {
  
  Double_t statev[5];
  
  statev[0]=fStateVPhi;
  statev[1]=fStateVZ;
  statev[2]=fStateVD;
  statev[3]=fStateVTgl;
  statev[4]=fStateVC;
  
  return statev;

} */


