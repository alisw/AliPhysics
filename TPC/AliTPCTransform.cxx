#include "AliTPCTransform.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "TMath.h"
#include "AliTPCExB.h"

/* To test it:
   cdb=AliCDBManager::Instance()
   cdb->SetDefaultStorage("local:///u/mmager/mycalib1")
   c=AliTPCcalibDB::Instance()
   c->SetRun(0)
   Double_t x[]={1.0,2.0,3.0}
   Int_t i[]={4}
   AliTPCTransform trafo
   trafo.Transform(x,i,0,1)
 */

AliTPCTransform::AliTPCTransform() {
  // Speed it up a bit!
  for (Int_t i=0;i<18;++i) {
    Double_t alpha=TMath::DegToRad()*(10.+20.*(i%18));
    fSins[i]=TMath::Sin(alpha);
    fCoss[i]=TMath::Cos(alpha);
  }
}

AliTPCTransform::~AliTPCTransform() {
}

void AliTPCTransform::Transform(Double_t *x,Int_t *i,UInt_t time,
				Int_t coordinateType) {
  // input: x[0] - pad
  //        x[1] - pad row
  //        x[2] - time in us
  //        i[0] - sector
  // output: x[0] - x (all in the rotated global coordinate frame)
  //         x[1] - y
  //         x[2] - z
  Int_t row=TMath::Nint(x[1]);
  Int_t pad=TMath::Nint(x[0]);
  Int_t sector=i[0];
  AliTPCcalibDB* const calib=AliTPCcalibDB::Instance();
  Double_t xx[3];

  //ugly:  calib->SetRun(time);

  // Time0
  //TODO:  x[2]-=calib->GetPadTime0()->GetCalROC(sector)->GetValue(row,pad);

  // Drift Velocity
  // (cm/us)
  // TODO: use a map or parametrisation!
  x[2]*=2.66;

  Pad2RotatedGlobal(pad,row,x);

  // Alignment
  //TODO:  calib->GetParameters()->GetClusterMatrix(sector)->LocalToMaster(x,xx);

  RotatedGlobal2Global(sector,x);

  // ExB
  calib->GetExB()->Correct(x,xx);

  Global2RotatedGlobal(sector,xx);

  x[0]=xx[0];x[1]=xx[1];x[2]=xx[2];
}

inline void AliTPCTransform::Pad2RotatedGlobal(Int_t pad,Int_t row,Double_t *x)
  const {
  Float_t tmp[3];
  AliTPCROC::Instance()->GetPositionLocal(0,row,pad,tmp);
  x[0]=tmp[0];
  x[1]=tmp[1];
}

//TODO rotation in the right direction?
inline void AliTPCTransform::RotatedGlobal2Global(Int_t sector,Double_t *x)
  const {
  Double_t cos,sin;
  GetCosAndSin(sector,cos,sin);
  Double_t tmp=x[0];
  x[0]= cos*tmp+sin*x[1];
  x[1]=-sin*tmp+cos*x[1];
}

inline void AliTPCTransform::Global2RotatedGlobal(Int_t sector,Double_t *x)
  const {
  Double_t cos,sin;
  GetCosAndSin(sector,cos,sin);
  Double_t tmp=x[0];
  x[0]= cos*tmp-sin*x[1];
  x[1]= sin*tmp+cos*x[1];
}

inline void AliTPCTransform::GetCosAndSin(Int_t sector,Double_t &cos,
					  Double_t &sin) const {
  cos=fCoss[sector%18];
  sin=fSins[sector%18];
}

ClassImp(AliTPCTransform)

