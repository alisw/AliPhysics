#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "AliTrackFitterRieman.h"

ClassImp(AliTrackFitterRieman)

AliTrackFitterRieman::AliTrackFitterRieman():
  AliTrackFitter()
{
  //
  // default constructor
  //
  fAlpha = 0.;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;
  fConv = kFALSE;
}


AliTrackFitterRieman::AliTrackFitterRieman(AliTrackPointArray *array, Bool_t owner):
  AliTrackFitter(array,owner)
{
  //
  // Constructor
  //
  fAlpha = 0.;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;
  fConv = kFALSE;
}

AliTrackFitterRieman::AliTrackFitterRieman(const AliTrackFitterRieman &rieman):
  AliTrackFitter(rieman)
{
  //
  // copy constructor
  //
  fAlpha = rieman.fAlpha;
  for (Int_t i=0;i<9;i++) fSumXY[i]  = rieman.fSumXY[i];
  for (Int_t i=0;i<9;i++) fSumXZ[i]  = rieman.fSumXZ[i];
  fConv = rieman.fConv;
}

//_____________________________________________________________________________
AliTrackFitterRieman &AliTrackFitterRieman::operator =(const AliTrackFitterRieman& rieman)
{
  // assignment operator
  //
  if(this==&rieman) return *this;
  ((AliTrackFitter *)this)->operator=(rieman);

  fAlpha = rieman.fAlpha;
  for (Int_t i=0;i<9;i++) fSumXY[i]  = rieman.fSumXY[i];
  for (Int_t i=0;i<9;i++) fSumXZ[i]  = rieman.fSumXZ[i];
  fConv = rieman.fConv;

  return *this;
}

AliTrackFitterRieman::~AliTrackFitterRieman()
{
  // destructor
  //
}

void AliTrackFitterRieman::Reset()
{
  // Reset the track parameters and
  // rieman sums
  AliTrackFitter::Reset();
  fAlpha = 0.;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;
  fConv =kFALSE;
}

Bool_t AliTrackFitterRieman::Fit(UShort_t volId,
				 AliTrackPointArray *pVolId, AliTrackPointArray *pTrack,
				 AliAlignObj::ELayerID layerRangeMin,
				 AliAlignObj::ELayerID layerRangeMax)
{
  // Fit the track points. The method takes as an input
  // the id (volid) of the volume to be skipped from fitting.
  // The following two parameters are used to define the
  // range of volumes to be used in the fitting
  // As a result two AliTrackPointArray's obects are filled.
  // The first one contains the space points with
  // volume id = volid. The second array of points represents
  // the track extrapolation corresponding to the space points
  // in the first array. The two arrays can be used to find
  // the residuals in the volid and consequently construct a
  // chi2 function to be minimized during the alignment
  // procedures. For the moment the track extrapolation is taken
  // as follows: in XY plane - at the CDA between track circle
  // and the space point; in Z - the track extrapolation on the Z
  // plane defined by the space point.

  pVolId = pTrack = 0x0;
  fConv = kFALSE;

  Int_t npoints = fPoints->GetNPoints();
  if (npoints < 3) return kFALSE;

  AliTrackPoint p;
  fPoints->GetPoint(p,0);
  fAlpha = TMath::ATan2(p.GetY(),p.GetX());
  Double_t sin = TMath::Sin(fAlpha);
  Double_t cos = TMath::Cos(fAlpha);

  Int_t npVolId = 0;
  Int_t npused = 0;
  Int_t *pindex = new Int_t[npoints];
  for (Int_t ipoint = 0; ipoint < npoints; ipoint++)
    {
      fPoints->GetPoint(p,ipoint);
      UShort_t iVolId = p.GetVolumeID();
      if (iVolId == volId) {
	pindex[npVolId] = ipoint;
	npVolId++;
      }
      if (iVolId < AliAlignObj::LayerToVolUID(layerRangeMin,0) ||
	  iVolId >= AliAlignObj::LayerToVolUID(layerRangeMax,0)) continue;
      Float_t x = p.GetX()*cos + p.GetY()*sin;
      Float_t y = p.GetY()*cos - p.GetX()*sin;
      AddPoint(x,y,p.GetZ(),1,1);
      npused++;
    }

  if (npused < 3) {
    delete [] pindex;
    return kFALSE;
  }

  Update();

  if (!fConv) {
    delete [] pindex;
    return kFALSE;
  }

  if ((fParams[0] == 0) ||
      ((-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1) <= 0)) {
    delete [] pindex;
    return kFALSE;
  }

  pVolId = new AliTrackPointArray(npVolId);
  pTrack = new AliTrackPointArray(npVolId);
  AliTrackPoint p2;
  for (Int_t ipoint = 0; ipoint < npVolId; ipoint++)
    {
      Int_t index = pindex[ipoint];
      fPoints->GetPoint(p,index);
      if (GetPCA(p,p2)) {
	pVolId->AddPoint(ipoint,&p);
	pTrack->AddPoint(ipoint,&p2);
      }
    }  

  delete [] pindex;

  return kTRUE;
}

void AliTrackFitterRieman::AddPoint(Float_t x, Float_t y, Float_t z, Float_t sy, Float_t sz)
{
  //
  //  Rieman update
  //
  //------------------------------------------------------
  // XY direction
  //
  //  (x-x0)^2+(y-y0)^2-R^2=0 ===>
  //
  //  (x^2+y^2 -2*x*x0 - 2*y*y0+ x0^2 -y0^2 -R^2 =0;  ==>
  //
  //   substitution t = 1/(x^2+y^2),   u = 2*x*t, y = 2*y*t,  D0 = R^2 - x0^2- y0^2
  //
  //  1 - u*x0 - v*y0 - t *D0 =0 ;  - linear equation
  //     
  //  next substition   a = 1/y0    b = -x0/y0   c = -D0/y0
  //
  //  final linear equation :   a + u*b +t*c - v =0;
  //
  // Minimization :
  //
  // sum( (a + ui*b +ti*c - vi)^2)/(sigmai)^2 = min;
  //
  // where sigmai is the error of  maesurement  (a + ui*b +ti*c - vi)
  //
  // neglecting error of xi, and supposing  xi>>yi    sigmai ~ sigmaVi ~ 2*sigmay*t  
  //
  //
  // XY part
  //
  Double_t  t  =  x*x+y*y;
  if (t<2) return;
  t            = 1./t;
  Double_t  u  =  2.*x*t;
  Double_t  v  =  2.*y*t;
  Double_t  error = 2.*sy*t;
  error *=error;
  Double_t weight = 1./error;
  fSumXY[0] +=weight;
  fSumXY[1] +=u*weight;      fSumXY[2]+=v*weight;  fSumXY[3]+=t*weight;
  fSumXY[4] +=u*u*weight;    fSumXY[5]+=t*t*weight;
  fSumXY[6] +=u*v*weight;    fSumXY[7]+=u*t*weight; fSumXY[8]+=v*t*weight;
  //
  // XZ part
  //
  weight = 1./sz;
  fSumXZ[0] +=weight;
  fSumXZ[1] +=weight*x;   fSumXZ[2] +=weight*x*x; fSumXZ[3] +=weight*x*x*x; fSumXZ[4] += weight*x*x*x*x;
  fSumXZ[5] +=weight*z;   fSumXZ[6] +=weight*x*z; fSumXZ[7] +=weight*x*x*z;
}

void AliTrackFitterRieman::Update(){
  //
  //  Rieman update
  //
  //
  for (Int_t i=0;i<6;i++)fParams[i]=0;
  Int_t conv=0;
  //
  // XY part
  //
  TMatrixDSym     smatrix(3);
  TMatrixD        sums(1,3);
  //
  //   smatrix(0,0) = s0; smatrix(1,1)=su2; smatrix(2,2)=st2;
  //   smatrix(0,1) = su; smatrix(0,2)=st; smatrix(1,2)=sut;
  //   sums(0,0)    = sv; sums(0,1)=suv; sums(0,2)=svt;

  smatrix(0,0) = fSumXY[0]; smatrix(1,1)=fSumXY[4]; smatrix(2,2)=fSumXY[5];
  smatrix(0,1) = fSumXY[1]; smatrix(0,2)=fSumXY[3]; smatrix(1,2)=fSumXY[7];
  sums(0,0)    = fSumXY[2]; sums(0,1)   =fSumXY[6]; sums(0,2)   =fSumXY[8];
  smatrix.Invert();
  if (smatrix.IsValid()){
    for (Int_t i=0;i<3;i++)
      for (Int_t j=0;j<=i;j++){
	(*fCov)(i,j)=smatrix(i,j);
      }
    TMatrixD  res = sums*smatrix;
    fParams[0] = res(0,0);
    fParams[1] = res(0,1);
    fParams[2] = res(0,2);
    conv++;
  }
  //
  // XZ part
  //
  TMatrixDSym     smatrixz(3);
  smatrixz(0,0) = fSumXZ[0]; smatrixz(0,1) = fSumXZ[1]; smatrixz(0,2) = fSumXZ[2];
  smatrixz(1,1) = fSumXZ[2]; smatrixz(1,2) = fSumXZ[3];
  smatrixz(2,2) = fSumXZ[4];
  smatrixz.Invert();
  if (smatrixz.IsValid()){
    sums(0,0)    = fSumXZ[5];
    sums(0,1)    = fSumXZ[6];
    sums(0,2)    = fSumXZ[7];
    TMatrixD res = sums*smatrixz;
    fParams[3] = res(0,0);
    fParams[4] = res(0,1);
    fParams[5] = res(0,2);
    for (Int_t i=0;i<3;i++)
      for (Int_t j=0;j<=i;j++){
	(*fCov)(i+2,j+2)=smatrixz(i,j);
    }
    conv++;
  }

  //  (x-x0)^2+(y-y0)^2-R^2=0 ===>
  //
  //  (x^2+y^2 -2*x*x0 - 2*y*y0+ x0^2 -y0^2 -R^2 =0;  ==>
  //   substitution t = 1/(x^2+y^2),   u = 2*x*t, y = 2*y*t,  D0 = R^2 - x0^2- y0^2
  //
  //  1 - u*x0 - v*y0 - t *D0 =0 ;  - linear equation
  //     
  //  next substition   a = 1/y0    b = -x0/y0   c = -D0/y0
  //  final linear equation :   a + u*b +t*c - v =0;
  //
  //  fParam[0]  = 1/y0
  //  fParam[1]  = -x0/y0
  //  fParam[2]  = - (R^2 - x0^2 - y0^2)/y0
  if (conv>1) fConv =kTRUE;
  else
    fConv=kFALSE;
}

Double_t AliTrackFitterRieman::GetYat(Double_t x){
  if (!fConv) return 0.;
  Double_t res2 = (x*fParams[0]+fParams[1]);
  res2*=res2;
  res2 = 1.-fParams[2]*fParams[0]+fParams[1]*fParams[1]-res2;
  if (res2<0) return 0;
  res2 = TMath::Sqrt(res2);
  res2 = (1-res2)/fParams[0];
  return res2;
}

Double_t AliTrackFitterRieman::GetDYat(Double_t x){
  if (!fConv) return 0.;
  Double_t x0 = -fParams[1]/fParams[0];
  if (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<0) return 0;
  Double_t Rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  if ( 1./(Rm1*Rm1)-(x-x0)*(x-x0)<=0) return 0;
  Double_t res = (x-x0)/TMath::Sqrt(1./(Rm1*Rm1)-(x-x0)*(x-x0));
  if (fParams[0]<0) res*=-1.;
  return res;
}



Double_t AliTrackFitterRieman::GetZat(Double_t x){
  if (!fConv) return 0.;
  Double_t x0 = -fParams[1]/fParams[0];
  if (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<=0) return 0;
  Double_t Rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  Double_t phi  = TMath::ASin((x-x0)*Rm1);
  Double_t phi0 = TMath::ASin((0.-x0)*Rm1);
  Double_t dphi = (phi-phi0);
  Double_t res = fParams[3]+fParams[4]*dphi/Rm1;
  return res;
}

Double_t AliTrackFitterRieman::GetDZat(Double_t x){
  if (!fConv) return 0.;
  Double_t x0 = -fParams[1]/fParams[0]; 
  if  (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<=0) return 0;
  Double_t Rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  Double_t res = fParams[4]/TMath::Cos(TMath::ASin((x-x0)*Rm1));
  return res;
}


Double_t AliTrackFitterRieman::GetC(){
  return fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1);
}

Bool_t AliTrackFitterRieman::GetXYZat(Double_t r, Float_t *xyz){
  if (!fConv) return kFALSE;
  Double_t res2 = (r*fParams[0]+fParams[1]);
  res2*=res2;
  res2 = 1.-fParams[2]*fParams[0]+fParams[1]*fParams[1]-res2;
  if (res2<0) return kFALSE;
  res2 = TMath::Sqrt(res2);
  res2 = (1-res2)/fParams[0];

  if (!fConv) return kFALSE;
  Double_t x0 = -fParams[1]/fParams[0];
  if (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<=0) return 0;
  Double_t Rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  Double_t phi  = TMath::ASin((r-x0)*Rm1);
  Double_t phi0 = TMath::ASin((0.-x0)*Rm1);
  Double_t dphi = (phi-phi0);
  Double_t res = fParams[3]+fParams[4]*dphi/Rm1;

  Double_t sin = TMath::Sin(fAlpha);
  Double_t cos = TMath::Cos(fAlpha);
  xyz[0] = r*cos - res2*sin;
  xyz[1] = res2*cos + r*sin;
  xyz[2] = res;

  return kTRUE;
}

Bool_t AliTrackFitterRieman::GetPCA(const AliTrackPoint &p, AliTrackPoint &p2) const
{
  // Get the closest to a given spacepoint track trajectory point
  // Look for details in the description of the Fit() method

  if (!fConv) return kFALSE;

  // First X and Y coordinates
  Double_t sin = TMath::Sin(fAlpha);
  Double_t cos = TMath::Cos(fAlpha);
  //  fParam[0]  = 1/y0
  //  fParam[1]  = -x0/y0
  //  fParam[2]  = - (R^2 - x0^2 - y0^2)/y0
  if (fParams[0] == 0) return kFALSE;
  Double_t x0 = -fParams[1]/fParams[0]*cos -         1./fParams[0]*sin;
  Double_t y0 =          1./fParams[0]*cos - fParams[1]/fParams[0]*sin;
  if ((-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1) <= 0) return kFALSE;
  Double_t R  = TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1)/
                fParams[0];

  Double_t x  = p.GetX();
  Double_t y  = p.GetY();
  Double_t dR = TMath::Sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
  Double_t xprime = TMath::Abs(R)*(x-x0)/dR + x0;
  Double_t yprime = TMath::Abs(R)*(y-y0)/dR + y0;

  // Now Z coordinate
  Double_t phi  = TMath::ASin((x-x0)/R);
  Double_t phi0 = TMath::ASin((0.-x0)/R);
  Double_t dphi = (phi-phi0);
  Double_t zprime = fParams[3]+fParams[4]*dphi*R;

  p2.SetXYZ(xprime,yprime,zprime);

  return kTRUE;
}
