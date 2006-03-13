#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "AliTrackFitterRieman.h"
#include "AliLog.h"

ClassImp(AliTrackFitterRieman)

AliTrackFitterRieman::AliTrackFitterRieman():
  AliTrackFitter()
{
  //
  // default constructor
  //
  fAlpha = 0.;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  fSumYY = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;
  fSumZZ = 0;
  fNUsed = 0;
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
  fSumYY = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;
  fSumZZ = 0;
  fNUsed = 0;
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
  fSumYY = rieman.fSumYY;
  for (Int_t i=0;i<9;i++) fSumXZ[i]  = rieman.fSumXZ[i];
  fSumZZ = rieman.fSumZZ;
  fNUsed = rieman.fNUsed;
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
  fSumYY = rieman.fSumYY;
  for (Int_t i=0;i<9;i++) fSumXZ[i]  = rieman.fSumXZ[i];
  fSumZZ = rieman.fSumZZ;
  fNUsed = rieman.fNUsed;
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
  fSumYY = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;
  fSumZZ = 0;
  fNUsed = 0;
  fConv =kFALSE;
}

Bool_t AliTrackFitterRieman::Fit(const TArrayI *volIds,const TArrayI *volIdsFit,
				 AliAlignObj::ELayerID layerRangeMin,
				 AliAlignObj::ELayerID layerRangeMax)
{
  // Fit the track points. The method takes as an input
  // the set of id's (volids) of the volumes in which
  // one wants to calculate the residuals.
  // The following parameters are used to define the
  // range of volumes to be used in the fitting
  // As a result two AliTrackPointArray's obects are filled.
  // The first one contains the space points with
  // volume id's from volids list. The second array of points represents
  // the track extrapolations corresponding to the space points
  // in the first array. The two arrays can be used to find
  // the residuals in the volids and consequently construct a
  // chi2 function to be minimized during the alignment
  // procedures. For the moment the track extrapolation is taken
  // at the space-point reference plane. The reference plane is
  // found using the covariance matrix of the point
  // (assuming sigma(x)=0 at the reference coordinate system.

  Reset();

  Int_t npoints = fPoints->GetNPoints();
  if (npoints < 3) return kFALSE;

  Bool_t isAlphaCalc = kFALSE;
  AliTrackPoint p,plocal;
//   fPoints->GetPoint(p,0);
//   fAlpha = TMath::ATan2(p.GetY(),p.GetX());

  Int_t npVolId = 0;
  fNUsed = 0;
  Int_t *pindex = new Int_t[npoints];
  fX  = new Float_t[npoints];
  fY  = new Float_t[npoints];
  fZ  = new Float_t[npoints];
  fSy = new Float_t[npoints];
  fSz = new Float_t[npoints];
  for (Int_t ipoint = 0; ipoint < npoints; ipoint++)
    {
      fPoints->GetPoint(p,ipoint);
      UShort_t iVolId = p.GetVolumeID();
      if (FindVolId(volIds,iVolId)) {
	pindex[npVolId] = ipoint;
	npVolId++;
      }
      if (volIdsFit != 0x0) {
	if (!FindVolId(volIdsFit,iVolId)) continue;
      }
      else {
	if (iVolId < AliAlignObj::LayerToVolUID(layerRangeMin,0) ||
	    iVolId > AliAlignObj::LayerToVolUID(layerRangeMax,
						AliAlignObj::LayerSize(layerRangeMax-
								       AliAlignObj::kFirstLayer))) continue;
      }
      if (!isAlphaCalc) {
	fAlpha = p.GetAngle();
	isAlphaCalc = kTRUE;
      }
      plocal = p.Rotate(fAlpha);
      AddPoint(plocal.GetX(),plocal.GetY(),plocal.GetZ(),
	       TMath::Sqrt(plocal.GetCov()[3]),TMath::Sqrt(plocal.GetCov()[5]));
      fNUsed++;
    }

  if (npVolId == 0 || fNUsed < 3) {
    delete [] pindex;
    delete [] fX;
    delete [] fY;
    delete [] fZ;
    delete [] fSy;
    delete [] fSz;
    return kFALSE;
  }

  Update();

  delete [] fX;
  delete [] fY;
  delete [] fZ;
  delete [] fSy;
  delete [] fSz;
 
  if (!fConv) {
    delete [] pindex;
    return kFALSE;
  }

  if ((fParams[0] == 0) ||
      ((-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1) <= 0)) {
    delete [] pindex;
    return kFALSE;
  }


  if (fNUsed < fMinNPoints) {
    delete [] pindex;
    return kFALSE;
  }

  fPVolId = new AliTrackPointArray(npVolId);
  fPTrack = new AliTrackPointArray(npVolId);
  AliTrackPoint p2;
  for (Int_t ipoint = 0; ipoint < npVolId; ipoint++)
    {
      Int_t index = pindex[ipoint];
      fPoints->GetPoint(p,index);
      if (GetPCA(p,p2)) {
	Float_t xyz[3],xyz2[3];
	p.GetXYZ(xyz); p2.GetXYZ(xyz2);
	//	printf("residuals %f %d %d %f %f %f %f %f %f\n",fChi2,fNUsed,fConv,xyz[0],xyz[1],xyz[2],xyz2[0]-xyz[0],xyz2[1]-xyz[1],xyz2[2]-xyz[2]);
	fPVolId->AddPoint(ipoint,&p);
	fPTrack->AddPoint(ipoint,&p2);
      }
    }  

  delete [] pindex;

  // debug info
//   Float_t chi2 = 0, chi22 = 0;
//   for (Int_t ipoint = 0; ipoint < npoints; ipoint++)
//     {
//       fPoints->GetPoint(p,ipoint);
//       UShort_t iVolId = p.GetVolumeID();
//       if (volIdFit != 0) {
// 	if (iVolId != volIdFit) continue;
//       }
//       else {
// 	if (iVolId < AliAlignObj::LayerToVolUID(layerRangeMin,0) ||
// 	    iVolId > AliAlignObj::LayerToVolUID(layerRangeMax,AliAlignObj::LayerSize(layerRangeMax-
// 										     AliAlignObj::kFirstLayer))) continue;
//       }
//       plocal = p.Rotate(fAlpha);
//       Float_t delta = (fParams[0]*(plocal.GetX()*plocal.GetX()+plocal.GetY()*plocal.GetY())+
// 		       2.*plocal.GetX()*fParams[1]+
// 		       fParams[2]-
// 		       2.*plocal.GetY())/
// 	              (2.*TMath::Sqrt(plocal.GetCov()[3]));
// //       Float_t delta2 = (fParams[3]+
// // 		       plocal.GetX()*fParams[4]+
// // 		       plocal.GetX()*plocal.GetX()*fParams[5]-
// // 		       plocal.GetZ())/
// // 	              (TMath::Sqrt(plocal.GetCov()[5]));
//       Double_t r = TMath::Sqrt(plocal.GetX()*plocal.GetX()+plocal.GetY()*plocal.GetY());
//       Double_t Rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
//       Float_t delta2 = (fParams[3]+
// 			r*fParams[4]+r*r*r*fParams[4]*Rm1*Rm1/24-
// 			plocal.GetZ())/
// 	               (TMath::Sqrt(plocal.GetCov()[5]));
//       chi2 += delta*delta;
//       chi22 += delta2*delta2;
//       //      printf("pulls %d %d %f %f\n",ipoint,iVolId,delta,delta2);
      
//     }
//   printf("My chi2 = %f + %f = %f\n",chi2,chi22,chi2+chi22);

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
  //   substitution t = 1/(x^2+y^2),   u = 2*x*t, v = 2*y*t,  D0 = R^2 - x0^2- y0^2
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
  fX[fNUsed] = x; fY[fNUsed]=y; fZ[fNUsed]=z; fSy[fNUsed]=sy; fSz[fNUsed]=sz;
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
  fSumYY += v*v*weight;
  //
  // XZ part
  //
  if (1) {
    weight = 1./(sz*sz);
//     fSumXZ[0] +=weight;
//     fSumXZ[1] +=weight*x;   fSumXZ[2] +=weight*x*x; fSumXZ[3] +=weight*x*x*x; fSumXZ[4] += weight*x*x*x*x;
//     fSumXZ[5] +=weight*z;   fSumXZ[6] +=weight*x*z; fSumXZ[7] +=weight*x*x*z;
    fSumZZ += z*z*weight;
  }
  else {
    weight = 1./(sz*sz);
    fSumXZ[0] +=weight;
    Double_t r = TMath::Sqrt(x*x+y*y);
    fSumXZ[1] +=weight*r;   fSumXZ[2] +=weight*r*r; fSumXZ[3] +=weight*z; fSumXZ[4] += weight*r*z;
    // Now the auxulary sums
    fSumXZ[5] +=weight*r*r*r/24; fSumXZ[6] +=weight*r*r*r*r/12; fSumXZ[7] +=weight*r*r*r*z/24;
    fSumXZ[8] +=weight*r*r*r*r*r*r/(24*24);
    fSumZZ += z*z*weight;
  }
}

void AliTrackFitterRieman::Update(){
  //
  //  Rieman update
  //
  //
  for (Int_t i=0;i<6;i++)fParams[i]=0;
  fChi2 = 0;
  fNdf = 0;
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
    TMatrixD  tmp = res*sums.T();
    fChi2 += fSumYY - tmp(0,0);
    fNdf  += fNUsed - 3;
    conv++;
  }
  //
  // XZ part
  //
  if (1) {
    Double_t x0 = -fParams[1]/fParams[0];
    Double_t Rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 

    for (Int_t i=0;i<fNUsed;i++){
      Double_t phi  = TMath::ASin((fX[i]-x0)*Rm1);
      Double_t phi0 = TMath::ASin((0.-x0)*Rm1);
      Double_t weight = 1/fSz[i];
      weight *=weight;
      Double_t dphi = (phi-phi0)/Rm1;
      fSumXZ[0] +=weight;
      fSumXZ[1] +=weight*dphi;
      fSumXZ[2] +=weight*dphi*dphi;
      fSumXZ[3] +=weight*fZ[i];
      fSumXZ[4] +=weight*dphi*fZ[i];
    }

    TMatrixDSym     smatrixz(2);
    smatrixz(0,0) = fSumXZ[0]; smatrixz(0,1) = fSumXZ[1]; smatrixz(1,1) = fSumXZ[2];
    smatrixz.Invert();
    TMatrixD        sumsxz(1,2);
    if (smatrixz.IsValid()){
      sumsxz(0,0)    = fSumXZ[3];
      sumsxz(0,1)    = fSumXZ[4];
      TMatrixD res = sumsxz*smatrixz;
      fParams[3] = res(0,0);
      fParams[4] = res(0,1);
      fParams[5] = 0;
      for (Int_t i=0;i<2;i++)
	for (Int_t j=0;j<=i;j++){
	  (*fCov)(i+3,j+3)=smatrixz(i,j);
	}
      TMatrixD  tmp = res*sumsxz.T();
      fChi2 += fSumZZ - tmp(0,0);
      fNdf  += fNUsed - 2;
      conv++;
    }
  }
  else {
    Double_t Rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
    fSumXZ[1] += fSumXZ[5]*Rm1*Rm1;
    fSumXZ[2] += fSumXZ[6]*Rm1*Rm1 + fSumXZ[8]*Rm1*Rm1*Rm1*Rm1;
    fSumXZ[4] += fSumXZ[7]*Rm1*Rm1;

    TMatrixDSym     smatrixz(2);
    smatrixz(0,0) = fSumXZ[0]; smatrixz(0,1) = fSumXZ[1]; smatrixz(1,1) = fSumXZ[2];
    smatrixz.Invert();
    TMatrixD        sumsxz(1,2);
    if (smatrixz.IsValid()){
      sumsxz(0,0)    = fSumXZ[3];
      sumsxz(0,1)    = fSumXZ[4];
      TMatrixD res = sumsxz*smatrixz;
      fParams[3] = res(0,0);
      fParams[4] = res(0,1);
      fParams[5] = 0;
      for (Int_t i=0;i<2;i++)
	for (Int_t j=0;j<=i;j++){
	  (*fCov)(i+3,j+3)=smatrixz(i,j);
	}
      TMatrixD  tmp = res*sumsxz.T();
      fChi2 += fSumZZ - tmp(0,0);
      fNdf  += fNUsed - 2;
      conv++;
    }
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

Double_t AliTrackFitterRieman::GetDYat(Double_t x) const {
  if (!fConv) return 0.;
  Double_t x0 = -fParams[1]/fParams[0];
  if (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<0) return 0;
  Double_t Rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  if ( 1./(Rm1*Rm1)-(x-x0)*(x-x0)<=0) return 0;
  Double_t res = (x-x0)/TMath::Sqrt(1./(Rm1*Rm1)-(x-x0)*(x-x0));
  if (fParams[0]<0) res*=-1.;
  return res;
}



Double_t AliTrackFitterRieman::GetZat(Double_t x) const {
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

Double_t AliTrackFitterRieman::GetDZat(Double_t x) const {
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
  // Track parameters in the global coordinate system
  Double_t x0 = -fParams[1]/fParams[0]*cos -         1./fParams[0]*sin;
  Double_t y0 =          1./fParams[0]*cos - fParams[1]/fParams[0]*sin;
  if ((-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1) <= 0) return kFALSE;
  Double_t R  = TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1)/
                fParams[0];

  // Define space-point refence plane
  Double_t alphap = p.GetAngle();
  Double_t sinp = TMath::Sin(alphap);
  Double_t cosp = TMath::Cos(alphap);
  Double_t x  = p.GetX()*cosp + p.GetY()*sinp;
  Double_t y  = p.GetY()*cosp - p.GetX()*sinp;
  Double_t x0p= x0*cosp + y0*sinp;
  Double_t y0p= y0*cosp - x0*sinp;
  if ((R*R - (x-x0p)*(x-x0p))<0) {
    AliWarning(Form("Track extrapolation failed ! (Track radius = %f, track circle x = %f, space-point x = %f, reference plane angle = %f\n",R,x0p,x,alphap));
    return kFALSE;
  }
  Double_t temp = TMath::Sqrt(R*R - (x-x0p)*(x-x0p));
  Double_t y1 = y0p + temp;
  Double_t y2 = y0p - temp;
  Double_t yprime = y1;
  if(TMath::Abs(y2-y) < TMath::Abs(y1-y)) yprime = y2;

  // Back to the global coordinate system
  Double_t xsecond = x*cosp - yprime*sinp;
  Double_t ysecond = yprime*cosp + x*sinp;

  // Now Z coordinate and track angles
  Double_t x2 = xsecond*cos + ysecond*sin;
  Double_t zsecond = GetZat(x2);
  Double_t dydx = GetDYat(x2);
  Double_t dzdx = GetDZat(x2);

  // Fill the cov matrix of the track extrapolation point
  Double_t cov[6] = {0,0,0,0,0,0};
  Double_t sigmax = 100*100.;
  cov[0] = sigmax;           cov[1] = sigmax*dydx;      cov[2] = sigmax*dzdx;
  cov[3] = sigmax*dydx*dydx; cov[4] = sigmax*dydx*dzdx;
  cov[5] = sigmax*dzdx*dzdx;

  Float_t  newcov[6];
  newcov[0] = cov[0]*cos*cos-
            2*cov[1]*sin*cos+
              cov[3]*sin*sin;
  newcov[1] = cov[1]*(cos*cos-sin*sin)-
             (cov[3]-cov[0])*sin*cos;
  newcov[2] = cov[2]*cos-
              cov[4]*sin;
  newcov[3] = cov[0]*sin*sin+
            2*cov[1]*sin*cos+
              cov[3]*cos*cos;
  newcov[4] = cov[4]*cos+
              cov[2]*sin;
  newcov[5] = cov[5];

  p2.SetXYZ(xsecond,ysecond,zsecond,newcov);

  return kTRUE;
}
