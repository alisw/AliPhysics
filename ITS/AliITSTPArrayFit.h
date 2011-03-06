/* Copyright(c) 2009-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIITSTPARRAYFIT_H
#define ALIITSTPARRAYFIT_H

///////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                           //
// The line is defined by equations (1)                                                      //
// a0*z+a1*x-a0*a1=0 and                                                                     //
// b0*z+b1*y-b0*b1=0                                                                         //
// where x,y,z are NOT the lab axes but z is the lab axis along which the track              //
// has the largest lever arm and x,y are the remaining 2 axis in                             //
// the order of fgkAxisID[z][0], fgkAxisID[z][1]                                             //
// The parameters are fParams[kA0,kB0,kA1,kB1] and the axis chosen as the independent        //
// var. is fParAxis (i.e. if fParAxis==kZ, then a0=ax,b0=bx, a1=ay,b1=by)                    //
//                                                                                           //
//                                                                                           //
// The helix is defined by the equations (2)                                                 //
// X(t) = (dr+R)*cos(phi0) - (R+sum{dRi})*cos(t+phi0) + sum{dRi*cos(phi0+ti)}                //
// Y(t) = (dr+R)*sin(phi0) - (R+sum{dRi})*sin(t+phi0) + sum{dRi*sin(phi0+ti)}                //
// Z(t) = dz - (R+sum{dRi})*t*tg(dip) + sum{dRi*ti}*tg(dip)                                  //
// where dRi is the change of the radius due to the ELoss at parameter ti                    //
//                                                                                           //
// Author: ruben.shahoyan@cern.ch                                                            //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////


#include <TObject.h>
#include <TMath.h>
#include <AliTrackPointArray.h>
class AliSymMatrix;
class AliLog;
class AliParamSolver;


class AliITSTPArrayFit : public TObject
{
 public:
  enum {kFitDoneBit=BIT(14),kCovInvBit=BIT(15),
	kCosmicsBit=BIT(16),kELossBit=BIT(17),
	kIgnoreCovBit=BIT(18),
	kMask=BIT(24)-1};
  enum {kXX=0,kXY=1,kXZ=2,kYX=kXY,kYY=3,kYZ=4,kZX=kXZ,kZY=kYZ,kZZ=5,kScl=6,kNCov};
  enum {kA0,kB0,kA1,kB1};                // line params
  enum {kD0,kPhi0,kR0,kDZ,kDip};         // helix params
  enum {kX,kY,kZ};
  enum {kMaxParam=6,kMaxParamSq = kMaxParam*(kMaxParam+1)/2};
  enum {kLrBeamPime, kLrSPD1,kLrSPD2, kLrShield1, kLrSDD1,kLrSDD2, kLrShield2, kLrSSD1,kLrSSD2,kMaxLrITS};
  //
 public:
  AliITSTPArrayFit();
  AliITSTPArrayFit(Int_t npoints);
  AliITSTPArrayFit(const AliITSTPArrayFit &fit);
  AliITSTPArrayFit& operator= (const AliITSTPArrayFit& src);
  virtual ~AliITSTPArrayFit();
  //
  void          AttachPoints(const AliTrackPointArray* points, Int_t pfirst=-1,Int_t plast=-1);
  Bool_t        SetFirstLast(Int_t pfirst=-1,Int_t plast=-1);
  AliTrackPointArray* GetPoints()                           const {return (AliTrackPointArray*)fkPoints;}
  //
  void          SetBz(Double_t bz)                                {fBz = bz;}
  Double_t      GetBz()                                     const {return fBz;}
  Bool_t        IsHelix()                                   const {return fParAxis<0;}
  Bool_t        IsFieldON()                                 const {return TMath::Abs(fBz)>1e-5;}
  Bool_t        IsTypeCosmics()                             const {return TestBit(kCosmicsBit);}
  Bool_t        IsTypeCollision()                           const {return !IsTypeCosmics();}
  Int_t         GetCharge()                                 const {return fCharge;}
  Int_t         GetSignQB()                                 const {return fBz<0 ? -fCharge:fCharge;}
  void          GetResiduals(Double_t *res, Int_t ipnt)     const;
  void          GetResiduals(Double_t *resPCA, const Double_t* xyz, const Double_t* covI=0, Double_t sclCovI=-1)  const;
  Double_t      GetPosition( Double_t *xyzPCA, const Double_t* xyz, const Double_t* covI=0, Double_t sclCovI=-1)  const;
  Double_t      GetPosition( Double_t *xyzPCA, const AliTrackPoint *pntCovInv,Bool_t useErr=kFALSE) const;
  void          GetResiduals(Double_t *xyzPCA, const AliTrackPoint *pntCovInv,Bool_t useErr=kFALSE) const;
  void          GetPosition(Double_t *xyz, Double_t t)      const;
  void          GetPosition(Double_t *xyz, Int_t pnt)       const;
  void          GetDirCos(Double_t *dircos, Double_t t)     const;
  Double_t      GetPCA2PlaneInfo(Double_t *xyz, Double_t *dir=0, Int_t axis=kY, Double_t axval=0) const;
  void          GetT0Info(Double_t *xyz, Double_t *dir=0)   const;
  Double_t      CalcChi2NDF()                               const;
  Double_t      GetChi2NDF()                                const {return fChi2NDF;}
  Double_t      GetParPCA(const double *xyz, const double *covI=0, Double_t sclCovI=-1)        const;
  Double_t      CalcParPCA(Int_t ipnt)                      const;
  Bool_t        CalcErrorMatrix();
  //
  void          GetDResDParamsLine (Double_t *dXYZdP, const Double_t *xyz, const Double_t *covI=0/*,Double_t sclCovI=-1*/) const;
  void          GetDResDParamsLine (Double_t *dXYZdP, Int_t ipnt) const;
  void          GetDResDParams(Double_t *dXYZdP, const Double_t *xyz, const Double_t *covI=0, Double_t sclCovI=-1);
  void          GetDResDParams(Double_t *dXYZdP, Int_t ipnt);
  //
  void          GetDResDPosLine (Double_t *dXYZdP,/*const Double_t *xyz,*/ const Double_t *covI=0/*,Double_t sclCovI=-1*/) const;
  void          GetDResDPosLine (Double_t *dXYZdP, Int_t ipnt) const;
  void          GetDResDPos(Double_t *dXYZdP, const Double_t *xyz, const Double_t *covI=0, Double_t sclCovI=-1) const;
  void          GetDResDPos(Double_t *dXYZdP, Int_t ipnt);
  //
  Double_t*     GetPoint(int ip)                            const;
  Bool_t        Converged()                                 const {return fIter<fMaxIter;}
  //
  Double_t      Fit(Int_t extQ=0, Double_t extPT=-1,Double_t extPTerr=0);
  Double_t      FitLine();
  Double_t      FitHelix(Int_t extQ=0, Double_t extPT=-1,Double_t extPTerr=0);
  Bool_t        FitLineCrude();
  Bool_t        FitHelixCrude(Int_t imposedQ=0);
  //
  Int_t         GetParAxis()                                const {return fParAxis;}
  Int_t         GetAxID(Int_t id)                           const {return fkAxID  ? fkAxID[id] : -1;}
  Int_t         GetAxCID(Int_t id)                          const {return fkAxCID ? fkAxCID[id] : -1;}
  Int_t         GetFirst()                                  const {return fPntFirst;}
  Int_t         GetLast()                                   const {return fPntLast;}
  //
  Int_t         GetNParams()                                const {return IsFieldON() ? 5:4;}
  Bool_t        InvertPointsCovMat();
  //
  Int_t*        GetElsId() const {return fElsId;}
  Double_t*     GetElsDR() const {return fElsDR;}
  //
  Double_t      GetCovIScale(Int_t ip)                      const {return ip<fNPBooked ? fCovI[ip*kNCov+kScl] : -1.0;}
  Double_t*     GetCovI(Int_t ip)                           const {return fCovI + kNCov*ip;}
  Double_t*     GetCovI()                                   const {return fCovI;}
  Double_t*     GetParams()                                 const {return (Double_t*)&fParams[0];}
  Double_t      GetParam(Int_t ip)                          const {return fParams[ip];}
  Double_t*     GetTs()                                     const {return (Double_t*)fCurT;}
  Double_t      GetT(Int_t ip)                              const {return fCurT[ip];}
  Double_t      GetLineOffset(Int_t axis)                   const;
  Double_t      GetLineSlope(Int_t axis)                    const;
  //
  Bool_t        IsSwitched2Line()                           const {return fSwitch2Line;}
  Bool_t        IsELossON()                                 const {return TestBit(kELossBit)&&IsFieldON();}
  Bool_t        IsFitDone()                                 const {return TestBit(kFitDoneBit);}
  Bool_t        IsCovInv()                                  const {return TestBit(kCovInvBit);}
  Bool_t        IsCovIgnored()                              const {return TestBit(kIgnoreCovBit);}
  Int_t         GetMaxIterations()                          const {return fMaxIter;}
  Int_t         GetNIterations()                            const {return fIter;}
  Double_t      GetEps()                                    const {return fEps;}
  Double_t      GetMass()                                   const {return fMass;}
  //
  Double_t      GetPt()                                     const;
  Double_t      GetP()                                      const;

  //
  void          Switch2Line(Bool_t v=kTRUE)                       {fSwitch2Line = v;}
  void          SetMaxRforHelix(Double_t r)                       {fMaxRforHelix = r>0 ? r : 1e9;}
  void          SetCharge(Int_t q=1)                              {fCharge = q<0 ? -1:1;}
  void          SetELossON(Bool_t v=kTRUE)                        {SetBit(kELossBit,v);}
  void          SetTypeCosmics(Bool_t v=kTRUE)                    {SetBit(kCosmicsBit,v);}
  void          SetTypeCollision(Bool_t v=kTRUE)                  {SetTypeCosmics(!v);}
  void          SetFitDone(Bool_t v=kTRUE)                        {SetBit(kFitDoneBit,v);}
  void          SetCovInv(Bool_t v=kTRUE)                         {SetBit(kCovInvBit,v);}
  void          SetIgnoreCov(Bool_t v=kTRUE)                      {SetBit(kIgnoreCovBit,v);}
  void          SetParAxis(Int_t ax);
  void          SetMaxIterations(Int_t n=20)                      {fMaxIter = n<2 ? 2:n;}
  void          SetEps(Double_t eps=1e-6)                         {fEps = eps<0 ? GetMachinePrec() : eps;}
  void          SetMass(Double_t m=0.13957)                       {fMass = m<5E-4 ? 5E-4 : m;}
  void          Reset();
  void          BuildMaterialLUT(Int_t ntri=3000);
  void          SetCovIScale(Int_t ip, Double_t scl=-1.0);
  void          ResetCovIScale(Double_t scl=-1.0)                 {for (int i=fNPBooked;i--;) SetCovIScale(i,scl);}
  //
  virtual void  Print(Option_t *opt="")                    const;
  //
  static void   GetNormal(Double_t *norm,const Float_t *covMat);
  //
 protected:
  void          InitAux();
  Int_t         ChoseParAxis()                                              const;
  Double_t      GetParPCALine(const Double_t *xyz, const Double_t *covI=0/*, Double_t sclCovI=-1*/)  const;
  Double_t      GetParPCAHelix(const Double_t *xyz, const Double_t *covI=0, Double_t sclCovI=-1) const;
  Double_t      GetParPCACircle(Double_t x, Double_t y)                     const;
  Double_t      GetHelixParAtR(Double_t r)                                  const;
  //
  void          GetDtDPosLine(Double_t *dtpos,/*const Double_t *xyz,*/  const Double_t *covI=0/*, Double_t sclCovI=-1*/)  const;
  Double_t      GetDtDParamsLine(Double_t *dtparam,const Double_t *xyz, const Double_t *covI=0/*, Double_t sclCovI=-1*/)  const;
  //
  Double_t      GetDRofELoss(Double_t t,Double_t cdip,Double_t rhoL,
			     const Double_t *normS, Double_t &p,Double_t &e) const;
  static Bool_t IsZero(Double_t v,Double_t threshold = 1e-16)     {return TMath::Abs(v)<threshold; }
  static Double_t      GetMachinePrec();
  //
 protected:
  const AliTrackPointArray *fkPoints;               // current points
  AliParamSolver* fParSol;                         // solver for parametric linearized systems
  //
  Double_t  fBz;                                   // magnetic field
  Int_t     fCharge;                               // track charge +1=+, -1=-
  Int_t     fPntFirst;                             // first point to fit
  Int_t     fPntLast;                              // last point to fit
  Int_t     fNPBooked;                             // number of points booked
  Int_t     fParAxis;                              // parameterization axis
  Double_t *fCovI;                                 //! inverted cov.matrix for each point
  Double_t  fParams[kMaxParam];                    // fitted params
  Double_t  fParamsCov[kMaxParamSq];               // fit cov matrix
  Double_t  fChi2NDF;                              // fit chi2/NDF
  Int_t     fMaxIter;                              // max number of iterations
  Int_t     fIter;                                 // real number of iterations
  Double_t  fEps;                                  // precision
  Double_t  fMass;                                 // assumed particle mass for ELoss Calculation
  Bool_t    fSwitch2Line;                          // decided to switch to line
  Double_t  fMaxRforHelix;                         // above this radius use straight line fit
  //
  const Int_t  *fkAxID;                            // axis IDs
  const Int_t  *fkAxCID;                           // axis combinations IDs
  //
  // internal storage
  Double_t *fCurT;                                 // track parameter for each point
  //
  // storage to account e-loss
  Int_t     fFirstPosT;                            // id of the first positive t index in fElsId
  Int_t     fNElsPnt;                              // number of e-loss layers seen by the track 
  Int_t    *fElsId;                                // index of increasing t-ordering in the fCurT
  Double_t *fElsDR;                                // delta_Radius for each e-loss layer
  //
  static       Double_t fgRhoLITS[kMaxLrITS];      // <rho*L> for each material layer
  static const Double_t fgkRLayITS[kMaxLrITS];     // radii of material layers
  static const Double_t fgkZSpanITS[kMaxLrITS];    // half Z span of the material layer
  static const Int_t    fgkPassivLrITS[3];         // list of passive layer enums
  static const Int_t    fgkActiveLrITS[6];         // list of active layer enums
  static const Double_t fgkAlmostZero;             // tiny double
  static const Double_t fgkCQConv;                 // R = PT/Bz/fgkCQConv with GeV,kGauss,cm
  static const Int_t    fgkAxisID[3][3];           // permutations of axis
  static const Int_t    fgkAxisCID[3][6];          // cov matrix elements for axis selection
  
  ClassDef(AliITSTPArrayFit,0);
};

//____________________________________________________
inline void AliITSTPArrayFit::GetPosition(Double_t *xyz, Int_t pnt) const 
{
  // track position at measured point pnt
  GetPosition(xyz,fCurT[pnt]);
}

//____________________________________________________
inline Double_t AliITSTPArrayFit::GetParPCA(const double *xyz, const double *covI, Double_t sclCovI) const
{
  // get parameter for the point with least weighted distance to the point
  if (IsFieldON()) return GetParPCAHelix(xyz,covI,sclCovI);
  else             return GetParPCALine(xyz,covI/*,sclCovI*/);
}

//____________________________________________________
inline Double_t* AliITSTPArrayFit::GetPoint(Int_t ip) const
{
  // get point xyz
  static double xyz[3];
  xyz[kX] = fkPoints->GetX()[ip];
  xyz[kY] = fkPoints->GetY()[ip];
  xyz[kZ] = fkPoints->GetZ()[ip];
  return &xyz[0];
}

//____________________________________________________
inline Double_t AliITSTPArrayFit::Fit(Int_t extQ,Double_t extPT,Double_t extPTerr)
{
  // perform the fit
  if (IsFieldON()) return FitHelix(extQ,extPT,extPTerr);
  else             return FitLine();
}

//____________________________________________________
inline void  AliITSTPArrayFit::SetCovIScale(Int_t ip, Double_t scl) 
{
  // rescale inverted error matrix of specific point
  if (ip>=fNPBooked) return;
  if (TMath::Abs(scl-GetCovIScale(ip))<1e-7) ResetBit(kFitDoneBit); 
  fCovI[ip*kNCov+kScl] = scl;
}

#endif

