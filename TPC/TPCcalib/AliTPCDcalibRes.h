#ifndef ALITPCCALIBRES_H
#define ALITPCCALIBRES_H
#include <TSystem.h>
#include <TNamed.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TVectorF.h>
#include <TVectorD.h>
#include <TString.h>
#include <TMath.h>
#include <TGeoGlobalMagField.h>
#include <TGrid.h>
#include <TNDArray.h>
#include <THn.h>
#include <TH1F.h>
#include <TF1.h>
#include <TEnv.h>
#include <TStopwatch.h>
#include <TGraphErrors.h>
#include <TGeoMatrix.h>
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliTPCcalibAlignInterpolation.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliMagF.h"
#include "AliSysInfo.h"
#include "TStatToolkit.h"
#include "AliSymMatrix.h"
#include "AliTPCChebCorr.h"


#define kMaxResid  10.0f    
#define kMaxTgSlp  2.0f

class AliTPCDcalibRes: public TNamed
{
 public:
  enum {kEpanechnikovKernel, kGaussianKernel};  // defined kernels
  enum {kNSect=18,kNSect2=2*kNSect,kNROC=4*kNSect,kNPadRows=159, kNRowIROC=63, kNRowOROC1=64, kNRowOROC2=32};
  enum {kAlignmentBugFixedBit = AliTPCcalibAlignInterpolation::kAlignmentBugFixedBit};
  enum {kExtractMode, kClosureTestMode};
  enum {kDistDone=BIT(0),kDispDone=BIT(1),kSmoothDone=BIT(2)};
  
  // the voxels are defined in following space
  enum {kVoxZ,   // Z/X sector coordinates
	kVoxF,   // y/x in sector coordinates
	kVoxX,   // sector X coordinate
	kVoxV,   // variable within the voxel (delta, stat, etc): last dimension of all THn histos
	kVoxHDim, kVoxDim=kVoxHDim-1};

  enum {kResX,kResY,kResZ,kResD,kResDim,kResDimG=kResDim-1}; // output dimensions
  //
  struct dts_t {  // struct for basic local residual
    Double32_t dy; //[-kMaxResid,kMaxResid,16] 
    Double32_t dz; //[-kMaxResid,kMaxResid,16]
    Double32_t tgSlp; //[kMaxTgSlp,kMaxTgSlp,16]
    UChar_t bvox[kVoxDim]; // voxel bin info: VoxF,kVoxX,kVoxZ
  };

  // structure for closure test residuals
  struct dtc_t
  {
    Int_t t;        // time stamp
    Double32_t dyR; //[-kMaxResid,kMaxResid,16] 
    Double32_t dzR; //[-kMaxResid,kMaxResid,16]
    Double32_t dyC; //[-kMaxResid,kMaxResid,16] 
    Double32_t dzC; //[-kMaxResid,kMaxResid,16]
    Double32_t q2pt;//[-kMaxQ2Pt,kMaxQ2Pt,16]
    Double32_t tgLam;//[-2.,2.,16]
    Double32_t tgSlp;//kMaxTgSlp,kMaxTgSlp,16]
    Float_t x;
    Float_t y;
    Float_t z;
    UChar_t bvox[kVoxDim]; // voxel bin info: kVoxQ,kVoxF,kVoxX,kVoxZ
    //
    dtc_t() {memset(this,0,sizeof(dtc_t));}
  };

  struct bres_t  {
    Float_t D[kResDim];      // values of extracted distortions
    Float_t E[kResDim];      // their errors
    Float_t DS[kResDim];     // smoothed residual
    Float_t DC[kResDim];     // Cheb parameterized residual
    Float_t stat[kVoxHDim];  // statistics: averages of each voxel dimension + entries
    UChar_t bvox[kVoxDim];   // voxel identifier, here the bvox[0] shows number of Q bins used for Y
    UChar_t bsec;            // sector ID (0-35)
    UChar_t flags;           // status flag
    //
    bres_t() {memset(this,0,sizeof(bres_t));}
  };
 
 public:

 AliTPCDcalibRes(int run=0,Long64_t tmin=0,Long64_t tmax=9999999999,const char* resList=0);
 virtual ~AliTPCDcalibRes();
  
 void ProcessFromDeltaTrees();
 void ProcessFromLocalBinnedTrees();
 void Save(const char* name=0);

 void Init();
 void CollectData(int mode = kExtractMode);
 void FillLocalResidualsTrees();
 void FillCorrectedResiduals();
 void ClosureTest();
 void CreateLocalResidualsTrees(int mode);
 void ProcessResiduals();
 void ProcessDispersions();
 void ProcessSectorResiduals(int is);
 void ProcessSectorDispersions(int is);
 void ProcessVoxelResiduals(int np, float* tg, float *dy, float *dz, bres_t& voxRes);
 void ProcessVoxelDispersions(int np, const float* tg, float *dy, bres_t& voxRes);
 //
 void InitGeom();
 THnF* CreateVoxelStatHisto(int sect);
 void    LoadVDrift();
 Float_t GetDriftCorrection(float z, float x, float phi, int rocID);
 Float_t tgpXY(float x, float y, float q2p, float bz);
  
 void WriteStatHistos();
 void LoadStatHistos();
 void WriteResTree();

 void  FixAlignmentBug(int sect, float q2pt, float bz, float& alp, float& x, float &z, float &deltaY, float &deltaZ);

 Bool_t ValidateTrack();
 Bool_t CompareToHelix(float *resHelixY, float *resHelixZ);

 int    CheckResiduals(Bool_t* kill,float &rmsLongMA);


 //------------------------------------ misc. stat. methods

 static void    FitCircle(int np, const float* x, const float* y, 
			  double &xc, double &yc, double &r, float* dy=0);
 static void    DiffToMA(int np, const float *y, const int winLR, float* diffMA);
 static int     DiffToLocLine(int np, const float* x, const float *y, const int nVoisin, float *diffY);
 static int     DiffToMedLine(int np, const float* x, const float *y, const int nVoisin, float *diffY);
 static float   RoFunc(int np, const float* x, const float* y, float b, float &aa);
 static Float_t SelKthMin(int k, int np, float* arr);
 static void    medFit(int np, const float* x, const float* y, float &a, float &b, float *err=0, float delI=0.f);
 static Int_t*  LTMUnbinnedSig(int np, const float *arr, TVectorF &params , Float_t sigTgt, Float_t minFrac=0.7, Bool_t sorted=kFALSE);
 static Float_t MAD2Sigma(int np, float* y);
 static Bool_t  FitPoly2(const float* x,const float* y, const float* w, int np, float *res, float *err);
 static Bool_t  FitPoly1(const float* x,const float* y, const float* w, int np, float *res, float *err);
 static Bool_t  GetTruncNormMuSig(double a, double b, double &mean, double &sig);
 static void    TruncNormMod(double a, double b, double mu0, double sig0, double &muCf, double &sigCf);
 static Double_t GetLogL(TH1F* histo, int bin0, int bin1, double &mu, double &sig, double &logL0);

 //------------------------------------
  

 Int_t   Smooth0(int isect);
 Bool_t  GetSmoothEstimate(int isect, float x, float p, float z, int which, float *res, float *deriv=0);
 Bool_t  GetSmoothEstimateDim(int isect, float x, float p, float z, int dim, float &res, float *deriv=0);
 void    SetKernelType(int tp=kEpanechnikovKernel, float bwX=2.5, float bwP=2.5, float bwZ=2.1, 
		       float scX=1.f,float scP=1.f,float scZ=1.f);
  
 void    CreateCorrectionObject();
 void    InitBinning();
 Int_t   GetXBin(float x);
 Int_t   GetRowID(float x);
  
  Bool_t  FindVoxelBin(int sectID, float x, float y, float z, UChar_t bin[kVoxHDim],float voxVars[kVoxHDim]);
  
  Int_t   GetXBinExact(float x);
  Float_t GetY2X(int ix, int iy);
  Float_t GetY2XLow(int ix, int iy);
  Float_t GetDY2X(int ix);
  Float_t GetDY2XI(int ix);
  Float_t GetX(int i);
  Float_t GetXLow(int i);
  Float_t GetDX(int i);
  Float_t GetDXI(int i);
  Int_t   GetY2XBinExact(float y2x, int ix);
  Int_t   GetY2XBin(float y2x, int ix);
  Int_t   GetZ2XBinExact(float z2x);
  Int_t   GetZ2XBin(float z2x);
  Float_t GetZ2XLow(int iz);
  Float_t GetZ2X(int iz);
  Float_t GetDZ2X();
  Float_t GetDZ2XI();
  void    FindVoxel(float x, float y2x, float z2x, int &ix,int &ip, int &iz);
  void    FindVoxel(float x, float y2x, float z2x, UChar_t &ix,UChar_t &ip, UChar_t &iz);
  void    GetVoxelCoordinates(int isec, int ix, int ip, int iz,float &x, float &p, float &z);
  Double_t GetKernelWeight(double *u2vec, int np) const;

  Long64_t    GetBin2Fill(const UChar_t binVox[kVoxDim], UShort_t bVal) const;
  UShort_t    GetVoxGBin(const UChar_t bvox[kVoxDim]) const;
  UShort_t    GetVoxGBin(int ix,int ip,int iz) const;
  void        GBin2Vox(UShort_t gbin, UChar_t bvox[kVoxDim]) const;
  //
  void     SetRun(int run)                       {fRun = run;}
  void     SetTMinMax(Long64_t tmin=0, Long64_t tmax=9999999999) {fTMin=tmin; fTMax=tmax;}
  void     SetNXBins(int n=kNPadRows)            {fNXBins = n;}
  void     SetNY2XBins(int n=15)                 {fNY2XBins = n;}
  void     SetNZ2XBins(int n=10)                 {fNZ2XBins = n;}
  void     SetMaxTracks(int n=10000000)          {fMaxTracks = n;}
  void     SetFixAligmentBug(Bool_t v=kTRUE)     {fFixAlignmentBug = v;}
  void     SetCacheLearnSize(int n=1)            {fLearnSize = n;}
  void     SetCacheInput(Int_t v=100)            {fCacheInp = v;}
  void     SetSwitchCache(Bool_t v=kFALSE)       {fSwitchCache = v;}
  void     SetApplyZt2Zc(Bool_t v=kTRUE)         {fApplyZt2Zc = v;}
  void     SetResidualList(const char* l)        {fResidualList = l;}
  void     SetOCDBPath(const char* l)            {fOCDBPath = l;}
  void     SetUseErrorInSmoothing(Bool_t v=kTRUE) {fUseErrInSmoothing = v;}
  void     SetNPrimTrackCuts(int n=600)          {fNPrimTracksCut = n;}
  void     SetMinEntriesVoxel(int n=15)          {fMinEntriesVoxel = n;}
  void     SetMinNClusters(int n=30)             {fMinNCl = n;}
  void     SetNVoisinMA(int n=3)                 {fNVoisinMA = n;}
  void     SetNVoisinMALong(int n=15)            {fNVoisinMALong = n;}
  void     SetMaxDevYHelix(float d=0.3)          {fMaxDevYHelix = d;}
  void     SetMaxDevZHelix(float d=0.3)          {fMaxDevZHelix = d;}
  void     SetMaxStdDevMA(float v=25.0)          {fMaxStdDevMA = v;}
  void     SetMaxRMSLong(float v=0.8)            {fMaxRMSLong = v;}
  void     SetMaxRejFrac(float v=0.15)           {fMaxRejFrac = v;}
  void     SetFilterOutliers(Bool_t v=kTRUE)     {fFilterOutliers = v;}

  void     SetMaxFitYErr2(float v=1.0)           {fMaxFitYErr2 = v;}
  void     SetMaxFitXErr2(float v=1.5)           {fMaxFitXErr2 = v;}

  Float_t  GetMaxFitYErr2()                 const {return fMaxFitYErr2;}
  Float_t  GetMaxFitXErr2()                 const {return fMaxFitXErr2;}

  Int_t    GetRun()                         const {return fRun;}
  Long64_t GetTMin()                        const {return fTMin;}
  Long64_t GetTMax()                        const {return fTMax;}  
  Int_t    GetNXBins()                      const {return fNXBins;}
  Int_t    GetNY2XBins()                    const {return fNY2XBins;}
  Int_t    GetNZ2XBins()                    const {return fNZ2XBins;}
  Int_t    GetMaxTracks()                   const {return fMaxTracks;}
  Int_t    GetCacheInput()                  const {return fCacheInp;}
  Int_t    GetCacheLearnSize()              const {return fLearnSize;}
  Int_t    GetNPrimTrackCuts()              const {return fNPrimTracksCut;}
  Int_t    GetMinEntriesVoxel()             const {return fMinEntriesVoxel;}
  Int_t    GetMinNClusters()                const {return fMinNCl;}
  Int_t    GetNVoisinMA()                   const {return fNVoisinMA;}
  Int_t    GetNVoisinMALong()               const {return fNVoisinMALong;}
  Float_t  GetMaxDevYHelix()                const {return fMaxDevYHelix;}
  Float_t  GetMaxDevZHelix()                const {return fMaxDevZHelix;}
  Float_t  GetMaxStdDevMA()                 const {return fMaxStdDevMA;}
  Float_t  GetMaxRMSLong()                  const {return fMaxRMSLong;}
  Float_t  GetMaxRejFrac()                  const {return fMaxRejFrac;}
  Bool_t   GetFilterOutliers()              const {return fFilterOutliers;}

  Bool_t   GetFixAlignmentBug()             const {return fFixAlignmentBug;}
  Bool_t   GetSwitchCache()                 const {return fSwitchCache;}
  Bool_t   GetApplyZt2Zc()                  const {return fApplyZt2Zc;}
  Bool_t   GetUseErrorInSmoothing()         const {return fUseErrInSmoothing;}
  
  const TString& GetOCDBPath()              const {return fOCDBPath;}
  const TString& GetReisdualList()          const {return fResidualList;}

  const AliTPCChebCorr* GetChebCorrObject() const {return fChebCorr;}


  static void SetUsedInstance(AliTPCDcalibRes* inst) {fgUsedInstance = inst;}
  static AliTPCDcalibRes* GetUsedInstance()          {return fgUsedInstance;}
  static float GetTPCRowX(int r)                     {return kTPCRowX[r];}
protected:
  //
  Bool_t   fInitDone;                               // init flag
  Bool_t   fUseErrInSmoothing;                      // weight kernel by point error
  Bool_t   fSwitchCache;                            // reset the cache when the reading mode is changing
  Bool_t   fFixAlignmentBug;                        // flag to apply the fix
  Bool_t   fApplyZt2Zc;                             // Apply fix for using Z_track instead of Z_cluster in the data

  // --------------------------------Chebyshev object creation 
  Int_t    fChebZSlicePerSide;                      // z partitions per side
  Int_t    fChebPhiSlicePerSector;                  // azimuthal partitions per sector
  Int_t    fNPCheb[kResDim][2];                     // cheb. nodes per slice

  Float_t  fChebPrecD[kResDim];                     // nominal precision per output dimension
  AliTPCChebCorr* fChebCorr;                        // final Chebyshev object

  // -------------------------------Task defintion
  Int_t    fRun;     // run numbet 
  Long64_t fTMin;    // time start
  Long64_t fTMax;    // time stop
  Int_t    fMaxTracks;  // max tracks to accept
  Int_t    fCacheInp;      // input trees cache in MB
  Int_t    fLearnSize;     // event to learn for the cache
  Float_t  fBz;            // B field
  Bool_t   fDeleteSectorTrees; // delete residuals trees once statistics tree is done
  TString  fResidualList;  // list of residuals tree
  TString  fOCDBPath;      // ocdb path

  // ------------------------------Selection/filtering cuts
  Int_t    fMinEntriesVoxel;         // min number of entries per voxel to consider
  Int_t    fNPrimTracksCut;          // of >0, cut on event multiplicity
  Float_t  fMinNCl;                  // min number of TPC clusters to consider
  Float_t  fMaxDevYHelix;            // max-min Y deviation of interpolating track from helix
  Float_t  fMaxDevZHelix;            // max-min Z deviation of interpolating track from helix
  Float_t  fNVoisinMA;               // N neighbours for moving average
  Float_t  fNVoisinMALong;           // max RMS of cleaned residuals wrt its fNVoisinMALong moving average
  Float_t  fMaxStdDevMA;             // max cluster N std.dev (Y^2+Z^2) wrt moving av. to accept
  Float_t  fMaxRMSLong;              // max RMS of cleaned residuals wrt its fNVoisinMALong moving average
  Float_t  fMaxRejFrac;              // max outlier clusters tagged to accept the track
  Bool_t   fFilterOutliers;          // reject outliers

  Float_t  fMaxFitYErr2;             // cut on median fit Y err^2
  Float_t  fMaxFitXErr2;             // cut on median fit X err^2

  // -------------------------------Binning
  Int_t    fNY2XBins;    // y/x bins per sector
  Int_t    fNZ2XBins;    // z/x bins per sector
  Int_t    fNXBins;      // n bins in radial dim.
  Int_t    fNXYBinsProd; // nx*ny bins
  Int_t    fNBins[kVoxDim]; // bins
  Bool_t   fUniformBins[kVoxDim]; // uniform binning? Currently only X may be non-uniform (per pad-row)
  

  Float_t  fDZ2X;            // Z2X bin size
  Float_t  fDX;            // X bin size
  Float_t  fDZ2XI;           // inverse Z2X bin size 
  Float_t  fDXI;           // inverse X bin size 

  Int_t    fNGVoxPerSector; // total number of geometrical voxels per sector

  Float_t  *fMaxY2X;        //[fNXBins] max Y/X at each X bin, account for dead zones
  Float_t  *fDY2X;          //[fNXBins] Y/X bin size at given X bin
  Float_t  *fDY2XI;         //[fNXBins] inverse of Y/X bin size at given X bin

  Long64_t fNBProdSt[kVoxHDim];   // aux arrays for fast bin calculation
  Int_t    fNBProdSectG[kVoxDim]; // aux info for fast bin index calculation in geom voxel space
  Int_t    fBins[kVoxDim];        // binning in voxel variables

  // ------------------------------Smoothing
  Int_t    fNMaxNeighb;        // max neighbours to loop for smoothing
  Int_t    fKernelType;        // kernel type
  Int_t    fStepKern[kVoxDim]; // N bins to consider with given kernel settings
  Float_t  fKernelWInv[kVoxDim];      // inverse kernel width in bins
  Float_t  fKernelScaleEdge[kVoxDim]; // optional scaling factors for kernel width on the edge
  // result of last kernel minimization: value and dV/dX,dV/dY,dV/dZ for each dim
  Double_t fLastSmoothingRes[kResDim*4];  

  // ------------------------------Selection Stats
  Int_t    fNTrSelTot;      // selected tracks
  Int_t    fNTrSelTotWO;    // would be selected w/o outliers rejection
  Int_t    fNReadCallTot;   // read calls from input trees
  Long64_t fNBytesReadTot;  // total bytes read


  // ------------------------------VDrift correction
  TVectorD     *fVDriftParam;
  TGraphErrors *fVDriftGraph;  
  Float_t      fCorrTime;   //! 

  // -----------------------------Results of processing
  bres_t *fSectGVoxRes[kNSect2];         //! [fNGVoxPerSector] sectors results for geometric voxel
  TTree* fStatTree;                      //! tree with voxels statistics
  TTree* fTmpTree[kNSect2];              //! IO tree per sector
  TFile* fTmpFile[kNSect2];              //! file for fTmpTree
  THnF*  fStatHist[kNSect2];             //! histos for statistics bins
  TNDArrayT<float> *fArrNDStat[kNSect2]; //! alias arrays for fast access to fStatHist

  TH1F* fHDelY;                          //! work histo for delta Y fits
  TH1F* fHDelZ;                          //! work histo for delta Z fits
  //
  // ----------------------------data exchange structures for trees and between routines
  dts_t fDTS;                            //! binned residuals
  dtc_t fDTC;                            //! corrected residuals for closure test
  //
  // ---------------------------track data-----------------------------------
  int   fTimeStamp;                       //! time stamp
  int   fNCl;                             //! number of clusters
  float fQ2Pt;                            //! fitted q2pt
  float fTgLam;                           //! fitted tgLambda
  float fArrPhi[kNPadRows];               //! cluster phi
  float fArrDY[kNPadRows];                //! cluster residual Y
  float fArrDZ[kNPadRows];                //! cluster residual Z
  float fArrX[kNPadRows];                 //! cluster X (row)
  float fArrYCl[kNPadRows];               //! cluster Y
  float fArrZCl[kNPadRows];               //! cluster Z
  float fArrYTr[kNPadRows];               //! ref track Y
  float fArrZTr[kNPadRows];               //! ref tracz Z
  float fArrTgSlp[kNPadRows];             //! track inclination at padrow
  int   fArrSectID[kNPadRows];            //! cluster sector id 
  //
  static AliTPCDcalibRes* fgUsedInstance; //! interface instance to use for parameterization
  //
  static const float kSecDPhi;
  static const float kMaxQ2Pt;
  //  static const float kMaxTgSlp;
  //  static const float kMaxResid; // max allowed residual  
  static const float kMinX;   // min X to cover
  static const float kMaxX;   // max X to cover
  static const float kMaxZ2X;   // max z/x
  static const float kZLim;   // endcap position
  static const char* kLocalResFileName;
  static const char* kClosureTestFileName;
  static const char* kStatOut;
  static const char* kResOut;
  static const char* kDriftFileName;
  static const float kDeadZone;  // dead zone on sector edges in cm
  static const float kInvalidR;  // to signal invalid R
  static const float kInvalidRes; // to signal invalid residual
  static const ULong64_t kMByte;
  static const Float_t kZeroK; // zero kernel weight

  static const char* kVoxName[];
  static const char* kResName[];
  
  static const Float_t kTPCRowX[]; // X of the pad-row
  static const Float_t kTPCRowDX[]; // pitch in X

  ClassDef(AliTPCDcalibRes,1);
};

//________________________________________________________________
inline Int_t AliTPCDcalibRes::GetXBinExact(float x) 
{
  // convert X to bin ID, following pad row widths
  if (fUniformBins[kVoxX]) {
    int ix = (x-kMinX)*fDXI;
    return (ix<0 || ix>=fNXBins) ? -2 : ix;
  }
  else return GetRowID(x);
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetY2X(int ix, int iy)
{
  // get Y2X bin center for ix,iy bin
  return (0.5f+iy)*fDY2X[ix] - fMaxY2X[ix];
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetY2XLow(int ix, int iy)
{
  // get Y2X bin low edge for ix,iy bin
  return iy*fDY2X[ix] - fMaxY2X[ix];
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetDY2X(int ix)
{
  // get Y2X bin size value for ix bin
  return fDY2X[ix];
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetDY2XI(int ix)
{
  // get Y2X inverse bin size  for ix bin
  return fDY2XI[ix];
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetX(int i)
{
  // low edge of i-th X bin
  return (fUniformBins[kVoxX]) ? kMinX+(0.5+i)*fDX : kTPCRowX[i];
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetXLow(int i)
{
  // low edge of i-th X bin
  return fUniformBins[kVoxX] ? kMinX+i*fDX : kTPCRowX[i] - 0.5*kTPCRowDX[i];
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetDX(int i)
{
  // width of i-th X bin
  return fUniformBins[kVoxX] ? fDX : kTPCRowDX[i];
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetDXI(int i)
{
  // inverse width of i-th X bin
  return (fUniformBins[kVoxX]) ? fDXI : 1.f/kTPCRowDX[i];
}

//________________________________________________________________
inline Int_t AliTPCDcalibRes::GetY2XBinExact(float y2x, int ix) 
{
  // get exact y2x bin at given x range
  float bf = ( y2x + fMaxY2X[ix] ) * GetDY2XI(ix);
  if (bf<0) return -1;
  else if (bf>=fNY2XBins) return fNY2XBins;
  return int(bf);
}

//________________________________________________________________
inline Int_t AliTPCDcalibRes::GetY2XBin(float y2x, int ix) 
{
  // get closest y2x bin at given x range
  int bf = ( y2x + fMaxY2X[ix] ) * GetDY2XI(ix);
  if (bf<0) bf = 0;
  else if (bf>=fNY2XBins) fNY2XBins-1;
  return bf;
}

//________________________________________________________________
inline Int_t AliTPCDcalibRes::GetZ2XBinExact(float z2x)
{
  // get exact z2x bin at given x range (z2x is positive for clusters not changing the side)
  float bz = z2x*GetDZ2XI();
  if (bz>=fNZ2XBins) return -1;
  if (bz<0) bz = 0; // to account for clusters which moved to wrong side
  return int(bz);
}

//________________________________________________________________
inline Int_t AliTPCDcalibRes::GetZ2XBin(float z2x) 
{
  // get closest z2x bin (z2x is positive for clusters not changing the side)
  int bz = z2x*GetDZ2XI();
  if (bz<0) bz = 0; // to account for clusters which moved to wrong side
  return bz<fNZ2XBins ? bz : fNZ2XBins-1;
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetZ2X(int iz)
{
  // get Z2X bin center for iz, !! always positive
  return (0.5f+iz)*GetDZ2X();
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetZ2XLow(int iz)
{
  // get Z2X bin low edge for iz !! bin positive
  return iz*GetDZ2X();
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetDZ2X()
{
  // get Z2X bin size value
  return fDZ2X;
}

//________________________________________________________________
inline Float_t AliTPCDcalibRes::GetDZ2XI()
{
  // get Z2X inverse bin size
  return fDZ2XI;
}

//_____________________________________
inline void AliTPCDcalibRes::FindVoxel(float x, float y2x, float z2x, int &ix, int &ip, int &iz)
{
  // calculate voxel center sector coordinates (wrt sector)
  ix = GetXBin(x);
  ip = GetY2XBin(y2x,ix);
  iz = GetZ2XBin(z2x);
  //
}

//_____________________________________
inline void AliTPCDcalibRes::FindVoxel(float x, float y2x, float z2x, UChar_t &ix, UChar_t &ip, UChar_t &iz)
{
  // calculate voxel center sector coordinates (wrt sector)
  ix = GetXBin(x);
  ip = GetY2XBin(y2x,ix);
  iz = GetZ2XBin(z2x);
  //
}

//_____________________________________
inline void AliTPCDcalibRes::GetVoxelCoordinates(int isec, int ix, int ip, int iz, float &x, float &p, float &z)
{
  // calculate voxel center sector coordinates (wrt sector)
  x = GetX(ix);
  p = GetY2X(ix,ip);
  z = GetZ2X(iz);
  if (isec>=kNSect) z = -z;
}


//_____________________________________________________
inline Long64_t AliTPCDcalibRes::GetBin2Fill(const UChar_t binVox[kVoxDim], UShort_t bVal)  const
{
  // TH4 bin calculation, bval is the last dimention binID
  ULong64_t binToFill = bVal+1; // 0 bin is undeflow
  for (int id=kVoxDim;id--;) binToFill += fNBProdSt[id]*(1+binVox[id]);
  return binToFill;
}

//_____________________________________________________
inline UShort_t AliTPCDcalibRes::GetVoxGBin(const UChar_t bvox[kVoxDim]) const
{
  // index of geometrix voxel 
  int binToFill = bvox[kVoxDim-1];
  for (int id=kVoxDim-1;id--;) binToFill += fNBProdSectG[id]*bvox[id];
  return binToFill;
}

//_____________________________________________________
inline UShort_t AliTPCDcalibRes::GetVoxGBin(int ix,int ip,int iz) const
{
  // index of geometrix voxel 
  UChar_t bvox[kVoxDim];
  bvox[kVoxX] = ix; bvox[kVoxF] = ip; bvox[kVoxZ] = iz;
  return GetVoxGBin(bvox);
}

//_____________________________________________________
inline void AliTPCDcalibRes::GBin2Vox(UShort_t gbin, UChar_t bvox[kVoxDim]) const
{
  // index to geometrix voxel 
  bvox[kVoxDim-1] = gbin%fNBins[kVoxDim-1];
  for (int id=kVoxDim-1;id--;) bvox[id] = (gbin/fNBProdSectG[id])%fNBins[id];
}


#endif
