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

class AliTPCDcalibRes: public TNamed
{
 public:
  enum {kEpanechnikovKernel, kGaussianKernel};  // defined kernels
  enum {kNSect=18,kNSect2=2*kNSect,kNROC=4*kNSect,kNPadRows=159, kNRowIROC=63, kNRowOROC1=64, kNRowOROC2=32};
  enum {kNQBins=4};
  enum {kAlignmentBugFixedBit = AliTPCcalibAlignInterpolation::kAlignmentBugFixedBit};
  enum {kExtractMode, kClosureTestMode};
  //
  // the voxels are defined in following space
  enum {kVoxQ,   // tg of track inclination wrt pad row, each voxel has its own range according to F,X,Y
	kVoxF,   // y/x in sector coordinates
	kVoxX,   // sector X coordinate
	kVoxZ,   // Z/X sector coordinates
	kVoxV,   // variable within the voxel (delta, stat, etc): last dimension of all THn histos
	kVoxHDim, kVoxDim=kVoxHDim-1};

  enum {kResX,kResY,kResZ,kResD,kResDim,kResDimG=kResDim-1}; // output dimensions

  // content of processed voxel
  enum {kEstNorm,kEstMean,kEstSig,kEstMax,  // statistics
	kEstMeanL,kEstMeanEL, kEstSigL,kEstSigEL, // LTM
	kEstNormG,kEstMeanG,kEstMeanEG,kEstSigG,kEstSigEG,kEstChi2G, // gaussian fit
	kNEstPar};

  struct dts_t {  // struct for basic residual
    //    UChar_t dy;   // Y residual
    //    UChar_t dz;   // Z residual
    Double32_t dy; //[-kMaxResid,kMaxResid,12] 
    Double32_t dz; //[-kMaxResid,kMaxResid,12] 
    UChar_t bvox[kVoxDim]; // voxel bin info: kVoxQ,kVoxF,kVoxX,kVoxZ
  };

 // structure for closure test residuals
 struct dtc_t
 {
   Int_t   t;    // time stamp
   Float_t dyR;   // raw Y residual
   Float_t dzR;   // raw Z residual
   Float_t dyC;   // corrected Y residual
   Float_t dzC;   // corrected Z residual
   Float_t q2pt;
   Float_t tgLam;
   Float_t tgSlp;
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
    Float_t stat[kVoxHDim];  // statistics: averages (weigted over Q bins) of each voxel dimension + entries
    UChar_t bvox[kVoxDim];   // voxel identifier, here the bvox[0] shows number of Q bins used for Y
    UChar_t bsec;            // sector ID (0-35)
    UChar_t smooth;          // smoother flag
    //
    bres_t() {memset(this,0,sizeof(bres_t));}
  };
  
  struct bstat_t {           // stat info on the voxel
    Float_t stat[kVoxHDim];  // statistics: averages of each voxel dimension + entries
    Float_t distY[kNEstPar]; // distortion estimators for Y
    Float_t distZ[kNEstPar]; // distortion estimators for Z
    UChar_t bvox[kVoxDim];   // voxel identifier
    UChar_t bsec;            // sector ID (0-35)
  };

  struct voxDef_t {          // to dumpe the voxel definition (within sector)
    UChar_t bvox[kVoxDim];   // voxel identifier
    Float_t vmin[kVoxDim];   // min boundary
    Float_t vmax[kVoxDim];   // max boundary
  };


public:

  AliTPCDcalibRes(int run=0,Long64_t tmin=0,Long64_t tmax=9999999999,const char* resList=0);
  virtual ~AliTPCDcalibRes();
  
  void ProcessFromDeltaTrees();
  void ProcessFromLocalBinnedTrees();
  void ProcessFromStatTree();
  void Save(const char* name=0);

  void Init();
  void CollectData(int mode = kExtractMode);
  void FillLocalResidualsTrees();
  void FillCorrectedResiduals();
  void ClosureTest();
  void CreateLocalResidualsTrees(int mode);
  void ProcessResiduals();
  void ProcessDispersions();
  void WriteVoxelDefinitions();
  void ProcessSectorResiduals(int is, bstat_t &voxStat);
  void ProcessSectorDispersions(int is);
  void ExtractVoxelDispersion(int is, const TNDArrayT<short>* harrYC, float maxGChi2=5);
  void ExtractVoxelData(bstat_t &stat,const TNDArrayT<short>* harrY,
			const TNDArrayT<short>* harrZ,const TNDArrayT<float>* harrStat);
  void ExtractDistortionsData(TH1F* histo, float est[kNEstPar], const UChar_t vox[kVoxDim], float minNorm=5.f, float fracLTM=0.7f);

  void InitGeom();
  THnF* CreateVoxelStatHisto(int sect);
  THn* CreateSectorResidualsHisto(int sect, int nbDelta,float range, const char* pref);

  void    LoadVDrift();
  Float_t GetDriftCorrection(float z, float x, float phi, int rocID);
  Float_t tgpXY(float x, float y, float q2p, float bz);
  
  void WriteStatHistos();
  void LoadStatHistos();
  void WriteResTree();

  Float_t ExtractResidualHisto(const TNDArrayT<short>* harr, const Long64_t bprod[kVoxHDim], 
			       const UChar_t vox[kVoxDim], TH1F* dest);

  Float_t ExtractResidualHisto(const TNDArrayT<short>* harr, const Long64_t bprod[kVoxHDim], 
			       const UChar_t voxMin[kVoxDim], const UChar_t voxMax[kVoxDim], TH1F* dest);

  TH1F* ExtractResidualHisto(int htype, int sect, const UChar_t vox[kVoxDim]);
  TH1F* ExtractResidualHisto(int htype, int sect, const UChar_t vox[kVoxDim], const UChar_t vox1[kVoxDim]);
  void  ExtractXYZDistortions();
  Bool_t ExtractVoxelXYZDistortions(const bstat_t voxIQ[kNQBins], bres_t &res, 
				    int minStat=20, float maxGChi2=5, int minYBinsOK=3);
  void  FixAlignmentBug(int sect, float q2pt, float bz, float& alp, float& x, float &z, float &deltaY, float &deltaZ);

  Bool_t ValidateTrack();
  Bool_t CompareToHelix(float *resHelixY, float *resHelixZ);

  int    CheckResiduals(Bool_t* kill,float &rmsLongMA);


//------------------------------------ misc. stat. methods

  void    FitCircle(int np, const float* x, const float* y, 
		    double &xc, double &yc, double &r, float* dy=0);
  void    DiffToMA(int np, const float *y, const int winLR, float* diffMA);
  int     DiffToLocLine(int np, const float* x, const float *y, const int nVoisin, float *diffY);
  int     DiffToMedLine(int np, const float* x, const float *y, const int nVoisin, float *diffY);
  float   RoFunc(int np, const float* x, const float* y, float b, float &aa);
  Float_t SelKthMin(int k, int np, float* arr);
  void    medFit(int np, const float* x, const float* y, float &a, float &b, float delI=0.f);
  Bool_t  FitPoly2(const float* x,const float* y, const float* w, int np, float *res, float *err);
  Bool_t  FitPoly1(const float* x,const float* y, const float* w, int np, float *res, float *err);
  Bool_t  GetTruncNormMuSig(double a, double b, double &mean, double &sig);
  void    TruncNormMod(double a, double b, double mu0, double sig0, double &muCf, double &sigCf);
  Double_t GetLogL(TH1F* histo, int bin0, int bin1, double &mu, double &sig, double &logL0);

//------------------------------------

  void    FillHoles(int isect, bres_t *sectData, const int fNBProdSectG[2], int minGoodPoints); // obsolete
  

  Int_t   Smooth0(int isect);
  Bool_t  GetSmoothEstimate(int isect, float x, float p, float z, int which, float *res, float *deriv=0);
  Bool_t  GetSmoothEstimateDim(int isect, float x, float p, float z, int dim, float &res, float *deriv=0);
  void    SetKernelType(int tp=kEpanechnikovKernel, float bwX=2.5, float bwP=2.5, float bwZ=2.1, 
	                float scX=1.f,float scP=1.f,float scZ=1.f);
  
  void    CreateCorrectionObject();
  void    InitBinning();
  Int_t   GetXBin(float x);
  Int_t   GetRowID(float x);
  
  Bool_t  FindVoxelBin(int sectID, float tgSlp, float q2pt, float x, float y, float z, UChar_t bin[kVoxHDim],float voxVars[kVoxHDim]);
  
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

  //  Int_t  GetQBin(float tgp, int binX, int binY);
  Int_t    GetQBin(float tgp);
  Long64_t GetBin2Fill(const Long64_t bprod[kVoxHDim],const UChar_t binVox[kVoxDim], UShort_t bVal);
  Int_t    GetVoxGBin(int ix, int ip, int iz);
  Int_t    GetVoxGBin(UChar_t bvox[kVoxDim]);

  //
  void     SetMaxDY(float m=6.0);
  void     SetMaxDZ(float m=6.0);
  void     SetRun(int run)                       {fRun = run;}
  void     SetTMinMax(Long64_t tmin=0, Long64_t tmax=9999999999) {fTMin=tmin; fTMax=tmax;}
  void     SetMaxQ2Pt(float v=3.0)               {fMaxQ2Pt = v;}
  void     SetMidQ2Pt(float v=1.22)              {fMidQ2Pt = v;}
  void     SetNXBins(int n=kNPadRows)            {fNXBins = n;}
  void     SetNY2XBins(int n=15)                 {fNY2XBins = n;}
  void     SetNZ2XBins(int n=10)                 {fNZ2XBins = n;}
  void     SetNDeltaBinsY(int n=120)             {fNDeltaYBins = n;}
  void     SetNDeltaBinsZ(int n=120)             {fNDeltaZBins = n;}
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
  void     SetMinNClusters(int n=30)             {fMinNCl = n;}
  void     SetNVoisinMA(int n=3)                 {fNVoisinMA = n;}
  void     SetNVoisinMALong(int n=15)            {fNVoisinMALong = n;}
  void     SetMaxDevYHelix(float d=0.3)          {fMaxDevYHelix = d;}
  void     SetMaxDevZHelix(float d=0.3)          {fMaxDevZHelix = d;}
  void     SetMaxStdDevMA(float v=25.0)          {fMaxStdDevMA = v;}
  void     SetMaxRMSLong(float v=0.8)            {fMaxRMSLong = v;}
  void     SetMaxRejFrac(float v=0.15)           {fMaxRejFrac = v;}
  void     SetFilterOutliers(Bool_t v=kTRUE)     {fFilterOutliers = v;}

  Float_t  GetMaxDY()                       const {return fMaxDY;}
  Float_t  GetMaxDZ()                       const {return fMaxDZ;}
  Int_t    GetRun()                         const {return fRun;}
  Long64_t GetTMin()                        const {return fTMin;}
  Long64_t GetTMax()                        const {return fTMax;}  
  Float_t  GetMaxQ2Pt()                     const {return fMaxQ2Pt;}
  Float_t  GetMidQ2Pt()                     const {return fMidQ2Pt;}
  Int_t    GetNXBins()                      const {return fNXBins;}
  Int_t    GetNY2XBins()                    const {return fNY2XBins;}
  Int_t    GetNZ2XBins()                    const {return fNZ2XBins;}
  Int_t    GetNDeltaBinsY()                 const {return fNDeltaYBins;}
  Int_t    GetNDeltaBinsZ()                 const {return fNDeltaZBins;}
  Int_t    GetMaxTracks()                   const {return fMaxTracks;}
  Int_t    GetCacheInput()                  const {return fCacheInp;}
  Int_t    GetCacheLearnSize()              const {return fLearnSize;}
  Int_t    GetNPrimTrackCuts()              const {return fNPrimTracksCut;}
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
  Int_t    fNPCheb[3][2];                           // cheb. nodes per slice

  Float_t  fChebPrecD[3];                           // nominal precision per output dimension
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


  // -------------------------------Binning
  Float_t  fMaxDY;   // max residual in Y
  Float_t  fMaxDZ;   // max residual in Z
  Float_t  fMaxQ2Pt; // max |q/pt|
  Float_t  fMidQ2Pt; // middle |q/pt| for slopes binning 
  Int_t    fNY2XBins;    // y/x bins per sector
  Int_t    fNZ2XBins;    // z/x bins per sector
  Int_t    fNXBins;      // n bins in radial dim.
  Int_t    fNXYBinsProd; // nx*ny bins
  Int_t    fNDeltaYBins; // n bins in Y residual space
  Int_t    fNDeltaZBins; // n bins in Z residual space
  Bool_t   fUniformBins[kVoxDim]; // uniform binning? Currently only X may be non-uniform (per pad-row)


  Float_t  fDZ2X;            // Z2X bin size
  Float_t  fDX;            // X bin size
  Float_t  fDZ2XI;           // inverse Z2X bin size 
  Float_t  fDXI;           // inverse X bin size 
  Float_t  fDeltaYbinI;    // inverse deltaY bin size
  Float_t  fDeltaZbinI;    // inverse deltaZ bin size

  Int_t    fNGVoxPerSector; // total number of geometrical voxels per sector (excluding Q binning)

  Float_t  *fMaxY2X;        //[fNXBins] max Y/X at each X bin, account for dead zones
  Float_t  *fDY2X;          //[fNXBins] Y/X bin size at given X bin
  Float_t  *fDY2XI;         //[fNXBins] inverse of Y/X bin size at given X bin
  // this is obsolete: we bin in q/pt, but convert it on the fly to tgSlp
  //  Float_t  *fBinMinQ;       //[fNXYBinsProd] min value of tg(inclination) at given X,Y bin
  //  Float_t  *fBinDQ;         //[fNXYBinsProd] tg(inclination) bin size at given X,Y bin
  //  Float_t  *fBinDQI;        //[fNXYBinsProd] inverse of tg(inclination) bin size at given X,Y bin
  Float_t  fQ2PTBound[kNQBins+1]; // q2pt bins boundaries

  Long64_t fNBProdSt[kVoxHDim]; // aux arrays for fast bin calculation
  Long64_t fNBProdDY[kVoxHDim];
  Long64_t fNBProdDZ[kVoxHDim];
  Int_t    fNBProdSectG[2];   // aux info for fast bin index calculation in geom voxel space


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
  static const float kMaxResid; // max allowed residual  
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
  static const char* kEstName[];
  
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
inline Int_t AliTPCDcalibRes::GetQBin(float q2pt)
{
  // get binID in track Q variable (tg of inclination) for given X,Y bin
  //
  float q2ptA = TMath::Abs(q2pt);
  if (q2ptA>=fMaxQ2Pt) return -1;
  for (int bin=kNQBins;bin--;) if (q2pt>fQ2PTBound[bin]) return bin;
}

//_____________________________________________________

/*
inline Int_t AliTPCDcalibRes::GetQBin(float tgp, int binX, int binY)
{
  // get binID in track Q variable (tg of inclination) for given X,Y bin
  //
  int id = binX*fNY2XBins+binY;
  float bf = (tgp - fBinMinQ[id])*fBinDQI[id];
  if (bf<0 || bf>kNQBins) return -1;
  return int(bf);
}
*/

//_____________________________________________________
inline Long64_t AliTPCDcalibRes::GetBin2Fill(const Long64_t bprod[kVoxHDim],
							const UChar_t binVox[kVoxDim], UShort_t bVal) 
{
  // TH5 bin calculation, bval is the last dimention binID
  ULong64_t binToFill = bVal+1; // 0 bin is undeflow
  for (int id=kVoxDim;id--;) binToFill += bprod[id]*(1+binVox[id]);
  return binToFill;
}

//_____________________________________________________
inline Int_t AliTPCDcalibRes::GetVoxGBin(int ix, int ip, int iz) 
{
  // index of geometrix voxel (no Q info)
  return iz+fNBProdSectG[1]*ip+fNBProdSectG[0]*ix;
}

//_____________________________________________________
inline Int_t AliTPCDcalibRes::GetVoxGBin(UChar_t bvox[kVoxDim]) 
{
  // index of geometrix voxel (no Q info)
  return bvox[kVoxZ]+fNBProdSectG[1]*bvox[kVoxF]+fNBProdSectG[0]*bvox[kVoxX];
}

//_____________________________________________________
inline void AliTPCDcalibRes::SetMaxDY(float v)
{
  // set max accepted residual
  if (v>kMaxResid) {
    AliWarningF("Redefining %.2f to max allowed %.2f",v,kMaxResid);
    v = kMaxResid;
  }
  fMaxDY = v;
}

//_____________________________________________________
inline void AliTPCDcalibRes::SetMaxDZ(float v)
{
  // set max accepted residual
  if (v>kMaxResid) {
    AliWarningF("Redefining %.2f to max allowed %.2f",v,kMaxResid);
    v = kMaxResid;
  }
  fMaxDZ = v;
}


#endif
