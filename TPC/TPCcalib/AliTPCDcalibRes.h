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
#include <TH2F.h>
#include <TProfile2D.h>
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
#include "AliTPCChebDist.h"


class AliTPCDcalibRes: public TNamed
{
 public:
  enum {kEpanechnikovKernel, kGaussianKernel};  // defined kernels
  enum {kNSect=18,kNSect2=2*kNSect,kNROC=4*kNSect,kNPadRows=159, kNRowIROC=63, kNRowOROC1=64, kNRowOROC2=32};
  enum {kAlignmentBugFixedBit = AliTPCcalibAlignInterpolation::kAlignmentBugFixedBit};
  enum {kVDriftCalibMode, kDistExtractMode, kDistClosureTestMode,kNCollectModes};
  enum {kDistDone=BIT(0),kDispDone=BIT(1),kSmoothDone=BIT(2),kMasked=BIT(7)};
  enum {kVDWithTOFBC=BIT(14),kVDWithoutTOFBC=BIT(15)};
  enum {kUseTRDonly,kUseTOFonly,kUseITSonly,kUseTRDorTOF,kNExtDetComb}; // which points to use
  enum {kSmtLinDim=4, kMaxSmtDim=7}; // max size of matrix for smoothing, for pol1 and pol2 options
  enum {kCtrITS,kCtrTRD,kCtrTOF,kCtrBC0,kCtrNtr,kCtrNbr}; // control branches for test stat
  enum {kLargeTimeDiff=1000*3600,kMaxOnStack=10000};  // misc consts

  // the voxels are defined in following space
  enum {kVoxZ,   // Z/X sector coordinates
	kVoxF,   // y/x in sector coordinates
	kVoxX,   // sector X coordinate
	kVoxV,   // variable within the voxel (delta, stat, etc): last dimension of all THn histos
	kVoxHDim, kVoxDim=kVoxHDim-1};

  enum {kResX,kResY,kResZ,kResD,kResDim,kResDimG=kResDim-1}; // output dimensions
  //
  struct dtv_t {  // struct for vdrift calibration loca residual
    Char_t  side;      // +/- for A/C
    Double32_t dz;     //[-40.,40.,15]   // +-kMaxResidZVD
    Double32_t drift;  //[-260.,260.,14] drift distance
    Double32_t ylab;   //[-260.,260.,14] y coordinate
    Double32_t r;      //[80.,250.,7] radius
    Int_t   t;         // time stamp
    //
    dtv_t() {memset(this,0,sizeof(dtv_t));}
  };
  //
  struct dts_t {  // struct for basic local residual
    Double32_t dy; //[-20.,20.,15] // [-kMaxResid,kMaxResid,14]
    Double32_t dz; //[-20.,20.,15] // [-kMaxResid,kMaxResid,14]
    Double32_t tgSlp; //[-2,2,14]  //[kMaxTgSlp,kMaxTgSlp,14]
    UChar_t bvox[kVoxDim]; // voxel bin info: VoxF,kVoxX,kVoxZ
    //
    dts_t() {memset(this,0,sizeof(dts_t));}
  };

  // structure for closure test residuals
  struct dtc_t
  {
    Int_t t;        // time stamp
    Double32_t dyR; //[-15.,15.,14] 
    Double32_t dzR; //[-15.,15.,14]
    Double32_t dyC; //[-15.,15.,14] 
    Double32_t dzC; //[-15.,15.,14]
    Double32_t q2pt;//[-3,3,14]
    Double32_t tgLam;//[-2.,2.,14]
    Double32_t tgSlp;//[-2,2,14] 
    Float_t x;  //[80,160,14]
    Float_t y; //[-50,50,14]
    Float_t z; //[-250,250,14]
    UChar_t bvox[kVoxDim]; // voxel bin info: kVoxQ,kVoxF,kVoxX,kVoxZ
    //
    dtc_t() {memset(this,0,sizeof(dtc_t));}
  };

  struct bres_t  {
    Float_t D[kResDim];      // values of extracted distortions
    Float_t E[kResDim];      // their errors
    Float_t DS[kResDim];     // smoothed residual
    Float_t DC[kResDim];     // Cheb parameterized residual
    Float_t EXYCorr;         // correlation between extracted X and Y
    Float_t dYSigMAD;        // MAD estimator of dY sigma (dispersion after slope removal)
    Float_t dZSigLTM;        // Z sigma from unbinned LTM estimator
    Float_t stat[kVoxHDim];  // statistics: averages of each voxel dimension + entries
    UChar_t bvox[kVoxDim];   // voxel identifier, here the bvox[0] shows number of Q bins used for Y
    UChar_t bsec;            // sector ID (0-35)
    UChar_t flags;           // status flag
    //
    bres_t() {memset(this,0,sizeof(bres_t));}
  };
 
  struct delta_t { // structure to organize the input from delta trees
    TVectorF *vecDYTRD,*vecDZTRD,*vecDYITS,*vecDZITS,*vecDYTOF,*vecDZTOF,*vecZ,*vecR,*vecSec,*vecPhi;
    AliExternalTrackParam* param;
    Double_t tofBC;
    Int_t nPrimTracks;
    Int_t timeStamp;
    UShort_t npValid;
    Char_t trdOK,tofOK,itsOK;
    //
    delta_t() {memset(this,0,sizeof(delta_t));}
  };
  
  struct dcTst_t {  // struct for distortion/correction check
    Double32_t xyz[3];
    Double32_t dxyz[4]; // [-20.,20.,15]
    Double32_t cxyz[4]; // [-20.,20.,15]
    UChar_t bsec;          // sector
    UChar_t bvox[kVoxDim]; // voxel bin info: VoxF,kVoxX,kVoxZ
    //
    dcTst_t() {memset(this,0,sizeof(dcTst_t));}
  };

 public:

  AliTPCDcalibRes(int run=0,Long64_t tmin=0,Long64_t tmax=9999999999,const char* resList=0);
  virtual ~AliTPCDcalibRes();
  
  void CalibrateVDrift();
  void FitDrift();
  TProfile2D* GetHistoVDTimeInt()    const {return fHVDTimeInt;}
  TH2F*       GetHistoVDTimeCorr()   const {return fHVDTimeCorr;}
  Int_t GetDeltaTVDrift()  const {return fDeltaTVD;}
  Int_t GetSigmaTVDrift()  const {return fSigmaTVD;}
  void  SetDeltaTVDrift(int v=120) {fDeltaTVD = v;}
  void  SetSigmaTVDrift(int v=600) {fSigmaTVD = v;}
  //
  void ProcessFromDeltaTrees();
  void ProcessFromLocalBinnedTrees();
  void ReProcessFromResVoxTree(const char* resTreeFile, Bool_t backup=kTRUE);
  void Save(const char* name=0);
  void SaveFitDrift();

  TTree* InitDeltaFile(const char* name, Bool_t connect=kTRUE, const char* treeName="delta");
  Bool_t EstimateChunkStatistics();
  Bool_t CollectDataStatistics();
  Bool_t EnoughStatForVDrift(int tstamp=-1,float maxHolesFrac=0.05);
  Int_t  ParseInputList();
  void CloseTreeFile(TTree* dtree);
  void Init();
  void CollectData(const int mode = kDistExtractMode);
  void FillLocalResidualsTrees();
  void FillCorrectedResiduals();
  void FillDriftResidualsTrees();
  void ClosureTest();
  void CreateLocalResidualsTrees(int mode);
  void CloseLocalResidualsTrees(int mode);
  void ProcessResiduals();
  void ReProcessResiduals();
  void ProcessDispersions();
  void ProcessSectorResiduals(int is);
  void ReProcessSectorResiduals(int is);
  void ProcessSectorDispersions(int is);
  void ProcessVoxelResiduals(int np, float* tg, float *dy, float *dz, bres_t& voxRes);
  void ProcessVoxelDispersions(int np, const float* tg, float *dy, bres_t& voxRes);
  Int_t ValidateVoxels(int isect);
  //
  void InitFieldGeom(Bool_t field=kTRUE,Bool_t geom=kTRUE);
  THnF* CreateVoxelStatHisto(int sect);
  void    LoadVDrift();
  void    MakeVDriftOCDB(TString targetOCDBstorage);
  Float_t GetDriftCorrection(float z, float x, float phi, int rocID);
  Float_t tgpXY(float x, float y, float q2p, float bz);

  TVectorD*     GetVDriftParam() const {return (TVectorD*)fVDriftParam;}
  TGraphErrors* GetVDriftGraph() const {return (TGraphErrors*)fVDriftGraph;}

  void WriteStatHistos();
  void LoadStatHistos();
  void WriteResTree();
  void WriteDistCorTestTree(int nx=250, int nphi2sect=30,int nzside=10);
  Bool_t LoadDataFromResTree(const char* resTreeFile);

  void  FixAlignmentBug(int sect, float q2pt, float bz, float& alp, float& x, float &z, float &deltaY, float &deltaZ);

  Bool_t ValidateTrack();
  Bool_t CompareToHelix(float *resHelixY, float *resHelixZ);

  int    CheckResiduals(Bool_t* mask,float &rmsLongMA);

  const char* GetVoxResFileName() const {return Form("%sTree.root",kResOut);}

  //------------------------------------ misc. stat. methods
  static Float_t FitPoly1Robust(int np, float* x, float* y, float* res, float* err, float ltmCut);
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
  static Bool_t  FitPoly2(const float* x,const float* y, const float* w, int np, float *res, float *err=0);
  static Bool_t  FitPoly1(const float* x,const float* y, const float* w, int np, float *res, float *err=0);
  static Bool_t  GetTruncNormMuSig(double a, double b, double &mean, double &sig);
  static void    TruncNormMod(double a, double b, double mu0, double sig0, double &muCf, double &sigCf);
  static Double_t GetLogL(TH1F* histo, int bin0, int bin1, double &mu, double &sig, double &logL0);
  static Bool_t  GetSmooth1D(double xQuery,double valerr[2],int np,const double* x,const double* y,const double* err,
			     double wKernel,int kType=kGaussianKernel,Bool_t usePoly2=kFALSE,
			     Bool_t xIncreasing=kTRUE, Float_t maxD2Range=3.0);
  static Bool_t  GradCorrCheb(const AliTPCChebCorr* cheb, int sect36, float x, float y, float z, 
			      float val[AliTPCDcalibRes::kResDim],
			      float grad[AliTPCDcalibRes::kResDimG][AliTPCDcalibRes::kResDim]);
  //------------------------------------
  

  Int_t   Smooth0(int isect);
  void    BringToActiveBoundary(int sect36, float xyz[3]) const;
  Bool_t  GetSmoothEstimate(int isect, float x, float p, float z, int which, float *res, float *deriv=0);
  Bool_t  GradValSmooth(int sect36, float x, float y, float z,float val[AliTPCDcalibRes::kResDim],
			float grad[AliTPCDcalibRes::kResDimG][AliTPCDcalibRes::kResDim], Bool_t bringToBoundary=kTRUE);
  Bool_t  InvertCorrection(int sect36,const float vecIn[AliTPCDcalibRes::kResDimG],float vecOut[AliTPCDcalibRes::kResDim],
			   int maxIt=15,float minDelta=50e-4);
  void    SetKernelType(int tp=kEpanechnikovKernel, float bwX=2.1, float bwP=2.1, float bwZ=1.7,float scX=1.f,float scP=1.f,float scZ=1.f);
  Bool_t  GetSmoothPol2(int i)                              const {return fSmoothPol2[i];}
  void    SetSmoothPol2(int i,Bool_t v=kTRUE)                     {fSmoothPol2[i] = v;}
  //  
  void    CreateCorrectionObject();
  void    CreateDistortionObject();
  void    InitBinning();
  Int_t   GetXBin(float x);
  static  Int_t  GetRowID(float x);
  
  Bool_t  FindVoxelBin(int sectID, float x, float y, float z, UChar_t bin[kVoxHDim],float voxVars[kVoxHDim]);
  TH1F*   GetEventsRateHisto()          const {return fEvRateH;}
  TH1F*   GetTracksRateHisto()          const {return fTracksRate;}
  TH1F*   GetTOFBCHisto()               const {return fTOFBCTestH;}
  void    SetTracksRateHisto(TH1F* h)  {fTracksRate = h;}   // needed for hacks
  TH1F*   ExtractTrackRate() const;
  Float_t GetValidFracXBin(int sect, int bin) const {return fValidFracXBin[sect][bin];}
  Int_t   GetNSmoothingFailed(int sect) const {return fNSmoothingFailedBins[sect];}
  Int_t   GetNTracksUsed()              const {return fNTrSelTot;}
  Int_t   GetNTracksWithOutliers()      const {return fNTrSelTotWO;}
  Int_t   GetMinTrackToUse()            const {return fMinTracksToUse;}
  Int_t   GetMinTracksPerVDBin()        const {return fMinTracksPerVDBin;}
  Int_t   GetNTestTracks()              const {return fNTestTracks;}
  Float_t GetTestStat(int row, int col) const {return fTestStat[row][col];}
  TGraph*  GetLumiGraph() const {return (TGraph*) fLumiCTPGraph;}
  Float_t  GetLumiCOG()   const {return fLumiCOG;}
  //
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
  static Double_t GetKernelWeight(double *u2vec, int np, int kernelType);

  Long64_t    GetBin2Fill(const UChar_t binVox[kVoxDim], UShort_t bVal) const;
  UShort_t    GetVoxGBin(const UChar_t bvox[kVoxDim]) const;
  UShort_t    GetVoxGBin(int ix,int ip,int iz) const;
  void        GBin2Vox(UShort_t gbin, UChar_t bvox[kVoxDim]) const;
  //
  void     SetRun(int run)                       {fRun = run;}
  void     SetTMinMax(Long64_t tmin=0, Long64_t tmax=9999999999);
  void     SetNXBins(int n=kNPadRows)            {fNXBins = n;}
  void     SetNY2XBins(int n=15)                 {fNY2XBins = n;}
  void     SetNZ2XBins(int n=5)                  {fNZ2XBins = n;}
  void     SetMaxTracks(int n=4000000)           {fMaxTracks = n;}
  void     SetNominalTimeBin(int nsec=2400)      {fNominalTimeBin = nsec;}
  void     SetNominalTimeBinPrec(float p=0.3)    {fNominalTimeBinPrec = p;}
  void     SetFixAligmentBug(Bool_t v=kTRUE)     {fFixAlignmentBug = v;}
  void     SetCacheLearnSize(int n=1)            {fLearnSize = n;}
  void     SetCacheInput(Int_t v=100)            {fCacheInp = v;}
  void     SetSwitchCache(Bool_t v=kFALSE)       {fSwitchCache = v;}
  void     SetApplyZt2Zc(Bool_t v=kTRUE)         {fApplyZt2Zc = v;}
  void     SetResidualList(const char* l)        {fResidualList = l;}
  void     SetOCDBPath(const char* l)            {fOCDBPath = l;}
  void     SetUseErrorInSmoothing(Bool_t v=kTRUE) {fUseErrInSmoothing = v;}
  void     SetNPrimTrackCuts(int n=400)          {fNPrimTracksCut = n;}
  void     SetMinTrackToUse(int n=600000)        {fMinTracksToUse = n;}
  void     SetMinTracksPerVDBin(int n=3000)      {fMinTracksPerVDBin = n;}
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
  void     SetMaxFitXErr2(float v=9.0)           {fMaxFitXErr2 = v;}
  void     SetMaxFitXYCorr(float v=0.95)         {fMaxFitXYCorr = v;}
  void     SetLTMCut(float v=0.75)               {fLTMCut = v;}

  Bool_t   GetXBinIgnored(int sect, int bin) const           {return fXBinIgnore[sect].TestBitNumber(bin);}
  void     SetXBinIgnored(int sect, int bin, Bool_t v=kTRUE) {fXBinIgnore[sect].SetBitNumber(bin,v);}
  //
  Float_t  GetMinValidVoxFracDrift()       const {return fMinValidVoxFracDrift;}
  Float_t  GetMaxSigY()                    const {return fMaxSigY;}
  Float_t  GetMaxSigZ()                    const {return fMaxSigZ;}
  Int_t    GetMaxBadXBinsToCover()         const {return fMaxBadXBinsToCover;}
  Int_t    GetMinGoodXBinsToCover()        const {return fMinGoodXBinsToCover;}
  Int_t    GetMaxBadRowsPerSector()        const {return fMaxBadRowsPerSector;}
  
  void     SetMinValidVoxFracDrift(float v=0.7)   {fMinValidVoxFracDrift = v;}
  void     SetMaxSigY(Float_t v=1.1)              {fMaxSigY = v;}
  void     SetMaxSigZ(Float_t v=0.7)              {fMaxSigZ = v;}  
  void     SetMaxBadXBinsToCover(int n=4)         {fMaxBadXBinsToCover = n;}
  void     SetMinGoodXBinsToCover(int n=2)        {fMinGoodXBinsToCover = n;}
  void     SetMaxBadRowsPerSector(float v=0.5)    {fMaxBadRowsPerSector = v;}
  //
  Bool_t   GetUseTOFBC()                   const {return fUseTOFBC;}
  Float_t  GetTOFBCMin()                   const {return fTOFBCMin;}
  Float_t  GetTOFBCMax()                   const {return fTOFBCMax;}
  void     SetUseTOFBC(Bool_t v)                 {fUseTOFBC = v;}
  void     SetTOFBCMin(Float_t v=-25.f)           {fTOFBCMin = v;}
  void     SetTOFBCMax(Float_t v=50.f)           {fTOFBCMax = v;}
  //
  Float_t  GetMaxFitYErr2()                 const {return fMaxFitYErr2;}
  Float_t  GetMaxFitXErr2()                 const {return fMaxFitXErr2;}
  Float_t  GetMaxFitXYCorr()                const {return fMaxFitXYCorr;}
  Float_t  GetLTMCut()                      const {return fLTMCut;}

  Int_t    GetRun()                         const {return fRun;}
  Long64_t GetTMin()                        const {return fTMin;}
  Long64_t GetTMax()                        const {return fTMax;}  
  Long64_t GetTMinGRP()                     const {return fTMinGRP;}
  Long64_t GetTMaxGRP()                     const {return fTMaxGRP;}  
  Long64_t GetTMinCTP()                     const {return fTMinCTP;}
  Long64_t GetTMaxCTP()                     const {return fTMaxCTP;}  
  Int_t    GetNXBins()                      const {return fNXBins;}
  Int_t    GetNY2XBins()                    const {return fNY2XBins;}
  Int_t    GetNZ2XBins()                    const {return fNZ2XBins;}
  Int_t    GetMaxTracks()                   const {return fMaxTracks;}
  Int_t    GetNominalTimeBin()              const {return fNominalTimeBin;}
  Float_t  GetNominalTimeBinPrec()          const {return fNominalTimeBinPrec;}
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
  Int_t    GetExternalDetectors()           const {return fExtDet;}
  void     SetExternalDetectors(int det=kUseTRDonly);

  //
  void     SetChebZSlicePerSide(int n=1)          {fChebZSlicePerSide = n;}
  void     SetChebPhiSlicePerSector(int n=1)      {fChebPhiSlicePerSector = n;}
  Int_t    GetChebZSlicePerSide()           const {return fChebZSlicePerSide;}
  Int_t    GetChebPhiSlicePerSector()       const {return fChebPhiSlicePerSector;}
  Int_t*   GetNPCheb(int dim)               const {return (Int_t*)fNPCheb[dim];}
  Float_t* GetChebPrecD()                   const {return (Float_t*)fChebPrecD;}
  //
  Bool_t   GetFixAlignmentBug()             const {return fFixAlignmentBug;}
  Bool_t   GetSwitchCache()                 const {return fSwitchCache;}
  Bool_t   GetApplyZt2Zc()                  const {return fApplyZt2Zc;}
  Bool_t   GetUseErrorInSmoothing()         const {return fUseErrInSmoothing;}

  Bool_t   GetCreateCorrection()            const {return fCreateCorrection;}
  Bool_t   GetCreateDistortion()            const {return fCreateDistortion;}
  void     SetCreateCorrection(Bool_t v=kTRUE)    {fCreateCorrection = v;}
  void     SetCreateDistortion(Bool_t v=kTRUE)    {fCreateDistortion = v;}
  
  const TString& GetOCDBPath()              const {return fOCDBPath;}
  const TString& GetReisdualList()          const {return fResidualList;}

  const AliTPCChebCorr* GetChebCorrObject() const {return fChebCorr;}
  const AliTPCChebDist* GetChebDistObject() const {return fChebDist;}

  static AliTPCDcalibRes* Load(const char* fname="alitpcdcalibres.root");
  static TTree*           LoadTree(const char* fname="voxelResTree.root",Int_t mStyle=7,Int_t mCol=kBlack);
  static void SetUsedInstance(AliTPCDcalibRes* inst) {fgUsedInstance = inst;}
  static AliTPCDcalibRes* GetUsedInstance()          {return fgUsedInstance;}
  static void SetUsedInstanceExt(AliTPCDcalibRes* inst) {fgUsedInstanceExt = inst;}
  static AliTPCDcalibRes* GetUsedInstanceExt()          {return fgUsedInstanceExt;}
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
  AliTPCChebCorr* fChebCorr;                        // final Chebyshev correction object
  AliTPCChebDist* fChebDist;                        // final Chebyshev distortion object
  // -------------------------------Task defintion
  Int_t    fRun;     // run number 
  Int_t    fExtDet;  // external detectors to use
  Long64_t fTMin;    // time start for timebin
  Long64_t fTMax;    // time stop for timebin
  Long64_t fTMinGRP;    // time start from GRP
  Long64_t fTMaxGRP;    // time stop from GRP
  Long64_t fTMinCTP;    // time start from CTP
  Long64_t fTMaxCTP;    // time stop from CTP
  Int_t    fDeltaTVD;      // vdrift T-bin size
  //
  Int_t    fSigmaTVD;      // vdrift smoothing window
  Int_t    fMaxTracks;     // max tracks to accept
  Int_t    fNominalTimeBin;// nominal time bin
  Float_t  fNominalTimeBinPrec; // allow deviation of time bin from nominal within this limit
  Int_t    fCacheInp;      // input trees cache in MB
  Int_t    fLearnSize;     // event to learn for the cache
  Float_t  fBz;            // B field
  Bool_t   fDeleteSectorTrees; // delete residuals trees once statistics tree is done
  TString  fResidualList;  // text file with list of residuals tree
  TObjArray* fInputChunks;  // list of input files used
  TString  fOCDBPath;      // ocdb path
  // ------------------------------Selection/filtering cuts
  Int_t    fMinTracksPerVDBin;       // min number of tracks per VDrift bin
  Int_t    fMinTracksToUse;          // produce warning if n tracks is too low
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
  Float_t  fTOFBCMin;                // min dTOF cut in ns if validation requested
  Float_t  fTOFBCMax;                // max dTOF cut in ns if validation requested
  Bool_t   fUseTOFBC;                // require TOF BC validation
  Bool_t   fFilterOutliers;          // reject outliers
  Bool_t   fFatalOnMissingDrift;     // if vdrift is needed by absent, produce fatal
  Bool_t   fCreateCorrection;        // request to create Cheb correction
  Bool_t   fCreateDistortion;        // request to create Cheb distortion
  Float_t  fMaxFitYErr2;             // cut on median fit Y err^2
  Float_t  fMaxFitXErr2;             // cut on median fit X err^2
  Float_t  fMaxFitXYCorr;            // cut on max correlation of X,Y errors in median fit
  Float_t  fLTMCut;                  // LTM cut for outliers suppression
  //
  //-------------------------------- voxel validation
  Float_t  fMaxSigY;                 // cut on sigY (after slop removal) MAD estimator
  Float_t  fMaxSigZ;                 // cut on sigZ LTM estimator
  Float_t  fMinValidVoxFracDrift;    // smooth/parameterize only Xbins with fraction of valid voxels above the threshold
  Int_t    fMaxBadXBinsToCover;      // do not extrapolate to more than this number of bad Xbins
  Int_t    fMinGoodXBinsToCover;     // requre at least this amount of consecutive good bins to parameterize
  Float_t  fMaxBadRowsPerSector;     // block whole sector once the fraction of bad Xbins exceeds

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
  Int_t    fKernelType;        // kernel type
  Int_t    fStepKern[kVoxDim]; // N bins to consider with given kernel settings
  Float_t  fKernelWInv[kVoxDim];      // inverse kernel width in bins
  Float_t  fKernelScaleEdge[kVoxDim]; // optional scaling factors for kernel width on the edge
  Bool_t   fSmoothPol2[kVoxDim];      // option for use pol1 or pol2 in each direction (no x-terms)
  // result of last kernel minimization: value and dV/dX,dV/dY,dV/dZ for each dim
  Double_t fLastSmoothingRes[kResDim*kMaxSmtDim];  //! results of last smoothing
  // ------------------------------Selection Stats
  Int_t    fNEvTot;         // total number of events
  Int_t    fNTrSelTot;      // selected tracks
  Int_t    fNTrSelTotWO;    // would be selected w/o outliers rejection
  Int_t    fNReadCallTot;   // read calls from input trees
  Long64_t fNBytesReadTot;  // total bytes read
  Int_t    fNTestTracks;                // number of control tracks in test stat
  Float_t  fEstTracksPerEvent;          // estimate for used tracks per selected event
  Float_t  fTestStat[kCtrNbr][kCtrNbr]; // control statistics (fraction to fNTestTracks) 
  TH1F*    fEvRateH;        // events per time histo
  TH1F*    fTracksRate;     // accepted tracks per second
  TH1F*    fTOFBCTestH;     // TOF BC (ns) for test fNTestTracks tracks
  TProfile2D* fHVDTimeInt;   // histo used for time integrated vdrift fit
  TH2F*    fHVDTimeCorr;     // histo used for vdrift time correction fit
  Float_t  fValidFracXBin[kNSect2][kNPadRows]; // fraction of voxels valid per padrow
  Int_t    fNSmoothingFailedBins[kNSect2];     // number of failed bins/sector, should be 0 to produce parameterization
  TBits    fXBinIgnore[kNSect2];    // flag to ignore Xbin
  Float_t  fLumiCOG;                // COG lumi for timebin
  TGraph*  fLumiCTPGraph;           // lumi graph from CTP
  // ------------------------------VDrift correction
  TVectorD     *fVDriftParam;
  TGraphErrors *fVDriftGraph;  
  Float_t       fCorrTime;   //! 
  // -----------------------------Results of processing
  bres_t *fSectGVoxRes[kNSect2];         //! [fNGVoxPerSector] sectors results for geometric voxel
  TTree* fStatTree;                      //! tree with voxels statistics
  TTree* fTmpTree[kNSect2];              //! IO tree per sector
  TFile* fTmpFile[kNSect2];              //! file for fTmpTree
  THnF*  fStatHist[kNSect2];             //! histos for statistics bins
  TNDArrayT<float> *fArrNDStat[kNSect2]; //! alias arrays for fast access to fStatHist
  //
  // ----------------------------data exchange structures for trees and between routines
  dts_t fDTS;                            //! binned residuals
  dtc_t fDTC;                            //! corrected residuals for closure test
  dtv_t fDTV;                            //! vdrift residuals
  delta_t fDeltaStr;                     //! input from delta tree
  //
  // ---------------------------track data-----------------------------------
  int   fTimeStamp;                       //! time stamp
  int   fNCl;                             //! number of clusters
  float fQ2Pt;                            //! fitted q2pt
  float fTgLam;                           //! fitted tgLambda
  float fArrPhi[kNPadRows];               //! cluster phi
  float fArrDY[kNPadRows];                //! cluster residual Y
  float fArrDZ[kNPadRows];                //! cluster residual Z
  float fArrR[kNPadRows];                 //! cluster R 
  float fArrX[kNPadRows];                 //! cluster X in sector coordinates (row) 
  float fArrYCl[kNPadRows];               //! cluster Y in sector coordinates
  float fArrZCl[kNPadRows];               //! cluster Z
  float fArrYTr[kNPadRows];               //! ref track Y in sector coordinates
  float fArrZTr[kNPadRows];               //! ref tracz Z in sector coordinates
  float fArrTgSlp[kNPadRows];             //! track inclination at padrow
  int   fArrSectID[kNPadRows];            //! cluster sector id 
  //
  static AliTPCDcalibRes* fgUsedInstance; //! interface instance to use for parameterization
  static AliTPCDcalibRes* fgUsedInstanceExt; //! interface to extra instance if average of 2 maps to be used
  //
  static const float kMaxResid;     // max range of distortions, must be <= than the double32_t range of dst_t
  static const float kMaxResidZVD;  // max range of distortions in VDrift calib, must be <= than the double32_t range of dtv_t
  static const float kMaxTgSlp;  // max range of tgSlope, must be <= than the double32_t range of dst_t
  static const float kSecDPhi;
  static const float kMaxQ2Pt;
  static const float kMaxGaussStdDev; // don't evaluate gaussian outside this number of StdDevs
  //  static const float kMaxTgSlp;
  //  static const float kMaxResid; // max allowed residual  
  static const float kMinX;   // min X to cover
  static const float kMaxX;   // max X to cover
  static const float kMaxZ2X;   // max z/x
  static const float kZLim[2];   // endcap positions
  static const char* kDriftResFileName;
  static const char* kLocalResFileName;
  static const char* kStatInfoFileName;
  static const char* kClosureTestFileName;
  static const char* kStatOut;
  static const char* kResOut;
  static const char* kDriftFileName;
  static const float kDeadZone;  // dead zone on sector edges in cm
  static const float kInvalidR;  // to signal invalid R
  static const float kInvalidRes; // to signal invalid residual
  static const ULong64_t kMByte;
  static const Float_t kZeroK; // zero kernel weight
  static const char* kControlBr[];
  static const char* kVoxName[];
  static const char* kResName[];
  static const char* kModeName[];
  static const char* kExtDetName[];
  
  static const Float_t kTPCRowX[]; // X of the pad-row
  static const Float_t kTPCRowDX[]; // pitch in X

  ClassDef(AliTPCDcalibRes,17);
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
  else if (bf>=fNY2XBins) bf = fNY2XBins-1;
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

//_____________________________________________
inline void AliTPCDcalibRes::SetTMinMax(Long64_t tmin, Long64_t tmax) {
  // set min max time
  fTMin = tmin;
  fTMax = tmax;
  if (fTMin>=fTMax) {
    fTMin = 0;
    fTMax = 9999999999;
  } 
}

#endif
