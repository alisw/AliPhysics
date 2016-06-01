#include "AliTPCDcalibRes.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliTriggerRunScalers.h"
#include "AliTriggerScalersRecord.h"
#include "AliDAQ.h"
#include <TKey.h>
#include <TF2.h>

using std::swap;

// this must be standalone f-n, since the signature is important for Chebyshev training
void trainCorr(int row, float* tzLoc, float* corrLoc);

// make sure these branches are always connected in InitDeltaFile
const char* AliTPCDcalibRes::kControlBr[kCtrNbr] = {"itsOK","trdOK","tofOK","tofBC","nPrimTracks"}; 
const char* AliTPCDcalibRes::kVoxName[AliTPCDcalibRes::kVoxHDim] = {"z2x","y2x","x","N"};
const char* AliTPCDcalibRes::kResName[AliTPCDcalibRes::kResDim] = {"dX","dY","dZ","Disp"};
const float  AliTPCDcalibRes::kMaxResid=20.0f;   
const float  AliTPCDcalibRes::kMaxResidZVD=40.0f;   
const float  AliTPCDcalibRes::kMaxTgSlp=2.0;

const float AliTPCDcalibRes::kSecDPhi = 20.f*TMath::DegToRad();
const float AliTPCDcalibRes::kMaxQ2Pt = 3.0f;
const float AliTPCDcalibRes::kMinX = 85.0f;
const float AliTPCDcalibRes::kMaxX = 246.0f;
const float AliTPCDcalibRes::kMaxZ2X = 1.0f;
const float AliTPCDcalibRes::kZLim[2] = {2.49725e+02,2.49698e+02};
const char* AliTPCDcalibRes::kDriftResFileName  = "tmpDriftTree";
const char* AliTPCDcalibRes::kLocalResFileName  = "tmpDeltaSect";
const char* AliTPCDcalibRes::kClosureTestFileName  = "closureTestSect";
const char* AliTPCDcalibRes::kStatOut      = "voxelStat";
const char* AliTPCDcalibRes::kResOut       = "voxelRes";
const char* AliTPCDcalibRes::kDriftFileName= "fitDrift";
const float AliTPCDcalibRes::kDeadZone = 1.5f;
const float AliTPCDcalibRes::kZeroK = 1e-6;
const float AliTPCDcalibRes::kInvalidR = 10.f;
const float AliTPCDcalibRes::kInvalidRes = -900.0f;
const float AliTPCDcalibRes::kMaxGaussStdDev = 5.0f;
const ULong64_t AliTPCDcalibRes::kMByte = 1024LL*1024LL;

const Float_t AliTPCDcalibRes::kTPCRowX[AliTPCDcalibRes::kNPadRows] = { // pad-row center X
  85.225, 85.975, 86.725, 87.475, 88.225, 88.975, 89.725, 90.475, 91.225, 91.975, 92.725, 93.475, 94.225, 94.975, 95.725,
  96.475, 97.225, 97.975, 98.725, 99.475,100.225,100.975,101.725,102.475,103.225,103.975,104.725,105.475,106.225,106.975,
  107.725,108.475,109.225,109.975,110.725,111.475,112.225,112.975,113.725,114.475,115.225,115.975,116.725,117.475,118.225,
  118.975,119.725,120.475,121.225,121.975,122.725,123.475,124.225,124.975,125.725,126.475,127.225,127.975,128.725,129.475,
  130.225,130.975,131.725,135.100,136.100,137.100,138.100,139.100,140.100,141.100,142.100,143.100,144.100,145.100,146.100,
  147.100,148.100,149.100,150.100,151.100,152.100,153.100,154.100,155.100,156.100,157.100,158.100,159.100,160.100,161.100,
  162.100,163.100,164.100,165.100,166.100,167.100,168.100,169.100,170.100,171.100,172.100,173.100,174.100,175.100,176.100,
  177.100,178.100,179.100,180.100,181.100,182.100,183.100,184.100,185.100,186.100,187.100,188.100,189.100,190.100,191.100,
  192.100,193.100,194.100,195.100,196.100,197.100,198.100,199.350,200.850,202.350,203.850,205.350,206.850,208.350,209.850,
  211.350,212.850,214.350,215.850,217.350,218.850,220.350,221.850,223.350,224.850,226.350,227.850,229.350,230.850,232.350,
  233.850,235.350,236.850,238.350,239.850,241.350,242.850,244.350,245.850
};
const Float_t AliTPCDcalibRes::kTPCRowDX[AliTPCDcalibRes::kNPadRows] = { // pad-row pitch in X
  0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,
  0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,
  0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,
  1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
  1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
  1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
  1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,
  1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500
};


AliTPCDcalibRes* AliTPCDcalibRes::fgUsedInstance = 0;


ClassImp(AliTPCDcalibRes)

//________________________________________
AliTPCDcalibRes::AliTPCDcalibRes(int run,Long64_t tmin,Long64_t tmax,const char* resList) : 
  fInitDone(kFALSE)
  ,fUseErrInSmoothing(kTRUE)
  ,fSwitchCache(kFALSE)
  ,fFixAlignmentBug(kTRUE)
  ,fApplyZt2Zc(kTRUE)

  ,fChebZSlicePerSide(1)
  ,fChebPhiSlicePerSector(1)
  ,fChebCorr(0)

  ,fRun(run)
  ,fExtDet(kUseTRDonly)
  ,fTMin(tmin)
  ,fTMax(tmax)
  ,fTMinGRP(0)
  ,fTMaxGRP(0)
  ,fTMinCTP(0)
  ,fTMaxCTP(0)
  ,fDeltaTVD(120)
  ,fSigmaTVD(600)
  ,fMaxTracks(4000000)
  ,fCacheInp(100)
  ,fLearnSize(1)
  ,fBz(0)
  ,fDeleteSectorTrees(kFALSE) // set to true for production
  ,fResidualList(resList)
  ,fInputChunks(0)
  ,fOCDBPath()

  ,fMinTracksToUse(600000)
  ,fMinEntriesVoxel(15)
  ,fNPrimTracksCut(400)
  ,fMinNCl(30)
  ,fMaxDevYHelix(0.3)
  ,fMaxDevZHelix(0.3) // !!! VDrift calib. screas up the Z fit, 0.3 w/o vdrift
  ,fNVoisinMA(3)
  ,fNVoisinMALong(15)
  ,fMaxStdDevMA(25.0)
  ,fMaxRMSLong(0.8)
  ,fMaxRejFrac(0.15)
  ,fTOFBCMin(-25.0)
  ,fTOFBCMax(50.0)
  ,fUseTOFBC(kFALSE)
  ,fFilterOutliers(kTRUE)
  ,fFatalOnMissingDrift(kTRUE)
  ,fMaxFitYErr2(1.0)
  ,fMaxFitXErr2(9.)
  ,fMaxFitXYCorr(0.95)
  ,fLTMCut(0.75)
  //
  ,fMaxSigY(1.1)
  ,fMaxSigZ(0.7)
  ,fMinValidVoxFracDrift(0.65)
  ,fMaxBadXBinsToCover(4)
  ,fMinGoodXBinsToCover(3)
  ,fMaxBadRowsPerSector(0.4)
  //
  ,fNY2XBins(15)
  ,fNZ2XBins(5)
  ,fNXBins(-1)
  ,fNXYBinsProd(0)
  ,fDZ2X(0)
  ,fDX(0)
  ,fDZ2XI(0)
  ,fDXI(0)
  ,fNGVoxPerSector(0)
  //
  ,fMaxY2X(0)
  ,fDY2X(0)
  ,fDY2XI(0)

  ,fKernelType(kGaussianKernel)

  ,fNTrSelTot(0)
  ,fNTrSelTotWO(0)
  ,fNReadCallTot(0)
  ,fNBytesReadTot(0)
  ,fNTestTracks(0)
  ,fTracksRate(0)
  ,fTOFBCTestH(0)
  ,fHVDTimeInt(0)
  ,fHVDTimeCorr(0)
  ,fVDriftParam(0)
  ,fVDriftGraph(0)
  ,fCorrTime(0)

  ,fStatTree(0)

  ,fDTS()
  ,fDTC()
  ,fDTV()
  ,fDeltaStr()

  ,fTimeStamp(0)
  ,fNCl(0)
  ,fQ2Pt(0)
  ,fTgLam(0)
{
  SetTMinMax(tmin,tmax);
  for (int i=0;i<kResDim;i++) {
    for (int j=0;j<2;j++) fNPCheb[i][j] = -1; //determine from job binning// 15;
    fChebPrecD[i] = 100e-4;
  }

  memset(fLastSmoothingRes,0,kResDim*4*sizeof(double));
  memset(fTestStat,0,kCtrNbr*kCtrNbr*sizeof(float));
  for (int i=kVoxDim;i--;) {
    fUniformBins[i] = kTRUE;
    fKernelScaleEdge[i] = 1.0f;
    fStepKern[i] = 1;
    fKernelWInv[i] = 0; // calculated later
    fSmoothPol2[i] = kFALSE;
  }
  fSmoothPol2[kVoxX] = kTRUE;
  fSmoothPol2[kVoxF] = kTRUE;
  //
  for (int i=kVoxHDim;i--;) fNBProdSt[i] = 0;
  for (int i=kVoxDim;i--;) fNBProdSectG[i] = 0;
  //
  for (int i=0;i<kNSect2;i++) {
    fSectGVoxRes[i] = 0;
    fTmpTree[i] = 0;
    fStatHist[i] = 0;
    fArrNDStat[i] = 0;
    fTmpFile[i] = 0;
    memset(fValidFracXBin[i],0,kNPadRows*sizeof(float));
    fNSmoothingFailedBins[i] = 0;
  }
  SetKernelType();
}

//________________________________________
AliTPCDcalibRes::~AliTPCDcalibRes() 
{
  // d-tor
  delete fChebCorr;
  delete[] fMaxY2X;
  delete[] fDY2X;
  delete[] fDY2XI;
  delete fVDriftParam;
  delete fVDriftGraph;
  for (int i=0;i<kNSect2;i++) {
    delete fSectGVoxRes[i];
    delete fStatHist[i];
  }
  delete fTracksRate;
  delete fTOFBCTestH;
  delete fHVDTimeInt;
  delete fHVDTimeCorr;
  delete fInputChunks;
}

//________________________________________
void AliTPCDcalibRes::SetExternalDetectors(int det)
{
  // set external detectos choice
  if (det<0 || det>=kNExtDetComb) {
    AliErrorF("Invalid external detector %d, allowed range 0:%d, see header file enum{...kNExtDetComb}",det,kNExtDetComb-1);
    return;
  }
  fExtDet = det;
}

//________________________________________
void AliTPCDcalibRes::CalibrateVDrift()
{
  // run VDrift calibration from residual trees. 
  // Parameters  time fDeltaTVD (binning) and  fSigmaTVD (smoothing) can be changed from their
  // default values by their setters
  TStopwatch sw;
  const float kSafeMargin = 0.25f;
  if (!fInitDone) Init();
  Long64_t tmn = fTMinCTP-fTMin>kLargeTimeDiff ? fTMinCTP : fTMin;
  Long64_t tmx = fTMax-fTMaxCTP>kLargeTimeDiff ? fTMaxCTP : fTMax;
  Int_t duration = tmx - tmn;
  int nTBins = duration/fDeltaTVD+1;
  //
  // check available statistics
  if (!EstimateStatistics()) AliFatal("Cannot process further");
  float fracMult = fTestStat[kCtrNtr][kCtrNtr];
  Bool_t useTOFBC = fUseTOFBC;
  float statEst = 0;
  if (useTOFBC) {
    if      (fExtDet==kUseTRDonly) statEst = fTestStat[kCtrBC0][kCtrTRD];
    else if (fExtDet==kUseTOFonly) statEst = fTestStat[kCtrBC0][kCtrTOF];
    else if (fExtDet==kUseITSonly) statEst = fTestStat[kCtrBC0][kCtrITS];
    else if (fExtDet==kUseTRDorTOF) statEst = fTestStat[kCtrBC0][kCtrTOF]; // contribution of TRD is negligable
  }
  else {
    if      (fExtDet==kUseTRDonly) statEst = fTestStat[kCtrTRD][kCtrTRD];
    else if (fExtDet==kUseTOFonly) statEst = fTestStat[kCtrTOF][kCtrTOF];
    else if (fExtDet==kUseITSonly) statEst = fTestStat[kCtrITS][kCtrITS];
    else if (fExtDet==kUseTRDorTOF) statEst = fTestStat[kCtrTRD][kCtrTRD]+
				      fTestStat[kCtrTOF][kCtrTOF]-fTestStat[kCtrTRD][kCtrTOF];
  }
  statEst *= fTestStat[kCtrNtr][kCtrNtr]; // loss due to the mult selection
  statEst *= fNTestTracks*fInputChunks->GetEntriesFast()*(1.-kSafeMargin); // safety margin for losses (outliers etc)
  //
  statEst *= TMath::Min(1.f,float(duration)/(fTMaxCTP-fTMinCTP));
  float statEstTBin = statEst/nTBins;
  // select tracks matching to time window and write compact local trees
  //
  CollectData(kVDriftCalibMode);
  //
  FitDrift(tmn,tmx,nTBins);

  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
  //
}

//________________________________________
void AliTPCDcalibRes::FitDrift(Int_t tmin,Int_t tmax, Int_t ntbins)
{
  // fit v.drift params
  const int kUsePoints0 = 1000000; // nominal number of points for 1st estimate
  const int kUsePoints1 = 3000000; // nominal number of points for time integrated estimate
  const int kNXBins = 150, kNYBins = 150, kNTCorrBins=100;
  const float kMaxTCorr = 0.02; // max time correction for vdrift
  const float kDiscardDriftEdge = 10.; // discard min and max drift values with this margin
  TStopwatch sw;  sw.Start();
  AliSysInfo::AddStamp("FitDrift",0,0,0,0);
  //
  if (!fInitDone) Init();
  //
  TString zresFileName = Form("%s.root",kDriftResFileName);
  TFile* zresFile = TFile::Open(zresFileName.Data());
  if (!zresFile) AliFatalF("file %s not found",zresFileName.Data());
  TString treeName = Form("resdrift");
  TTree *zresTree = (TTree*) zresFile->Get(treeName.Data());
  if (!zresTree) AliFatalF("tree %s is not found in file %s",treeName.Data(),zresFileName.Data());
  //
  dtv_t *dtvP = &fDTV; 
  zresTree->SetBranchAddress("dtv",&dtvP);
  int entries = zresTree->GetEntries();
  if (!entries) AliFatalF("tree %s in %s has no entries",treeName.Data(),zresFileName.Data());
  //
  float prescale = kUsePoints0/float(entries);
  int npUse = entries*prescale;
  //
  float *xArr = new Float_t[npUse];
  float *zArr = new Float_t[npUse];
  int nAcc = 0;
  float dmaxCut = (kZLim[0]+kZLim[1])/2. - kDiscardDriftEdge;
  float res[2],err[3];
  for (int i=0;i<entries;i++) {
    if (gRandom->Rndm()>prescale) continue;
    zresTree->GetEntry(i);
    if (fDTV.drift<kDiscardDriftEdge || fDTV.drift>dmaxCut) continue;
    xArr[nAcc] = fDTV.drift;
    zArr[nAcc] = fDTV.side>0 ? fDTV.dz : -fDTV.dz;
    nAcc++;
  }
  AliInfoF("Will use %d points out of %d for outliers rejection estimate",nAcc,entries);
  //
  float sigMAD = AliTPCDcalibRes::FitPoly1Robust(nAcc,xArr,zArr,res,err,0.95);
  if (sigMAD<0) AliFatal("Unbinned fit failed");
  AliInfoF("Will clean outliers outside of |dz*side-(%.3e+drift*(%.3e))|<3*%.3f band",res[0],res[1],sigMAD);
  delete[] xArr;
  delete[] zArr;
  //
  const float outlCut = sigMAD*3.0f;
  delete fHVDTimeInt;
  fHVDTimeInt = new TProfile2D("driftTInt","",kNXBins,-dmaxCut,dmaxCut,kNYBins,-kMaxX,kMaxX);
  fHVDTimeInt->SetXTitle("side*drift");
  fHVDTimeInt->SetYTitle("ylab*side*drift/zmax");  
  fHVDTimeInt->SetZTitle("#deltaz*side");  
  fHVDTimeInt->SetDirectory(0);
  prescale = kUsePoints1/float(entries);
  npUse = entries*prescale;
  AliInfoF("Will use %d points out of %d for time-integrated estimate",npUse,entries);
  for (int i=0;i<entries;i++) {
    if (gRandom->Rndm()>prescale) continue;
    zresTree->GetEntry(i);
    if (fDTV.drift<kDiscardDriftEdge || fDTV.drift>dmaxCut) continue;
    float dz = fDTV.side>0 ? fDTV.dz : -fDTV.dz;
    if (TMath::Abs(dz - (res[0]+res[1]*fDTV.drift)) > outlCut) continue;
    float sdrift = fDTV.side>0 ? fDTV.drift : -fDTV.drift;
    fHVDTimeInt->Fill(sdrift, fDTV.ylab*sdrift/kZLim[fDTV.side<0], dz);
  }
  TF2* ftd = new TF2("ftd","[0]+[1]*sign(x)+[2]*sign(x)*y + [3]*sign(x)*x",-250,250,-250,250);
  ftd->SetParameters(res[0],0.,0.,res[1]); // initial values from unbinned fit
  AliInfoF("Fitting time-integrated vdrift params by %s",ftd->GetTitle());
  TFitResultPtr rf = fHVDTimeInt->Fit(ftd,"0S");
  int ndf = rf->Ndf();
  float chi2 = rf->Chi2();
  AliInfoF("Fit chi2: %f per %d DOFs -> %f",chi2,ndf,ndf>0 ? chi2/ndf : -1.f);
  delete fVDriftParam;
  fVDriftParam = new TVectorD(4);
  for (int i=4;i--;) (*fVDriftParam)[i] = ftd->GetParameter(i);
  //
  // time correction
  delete fHVDTimeCorr;
  fHVDTimeCorr = new TH2F("driftTCorr","drift time correction",ntbins,tmin,tmax,kNTCorrBins,-kMaxTCorr,kMaxTCorr);
  fHVDTimeCorr->SetXTitle("time");
  fHVDTimeCorr->SetYTitle("time");  
  //
  double *vpar = fVDriftParam->GetMatrixArray();
  for (int i=0;i<entries;i++) {
    zresTree->GetEntry(i);
    if (fDTV.drift<kDiscardDriftEdge || fDTV.drift>dmaxCut) continue;
    float dz = fDTV.side>0 ? fDTV.dz : -fDTV.dz;
    float sdrift = fDTV.side>0 ? fDTV.drift : -fDTV.drift;
    float d2z = fDTV.drift/kZLim[fDTV.side<0];
    Double_t expected = vpar[0]+vpar[1]*fDTV.side + vpar[2]*fDTV.ylab*d2z + vpar[3]*fDTV.drift;
    dz -= expected;
    if (TMath::Abs(dz) > outlCut) continue;
    fHVDTimeCorr->Fill(fDTV.t,  dz/fDTV.drift, d2z); // ?? why this weight 
  }
  //
  sw.Stop(); 
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
  AliSysInfo::AddStamp("FitDrift",1,0,0,0);
}

//________________________________________
void AliTPCDcalibRes::ProcessFromDeltaTrees()
{
  // process from residual trees
  TStopwatch sw;
  // select tracks matching to time window and write compact local trees
  EstimateStatistics();
  //
  CollectData(kDistExtractMode);
  //
  if (fNTrSelTot<fMinTracksToUse) {
    AliErrorF("Low statistics: number of contributing tracks %d, min.requested %d",fNTrSelTot,fMinTracksToUse);
    TString stopOnLosStat = gSystem->Getenv("stopOnLowStat");
    if (!stopOnLosStat.IsNull()) {
      AliInfo("Stop on low statistics requested: abandoning map creation");
      exit(1);
    }
    else {
      AliInfo("No stop on low statistics requested: starting map creation");
    }
  }
  ProcessFromLocalBinnedTrees();
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
  //
}

//________________________________________
void AliTPCDcalibRes::ProcessFromLocalBinnedTrees()
{
  // process starting from local binned trees created by CollectData(kDistExtractMode)
  TStopwatch sw;
  sw.Start();

  // do per-sector projections and fits
  ProcessResiduals();
  //
  //  ProcessDispersions();
  //
  CreateCorrectionObject();
  //
  WriteResTree();
  //
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
}

//________________________________________
void AliTPCDcalibRes::ReProcessFromResVoxTree(const char* resTreeFile, Bool_t backup)
{
  // reprocess starting from the raw data filled from existing resVox tree
  TStopwatch sw;
  sw.Start();
  if (!LoadResTree(resTreeFile)) return;
  ReProcessResiduals();
  //
  if (fChebCorr) delete fChebCorr; fChebCorr = 0;
  //
  CreateCorrectionObject();
  //
  if (backup) { 
    TString inps = resTreeFile;
    if (inps == GetVoxResFileName()) {
      TString inpsb = resTreeFile;
      inpsb.ReplaceAll(".root","_1.root");
      rename(inps.Data(),inpsb.Data());
      AliInfoF("Input file %s backed up to %s",inps.Data(),inpsb.Data());
    }
  }
  WriteResTree();
  //
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
}

//________________________________________
AliTPCDcalibRes* AliTPCDcalibRes::Load(const char* fname)
{
  // load AliTPCDcalibRes object from input file
  TFile* fl = TFile::Open(fname);
  if (!fl) {AliErrorClassF("Failed to open %s",fname); return 0;}
  TList* lstk = fl->GetListOfKeys();
  TIter next(lstk);
  TKey* key = 0;
  AliTPCDcalibRes* res = 0;
  while (key=(TKey*)next()) {
    TString keyt = key->GetTitle();
    if (keyt == "AliTPCDcalibRes") {res = (AliTPCDcalibRes*)fl->Get(key->GetName()); break;}
  }
  if (!res) AliErrorClassF("Did not find AliTPCDcalibRes object in %s",fname);
  return res;
  //
}

//________________________________________
void AliTPCDcalibRes::Save(const char* name)
{
  // save itself
  TString names = name;
  if (names.IsNull()) {
    //    names = Form("%s_run%d_%lld_%lld.root",IsA()->GetName(),fRun,fTMin,fTMax);
    names = Form("%s.root",IsA()->GetName());
    names.ToLower();
  }
  TFile* flout = TFile::Open(names.Data(),"recreate");
  this->Write("",TObject::kOverwrite);
  flout->Close();
  delete flout;
  AliInfoF("Saved itself to %s",names.Data());
  //
}

//==================================================================================
void AliTPCDcalibRes::Init()
{          
  // do the initialization once
  const int kMaxResBins = 0xff;
  AliSysInfo::AddStamp("ProjStart",0,0,0,0);
  if (fInitDone) {AliInfo("Init already done"); return;}
  if (fRun<1) {
    int run = TString(gSystem->Getenv("runNumber")).Atoi();
    if (run<1) AliFatal("Run number is neither set nor provided via runNumber env.var");
    SetRun(run);
  }
  //
  AliCDBManager* man = AliCDBManager::Instance();
  if (fOCDBPath.IsNull()) fOCDBPath = "raw://";
  if (!man->IsDefaultStorageSet()) man->SetDefaultStorage(fOCDBPath);
  if (man->GetRun()!=fRun) man->SetRun(fRun); 
  //
  // memorize GRP time
  AliGRPObject* grp = (AliGRPObject*)man->Get(AliCDBPath("GRP/GRP/Data"))->GetObject();
  fTMinGRP = grp->GetTimeStart();
  fTMaxGRP = grp->GetTimeEnd();
  //
  if (fUseTOFBC) { // check if TOF is present
    Int_t activeDetectors = grp->GetDetectorMask();
    if (!(activeDetectors&AliDAQ::DetectorPattern("TOF"))) {
      AliWarning("Disabling TOF BC validation since TOF is not in the GRP");
      fUseTOFBC = kFALSE;
    }
    if (fTOFBCMin>=fTOFBCMax) {
      AliWarningF("Disabling TOF BC validation: inconsistent cuts %.3f:%.3f",fTOFBCMin,fTOFBCMax);
      fUseTOFBC = kFALSE;      
    }
  }
  // get real data time from CTP scalers
  AliTriggerRunScalers* scalers = (AliTriggerRunScalers*)man->Get(AliCDBPath("GRP/CTP/Scalers"))->GetObject();
  Int_t nEntries = scalers->GetScalersRecords()->GetEntriesFast();
  fTMinCTP = scalers->GetScalersRecord(0         )->GetTimeStamp()->GetSeconds();
  fTMaxCTP = scalers->GetScalersRecord(nEntries-1)->GetTimeStamp()->GetSeconds();
  //
  // init histo for track rate
  fTracksRate = new TH1F("TracksRate","TracksRate", 1+fTMaxCTP-fTMinCTP, -0.5+fTMinCTP,0.5+fTMaxCTP);
  fTracksRate->SetDirectory(0);
  //
  InitGeom();
  SetName(Form("run%d_%lld_%lld",fRun,fTMin,fTMax));
  SetTitle(IsA()->GetName());
  //
  // define boundaries
  InitBinning();
  //
  LoadVDrift(); //!!!
  //
  // prepare aux info for stat and residuals histo bin calculation, see doc of TNDArray bin calculation
  fNBProdSt[kVoxHDim-1] = 1;
  fNBProdSectG[kVoxDim-1] = 1;  
  for (int i=kVoxHDim-1;i--;) {   // +2 to account for under/over-flows
    fNBProdSt[i] = fNBProdSt[i+1]*(2 + ((i==kVoxDim-1) ? kVoxHDim     : fNBins[i+1]));
  }
  for (int i=kVoxDim-1;i--;) {
    fNBProdSectG[i] = fNBProdSectG[i+1]*fNBins[i+1];
  }
  //
  AliSysInfo::AddStamp("Init",0,0,0,0);
  //
  fInitDone = kTRUE;
}

//_____________________________________________________
void AliTPCDcalibRes::CloseDeltaFile(TTree* dtree)
{
  // close input delta chunk
  TFile* fl = 0;
  TDirectory* dir = dtree->GetDirectory();
  if (dir) fl = dir->GetFile();
  delete dtree;
  if (fl) fl->Close();
  delete fl;
}

//_____________________________________________________
TTree* AliTPCDcalibRes::InitDeltaFile(const char* name, Bool_t connect, const char* treeName) 
{
  // init residuals delta file, attach necessary branches
  // 
  static delta_t *delta = &fDeltaStr;
  TString fileNameString(name);
  if (fileNameString.Contains("alien://") && (!gGrid || (gGrid && !gGrid->IsConnected()))) TGrid::Connect("alien://");
  TFile* file = TFile::Open(fileNameString.Data());
  if (!file) {
    AliErrorF("Cannot open file %s",fileNameString.Data());
    return 0;
  }
  TTree* tree = (TTree*)file->Get(treeName);
  if (!tree) {
    AliErrorF("No tree %s in %s",treeName,fileNameString.Data());
    delete file; file = 0;
    return 0;
  }
  //
  if (!connect) return tree;
  //
  Bool_t needTRD = fExtDet==kUseTRDorTOF || fExtDet==kUseTRDonly;
  Bool_t needTOF = fExtDet==kUseTRDorTOF || fExtDet==kUseTOFonly;
  //
  tree->SetCacheLearnEntries(fLearnSize);
  tree->SetCacheSize(0);
  tree->SetCacheSize(fCacheInp*kMByte);
  //
  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("timeStamp",kTRUE);
  tree->SetBranchStatus("itsOK",kTRUE);
  tree->SetBranchStatus("trdOK",kTRUE);
  tree->SetBranchStatus("tofOK",kTRUE);
  tree->SetBranchStatus("vecR.",kTRUE);
  tree->SetBranchStatus("vecSec.",kTRUE);
  tree->SetBranchStatus("vecPhi.",kTRUE);
  tree->SetBranchStatus("vecZ.",kTRUE);
  tree->SetBranchStatus("track.*",kTRUE);      
  tree->SetBranchStatus("npValid",kTRUE);
  tree->SetBranchStatus("its0.",kTRUE);
  tree->SetBranchStatus("its1.",kTRUE);
  tree->SetBranchStatus("tofBC",kTRUE);
  tree->SetBranchStatus("nPrimTracks",kTRUE);
  //
  tree->SetBranchAddress("timeStamp",&fDeltaStr.timeStamp);
  tree->SetBranchAddress("itsOK",&fDeltaStr.itsOK);
  tree->SetBranchAddress("trdOK",&fDeltaStr.trdOK);
  tree->SetBranchAddress("tofOK",&fDeltaStr.tofOK);
  tree->SetBranchAddress("vecR.",&fDeltaStr.vecR);
  tree->SetBranchAddress("vecSec.",&fDeltaStr.vecSec);
  tree->SetBranchAddress("vecPhi.",&fDeltaStr.vecPhi);
  tree->SetBranchAddress("vecZ.",&fDeltaStr.vecZ);
  tree->SetBranchAddress("track.",&fDeltaStr.param);
  tree->SetBranchAddress("npValid",&fDeltaStr.npValid);
  tree->SetBranchAddress("its0.",&fDeltaStr.vecDYITS);
  tree->SetBranchAddress("its1.",&fDeltaStr.vecDZITS);
  tree->SetBranchAddress("tofBC",&fDeltaStr.tofBC);
  tree->SetBranchAddress("nPrimTracks",&fDeltaStr.nPrimTracks);
  //
  if (needTRD) {
    tree->SetBranchStatus("trd0.",kTRUE);
    tree->SetBranchStatus("trd1.",kTRUE);
    tree->SetBranchAddress("trd0.",&fDeltaStr.vecDYTRD);
    tree->SetBranchAddress("trd1.",&fDeltaStr.vecDZTRD);
  }
  if (needTOF) {
    tree->SetBranchStatus("tof0.",kTRUE);
    tree->SetBranchStatus("tof1.",kTRUE);
    tree->SetBranchAddress("tof0.",&fDeltaStr.vecDYTOF);
    tree->SetBranchAddress("tof1.",&fDeltaStr.vecDZTOF);
  }
  //
  tree->GetEntry(0);
  //
  return tree;
}

//_____________________________________________________
Int_t AliTPCDcalibRes::ParseInputList()
{
  // convert text file to array of file names
  if (fInputChunks) delete fInputChunks;
  if (fResidualList.IsNull()) {
    AliError("Input residuals list is not assigned");
    return 0;
  }
  TString  chunkList = gSystem->GetFromPipe(TString::Format("cat %s",fResidualList.Data()).Data());
  fInputChunks = chunkList.Tokenize("\n");  
  return fInputChunks->GetEntriesFast();
}

//_____________________________________________________
Bool_t AliTPCDcalibRes::EstimateStatistics()
{
  // make rough estimate of statistics with different options
  AliInfo("Performing rough statistics check");
  if (!fInitDone) Init();
  if (!AliGeomManager::GetGeometry()) InitGeom(); // in case started from saved object
  // 
  // pick 1st chunk
  int nChunks = (!fInputChunks) ? ParseInputList() : fInputChunks->GetEntriesFast();
  //
  TTree *tree = 0;
  int chunk=0; // read 1st accessible chunk
  while ( !(tree=InitDeltaFile(fInputChunks->At(chunk)->GetName())) && ++chunk<nChunks) {}
  if (!tree) {
    AliError("No data chunk was accessible");
    return kFALSE;
  }
  //
  // make sure these branches are always connected in InitDeltaFile
  const char* kNeedBr[]={"itsOK","trdOK","tofOK","tofBC","nPrimTracks"}; 
  TBranch* br[kCtrNbr];
  for (int i=0;i<kCtrNbr;i++) {
    br[i] = tree->GetBranch(kControlBr[i]);
    if (!br[i]) AliFatalF("Control branch %s is not in the delta tree",kControlBr[i]);
    if (!br[i]->GetAddress()) AliFatalF("Control branch %s address is not set",kControlBr[i]);
    if (!tree->GetBranchStatus(kControlBr[i])) AliFatalF("Control branch %s is not active",kControlBr[i]);
  }	       
  fNTestTracks = tree->GetEntries();
  memset(fTestStat,0,kCtrNbr*kCtrNbr*sizeof(float));
  if (fTOFBCTestH) delete fTOFBCTestH;
  fTOFBCTestH = new TH1F("TOFBCtest",Form("TOF BC for %d tracks",fNTestTracks),1000,-500.,500.);
  fTOFBCTestH->SetDirectory(0);
  Bool_t condOK[kCtrNbr];
  for (int itr=0;itr<fNTestTracks;itr++) {
    for (int ib=kCtrNbr;ib--;) br[ib]->GetEntry(itr);
    condOK[kCtrITS] = fDeltaStr.itsOK;
    condOK[kCtrTRD] = fDeltaStr.trdOK;
    condOK[kCtrTOF] = fDeltaStr.tofOK;
    condOK[kCtrBC0] = fDeltaStr.tofBC>fTOFBCMin && fDeltaStr.tofBC<=fTOFBCMax;
    condOK[kCtrNtr] = fNPrimTracksCut<0 || fDeltaStr.nPrimTracks<=fNPrimTracksCut;
    for (int i=kCtrNbr;i--;) for (int j=kCtrNbr;j--;) if (condOK[i]&&condOK[j]) fTestStat[i][j]++;
    fTOFBCTestH->Fill(fDeltaStr.tofBC);
  }
  if (fNTestTracks)  for (int i=kCtrNbr;i--;)  for (int j=kCtrNbr;j--;) fTestStat[i][j] /= fNTestTracks;
  AliInfoF("Accepted statistics wrt %d tracks in chunk id=%d (out of %d)",fNTestTracks,chunk,nChunks);
  printf(" %11s",""); for (int i=0;i<kCtrNbr;i++) printf("  %11s",kControlBr[i]); printf("\n");
  for (int i=0;i<kCtrNbr;i++) {
    printf("*%11s",kControlBr[i]);
    for (int j=0;j<kCtrNbr;j++) printf("  %11.5f",fTestStat[i][j]);
    printf("\n");
  }
  //
  CloseDeltaFile(tree);
  return kTRUE;
}

//_____________________________________________________
void AliTPCDcalibRes::CollectData(const int mode) 
{
  const float kEps = 1e-6;
  const float q2ptIniTolerance = 1.5;
  if (!fInitDone) Init();
  if (!AliGeomManager::GetGeometry()) InitGeom(); // in case started from saved object
  //
  TStopwatch swTot;
  swTot.Start();
  fNTrSelTot = 0;
  fNTrSelTotWO = 0;
  fNReadCallTot = 0;
  fNBytesReadTot = 0;
  //
  float maxAbsResid = kMaxResid - kEps; // discard residuals exceeding this
  Bool_t correctVDrift = kTRUE;
  if (mode==kVDriftCalibMode) {
    AliInfo("VDrift calibration mode: drift correction disabled");
    maxAbsResid = kMaxResidZVD - kEps; // vdrift calibration may see larger residuals
    correctVDrift = kFALSE;
  }
  else {
    if (!fVDriftGraph || !fVDriftParam) {
      AliErrorF("We are in mode %d but DriftGraph=%p or VDriftParam=%p missing",
		mode,fVDriftGraph,fVDriftParam);
      if (fFatalOnMissingDrift) AliFatal("Aborting...");
    }
  }
  CreateLocalResidualsTrees(mode);
  //
  // if cheb object is present, do on-the-fly init to attach internal structures
  if (fChebCorr) fChebCorr->Init();
  // prepare input tree
  int nChunks = (!fInputChunks) ? ParseInputList() : fInputChunks->GetEntriesFast();
  //
  AliSysInfo::AddStamp("ProjInit",0,0,0,0);
  //
  for (int ichunk=0;ichunk<nChunks;ichunk++) {
    //
    int ntrSelChunkWO=0, ntrSelChunk=0,nReadCallsChunk=0,nBytesReadChunk=0;
    //
    TStopwatch swc;
    swc.Start();
    TString deltaFName = fInputChunks->At(ichunk)->GetName();
    TTree *tree = InitDeltaFile(deltaFName.Data());
    if (!tree) continue;
    //
    TBranch* brTime = tree->GetBranch("timeStamp");
    TBranch* brTRDOK = tree->GetBranch("trdOK");
    TBranch* brTOFOK = tree->GetBranch("tofOK");
    TBranch* brITSOK = tree->GetBranch("itsOK");
    TBranch* brTOFBC = 0;
    if (fUseTOFBC) brTOFBC = tree->GetBranch("tofBC");
    //
    int nTracks = tree->GetEntries();
    AliInfoF("Processing %d tracks of %s",nTracks,deltaFName.Data());

    float residHelixY[kNPadRows],residHelixZ[kNPadRows];
    //
    // reset the cache when swithching between the timeStamp and Event read modes
    Bool_t lastReadMatched = kFALSE; 
    for (int itr=0;itr<nTracks;itr++) {
      nBytesReadChunk += brTime->GetEntry(itr);
      fTimeStamp = fDeltaStr.timeStamp;
      if (fTimeStamp<fTMin  || fTimeStamp>fTMax) {
	if (lastReadMatched && fSwitchCache) { // reset the cache
	  tree->SetCacheSize(0);
	  tree->SetCacheSize(fCacheInp*kMByte);
	  lastReadMatched = kFALSE;
	}
	continue;	
      }
      //
      brITSOK->GetEntry(itr);
      brTRDOK->GetEntry(itr);
      brTOFOK->GetEntry(itr);
      if (!fDeltaStr.itsOK) continue;     
      if (!fDeltaStr.trdOK && fExtDet==kUseTRDonly) continue;
      if (!fDeltaStr.tofOK && fExtDet==kUseTOFonly) continue;
      if (!fDeltaStr.tofOK && !fDeltaStr.trdOK && fExtDet!=kUseITSonly) continue;
      //
      if (brTOFBC && brTOFBC->GetEntry(itr) && (fDeltaStr.tofBC<fTOFBCMin || fDeltaStr.tofBC>fTOFBCMax)) continue;      
      //
      if (!lastReadMatched && fSwitchCache) { // reset the cache before switching to event reading mode
	tree->SetCacheSize(0);
	tree->SetCacheSize(fCacheInp*kMByte);
      }
      lastReadMatched = kTRUE;
      nBytesReadChunk += tree->GetEntry(itr);
      if (fNPrimTracksCut>0 && fDeltaStr.nPrimTracks>fNPrimTracksCut) continue;
      //
      fQ2Pt = fDeltaStr.param->GetParameter()[4];
      fTgLam = fDeltaStr.param->GetParameter()[3];
      if (TMath::Abs(fQ2Pt)>kMaxQ2Pt*q2ptIniTolerance) continue;
      //
      const Float_t *vSec= fDeltaStr.vecSec->GetMatrixArray();
      const Float_t *vPhi= fDeltaStr.vecPhi->GetMatrixArray();
      const Float_t *vR  = fDeltaStr.vecR->GetMatrixArray();
      const Float_t *vZ  = fDeltaStr.vecZ->GetMatrixArray();
      const Float_t *vDYITS = fDeltaStr.vecDYITS->GetMatrixArray();
      const Float_t *vDZITS = fDeltaStr.vecDZITS->GetMatrixArray();
      //
      const Float_t *vDY=0,*vDZ = 0;
      if (fExtDet==kUseTRDonly || (fExtDet==kUseTRDorTOF && fDeltaStr.trdOK)) {
	vDY = fDeltaStr.vecDYTRD->GetMatrixArray();
	vDZ = fDeltaStr.vecDZTRD->GetMatrixArray();
      }
      else if (fExtDet==kUseITSonly) {  // ignore other detectos
	vDY = fDeltaStr.vecDYITS->GetMatrixArray();
	vDZ = fDeltaStr.vecDZITS->GetMatrixArray();
      }
      else { // only TOF
	vDY = fDeltaStr.vecDYTOF->GetMatrixArray();
	vDZ = fDeltaStr.vecDZTOF->GetMatrixArray();	
      }
      //
      fCorrTime = (correctVDrift && fVDriftGraph) ? fVDriftGraph->Eval(fTimeStamp):0; // for VDrift correction
      //
      fNCl = 0;
      // 1st iteration: collect data in cluster frame
      for (int ip=0;ip<fDeltaStr.npValid;ip++) { // 1st fill selected track data to buffer for eventual outlier rejection
	if (vR[ip]<kInvalidR || vDY[ip]<kInvalidRes || vDYITS[ip]<kInvalidRes) continue;
	//
	fArrX[fNCl]   = -1;
	fArrR[fNCl]   = vR[ip];  // X (R) is the same for cluster and track
	fArrZTr[fNCl] = vZ[ip];  // Z of ITS track was stored!!
	fArrDY[fNCl]  = vDY[ip]; // this is also the track coordinate in cluster frame
	fArrDZ[fNCl]  = vDZ[ip];
	fArrPhi[fNCl] = vPhi[ip];
	int rocID = TMath::Nint(vSec[ip]);
	//
	// !!! fArrZTr corresponds to ITS track Z, we need that of TRD-ITS
	fArrZTr[fNCl] += fArrDZ[fNCl] - vDZITS[ip]; // recover ITS-TRD track position from ITS and deltas
	
	if (fFixAlignmentBug && !fDeltaStr.param->TestBit(kAlignmentBugFixedBit)) {
	  FixAlignmentBug(rocID, fQ2Pt, fBz, fArrPhi[fNCl], fArrR[fNCl], fArrZTr[fNCl], fArrDY[fNCl],fArrDZ[fNCl]);
	}
	if (fArrPhi[fNCl]<0) fArrPhi[fNCl] += 2.*TMath::Pi();
	//
	// correct for drift velocity calibration if needed
	if (correctVDrift) fArrDZ[fNCl] += GetDriftCorrection(fArrZTr[fNCl],fArrR[fNCl],fArrPhi[fNCl],rocID);
	//
	fArrSectID[fNCl] = rocID%kNSect2; // 0-36 for sectors from A0 to C17
	//
	fNCl++;
      }
      // fit track coordinates by helix to get interpolated track q/pt: 
      // more precise than the distorted TPC q/pt
      if (fNCl<fMinNCl) continue;
      //
      ntrSelChunkWO++;
      //
      float q2ptTPC = fQ2Pt;
      Bool_t resH = CompareToHelix(residHelixY,residHelixZ);
      //
      if (fFilterOutliers && !resH) continue; // too strong deviation to helix, discard track
      if (TMath::Abs(fQ2Pt)>kMaxQ2Pt) continue; // now we have more precise estimate of q/pt
      //
      // 2nd iteration: convert everything to sector frame
      // *****************************************************************
      //
      // All these manipulations are needed because the ResidualTree is stored
      // in cluster frame, while we need the sector frame
      //
      // *****************************************************************
      int nc0 = fNCl; 
      fNCl = 0;
      for (int ip=0;ip<nc0;ip++) {
	int side = ((fArrSectID[ip] /kNSect)&0x1);
	float sna = TMath::Sin(fArrPhi[ip]-(0.5f +fArrSectID[ip]%kNSect)*kSecDPhi);
	float csa = TMath::Sqrt((1.f-sna)*(1.f+sna));
	//
	// by using propagation in cluster frame in AliTPCcalibAlignInterpolation::Process,
	// the X of the track is evaluated not at the pad-row x=r*csa but at x=r*sca-dy*sna
	double xrow = fArrR[ip]*csa;
	double dx   = fArrDY[ip]*sna;
	double xtr  = xrow - dx;
	double ycl  = fArrR[ip]*sna;      // cluster Y in the sector frame
	double ytr  = ycl + fArrDY[ip]*csa; // track Y in the sector frame at x=xtr is 
	//
	double ztr  = fArrZTr[ip];          // Z of the track at x=xtr
	double zcl  = ztr - fArrDZ[ip];     // and the Z of the cluster is Ztr-deltaZ
	//
	// Now we need to take the track to real pad-row X
	// use linear extrapolation:
	float tgs = fArrTgSlp[ip];
	if (TMath::Abs(tgs)>kMaxTgSlp) continue;
	ytr += dx*tgs;
	double csXtrInv = TMath::Sqrt(1.+tgs*tgs); // (inverse cosine of track angle)
	ztr += dx*fTgLam*csXtrInv;
	//
	// assign to arrays and recalculate residuals
	fArrX[fNCl]   = xrow;
	fArrYTr[fNCl] = ytr;
	fArrZTr[fNCl] = ztr;
	//
	fArrYCl[fNCl] = ycl;
	fArrZCl[fNCl] = zcl;
	fArrDY[fNCl]  = ytr - ycl;
	fArrDZ[fNCl]  = ztr - zcl;
	//
	// we don't want under/overflows
	if (TMath::Abs(fArrDY[fNCl])>maxAbsResid) continue;
	if (TMath::Abs(fArrDZ[fNCl])>maxAbsResid) continue;
	//
	if (fArrX[fNCl]<kMinX || fArrX[fNCl]>kMaxX) continue;
	if (TMath::Abs(fArrZCl[fNCl])>kZLim[side]) continue;;
	//
	// End of manipulations to go to the sector frame
	//
	fNCl++;
      }

      if (fFilterOutliers && !ValidateTrack()) continue;

      ntrSelChunk++;
      
      switch(mode) {
      case kVDriftCalibMode:     FillDriftResidualsTrees(); break;
      case kDistExtractMode:     FillLocalResidualsTrees(); break;
      case kDistClosureTestMode: FillCorrectedResiduals();  break;
      default: AliFatalF("Uknown mode %d",mode);
      };
    } // loop over tracks
    //
    swc.Stop();
    TFile* chunkFile = tree->GetDirectory()?tree->GetDirectory()->GetFile():0;
    nReadCallsChunk =  chunkFile ? chunkFile->GetReadCalls():0;
    AliInfoF("Chunk%3d: selected %d tracks (%d with outliers) from chunk %d | %.1f MB read in %d read calls",
	     ichunk,ntrSelChunk,ntrSelChunkWO, ichunk,float(nBytesReadChunk)/kMByte,nReadCallsChunk); swc.Print();
    fNTrSelTot += ntrSelChunk;
    fNTrSelTotWO += ntrSelChunkWO;
    fNReadCallTot += nReadCallsChunk;
    fNBytesReadTot += nBytesReadChunk;
    //
    CloseDeltaFile(tree);
    AliSysInfo::AddStamp("ProjTreeLoc", ichunk ,fNTrSelTot,fNTrSelTot,fNReadCallTot );
    //
    if (fNTrSelTot > fMaxTracks) {
      AliInfo("Max number of tracks exceeded");
      break;
    }
    //
  } // loop over chunks
  //
  // write/close local trees
  CloseLocalResidualsTrees(mode);
  //
  AliInfoF("Summary: selected %d tracks (%d with outliers) | %.1f MB read in %d read calls",
	   fNTrSelTot,fNTrSelTotWO,float(fNBytesReadTot)/kMByte,fNReadCallTot); 
  swTot.Print();

  AliSysInfo::AddStamp("ProjTreeLocSave");

  if (mode==kDistExtractMode) WriteStatHistos();
  //
}

//________________________________________________
void AliTPCDcalibRes::FillDriftResidualsTrees()
{
  // fill local trees for vdrift calibration
  fDTV.t = fTimeStamp;
  for (int icl=fNCl;icl--;) {
    if (fArrR[icl]<kInvalidR) continue; // rejected outlier
    Bool_t isCside = ((fArrSectID[icl]/kNSect)&0x1);
    fDTV.side  = isCside ? -1:1;
    fDTV.dz    = fArrDZ[icl];
    fDTV.drift = kZLim[isCside] - fDTV.side*fArrZCl[icl];
    fDTV.ylab  = fArrR[icl]*TMath::Sin(fArrPhi[icl]);
    // 
    fTmpTree[0]->Fill();
  }
  if (fTracksRate) fTracksRate->Fill(fTimeStamp); // register track time
  //
}

//________________________________________________
void AliTPCDcalibRes::FillLocalResidualsTrees()
{
  // fill local trees with binned data
  float voxVars[kVoxHDim]={0}; // voxel variables (unbinned)
  for (int icl=fNCl;icl--;) {
    if (fArrX[icl]<kInvalidR) continue; // rejected outlier
    int sectID = fArrSectID[icl]; // 0-35 numbering
    // 
    // calculate voxel variables and bins
    // 
    if (!FindVoxelBin(sectID, fArrX[icl], fArrYCl[icl], fArrZCl[icl], fDTS.bvox, voxVars)) continue;    
    fDTS.dy   = fArrDY[icl];
    fDTS.dz   = fArrDZ[icl];
    fDTS.tgSlp = fArrTgSlp[icl];
    //
    fTmpTree[sectID]->Fill();
    //
    // fill statistics on distribution within the voxel, last dimension, kVoxV is for Nentries
    ULong64_t binToFill = GetBin2Fill(fDTS.bvox,kVoxV); // bin of sector stat histo
    float &binEntries = fArrNDStat[sectID]->At(binToFill); // entries in the voxel
    float oldEntries  = binEntries++;
    float norm        = 1.f/binEntries;
    for (int iv=kVoxDim;iv--;) {
      float &mean = fArrNDStat[sectID]->At(binToFill+iv-kVoxV);
      mean = ( mean*oldEntries + voxVars[iv]) * norm; // account new bin entry in averages calculation
    }
    //
  } // loop over clusters
  //
  if (fTracksRate) fTracksRate->Fill(fTimeStamp); // register track time
  //
}

//________________________________________________
void AliTPCDcalibRes::FillCorrectedResiduals()
{
  // fill local trees result of closure test: corrected distortions
  
  float voxVars[kVoxHDim]={0}; // voxel variables (unbinned)
  fDTC.t = fTimeStamp;
  fDTC.q2pt   = fQ2Pt;
  fDTC.tgLam  = fTgLam;
  //
  for (int icl=fNCl;icl--;) {
    if (fArrX[icl]<kInvalidR) continue; // rejected outlier
    int sectID = fArrSectID[icl]; // 0-35 numbering
    // 
    // extract correction
    // calculate voxel variables and bins
    if (!FindVoxelBin(sectID,fArrX[icl], fArrYCl[icl], fArrZCl[icl], fDTC.bvox, voxVars)) continue;    
    int row159 = GetRowID(fArrX[icl]);
    if (row159<0) continue;
    float corr[3];

    fChebCorr->Eval(sectID, row159, fArrYCl[icl]/fArrX[icl], fArrZCl[icl]/fArrX[icl], corr);
    // 
    fDTC.dyR = fArrDY[icl];
    fDTC.dzR = fArrDZ[icl];

    fDTC.dyC = fArrDY[icl] - (corr[kResY]-corr[kResX]*fArrTgSlp[icl]);
    fDTC.dzC = fArrDZ[icl] - (corr[kResZ]-corr[kResX]*fTgLam); // we evaluate at pad-row

    fDTC.tgSlp  = fArrTgSlp[icl];
    fDTC.x      = fArrX[icl];
    fDTC.y      = fArrYCl[icl];
    fDTC.z      = fArrZCl[icl];
    //
    fTmpTree[sectID]->Fill();
    //
  } // loop over clusters
}

//________________________________________________
void AliTPCDcalibRes::CreateLocalResidualsTrees(int mode)
{
  // temporary trees for local delta's storage
  //
  static dts_t *dtsP = &fDTS;
  static dtc_t *dtcP = &fDTC;
  static dtv_t *dtvP = &fDTV;
  TString namef;
  if (mode==kVDriftCalibMode) {
    namef = Form("%s.root",kDriftResFileName);
    fTmpFile[0] = TFile::Open(namef.Data(),"recreate");
    fTmpTree[0] = new TTree("resdrift","");
    fTmpTree[0]->Branch("dtv", &dtvP);
  }
  else if (mode==kDistExtractMode||mode==kDistClosureTestMode) {    
    for (int is=0;is<kNSect2;is++) {
      if      (mode==kDistExtractMode)         namef = Form("%s%d.root",kLocalResFileName,is);
      else /*if (mode==kDistClosureTestMode)*/ namef = Form("%s%d.root",kClosureTestFileName,is);
      fTmpFile[is] = TFile::Open(namef.Data(),"recreate");
      fTmpTree[is] = new TTree(Form("ts%d",is),"");
      //
      if (mode==kDistExtractMode) {
	fTmpTree[is]->Branch("dts", &dtsP);
	//fTmpTree[is]->SetAutoFlush(150000);
	//
	fStatHist[is] = CreateVoxelStatHisto(is);
	fArrNDStat[is] = (TNDArrayT<float>*)&fStatHist[is]->GetArray();
      }
      else if (mode==kDistClosureTestMode) {
	fTmpTree[is]->Branch("dtc", &dtcP);
      }
    }
  }
  else AliFatalF("Unknown mode %d",mode);
  //
}

//________________________________________________
void AliTPCDcalibRes::CloseLocalResidualsTrees(int /*mode*/)
{
  // close trees for local delta's storage
  //
  for (int is=0;is<kNSect2;is++) {
    if (!fTmpFile[is]) continue;
    fTmpFile[is]->cd();
    fTmpTree[is]->Write("", TObject::kOverwrite);
    delete fTmpTree[is];
    fTmpTree[is] = 0;
    fTmpFile[is]->Close();
    delete fTmpFile[is];
    fTmpFile[is] = 0;
  }
  //
}

//__________________________________________________________________________________
Bool_t AliTPCDcalibRes::CompareToHelix(float *resHelixY, float *resHelixZ)
{
  // compare track to helix, refit q/pt and tgLambda and build array of tg(slope) at pad-rows
  const double kEps = 1e-12;
  float xlab[kNPadRows],ylab[kNPadRows],spath[kNPadRows]; // lab X,Y rotated to for sectort of 1st cluster
  // fill lab coordinates
  float crv = TMath::Abs(fQ2Pt*fBz*0.299792458e-3f), cs,sn;
  int sectPrev=-1,sect0 = fArrSectID[0]%kNSect; // align to the sector of 1st point
  float phiSect = (sect0+0.5)*20*TMath::DegToRad();
  double sna = TMath::Sin(phiSect), csa = TMath::Cos(phiSect);
  //
  spath[0] = 0.f;
  for (int ip=0;ip<fNCl;ip++) {
    cs = TMath::Cos(fArrPhi[ip]-phiSect);
    sn = TMath::Sin(fArrPhi[ip]-phiSect);
    xlab[ip] = fArrR[ip]*cs - fArrDY[ip]*sn;
    ylab[ip] = fArrDY[ip]*cs + fArrR[ip]*sn;
    if (ip) {
      float dx = xlab[ip]-xlab[ip-1];
      float dy = ylab[ip]-ylab[ip-1];
      float ds2 = dx*dx+dy*dy;
      float ds  = TMath::Sqrt(ds2); // circular path
      if (ds*crv>0.05) { 
	// account for the arc-chord difference as 1st 2 terms of asin expansion	
	ds *= (1.f+ds2*crv*crv/24.f);
      }
      spath[ip] = spath[ip-1]+ds;
    }
  }
  double xcSec=0,ycSec=0,xc=0,yc=0,r=0;
  FitCircle(fNCl,xlab,ylab,xcSec,ycSec,r,resHelixY);
  // determine qurvature
  float phi0 = TMath::ATan2(ylab[0],xlab[0]);
  if (phi0<0) phi0 += TMath::Pi()*2;
  float phi1 = TMath::ATan2(ylab[fNCl-1],xlab[fNCl-1]);
  if (phi1<0) phi1 += TMath::Pi()*2;
  float dphi = phi1-phi0;
  int curvSign = 1;
  if (dphi>0) {
    if (dphi<TMath::Pi()) curvSign = -1; // clockwise, no 2pi-0 crossing
  }
  else if (dphi<-TMath::Pi()) curvSign = -1; // clockwise, 2pi-0 crossing
  //
  fQ2Pt = curvSign/(r*fBz*0.299792458e-3f);
  //
  // calculate circle coordinates in the lab frame
  xc = xcSec*csa - ycSec*sna;
  yc = ycSec*csa + xcSec*sna;
  //
  float pol1z[2],pol1zE[4] ;
  Bool_t resfZ = FitPoly1(spath, fArrZTr, 0, fNCl, pol1z, pol1zE);
  //
  fTgLam = pol1z[1]; // new tg. lambda
  // extract deviations wrt helical fit and fill track slopes in sector frame
  float hmnY=1e9,hmxY=-1e9,hmnZ=1e9,hmxZ=-1e9;

  for (int ip=0;ip<fNCl;ip++) {
    float val = fArrZTr[ip] - (pol1z[0]+spath[ip]*pol1z[1]);
    resHelixZ[ip] = val;
    if (val<hmnZ) hmnZ = val;
    if (val>hmxZ) hmxZ = val;
    //    
    val = resHelixY[ip];
    if (val<hmnY) hmnY = val;
    if (val>hmxY) hmxY = val;
    //  
    int sect = fArrSectID[ip]%kNSect;
    if (sect!=sect0) {
      sect0 = sect;
      phiSect = (0.5f + sect)*kSecDPhi;
      sna = TMath::Sin(phiSect);
      csa = TMath::Cos(phiSect);
      xcSec = xc*csa + yc*sna; // recalculate circle center in the sector frame
    }
    // find intersection of the circle with the padrow
    // 1) equation of circle in lab: x=xc+r*cos(tau), y=yc+r*sin(tau)
    // 2) equation of circle in sector frame: 
    //    x=xc'+R*cos(tau-alpSect), y=yc'+R*sin(tau-alpSect)
    //    with xc'=xc*cos(-alp)-yc*sin(-alp); yc'=yc*cos(-alp)+xc*sin(-alp)
    // The circle and padrow at X cross at cos(tau) = (X-xc*csa+yc*sna)/R
    // Hence the derivative of y vs x in sector frame:
    cs = TMath::Cos(fArrPhi[ip]-phiSect);
    double xRow = fArrR[ip]*cs; 
    double cstalp = (xRow - xcSec)/r;
    if (TMath::Abs(cstalp)>1.-kEps) { // track cannot reach this padrow
      cstalp = TMath::Sign(1.-kEps,cstalp);
    }
    // and the slope in sector frame is +-1/tg(acos(cstalp)) = +-cstalp/sqrt(1-cstalp^2)
    // The sign is defined by the fact that in B+ the slope of q- should increase with X.
    // Since the derivative of cstalp/sqrt(1-cstalp^2) on X is positive, just look on qB
    fArrTgSlp[ip] = cstalp/TMath::Sqrt((1.-cstalp)*(1.+cstalp));
    if (fQ2Pt*fBz>0) fArrTgSlp[ip] = -fArrTgSlp[ip];
  }
  //
  //  if (TMath::Abs(hmxY-hmnY)>fMaxDevYHelix || TMath::Abs(hmxZ-hmnZ)>fMaxDevZHelix)
  //    printf("MinMax%d: %e %e %e %e\n",evID,hmnY,hmxY,hmnZ,hmxZ);
  return TMath::Abs(hmxY-hmnY)<fMaxDevYHelix && TMath::Abs(hmxZ-hmnZ)<fMaxDevZHelix;
}

//________________________________________________
void AliTPCDcalibRes::ClosureTest()
{
  // correct distortions
  TStopwatch sw;
  sw.Start();
  if (!fChebCorr) {
    AliError("Chebyshev correction object was not created, cannot run closure test");
    return;
  }
  CollectData(kDistClosureTestMode);
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
  //
}

//________________________________________________
void AliTPCDcalibRes::ProcessResiduals()
{
  // project local trees, extract distortions
  if (!fInitDone) Init(); //{AliError("Init not done"); return;}
  LoadStatHistos();
  AliSysInfo::AddStamp("ProcResid",0,0,0,0);
  //
  for (int is=0;is<kNSect2;is++) ProcessSectorResiduals(is);
  //
  AliSysInfo::AddStamp("ProcResid",1,0,0,0);
  //
}

//________________________________________________
void AliTPCDcalibRes::ProcessDispersions()
{
  // extract distortions of corrected Y residuals ||| DEPRECATED
  if (!fInitDone) Init(); //{AliError("Init not done"); return;}
  //
  LoadStatHistos();
  for (int is=0;is<kNSect2;is++) {
    ProcessSectorDispersions(is);
    if (fDeleteSectorTrees) {
      TString sectFileName = Form("%s%d.root",kLocalResFileName,is);
      AliInfoF("Deleting %s",sectFileName.Data());
      unlink(sectFileName.Data());
    }
  }
  //
}

//________________________________________________
void AliTPCDcalibRes::ProcessSectorDispersions(int is)
{
  // extract dispersion of corrected residuals ||| DEPRECATED
  const float kEps = 1e-6;

  if (!fInitDone) {AliError("Init not done"); return;}
  TStopwatch sw;  sw.Start();
  AliSysInfo::AddStamp("ProcessSectorDispersions",is,0,0,0);

  TString sectFileName = Form("%s%d.root",kLocalResFileName,is);
  TFile* sectFile = TFile::Open(sectFileName.Data());
  if (!sectFile) AliFatalF("file %s not found",sectFileName.Data());
  TString treeName = Form("ts%d",is);
  TTree *sectTree = (TTree*) sectFile->Get(treeName.Data());
  if (!sectTree) AliFatalF("tree %s is not found in file %s",treeName.Data(),sectFileName.Data());
  //
  dts_t *dtsP = &fDTS; 
  sectTree->SetBranchAddress("dts",&dtsP);
  int npoints = sectTree->GetEntries();
  if (!npoints) {
    AliWarningF("No entries for sector %d",is);
    delete sectTree;
    sectFile->Close(); // to reconsider: reuse the file
    delete sectFile;
    return;
  }
  Short_t *resYArr = new Short_t[npoints];
  Short_t *tgslArr = new Short_t[npoints];
  UShort_t *binArr = new UShort_t[npoints];
  Int_t* index = new Int_t[npoints];
  TArrayF dya(1000),tga(1000);//, dza(1000),
  float *dy = dya.GetArray(), *tg = tga.GetArray();//, *dz = dza.GetArray();
  int nacc = 0;
  bres_t* sectData = fSectGVoxRes[is];
  for (int ie=0;ie<npoints;ie++) {
    sectTree->GetEntry(ie);
    if (TMath::Abs(fDTS.tgSlp)>=kMaxTgSlp) continue;
    resYArr[nacc] = Short_t(fDTS.dy*0x7fff/kMaxResid);
    tgslArr[nacc] = Short_t(fDTS.tgSlp*0x7fff/kMaxTgSlp);
    binArr[nacc] = GetVoxGBin(fDTS.bvox);
    nacc++;
  }
  TMath::Sort(nacc, binArr, index, kFALSE); // sort in voxel increasing order
  UShort_t curBin = 0xffff;
  UChar_t bvox[kVoxDim];
  int nproc = 0, npBin = 0;
  while (nproc<nacc) {
    int ip = index[nproc++];
    if (curBin!=binArr[ip]) {
      if (npBin) {
	bres_t& resVox = sectData[curBin];
	GBin2Vox(curBin,resVox.bvox);  // parse voxel
	ProcessVoxelDispersions(npBin,tg,dy,resVox);	
      }
      curBin = binArr[ip];
      npBin = 0;
    }
    if (npBin==dya.GetSize()) {
      dya.Set(100+npBin); dy = dya.GetArray();
      tga.Set(100+npBin); tg = tga.GetArray();
    }
    dy[npBin] = resYArr[ip]*kMaxResid/0x7fff;
    tg[npBin] = tgslArr[ip]*kMaxTgSlp/0x7fff;
    npBin++;
  }
  if (npBin) {
    bres_t& resVox = sectData[curBin];
    GBin2Vox(curBin,resVox.bvox);  // parse voxel
    ProcessVoxelDispersions(npBin,tg,dy,resVox);
  }
  //
  delete[] binArr;
  delete[] resYArr;
  delete[] tgslArr;
  delete[] index;
  //
  delete sectTree;
  sectFile->Close(); // to reconsider: reuse the file
  delete sectFile;
  //
  // now smooth the dispersion
  for (bvox[kVoxX]=0;bvox[kVoxX]<fNXBins;bvox[kVoxX]++) { 
    if (GetXBinIgnored(is,bvox[kVoxX])) continue;
    for (bvox[kVoxZ]=0;bvox[kVoxZ]<fNZ2XBins;bvox[kVoxZ]++) {
      for (bvox[kVoxF]=0;bvox[kVoxF]<fNY2XBins;bvox[kVoxF]++) {
	int binGlo = GetVoxGBin(bvox);
	bres_t *voxRes = &sectData[binGlo];
	Bool_t res = GetSmoothEstimate(is,voxRes->stat[kVoxX],voxRes->stat[kVoxF],voxRes->stat[kVoxZ],
				       BIT(kResD), voxRes->DS);
      }
    }
  }
  //
  sw.Stop(); 
  AliInfoF("Sector%2d | timing: real: %.3f cpu: %.3f",is, sw.RealTime(), sw.CpuTime());
  AliSysInfo::AddStamp("ProcessSectorDispersions",1,0,0,0);
  //
}

//_________________________________________________
void AliTPCDcalibRes::ProcessSectorResiduals(int is)
{
  // process residuals for single sector staring from local binned per-sector trees
  //
  const int kMaxPnt = 30000000; // max points per sector to accept
  TStopwatch sw;  sw.Start();
  AliSysInfo::AddStamp("ProcSectRes",is,0,0,0);
  //
  fNSmoothingFailedBins[is] = 0;
  //
  TString sectFileName = Form("%s%d.root",kLocalResFileName,is);
  TFile* sectFile = TFile::Open(sectFileName.Data());
  if (!sectFile) AliFatalF("file %s not found",sectFileName.Data());
  TString treeName = Form("ts%d",is);
  TTree *sectTree = (TTree*) sectFile->Get(treeName.Data());
  if (!sectTree) AliFatalF("tree %s is not found in file %s",treeName.Data(),sectFileName.Data());
  //
  if (fSectGVoxRes[is]) delete[] fSectGVoxRes[is];
  fSectGVoxRes[is] = new bres_t[fNGVoxPerSector]; // here we keep main result
  bres_t*  sectData = fSectGVoxRes[is];
  // by default set the COG estimates to bin center
  for (int ix=0;ix<fNXBins;ix++) {
    for (int ip=0;ip<fNY2XBins;ip++) {
      for (int iz=0;iz<fNZ2XBins;iz++) {  // extract line in z
	int binGlo = GetVoxGBin(ix,ip,iz);
	bres_t &resVox = sectData[binGlo];
	resVox.bvox[kVoxX] = ix;
	resVox.bvox[kVoxF] = ip;
	resVox.bvox[kVoxZ] = iz;	
	resVox.bsec = is;
	GetVoxelCoordinates(resVox.bsec,resVox.bvox[kVoxX],resVox.bvox[kVoxF],resVox.bvox[kVoxZ],
			    resVox.stat[kVoxX],resVox.stat[kVoxF],resVox.stat[kVoxZ]);
      }
    } 
  }
  //
  dts_t *dtsP = &fDTS; 
  sectTree->SetBranchAddress("dts",&dtsP);
  int npoints = sectTree->GetEntries();
  if (!npoints) {
    AliWarningF("No entries for sector %d, masking all rows",is);
    for (int ix=fNXBins;ix--;) SetXBinIgnored(is,ix);
    delete sectTree;
    sectFile->Close(); // to reconsider: reuse the file
    delete sectFile;
    return;
  }
  if (npoints>kMaxPnt) npoints = kMaxPnt;
  sw.Stop();
  AliInfoF("Sector%2d. Extracted %d points of unbinned data. Timing: real: %.3f cpu: %.3f",
	   is, npoints, sw.RealTime(), sw.CpuTime());
  sw.Start(kFALSE);
  //
  Short_t *resYArr = new Short_t[npoints];
  Short_t *resZArr = new Short_t[npoints];
  Short_t *tgslArr = new Short_t[npoints];
  UShort_t *binArr = new UShort_t[npoints];
  Int_t* index = new Int_t[npoints];
  TArrayF dya(1000),dza(1000),tga(1000);
  float *dy = dya.GetArray(), *dz = dza.GetArray(), *tg = tga.GetArray();
  int nacc = 0;
  for (int ie=0;ie<npoints;ie++) {
    sectTree->GetEntry(ie);
    if (TMath::Abs(fDTS.tgSlp)>=kMaxTgSlp) continue;
    resYArr[nacc] = Short_t(fDTS.dy*0x7fff/kMaxResid);
    resZArr[nacc] = Short_t(fDTS.dz*0x7fff/kMaxResid);
    tgslArr[nacc] = Short_t(fDTS.tgSlp*0x7fff/kMaxTgSlp);
    binArr[nacc] = GetVoxGBin(fDTS.bvox);
    nacc++;
  }
  //
  delete sectTree;
  sectFile->Close(); // to reconsider: reuse the file
  delete sectFile;
  //
  TMath::Sort(nacc, binArr, index, kFALSE); // sort in voxel increasing order
  UShort_t curBin = 0xffff;
  UChar_t bvox[kVoxDim];
  int nproc = 0, npBin = 0;
  while (nproc<nacc) {
    int ip = index[nproc++];
    if (curBin!=binArr[ip]) {
      if (npBin) {
	bres_t& resVox = sectData[curBin];
	GBin2Vox(curBin,resVox.bvox);  // parse voxel
	ProcessVoxelResiduals(npBin,tg,dy,dz,resVox);	
      }
      curBin = binArr[ip];
      npBin = 0;
    }
    if (npBin==dya.GetSize()) {
      dya.Set(100+npBin); dy = dya.GetArray();
      dza.Set(100+npBin); dz = dza.GetArray();
      tga.Set(100+npBin); tg = tga.GetArray();
    }
    dy[npBin] = resYArr[ip]*kMaxResid/0x7fff;
    dz[npBin] = resZArr[ip]*kMaxResid/0x7fff;
    tg[npBin] = tgslArr[ip]*kMaxTgSlp/0x7fff;
    npBin++;
  }
  if (npBin) {
    bres_t &resVox = sectData[curBin];
    GBin2Vox(curBin,resVox.bvox);  // parse voxel
    ProcessVoxelResiduals(npBin,tg,dy,dz,resVox);
  }
  //
  sw.Stop();
  AliInfoF("Sector%2d. Extracted residuals. Timing: real: %.3f cpu: %.3f",
	   is, sw.RealTime(), sw.CpuTime());
  sw.Start(kFALSE);
  
  int nrowOK = ValidateVoxels(is);
  if (!nrowOK) AliWarningF("Sector%2d: all X-bins disabled, abandon smoothing",is);
  else Smooth0(is);
  //
  sw.Stop();
  AliInfoF("Sector%2d. Smoothed residuals. Timing: real: %.3f cpu: %.3f",
	   is, sw.RealTime(), sw.CpuTime());
  sw.Start(kFALSE);

  // now process dispersions
  curBin = 0xffff;

  nproc = 0;
  npBin = 0;
  while (nproc<nacc) {
    int ip = index[nproc++];
    if (curBin!=binArr[ip]) {
      if (npBin) {
	bres_t& resVox = sectData[curBin];
	GBin2Vox(curBin,resVox.bvox);  // parse voxel
	if (!GetXBinIgnored(is,resVox.bvox[kVoxX])) ProcessVoxelDispersions(npBin,tg,dy,resVox);	
      }
      curBin = binArr[ip];
      npBin = 0;
    }
    if (npBin==dya.GetSize()) {
      dya.Set(100+npBin); dy = dya.GetArray();
      tga.Set(100+npBin); tg = tga.GetArray();
    }
    dy[npBin] = resYArr[ip]*kMaxResid/0x7fff;
    tg[npBin] = tgslArr[ip]*kMaxTgSlp/0x7fff;
    npBin++;
  }
  if (npBin) {
    bres_t& resVox = sectData[curBin];
    GBin2Vox(curBin,resVox.bvox);  // parse voxel
    if (!GetXBinIgnored(is,resVox.bvox[kVoxX])) ProcessVoxelDispersions(npBin,tg,dy,resVox);
  }
  //
  // now smooth the dispersion
  for (bvox[kVoxX]=0;bvox[kVoxX]<fNXBins;bvox[kVoxX]++) { 
    if (GetXBinIgnored(is,bvox[kVoxX])) continue;
    for (bvox[kVoxZ]=0;bvox[kVoxZ]<fNZ2XBins;bvox[kVoxZ]++) {
      for (bvox[kVoxF]=0;bvox[kVoxF]<fNY2XBins;bvox[kVoxF]++) {
	int binGlo = GetVoxGBin(bvox);
	bres_t *voxRes = &sectData[binGlo];
	Bool_t res = GetSmoothEstimate(is,voxRes->stat[kVoxX],voxRes->stat[kVoxF],voxRes->stat[kVoxZ],
				       BIT(kResD), voxRes->DS);
      }
    }
  }

  delete[] binArr;
  delete[] resYArr;
  delete[] resZArr;
  delete[] tgslArr;
  delete[] index;
  //
  sw.Stop(); 
  AliInfoF("Sector%2d. Processed dispersion. Timing: real: %.3f cpu: %.3f",is, sw.RealTime(), sw.CpuTime());
  AliSysInfo::AddStamp("ProcSectRes",is,1,0,0);
  //
}
//________________________________________________
void AliTPCDcalibRes::ReProcessResiduals()
{
  // reprocess residuals using raw voxel info filled from existing resVoxTree
  // The raw data is already loaded from the tree
  AliSysInfo::AddStamp("ReProcResid",0,0,0,0);
  for (int is=0;is<kNSect2;is++) ReProcessSectorResiduals(is);
  AliSysInfo::AddStamp("ReProcResid",1,0,0,0);
}

//_________________________________________________
void AliTPCDcalibRes::ReProcessSectorResiduals(int is)
{
  // Reprocess residuals for single sector filled from existing resVoxTree
  // The raw data is already loaded from the tree
  //
  TStopwatch sw;  sw.Start();
  AliSysInfo::AddStamp("RProcSectRes",is,0,0,0);
  //
  fNSmoothingFailedBins[is] = 0;
  bres_t*  sectData = fSectGVoxRes[is];
  if (!sectData) AliFatalF("No SectGVoxRes data for sector %d",is);
  //
  int nrowOK = ValidateVoxels(is);
  if (!nrowOK) AliWarningF("Sector%2d: all X-bins disabled, abandon smoothing",is);
  else Smooth0(is);
  //
  UChar_t bvox[kVoxDim];
  // now smooth the dispersion
  for (bvox[kVoxX]=0;bvox[kVoxX]<fNXBins;bvox[kVoxX]++) { 
    if (GetXBinIgnored(is,bvox[kVoxX])) continue;
    for (bvox[kVoxZ]=0;bvox[kVoxZ]<fNZ2XBins;bvox[kVoxZ]++) {
      for (bvox[kVoxF]=0;bvox[kVoxF]<fNY2XBins;bvox[kVoxF]++) {
	int binGlo = GetVoxGBin(bvox);
	bres_t *voxRes = &sectData[binGlo];
	Bool_t res = GetSmoothEstimate(is,voxRes->stat[kVoxX],voxRes->stat[kVoxF],voxRes->stat[kVoxZ],
				       BIT(kResD), voxRes->DS);
      }
    }
  } 
  sw.Stop(); 
  AliInfoF("Sector%2d. Processed dispersion. Timing: real: %.3f cpu: %.3f",is, sw.RealTime(), sw.CpuTime());
  AliSysInfo::AddStamp("ReProcSectRes",is,1,0,0);
  //
}

//_________________________________________________________
Float_t AliTPCDcalibRes::FitPoly1Robust(int np, float* x, float* y, float* res, float* err, float ltmCut)
{
  // robust pol1 fit, modifies input arrays order
  res[0] = res[1] = 0.f;
  if (np<2) return -1;
  TVectorF yres(7);
  int *indY =  TStatToolkit::LTMUnbinned(np,y,yres,ltmCut);
  if (!indY) return -1;
  // rearrange used events in increasing order
  TStatToolkit::Reorder(np,y,indY);
  TStatToolkit::Reorder(np,x,indY);
  //
  // 1st fit to get crude slope
  int npuse = TMath::Nint(yres[0]);
  int offs =  TMath::Nint(yres[5]);
  // use only entries selected by LTM for the fit
  float a,b;
  AliTPCDcalibRes::medFit(npuse, x+offs, y+offs, a, b, err);
  //
  // don't abuse stack
  float *ycmHeap=0,ycmStack[np<kMaxOnStack ? np:1],*ycm=np<kMaxOnStack ? &ycmStack[0] : (ycmHeap=new float[np]);
  int   *indcmHeap=0,indcmStack[np<kMaxOnStack ? np:1],*indcm=np<kMaxOnStack ? &indcmStack[0] : (indcmHeap=new int[np]);
  //  
  for (int i=np;i--;) ycm[i] = y[i]-(a+b*x[i]);
  TMath::Sort(np,ycm,indcm,kFALSE);
  TStatToolkit::Reorder(np,ycm,indcm);
  TStatToolkit::Reorder(np,y,indcm); // we must keep the same order
  TStatToolkit::Reorder(np,x,indcm);
  //
  // robust estimate of sigma after crude slope correction
  float sigMAD = AliTPCDcalibRes::MAD2Sigma(npuse,ycm+offs);
  // find LTM estimate matching to sigMAD, keaping at least given fraction
  indY = AliTPCDcalibRes::LTMUnbinnedSig(np, ycm, yres, sigMAD,0.5,kTRUE);
  delete[] ycmHeap;
  delete[] indcmHeap;
  //
  if (!indY) return -1;
  // final fit
  npuse = TMath::Nint(yres[0]);
  offs =  TMath::Nint(yres[5]);
  AliTPCDcalibRes::medFit(npuse, x+offs, y+offs, a,b, err);
  res[0] = a;
  res[1] = b;
  return sigMAD;
}

//_________________________________________________
void AliTPCDcalibRes::ProcessVoxelResiduals(int np, float* tg, float *dy, float *dz, bres_t& voxRes)
{
  // extract X,Y,Z distortions of the voxel
  if (np<fMinEntriesVoxel) return;
  TVectorF zres(7);
  voxRes.flags = 0;
  if (!TStatToolkit::LTMUnbinned(np,dz,zres,fLTMCut)) return; 
  //
  float ab[2],err[3];
  float sigMAD = FitPoly1Robust(np,tg,dy,ab,err,fLTMCut);
  if (sigMAD<0) return;
  float corrErr = err[0]*err[2];
  corrErr = corrErr>0 ? err[1]/TMath::Sqrt(corrErr) : -999;
  //printf("N:%3d A:%+e B:%+e / %+e %+e %+e | %+e %+e / %+e %+e\n",np,a,b,err[0],err[1],err[2], zres[1],zres[2], zres[3],zres[4]);
  //
  voxRes.D[kResX] = -ab[1];
  voxRes.D[kResY] = ab[0];
  voxRes.D[kResZ] = zres[1];
  voxRes.E[kResX] = TMath::Sqrt(err[2]);
  voxRes.E[kResY] = TMath::Sqrt(err[0]);
  voxRes.E[kResZ] = zres[4];
  voxRes.EXYCorr  = corrErr;
  voxRes.D[kResD] = voxRes.dYSigMAD = sigMAD; // later will be overriden by real dispersion
  voxRes.dZSigLTM = zres[2];
  //
  // store the statistics
  ULong64_t binStat = GetBin2Fill(voxRes.bvox,kVoxV);
  voxRes.stat[kVoxV] = fArrNDStat[voxRes.bsec]->At(binStat);
  for (int iv=kVoxDim;iv--;) voxRes.stat[iv] = fArrNDStat[voxRes.bsec]->At(binStat+iv-kVoxV);
  //
  voxRes.flags |= kDistDone;
}

//_________________________________________________
void AliTPCDcalibRes::ProcessVoxelDispersions(int np, const float* tg, float *dy, bres_t& voxRes)
{
  // extract Y (Z ignored at the moment) dispersions of the voxel
  // correct Y distortions
  if (np<2) return;
  for (int i=np;i--;) dy[i] -= voxRes.DS[kResY] - voxRes.DS[kResX]*tg[i];
  voxRes.D[kResD] = MAD2Sigma(np,dy);
  voxRes.E[kResD] = voxRes.D[kResD]/TMath::Sqrt(2.*np); // a la gaussian RMS error, this is very crude
  voxRes.flags |= kDispDone;
  //
}

//_____________________________________________________________________
Double_t AliTPCDcalibRes::GetLogL(TH1F* histo, int bin0, int bin1, double &mu, double &sig, double &logL0)
{
  // Calculate log likelihood of normal distribution for the histo between boundaries 
  // bin0 and bin for given mu and sigma assumption. Exact Poisson statistics is assumed
  // Also the approximate "reference" log-likelihood logLO is calculated in the following way:
  // if the Poisson prob. for given bin is "m", then the logL0 gets contribution for this bin
  // log("reference probability"), which I define as a geometric mean of probabilities for
  // observing [M] and [M]+1 entries if m is large: 
  // P_ref = exp(-m)/m! M^m sqrt( m/(M+1) )
  // -> ln(P_ref) = -(1/2+M)*ln(M/m) + M-m + 1/2 ln(M+1) - 1/2 ln(2 pi) -> -1/2 ln(2 pi m) for m>>1 (take m>5)
  //                (precise up to 1/2 ln(2 pi m) term of Stirling formula
  // or           = -m + m*log(m) - log( Gamma(1.+m) )                                     for m<~1
  // 
  // integral
  const double kNuLarge = 5.0, kMinSig2BinH = 0.01;
  double dxh = 0.5*histo->GetBinWidth(1);
  if ((sig/dxh)<kMinSig2BinH) {
    AliWarningClassF("Too small sigma %.4e is provided for bin width %.4e",sig,dxh);
    logL0 = -1;
    return -1e9;
  }
  double sum=0, sum1=0, sum2=0;
  
  for (int ib=bin0;ib<=bin1;ib++) {
    double w = histo->GetBinContent(ib);
    double x = histo->GetBinCenter(ib);
    sum += w;
    sum1 += w*x;
    sum2 += w*x*x;
  }  
  //
  double xb0 = histo->GetBinCenter(bin0)-dxh;
  double xb1 = histo->GetBinCenter(bin1)+dxh;
  //
  if (sum<1e-6) {logL0 = -1e6; return -1e9;}
  mu = sum1/sum;
  sig = sum2/sum - mu*mu;
  sig = TMath::Max(sig>0 ? TMath::Sqrt(sig) : 0.f, dxh/TMath::Sqrt(3)); // don't allow too small sigma
  
  //printf("Sample mu : %e sig: %e in %e %e\n",mu,sig,xb0,xb1);

  // estimated sig, mu are from the truncated sample, try to recover the truth
  GetTruncNormMuSig(xb0,xb1, mu, sig);
  //
  xb0 -= mu;
  xb1 -= mu;
  double sqri2 = 1./(TMath::Sqrt(2.)*sig);
  // normalization constant
  double norm = 2.*sum / (TMath::Erf(xb1*sqri2) - TMath::Erf(xb0*sqri2));
  //
  //  printf("Norm: %e\n",norm);
  // likelihood
  double logL = 0;
  logL0 = 0;
  const double kMinExp = 1e-100;
  for (int i=bin0;i<=bin1;i++) {
    double x = histo->GetBinCenter(i)-mu;
    double w = histo->GetBinContent(i);
    xb0 = x-dxh;
    xb1 = x+dxh;
    // bin expectation: normal integral within the bin
    double nu = 0.5*norm*(TMath::Erf(xb1*sqri2) - TMath::Erf(xb0*sqri2));  
    if (nu<kMinExp) nu = kMinExp;
    double logNFac = w<100 ? TMath::Log(TMath::Factorial(w)) : w*TMath::Log(w)-w + TMath::Log( sqrt(2*TMath::Pi()*w));
    double logNu = TMath::Log(nu);
    double logc = -nu + w*logNu - logNFac;  // contribution of this bin to log-likelihood
    logL += logc;
    // now get the reference contribution
    double logc0 = 0;
    if (nu>kNuLarge) logc0 = -0.5*TMath::Log(2.*TMath::Pi()*nu);
    else {
      logc0 = -nu + nu*logNu - TMath::Log( TMath::Gamma(1.+nu) );
    }
    logL0 += logc0;  // reference LL update
    //printf("b: %d x:%+.2e nstd:%+.2e Exp:%e Obs:%e logc: %e logc0: %e\n",i,x,(x-mu)/sig, nu,w,logc, logc0);

  }
  //  printf("LogL: %e LogL0: %e\n",logL,logL0);
  //
  return logL;
}

//___________________________________________________________________________
void AliTPCDcalibRes::TruncNormMod(double a, double b, double mu0, double sig0, double &muCf, double &sigCf)
{
  // calculate truncated mean and sigma of normal distribution as 
  // mu_tr  = mu0 + sig0*muCf
  // sig_tr = sig0 * sigCf
  //
  const double sqrt2PiI = 1./TMath::Sqrt(TMath::Pi()*2.);
  double sigI = 1./(sig0*TMath::Sqrt(2.));
  double ra = (a-mu0)*sigI, rb = (b-mu0)*sigI;
  double ra2 = ra*ra, rb2 = rb*rb;
  double af = ra2<100 ? sqrt2PiI*TMath::Exp(-ra2) : 0;
  double bf = rb2<100 ? sqrt2PiI*TMath::Exp(-rb2) : 0;
  //  double aF = 0.5*(1.+TMath::Erf(ra)), bF = 0.5*(1.+TMath::Erf(rb)), deltaF = bF-aF
  double deltaF = 0.5*( TMath::Erf(rb) - TMath::Erf(ra) );
  double deltaf = af - bf;
  muCf = deltaf / deltaF;
  sigCf = 1./TMath::Sqrt(1. + TMath::Sqrt(2)*(ra*af-rb*bf)/deltaF - muCf*muCf); 
  //
}

//_____________________________________________________
Bool_t AliTPCDcalibRes::GetTruncNormMuSig(double a, double b, double &mean, double &sig)
{
  // get estimate of real mu and sigma of normal distribution provided
  // the mean and rms of sample truncated between a and b
  const double kMinWindow=1e-2,kEpsRMS = 1e-4, kEpsMu = 1e-4;
  const int kMaxIter = 200;
  //
  if (sig<1e-12) {
    AliWarningClassF("Input sigma %e is too small",sig);
    return kFALSE;
  }
  if ( (b-a)/sig<kMinWindow ) {
    AliWarningClassF("Truncation window %e-%e is %e sigma only",a,b,(b-a)/sig);
    return kFALSE;
  }
  //
  double sig0=sig, mean0=mean; // initial values
  // for protection, don't allow the sigma to grow above a factor of the flat distribution
  double sigMax = 1.2*(b-a)/TMath::Sqrt(12.);
  //
  double m = mean, s = sig;
  for (int i=0;i<kMaxIter;i++) {
    double sclRMS,sclMU;
    TruncNormMod(a,b,m,s, sclMU,sclRMS);
    //
    s = sig * sclRMS;
    double mPrev = m, sPrev = s;
    m = mean - sclMU*s;
    //printf("%d -> M: %e S: %e\n",i,m, s);
    if ( s>sigMax) {

      //      printf("Iteration took sigma to twice of the flat distribution for "
      //	     "mu0=%+.3e sig0=%.3e in %+.3e:%+.3e interval\n",mean0,sig0, a,b);
      if (TMath::Abs(m-mean0)>sig0) {
	//	printf("Abandoning and returning input sigma and mean\n");
	m = mean0; s = sig0;
      }
      break;
    }
    if (TMath::Abs(1.-sPrev/s)<kEpsRMS && TMath::Abs(m-mPrev)<kEpsMu ) break;
  }
  //
  mean = m;
  sig  = s;

  return kTRUE;
}

//_________________________________________________
void AliTPCDcalibRes::InitGeom()
{
  // init geometry and field
  // this requires the field and the geometry ...
  if (fRun<1) AliFatal("Run number is not provided");
  Bool_t geomOK = AliGeomManager::GetGeometry() != 0;
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!geomOK || !fld) { // need to setup ocdb?
    AliCDBManager* man = AliCDBManager::Instance();
    if (fOCDBPath.IsNull()) fOCDBPath = "raw://";
    if (!man->IsDefaultStorageSet()) man->SetDefaultStorage(fOCDBPath);
    if (man->GetRun()!=fRun) man->SetRun(fRun);
  }
  if (!geomOK) {
    AliGeomManager::LoadGeometry();
    AliGeomManager::ApplyAlignObjsFromCDB("TPC");
  }
  if (!fld) {
    AliGRPManager grpMan;
    grpMan.ReadGRPEntry();
    grpMan.SetMagField();
    fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  }
  fBz = fld->SolenoidField();
}


//___________________________________________________________________________
THnF* AliTPCDcalibRes::CreateVoxelStatHisto(int sect)
{
  // prepare histogram to store the means of distributions within the voxel

  // create binning for voxels statistics histograms
  Int_t    voxNBins[kVoxHDim];
  Double_t voxBinMin[kVoxHDim],voxBinMax[kVoxHDim];
  TString  voxAxisName[kVoxHDim];

  voxAxisName[kVoxF] = "Y2X_Bin";
  voxNBins[kVoxF]    = fNY2XBins;
  voxBinMin[kVoxF]   = 0;
  voxBinMax[kVoxF]   = fNY2XBins;
  //
  voxAxisName[kVoxX]   = "X_Bin";
  voxNBins[kVoxX]      = fNXBins;
  voxBinMin[kVoxX]     = 0;
  voxBinMax[kVoxX]     = fNXBins;
  //
  voxAxisName[kVoxZ]   = "Z2X_Bin";
  voxNBins[kVoxZ]      = fNZ2XBins;
  voxBinMin[kVoxZ]     = 0;
  voxBinMax[kVoxZ]     = fNZ2XBins;
  //
  voxAxisName[kVoxV] = "Stat_Bin";
  voxNBins[kVoxV]    = kVoxHDim;
  voxBinMin[kVoxV]   = 0;
  voxBinMax[kVoxV]   = kVoxHDim;
  //
  THnF* h = new THnF(Form("hs%d",sect),"",kVoxHDim,voxNBins,voxBinMin,voxBinMax);
  for (int i=0;i<kVoxHDim;i++) h->GetAxis(i)->SetName(voxAxisName[i].Data());
  h->SetEntries(1); // otherwise drawing does not work well
  return h;
}

//===============================================================
//
//                   TRACK VALIDATION >>>>>>>>>>>>>>>>>>>>>>>>>>>
//

//__________________________________________________________________________________
Bool_t AliTPCDcalibRes::ValidateTrack()
{
  // if (nCl<fMinNCl) return kFALSE;
 if (fNCl<fNVoisinMALong) return kFALSE;

  Bool_t rejCl[kNPadRows];
  float rmsLong = 0.f;
  int nRej = CheckResiduals(rejCl, rmsLong);
  if (float(nRej)/fNCl > fMaxRejFrac) return kFALSE;
  if (rmsLong>fMaxRMSLong) return kFALSE;
  //
  // flag outliers
  for (int i=fNCl;i--;) if (rejCl[i]) fArrR[i] = fArrX[i] = -1;

  return kTRUE;
}

//______________________________________________
void AliTPCDcalibRes::FixAlignmentBug(int sect, float q2pt, float bz, float& alp, 
				      float& x, float &z, float &deltaY, float &deltaZ)
{
  // fix alignment bug: https://alice.its.cern.ch/jira/browse/ATO-339?focusedCommentId=170850&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-170850
  //
  // NOTE: deltaZ in the buggy code is calculated as Ztrack_with_bug - Zcluster_w/o_bug
  static TGeoHMatrix *mCache[72] = {0};
  if (sect<0||sect>=72) {
    AliErrorF("Invalid sector %d",sect);
    return;
  }
  int lr = sect/36 ? (AliGeomManager::kTPC2) : (AliGeomManager::kTPC1);
  TGeoHMatrix* mgt = mCache[sect];
  if (!mgt) {
    int volID = AliGeomManager::LayerToVolUIDSafe(lr,sect%36);
    mgt = new TGeoHMatrix(*AliGeomManager::GetTracking2LocalMatrix(volID));
    mgt->MultiplyLeft(AliGeomManager::GetMatrix(volID));
    mCache[sect] = mgt;
    AliInfoF("Caching matrix for sector %d",sect);
  }  
  double alpSect = ((sect%18)+0.5)*20.*TMath::DegToRad();

  // cluster in its proper alpha frame with alignment bug
  double xyzClUse[3] = {x,0,z}; // this is what we read from the residual tree
  double xyzTrUse[3] = {x, deltaY, z}; // track in bad cluster frame
  //
  // recover cluster Z position by adding deltaZ
  double zClSave = xyzClUse[2] -= deltaZ;  // here the cluster is not affected by Z alignment component of the bug!
  static AliExternalTrackParam trDummy;
  trDummy.Local2GlobalPosition(xyzClUse,alp); // misaligned cluster in global frame
  double xyz0[3]={xyzClUse[0],xyzClUse[1],xyzClUse[2]};
  mgt->MasterToLocal(xyz0,xyzClUse);
  // we got ideal cluster in the sector tracking frame, but now the Z is wrong, since it was not affected by the bug!!!
  //
  xyzClUse[2] = zClSave;

  // go to ideal cluster frame
  trDummy.Local2GlobalPosition(xyzClUse,alpSect); // ideal global
  double alpFix = TMath::ATan2(xyzClUse[1],xyzClUse[0]);    // fixed cluster phi
  trDummy.Global2LocalPosition(xyzClUse,alpFix);     // fixed cluster in in its frame
  //
  trDummy.Local2GlobalPosition(xyzTrUse,alp); // track in global frame
  trDummy.Global2LocalPosition(xyzTrUse,alpFix); // track in cluster frame
  alp = alpFix;
  //
  double dx = xyzTrUse[0] - xyzClUse[0]; // x might not be the same after alignment fix
  // deduce track slopes assuming it comes from the vertex
  double tgphi = tgpXY(xyzClUse[0],xyzTrUse[1]/xyzClUse[0],q2pt,bz);
  xyzTrUse[1] -= dx*tgphi;
  xyzTrUse[2] -= dx*xyzClUse[2]/xyzClUse[0]; // z2x
  //
  x = xyzClUse[0];
  z = xyzTrUse[2]; // we still use track Z as a reference ...
  deltaY = xyzTrUse[1]-xyzClUse[1];
  deltaZ = xyzTrUse[2]-xyzClUse[2];
  //
}


//_______________________________________________________________
int AliTPCDcalibRes::CheckResiduals(Bool_t* mask, float &rmsLongMA)
{

  int ip0=0,ip1;
  int sec0 = fArrSectID[ip0];
  int npLast = fNCl-1;
  //
  const int nMinAcc = 30;
  float yDiffLL[kNPadRows] = {0.f};
  float zDiffLL[kNPadRows] = {0.f};
  float absDevY[kNPadRows] = {0.f};
  float absDevZ[kNPadRows] = {0.f};
  
  rmsLongMA = 0.f;

  memset(mask,0,fNCl*sizeof(Bool_t));
  for (int i=0;i<fNCl;i++) {
    if (fArrSectID[i]==sec0 && i<npLast) continue;
    //
    // sector change or end of input reached
    // run estimators for the points in the same sector
    int npSec = i-ip0;
    if (i==npLast) npSec++;
    //
    DiffToLocLine(npSec, fArrX+ip0, fArrDY+ip0, fNVoisinMA, yDiffLL+ip0);
    DiffToLocLine(npSec, fArrX+ip0, fArrDZ+ip0, fNVoisinMA, zDiffLL+ip0);
    //    DiffToMA(npSec, fArrX+ip0, fArrDY+ip0, fNVoisinMA, yDiffLL+ip0);
    //    DiffToMA(npSec, fArrX+ip0, fArrDZ+ip0, fNVoisinMA, zDiffLL+ip0);
    //
    ip0 = i;
    sec0 = fArrSectID[ip0];
  }
  // store abs deviations
  int naccY=0,naccZ=0;
  for (int i=fNCl;i--;) {
    if (yDiffLL[i]) absDevY[naccY++] = TMath::Abs(yDiffLL[i]);
    if (zDiffLL[i]) absDevZ[naccZ++] = TMath::Abs(zDiffLL[i]);
  }
  //
  // estimate rms on 90% smallest deviations
  int kmnY = 0.9*naccY,kmnZ = 0.9*naccZ;
  if (naccY<nMinAcc || naccZ<nMinAcc) { // mask all
    for (int i=fNCl;i--;) mask[i] = kTRUE;
    return fNCl;
  }

  SelKthMin(kmnY, naccY, absDevY);
  SelKthMin(kmnZ, naccZ, absDevZ);
  float rmsKY=0,rmsKZ=0;
  for (int i=kmnY;i--;) rmsKY += absDevY[i]*absDevY[i];
  for (int i=kmnZ;i--;) rmsKZ += absDevZ[i]*absDevZ[i];
  rmsKY = TMath::Sqrt(rmsKY/kmnY);
  rmsKZ = TMath::Sqrt(rmsKZ/kmnZ);
  //
  if (rmsKY<1e-6 || rmsKZ<1e-6) {
    AliWarningF("Too small RMS: %f %f",rmsKY,rmsKZ);
    for (int i=fNCl;i--;) mask[i] = kTRUE;
    return fNCl;
  }
  //
  //  printf("RMSY %d min of %d: %f | RMSZ %d min of %d: %f\n",kmnY,naccY,rmsKY, kmnZ,naccZ,rmsKZ);
  //
  //
  float rmsKYI = 1./rmsKY;
  float rmsKZI = 1./rmsKZ;
  int nMask=0, nacc = 0;
  float yacc[kNPadRows],yDiffLong[kNPadRows];
  for (int ip=0;ip<fNCl;ip++) {

    yDiffLL[ip] *= rmsKYI;
    zDiffLL[ip] *= rmsKZI;
    float dy = yDiffLL[ip], dz = zDiffLL[ip];
    if (dy*dy+dz*dz>fMaxStdDevMA) {
      mask[ip] = kTRUE;
      nMask++;
    }
    else yacc[nacc++] = fArrDY[ip];
  }
  // rms to long-range moving average wrt surviving clusters
  if (nacc>fNVoisinMALong) {
    DiffToMA(nacc, yacc, fNVoisinMALong, yDiffLong);
    float av=0,rms=0;
    for (int i=0;i<nacc;i++) {
      av += yDiffLong[i];
      rms  += yDiffLong[i]*yDiffLong[i];
    }
    av /= nacc;
    rmsLongMA = rms/nacc - av*av;
    rmsLongMA = rmsLongMA>0 ? TMath::Sqrt(rmsLongMA) : 0.f;
  }
  return nMask;
  //
}

//____________________________________________________________________
void AliTPCDcalibRes::FitCircle(int np, const float* x, const float* y, 
				double &xc, double &yc, double &r, float* dy)
{
  // fit points to circle, if dy!=0, fill residuals
  double x0=0.,y0=0.;
  for (int i=np;i--;) {
    x0 += x[i];
    y0 += y[i];
  } 
  x0 /= np;
  y0 /= np;
  double su2=0,sv2=0,suv=0,su3=0,sv3=0,su2v=0,suv2=0;
  for (int i=np;i--;) {
    double ui = x[i]-x0, ui2 = ui*ui;
    double vi = y[i]-y0, vi2 = vi*vi;
    suv += ui*vi;
    su2 += ui2;
    sv2 += vi2;
    su3 += ui2*ui;
    sv3 += vi2*vi;
    su2v += ui2*vi;
    suv2 += ui*vi2;
  } 
  double rhsU = 0.5*(su3+suv2), rhsV = 0.5*(sv3+su2v);
  double det = su2*sv2-suv*suv;
  double uc  = (rhsU*sv2 - rhsV*suv)/det;
  double vc  = (su2*rhsV - suv*rhsU)/det;
  double r2  = uc*uc + vc*vc + (su2+sv2)/np;
  xc = uc + x0;
  yc = vc + y0;
  //
  if (dy) {
    for (int i=np;i--;) {
      double dx = x[i]-xc;
      double dxr = r2 - dx*dx;
      double ys = dxr>0 ? TMath::Sqrt(dxr) : 0;
      double dy0 = y[i]-yc;
      double dysp = dy0-ys;
      double dysm = dy0+ys;
      dy[i] = TMath::Abs(dysp)<TMath::Abs(dysm) ? dysp : dysm;
    }
  }
  r = TMath::Sqrt(r2);
}

//_________________________________________
void AliTPCDcalibRes::DiffToMA(int np, const float *y, const int winLR, float* diffMA)
{
  // difference to moving average, excluding central element
  //
  double arrSumO[kNPadRows+1], *arrSum=arrSumO+1;
  arrSum[-1] = 0.;
  for (int ip=0;ip<np;ip++) arrSum[ip] = arrSum[ip-1]+y[ip];
  for (int ip=0;ip<np;ip++) {
    diffMA[ip] = 0;
    int ipmn = ip-winLR;
    int ipmx = ip+winLR;
    if (ipmn<0)   ipmn=0;
    if (ipmx>=np) ipmx=np-1;
    int nrm = (ipmx-ipmn);
    if (nrm<winLR) continue;
    double ma = (arrSum[ipmx]-arrSum[ipmn-1] - (arrSum[ip]-arrSum[ip-1]))/nrm;
    diffMA[ip] = y[ip] - ma;
  }

}


//_______________________________________________
int AliTPCDcalibRes::DiffToLocLine(int np, const float* x, const float *y, const int nVoisin, float *diffY)
{
  // calculate the difference between the point and linear extrapolation from neigbourhood
  const float kEps = 1e-9;
  double sumX1b[kNPadRows+1],sumX2b[kNPadRows+1],sumYb[kNPadRows+1],sumYXb[kNPadRows+1];
  double *sumX1 = sumX1b+1, *sumX2 = sumX2b+1, *sumY0 = sumYb+1, *sumYX = sumYXb+1;
  //
  sumX1[-1]=0.f;
  sumX2[-1]=0.f;
  sumY0[-1]=0.f;
  sumYX[-1]=0.f;
  // accumulate matrix elements for whole array, element -1 is at 0
  double sX1=0., sX2=0., sY0=0., sYX=0.;

  for (int ip=0;ip<np;ip++) {
    sumX1[ip] = sX1 += x[ip];
    sumX2[ip] = sX2 += x[ip]*x[ip];
    sumY0[ip] = sY0 += y[ip];
    sumYX[ip] = sYX += y[ip]*x[ip];
    diffY[ip]  = 0.0f;
  }
  // 
  int nAcc = 0;
  for (int ip=0;ip<np;ip++) {
    float &yEst = diffY[ip];
    // estimate from the left
    int ip0 = ip-nVoisin;
    int ip1 = ip+nVoisin;
    if (ip0<0)   ip0=0;
    if (ip1>=np) ip1=np-1;
    int nrm = (ip1-ip0);
    if (nrm<nVoisin) continue;
    int ip0m = ip0-1, ipm = ip-1;
    // extract sum from ip0 to ip1 from cumulant, S00=nrm, excluding tested point
    sX1 = sumX1[ip1] - sumX1[ip0m] - (sumX1[ip]-sumX1[ipm]); // S01
    sX2 = sumX2[ip1] - sumX2[ip0m] - (sumX2[ip]-sumX2[ipm]); // S11
    sY0 = sumY0[ip1] - sumY0[ip0m] - (sumY0[ip]-sumY0[ipm]); // RHS0
    sYX = sumYX[ip1] - sumYX[ip0m] - (sumYX[ip]-sumYX[ipm]); // RHS1
    double det = nrm*sX2 - sX1*sX1;
    if (det<kEps) continue;
    double detI = 1./det;
    // yLEst = offs + slop*x[ip] 
    // with offs=(sY0*sX2-sYX*sX1)/det and slop=(nVoisin*sYX-sY0*sX1)/det
    // inverse err^2 = 1/(errOffs+x^2*errSlop+2*x*errSlpOffs) with
    // errOffs = S11/det, errSlp=S00/det, errSlpOffs=-S01/det
    yEst = y[ip]-((sY0*sX2-sYX*sX1) + (nrm*sYX-sY0*sX1)*x[ip])*detI;
    nAcc++;
    //    
  }
  return nAcc;
  //
}

//====================================================================
int AliTPCDcalibRes::DiffToMedLine(int np, const float* x, const float *y, const int nVoisin, float *diffY)
{
  int nAcc = 0;
  float offs=0,slp=0;
  float buff[kNPadRows];
  for (int ip=0;ip<np;ip++) {
    float &yEst = diffY[ip];
    yEst = 0;
    // estimate from the left
    int ip0 = ip-nVoisin;
    int ip1 = ip+nVoisin;
    if (ip0<0)   ip0=0;
    if (ip1>=np) ip1=np-1;
    int nrm = (ip1-ip0+1);
    if (nrm<nVoisin) continue;
    int ip0m = ip0-1, ipm = ip-1;
    const float *arrx = x+ip0, *arry = y+ip0;
    medFit(nrm, arrx, arry, offs,slp);
    /*
    float asum = 0;
    for (int i=nrm;i--;) {
      buff[i] = arry[i] - (offs + slp*arrx[i]);
      if (i+ip0!=ip) asum += TMath::Abs(buff[i]);
    }
    asum /= nrm-1;
    yEst = buff[ip-ip0]/asum;
    */
    yEst = y[ip] - (offs + slp*x[ip]);
    nAcc++;
  }
  return nAcc;
}

//_________________________________________________________
Int_t* AliTPCDcalibRes::LTMUnbinnedSig(int np, const float *arr, TVectorF &params , Float_t sigTgt, Float_t minFrac, Bool_t sorted)
{
  //
  // LTM : Trimmed keeping at most minFrac of unbinned array to reach targer sigma
  // 
  // Robust statistic to estimate properties of the distribution
  // To handle binning error special treatment
  // for definition of unbinned data see:
  //     http://en.wikipedia.org/w/index.php?title=Trimmed_estimator&oldid=582847999
  //
  // Function parameters:
  //     np      - number of points in the array
  //     arr     - data array (unsorted)
  //     params  - vector with parameters
  //             - 0 - area
  //             - 1 - mean
  //             - 2 - rms 
  //             - 3 - error estimate of mean
  //             - 4 - error estimate of RMS
  //             - 5 - first accepted element (of sorted array)
  //             - 6 - last accepted  element (of sorted array)
  //
  // On success returns index of sorted events 
  //
  static int *index = 0, book = 0;
  static double* w = 0;
  params[0] = 0.0f;
  if (book<np) {
    delete[] index;
    book = np;
    index = new int[book];
    delete[] w;
    w = new double[book+book];
  }
  //
  double *wx1 = w, *wx2 = wx1+np;
  if (!sorted) TMath::Sort(np,arr,index,kFALSE); // sort in increasing order
  else for (int i=0;i<np;i++) index[i]=i;
  // build cumulants
  double sum1=0.0,sum2=0.0;
  for (int i=0;i<np;i++) {
    double x = arr[index[i]];
    wx1[i] = (sum1+=x);
    wx2[i] = (sum2+=x*x);
  }
  //
  int keepMax = np;
  int keepMin = minFrac*np;
  if (keepMin>keepMax) keepMin = keepMax;
  //
  float sig2Tgt = sigTgt*sigTgt;
  while (1) {
    double minRMS = sum2+1e6;
    int keepN = (keepMax+keepMin)>>1;
    if (keepN<2) return 0;
    //
    params[0] = keepN;
    int limI = np - keepN+1;
    for (int i=0;i<limI;i++) {
      int limJ = i+keepN-1;
      Double_t sum1 = wx1[limJ] - (i ? wx1[i-1] : 0.0);
      Double_t sum2 = wx2[limJ] - (i ? wx2[i-1] : 0.0);
      double mean = sum1/keepN;
      double rms2 = sum2/keepN - mean*mean;
      if (rms2>minRMS) continue;
      minRMS = rms2;
      params[1] = mean;
      params[2] = rms2;
      params[5] = i;
      params[6] = limJ;
    }
    if (minRMS<sig2Tgt) keepMin = keepN;
    else                keepMax = keepN;
    if (keepMin>=keepMax-1) break;
  }
  //
  if (!params[0]) return 0;
  params[2] = TMath::Sqrt(params[2]);
  params[3] = params[2]/TMath::Sqrt(params[0]); // error on mean
  params[4] = params[3]/TMath::Sqrt(2.0); // error on RMS
  return index;
}

//___________________________________________________________________
float AliTPCDcalibRes::MAD2Sigma(int np, float* y)
{
  // Sigma calculated from median absolute deviations, https://en.wikipedia.org/wiki/Median_absolute_deviation
  // the input array is not modified
  if (np<2) return 0;
  int nph = np>>1;
  if (nph&0x1) nph -= 1;
  // don't abuse stack
  float *ycHeap=0, ycStack[np<kMaxOnStack ? np:1],*yc=np<kMaxOnStack ? &ycStack[0] : (ycHeap = new float[np]);
  memcpy(yc,y,np*sizeof(float));
  float median = (np&0x1) ? SelKthMin(nph,np,yc) : 0.5f*(SelKthMin(nph-1,np,yc)+SelKthMin(nph,np,yc));
  // build abs differences to median
  for (int i=np;i--;) yc[i] = TMath::Abs(yc[i]-median);
  // now get median of abs deviations
  median = (np&0x1) ? SelKthMin(nph,np,yc) : 0.5f*(SelKthMin(nph-1,np,yc)+SelKthMin(nph,np,yc));
  delete[] ycHeap; // if any...
  return median*1.4826; // convert to Gaussian sigma
}

//___________________________________________________________________
void AliTPCDcalibRes::medFit(int np, const float* x, const float* y, float &a, float &b, float* err,float delI)
{
  // Median linear fit: minimizes abs residuals instead of squared ones
  // Adapted from "Numerical Recipes in C"
  float aa,bb,b1,b2,f,f1,f2,sigb,chisq=0.0f;
  if (np<2) {
    a = b = 0.0;
    if (err) {
      err[0] = err[1] = err[2] = 999.;
    }
    return;
  }
  if (!delI) {
    float sx=0.0f,sxx=0.0f,sy=0.0f,sxy=0.0f,del;
    //
    for (int j=np;j--;) { sx += x[j]; sxx += x[j]*x[j];}
    del = np*sxx-sx*sx;
    //
    for (int j=np;j--;) { sy += y[j]; sxy += x[j]*y[j];}
    //
    delI = 1./del;
    aa = (sxx*sy-sx*sxy)*delI;
    bb = (np*sxy-sx*sy)*delI;
    if (err) {
      err[0] = sxx*delI;
      err[1] = sx*delI;
      err[2] = np*delI;
    }
  }
  else { // initial values provided
    aa = a;
    bb = b;
  }
  //
  for (int j=np;j--;) {
    float temp = y[j]-(aa+bb*x[j]);
    chisq += temp*temp;
  }
  //
  sigb = TMath::Sqrt(chisq*delI);
  b1=bb;
  f1 = RoFunc(np,x,y,b1,aa);
  if (sigb>0) {
    b2 = bb+TMath::Sign(float(3.0f*sigb),f1);
    f2 = RoFunc(np,x,y,b2,aa);
    if (f1==f2) {
      a = aa;
      b = bb;
      return;
    }
    while (f1*f2 > 0.0f) { // bracketing
      bb = b2 + 1.6f*(b2-b1);
      b1 = b2;
      f1 = f2;
      b2 = bb;
      f2 = RoFunc(np,x,y,b2,aa);
    }
    sigb = 0.01*sigb;
    while (fabs(b2-b1)>sigb) {
      bb = b1 + 0.5f*(b2-b1);
      if (bb==b1 || bb==b2) break;
      f = RoFunc(np,x,y,bb,aa);
      if (f*f1 >= 0.0f) {
	f1=f;
	b1=bb;
      } 
      else {
	f2 = f;
	b2 = bb;
      }
    }
  }
  a = aa;
  b = bb;
  //
}

//___________________________________________________________________
float AliTPCDcalibRes::RoFunc(int np, const float* x, const float* y, float b, float &aa)
{
  const float kEPS = 1.0e-7f; 
  static float* arrTmp = 0;
  static int nBook = 0;
  if (np>nBook) { // make sure the buffer is ok
    nBook = np;
    delete[] arrTmp;
    arrTmp = new float[nBook];
  }
  float d,sum=0.0f;
  for (int j=np;j--;) arrTmp[j] = y[j]-b*x[j];
  //
  int nph = np>>1;
  if (np<20) {  // it is faster to do insertion sort 
    for (int i=1;i<np;i++) {
      float v = arrTmp[i];
      int j;
      for (j=i;j--;) if (arrTmp[j]>v) arrTmp[j+1]=arrTmp[j]; else break;
      arrTmp[j+1] = v;
    }
    aa = (np&0x1) ? arrTmp[nph] : 0.5f*(arrTmp[nph-1]+arrTmp[nph]);
  }
  else {
    aa = (np&0x1) ? SelKthMin(nph,np,arrTmp) : 
      0.5f*(SelKthMin(nph-1,np,arrTmp)+SelKthMin(nph,np,arrTmp));
  }
  for (int j=np;j--;) {
    d = y[j] - (b*x[j] + aa);
    if (y[j] != 0.0f) d /= TMath::Abs(y[j]);
    if (TMath::Abs(d) > kEPS) sum += (d >= 0.0f ? x[j] : -x[j]);
  }
  return sum;
}

//___________________________________________________________________
Float_t AliTPCDcalibRes::SelKthMin(int k, int np, float* arr)
{
  // Returns the k th smallest value in the array. The input array will be rearranged
  // to have this value in location arr[k] , with all smaller elements moved before it
  // (in arbitrary order) and all larger elements after (also in arbitrary order).
  // From Numerical Recipes in C++

  int i,ir,j,l,mid;
  float a;
  l=0;ir=np-1;
  for (;;) {
    if (ir<=l+1) {
      if (ir==l+1 && arr[ir]<arr[l]) swap(arr[l],arr[ir]);
      return arr[k];
    } 
    else {
      int mid = (l+ir)>>1, i=l+1;
      swap(arr[mid],arr[i]);
      if (arr[i]>arr[ir]) swap(arr[i],arr[ir]);
      if (arr[l]>arr[ir]) swap(arr[l]  ,arr[ir]);
      if (arr[i]>arr[l])  swap(arr[i],arr[l]);
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i) break;
	swap(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      if (j>=k) ir=j-1;
      if (j<=k) l=i;
    }
  }
}

//===============================================================


//_________________________________________
void AliTPCDcalibRes::LoadVDrift()
{
  // load vdrift params
  fVDriftGraph = 0;
  fVDriftParam = 0;
  TString vdname = Form("%s.root",kDriftFileName);
  TFile *fdrift = 0;
  if (gSystem->AccessPathName(vdname.Data()) || !(fdrift=TFile::Open(vdname.Data())) ) {
    AliWarningF("vdrift file %s not accessible",vdname.Data());
    return;
  }
  TTree * tree = (TTree*)fdrift->Get("fitTimeStat");
  if (!tree) AliWarningF("tree fitTimeStat not avaliable in file %s.root",kDriftFileName);
  else {      
    tree->SetBranchAddress("grTRDReg.",&fVDriftGraph);
    tree->SetBranchAddress("paramRobust.",&fVDriftParam);
    tree->GetEntry(0);
    if (fVDriftGraph==NULL || fVDriftGraph->GetN()<=0) {
      AliWarning("ITS/TRD drift calibration not availalble. Trying ITS/TOF");
      tree->SetBranchAddress("grTOFReg.",&fVDriftGraph);
      tree->GetEntry(0);
    }
  }
  delete tree;
  delete fdrift;
}

//_________________________________________
float AliTPCDcalibRes::GetDriftCorrection(float z, float x, float phi, int rocID)
{
  // apply vdrift correction
  if (!fVDriftParam) return 0.f;
  int side = ((rocID/kNSect)&0x1) ? -1:1; // C:A
  float drift = side>0 ? kZLim[0]-z : z+kZLim[1];
  float gy    = TMath::Sin(phi)*x;
  Double_t pvecFit[3];
  pvecFit[0]= side;             // z shift (cm)
  pvecFit[1]= drift*gy/kZLim[side<0];   // global y gradient
  pvecFit[2]= drift;            // drift length
  float expected = (fVDriftParam==NULL) ? 0:
    (*fVDriftParam)[0]+
    (*fVDriftParam)[1]*pvecFit[0]+
    (*fVDriftParam)[2]*pvecFit[1]+
    (*fVDriftParam)[3]*pvecFit[2];
  return -side*(expected+fCorrTime*drift);
  //
}

//_____________________________________________________
float AliTPCDcalibRes::tgpXY(float x, float y, float q2p, float bz)
{
  // get the tg of primary track inclination wrt padrow given
  // that it was registered at X,Y sector coordinates
  float c = q2p*bz*(-0.299792458e-3f);
  if (TMath::Abs(c)<1e-9) return y/x;
  float r2 = x*x+y*y;
  float det = 4./r2 - c*c;
  float snp  = 0;
  if (det<0) {
    snp = TMath::Sign(-0.8f,c);
    AliWarningF("track of q2p=%f cannot reach x:%f y:%f",q2p,x,y);
  }
  else {
    snp = 0.5f*(y*TMath::Sqrt(det)-c*x); // snp at vertex
    snp += x*c;  // snp at x,y
  }
  return snp/TMath::Sqrt((1.f-snp)*(1.f+snp));
}


//=========================================================================
//
//   IO related methods
//
//=========================================================================
//
//___________________________________________________________________
void AliTPCDcalibRes::WriteStatHistos()
{
  // write stat histos
  TString statOutName = Form("%s.root",kStatOut);
  TFile* statOutFile = TFile::Open(statOutName.Data(),"recreate");
  for (int is=0;is<kNSect2;is++) fStatHist[is]->Write("", TObject::kOverwrite);
  statOutFile->Close();
  delete statOutFile;
  //
}

//___________________________________________________________________
void AliTPCDcalibRes::LoadStatHistos()
{
  // load bin stat histos
  if (fStatHist[0]) return; // histos are in memory
  TString statOutName = Form("%s.root",kStatOut);
  TFile* statOutFile = TFile::Open(statOutName.Data());
  if (!statOutFile) AliFatalF("LoadStatHistos: failed to read file %s",statOutName.Data()); 
  for (int is=0;is<kNSect2;is++) {
    fStatHist[is] = (THnF*) statOutFile->Get(Form("hs%d",is));
    if (!fStatHist[is]) AliFatalF("LoadStatHistos: failed to read secto %d histo from %s",is,statOutName.Data());
    fArrNDStat[is] = (TNDArrayT<float>*)&fStatHist[is]->GetArray();
  }
  statOutFile->Close();
  delete statOutFile;
  //
}

//___________________________________________________________________
void AliTPCDcalibRes::WriteResTree()
{
  // output file for results tree
  TStopwatch sw;
  sw.Start();
  bres_t voxRes, *voxResP=&voxRes;

  AliSysInfo::AddStamp("ResTree",0,0,0,0);
  if (fChebCorr) fChebCorr->Init();
  TFile* flOut = new TFile(GetVoxResFileName(),"recreate");
  TTree* resTree = new TTree("voxRes","final distortions, see GetListOfAliases");
  resTree->Branch("res", &voxRes);
  for (int i=0;i<kVoxDim;i++) {
    resTree->SetAlias(kVoxName[i],Form("bvox[%d]",i));
    resTree->SetAlias(Form("%sAV",kVoxName[i]),Form("stat[%d]",i));
  }
  resTree->SetAlias(Form("%sAV",kVoxName[kVoxV]),Form("stat[%d]",kVoxV));
  for (int i=0;i<kResDim;i++) {
    resTree->SetAlias(kResName[i],Form("D[%d]",i));
    resTree->SetAlias(Form("%sS",kResName[i]),Form("DS[%d]",i));
    resTree->SetAlias(Form("%sC",kResName[i]),Form("DC[%d]",i));
    resTree->SetAlias(Form("%sE",kResName[i]),Form("E[%d]",i));
  }
  //
  resTree->SetAlias("fitOK",Form("(flags&0x%x)==0x%x",kDistDone,kDistDone));
  resTree->SetAlias("dispOK",Form("(flags&0x%x)==0x%x",kDispDone,kDispDone));
  resTree->SetAlias("smtOK",Form("(flags&0x%x)==0x%x",kSmoothDone,kSmoothDone));
  resTree->SetAlias("masked",Form("(flags&0x%x)==0x%x",kMasked,kMasked));
  //
  for (int is=0;is<kNSect2;is++) { 

    bres_t* sectData = fSectGVoxRes[is];
    if (!sectData) {
      AliWarningF("No processed data for sector %d",is);
      continue;
    }

    // now store the data for the sector
    for (int iz=0;iz<fNZ2XBins;iz++) {
      for (int ix=0;ix<fNXBins;ix++) { 
	for (int ip=0;ip<fNY2XBins;ip++) {
	  int binGlo = GetVoxGBin(ix,ip,iz);
	  bres_t *voxel = &sectData[binGlo];
	  if (fChebCorr) {
	    int row = GetRowID(voxel->stat[kVoxX]);
	    int roc = is;
	    if (row>=kNRowIROC) {roc += kNSect2; row -= kNRowIROC;}
	    fChebCorr->Eval(roc,row, voxel->stat[kVoxF],voxel->stat[kVoxZ],voxel->DC);
	  }
	  memcpy(&voxRes,voxel,sizeof(bres_t)); // store in the sector data array
	  resTree->Fill();
	}
      }
    }    
    //
  } // end of sector loop
  //
  flOut->cd();
  resTree->Write("", TObject::kOverwrite);
  delete resTree;
  //
  if (fChebCorr) {
    fChebCorr->Write();
  }
  //
  flOut->Close();
  delete flOut;
  //
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
  AliSysInfo::AddStamp("ResTree",1,0,0,0);

}

//___________________________________________________________________
Bool_t AliTPCDcalibRes::LoadResTree(const char* resTreeFile)
{
  // Fill voxels info from existing resVox tree for reprocessing
  TStopwatch sw;
  sw.Start();
  bres_t voxRes, *voxResP=&voxRes;
  //
  TFile* flIn = new TFile(resTreeFile);
  if (!flIn) {AliErrorF("Failed to open %s",resTreeFile); return kFALSE;}
  TTree* resTree = (TTree*) flIn->Get("voxRes");
  if (!resTree) {AliErrorF("Failed to extract resTree from %s",resTreeFile); delete flIn; return kFALSE;}
  resTree->SetBranchAddress("res", &voxResP);
  //
  int nent = resTree->GetEntries();
  if (nent != kNSect2*fNGVoxPerSector) AliFatalF("resTree from %s has %d voxels per sector, this object: %d",
						 resTreeFile,nent/kNSect2,fNGVoxPerSector);
  //
  for (int is=0;is<kNSect2;is++) { 
    if (fSectGVoxRes[is]) delete[] fSectGVoxRes[is];
    bres_t* sectData = fSectGVoxRes[is] = new bres_t[fNGVoxPerSector];
  }
  for (int ient=0;ient<nent;ient++) {
    //
    resTree->GetEntry(ient);
    bres_t* sectData = fSectGVoxRes[voxRes.bsec];
    int binGlo = GetVoxGBin(voxRes.bvox);
    bres_t *voxel = &sectData[binGlo];
    memcpy(voxel,&voxRes,sizeof(bres_t));
    // the X distortion contribution was already subtracted from the Y,Z components, restore it
    voxel->D[kResZ]  -= voxel->stat[kVoxZ]*voxel->DS[kResX];
    // reset smoothed params
    voxel->flags &= ~(kSmoothDone|kMasked);
    for (int ir=kResDim;ir--;) voxel->DS[ir]=voxel->DC[ir]=0.f;
  } // end of sector loop
  //
  delete resTree;
  flIn->Close();
  delete flIn;
  //
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
  //
  return kTRUE;
}


//=========================================================================
//
//   fitting/smoothing related methods
//
//=========================================================================

//_____________________________________________________
Bool_t AliTPCDcalibRes::FitPoly2(const float* x,const float* y, const float* w, int np, float *res, float *err)
{
  // poly2 fitter
  if (np<3) return kFALSE; // no enough points
  double sumW[5]={0},sumY[3]={0};
  for (int ip=np;ip--;) {
    double ww = w ? w[ip] : 1.0;
    sumW[0] += ww;
    sumY[0] += ww*y[ip];
    sumW[1] += (ww*=x[ip]); 
    sumY[1] += ww*y[ip];
    sumW[2] += (ww*=x[ip]); 
    sumY[2] += ww*y[ip];
    sumW[3] += (ww*=x[ip]); 
    sumW[4] += (ww*=x[ip]); 
  }
  double min00 = (sumW[2]*sumW[4]-sumW[3]*sumW[3]);
  double min01 = (sumW[1]*sumW[4]-sumW[3]*sumW[2]);
  double min02 = (sumW[1]*sumW[3]-sumW[2]*sumW[2]);
  double min11 = (sumW[0]*sumW[4]-sumW[2]*sumW[2]);
  double min12 = (sumW[0]*sumW[3]-sumW[2]*sumW[1]);
  double min22 = (sumW[0]*sumW[2]-sumW[1]*sumW[1]);
  double det  = sumW[0]*min00
    -           sumW[1]*min01
    +           sumW[2]*min02;
  if (TMath::Abs(det)<1e-12) return kFALSE;
  double detI = 1./det;
  double det0 = sumY[0]*min00
    -           sumW[1]*(sumY[1]*sumW[4]-sumW[3]*sumY[2]) 
    +           sumW[2]*(sumY[1]*sumW[3]-sumW[2]*sumY[2]);
  double det1 = sumW[0]*(sumY[1]*sumW[4]-sumW[3]*sumY[2]) 
    -           sumY[0]*min01
    +           sumW[2]*(sumW[1]*sumY[2]-sumY[1]*sumW[2]);
  double det2 = sumW[0]*(sumW[2]*sumY[2]-sumY[1]*sumW[3]) 
    -           sumW[1]*(sumW[1]*sumY[2]-sumY[1]*sumW[2]) 
    +           sumY[0]*min02;
  res[0] = det0*detI;
  res[1] = det1*detI;
  res[2] = det2*detI;
  //
  err[0] = min00*detI; // e00
  err[1] =-min01*detI; // e10
  err[2] = min11*detI; // e11
  err[3] = min02*detI; // e20
  err[4] =-min12*detI; // e21
  err[5] = min22*detI; // e21
  //
  return kTRUE;
}

//_____________________________________________________
Bool_t AliTPCDcalibRes::FitPoly1(const float* x,const float* y, const float* w, int np, float *res, float *err)
{
  // poly1 fitter
  if (np<2) return kFALSE; // no enough points
  double sumW[3]={0},sumY[2]={0};
  for (int ip=np;ip--;) {
    double ww = w ? w[ip]:1.0;
    sumW[0] += ww;
    sumY[0] += ww*y[ip];
    sumW[1] += (ww*=x[ip]); 
    sumY[1] += ww*y[ip];
    sumW[2] += (ww*=x[ip]);
  }
  double det  = sumW[0]*sumW[2] - sumW[1]*sumW[1];
  if (TMath::Abs(det)<1e-12) return kFALSE;
  double detI = 1./det;
  double det0 = sumY[0]*sumW[2] - sumW[1]*sumY[1];
  double det1 = sumW[0]*sumY[1] - sumY[0]*sumW[1];
  res[0] = det0*detI;
  res[1] = det1*detI;
  //
  err[0] = sumW[2]*detI; // e00
  err[1] =-sumW[1]*detI; // e10
  err[2] = sumW[0]*detI; // e11
  //
  return kTRUE;
}

//________________________________
Int_t AliTPCDcalibRes::ValidateVoxels(int isect)
{
  // apply voxel validation cuts, calculate number of low stat or masked voxels
  int cntMasked=0, cntInvalid=0;
  //
  fXBinIgnore[isect].ResetAllBits();
  bres_t* sectData = fSectGVoxRes[isect];
  for (int ix=0;ix<fNXBins;ix++) {
    int nvalidXBin = 0;
    for (int ip=0;ip<fNY2XBins;ip++) {
      for (int iz=0;iz<fNZ2XBins;iz++) {  // extract line in z
	int binGlo = GetVoxGBin(ix,ip,iz);
	bres_t *voxRes = &sectData[binGlo];
	Bool_t voxOK = voxRes->flags&kDistDone && !(voxRes->flags&kMasked);
	if (voxOK) {
	  // check fit errors
	  if (voxRes->E[kResY]*voxRes->E[kResY]>fMaxFitYErr2 ||
	      voxRes->E[kResX]*voxRes->E[kResX]>fMaxFitXErr2 ||
	      TMath::Abs(voxRes->EXYCorr)>fMaxFitXYCorr) voxOK = kFALSE;
	  // check raw distributions sigmas
	  if (voxRes->dYSigMAD > fMaxSigY || voxRes->dZSigLTM > fMaxSigZ) voxOK = kFALSE;
	  if (!voxOK) cntMasked++;
	}
	if (voxOK) {
	  nvalidXBin++;
	  voxRes->flags |= kDistDone;
	}
	else {
	  cntInvalid++;
	  voxRes->flags |= kMasked;
	}
      } // loop over Z
    } // loop over Y
    //
    fValidFracXBin[isect][ix] = Float_t(nvalidXBin)/(fNY2XBins*fNZ2XBins);
    if (fValidFracXBin[isect][ix]<fMinValidVoxFracDrift) {
      AliWarningF("Sector%2d: Xbin%3d has %4.1f%% of voxels valid (%d out of %d)",
		  isect,ix,100*fValidFracXBin[isect][ix],nvalidXBin,fNY2XBins*fNZ2XBins);
    }
    //
  } // loop over X
  //
  // mask X-bins which cannot be smoothed
  int nValidInPatch = 0;
  // 1st loop: find bad regions
  short nbadReg=0,badStart[kNPadRows],badEnd[kNPadRows];
  Bool_t prevBad = kFALSE;
  float fracBadRows = 0;
  for (int ix=0;ix<fNXBins;ix++) {
    if (fValidFracXBin[isect][ix]<fMinValidVoxFracDrift) {
      fracBadRows++;
      if (prevBad) badEnd[nbadReg] = ix;
      else {
	badStart[nbadReg] = badEnd[nbadReg] = ix;
	prevBad = kTRUE; 
      }
    }
    else {
      if (prevBad) {
	prevBad = kFALSE;
	nbadReg++;
      }
    }
  }
  if (prevBad) nbadReg++; // last bad region was not closed
  //
  fracBadRows /= fNXBins;
  if (fracBadRows>fMaxBadRowsPerSector) {
    AliWarningF("Sector%2d: Fraction of bad X-bins %.3f > %.3f: masking whole sector",
		isect,fracBadRows,fMaxBadRowsPerSector);
    for (int ix=0;ix<fNXBins;ix++) SetXBinIgnored(isect,ix);
  }
  else {
  //
  // 2nd loop: disable those regions which cannot be smoothed
    for (int ibad=0;ibad<nbadReg;ibad++) {
      short badSize=badEnd[ibad]-badStart[ibad]+1,badSizeNext=ibad<(nbadReg-1) ? badEnd[ibad]-badStart[ibad]+1 : 0;
      // disable too large bad patches
      if (badSize>fMaxBadXBinsToCover) for (int i=0;i<badSize;i++) SetXBinIgnored(isect,badStart[ibad]+i);
      // disable too small isolated good patches
      if (badSizeNext>fMaxBadXBinsToCover && (badStart[ibad+1]-badEnd[ibad]-1)<fMinGoodXBinsToCover) {
	for (int i=badEnd[ibad]+1;i<badStart[ibad+1];i++) SetXBinIgnored(isect,i);
      }
    }
    if (nbadReg) {
      int ib=0;
      // is 1st good patch too small?
      if (GetXBinIgnored(isect,badStart[0]) && badStart[ib]<fMinGoodXBinsToCover) {
	for (int i=0;i<badStart[ib];i++) SetXBinIgnored(isect,i);
      }
      // last good patch is too small?
      ib = nbadReg-1;
      if (GetXBinIgnored(isect,badStart[ib]) && (fNXBins-badEnd[ib]-1)<fMinGoodXBinsToCover) {
	for (int i=badEnd[ib]+1;i<fNXBins;i++) SetXBinIgnored(isect,i);
      }
    }
    //
  }
  //
  int nMaskedRows = fXBinIgnore[isect].CountBits();
  AliInfoF("Sector%2d: Voxel stat: Masked: %5d(%07.3f%%) Invalid: %5d(%07.3f%%)  -> Masked %3d rows out of %3d",isect,
	   cntMasked, 100*float(cntMasked)/fNGVoxPerSector,
	   cntInvalid,100*float(cntInvalid)/fNGVoxPerSector,
	   nMaskedRows,fNXBins);
  //
  return fNXBins-nMaskedRows;
}

//________________________________
Int_t AliTPCDcalibRes::Smooth0(int isect)
{
  // apply linear regression kernel smoother 
  int cnt = 0;
  bres_t* sectData = fSectGVoxRes[isect];
  for (int ix=0;ix<fNXBins;ix++) {
    if (GetXBinIgnored(isect,ix)) continue;
    for (int ip=0;ip<fNY2XBins;ip++) {
      for (int iz=0;iz<fNZ2XBins;iz++) {  // extract line in z
	int binGlo = GetVoxGBin(ix,ip,iz);
	bres_t *vox = &sectData[binGlo];
	vox->flags &= ~kSmoothDone;
	Bool_t res = GetSmoothEstimate(vox->bsec,vox->stat[kVoxX],vox->stat[kVoxF],vox->stat[kVoxZ],
				       BIT(kResX)|BIT(kResY)|BIT(kResZ), // at this moment we cannot smooth dispersion
				       vox->DS);
	if (!res) fNSmoothingFailedBins[isect]++;
	else vox->flags |= kSmoothDone;
      }
    }
  }
  // now subtract the dX contribution to DZ
  for (int ix=0;ix<fNXBins;ix++) {
    if (GetXBinIgnored(isect,ix)) continue;
    for (int ip=0;ip<fNY2XBins;ip++) {
      for (int iz=0;iz<fNZ2XBins;iz++) {
	int binGlo = GetVoxGBin(ix,ip,iz);
	bres_t *vox = &sectData[binGlo];
	if (!(vox->flags&kSmoothDone)) continue;
	vox->DS[kResZ] += vox->stat[kVoxZ]*vox->DS[kResX]; // remove slope*dx contribution from DZ
	vox->D[kResZ]  += vox->stat[kVoxZ]*vox->DS[kResX]; // remove slope*dx contribution from DZ
      }
    }
  }
  //
  return cnt;
}

//________________________________________________________________
Bool_t AliTPCDcalibRes::GetSmoothEstimate(int isect, float x, float p, float z, Int_t which, float *res, float *deriv)
{
  // get smooth estimate of distortions mentioned in "which" bit pattern for point in sector coordinates (x,y/x,z/x)
  // smoothing results also saved in the fLastSmoothingRes (allow derivative calculation)
  //
  int minPointsDir[kVoxDim]={0}; // min number of points per direction
  //
  const float kTrialStep = 0.5;
  Bool_t doDim[kResDim] = {kFALSE};
  for (int i=0;i<kResDim;i++) {
    doDim[i] = (which&(0x1<<i))>0;
    if (doDim[i]) res[i] = 0;
  }
  //
  // extimate smoothing matrix size and min number of points
  int matSize = kSmtLinDim;
  for (int i=0;i<kVoxDim;i++) {
    minPointsDir[i] = 3; // for pol1 smoothing require at least 3 points 
    if (fSmoothPol2[i]) {
      minPointsDir[i]++;
      matSize++;
    }
  } 
  double cmat[kResDim][kMaxSmtDim*(kMaxSmtDim+1)/2];
  static int maxNeighb = 10*10*10;
  static bres_t **currClus = new bres_t*[maxNeighb];
  static float* currCache = new float[maxNeighb*kVoxHDim];
  //
  //loop over neighbours which can contribute
  //
  //
  int ix0,ip0,iz0;
  FindVoxel(x,p, isect<kNSect ? z : -z, ix0,ip0,iz0); // find nearest voxel
  bres_t* sectData = fSectGVoxRes[isect];
  int binCen = GetVoxGBin(ix0,ip0,iz0);  // global bin of nearest voxel
  bres_t* voxCen = &sectData[binCen]; // nearest voxel
  //
  int maxTrials[kVoxDim];
  maxTrials[kVoxZ] = fNBins[kVoxZ]/2;
  maxTrials[kVoxF] = fNBins[kVoxF]/2;
  maxTrials[kVoxX] = fMaxBadXBinsToCover*2;

  int trial[kVoxDim]={0};
  while(1)  {
    //
    memset(fLastSmoothingRes,0,kResDim*kMaxSmtDim*sizeof(double));
    memset(cmat,0,kResDim*kMaxSmtDim*(kMaxSmtDim+1)/2*sizeof(double));
    //
    int nbOK=0; // accounted neighbours
    //
    float stepX = fStepKern[kVoxX]*(1. + kTrialStep*trial[kVoxX]);
    float stepF = fStepKern[kVoxF]*(1. + kTrialStep*trial[kVoxF]);
    float stepZ = fStepKern[kVoxZ]*(1. + kTrialStep*trial[kVoxZ]);
    //
    if (!(voxCen->flags&kDistDone) || (voxCen->flags&kMasked) || GetXBinIgnored(isect,ix0)) { 
      // closest voxel has no data, increase smoothing step
      stepX+=kTrialStep*fStepKern[kVoxX];
      stepF+=kTrialStep*fStepKern[kVoxF];
      stepZ+=kTrialStep*fStepKern[kVoxZ];
    }
    //
    // effective kernel widths accounting for the increased bandwidth at the edges and missing data
    float kWXI = GetDXI(ix0)  *fKernelWInv[kVoxX]*fStepKern[kVoxX]/stepX;
    float kWFI = GetDY2XI(ix0)*fKernelWInv[kVoxF]*fStepKern[kVoxF]/stepF;
    float kWZI = GetDZ2XI()   *fKernelWInv[kVoxZ]*fStepKern[kVoxZ]/stepZ;
    int istepX = TMath::Nint(stepX+0.5);
    int istepF = TMath::Nint(stepF+0.5);
    int istepZ = TMath::Nint(stepZ+0.5);
    // for edge bins increase kernel size and neighbours search
    int ixMn=ix0-istepX,ixMx=ix0+istepX;
    if (ixMn<0) {
      ixMn = 0;
      ixMx = TMath::Min(TMath::Nint(ix0+stepX*fKernelScaleEdge[kVoxX]),fNXBins-1);
      kWXI /= fKernelScaleEdge[kVoxX];
    }
    if (ixMx>=fNXBins) {
      ixMx = fNXBins-1;
      ixMn = TMath::Max(TMath::Nint(ix0-stepX*fKernelScaleEdge[kVoxX]),0);
      kWXI /= fKernelScaleEdge[kVoxX];
    }
    //
    int ipMn=ip0-istepF,ipMx=ip0+istepF;
    if (ipMn<0) {
      ipMn = 0;
      ipMx = TMath::Min(TMath::Nint(ip0+stepF*fKernelScaleEdge[kVoxF]),fNY2XBins-1);
      kWFI /= fKernelScaleEdge[kVoxF];
    }
    if (ipMx>=fNY2XBins) {
      ipMx = fNY2XBins-1; 
      ipMn = TMath::Max(TMath::Nint(ip0-stepF*fKernelScaleEdge[kVoxF]),0);
      kWFI /= fKernelScaleEdge[kVoxF];
    }
    //
    int izMn=iz0-istepZ,izMx=iz0+istepZ;
    if (izMn<0) {
      izMn = 0;
      izMx = TMath::Min(TMath::Nint(iz0+stepZ*fKernelScaleEdge[kVoxZ]),fNZ2XBins-1);
      kWZI /= fKernelScaleEdge[kVoxZ];
    }
    if (izMx>=fNZ2XBins) {
      izMx = fNZ2XBins-1;
      izMn = TMath::Max(TMath::Nint(iz0-stepZ*fKernelScaleEdge[kVoxZ]),0);
      kWZI /= fKernelScaleEdge[kVoxZ];
    }
    //
    int nOccX[ixMx-ixMn+1],nOccF[ipMx-ipMn+1],nOccZ[izMx-izMn+1];
    //
    // check if cache arrays should be expanded
    int nbCheck = (ixMx-ixMn+1)*(ipMx-ipMn+1)*(izMx-izMn+1);
    if (nbCheck>=maxNeighb) { // need to expand caches
      int mxNb = nbCheck+100;
      delete[] currClus;
      delete[] currCache;
      currClus = new bres_t*[mxNb];
      currCache= new float[mxNb*kVoxHDim];
      maxNeighb = mxNb;
    }
    //
    for (int i=ixMx-ixMn+1;i--;) nOccX[i]=0;
    for (int i=ipMx-ipMn+1;i--;) nOccF[i]=0;
    for (int i=izMx-izMn+1;i--;) nOccZ[i]=0;
    double u2Vec[3];
    //1st loop, check presence of enough points, cache precalculated values
    float *cacheVal = currCache;
    for (int ix=ixMn;ix<=ixMx;ix++) {
      for (int ip=ipMn;ip<=ipMx;ip++) {
	for (int iz=izMn;iz<=izMx;iz++) {
	  //
	  int binNb = GetVoxGBin(ix,ip,iz);  // global bin
	  bres_t* voxNb = &sectData[binNb];
	  if (!(voxNb->flags&kDistDone) || (voxNb->flags&kMasked) || GetXBinIgnored(isect,ix)) continue; // skip voxels w/o data
	  // estimate weighted distance
	  float dx = voxNb->stat[kVoxX]-x;
	  float df = voxNb->stat[kVoxF]-p;
	  float dz = voxNb->stat[kVoxZ]-z;
	  float dxw = dx*kWXI, dfw = df*kWFI, dzw = dz*kWZI;
	  u2Vec[0] = dxw*dxw;
	  u2Vec[1] = dfw*dfw;
	  u2Vec[2] = dzw*dzw;
	  double kernW = GetKernelWeight(u2Vec,3);
	  if (kernW<kZeroK) continue; 
	  // new point is validated	  
	  nOccX[ix-ixMn]++;
	  nOccF[ip-ipMn]++;
	  nOccZ[iz-izMn]++;
	  currClus[nbOK] = voxNb;
	  cacheVal[kVoxX] = dx;
	  cacheVal[kVoxF] = df;
	  cacheVal[kVoxZ] = dz;
	  cacheVal[kVoxV] = kernW;
	  cacheVal += kVoxHDim;
	  nbOK++;
	}
      }
    }
    //
    // check if we have enough points in every dimension
    int np[kVoxDim]={0};
    for (int i=ixMx-ixMn+1;i--;) if (nOccX[i]) np[kVoxX]++; 
    for (int i=ipMx-ipMn+1;i--;) if (nOccF[i]) np[kVoxF]++;
    for (int i=izMx-izMn+1;i--;) if (nOccZ[i]) np[kVoxZ]++;
    Bool_t enoughPoints = kTRUE, incrDone[kVoxDim] = {0};
    for (int i=0;i<kVoxDim;i++) {
      if (np[i]<minPointsDir[i]) { // need to extend smoothing neighborhood
	enoughPoints=kFALSE;
	if (trial[i]<maxTrials[i] && !incrDone[i]) { //try to increment only missing direction
	  trial[i]++; incrDone[i]=kTRUE;
	} 
	else if (trial[i]==maxTrials[i]) { // cannot increment missing direction, try others
	  for (int j=kVoxDim;j--;) {
	    if (i!=j && trial[j]<maxTrials[j] && !incrDone[j]) {
	      trial[j]++; incrDone[j]=kTRUE;
	    }
	  }
	}
      }
    }
    if (!enoughPoints) {
      if (!(incrDone[kVoxX]||incrDone[kVoxF]||incrDone[kVoxZ])) {
	AliErrorF("Voxel Z:%d F:%d X:%d Trials limit reached: Z:%d F:%d X:%d",
		  voxCen->bvox[kVoxZ],voxCen->bvox[kVoxF],voxCen->bvox[kVoxX],
		  trial[kVoxZ],trial[kVoxF],trial[kVoxX]);
	return kFALSE;
      }
      /*
	AliWarningF("Sector:%2d x=%.3f y/x=%.3f z/x=%.3f (iX:%d iY2X:%d iZ2X:%d)\n"
	"not enough neighbours (need min %d) %d %d %d (tot: %d) | Steps: %.1f %.1f %.1f\n"
	"trying to increase filter bandwidth (trialXFZ:%d %d %d)\n",
	isect,x,p,z,ix0,ip0,iz0,2,np[kVoxX],np[kVoxF],np[kVoxZ],nbOK,stepX,stepF,stepZ,
	trial[kVoxX],trial[kVoxF],trial[kVoxZ]);
      */
      continue;
    }
    //
    // now fill the matrices and solve 
    cacheVal = currCache;
    for (int ib=0;ib<nbOK;ib++) {
      double kernW = cacheVal[kVoxV];
      double dx=cacheVal[kVoxX], df=cacheVal[kVoxF], dz=cacheVal[kVoxZ],dx2=dx*dx, df2=df*df, dz2 = dz*dz;
      cacheVal += kVoxHDim;
      const bres_t* voxNb = currClus[ib];
      for (int id=0;id<kResDim;id++) {
	if (!doDim[id]) continue;
	double kernWD = kernW;
	if (fUseErrInSmoothing) kernWD /= (voxNb->E[id]*voxNb->E[id]); // apart from the kernel value, account for the point error
	double *cmatD = cmat[id];
	double *rhsD = &fLastSmoothingRes[id*kMaxSmtDim];
	//
	double kernWDx=kernWD*dx, kernWDf=kernWD*df, kernWDz=kernWD*dz;
	double kernWDx2=kernWDx*dx, kernWDf2=kernWDf*df, kernWDz2=kernWDz*dz;
	//
	// linear part
	int el=-1,elR=-1;
	cmatD[++el] += kernWD;
	rhsD[++elR] += kernWD*voxNb->D[id];
	//
	cmatD[++el] += kernWDx;   cmatD[++el] += kernWDx2;
	rhsD[++elR] += kernWDx*voxNb->D[id];
	//
	cmatD[++el] += kernWDf;   cmatD[++el] += kernWDx*df;  cmatD[++el] += kernWDf2;
	rhsD[++elR] += kernWDf*voxNb->D[id];
	//
	cmatD[++el] += kernWDz;   cmatD[++el] += kernWDx*dz;  cmatD[++el] += kernWDf*dz;   cmatD[++el] += kernWDz2;
	rhsD[++elR] += kernWDz*voxNb->D[id];	      
	//
	// check if quadratic part is needed
	if (fSmoothPol2[kVoxX]) {
	  cmatD[++el] += kernWDx2;   cmatD[++el] += kernWDx2*dx; cmatD[++el] += kernWDf*dx2; cmatD[++el] += kernWDz*dx2; cmatD[++el] += kernWDx2*dx2;
	  rhsD[++elR] += kernWDx2*voxNb->D[id];
	}
	if (fSmoothPol2[kVoxF]) {
	  cmatD[++el] += kernWDf2;   cmatD[++el] += kernWDx*df2; cmatD[++el] += kernWDf*df2; cmatD[++el] += kernWDz*df2; cmatD[++el] += kernWDx2*df2; cmatD[++el] += kernWDf2*df2;
	  rhsD[++elR] += kernWDf2*voxNb->D[id];
	}
	if (fSmoothPol2[kVoxZ]) {
	  cmatD[++el] += kernWDz2;   cmatD[++el] += kernWDx*dz2; cmatD[++el] += kernWDf*dz2; cmatD[++el] += kernWDz*dz2; cmatD[++el] += kernWDx2*dz2; cmatD[++el] += kernWDf2*dz2; cmatD[++el] += kernWDz2*dz2;
	  rhsD[++elR] += kernWDz2*voxNb->D[id];
	}
      }
    }
    //
    Bool_t fitRes = kTRUE;
    //
    // solve system of linear equations
    AliSymMatrix mat(matSize);
    for (int id=0;id<kResDim;id++) {
      if (!doDim[id]) continue;
      mat.Reset();
      double *cmatD = cmat[id];
      double *rhsD = &fLastSmoothingRes[id*kMaxSmtDim];
      int el=-1,elR=-1,row=-1;
      mat(++row,0) = cmatD[++el];
      mat(++row,0) = cmatD[++el];   mat(row,1) = cmatD[++el];
      mat(++row,0) = cmatD[++el];   mat(row,1) = cmatD[++el];  mat(row,2) = cmatD[++el]; 
      mat(++row,0) = cmatD[++el];   mat(row,1) = cmatD[++el];  mat(row,2) = cmatD[++el];  mat(row,3) = cmatD[++el];
      // pol2 elements if needed
      if (fSmoothPol2[kVoxX]) {
	const int colLim = (++row)+1;
	for (int col=0;col<colLim;col++) mat(row,col) = cmatD[++el];
      }
      if (fSmoothPol2[kVoxF]) {
	const int colLim = (++row)+1;
	for (int col=0;col<colLim;col++) mat(row,col) = cmatD[++el];
      }
      if (fSmoothPol2[kVoxZ]) {
	const int colLim = (++row)+1;
	for (int col=0;col<colLim;col++) mat(row,col) = cmatD[++el];
      }
      //
      fitRes &= mat.SolveChol(rhsD);
      if (!fitRes) {
	for (int i=kVoxDim;i--;) trial[i]++;
	AliWarningF("Sector:%2d x=%.3f y/x=%.3f z/x=%.3f (iX:%d iY2X:%d iZ2X:%d)\n"
		    "neighbours range used %d %d %d (tot: %d) | Steps: %.1f %.1f %.1f\n"
		    "Solution for smoothing Failed, trying to increase filter bandwidth (trialXFZ: %d %d %d)",
		    isect,x,p,z,ix0,ip0,iz0,np[kVoxX],np[kVoxF],np[kVoxZ],nbOK,stepX,stepF,stepZ,trial[kVoxX],trial[kVoxF],trial[kVoxZ]);
	continue;
      }
      res[id] = rhsD[0];
      if (deriv) for (int j=0;j<3;j++) deriv[id*3 +j] = rhsD[j+1]; // ignore eventual pol2 term
    }
    //
    break; // success
  } // end of loop over allowed trials
  
  return kTRUE;

}


//_____________________________________
Double_t AliTPCDcalibRes::GetKernelWeight(double* u2vec,int np) const
{
  double w = 1;
  if (fKernelType == kEpanechnikovKernel) {
    for (int i=np;i--;) {
      if (u2vec[i]>1) return 0.;
      w *= 3./4.*(1.-u2vec[i]);
    }
  }
  else if (fKernelType == kGaussianKernel) {
    double u2 = 0;
    for (int i=np;i--;) u2 += u2vec[i];
    w = u2<kMaxGaussStdDev*kMaxGaussStdDev*np ? TMath::Exp(-u2)/TMath::Sqrt(2.*TMath::Pi()) : 0;
  }
  else {
    AliFatalF("Kernel type %d is not defined",fKernelType);
  }
  return w;
}


//_____________________________________
void AliTPCDcalibRes::SetKernelType(int tp, float bwX, float bwP, float bwZ, float scX,float scP,float scZ)
{
  // set kernel type and widths in terms of binning in X,Y/X and Z/X, define aux variables
  fKernelType = tp;
  //
  fKernelScaleEdge[kVoxX] = scX;
  fKernelScaleEdge[kVoxF] = scP;
  fKernelScaleEdge[kVoxZ] = scZ;

  fKernelWInv[kVoxX] = bwX>0 ? 1./bwX : 1.;
  fKernelWInv[kVoxF] = bwP>0 ? 1./bwP : 1.;
  fKernelWInv[kVoxZ] = bwZ>0 ? 1./bwZ : 1.;

  if (fKernelType == kEpanechnikovKernel) { // bandwidth 1
    fStepKern[kVoxX] = TMath::Nint(bwX+0.5);
    fStepKern[kVoxF] = TMath::Nint(bwP+0.5);
    fStepKern[kVoxZ] = TMath::Nint(bwZ+0.5);    
  }
  else if (kGaussianKernel) {  // look in ~5 sigma
    fStepKern[kVoxX] = TMath::Nint(bwX*5.+0.5);
    fStepKern[kVoxF] = TMath::Nint(bwP*5.+0.5);
    fStepKern[kVoxZ] = TMath::Nint(bwZ*5.+0.5);
  }
  else {
    AliFatalF("Kernel type %d is not defined",fKernelType);
  }
  for (int i=kVoxDim;i--;) if (fStepKern[i]<1) fStepKern[i] = 1;
  //
}


//=========================================================================

//_____________________________________________
void AliTPCDcalibRes::CreateCorrectionObject()
{
  // create correction object for given time slice

  AliSysInfo::AddStamp("CreateCorrectionObject",0,0,0,0);
  //  
  // check if there are failures
  Bool_t create = kTRUE;
  for (int i=0;i<kNSect2;i++) {
    if (fNSmoothingFailedBins[i]>0) {
      AliErrorF("%d failed voxels in sector %d",fNSmoothingFailedBins[i],i); 
      create = kFALSE;
    }
  }
  //
  if (!create) {
    AliError("ATTENTION: MAP WILL NOT BE CREATED");
    return;
  }
  //
  TString name = Form("run%d_%lld_%lld",fRun,fTMin,fTMax);
  fChebCorr = new AliTPCChebCorr(name.Data(),name.Data(),
				 fChebPhiSlicePerSector,fChebZSlicePerSide,1.0f);
  fChebCorr->SetUseFloatPrec(kFALSE);
  fChebCorr->SetRun(fRun);
  fChebCorr->SetTimeStampStart(fTMin);
  fChebCorr->SetTimeStampEnd(fTMax);
  fChebCorr->SetTimeDependent(kFALSE);
  fChebCorr->SetUseZ2R(kTRUE);
  //
  if      (fBz> 0.01) fChebCorr->SetFieldType(AliTPCChebCorr::kFieldPos);
  else if (fBz<-0.01) fChebCorr->SetFieldType(AliTPCChebCorr::kFieldNeg);
  else                fChebCorr->SetFieldType(AliTPCChebCorr::kFieldZero);
  // Note: to create universal map, set manually SetFieldType(AliTPCChebCorr::kFieldAny)

  SetUsedInstance(this);
  int npCheb[kResDim][2];
  for (int i=0;i<kResDim;i++) { // do we need to determine N nodes automatically?
    int nbFauto = TMath::Max(int(fNY2XBins*1.2),fNY2XBins+3);
    int nbZauto = TMath::Max(int(fNZ2XBins*1.2),fNZ2XBins+3);
    npCheb[i][0] = fNPCheb[i][0];
    npCheb[i][1] = fNPCheb[i][1];
    if (npCheb[i][0]<1) {
      npCheb[i][0] = nbFauto; // 1st dimension: sector coordinate y/x
      AliInfoF("Nnodes for Cheb.%4s segmentation in %4s is set to %2d",kResName[i],kVoxName[kVoxF],nbFauto);
    }
    if (npCheb[i][1]<1) {
      npCheb[i][1] = nbZauto; // 2nd dimension: z/x
      AliInfoF("Nnodes for Cheb.%4s segmentation in %4s is set to %2d",kResName[i],kVoxName[kVoxZ],nbZauto);
    }
  }
  fChebCorr->Parameterize(trainCorr,kResDim,npCheb,fChebPrecD);
  //
  // register tracks rate for lumi weighting
  fChebCorr->SetTracksRate(ExtractTrackRate());
  //
  AliSysInfo::AddStamp("CreateCorrectionObject",1,0,0,0);
}

//________________________________________________________________
TH1F* AliTPCDcalibRes::ExtractTrackRate() const
{
  // create histo with used tracks per timestamp
  TH1F* hTr = 0;
  if (fTracksRate) { 
    const float *tarr = fTracksRate->GetArray();
    int nb = fTracksRate->GetNbinsX();
    int bin0=1,bin1=nb;
    while(tarr[bin0]<0.5 && bin0<=nb) bin0++; // find first significant time bin 
    while(tarr[bin1]<0.5 && bin1>=bin0) bin1--; // find last significant time bin
    nb = bin1 - bin0 + 1;
    Long64_t tmn = (Long64_t)fTracksRate->GetBinCenter(bin0);
    Long64_t tmx = (Long64_t)fTracksRate->GetBinCenter(bin1);
    hTr = new TH1F(Form("TrackRate%lld_%lld",tmn,tmx),"TracksRate",nb,
		   fTracksRate->GetBinLowEdge(bin0),fTracksRate->GetBinLowEdge(bin1+1));
    for (int ib=0;ib<nb;ib++) {
      int ibh = ib+bin0;
      hTr->SetBinContent(ib+1,tarr[ibh]);
    }
  }
  else AliError("TracksRate accumulation histo was not initialized");
  return hTr;
}

//________________________________________________________________
void AliTPCDcalibRes::InitBinning()
{
  // initialize binning structures
  //
  // X binning
  if (fNXBins>0 && fNXBins<kNPadRows) {
    AliInfoF("X-binning: uniform %d bins from %.2f to %.2f",fNXBins,kMinX,kMaxX);
    fDXI         = fNXBins/(kMaxX-kMinX);
    fDX          = 1.0f/fDXI;
    fUniformBins[kVoxX] = kTRUE;
  }
  else {
    fNXBins = kNPadRows;
    AliInfo("X-binning: bin per pad-row");
    fUniformBins[kVoxX] = kFALSE;
    fDX = kTPCRowDX[0];
    fDXI = 1.f/fDX; // should not be used
  }
  //
  // Y binning
  if (fNXBins<1) AliFatal("X bins must be initialized first");
  fMaxY2X = new Float_t[fNXBins];        // max Y/X at each X bin, account for dead zones
  fDY2XI  = new Float_t[fNXBins];        // inverse of Y/X bin size at given X bin
  fDY2X   = new Float_t[fNXBins];        // Y/X bin size at given X bin
  //
  const float kMaxY2X = TMath::Tan(0.5f*kSecDPhi);

  for (int ix=0;ix<fNXBins;ix++) {
    float x = GetX(ix);
    fMaxY2X[ix] = kMaxY2X - kDeadZone/x;
    fDY2XI[ix] = fNY2XBins / (2.f*fMaxY2X[ix]);
    fDY2X[ix] = 1.f/fDY2XI[ix];
  }
  //
  fDZ2XI = fNZ2XBins/kMaxZ2X;
  fDZ2X  = 1.0f/fDZ2XI;
  //
  fNBins[kVoxX] = fNXBins;
  fNBins[kVoxF] = fNY2XBins;
  fNBins[kVoxZ] = fNZ2XBins;

  fNGVoxPerSector = fNY2XBins*fNZ2XBins*fNXBins;

}

//________________________________________________________________
Int_t AliTPCDcalibRes::GetXBin(float x) 
{
  // convert X to bin ID, following pad row widths
  if (fUniformBins[kVoxX]) {
    int ix = (x-kMinX)*fDXI;
    if (ix<0) return 0;
    else if (ix>=fNXBins) return fNXBins-1;
    return ix;
  }
  else {
    int ix;
    if (x<kTPCRowX[kNRowIROC-1]+0.5*kTPCRowDX[kNRowIROC-1]) {     // uniform pad size in IROC
      ix = (x-(kTPCRowX[0]-kTPCRowDX[0]*0.5))/kTPCRowDX[0];
      if (ix<0) ix = 0;
    }
    // uniform pad size in OROC2
    else if ( x>= kTPCRowX[kNRowIROC+kNRowOROC1]-0.5*kTPCRowDX[kNRowIROC+kNRowOROC1] ) {
      ix = (x-(kTPCRowX[kNRowIROC+kNRowOROC1]-0.5*kTPCRowDX[kNRowIROC+kNRowOROC1]))/kTPCRowDX[kNPadRows-1]
	+ kNRowIROC + kNRowOROC1;
      if (ix>=kNPadRows) ix = kNPadRows-1;
    }
    else { // uniform pad size in OROC1
      ix = (x-(kTPCRowX[kNRowIROC]-0.5*kTPCRowDX[kNRowIROC]))/kTPCRowDX[kNRowIROC] + kNRowIROC;
      if (ix<kNRowIROC) { // we go between IROC and OROC?
	if (x> 0.5*(kTPCRowX[kNRowIROC-1]+kTPCRowX[kNRowIROC]))  ix = kNRowIROC; // 1st OROC1 row
	else ix = kNRowIROC-1;
      }
      
    }
    return ix;
    // 
  }
}

//________________________________________________________________
Int_t AliTPCDcalibRes::GetRowID(float x)
{
  // return row ID
  int ix;
  if (x<kTPCRowX[kNRowIROC-1]+0.5*kTPCRowDX[kNRowIROC-1]) {     // uniform pad size in IROC
    ix = (x-(kTPCRowX[0]-kTPCRowDX[0]*0.5))/kTPCRowDX[0];
    if (ix<0) ix = -1;
  }
  // uniform pad size in OROC2
  else if ( x>= kTPCRowX[kNRowIROC+kNRowOROC1]-0.5*kTPCRowDX[kNRowIROC+kNRowOROC1] ) {
    ix = (x-(kTPCRowX[kNRowIROC+kNRowOROC1]-0.5*kTPCRowDX[kNRowIROC+kNRowOROC1]))/kTPCRowDX[kNPadRows-1]
      + kNRowIROC + kNRowOROC1;
    if (ix>=kNPadRows) ix = -2;
  }
  else { // uniform pad size in OROC1
    ix = (x-(kTPCRowX[kNRowIROC]-0.5*kTPCRowDX[kNRowIROC]))/kTPCRowDX[kNRowIROC] + kNRowIROC;
    if (ix<kNRowIROC) { // we go between IROC and OROC?
      ix = -3;
    }
    
  }
  return ix;
}

//_____________________________________________________
Bool_t AliTPCDcalibRes::FindVoxelBin(int sectID, float x, float y, float z, UChar_t bin[kVoxHDim],float voxVars[kVoxHDim])
{
  // define voxel variables and bin
  //  
  //
  // track Z/X bin
  voxVars[kVoxZ] = z/x;
  if (TMath::Abs(voxVars[kVoxZ])>kMaxZ2X) return kFALSE;
  int binZ       = GetZ2XBinExact( sectID<kNSect ? voxVars[kVoxZ] : -voxVars[kVoxZ] );
  if (binZ<0) return kFALSE;
  bin[kVoxZ] = binZ;
  //
  // track X bin
  voxVars[kVoxX] = x;
  int binX = GetXBinExact(x); // exact matching to pad-rows
  if (binX<0 || binX>=fNXBins) return kFALSE;
  bin[kVoxX] = binX;
  //
  // track Y/X bin accounting for sector edges dead zones
  voxVars[kVoxF] = y/x;
  int binF = GetY2XBinExact(voxVars[kVoxF],binX);
  if (binF<0||binF>=fNY2XBins) return kFALSE;
  bin[kVoxF] = binF;
  //
  return kTRUE;
}

//=====================================================================
///////////////////////////////////////////////////
//
// This is a special non-meber f-n for cheb. learning
//
void trainCorr(int row, float* tzLoc, float* corrLoc)
{
  // Cheb. object training f-n: compute correction for the point
  //
  // xtzLoc: y2x and z2x in sector frame
  // corrLoc: x, y, z corrections in sector frame
  // when called with pointer at 0, row will set the sector/side 
  // (row should be sector in 0-35 or 0-71 format)
  static int sector=0;
  if (!tzLoc || !corrLoc) {
    sector = row%AliTPCDcalibRes::kNSect2; 
    printf("training Sector%d\n",sector);
    return;
  }
  //
  AliTPCDcalibRes* calib = AliTPCDcalibRes::GetUsedInstance();

  float x = AliTPCDcalibRes::GetTPCRowX(row);
  float dist[AliTPCDcalibRes::kResDim] = {0};
  int xbin = calib->GetNXBins()==AliTPCDcalibRes::kNPadRows ? row : calib->GetXBin(x);
  if (calib->GetXBinIgnored(sector,xbin)) return;
  float y2x = tzLoc[0];
  float z2x = tzLoc[1];
  //
  Bool_t res = calib->GetSmoothEstimate(sector, x, y2x, z2x, 0xff, dist);
  if (!res) { printf("Failed to evaluate smooth distortion\n"); exit(1); }

  /*
  // Marian stored Z track coordinate instead of cluster one, need to correct for this
  if (fApplyZt2Zc) {
    float deriv[AliTPCDcalibRes::kResDim*3];
    const double inversionEps = 20e-4; // when inverting, stop Newton-Raphson iterations at this eps
    const int    inversionMaxIt = 3; // when inverting, stop Newton-Raphson after some numbers of iterations
    double change = 0, xInv = 1./x;
    float zc = z2x*x;
    float zt = zc + dist[kResZ]; // 1st guess on true zt at measured zc
    //
    // use Newton-Raphson method for NDLocal inversion to get zt = F(zc) from 
    // dz == zt - zc = NDLoc(zt)   ->  zc - [zt - NDLoc(zt)]=0 
    // ->  zt_{i+1} = zt_i - (zc - [zt - NDLoct(zt_i)])/(-d[zt_i-NDLoc(zt_i)]/dz)
    //
    int it = 0;
    do {
      z2x = zt*xInv;
      res = GetSmoothEstimate(sector, x, y2x, z2x, dist, deriv);
      if (!res) {printf("Failed to evaluate smooth distortion\n");exit(1);}
      double bot = 1. - deriv[kResZ*3+2]*xInv;  // dF(zt_i)/dz
      if (TMath::Abs(bot)<1e-6) break;
      double top = zc - (zt - dist[kResZ]); // zc - F(zt_i) 
      double change = top/bot;
      //  printf("It %d Eps:%+e, Zc:%+e, Zt:%+e Ztn:%+e DZ:%+e | dztmp: :%+e DD:%+e\n",
      //     it,change,zc,zt,zt+change,dz,dztmp,deriv[kndZ2R]*xInv);
      zt += change;
      it++;
    } while(it<inversionMaxIt && TMath::Abs(change)>inversionEps);
    //
    // now query at fixed Z2X
    z2x = zt*xInv;
    res = GetSmoothEstimate(sector, x, y2x, z2x, dist);
    if (!res) {printf("Failed to evaluate smooth distortion\n");exit(1);}
  }
  */

  for (int i=0;i<AliTPCDcalibRes::kResDim;i++) corrLoc[i] = dist[i];
  //
}
//======================================================================================
