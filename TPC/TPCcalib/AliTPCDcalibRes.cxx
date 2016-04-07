#include "AliTPCDcalibRes.h"

// this must be standalone f-n, since the signature is important for Chebyshev training
void trainCorr(int row, float* tzLoc, float* corrLoc);


const char* AliTPCDcalibRes::kVoxName[AliTPCDcalibRes::kVoxDim] = {"tgSlp","y2x","x","z2x"};
const char* AliTPCDcalibRes::kResName[AliTPCDcalibRes::kResDim] = {"dX","dY","dZ","Disp"};
const char* AliTPCDcalibRes::kEstName[AliTPCDcalibRes::kNEstPar] = {
  "Nrm","Mean","Sig","Max",
  "MeanL","MeanEL","SigL","SigEL",
  "NormG","MeanG","MeanEG","SigG","SigEG","Chi2G"};

const float AliTPCDcalibRes::kSecDPhi = 20.f*TMath::DegToRad();
const float AliTPCDcalibRes::kMaxResid = 7.0f;
const float AliTPCDcalibRes::kMinX = 85.0f;
const float AliTPCDcalibRes::kMaxX = 246.0f;
const float AliTPCDcalibRes::kMaxZ2X = 1.0f;
const float AliTPCDcalibRes::kZLim = 250.0f;
const char* AliTPCDcalibRes::kLocalResFileName  = "tmpDeltaSect";
const char* AliTPCDcalibRes::kClosureTestFileName  = "closureTestSect";
const char* AliTPCDcalibRes::kStatOut      = "voxelStat";
const char* AliTPCDcalibRes::kResOut       = "voxelRes";
const char* AliTPCDcalibRes::kDriftFileName= "fitDrift";
const float AliTPCDcalibRes::kDeadZone = 1.5;
const float AliTPCDcalibRes::kZeroK = 1e-6;
const float AliTPCDcalibRes::kInvalidR = 10.f;
const float AliTPCDcalibRes::kInvalidRes = -900;
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
  ,fChebPhiSlicePerSector(2)
  ,fChebCorr(0)

  ,fRun(run)
  ,fTMin(tmin)
  ,fTMax(tmax)
  ,fMaxTracks(9999999)
  ,fCacheInp(100)
  ,fLearnSize(1)
  ,fBz(0)
  ,fDeleteSectorTrees(kFALSE) // set to true for production
  ,fResidualList(resList)
  ,fOCDBPath()

  ,fNPrimTracksCut(600)
  ,fMinNCl(30)
  ,fMaxDevYHelix(0.3)
  ,fMaxDevZHelix(0.3) // !!! VDrift calib. screas up the Z fit, 0.3 w/o vdrift
  ,fNVoisinMA(3)
  ,fNVoisinMALong(15)
  ,fMaxStdDevMA(25.0)
  ,fMaxRMSLong(0.8)
  ,fMaxRejFrac(0.15)
  ,fFilterOutliers(kTRUE) 

  ,fMaxDY(6.f)
  ,fMaxDZ(6.f)
  ,fMaxQ2Pt(3.f)
  ,fMidQ2Pt(1.22f)
  ,fNY2XBins(15)
  ,fNZ2XBins(10)
  ,fNXBins(-1)
  ,fNXYBinsProd(0)
  ,fNDeltaYBins(120)
  ,fNDeltaZBins(120)  
  ,fDZ2X(0)
  ,fDX(0)
  ,fDZ2XI(0)
  ,fDXI(0)
  ,fDeltaYbinI(0)
  ,fDeltaZbinI(0)
  ,fNGVoxPerSector(0)
  //
  ,fMaxY2X(0)
  ,fDY2X(0)
  ,fDY2XI(0)
//  ,fBinMinQ(0) // obsolete binning
//  ,fBinDQ(0)
//  ,fBinDQI(0)

  ,fNMaxNeighb(0)
  ,fKernelType(kGaussianKernel)
  

  ,fNTrSelTot(0)
  ,fNTrSelTotWO(0)
  ,fNReadCallTot(0)
  ,fNBytesReadTot(0)

  ,fVDriftParam(0)
  ,fVDriftGraph(0)
  ,fCorrTime(0)

  ,fStatTree(0)
  ,fHDelY(0)
  ,fHDelZ(0)

  ,fDTS()
  ,fDTC()

  ,fTimeStamp(0)
  ,fNCl(0)
  ,fQ2Pt(0)
  ,fTgLam(0)

{
  for (int i=0;i<kResDim;i++) {
    for (int j=0;j<2;j++) fNPCheb[i][j] = 15;
    fChebPrecD[i] = 100e-4;
  }

  memset(fLastSmoothingRes,0,kResDim*4*sizeof(double));

  for (int i=kVoxDim;i--;) {
    fUniformBins[i] = kTRUE;
    fKernelScaleEdge[i] = 1.0f;
    fStepKern[i] = 1;
    fKernelWInv[i] = 0; // calculated later
  }

  for (int i=kVoxHDim;i--;) {
    fNBProdSt[i] = 0;
    fNBProdDY[i] = 0;
    fNBProdDZ[0] = 0;
  }
  fNBProdSectG[0] = fNBProdSectG[1] = 0;
  //
  for (int i=0;i<kNSect2;i++) {
    fSectGVoxRes[i] = 0;
    fTmpTree[i] = 0;
    fStatHist[i] = 0;
    fArrNDStat[i] = 0;
    fTmpFile[i] = 0;
  }
  for (int i=kNQBins+1;i--;) fQ2PTBound[i] = 0;
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
  //  delete[] fBinMinQ;
  //  delete[] fBinDQ;
  //  delete[] fBinDQI;
  delete fVDriftParam;
  delete fVDriftGraph;
  delete fHDelY;
  delete fHDelZ;
  for (int i=0;i<kNSect2;i++) {
    delete fSectGVoxRes[i];
    delete fStatHist[i];
  }
}

//________________________________________
void AliTPCDcalibRes::ProcessFromDeltaTrees()
{
  // process from residual trees
  TStopwatch sw;
  // select tracks matching to time window and write compact local trees
  CollectData(kExtractMode);
  //
  ProcessFromLocalBinnedTrees();
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
  //
}

//________________________________________
void AliTPCDcalibRes::ProcessFromLocalBinnedTrees()
{
  // process starting from local binned trees created by CollectData(kExtractMode)
  TStopwatch sw;
  sw.Start();

  // do per-sector projections and fits
  ProcessResiduals();
  // store treee with voxels definitions
  WriteVoxelDefinitions();
  //
  ProcessFromStatTree();
  //
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
}

//________________________________________
void AliTPCDcalibRes::ProcessFromStatTree()
{
  // process starting from the voxels statistics tree created by the ProcessFromLocalBinnedTrees
  TStopwatch sw;
  sw.Start();

  ExtractXYZDistortions();
  //
  ProcessDispersions();
  //
  CreateCorrectionObject();
  //
  WriteResTree();
  //
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
}

//________________________________________
void AliTPCDcalibRes::Save(const char* name)
{
  // save itself
  TString names = name;
  if (names.IsNull()) {
    names = Form("%s_run%d_%lld_%lld.root",IsA()->GetName(),fRun,fTMin,fTMax);
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
  InitGeom();
  SetName(Form("run%d_%lld_%lld",fRun,fTMin,fTMax));
  SetTitle(IsA()->GetName());
  //
  if (fMidQ2Pt<0) fMidQ2Pt = fMaxQ2Pt/2.f;
  //  
  if (fNDeltaYBins>kMaxResBins) {
    AliErrorF("N DeltaY bins %d exceeds max allowed, setting to %d",fNDeltaYBins,kMaxResBins);
    fNDeltaYBins = kMaxResBins;
  }
  if (fNDeltaZBins>kMaxResBins) {
    AliErrorF("N DeltaZ bins %d exceeds max allowed, setting to %d",fNDeltaZBins,kMaxResBins);
    fNDeltaZBins = kMaxResBins;
  }
  // define boundaries
  InitBinning();
  //
  LoadVDrift(); //!!!
  //
  // prepare aux info for stat and residuals histo bin calculation, see doc of TNDArray bin calculation
  fNBProdDY[kVoxHDim-1] = 1;
  fNBProdDZ[kVoxHDim-1] = 1;
  fNBProdSt[kVoxHDim-1] = 1;
  int nbh[kVoxDim];
  nbh[kVoxQ] = kNQBins;
  nbh[kVoxF] = fNY2XBins;
  nbh[kVoxX] = fNXBins;
  nbh[kVoxZ] = fNZ2XBins;
  for (int i=kVoxDim;i--;) {   // +2 to account for under/over-flows
    fNBProdSt[i] = fNBProdSt[i+1]*(2 + ((i==kVoxDim-1) ? kVoxHDim     : nbh[i+1])); 
    fNBProdDY[i] = fNBProdDY[i+1]*(2 + ((i==kVoxDim-1) ? fNDeltaYBins : nbh[i+1]));
    fNBProdDZ[i] = fNBProdDZ[i+1]*(2 + ((i==kVoxDim-1) ? fNDeltaZBins : nbh[i+1]));
  }
  //
  AliSysInfo::AddStamp("Init",0,0,0,0);
  //
  fInitDone = kTRUE;
}

//_____________________________________________________
void AliTPCDcalibRes::CollectData(int mode) 
{
  const float kEps = 1e-6;
  const float q2ptIniTolerance = 1.5;
  if (!fInitDone) Init();
  //  gEnv->SetValue("TFile.AsyncPrefetching", 1);
  TVectorF *vecDY=0,*vecDZ=0,*vecZ=0,*vecR=0,*vecSec=0,*vecPhi=0, *vecDYITS=0,*vecDZITS=0;
  UShort_t npValid = 0;
  Int_t nPrimTracks = 0;
  Char_t trdOK=0;
  AliExternalTrackParam* param = 0;
  //
  TStopwatch swTot;
  swTot.Start();
  fNTrSelTot = 0;
  fNTrSelTotWO = 0;
  fNReadCallTot = 0;
  fNBytesReadTot = 0;
  //
  CreateLocalResidualsTrees(mode);
  //
  // prepare input tree
  TString  chunkList = gSystem->GetFromPipe(TString::Format("cat %s",fResidualList.Data()).Data());
  TObjArray *chunkArray= chunkList.Tokenize("\n");  
  Int_t nChunks = chunkArray->GetEntriesFast();  
  //
  AliSysInfo::AddStamp("ProjInit",0,0,0,0);

  for (int ichunk=0;ichunk<nChunks;ichunk++) {
    //
    int ntrSelChunkWO=0, ntrSelChunk=0,nReadCallsChunk=0,nBytesReadChunk=0;
    //
    TStopwatch swc;
    swc.Start();
    TString fileNameString(chunkArray->At(ichunk)->GetName());
    if (fileNameString.Contains("alien://") && (!gGrid || (gGrid && !gGrid->IsConnected()))) TGrid::Connect("alien://");
    TFile *chunkFile = TFile::Open(fileNameString.Data());
    if (!chunkFile) continue;
    TTree *tree = (TTree*)chunkFile->Get("delta");
    if (!tree) {AliWarningF("No delta tree in %s",fileNameString.Data());continue;}
    tree->SetCacheLearnEntries(fLearnSize);
    tree->SetCacheSize(0);
    tree->SetCacheSize(fCacheInp*kMByte);
    //
    tree->SetBranchStatus("*",kFALSE);
    if (fNPrimTracksCut>0) tree->SetBranchStatus("nPrimTracks",kTRUE);
    tree->SetBranchStatus("timeStamp",kTRUE);
    tree->SetBranchStatus("trdOK",kTRUE);
    tree->SetBranchStatus("vecR.",kTRUE);
    tree->SetBranchStatus("vecSec.",kTRUE);
    tree->SetBranchStatus("vecPhi.",kTRUE);
    tree->SetBranchStatus("vecZ.",kTRUE);
    tree->SetBranchStatus("track.*",kTRUE);      
    tree->SetBranchStatus("npValid",kTRUE);
    tree->SetBranchStatus("trd0.",kTRUE);
    tree->SetBranchStatus("trd1.",kTRUE);
    tree->SetBranchStatus("its0.",kTRUE);
    tree->SetBranchStatus("its1.",kTRUE);
    //
    tree->SetBranchAddress("timeStamp",&fTimeStamp);
    tree->SetBranchAddress("trdOK",&trdOK);
    if (fNPrimTracksCut>0) tree->SetBranchAddress("nPrimTracks",&nPrimTracks);
    tree->SetBranchAddress("vecR.",&vecR);
    tree->SetBranchAddress("vecSec.",&vecSec);
    tree->SetBranchAddress("vecPhi.",&vecPhi);
    tree->SetBranchAddress("vecZ.",&vecZ);
    tree->SetBranchAddress("track.",&param);
    tree->SetBranchAddress("npValid",&npValid);
    tree->SetBranchAddress("trd0.",&vecDY);
    tree->SetBranchAddress("trd1.",&vecDZ);
    tree->SetBranchAddress("its0.",&vecDYITS);
    tree->SetBranchAddress("its1.",&vecDZITS);
    //
    tree->GetEntry(0);

    TBranch* brTime = tree->GetBranch("timeStamp");
    TBranch* brTRDOK = tree->GetBranch("trdOK");
    //
    int nTracks = tree->GetEntries();
    AliInfoF("Processing %d tracks of %s",nTracks,fileNameString.Data());

    float residHelixY[kNPadRows],residHelixZ[kNPadRows];
    //
    // reset the cache when swithching between the timeStamp and Event read modes
    Bool_t lastReadMatched = kFALSE; 
    for (int itr=0;itr<nTracks;itr++) {
      nBytesReadChunk += brTime->GetEntry(itr);
      if (fTimeStamp<fTMin  || fTimeStamp>fTMax) {
	if (lastReadMatched && fSwitchCache) { // reset the cache
	  tree->SetCacheSize(0);
	  tree->SetCacheSize(fCacheInp*kMByte);
	  lastReadMatched = kFALSE;
	}
	continue;	
      }
      //
      brTRDOK->GetEntry(itr);
      if (!trdOK) continue;
      //
      if (!lastReadMatched && fSwitchCache) { // reset the cache before switching to event reading mode
	tree->SetCacheSize(0);
	tree->SetCacheSize(fCacheInp*kMByte);
      }
      lastReadMatched = kTRUE;
      nBytesReadChunk += tree->GetEntry(itr);
      if (fNPrimTracksCut>0 && nPrimTracks>fNPrimTracksCut) continue;
      //
      fQ2Pt = param->GetParameter()[4];
      fTgLam = param->GetParameter()[3];
      if (TMath::Abs(fQ2Pt)>fMaxQ2Pt*q2ptIniTolerance) continue;
      //
      const Float_t *vSec= vecSec->GetMatrixArray();
      const Float_t *vPhi= vecPhi->GetMatrixArray();
      const Float_t *vR  = vecR->GetMatrixArray();
      const Float_t *vZ  = vecZ->GetMatrixArray();
      const Float_t *vDY = vecDY->GetMatrixArray();
      const Float_t *vDZ = vecDZ->GetMatrixArray();
      const Float_t *vDYITS = vecDYITS->GetMatrixArray();
      const Float_t *vDZITS = vecDZITS->GetMatrixArray();
      //
      fCorrTime = (fVDriftGraph!=NULL) ? fVDriftGraph->Eval(fTimeStamp):0; // for VDrift correction
      //
      fNCl = 0;
      // 1st iteration: collect data in cluster frame
      for (int ip=0;ip<npValid;ip++) { // 1st fill selected track data to buffer for eventual outlier rejection
	if (vR[ip]<kInvalidR || vDY[ip]<kInvalidRes || vDYITS[ip]<kInvalidRes) continue;
	//
	fArrX[fNCl]   = vR[ip];  // X (R) is the same for cluster and track
	fArrZTr[fNCl] = vZ[ip];  // Z of ITS track was stored!!
	fArrDY[fNCl]  = vDY[ip]; // this is also the track coordinate in cluster frame
	fArrDZ[fNCl]  = vDZ[ip];
	fArrPhi[fNCl] = vPhi[ip];
	int rocID = TMath::Nint(vSec[ip]);
	//
	// !!! fArrZTr corresponds to ITS track Z, we need that of TRD-ITS
	fArrZTr[fNCl] += fArrDZ[fNCl] - vDZITS[ip]; // recover ITS-TRD track position from ITS and deltas
	
	if (fFixAlignmentBug && !param->TestBit(kAlignmentBugFixedBit)) {
	  FixAlignmentBug(rocID, fQ2Pt, fBz, fArrPhi[fNCl], fArrX[fNCl], fArrZTr[fNCl], fArrDY[fNCl],fArrDZ[fNCl]);
	}
	if (fArrPhi[fNCl]<0) fArrPhi[fNCl] += 2.*TMath::Pi();
	//
	// calculate drift velocity calibration if available
	float dzDrift = GetDriftCorrection(fArrZTr[fNCl],fArrX[fNCl],fArrPhi[fNCl],rocID);
	// apply drift velocity calibration if available
	fArrDZ[fNCl] += dzDrift;
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
      if (TMath::Abs(fQ2Pt)>fMaxQ2Pt) continue; // now we have more precise estimate of q/pt
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
	
	float sna = TMath::Sin(fArrPhi[ip]-(0.5f +fArrSectID[ip]%kNSect)*kSecDPhi);
	float csa = TMath::Sqrt((1.f-sna)*(1.f+sna));
	//
	// by using propagation in cluster frame in AliTPCcalibAlignInterpolation::Process,
	// the X of the track is evaluated not at the pad-row x=r*csa but at x=r*sca-dy*sna
	double xrow = fArrX[ip]*csa;
	double dx   = fArrDY[ip]*sna;
	double xtr  = xrow - dx;
	double ycl  = fArrX[ip]*sna;      // cluster Y in the sector frame
	double ytr  = ycl + fArrDY[ip]*csa; // track Y in the sector frame at x=xtr is 
	//
	double ztr  = fArrZTr[ip];          // Z of the track at x=xtr
	double zcl  = ztr - fArrDZ[ip];     // and the Z of the cluster is Ztr-deltaZ
	//
	// Now we need to take the track to real pad-row X
	// use linear extrapolation:
	float tgs = fArrTgSlp[ip];
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
	if (TMath::Abs(fArrDY[fNCl])>fMaxDY-kEps) continue;
	if (TMath::Abs(fArrDZ[fNCl])>fMaxDZ-kEps) continue;
	//
	if (fArrX[fNCl]<kMinX || fArrX[fNCl]>kMaxX) continue;
	if (TMath::Abs(fArrZCl[fNCl])>kZLim) continue;;
	//
	// End of manipulations to go to the sector frame
	//
	fNCl++;
      }

      if (fFilterOutliers && !ValidateTrack()) continue;

      ntrSelChunk++;

      if (mode==kExtractMode) {
	FillLocalResidualsTrees();
      }
      else if (mode==kClosureTestMode) {
	FillCorrectedResiduals();
      }
    } // loop over tracks
    //
    swc.Stop();
    nReadCallsChunk =  chunkFile->GetReadCalls();
    AliInfoF("Chunk%3d: selected %d tracks (%d with outliers) from chunk %d | %.1f MB read in %d read calls",
	     ichunk,ntrSelChunk,ntrSelChunkWO, ichunk,float(nBytesReadChunk)/kMByte,nReadCallsChunk); swc.Print();
    fNTrSelTot += ntrSelChunk;
    fNTrSelTotWO += ntrSelChunkWO;
    fNReadCallTot += nReadCallsChunk;
    fNBytesReadTot += nBytesReadChunk;
    //
    delete tree;
    chunkFile->Close();
    delete chunkFile;
    //
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
  for (int is=0;is<kNSect2;is++) {
    fTmpFile[is]->cd();
    fTmpTree[is]->Write("", TObject::kOverwrite);
    delete fTmpTree[is];
    fTmpTree[is] = 0;
    fTmpFile[is]->Close();
    delete fTmpFile[is];
    fTmpFile[is] = 0;
  }
  //
  AliInfoF("Summary: selected %d tracks (%d with outliers) | %.1f MB read in %d read calls",
	   fNTrSelTot,fNTrSelTotWO,float(fNBytesReadTot)/kMByte,fNReadCallTot); 
  swTot.Print();

  AliSysInfo::AddStamp("ProjTreeLocSave");

  if (mode==kExtractMode) WriteStatHistos();
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
    if (!FindVoxelBin(sectID, fArrTgSlp[icl], fQ2Pt, fArrX[icl], fArrYCl[icl], fArrZCl[icl], fDTS.bvox, voxVars)) continue;    
    fDTS.dy   = fArrDY[icl];
    fDTS.dz   = fArrDZ[icl];
    //    fDTS.dy   = (fArrDY[icl]+fMaxDY)*fDeltaYbinI;
    //    fDTS.dz   = (fArrDZ[icl]+fMaxDZ)*fDeltaZbinI;
    //
    fTmpTree[sectID]->Fill();
    //
    // fill statistics on distribution within the voxel, last dimension, kVoxV is for Nentries
    ULong64_t binToFill = GetBin2Fill(fNBProdSt,fDTS.bvox,kVoxV); // bin of sector stat histo
    float &binEntries = fArrNDStat[sectID]->At(binToFill); // entries in the voxel
    float oldEntries  = binEntries++;
    float norm        = 1.f/binEntries;
    for (int iv=kVoxDim;iv--;) {
      float &mean = fArrNDStat[sectID]->At(binToFill+iv-kVoxV);
      mean = ( mean*oldEntries + voxVars[iv]) * norm; // account new bin entry in averages calculation
    }
    //
  } // loop over clusters
}

//________________________________________________
void AliTPCDcalibRes::FillCorrectedResiduals()
{
  // fill local trees result of closure test: corrected distortions
  
  float voxVars[kVoxHDim]={0}; // voxel variables (unbinned)
  for (int icl=fNCl;icl--;) {
    if (fArrX[icl]<kInvalidR) continue; // rejected outlier
    int sectID = fArrSectID[icl]; // 0-35 numbering
    // 
    // extract correction
    // calculate voxel variables and bins
    if (!FindVoxelBin(sectID, fArrTgSlp[icl], fQ2Pt, fArrX[icl], fArrYCl[icl], fArrZCl[icl], fDTC.bvox, voxVars)) continue;    
    int row159 = GetRowID(fArrX[icl]);
    if (row159<0) continue;
    float corr[3];

    fChebCorr->Eval(sectID, row159, fArrYCl[icl]/fArrX[icl], fArrZCl[icl]/fArrX[icl], corr);
    // 
    fDTC.t   = fTimeStamp;
    fDTC.dyR = fArrDY[icl];
    fDTC.dzR = fArrDZ[icl];

    fDTC.dyC = fArrDY[icl] - (corr[kResY]-corr[kResX]*fArrTgSlp[icl]);
    fDTC.dzC = fArrDZ[icl] - (corr[kResZ]+corr[kResX]*fTgLam);

    fDTC.q2pt   = fQ2Pt;
    fDTC.tgLam  = fTgLam;
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
  TString namef;
  for (int is=0;is<kNSect2;is++) {
    if      (mode==kExtractMode)     namef = Form("%s%d.root",kLocalResFileName,is);
    else if (mode==kClosureTestMode) namef = Form("%s%d.root",kClosureTestFileName,is);
    else AliFatalF("unknown mode: %d",mode);
    fTmpFile[is] = TFile::Open(namef.Data(),"recreate");
    fTmpTree[is] = new TTree(Form("ts%d",is),"");
    //
    if (mode==kExtractMode) {
      fTmpTree[is]->Branch("dts",&dtsP);
      //fTmpTree[is]->SetAutoFlush(150000);
      //
      fStatHist[is] = CreateVoxelStatHisto(is);
      fArrNDStat[is] = (TNDArrayT<float>*)&fStatHist[is]->GetArray();
    }
    else if (mode==kClosureTestMode) {
      fTmpTree[is]->Branch("dtc",&dtcP);
    }
  }
}

//__________________________________________________________________________________
Bool_t AliTPCDcalibRes::CompareToHelix(float *resHelixY, float *resHelixZ)
{
  // compare track to helix, refit q/pt and tgLambda and build array of tg(slope) at pad-rows
  const double kEps = 1e-12;
  float xlab[kNPadRows],ylab[kNPadRows],spath[kNPadRows];
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
    xlab[ip] = fArrX[ip]*cs - fArrDY[ip]*sn;
    ylab[ip] = fArrDY[ip]*cs + fArrX[ip]*sn;
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
    double xRow = fArrX[ip]*cs; 
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
  CollectData(kClosureTestMode);
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
  bstat_t voxStat, *statP = &voxStat;
  AliSysInfo::AddStamp("ProcResid",0,0,0,0);
  TFile* flOut = new TFile(Form("%sTree.root",kStatOut),"recreate");
  fStatTree = new TTree("voxStat","");
  fStatTree->Branch("bins",&statP);
  //
  for (int i=0;i<kNEstPar;i++) {
    fStatTree->SetAlias(Form("Y%s",kEstName[i]),Form("distY[%d]",i));
    fStatTree->SetAlias(Form("Z%s",kEstName[i]),Form("distZ[%d]",i));
  }
  for (int i=0;i<kVoxDim;i++) {
    fStatTree->SetAlias(kVoxName[i],Form("bvox[%d]",i));
    fStatTree->SetAlias(Form("%sAV",kVoxName[i]),Form("stat[%d]",i));
  }
  fStatTree->SetAlias("N",Form("stat[%d]",kVoxV)); // entries
  //
  for (int is=0;is<kNSect2;is++) {
    ProcessSectorResiduals(is, voxStat);
    AliSysInfo::AddStamp("ProjResid",is);
    //
    // don't delete here, they are needed for dispersions extraction
    //if (fDeleteSectorTrees) {
    //  TString sectFileName = Form("%s%d.root",kLocalResFileName,is);
    //  AliInfoF("Deleting %s",sectFileName.Data());
    //  unlink(sectFileName.Data());
    // }
  }
  //
  flOut->cd();
  fStatTree->Write("", TObject::kOverwrite);
  delete fStatTree;
  fStatTree = 0;
  flOut->Close();
  delete flOut;
  //
  AliSysInfo::AddStamp("ProcResid",1,0,0,0);
  //
}

//________________________________________________
void AliTPCDcalibRes::ProcessDispersions()
{
  // extract distortions of corrected Y residuals
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
  // extract dispersion of corrected residuals
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
  //
  THnS *hisCY = (THnS*)CreateSectorResidualsHisto(is, fNDeltaYBins, fMaxDY, kResName[kResD]);
  TNDArrayT<short>* ndArrYC = (TNDArrayT<short>*)&hisCY->GetArray();
  bres_t* sectData = fSectGVoxRes[is]; // raw and smoothed results
  int npAcc = 0;
  for (int ip=0;ip<npoints;ip++) {
    sectTree->GetEntry(ip);
    float dy = fDTS.dy; 
    //   float dz = fDTS.dz; // we use only dy to extract dispersion
    //
    // extract smoothed distortions for this voxel
    Long64_t binGlo = GetVoxGBin(fDTS.bvox[kVoxX],fDTS.bvox[kVoxF],fDTS.bvox[kVoxZ]);
    bres_t* voxRes = &sectData[binGlo];
    //
    // extract mean q/pt for this voxel
    binGlo = GetBin2Fill(fNBProdSt,fDTS.bvox,kVoxQ);
    float tgSlp = fArrNDStat[is]->At(binGlo);
    //
    // correct 
    dy -= voxRes->DS[kResY] - voxRes->DS[kResX]*tgSlp;
    //
    // don't allow for over/underflow
    if (TMath::Abs(dy)>fMaxDY-kEps) continue;
    //
    UShort_t  binDY = (dy+fMaxDY)*fDeltaYbinI;
    binGlo = GetBin2Fill(fNBProdDY,fDTS.bvox,binDY);
    ndArrYC->At(binGlo)++;
    npAcc++;
    //
  }
  //
  hisCY->SetEntries(npAcc);
  //
  ExtractVoxelDispersion(is, ndArrYC);
  //
  // at the moment save corrected histos ...
  TFile* flOut = TFile::Open(Form("residualSect%d.root",is),"update"); // RS shall we dump all sector histos to 1 file?
  hisCY->Write("", TObject::kOverwrite);
  flOut->Close(); 
  //
  delete sectTree;
  sectFile->Close();
  delete sectFile;
  //
  AliInfoF("Sector %2d | timing: real: %.3f cpu: %.3f",is, sw.RealTime(), sw.CpuTime());
  AliSysInfo::AddStamp("ProcessSectorDispersions",1,0,0,0);
  //
}

//_________________________________
void  AliTPCDcalibRes::WriteVoxelDefinitions()
{
  // Store voxel boundaries
  if (!fInitDone) {AliError("Init not done"); return;}

  voxDef_t vdef, *vdefP = &vdef;
  //
  TFile* flOut = new TFile(Form("%sTree.root",kStatOut),"update");
  //
  TTree* trDef = new TTree("voxDef","Voxel Boundaries definition");
  trDef->Branch("vDef",&vdef);
  for (int ix=0;ix<fNXBins;ix++) {
    vdef.bvox[kVoxX] = ix;
    vdef.vmin[kVoxX] = GetXLow(ix);
    vdef.vmax[kVoxX] = vdef.vmin[kVoxX] + GetDX(ix);
    //
    for (int ip=0;ip<fNY2XBins;ip++) {
      vdef.bvox[kVoxF] = ip;
      vdef.vmin[kVoxF] = GetY2XLow(ix,ip);
      vdef.vmax[kVoxF] = vdef.vmin[kVoxF] + GetDY2X(ix);
      //
      for (int iq=0;iq<kNQBins;iq++) {
	vdef.bvox[kVoxQ] = iq;
	int idxy = ix*fNY2XBins + ip;
	float xc = 0.5*(vdef.vmin[kVoxX]+vdef.vmax[kVoxX]);
	float yc = 0.5*(vdef.vmin[kVoxF]+vdef.vmax[kVoxF])*xc;
	float tgp0 = tgpXY(xc,yc, 0.5*(vdef.vmin[kVoxX]+vdef.vmax[kVoxX],fQ2PTBound[iq]),fBz);
	float tgp1 = tgpXY(xc,yc, 0.5*(vdef.vmin[kVoxX]+vdef.vmax[kVoxX],fQ2PTBound[iq+1]),fBz);
	vdef.vmin[kVoxQ] = tgp0<tgp1 ? tgp0:tgp1;
	vdef.vmax[kVoxQ] = tgp0<tgp1 ? tgp1:tgp0;
	/*
	if (fBinDQI[idxy]>0) {
	  vdef.vmin[kVoxQ] = fBinMinQ[idxy] + iq*fBinDQ[idxy];
	  vdef.vmax[kVoxQ] = vdef.vmin[kVoxQ] + fBinDQ[idxy];
	}
	else {
	  vdef.vmin[kVoxQ] = -1;
	  vdef.vmax[kVoxQ] =  1;	  
	}
	*/	
	//
	for (int iz=0;iz<fNZ2XBins;iz++) {
	  vdef.bvox[kVoxZ] = iz;
	  vdef.vmin[kVoxZ] = GetZ2XLow(iz);
	  vdef.vmax[kVoxZ] = vdef.vmin[kVoxZ] + GetDZ2X();
	  //
	  trDef->Fill();
	}
      }
    }
  }
  //
  for (int i=0;i<kVoxDim;i++) {
    trDef->SetAlias(kVoxName[i],Form("bvox[%d]",i));
    trDef->SetAlias(Form("%sMN",kVoxName[i]),Form("vmin[%d]",i));
    trDef->SetAlias(Form("%sMX",kVoxName[i]),Form("vmax[%d]",i));
  }
  //
  trDef->Write("", TObject::kOverwrite);
  delete trDef;
  //
  flOut->Close();
  delete flOut;
  //
  AliSysInfo::AddStamp("WriteVoxDef",0,0,0,0);
  //
}


//_________________________________________________
void AliTPCDcalibRes::ProcessSectorResiduals(int is, bstat_t &voxStat)
{
  // process residuals for single sector and store in the tree
  //
  TStopwatch sw;  sw.Start();
  AliSysInfo::AddStamp("ProcSectRes",is,0,0,0);
  //
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
  //
  THnS *hisY = (THnS*)CreateSectorResidualsHisto(is, fNDeltaYBins, fMaxDY, kResName[kResY]);
  THnS *hisZ = (THnS*)CreateSectorResidualsHisto(is, fNDeltaZBins, fMaxDZ, kResName[kResZ]);
  TNDArrayT<short>* ndArrY = (TNDArrayT<short>*)&hisY->GetArray();
  TNDArrayT<short>* ndArrZ = (TNDArrayT<short>*)&hisZ->GetArray();
  //
  for (int ip=0;ip<npoints;ip++) {
    sectTree->GetEntry(ip);
    UShort_t  binDY = (fDTS.dy+fMaxDY)*fDeltaYbinI;
    UShort_t  binDZ = (fDTS.dz+fMaxDZ)*fDeltaZbinI;
    //    UShort_t  binDY = fDTS.dy;
    //    UShort_t  binDZ = fDTS.dz;
    ULong64_t binToFillY = GetBin2Fill(fNBProdDY,fDTS.bvox,binDY);
    ULong64_t binToFillZ = GetBin2Fill(fNBProdDZ,fDTS.bvox,binDZ);
    //
    ndArrY->At(binToFillY)++;
    ndArrZ->At(binToFillZ)++;
    //
  }
  //
  hisY->SetEntries(npoints);
  hisZ->SetEntries(npoints);
  //
  delete sectTree;
  sectFile->Close(); // to reconsider: reuse the file
  delete sectFile;
  //
  AliSysInfo::AddStamp("ProjSectRes",is,0,0,0);
  //
  // extract full info about the voxel and write in the tree
  voxStat.bsec = is;
  ExtractVoxelData(voxStat, ndArrY, ndArrZ, fArrNDStat[is]);

  // at the moment save histos ...
  TFile* flOut = TFile::Open(Form("residualSect%d.root",is),"update"); // RS shall we dump all sector histos to 1 file?
  hisY->Write("", TObject::kOverwrite);
  hisZ->Write("", TObject::kOverwrite);
  flOut->Close();
  //
  delete hisY;
  delete hisZ;
  //
  sw.Stop(); 
  AliInfoF("Sector %2d | timing: real: %.3f cpu: %.3f",is, sw.RealTime(), sw.CpuTime());
  AliSysInfo::AddStamp("ProcSectRes",is,1,0,0);
  //
}

//_________________________________________________
void AliTPCDcalibRes::ExtractVoxelData(bstat_t &stat, 
				       const TNDArrayT<short>* harrY, 
				       const TNDArrayT<short>* harrZ,
				       const TNDArrayT<float>* harrStat)
{
  // Extract distortion estimators from each voxel histo
  if (!fHDelY) {fHDelY = new TH1F("dy","dy",fNDeltaYBins,-fMaxDY,fMaxDY); fHDelY->SetDirectory(0);}
  if (!fHDelZ) {fHDelZ = new TH1F("dz","dz",fNDeltaZBins,-fMaxDZ,fMaxDZ); fHDelZ->SetDirectory(0);}
  //
  UChar_t bvox1[kVoxDim];
  for (stat.bvox[kVoxZ]=0;stat.bvox[kVoxZ]<fNZ2XBins;stat.bvox[kVoxZ]++) {
    for (stat.bvox[kVoxX]=0;stat.bvox[kVoxX]<fNXBins;stat.bvox[kVoxX]++) { 
      for (stat.bvox[kVoxF]=0;stat.bvox[kVoxF]<fNY2XBins;stat.bvox[kVoxF]++) {
	//
	// in Z we integrate over all Q bins, since no dependence is expected
	stat.bvox[kVoxQ]=0;
	for (int i=kVoxDim;i--;) bvox1[i] = stat.bvox[i];
	bvox1[kVoxQ] = kNQBins-1;
	ExtractResidualHisto(harrZ,fNBProdDZ,stat.bvox,bvox1,fHDelZ); // integrate over Q bins
	ExtractDistortionsData(fHDelZ, stat.distZ, stat.bvox);
	//	
	for (stat.bvox[kVoxQ]=0;stat.bvox[kVoxQ]<kNQBins;stat.bvox[kVoxQ]++) {
	  //
	  ExtractResidualHisto(harrY,fNBProdDY,stat.bvox,fHDelY);
	  ExtractDistortionsData(fHDelY, stat.distY,stat.bvox);
	  //
	  // extract voxel statistics: COG for each dimension
	  Long64_t bglo = GetBin2Fill(fNBProdSt,stat.bvox,kVoxV);
	  for (int i=0;i<kVoxHDim;i++) stat.stat[i] = harrStat->At(bglo+i-kVoxV);
	  fStatTree->Fill();
	}
      }
    }
  }
  //
}

//_________________________________________________
void AliTPCDcalibRes::ExtractVoxelDispersion(int is, const TNDArrayT<short>* harrYC,
					     float maxGChi2)
{
  // Extract distortion estimators from each voxel histo
  const float kMaxSigG2L = 5.0f; // max ratio between Gaussian sigma and RMS LTM
  const float kMinNormG2M = 0.2f; // min ratio between Gaussian amplitude and max Value
  const float kZeroSigma = 1e-4; 
  //
  if (!fHDelY) {fHDelY = new TH1F("dy","dy",fNDeltaYBins,-fMaxDY,fMaxDY); fHDelY->SetDirectory(0);}
  //
  bres_t* sectData = fSectGVoxRes[is]; // raw and smoothed results
  float estDisp[kNEstPar] = {0.0f};
  UChar_t bvox[kVoxDim],bvox1[kVoxDim];
  // 1st extract raw dispersion
  for (bvox[kVoxZ]=0;bvox[kVoxZ]<fNZ2XBins;bvox[kVoxZ]++) {
    for (bvox[kVoxX]=0;bvox[kVoxX]<fNXBins;bvox[kVoxX]++) { 
      for (bvox[kVoxF]=0;bvox[kVoxF]<fNY2XBins;bvox[kVoxF]++) {
	//
	// we integrate over all Q bins, since no dependence is expected
	bvox[kVoxQ]=0;
	for (int i=kVoxDim;i--;) bvox1[i] = bvox[i];
	bvox1[kVoxQ] = kNQBins-1;
	ExtractResidualHisto(harrYC,fNBProdDY,bvox,bvox1,fHDelY); // integrate over Q bins
	ExtractDistortionsData(fHDelY, estDisp , bvox);
	//
	Long64_t binGlo = GetVoxGBin(bvox[kVoxX],bvox[kVoxF],bvox[kVoxZ]);
	bres_t* voxRes = &sectData[binGlo]; // destinate voxel to store the dispersion
	//
	Bool_t okG=kFALSE,okL;
	okG = okL = estDisp[kEstSigL]>kZeroSigma && estDisp[kEstMeanEL]>0;
	voxRes->D[kResD] = 0;
	voxRes->E[kResD] = -1;
	//
	if (okG && 
	    estDisp[kEstNormG]<kMinNormG2M*estDisp[kEstMax] || // gaussian norm should not be negligible
	    estDisp[kEstChi2G]>maxGChi2 || estDisp[kEstSigG]<kZeroSigma ||
	    estDisp[kEstSigG]>kMaxSigG2L) okG = kFALSE;

	if (okL) { 
	  voxRes->D[kResD] = okG ? estDisp[kEstSigG]  : estDisp[kEstSigL];
	  voxRes->E[kResD] = okG ? estDisp[kEstSigEG]  : estDisp[kEstSigEL];
	  if (!voxRes->E[kResD]) voxRes->E[kResD] = voxRes->D[kResD]/TMath::Sqrt(estDisp[kEstNorm]); // LTM does not provide the error
	}
      } // loop over voxF
    } // loop over voxX
  } // loop over voxZ
  //
  // now smooth the dispersion
  for (bvox[kVoxZ]=0;bvox[kVoxZ]<fNZ2XBins;bvox[kVoxZ]++) {
    for (bvox[kVoxX]=0;bvox[kVoxX]<fNXBins;bvox[kVoxX]++) { 
      for (bvox[kVoxF]=0;bvox[kVoxF]<fNY2XBins;bvox[kVoxF]++) {
	int binGlo = GetVoxGBin(bvox[kVoxX],bvox[kVoxF],bvox[kVoxZ]);
	bres_t *voxRes = &sectData[binGlo];
	Bool_t res = GetSmoothEstimateDim(is,voxRes->stat[kVoxX],voxRes->stat[kVoxF],voxRes->stat[kVoxZ],
				       int(kResD), voxRes->DS[kResD]);
      }
    }
  }
}

//______________________________________________________________________________
void AliTPCDcalibRes::ExtractDistortionsData(TH1F* histo, float est[kNEstPar], const UChar_t vox[kVoxDim], float minNorm, float fracLTM)
{
  const float kMinEntries=30;
  static TF1 fgaus("fgaus","gaus",-10,10);
  float nrm=0,mean=0,mom2=0,rms=0,maxVal=0;
  //
  memset(est,0,kNEstPar*sizeof(float));
  //
  float *w = (float*) histo->GetArray();
  w++; // skip underflows
  int nb = histo->GetNbinsX();
  float x = histo->GetXaxis()->GetXmax();
  float dx = (x+x)/nb;
  float dxh = 0.5*dx;
  x -= dxh;
  for (int ip=nb;ip--;) {
    nrm  += w[ip];
    mean += x*w[ip];
    mom2 += x*x*w[ip];
    x -= dx;
    if (maxVal<w[ip]) maxVal = w[ip];
  }
  if (nrm>0) {
    mean /= nrm;
    mom2 /= nrm;
    rms = mom2 - mean*mean;
    rms = rms>0 ? TMath::Sqrt(rms):0;
  }
  est[kEstNorm] = nrm;
  est[kEstMean] = mean;
  est[kEstSig]  = rms;
  est[kEstMax]  = maxVal;
  if (nrm<minNorm) return;
  //
  float bwsig = histo->GetBinWidth(1)/TMath::Sqrt(12);
  const int kNLTMTests = 11;
  const float kLTMTests[kNLTMTests]={1.00,0.95,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.50};
  double ltmMuEst[kNLTMTests], ltmSigEst[kNLTMTests];
  double logLArr[kNLTMTests],logL0Arr[kNLTMTests],logLDifArr[kNLTMTests];
  double sigEstArr[kNLTMTests],muEstArr[kNLTMTests];
  int bminArr[kNLTMTests],bmaxArr[kNLTMTests];
  TVectorF vecLTM(10);
  //
  // store reference LTM
  TStatToolkit::LTMHisto(histo, vecLTM, fracLTM); 
  if (!vecLTM[3]) return;
  est[kEstMeanL]  = vecLTM[1];
  est[kEstSigL]   = TMath::Max(vecLTM[2],bwsig);
  est[kEstMeanEL] = vecLTM[3];
  est[kEstSigEL]  = vecLTM[4];
  //
  est[kEstChi2G] = 999.0f;
  //
  if (nrm<kMinEntries) return; // don't do fit if there are no enough entries
  //
  int nltmAcc = 0;
  int bminPrev = -1, bmaxPrev = -1;
  for (int iltm=0;iltm<kNLTMTests;iltm++) {
    TStatToolkit::LTMHisto(histo, vecLTM, kLTMTests[iltm]);
    int bmin = int(vecLTM[5]), bmax = int(vecLTM[6]);
    if (bmin==bminPrev && bmax==bmaxPrev) continue; // same range
    bminPrev = bminArr[nltmAcc] = bmin;
    bmaxPrev = bmaxArr[nltmAcc] = bmax;
    if (bmax-bmin<2) continue; // don't use too narow window
    // skip empty bins from edges
    while (!histo->GetBinContent(bmin)) bmin++;
    while (!histo->GetBinContent(bmax)) bmax--;
    double muEst = ltmMuEst[nltmAcc]  = vecLTM[1];
    double sigEst = ltmSigEst[nltmAcc] = TMath::Max(vecLTM[2],bwsig);
    if (sigEst<bwsig) sigEst = bwsig;
    //
    // extract non-truncated estimators and sample and reference log-likelihoods
    logLArr[nltmAcc] = GetLogL(histo,bmin,bmax,muEst,sigEst,logL0Arr[nltmAcc]);
    if (logLArr[nltmAcc]<-1e8) {
      AliWarningF("Failure for LTM_%.3f in voxel Q:%d F:%d X:%d Z:%d",
		  kLTMTests[iltm],vox[kVoxQ],vox[kVoxF],vox[kVoxX],vox[kVoxZ]);
    }   
    sigEstArr[nltmAcc] = sigEst;
    muEstArr[nltmAcc]  = muEst;
    //
    nltmAcc++;
  }
  // select best cut and fit
  const int kNCutsLL=4;
  const float kCutsLL[kNCutsLL] = {3.,5.,8.,10.};
  int ind = -1;
  for (int ict=0;ict<kNCutsLL;ict++) {
    for (int ift=1;ift<nltmAcc-1;ift++) {
      double logLRat = logL0Arr[ift] - logLArr[ift];
      double logLRatNxt = logL0Arr[ift+1] - logLArr[ift+1];
      if (logLRat<kCutsLL[ict] && logLRat>logLRatNxt) {
	float rngMin = histo->GetBinCenter(bminArr[ift]) - dxh;
	float rngMax = histo->GetBinCenter(bmaxArr[ift]) + dxh;

	fgaus.SetParameters(maxVal ,muEstArr[ift], sigEstArr[ift]);
	TFitResultPtr fitPtr = histo->Fit(&fgaus,"qnrLS","",rngMin,rngMax);	
	TFitResult * result = fitPtr.Get();
	float estMG=0,estSG=0,estMGE=0,estSGE=0,chi2=0;
	if (result!=NULL) {
	  est[kEstMeanG] = fgaus.GetParameter(1);
	  if (TMath::Abs(est[kEstMeanG]-ltmMuEst[ift])>ltmSigEst[ift]) continue;
	  est[kEstNormG] = fgaus.GetParameter(0);
	  est[kEstSigG]  = fgaus.GetParameter(2);
	  est[kEstMeanEG] = fgaus.GetParError(1);
	  est[kEstSigEG]  = fgaus.GetParError(2);
	  int npf = fgaus.GetNumberFreeParameters();
	  est[kEstChi2G] = npf>0 ? fgaus.GetChisquare()/npf : 0;
	  //
	  ind = ift;
	  break;
	}	
	//
      }
    }
    if (ind>-1) break;
  }
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
    AliWarningF("Too small sigma %.4e is provided for bin width %.4e",sig,dxh);
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
    AliWarningF("Input sigma %e is too small",sig);
    return kFALSE;
  }
  if ( (b-a)/sig<kMinWindow ) {
    AliWarningF("Truncation window %e-%e is %e sigma only",a,b,(b-a)/sig);
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

  voxAxisName[kVoxQ] = "tgphi_Bin";
  voxNBins[kVoxQ]    = kNQBins;
  voxBinMin[kVoxQ]   = 0;
  voxBinMax[kVoxQ]   = kNQBins;
  //
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

//___________________________________________________________________________
THn* AliTPCDcalibRes::CreateSectorResidualsHisto(int sect, int nbDelta,float range, const char* pref)
{
  // prepare histogram to store the residuals within the sector

  // create binning for voxels and residuals
  Int_t voxNBins[kVoxHDim];
  Double_t voxBinMin[kVoxHDim],voxBinMax[kVoxHDim];
  TString  voxAxisName[kVoxHDim];

  voxAxisName[kVoxQ] = "tgphi_Bin";
  voxNBins[kVoxQ]    = kNQBins;
  voxBinMin[kVoxQ]   = 0;
  voxBinMax[kVoxQ]   = kNQBins;
  //
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
  //
  voxAxisName[kVoxV]   = pref;
  voxNBins[kVoxV]      = nbDelta;
  voxBinMin[kVoxV]     =-range;
  voxBinMax[kVoxV]     = range;
  //
  THnS* h = new THnS(Form("delta%d_%s",sect,pref),"",kVoxHDim,voxNBins,voxBinMin,voxBinMax);
  for (int i=0;i<kVoxHDim;i++) {
    h->GetAxis(i)->SetName(voxAxisName[i].Data());
    h->GetAxis(i)->SetTitle(voxAxisName[i].Data());
  }
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
  for (int i=fNCl;i--;) if (rejCl[i]) fArrX[i] = -1;

  return kTRUE;
}

//______________________________________________
void AliTPCDcalibRes::FixAlignmentBug(int sect, float q2pt, float bz, float& alp, 
				      float& x, float &z, float &deltaY, float &deltaZ)
{
  // fix alignment bug: https://alice.its.cern.ch/jira/browse/ATO-339?focusedCommentId=170850&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-170850
  //
  // alp, x, z correspond to 
  //
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

  // cluster in its proper alpha frame with alignment bug, Z trackITS is used !!! 
  double xyzClUse[3] = {x,0,z}; // this is what we read from the residual tree, ITS Z only is stored
  double xyzTrUse[3] = {x, deltaY, z}; // track in bad cluster frame
  //
  // recover cluster Z position by adding deltaZ, this is approximate, since ITS track Z was used...
  xyzClUse[2] -= deltaZ;
  static AliExternalTrackParam trDummy;
  trDummy.Local2GlobalPosition(xyzClUse,alp); // misaligned cluster in global frame
  double xyz0[3]={xyzClUse[0],xyzClUse[1],xyzClUse[2]};
  mgt->MasterToLocal(xyz0,xyzClUse);
  // we got ideal cluster in the sector tracking frame, 
  //
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
int AliTPCDcalibRes::CheckResiduals(Bool_t* kill, float &rmsLongMA)
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

  memset(kill,0,fNCl*sizeof(Bool_t));
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
  if (naccY<nMinAcc || naccZ<nMinAcc) { // kill all
    for (int i=fNCl;i--;) kill[i] = kTRUE;
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
    for (int i=fNCl;i--;) kill[i] = kTRUE;
    return fNCl;
  }
  //
  //  printf("RMSY %d min of %d: %f | RMSZ %d min of %d: %f\n",kmnY,naccY,rmsKY, kmnZ,naccZ,rmsKZ);
  //
  //
  float rmsKYI = 1./rmsKY;
  float rmsKZI = 1./rmsKZ;
  int nKill=0, nacc = 0;
  float yacc[kNPadRows],yDiffLong[kNPadRows];
  for (int ip=0;ip<fNCl;ip++) {

    yDiffLL[ip] *= rmsKYI;
    zDiffLL[ip] *= rmsKZI;
    float dy = yDiffLL[ip], dz = zDiffLL[ip];
    if (dy*dy+dz*dz>fMaxStdDevMA) {
      kill[ip] = kTRUE;
      nKill++;
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
  return nKill;
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

//___________________________________________________________________
void AliTPCDcalibRes::medFit(int np, const float* x, const float* y, float &a, float &b, float delI)
{
  // Median linear fit: minimizes abs residuals instead of squared ones
  // Adapted from "Numerical Recipes in C"
  float aa,bb,b1,b2,f,f1,f2,sigb,chisq=0.0f;
  if (!delI) {
    float sx=0.0f,sxx=0.0f,sy=0.0f,sxy=0.0f,del;
    //
    for (int j=np;j--;) { sx += x[j]; sxx += x[j]*x[j];}
    del = np*sxx-sx*sx;
    //
    for (int j=np;j--;) { sy += y[j]; sxy += x[j]*y[j];}
    //
    float delI = 1./del;
    aa = (sxx*sy-sx*sxy)*delI;
    bb = (np*sxy-sx*sy)*delI;
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
  TFile *fdrift = TFile::Open(Form("%s.root",kDriftFileName));
  if (fdrift) {
    TTree * tree = (TTree*)fdrift->Get("fitTimeStat");
    if (tree==NULL) {
      ::Error("LoadDriftCalibration FAILED", "tree fitTimeStat not avaliable in file %s.root",kDriftFileName);
    }
    else {      
      tree->SetBranchAddress("grTRDReg.",&fVDriftGraph);
      tree->SetBranchAddress("paramRobust.",&fVDriftParam);
      tree->GetEntry(0);
      if (fVDriftGraph==NULL || fVDriftGraph->GetN()<=0) {
	AliInfo("ITS/TRD drift calibration not availalble. Trying ITS/TOF");
	tree->SetBranchAddress("grTOFReg.",&fVDriftGraph);
	tree->GetEntry(0);
      }
      /*
      else {
	::Info("LoadDriftCalibration", "tree fitTimeStat not avaliable in file %s.root",kDriftFileName);
      }
      */
    }
    delete tree;
  }
  else {
    ::Error("LoadDriftCalibration FAILED", "fitDrift.root not present");
  }
  if (fdrift) fdrift->Close();
  delete fdrift;
}

//_________________________________________
float AliTPCDcalibRes::GetDriftCorrection(float z, float x, float phi, int rocID)
{
  // apply vdrift correction
  int side = ((rocID/kNSect)&0x1) ? -1:1; // C:A
  float drift = side>0 ? kZLim-z : z+kZLim;
  float gy    = TMath::Sin(phi)*x;
  Double_t pvecFit[3];
  pvecFit[0]= side;             // z shift (cm)
  pvecFit[1]= drift*gy/kZLim;   // global y gradient
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

  TFile* flOut = new TFile(Form("%sTree.root",kResOut),"recreate");
  TTree* resTree = new TTree("voxRes","final distortions, see GetListOfAliases");
  resTree->Branch("res",&voxRes);
  for (int i=0;i<kVoxDim;i++) {
    resTree->SetAlias(kVoxName[i],Form("bvox[%d]",i));
    resTree->SetAlias(Form("%sAV",kVoxName[i]),Form("stat[%d]",i));
  }
  for (int i=0;i<kResDim;i++) {
    resTree->SetAlias(kResName[i],Form("D[%d]",i));
    resTree->SetAlias(Form("%sE",kResName[i]),Form("E[%d]",i));
  }
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



//=========================================================================
//
//   voxels processing related methods
//
//=========================================================================

//___________________________________________________________________________________
float AliTPCDcalibRes::ExtractResidualHisto(const TNDArrayT<short>* harr, const Long64_t bprod[kVoxHDim], 
						       const UChar_t vox[kVoxDim], TH1F* dest)
{
  // extract residuals for the voxel from THn to 1D histo
  // Here we use shortcut to fill the histo, don't modify neither this part, not the histo def.
  float* arr = (float*)dest->GetArray();
  int nb = dest->GetNbinsX();
  arr++; // skip the underflows bin
  ULong64_t binGlo = GetBin2Fill(bprod,vox,0);
  float sum = 0;
  for (int i=0;i<nb;i++) sum += arr[i] = harr->At(binGlo++);
  dest->SetEntries(sum);
  return sum;
}

//___________________________________________________________________________________
float AliTPCDcalibRes::ExtractResidualHisto(const TNDArrayT<short>* harr, const Long64_t bprod[kVoxHDim], 
						       const UChar_t voxMin[kVoxDim], const UChar_t voxMax[kVoxDim], TH1F* dest)
{
  // extract residuals for the voxel from THn to 1D histo
  // Here we use shortcut to fill the histo, don't modify neither this part, not the histo def.
  float* arr = (float*)dest->GetArray();
  int nb = dest->GetNbinsX();
  memset(arr,0,sizeof(float)*(nb+2));
  arr++; // skip the underflows bin
  float sum = 0;
  UChar_t vox[kVoxDim];
  
  for (vox[kVoxQ]=voxMin[kVoxQ];vox[kVoxQ]<=voxMax[kVoxQ];vox[kVoxQ]++) 
    for (vox[kVoxF]=voxMin[kVoxF];vox[kVoxF]<=voxMax[kVoxF];vox[kVoxF]++) 
      for (vox[kVoxX]=voxMin[kVoxX];vox[kVoxX]<=voxMax[kVoxX];vox[kVoxX]++) 
	for (vox[kVoxZ]=voxMin[kVoxZ];vox[kVoxZ]<=voxMax[kVoxZ];vox[kVoxZ]++) {
	  ULong64_t binGlo = GetBin2Fill(bprod,vox,0);
	  for (int i=0;i<nb;i++) arr[i] += harr->At(binGlo++);
	}
  //
  for (int i=0;i<nb;i++) sum += arr[i];
  dest->SetEntries(sum);
  return sum;
}

//___________________________________________________________________________________
TH1F* AliTPCDcalibRes::ExtractResidualHisto(int htype, int sect, const UChar_t vox[kVoxDim])
{
  // create Y or Z residuals for the voxel in sector sect,  from THn to 1D histo
  TString fln = Form("residualSect%d.root",sect);
  TFile* fl = TFile::Open(fln.Data());
  if (!fl) {printf("Failed to open %s\n",fln.Data()); return 0;}
  TString hname = Form("delta%d_%s",sect,kResName[htype]);
  THnS* hn = (THnS*)fl->Get(hname.Data());
  if (!hn) {printf("Did not find histo %s in %s\n",hname.Data(),fln.Data());return 0;}
  TString ttl = Form("d%s_S%d_vox%d_%d_%d_%d",kResName[htype] ,sect,vox[kVoxQ],vox[kVoxF],vox[kVoxX],vox[kVoxZ]);
  Bool_t typeY = htype==kResY||htype==kResD;
  TH1F* h1 = new TH1F(ttl.Data(),ttl.Data(),
		      typeY ? fNDeltaYBins : fNDeltaZBins ,
		      typeY ? -fMaxDY : -fMaxDZ,
		      typeY ?  fMaxDY :  fMaxDZ);
  h1->SetDirectory(0);
  const TNDArrayT<short>* harr = (TNDArrayT<short>*)&hn->GetArray();
  ExtractResidualHisto(harr, typeY ? fNBProdDY : fNBProdDZ, vox, h1);
  delete hn;
  fl->Close();
  delete fl;
  printf("Don't forget to delete (TH1F*)%p %s\n",h1,h1->GetName());
  return h1;
}

//___________________________________________________________________________________
TH1F* AliTPCDcalibRes::ExtractResidualHisto(Int_t htype, int sect, 
					    const UChar_t vox[kVoxDim], const UChar_t vox1[kVoxDim])
{
  // create Y or Z residuals for the voxel in sector sect,  from THn to 1D histo
  TString fln = Form("residualSect%d.root",sect);
  TFile* fl = TFile::Open(fln.Data());
  if (!fl) {printf("Failed to open %s\n",fln.Data()); return 0;}
  TString hname = Form("delta%d_%s",sect,kResName[htype]);
  THnS* hn = (THnS*)fl->Get(hname.Data());
  if (!hn) {printf("Did not find histo %s in %s\n",hname.Data(),fln.Data());return 0;}
  Bool_t typeY = htype==kResY||htype==kResD;
  TString ttl = Form("d%s_S%d_vox%d_%d_%d_%d__%d_%d_%d_%d", kResName[htype], sect,
		     vox[kVoxQ],vox[kVoxF],vox[kVoxX],vox[kVoxZ],
		     vox1[kVoxQ],vox1[kVoxF],vox1[kVoxX],vox1[kVoxZ]);
  TH1F* h1 = new TH1F(ttl.Data(),ttl.Data(),
		      typeY ? fNDeltaYBins : fNDeltaZBins ,
		      typeY ? -fMaxDY : -fMaxDZ,
		      typeY ?  fMaxDY :  fMaxDZ);
  h1->SetDirectory(0);
  const TNDArrayT<short>* harr = (TNDArrayT<short>*)&hn->GetArray();
  ExtractResidualHisto(harr, typeY ? fNBProdDY : fNBProdDZ, vox, vox1, h1);
  delete hn;
  fl->Close();
  delete fl;
  printf("Don't forget to delete (TH1F*)%p %s\n",h1,h1->GetName());
  return h1;
}

//__________________________________________________________________
void AliTPCDcalibRes::ExtractXYZDistortions()
{
  if (!fInitDone) Init(); //{AliError("Init not done"); return;}
  TStopwatch sw;
  LoadStatHistos();
  // extract XYZ distortions from fitted Y,Z residuals vs Q variable
  bstat_t voxStat, *statP = &voxStat;
  bstat_t voxIQ[kNQBins];
  bres_t voxRes, *voxResP=&voxRes;
  //
  AliSysInfo::AddStamp("ExtXYZ",0,0,0,0);
  TFile* flStat = 0;
  if (!fStatTree) {
    TString fname = Form("%sTree.root",kStatOut);
    flStat = new TFile(fname.Data());
    if (!flStat) AliFatalF("file %s not found",fname.Data());
    fStatTree = (TTree*)flStat->Get("voxStat");
    if (!fStatTree) AliFatalF("voxStat tree not found in %s",fname.Data());
  }
  fStatTree->SetBranchAddress("bins",&statP);  
  //
  int ent = 0;

  // 1st loop over good voxels, extract X distortions
  for (int is=0;is<kNSect2;is++) { 
    
    bres_t* sectData = fSectGVoxRes[is] = new bres_t[fNGVoxPerSector]; // here we keep main result
    voxRes.bsec = is;
    float cntGood = 0;
    for (voxRes.bvox[kVoxZ]=0;voxRes.bvox[kVoxZ]<fNZ2XBins;voxRes.bvox[kVoxZ]++) {
      for (voxRes.bvox[kVoxX]=0;voxRes.bvox[kVoxX]<fNXBins;voxRes.bvox[kVoxX]++) { 
	for (voxRes.bvox[kVoxF]=0;voxRes.bvox[kVoxF]<fNY2XBins;voxRes.bvox[kVoxF]++) {
	  for (int iq=0;iq<kNQBins;iq++) {
	    fStatTree->GetEntry(ent++);
	    //
	    // check
	    if (voxStat.bvox[kVoxZ]!=voxRes.bvox[kVoxZ] || voxStat.bvox[kVoxX]!=voxRes.bvox[kVoxX] ||
		voxStat.bvox[kVoxF]!=voxRes.bvox[kVoxF] || voxStat.bvox[kVoxQ]!=iq) {
	      AliErrorF("voxel QFXZ : Expected %d %2d %3d %2d | Read %d %2d %3d %2d",iq,
		     voxRes.bvox[kVoxF],voxRes.bvox[kVoxX],voxRes.bvox[kVoxZ],
		     voxStat.bvox[kVoxQ],voxStat.bvox[kVoxF],voxStat.bvox[kVoxX], voxStat.bvox[kVoxZ]);
	      AliFatalF("Mismatch between expected and obtained voxel %d",ent-1);
	    }
	    memcpy(&voxIQ[iq],&voxStat,sizeof(bstat_t)); // save for the analysis vs Q
	  }
	  ExtractVoxelXYZDistortions(voxIQ,voxRes);      // extract residuals deconvoluted for the slopes etc
	  cntGood += voxRes.bvox[kVoxQ]>0;
	  int binGlo = GetVoxGBin(voxRes.bvox);
	  memcpy(&sectData[binGlo],&voxRes,sizeof(bres_t)); // store in the sector data array
	  //
	}
      }
    } 
    //
    int cntSmooth = Smooth0(is); // smooth sector data
    //FillHoles(is, sectData, fNBProdSectG);

    AliInfoF("Sector%2d: voxels with data %6d (%4.1f%%) smoothed %6d (%4.1f%%) of %d",is,int(cntGood),
	   cntGood/fNGVoxPerSector*100.,cntSmooth,float(cntSmooth)/fNGVoxPerSector*100.,fNGVoxPerSector);
  }

  delete fStatTree;
  fStatTree = 0;
  flStat->Close();
  delete flStat;
  //
  sw.Stop();
  AliInfoF("timing: real: %.3f cpu: %.3f",sw.RealTime(), sw.CpuTime());
}

//_____________________________________________
Bool_t AliTPCDcalibRes::ExtractVoxelXYZDistortions(const bstat_t voxIQ[kNQBins], 
						   bres_t &res, int minStat, 
						   float maxGChi2, int minYBinsOK)
{
  // extract XYZ distortions from voxel fitted Y,Z residuals vs Q variable
  //
  const float kMaxSigG2L = 5.0f; // max ratio between Gaussian sigma and RMS LTM
  const float kMinNormG2M = 0.2f; // min ratio between Gaussian amplitude and max Value
  const float kZeroSigma = 1e-4; 
  Bool_t okG=kFALSE,okL;
  //
  int nyOK = 0;
  float av[kNQBins],meas[kNQBins],wgh[kNQBins],resFit[2],errFit[3];
  //
  for (int i=kResDim;i--;) res.D[i] = res.DS[i] = res.DC[i] = res.E[i] = 0;
  for (int i=kVoxHDim;i--;) res.stat[i] = 0;
  //
  for (int iq=kNQBins;iq--;) {
    const bstat_t &vox = voxIQ[iq];
    float ent = vox.stat[kVoxV];
    if (ent<minStat) continue;
    //
    okG = okL = vox.distY[kEstSigL]>kZeroSigma && vox.distY[kEstMeanEL]>0;
    //
    if (okG && 
	vox.distY[kEstNormG]<kMinNormG2M*vox.distY[kEstMax] || // gaussian norm should not be negligible
	vox.distY[kEstChi2G]>maxGChi2 || vox.distY[kEstSigG]<kZeroSigma ||
	vox.distY[kEstSigG]>kMaxSigG2L) okG = kFALSE;
    //
    // assume that measured Y resydual dy is related to real residuals DY and DX as
    // dy = DY - DX*tg(slope) 
    // where the slope is average track inclination angle at the pad-row
    // For DZ calculate simple weighted mean
    //
    if (okL) { // collect data for linear fit
      av[nyOK] = vox.stat[kVoxQ]; // mean value of tg(slope) for this Q bin
      meas[nyOK] = okG ? vox.distY[kEstMeanG]  : vox.distY[kEstMeanL];
      float wy   = okG ? vox.distY[kEstMeanEG] : vox.distY[kEstMeanEL];
      wgh[nyOK]  = 1./(wy*wy);
      //
      float st = vox.stat[kVoxV]; // statistics of the voxel
      res.stat[kVoxV] += st;
      for (int i=kVoxDim;i--;) res.stat[i] += vox.stat[i]*st;
      nyOK++;
    }
    //
  }
  //
  if (res.stat[kVoxV]>0) {
    float stI = 1.0f/res.stat[kVoxV];
    for (int i=kVoxDim;i--;) res.stat[i] *= stI; // average of each voxel dimension for selected bins    
  }
  if (nyOK>=minYBinsOK && FitPoly1(av,meas,wgh,nyOK,resFit,errFit)) {
    //
    res.D[kResY] = resFit[0];
    res.D[kResX] =-resFit[1];
    res.E[kResY] = errFit[0];
    res.E[kResX] = errFit[2];
    res.bvox[kVoxQ] = nyOK; // number of points used
  }
  else { // estimation impossible
    res.bvox[kVoxQ] = 0;
    // set coordinates to bin center
    GetVoxelCoordinates(res.bsec,res.bvox[kVoxX],res.bvox[kVoxF],res.bvox[kVoxZ],
			res.stat[kVoxX],res.stat[kVoxF],res.stat[kVoxZ]);
  }
  //
  // Z fits were integrated over Q, use just 1st bin (they are all the same), correcting for X shift
  // as measured DZ -> dZ + DX*<Z/X>
  const bstat_t &vox0 = voxIQ[0];
  okG = okL = vox0.distZ[kEstSigL]>kZeroSigma;  
  if (okG && 
      vox0.distZ[kEstNormG]<kMinNormG2M*vox0.distZ[kEstMax] || // gaussian norm should not be negligible
      vox0.distZ[kEstChi2G]>maxGChi2 || vox0.distZ[kEstSigG]<kZeroSigma ||
      vox0.distZ[kEstSigG]>kMaxSigG2L) okG = kFALSE;
    //
  if (okL) {
    res.D[kResZ] = (okG ? vox0.distZ[kEstMeanG]  : vox0.distZ[kEstMeanL]) + res.stat[kVoxZ]*res.D[kResX];
    res.E[kResZ] = okG ? vox0.distZ[kEstMeanEG] : vox0.distZ[kEstMeanEL];
  }
  //
  for (int i=0;i<kResDimG;i++) res.E[i] = res.E[i]>0 ? TMath::Sqrt(res.E[i]) : 0;
  //
}


//=========================================================================
//
//   fitting/smoothing related methods
//
//=========================================================================
//________________________________
void AliTPCDcalibRes::FillHoles(int isect, bres_t *sectData, const int fNBProdSectG[2], int minGoodPoints)
{
  /// RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR DELETE?
  // fill holes within 1 sector
  // scan transverse plane and extra/interpolate in Z/X 
  const float kZeroError=1e-6, kDummyError = 1.0;
  float resFit0[3],errFit0[6],resFit1[3],errFit1[6];
  float dz = GetDZ2X();
  int maxDim = TMath::Max(fNXBins,TMath::Max(fNY2XBins,fNZ2XBins));
  bres_t **currLine = new bres_t*[maxDim];
  Bool_t missY[maxDim],missZ[maxDim];
  Float_t val0[maxDim],pos0[maxDim],wgh0[maxDim];
  Float_t val1[maxDim],pos1[maxDim],wgh1[maxDim];
  AliInfoF("FillHoles Sector%d",isect);
  for (int ix=0;ix<fNXBins;ix++) {
    for (int ip=0;ip<fNY2XBins;ip++) {
      int nmissY=0,nmissZ=0;

      for (int iz=0;iz<fNZ2XBins;iz++) {  // extract line in z
	int binGlo = iz + fNBProdSectG[1]*ip + fNBProdSectG[0]*ix; // global bin
	currLine[iz] = &sectData[binGlo];
	missY[iz] = missZ[iz] = kFALSE;
	if (currLine[iz]->E[kResY]<kZeroError) {
	  missY[iz] = kTRUE;
	  nmissY++;
	}
	if (currLine[iz]->E[kResZ]<kZeroError) {
	  missZ[iz] = kTRUE;
	  nmissZ++;
	}
      }
      // 
      if (nmissY && (fNZ2XBins-nmissY)>=minGoodPoints) { // recover points in Y and X
	int npGood = 0;
	for (int iz=0;iz<fNZ2XBins;iz++) {
	  if (missY[iz]) continue;
	  val0[npGood] = currLine[iz]->D[kResY];
	  pos0[npGood] = currLine[iz]->stat[kVoxZ]; // average Z position
	  wgh0[npGood] = 1./(currLine[iz]->E[kResY]*currLine[iz]->E[kResY]);
	  // X distortion goes with Y
	  val1[npGood] = currLine[iz]->D[kResX];
	  pos1[npGood] = currLine[iz]->stat[kVoxZ]; // average Z position
	  wgh1[npGood] = 1./(currLine[iz]->E[kResX]*currLine[iz]->E[kResX]);
	  npGood++;
	}
	Bool_t res = FitPoly2(pos0,val0,wgh0,npGood, resFit0, errFit0)
	  &&         FitPoly2(pos1,val1,wgh1,npGood, resFit1, errFit1);
	if (res) {
	  for (int iz=0;iz<fNZ2XBins;iz++) {
	    if (!missY[iz]) continue;
	    // evaluate in the center of the bin
	    double z = (isect>=kNSect ? -1.0f:1.0f)*(iz+0.5)*dz, z2=z*z, z3=z2*z, z4=z3*z;
	    currLine[iz]->D[kResY] = resFit0[0]+z*resFit0[1]+z2*resFit0[2];
	    currLine[iz]->D[kResX] = resFit1[0]+z*resFit1[1]+z2*resFit1[2];
	    double evErr0= errFit0[0] + errFit0[2]*z2 + errFit0[5]*z4
	      +            2.*(errFit0[1]*z + errFit0[3]*z2 + errFit0[4]*z3);
	    double evErr1= errFit1[0] + errFit1[2]*z2 + errFit1[5]*z4
	      +            2.*(errFit1[1]*z + errFit1[3]*z2 + errFit1[4]*z3);
	    //
	    currLine[iz]->E[kResY] = evErr0>0 ? TMath::Sqrt(evErr0) : kDummyError;
	    currLine[iz]->E[kResX] = evErr1>0 ? TMath::Sqrt(evErr1) : kDummyError;
	    //		
	  }
	  AliInfoF("Sect%2d bX=%3d bF=%3d DY vs Z: filled %d holes using %d values",isect,ix,ip, nmissY,npGood);
	}
	else printf("Sect%2d bX=%3d bF=%3d DY vs Z: FAILED to fill %d holes using %d values",isect,ix,ip, nmissY,npGood);	
      }
      //
      //
      if (nmissZ && (fNZ2XBins-nmissZ)>=minGoodPoints) { // recover points in Z
	int npGood = 0;
	for (int iz=0;iz<fNZ2XBins;iz++) {
	  if (missY[iz]) continue;
	  val0[npGood] = currLine[iz]->D[kResZ];
	  pos0[npGood] = currLine[iz]->stat[kVoxZ]; // average Z position
	  wgh0[npGood] = 1./(currLine[iz]->E[kResZ]*currLine[iz]->E[kResZ]);
	  npGood++;
	}
	Bool_t res = FitPoly2(pos0,val0,wgh0,npGood, resFit0, errFit0);
	if (res) {
	  for (int iz=0;iz<fNZ2XBins;iz++) {
	    if (!missZ[iz]) continue;
	    // evaluate in the center of the bin
	    double z = (isect>=kNSect ? -1.0f:1.0f)*(iz+0.5)*dz, z2=z*z, z3=z2*z, z4=z3*z;
	    currLine[iz]->D[kResZ] = resFit0[0]+z*resFit0[1]+z2*resFit0[2];
	    double evErr = errFit0[0] + errFit0[2]*z2 + errFit0[5]*z4
	      +            2.*(errFit0[1]*z + errFit0[3]*z2 + errFit0[4]*z3);
	    currLine[iz]->E[kResZ] = evErr>0 ? TMath::Sqrt(evErr) : kDummyError;
	  }
	  printf("Sect%2d bX=%3d bF=%3d DZ vs Z: filled %d holes using %d values\n",isect,ix,ip, nmissZ,npGood);
	}
	else printf("Sect%2d bX=%3d bF=%3d DZ vs Z: FAILED to fill %d holes using %d values\n",isect,ix,ip, nmissY,npGood);	
      }
      //
    } // loop in phi bins
  } // loop in x bins
  delete[] currLine;
  //
}

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
Int_t AliTPCDcalibRes::Smooth0(int isect)
{
  // apply linear regression kernel smoother 
  int cnt = 0;
  bres_t* sectData = fSectGVoxRes[isect];
  for (int ix=0;ix<fNXBins;ix++) {
    for (int ip=0;ip<fNY2XBins;ip++) {
      for (int iz=0;iz<fNZ2XBins;iz++) {  // extract line in z
	int binGlo = GetVoxGBin(ix,ip,iz);
	bres_t *vox = &sectData[binGlo];
	Bool_t res = GetSmoothEstimate(vox->bsec,vox->stat[kVoxX],vox->stat[kVoxF],vox->stat[kVoxZ],
				       BIT(kResX)|BIT(kResY)|BIT(kResZ), // at this moment we cannot smooth dispersion
				       vox->DS);
	vox->smooth = res;
	if (res) cnt++;
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
  const int kMinPointsTot = 4; // we fit 12 paremeters, each point provides 3 values
  const int kMaxTrials = 5; // max allowed iterations if neighbours are missing
  const float kTrialStep = 0.5;
  
  Bool_t doDim[kResDim] = {kFALSE};
  for (int i=0;i<kResDim;i++) {
    doDim[i] = (which&(0x1<<i))>0;
    if (doDim[i]) res[i] = 0;
  }
  //
  enum {kM00,
	kM10,kM11,
	kM20,kM21,kM22,
	kM30,kM31,kM32,kM33,kMN};
  double cmat[kResDim][kMN];

  //loop over neighbours which can contribute
  //
  bres_t *currClus[fNMaxNeighb];
  double *rhsX=&fLastSmoothingRes[0],*rhsY=&fLastSmoothingRes[4],*rhsZ=&fLastSmoothingRes[8],*rhsD=&fLastSmoothingRes[12];
  //
  int ix0,ip0,iz0;
  FindVoxel(x,p, isect<kNSect ? z : -z, ix0,ip0,iz0); // find nearest voxel
  bres_t* sectData = fSectGVoxRes[isect];
  int binCen = GetVoxGBin(ix0,ip0,iz0);  // global bin of nearest voxel
  bres_t* voxCen = &sectData[binCen]; // nearest voxel
  //
  int trial = 0, nbOK = 0;
  while(1)  {
    //
    memset(fLastSmoothingRes,0,kResDim*4*sizeof(double));
    if (trial>kMaxTrials) {printf("Trials limit reached\n"); return kFALSE;}

    memset(cmat,0,kResDim*10*sizeof(double));
    //
    nbOK=0; // accounted neighbours
    //
    float stepX = fStepKern[kVoxX]*(1. + kTrialStep*trial);
    float stepF = fStepKern[kVoxF]*(1. + kTrialStep*trial);
    float stepZ = fStepKern[kVoxZ]*(1. + kTrialStep*trial);
    //
    if (!voxCen->bvox[kVoxQ]) { // closest voxel has no data, increase smoothing step
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
    for (int i=ixMx-ixMn+1;i--;) nOccX[i]=0;
    for (int i=ipMx-ipMn+1;i--;) nOccF[i]=0;
    for (int i=izMx-izMn+1;i--;) nOccZ[i]=0;
    double u2Vec[3];
    for (int ix=ixMn;ix<=ixMx;ix++) {
      for (int ip=ipMn;ip<=ipMx;ip++) {
	for (int iz=izMn;iz<=izMx;iz++) {
	  //
	  int binNb = GetVoxGBin(ix,ip,iz);  // global bin
	  bres_t* voxNb = &sectData[binNb];
	  if (!voxNb->bvox[kVoxQ]) continue; // skip voxels w/o data
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
	  nOccX[ix-ixMn]++;
	  nOccF[ip-ipMn]++;
	  nOccZ[iz-izMn]++;
	  //
	  for (int id=0;id<kResDim;id++) {
	    if (!doDim[id]) continue;
	    double kernWD = kernW;
	    if (fUseErrInSmoothing) kernWD /= (voxNb->E[id]*voxNb->E[id]); // apart from the kernel value, account for the point error
	    double *cmatD = cmat[id];
	    cmatD[kM00] += kernWD;
	    cmatD[kM10] += kernWD*dx;   cmatD[kM11] += kernWD*dx*dx;
	    cmatD[kM20] += kernWD*df;   cmatD[kM21] += kernWD*dx*df;  cmatD[kM22] += kernWD*df*df;
	    cmatD[kM30] += kernWD*dz;   cmatD[kM31] += kernWD*dx*dz;  cmatD[kM32] += kernWD*df*dz;   cmatD[kM33] += kernWD*dz*dz;
	    double *rhsD = &fLastSmoothingRes[id*4];
	    rhsD[0] += kernWD*voxNb->D[id];
	    rhsD[1] += kernWD*voxNb->D[id]*dx;
	    rhsD[2] += kernWD*voxNb->D[id]*df;
	    rhsD[3] += kernWD*voxNb->D[id]*dz;	      
	  }
	  //
	  currClus[nbOK] = voxNb;
	  nbOK++;
	}
      }
    }
  
    // check if we have enough points in every dimension
    int npx=0,npp=0,npz=0;
    for (int i=ixMx-ixMn+1;i--;) if (nOccX[i]) npx++; 
    for (int i=ipMx-ipMn+1;i--;) if (nOccF[i]) npp++;
    for (int i=izMx-izMn+1;i--;) if (nOccZ[i]) npz++;
    if (npx<2 || npp<2 || npz<2 || nbOK<kMinPointsTot) {
      trial++;
      AliWarningF("Sector:%2d x=%.3f y/x=%.3f z/x=%.3f (iX:%d iY2X:%d iZ2X:%d)\n"
		  "not enough neighbours (need min %d) %d %d %d (tot: %d) | Steps: %.1f %.1f %.1f\n"
		  "trying to increase filter bandwidth (trial%d)\n",
		  isect,x,p,z,ix0,ip0,iz0,2,npx,npp,npz,nbOK,stepX,stepF,stepZ,trial);
      continue;
    }
    //
    Bool_t fitRes = kTRUE;
    //
    // solve system of linear equations
    AliSymMatrix mat(4);
    for (int id=0;id<kResDim;id++) {
      if (!doDim[id]) continue;
      mat.Reset();
      double *cmatD = cmat[id];
      double *rhsD = &fLastSmoothingRes[id*4];
      mat(0,0) = cmatD[kM00];
      mat(1,0) = cmatD[kM10];   mat(1,1) = cmatD[kM11];
      mat(2,0) = cmatD[kM20];   mat(2,1) = cmatD[kM21];  mat(2,2) = cmatD[kM22]; 
      mat(3,0) = cmatD[kM30];   mat(3,1) = cmatD[kM31];  mat(3,2) = cmatD[kM32];  mat(3,3) = cmatD[kM33];
      fitRes &= mat.SolveChol(rhsD);
      if (!fitRes) {
	trial++;
	AliWarningF("Sector:%2d x=%.3f y/x=%.3f z/x=%.3f (iX:%d iY2X:%d iZ2X:%d)\n"
		    "neighbours range used %d %d %d (tot: %d) | Steps: %.1f %.1f %.1f\n"
		    "Solution for smoothing Failed, trying to increase filter bandwidth (trial%d)",
		    isect,x,p,z,ix0,ip0,iz0,npx,npp,npz,nbOK,stepX,stepF,stepZ,trial);
	continue;
      }
      res[id] = rhsD[0];
      if (deriv) for (int j=0;j<3;j++) deriv[id*3 +j] = rhsD[j+1];
    }
    //
    break; // success
  } // end of loop over allowed trials

  return kTRUE;

}


//________________________________________________________________
Bool_t AliTPCDcalibRes::GetSmoothEstimateDim(int isect, float x, float p, float z, int dim, 
					     float& res, float *deriv)
{
  // get smooth estimate of single dimension dim for point in sector coordinates (x,y/x,z/x)
  // smoothing results also saved in the fLastSmoothingRes (allow derivative calculation)
  //
  const int kMinPointsTot = 4; // we fit 12 paremeters, each point provides 3 values
  const int kMaxTrials = 5; // max allowed iterations if neighbours are missing
  const float kTrialStep = 0.5;

  res = 0;
  //
  double cmat[10];
  double &m00=cmat[0], 
    &m10=cmat[1], &m11=cmat[2], 
    &m20=cmat[3], &m21=cmat[4], &m22=cmat[5], 
    &m30=cmat[6], &m31=cmat[7], &m32=cmat[8], &m33=cmat[9];

  //loop over neighbours which can contribute
  //
  bres_t *currClus[fNMaxNeighb];
  double *rhs = &fLastSmoothingRes[dim*4];
  //
  int ix0,ip0,iz0;
  FindVoxel(x,p, isect<kNSect ? z : -z, ix0,ip0,iz0); // find nearest voxel
  bres_t* sectData = fSectGVoxRes[isect];
  int binCen = GetVoxGBin(ix0,ip0,iz0);  // global bin of nearest voxel
  bres_t* voxCen = &sectData[binCen]; // nearest voxel
  //
  int trial = 0, nbOK = 0;
  while(1)  {
    //
    memset(rhs,0,4*sizeof(double));
    if (trial>kMaxTrials) {printf("Trials limit reached\n"); return kFALSE;}
    memset(cmat,0,10*sizeof(double));
    nbOK=0; // accounted neighbours
    //
    float stepX = fStepKern[kVoxX]*(1. + kTrialStep*trial);
    float stepF = fStepKern[kVoxF]*(1. + kTrialStep*trial);
    float stepZ = fStepKern[kVoxZ]*(1. + kTrialStep*trial);
    //
    if (!voxCen->bvox[kVoxQ]) { // closest voxel has no data, increase smoothing step
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
    for (int i=ixMx-ixMn+1;i--;) nOccX[i]=0;
    for (int i=ipMx-ipMn+1;i--;) nOccF[i]=0;
    for (int i=izMx-izMn+1;i--;) nOccZ[i]=0;
    double u2Vec[3];
    for (int ix=ixMn;ix<=ixMx;ix++) {
      for (int ip=ipMn;ip<=ipMx;ip++) {
	for (int iz=izMn;iz<=izMx;iz++) {
	  //
	  int binNb = GetVoxGBin(ix,ip,iz);  // global bin
	  bres_t* voxNb = &sectData[binNb];
	  if (!voxNb->bvox[kVoxQ]) continue; // skip voxels w/o data
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
	  nOccX[ix-ixMn]++;
	  nOccF[ip-ipMn]++;
	  nOccZ[iz-izMn]++;
	  //
	  // apart from the kernel value, we may account for the point error
	  if (fUseErrInSmoothing) kernW /= (voxNb->E[dim]*voxNb->E[dim]);
	  //
	  m00 += kernW;
	  m10 += kernW*dx;   m11 += kernW*dx*dx;
	  m20 += kernW*df;   m21 += kernW*dx*df;  m22 += kernW*df*df;
	  m30 += kernW*dz;   m31 += kernW*dx*dz;  m32 += kernW*df*dz;   m33 += kernW*dz*dz;
	  //
	  rhs[0] += kernW*voxNb->D[dim];
	  rhs[1] += kernW*voxNb->D[dim]*dx;
	  rhs[2] += kernW*voxNb->D[dim]*df;
	  rhs[3] += kernW*voxNb->D[dim]*dz;
	  //
	  currClus[nbOK] = voxNb;
	  nbOK++;
	}
      }
    }
    //  
    // check if we have enough points in every dimension
    int npx=0,npp=0,npz=0;
    for (int i=ixMx-ixMn+1;i--;) if (nOccX[i]) npx++; 
    for (int i=ipMx-ipMn+1;i--;) if (nOccF[i]) npp++;
    for (int i=izMx-izMn+1;i--;) if (nOccZ[i]) npz++;
    if (npx<2 || npp<2 || npz<2 || nbOK<kMinPointsTot) {
      trial++;
      AliWarningF("Sector:%2d Dim%d x=%.3f y/x=%.3f z/x=%.3f (iX:%d iY2X:%d iZ2X:%d)\n"
		  "not enough neighbours (need min %d) %d %d %d (tot: %d) | Steps: %.1f %.1f %.1f\n"
		  "trying to increase filter bandwidth (trial%d)\n",
		  isect,dim,x,p,z,ix0,ip0,iz0,2,npx,npp,npz,nbOK,stepX,stepF,stepZ,trial);
      continue;
    }
    //
    Bool_t fitRes = kTRUE;
    //
    // solve system of linear equations
    AliSymMatrix mat(4);
    mat(0,0) = m00;
    mat(1,0) = m10;   mat(1,1) = m11;
    mat(2,0) = m20;   mat(2,1) = m21;  mat(2,2) = m22; 
    mat(3,0) = m30;   mat(3,1) = m31;  mat(3,2) = m32;  mat(3,3) = m33;
    fitRes &= mat.SolveChol(rhs);
    if (!fitRes) {
      trial++;
      AliWarningF("Sector:%2d Dim%d x=%.3f y/x=%.3f z/x=%.3f (iX:%d iY2X:%d iZ2X:%d)\n"
		  "neighbours range used %d %d %d (tot: %d) | Steps: %.1f %.1f %.1f\n"
		  "Solution for smoothing Failed, trying to increase filter bandwidth (trial%d)",
		  isect,dim,x,p,z,ix0,ip0,iz0,npx,npp,npz,nbOK,stepX,stepF,stepZ,trial);
      continue;
    }
    //
    break; // success
  } // end of loop over allowed trials
  res = rhs[0];
  //
  if (deriv) for (int j=0;j<3;j++) deriv[j] = rhs[j+1]; // derivatives are requested

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
    w = u2<25*3 ? TMath::Exp(-u2)/TMath::Sqrt(2.*TMath::Pi()) : 0;
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
  fNMaxNeighb = 2*(2*fStepKern[kVoxX]+1)*(2*fStepKern[kVoxF]+1)*(2*fStepKern[kVoxZ]+1);
}


//=========================================================================

//_____________________________________________
void AliTPCDcalibRes::CreateCorrectionObject()
{
  // create correction object for given time slice

  AliSysInfo::AddStamp("CreateCorrectionObject",0,0,0,0);
  
  TString name = Form("run%d_%lld_%lld",fRun,fTMin,fTMax);
  fChebCorr = new AliTPCChebCorr(name.Data(),name.Data(),
				 fChebPhiSlicePerSector,fChebZSlicePerSide,1.0f);
  fChebCorr->SetUseFloatPrec(kFALSE);
  fChebCorr->SetTimeStampStart(fTMin);
  fChebCorr->SetTimeStampEnd(fTMax);
  fChebCorr->SetTimeDependent(kFALSE);
  fChebCorr->SetUseZ2R(kTRUE);
  //
  SetUsedInstance(this);
  fChebCorr->Parameterize(trainCorr,kResDim,fNPCheb,fChebPrecD);
  //
  AliSysInfo::AddStamp("CreateCorrectionObject",1,0,0,0);
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
  fQ2PTBound[0] = -fMaxQ2Pt;
  fQ2PTBound[1] = -fMidQ2Pt;
  fQ2PTBound[2] =  0;
  fQ2PTBound[3] =  fMidQ2Pt;
  fQ2PTBound[4] =  fMaxQ2Pt;
  //
  // fNXYBinsProd = fNXBins*fNY2XBins;
  // fBinMinQ = new Float_t[fNXYBinsProd];
  // fBinDQI  = new Float_t[fNXYBinsProd];
  // fBinDQ   = new Float_t[fNXYBinsProd];
  //
  const float kMaxY2X = TMath::Tan(0.5f*kSecDPhi);

  for (int ix=0;ix<fNXBins;ix++) {
    float x = GetX(ix);
    fMaxY2X[ix] = kMaxY2X - kDeadZone/x;
    fDY2XI[ix] = fNY2XBins / (2.f*fMaxY2X[ix]);
    fDY2X[ix] = 1.f/fDY2XI[ix];
    // for (int iy=0;iy<fNY2XBins;iy++) {
    //   float y = GetY2X(ix,iy)*x;
    //   float tgMn = tgpXY(x,y,-fMaxQ2Pt,fBz);
    //   float tgMx = tgpXY(x,y, fMaxQ2Pt,fBz);
    //   if (tgMn>tgMx) swap(tgMn,tgMx);
    //   int ixy = ix*fNY2XBins + iy;
    //   fBinMinQ[ixy] = TMath::Abs(fBz)>0.01 ? tgMn : -0.5;
    //   fBinDQ[ixy]   = TMath::Abs(fBz)>0.01 ? (tgMx-tgMn)/kNQBins : 1.;
    //   fBinDQI[ixy]  = 1./fBinDQ[ixy];
    // }
  }
  //
  fDZ2XI = fNZ2XBins/kMaxZ2X;
  fDZ2X  = 1.0f/fDZ2XI;
  //
  // inverse bin sizes for residuals
  fDeltaYbinI  = fNDeltaYBins/(2.0f*fMaxDY);
  fDeltaZbinI  = fNDeltaZBins/(2.0f*fMaxDZ);
  //


  fNGVoxPerSector = fNY2XBins*fNZ2XBins*fNXBins;
  fNBProdSectG[0] = fNY2XBins*fNZ2XBins;
  fNBProdSectG[1] = fNZ2XBins;

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
Bool_t AliTPCDcalibRes::FindVoxelBin(int sectID, float tgsl, float q2pt, float x, float y, float z, UChar_t bin[kVoxHDim],float voxVars[kVoxHDim])
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
  // track inclination at pad row
  // the one calculated from q/pt is less precise
  voxVars[kVoxQ] = tgsl; // tgpXY(x,y,q2pt,fBz);
  int binQ = GetQBin(q2pt);  //GetQBin(q2pt,binX,binY)
  if (binQ<0) return kFALSE;
  bin[kVoxQ] = binQ;
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
    printf("training Sector %d\n",sector);
    return;
  }
  //
  float x = AliTPCDcalibRes::GetTPCRowX(row);
  float dist[AliTPCDcalibRes::kResDim], deriv[AliTPCDcalibRes::kResDim*3];

  float y2x = tzLoc[0];
  float z2x = tzLoc[1];
  //
  Bool_t res = AliTPCDcalibRes::GetUsedInstance()->GetSmoothEstimate(sector, x, y2x, z2x, 
								     0xff, dist);
  if (!res) { printf("Failed to evaluate smooth distortion\n"); exit(1); }

  /*
  // Marian stored Z track coordinate instead of cluster one, need to correct for this
  if (fApplyZt2Zc) {
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
