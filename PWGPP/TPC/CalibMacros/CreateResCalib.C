void CreateResCalib(int run=245231
		    ,Long64_t tmin=0
		    ,Long64_t tmax= 9999999999
		    ,const char * resList="lst.txt"
		    ,const char * cdb="raw://"
		    ,Bool_t useTOFBC = kFALSE
	  )
{
  //
  AliTPCDcalibRes* clb = new AliTPCDcalibRes(run, tmin, tmax, resList);
  clb->SetOCDBPath(cdb);
  clb->SetUseTOFBC(useTOFBC);
  //  
  CheckResCalibEnvVars(clb);
  //
  clb->ProcessFromDeltaTrees();
  clb->Save();
}


//__________________________________________
void CheckResCalibEnvVars(AliTPCDcalibRes* clb)
{
  TString envs;
  // if < -999 use default values
  int nBinsZ=-1000,nBinsY=-1000,nMaxTracks=-1000,nMinTracks=-1000;
  int kernelType=-1000;
  int xSmtPol2=-1000,ySmtPol2=-1000,zSmtPol2=-1000;
  float kernelWX=2.1,kernelWY=2.1,kernelWZ=1.7;
  int detToUse = AliTPCDcalibRes::kUseTRDonly;
  //

  // detectors to use >>>>>>>>>>>>>>>>
  envs = gSystem->Getenv("distUseDet");
  if (envs.IsDigit()) detToUse = envs.Atoi();
  //
  // detectors to use <<<<<<<<<<<<<<<<


  // binning >>>>>>>>>>>>>>>>>>>>>>>>>
  envs = gSystem->Getenv("distNBinsZ");
  if (envs.IsDigit()) nBinsZ = envs.Atoi();
  //
  envs = gSystem->Getenv("distNBinsY");
  if (envs.IsDigit()) nBinsY = envs.Atoi();
  //
  // binning <<<<<<<<<<<<<<<<<<<<<<<<<
  //
  // stat limits >>>>>>>>>>>>>>>>>>>>>
  envs = gSystem->Getenv("distMaxTracks");
  if (envs.IsDigit()) nMaxTracks = envs.Atoi();
  //
  envs = gSystem->Getenv("distMinTracks");
  if (envs.IsDigit()) nMinTracks = envs.Atoi();
  //
  // stat limits <<<<<<<<<<<<<<<<<<<<<
  //
  // Kernel Smoother settings >>>>>>>>
  envs = gSystem->Getenv("distKernelType");
  if (envs.IsDigit()) kernelType = envs.Atoi();  
  //
  envs = gSystem->Getenv("distKernelWX");
  if (envs.IsFloat()) kernelWX = envs.Atof();  
  //
  envs = gSystem->Getenv("distKernelWY");
  if (envs.IsFloat()) kernelWY = envs.Atof();  
  //
  envs = gSystem->Getenv("distKernelWZ");
  if (envs.IsFloat()) kernelWZ = envs.Atof();  
  //
  envs = gSystem->Getenv("distKernelPol2X");
  if (envs.IsDigit()) xSmtPol2 = envs.Atoi();  
  //
  envs = gSystem->Getenv("distKernelPol2Y");
  if (envs.IsDigit()) ySmtPol2 = envs.Atoi();  
  //
  envs = gSystem->Getenv("distKernelPol2Z");
  if (envs.IsDigit()) zSmtPol2 = envs.Atoi();  
  //
  // Kernel Smoother settings <<<<<<<<
  //
  // set values

  if (detToUse>=0 && detToUse<AliTPCDcalibRes::kNExtDetComb) {
    ::Info("CreateResCalib","SetExternalDetectors %d",detToUse);
    clb->SetExternalDetectors(detToUse);
  }
  else {
    ::Info("CreateResCalib","Invalid value for distUseDet=%d, keep default %d",detToUse,clb->GetExternalDetectors());
  }

  if (nMaxTracks>0) {
    ::Info("CreateResCalib","SetMaxTracks %d",nMaxTracks);
    clb->SetMaxTracks(nMaxTracks);
  }
  if (nMinTracks>0) {
    ::Info("CreateResCalib","SetMinTracks %d",nMinTracks);
    clb->SetMaxTracks(nMinTracks);
  }
  if (nBinsZ>0) {
    ::Info("CreateResCalib","SetNZ2XBins %d",nBinsZ);
    clb->SetNZ2XBins(nBinsZ);
  }
  if (nBinsY>0) {
    ::Info("CreateResCalib","SetNY2XBins %d",nBinsY);
    clb->SetNY2XBins(nBinsY);
  }
  //
  if (kernelType>=0) {
    if (kernelType!=AliTPCDcalibRes::kEpanechnikovKernel &&
	kernelType!=AliTPCDcalibRes::kGaussianKernel) {
      ::Error("CreateResCalib","Wrong kernel type %d",kernelType);
      exit(1);
    }
    ::Info("CreateResCalib","SetKernelType(%d,%.2f,%.2f,%.2f)",
	   kernelType,kernelWX,kernelWY,kernelWZ);
    clb->SetKernelType(kernelType,kernelWX,kernelWY,kernelWZ);
  }
  //
  if (xSmtPol2>=0) {
    ::Info("CreateResCalib","SetSmoothPol2(%d,%d)",AliTPCDcalibRes::kVoxX,xSmtPol2>0);
    clb->SetSmoothPol2(AliTPCDcalibRes::kVoxX,xSmtPol2>0);
  }
  //
  if (ySmtPol2>=0) {
    ::Info("CreateResCalib","SetSmoothPol2(%d,%d)",AliTPCDcalibRes::kVoxF,ySmtPol2>0);
    clb->SetSmoothPol2(AliTPCDcalibRes::kVoxF,ySmtPol2>0);
  }
  //
  if (xSmtPol2>=0) {
    ::Info("CreateResCalib","SetSmoothPol2(%d,%d)",AliTPCDcalibRes::kVoxZ,zSmtPol2>0);
    clb->SetSmoothPol2(AliTPCDcalibRes::kVoxZ,zSmtPol2>0);
  }
  //
  TString useTOFBCStr = gSystem->Getenv("useTOFBC");    
  if (!useTOFBCStr.IsNull()) {
    useTOFBCStr.ToLower();
    Bool_t useTOFBC=clb->GetUseTOFBC();
    if      (useTOFBCStr.Contains("true")) useTOFBC = kTRUE;
    else if (useTOFBCStr.Contains("false")) useTOFBC = kFALSE;
    else {
      ::Error("CreateResCalib","Wrong useTOFBC = %s",useTOFBCStr.Data());
      exit(1);
    }
    //
    ::Info("CreateResCalib","SetUseTOFBC %d",useTOFBC);
    clb->SetUseTOFBC(useTOFBC);
  }
  //
}
