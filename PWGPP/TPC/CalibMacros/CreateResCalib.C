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
  float kernelWX=2.5,kernelWY=2.5,kernelWZ=2.1;
  //
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
  // Kernel Smoother settings <<<<<<<<
  //
  // set values
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
  TString useTOFBCStr = gSystem->Getenv("useTOFBC");    
  if (!useTOFBCStr.IsNull()) {
    useTOFBCStr.ToLower();
    Bool_t useTOFBC=clb->GetUseTOFBC();
    if      (useTOFBCStr.Contains("true")) 
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
