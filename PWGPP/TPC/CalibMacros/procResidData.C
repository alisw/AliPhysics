#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliTPCDcalibRes.h"
#include <TString.h>
#include <TSystem.h>
#endif

AliTPCDcalibRes*  CreateSetCalib(int run,int tmin=0,int tmax=0x7fffffff,const char* inp="lst.txt");
AliTPCDcalibRes* clb = 0;

//================================================================================================
void procResidData(int mode,                       // processing mode
		   int run,                        // run number
		   int tmin=0,                     // start time
		   int tmax=0x7ffffff,             // end time
		   const char* inp="lst.txt",       // input residuals list
		   Bool_t completeMissing=kFALSE   // if info for given stage is missing, create it
		   ) 
{
  clb = 0;
  //
  if (mode==0) {   // vdrift and time bins definition
    printf("VDrift extraction and time bins defintion for run %d\n",run);
    clb = CreateSetCalib(run,tmin,tmax,inp);
    clb->CalibrateVDrift();
    clb->Save();
    return;
  }
  //
  if (mode==1) {  // write vdrift OCDB object from existing processor object 
    printf("Creation of VDrift OCDB object run %d\n",run);
    if (!clb) clb = AliTPCDcalibRes::Load();
    if (!clb) {
      printf("did not find preprocessed object to extract VDrift\n");
      if (!completeMissing) return;
      printf("autocompletion mode ON: creating object and extractin vdrift");
      clb = CreateSetCalib(run,tmin,tmax,inp);
      clb->CalibrateVDrift();
      clb->Save();
    }
    clb->MakeVDriftOCDB(gSystem->Getenv("targetOCDBDir"));
    return;
  }
  //
  if (mode==2) {  // creation of distortion maps
    printf("Creation of distortion map for run %d\n",run);
    clb = AliTPCDcalibRes::Load();
    if (!clb) {
      printf("did not find preprocessed object to extract maps\n");
      if (!completeMissing) return;
      printf("autocompletion mode ON: creating object and extracting vdrift and maps");
      clb = CreateSetCalib(run,tmin,tmax,inp);
      clb->CalibrateVDrift();
      clb->Save();      
    }
    clb->SetTMinMax(tmin,tmax);
    clb->ProcessFromDeltaTrees();
    clb->Save();
    return;
  }
  //
  if (mode==3) { // closure test
    printf("Performing closure test for run %d\n",run);
    clb = AliTPCDcalibRes::Load();
    if (!clb) {
      printf("did not find preprocessed object to perform closure test\n");
      if (!completeMissing) return;
      printf("autocompletion mode ON: creating object and extracting vdrift and maps");
      clb = CreateSetCalib(run,tmin,tmax,inp);
      clb->CalibrateVDrift();
      clb->ProcessFromDeltaTrees();
      clb->Save();      
    }
    TString envs = gSystem->Getenv("distNTracksClosureTest");
    if (envs.IsDigit()) clb->SetMaxTracks(envs.Atoi());
    else                clb->SetMaxTracks(clb->GetMinTrackToUse());
    clb->ClosureTest();
    return;
  }
  //
  printf("Mode %d is not supported\n",mode);
  //
}


//________________________________________________________________________
AliTPCDcalibRes*  CreateSetCalib(int run,int tmin,int tmax,const char* inp)
{
  // create calibration object and init with provided env.variables
  AliTPCDcalibRes* clbn = new AliTPCDcalibRes(run);
  clbn->SetTMinMax(tmin,tmax);
  clbn->SetResidualList(inp);
  //
  TString envs;

  // detectors to use >>>>>>>>>>>>>>>>
  envs = gSystem->Getenv("distUseDet");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetMinTracks %s",envs.Data());
    clbn->SetExternalDetectors(envs.Atoi());
  }
  //
  envs = gSystem->Getenv("useTOFBC");    
  if (!envs.IsNull()) {
    envs.ToLower();
    Bool_t useTOFBC=clbn->GetUseTOFBC();
    if      (envs.Contains("true")) useTOFBC = kTRUE;
    else if (envs.Contains("false")) useTOFBC = kFALSE;
    else {
      ::Error("CreateResCalib","Wrong useTOFBC = %s",envs.Data());
      exit(1);
    }
    //
    ::Info("CreateResCalib","SetUseTOFBC %s",envs.Data());
    clbn->SetUseTOFBC(useTOFBC);
  }
  // detectors to use <<<<<<<<<<<<<<<<


  // binning >>>>>>>>>>>>>>>>>>>>>>>>>
  envs = gSystem->Getenv("distNBinsZ");;
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetNZ2XBins %s",envs.Data());
    clbn->SetNZ2XBins(envs.Atoi());
  }
  //
  envs = gSystem->Getenv("distNBinsY");;
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetNY2XBins %s",envs.Data());
    clbn->SetNZ2XBins(envs.Atoi());
  }
  //
  // binning <<<<<<<<<<<<<<<<<<<<<<<<<
  //
  // stat limits >>>>>>>>>>>>>>>>>>>>>
  envs = gSystem->Getenv("distTimeBin");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetNominalTimeBin %s",envs.Data());
    clbn->SetNominalTimeBin(envs.Atoi());
  }
  //
  envs = gSystem->Getenv("distTimeBinPrec");
  if (envs.IsFloat()) {
    ::Info("CreateResCalib","SetNominalTimeBinPrec %s",envs.Data());
    clbn->SetNominalTimeBinPrec(envs.Atof());
  }
  //
  envs = gSystem->Getenv("distMaxTracks");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetMaxTracks %s",envs.Data());
    clbn->SetMaxTracks(envs.Atoi());
  }
  //
  envs = gSystem->Getenv("distMinTracks");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetMinTrackToUse %s",envs.Data());
    clbn->SetMinTrackToUse(envs.Atoi());
  }
  //
  // stat limits <<<<<<<<<<<<<<<<<<<<<
  //
  // Kernel Smoother settings >>>>>>>>
  int kernelType=-1000,xSmtPol2=-1000,ySmtPol2=-1000,zSmtPol2=-1000;
  float kernelWX=2.1,kernelWY=2.1,kernelWZ=1.7;
  //
  envs = gSystem->Getenv("distKernelType");
  if (envs.IsDigit()) {kernelType = envs.Atoi();}
  envs = gSystem->Getenv("distKernelWX");
  if (envs.IsFloat()) kernelWX = envs.Atof();  
  envs = gSystem->Getenv("distKernelWY");
  if (envs.IsFloat()) kernelWY = envs.Atof();  
  envs = gSystem->Getenv("distKernelWZ");
  if (envs.IsFloat()) kernelWZ = envs.Atof();  
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
  envs = gSystem->Getenv("distKernelPol2X");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetSmoothPol2X %s",envs.Data());
    clb->SetSmoothPol2(AliTPCDcalibRes::kVoxX,envs.Atoi()>0);
  }
  //
  envs = gSystem->Getenv("distKernelPol2F");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetSmoothPol2Y %s",envs.Data());
    clb->SetSmoothPol2(AliTPCDcalibRes::kVoxF,envs.Atoi()>0);
  }
  //
  envs = gSystem->Getenv("distKernelPol2Z");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetSmoothPol2Z %s",envs.Data());
    clb->SetSmoothPol2(AliTPCDcalibRes::kVoxZ,envs.Atoi()>0);
  }
  //
  // Kernel Smoother settings <<<<<<<<
  //
  // settings for VDrift extraction >>
  //
  envs = gSystem->Getenv("driftDeltaT");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetDeltaTVDrift %s",envs.Data());
    clbn->SetDeltaTVDrift(envs.Atoi());
  }
  //
  envs = gSystem->Getenv("driftSigmaT");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetSigmaTVDrift %s",envs.Data());
    clbn->SetSigmaTVDrift(envs.Atoi());
  }
  //
  envs = gSystem->Getenv("distTimeBin");
  if (envs.IsDigit()) {
    ::Info("CreateResCalib","SetNominalTimeBin %s",envs.Data());
    clbn->SetNominalTimeBin(envs.Atoi());
  }
  //
  envs = gSystem->Getenv("distTimeBinPrec");
  if (envs.IsFloat()) {
    ::Info("CreateResCalib","SetNominalTimeBinPrec %s",envs.Data());
    clbn->SetNominalTimeBinPrec(envs.Atof());
  }
  //
  return clbn;
}
