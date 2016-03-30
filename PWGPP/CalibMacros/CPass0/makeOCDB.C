/*
  macro to extract the OCDB entries

  input: CalibObjects.root
  ouput: TimeGain and TimeVdrift calibration objects for TPC and TRD

  Example:
  .L $ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/makeOCDB.C
  makeOCDB("105160");

*/

void printCalibStat(Int_t run, const char * fname,  TTreeSRedirector * pcstream);

void makeOCDB(Int_t runNumber, TString  targetOCDBstorage="", TString sourceOCDBstorage="raw://", Int_t detectorBitsQualityFlag = -1)
{
  //
  // extract OCDB entries for detectors participating in the calibration for the current run
  //

  gROOT->Macro("$ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/LoadLibraries.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/ConfigCalibTrain.C");

  // switch off log info
  AliLog::SetClassDebugLevel("AliESDEvent",0);

  // config GRP
  printf("runNumber from runCalibTrain = %d\n",runNumber);
  ConfigCalibTrain(runNumber, sourceOCDBstorage.Data());

  // check the presence of the detectors
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());
  if (!grpData) {printf("Failed to get GRP data for run",runNumber); return;}
  Int_t activeDetectors = grpData->GetDetectorMask();
  TString detStr = AliDAQ::ListOfTriggeredDetectors(activeDetectors);
  printf("Detectors in the data:\n%s\n",detStr.Data());
  TString LHCperiod = grpData->GetLHCPeriod();
  Bool_t isLHC10 =  LHCperiod.Contains("LHC10");
  Bool_t isLHC11 =  LHCperiod.Contains("LHC11");
  Bool_t isLHC12 =  LHCperiod.Contains("LHC12");
  Bool_t isLHC13 =  LHCperiod.Contains("LHC13");
  Bool_t isLHC13b =  LHCperiod.Contains("LHC13b");
  Bool_t isLHC13c =  LHCperiod.Contains("LHC13c");
  printf("LHCperiod:%s\n isLHC10:%d isLHC11:%d isLHC12:%d isLHC13:%d isLHC13b:%d isLHC13c:%d\n",LHCperiod.Data(),(Int_t)isLHC10,(Int_t)isLHC11,(Int_t)isLHC12,(Int_t)isLHC13,(Int_t)isLHC13b,(Int_t)isLHC13c);

  // Steering Tasks - set output storage
  // DefaultStorage set already before - in ConfigCalibTrain.C

  // Setting the mirror SEs for the default storage
  TString mirrorsStr("ALICE::CERN::OCDB,ALICE::FZK::SE,ALICE::LLNL::SE");
  AliCDBManager::Instance()->SetMirrorSEs(mirrorsStr.Data());
  printf("List of mirror SEs set to: \"%s\"\n",mirrorsStr.Data());

  // activate target OCDB storage
  AliCDBStorage* targetStorage = 0x0;
  if (targetOCDBstorage.Length()==0) {
    targetOCDBstorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    targetStorage = AliCDBManager::Instance()->GetStorage(targetOCDBstorage.Data());
  }
  else if (targetOCDBstorage.CompareTo("same",TString::kIgnoreCase) == 0 ){
    targetStorage = AliCDBManager::Instance()->GetDefaultStorage();
  }
  else {
    targetStorage = AliCDBManager::Instance()->GetStorage(targetOCDBstorage.Data());
  }
  printf("** targetOCDBstorage: \"%s\"\n",targetOCDBstorage.Data());

  // specific storage for TPC/Calib/Correction entry
  if (gSystem->AccessPathName("TPC", kFileExists)==0) {  
    AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Correction","local://");
  }

  // Magnetic field
  AliMagF* fld = TGeoGlobalMagField::Instance()->GetField();
  Double_t bz = fld->SolenoidField();
  Bool_t isMagFieldON = kTRUE;
  if (TMath::Abs(bz)>0) {
    printf("Mag field is %f --> ON\n", bz);
  }
  else {
    isMagFieldON = kFALSE;
    printf("Mag field is %f --> OFF\n", bz);
  }

  // Quality flags
  Bool_t TPC_qf = kTRUE;
  Bool_t TOF_qf = kTRUE;
  Bool_t TRD_qf = kTRUE;
  Bool_t T0_qf  = kTRUE;
  Bool_t SDD_qf = kTRUE;
  Bool_t SPD_qf = kTRUE;

  /*
    // RS Commenting this to sync with working version from alidaq
  if (detectorBitsQualityFlag != -1){
    TPC_qf = ((detectorBitsQualityFlag & AliDAQ::kTPC_QF) == AliDAQ::kTPC_QF)? kTRUE : kFALSE;
    TOF_qf = ((detectorBitsQualityFlag & AliDAQ::kTOF_QF) == AliDAQ::kTOF_QF)? kTRUE : kFALSE;
    TRD_qf = ((detectorBitsQualityFlag & AliDAQ::kTRD_QF) == AliDAQ::kTRD_QF)? kTRUE : kFALSE;
    T0_qf  = ((detectorBitsQualityFlag & AliDAQ::kT0_QF)  == AliDAQ::kT0_QF)?  kTRUE : kFALSE;
    SDD_qf = ((detectorBitsQualityFlag & AliDAQ::kSDD_QF) == AliDAQ::kSDD_QF)? kTRUE : kFALSE;
    SPD_qf = ((detectorBitsQualityFlag & AliDAQ::kSPD_QF) == AliDAQ::kSPD_QF)? kTRUE : kFALSE;
  } 
  */   

  Printf("Quality flags: detectorBitsQualityFlag = %d, TPC = %d, TOF = %d, TRD = %d, T0 = %d, SDD = %d, SPD = %d", detectorBitsQualityFlag, (Int_t)TPC_qf, (Int_t)TOF_qf, (Int_t)TRD_qf, (Int_t)T0_qf, (Int_t)SDD_qf, (Int_t)SPD_qf);

  // ===========================================================================
  // ===| TPC part |============================================================
  //
  AliTPCPreprocessorOffline *procesTPC = 0;
  if (detStr.Contains("TPC") && TPC_qf){
    Printf("\n******* Calibrating TPC *******");

    // ===| set up residual storage |===========================================
    TString targetStorageResidual="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";

    // --- check for overwrites
    const TString targetStorageResidualEnv(gSystem->Getenv("targetStorageResidual"));
    if (!targetStorageResidualEnv.IsNull())  targetStorageResidual=targetStorageResidualEnv;
    AliCDBStorage *residualStorage = AliCDBManager::Instance()->GetStorage(targetStorageResidual.Data());

    // ===| set up TPC calibration object |=====================================
    procesTPC = new AliTPCPreprocessorOffline;

    // ---| set up gain calibratin type |---------------------------------------
    //
    // default is Full Calibration in CPass0
    // will be overwritte by mergeMakeOCDB.byComponent.perStage.sh
    // NOTE: This must be consistent to the settings in runCPass*.sh (runCalibTrain.C)
    //
    procesTPC->SetGainCalibrationType(AliTPCPreprocessorOffline::kFullGainCalib);

    // --- check for overwrites from environment variable
    //
    const TString sGainTypeFromEnv(gSystem->Getenv("TPC_CPass0_GainCalibType"));
    if (!sGainTypeFromEnv.IsNull()) {
      if (!procesTPC->SetGainCalibrationType(sGainTypeFromEnv)) {
        ::Fatal("makeOCDB","Could not set up gain calibration type from environment variable TPC_CPass0_GainCalibType: %s",sGainTypeFromEnv.Data());
      }

      ::Info("makeOCDB","Setting gain calibration type from environment variable TPC_CPass0_GainCalibType: %d", Int_t(procesTPC->GetGainCalibrationType()));
    }

    // ---| switch on parameter validation |------------------------------------
    procesTPC->SetTimeGainRange(0.5,5.0);
    procesTPC->SetMaxVDriftCorr(0.2); 
    //procesTPC->SetMinTracksVdrift(100000);
    procesTPC->SwitchOnValidation();

    // ===| Run calibration |===================================================
    //
    // ---| Make time gain calibration |----------------------------------------
    if (isMagFieldON) procesTPC->CalibTimeGain("CalibObjects.root", runNumber, runNumber, targetStorage, residualStorage);

    // ---| Make vdrift calibration |-------------------------------------------
    procesTPC->CalibTimeVdrift("CalibObjects.root",runNumber,runNumber,targetStorage);
  }
  else {
    Printf("\n******* NOT Calibrating TPC: detStr = %s, TPC_qf = %d *******", detStr.Data(), (Int_t)TPC_qf);
  }

  // ===========================================================================
  // ===| TOF part |============================================================
  //
  AliTOFAnalysisTaskCalibPass0 *procesTOF=0;
  if (detStr.Contains("TOF") && detStr.Contains("TPC") && TOF_qf){
    procesTOF = new AliTOFAnalysisTaskCalibPass0;
    Printf("\n******* Calibrating TOF *******");
    if (isMagFieldON) procesTOF->ProcessOutput("CalibObjects.root", targetStorage);
    else {
      printf("Not calibrating TOF in case of mag field OFF\n");
    }
  }
  else {
    Printf("\n******* NOT Calibrating TOF: detStr = %s, TOF_qf = %d *******", detStr.Data(), (Int_t)TOF_qf);
  }

  // T0 part
  AliT0PreprocessorOffline *procesT0= 0;
  if (detStr.Contains("T0") && T0_qf) {
    Printf("\n******* Calibrating T0 *******");
    // Make  calibration of channels offset
    procesT0 = new AliT0PreprocessorOffline;
    if(isLHC10)
      procesT0->CalibOffsetChannels("CalibObjects.root",runNumber, runNumber, targetStorage);
    else 
      procesT0->Process("CalibObjects.root",runNumber, runNumber, targetStorage);
  }
  else {
    Printf("\n******* NOT Calibrating T0: detStr = %s, T0_qf = %d *******", detStr.Data(), (Int_t)T0_qf);
  }

  //TRD part
  AliTRDPreprocessorOffline *procesTRD = 0;
  if (detStr.Contains("TRD") && detStr.Contains("TPC") && TRD_qf){
    Printf("\n******* Calibrating TRD *******");
    procesTRD = new  AliTRDPreprocessorOffline;
    if(isLHC10 || isLHC13b || isLHC13c) procesTRD->SetSwitchOnChamberStatus(kFALSE);
    procesTRD->SetLinearFitForVdrift(kTRUE);
    procesTRD->SetMinStatsVdriftT0PH(600*10);
    procesTRD->SetMinStatsVdriftLinear(50);
    procesTRD->SetMinStatsGain(600);
    procesTRD->SetLimitValidateNoData(100);
    procesTRD->SetLimitValidateBadCalib(90);
    procesTRD->SetMinTimeOffsetValidate(-2.1);
    procesTRD->SetAlternativeDriftVelocityFit(kTRUE);
    if((!isLHC10) && (!isLHC11) && (!isLHC12) && (!isLHC13)) {
      printf("Run II\n");
      procesTRD->SetT0Shift1(0.2524132);// release the condition on the first bin and last bins
    }
    procesTRD->Init("CalibObjects.root");
    Int_t versionVdriftUsed = procesTRD->GetVersionVdriftUsed();
    Int_t subversionVdriftUsed = procesTRD->GetSubVersionVdriftUsed();
    Int_t versionGainUsed = procesTRD->GetVersionGainUsed();
    Int_t subversionGainUsed = procesTRD->GetSubVersionGainUsed();
    Int_t versionExBUsed = procesTRD->GetVersionExBUsed();
    Int_t subversionExBUsed = procesTRD->GetSubVersionExBUsed();
    printf("version and subversion vdrift %d and %d\n",versionVdriftUsed,subversionVdriftUsed);
    printf("version and subversion gain %d and %d\n",versionGainUsed,subversionGainUsed);
    printf("version and subversion exb %d and %d\n",versionExBUsed,subversionExBUsed);
    procesTRD->Process("CalibObjects.root",runNumber,runNumber,targetStorage);
  }
  else {
    Printf("\n******* NOT Calibrating TRD: detStr = %s, TRD_qf = %d *******", detStr.Data(), (Int_t)TRD_qf);
  }
  
  TF1 *gsf = (TF1 *)gROOT->GetFunction("gaus");
  if (gsf) for (int i=gsf->GetNpar();i--;) gsf->SetParError(i,0); // reset errors from previous fits
  
  //Mean Vertex
  AliMeanVertexPreprocessorOffline * procesMeanVtx=0;
  if (detStr.Contains("ITSSPD") && SPD_qf) {
    Printf("\n******* Calibrating MeanVertex *******");
    procesMeanVtx = new AliMeanVertexPreprocessorOffline;
    procesMeanVtx->ProcessOutput("CalibObjects.root", targetStorage, runNumber);
  }
  else {
    Printf("\n******* NOT Calibrating MeanVertex: detStr = %s, SPD_qf = %d *******", detStr.Data(), (Int_t)SPD_qf);
  }

  //
  // Print calibration status into the stdout
  //
  Int_t trdStatus = (procesTRD) ?  procesTRD->GetStatus():0;
  Int_t tofStatus = (procesTOF) ?  procesTOF->GetStatus():0;
  Int_t t0Status  = (procesT0)  ?  procesT0->GetStatus():0;
  Int_t tpcStatus = (procesTPC) ?  procesTPC->GetStatus():0;
  Int_t meanVtxStatus = (procesMeanVtx) ? procesMeanVtx->GetStatus():0;
  //
  printf("\n");
  printf("******* CPass0 calibration status *******\n");
  printf("TRD calibration status=%d\n",trdStatus);
  printf("TOF calibration status=%d\n",tofStatus);
  printf("TPC calibration status=%d\n",tpcStatus);
  printf("T0  calibration status=%d\n",t0Status);
  printf("MeanVertex  calibration status=%d\n",meanVtxStatus);
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("cpassStat.root","recreate");
  printCalibStat(runNumber, "CalibObjects.root",pcstream);
  delete pcstream;
  return;
}


// function to print statistics used to calibrate the various detectors

void printCalibStat(Int_t run, const char * fname,  TTreeSRedirector * pcstream){

  //
  // Dump the statistical information about all histograms in the calibration files 
  //    into the statistical tree, print on the screen (log files) as well 
  //
  //
  // 1. Default dump for all histograms
  //    Information to dump:
  //    stat =Entries, Mean, MeanError,  RMS, MaxBin
  //    Branch naming convention:
  //    <detName>_<hisName><statName>
  //
  // 2. Detector statistical information  - to be implemented by expert
  //                                      - First version implemented by MI 
  //  
  // 

  TFile *fin = TFile::Open(fname);
  if (!fin) return;
  const Double_t kMaxHis=10000;
  
  TList * keyList = fin->GetListOfKeys();
  Int_t nkeys=keyList->GetEntries();
  Double_t *hisEntries = new Double_t[kMaxHis];
  Double_t *hisMean = new Double_t[kMaxHis];
  Double_t *hisMeanError = new Double_t[kMaxHis];
  Double_t *hisRMS = new Double_t[kMaxHis];
  Double_t *hisMaxBin = new Double_t[kMaxHis];
  Int_t counter=0;
  
  if (pcstream) (*pcstream)<<"calibStatAll"<<"run="<<run;
  for (Int_t ikey=0; ikey<nkeys; ikey++){
    TObject * object = fin->Get(keyList->At(ikey)->GetName());
    if (!object) continue;
    if (object->InheritsFrom("TCollection")==0) continue;
    TSeqCollection *collection  = (TSeqCollection*)object; 
    Int_t nentries= collection->GetEntries();
    for (Int_t ihis=0; ihis<nentries; ihis++){
      TObject * ohis = collection->At(ihis);
      if (!ohis) continue;
      if (ohis->InheritsFrom("TH1")==0) continue;
      TH1* phis = (TH1*)ohis;
      hisEntries[counter]=phis->GetEntries();	
      Int_t idim=1;
      if (ohis->InheritsFrom("TH2")) idim=2;
      if (ohis->InheritsFrom("TH3")) idim=3;        
      hisMean[counter]=phis->GetMean(idim);	
      hisMeanError[counter]=phis->GetMeanError(idim);	
      hisRMS[counter]=phis->GetRMS(idim);	
      hisMaxBin[counter]=phis->GetBinCenter(phis->GetMaximumBin());	
      if (pcstream) (*pcstream)<<"calibStatAll"<<
		      Form("%s_%sEntries=",keyList->At(ikey)->GetName(), phis->GetName())<<hisEntries[counter]<<	
		      Form("%s_%sMean=",keyList->At(ikey)->GetName(), phis->GetName())<<hisMean[counter]<<	
		      Form("%s_%sMeanError=",keyList->At(ikey)->GetName(), phis->GetName())<<hisMeanError[counter]<<	
		      Form("%s_%sRMS=",keyList->At(ikey)->GetName(), phis->GetName())<<hisRMS[counter]<<	
		      Form("%s_%sMaxBin=",keyList->At(ikey)->GetName(), phis->GetName())<<hisMaxBin[counter];	
      //printf("Histo:\t%s_%s\t%f\t%d\n",keyList->At(ikey)->GetName(), phis->GetName(), hisEntries[counter],idim);
      counter++;
    }
    delete object;
  }    
  
  //
  // Expert dump example (MI first iteration):
  //
  // 0.)  TOF dump
  //

  Int_t tofEvents=0;
  Int_t tofTracks=0;
  TList * TOFCalib = (TList*)fin->Get("TOFHistos");      
  if (TOFCalib) {
    TH1 *histoEvents = (TH1*)TOFCalib->FindObject("hHistoVertexTimestamp");
    TH1 *histoTracks = (TH1*)TOFCalib->FindObject("hHistoDeltatTimestamp");
    if (histoEvents && histoTracks){
      tofEvents = TMath::Nint(histoEvents->GetEntries());
      tofTracks = TMath::Nint(histoTracks->GetEntries());
    }
    delete TOFCalib;
  }
  printf("Monalisa TOFevents\t%d\n",tofEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TOFevents="<<tofEvents;
  printf("Monalisa TOFtracks\t%d\n",tofTracks);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TOFtracks="<<tofTracks;

  //
  // 1.)  TPC  dump - usefull events/tracks  for the calibration
  //
  Int_t tpcEvents=0;
  Int_t tpcTracks=0;
  TObject* obj = dynamic_cast<TObject*>(fin->Get("TPCCalib"));
  TObjArray* array = dynamic_cast<TObjArray*>(obj);
  TDirectory* dir = dynamic_cast<TDirectory*>(obj);
  AliTPCcalibTime  * calibTime = NULL;
  if (dir) {
    calibTime = dynamic_cast<AliTPCcalibTime*>(dir->Get("calibTime"));
  }
  else if (array){
    calibTime = (AliTPCcalibTime *)array->FindObject("calibTime");
  }
  if (calibTime) {
      tpcEvents = TMath::Nint(calibTime->GetTPCVertexHisto(0)->GetEntries());
      tpcTracks = TMath::Nint(calibTime->GetResHistoTPCITS(0)->GetEntries());
  }
  printf("Monalisa TPCevents\t%d\n",tpcEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TPCevents="<<tpcEvents;
  printf("Monalisa TPCtracks\t%d\n",tpcTracks);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TPCtracks="<<tpcTracks;

  //
  // 2. TRD dump 
  //
  Int_t trdEvents=0;
  Int_t trdTracks=0;
  TList * TRDCalib = (TList*)fin->Get("TRDCalib");      
  if (TRDCalib) {
    TH1  *histoEvents = (TH1*)TRDCalib->FindObject("NEventsInput_AliTRDCalibTask");
    TH1  *histoTracks = (TH1*)TRDCalib->FindObject("AbsoluteGain_AliTRDCalibTask");
    if (histoEvents && histoTracks){
      trdEvents= TMath::Nint(histoEvents->GetEntries());
      trdTracks= TMath::Nint(histoTracks->GetEntries());
    }
    delete TRDCalib;
  }
  printf("Monalisa TRDevents\t%d\n",trdEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TRDevents="<<trdEvents;
  printf("Monalisa TRDtracks\t%d\n",trdTracks);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TRDtracks="<<trdTracks;

  //
  // 3. T0 dump 
  //
  Int_t T0Events=0;
  TList * T0Calib = (TList*)fin->Get("T0Calib");      
  if (T0Calib) {
    TH1  *histoEvents = (TH1*) T0Calib->FindObject("fTzeroORAplusORC");
    if (histoEvents){
      T0Events= TMath::Nint(histoEvents->GetEntries());
    }
    delete T0Calib;
  }
  printf("Monalisa T0events\t%d\n",T0Events);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"T0events="<<T0Events;

  //
  // 4. Mean vertex -   dump 
    Int_t meanVertexEvents=0;
  TList * meanVertexCalib = (TList*)fin->Get("MeanVertex");      
  if (meanVertexCalib) {
    TH1  *histoEvents = (TH1*) meanVertexCalib->FindObject("hTRKVertexX");
    if (histoEvents){
      meanVertexEvents = TMath::Nint(histoEvents->GetEntries());
    }
    delete meanVertexCalib;
  }
  printf("Monalisa MeanVertexevents\t%d\n",meanVertexEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"MeanVertexevents="<<meanVertexEvents;

  //
  // 5. SDD dump 
  //
  Int_t sddEvents=0;
  Int_t sddTracks=0;
  TList * SDDCalib = (TList*)fin->Get("clistSDDCalib");      
  if (SDDCalib) {
    TH1  *histoEvents = (TH1*) SDDCalib->FindObject("hNEvents");
    if (histoEvents ){
      sddEvents = TMath::Nint(histoEvents->GetBinContent(4));
      sddTracks = TMath::Nint(histoEvents->GetBinContent(5));
    }
    delete SDDCalib;
  }
  printf("Monalisa SDDevents\t%d\n",sddEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"SDDevents="<<sddEvents;
  printf("Monalisa SDDtracks\t%d\n",sddTracks);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"SDDtracks="<<sddTracks;

  //
  if (pcstream) (*pcstream)<<"calibStatAll"<<"\n";
  delete fin;

}

