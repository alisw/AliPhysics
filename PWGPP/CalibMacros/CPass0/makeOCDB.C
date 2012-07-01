/*
  macro to extract the OCDB entries

  input: CalibObjects.root
  ouput: TimeGain and TimeVdrift calibration objects for TPC and TRD

  Example:
  .L $ALICE_ROOT/PWGPP/CalibMacros/CPass0/makeOCDB.C
  makeOCDB("105160");

*/

void printCalibStat(Int_t run, const char * fname,  TTreeSRedirector * pcstream);


void makeOCDB(Int_t runNumber, TString  ocdbStorage="", TString defaultOCDBstorage="raw://")
{
  //
  // extract TPC OCDB entries
  //
  gROOT->Macro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/LoadLibraries.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/ConfigCalibTrain.C");

  // switch off log info
  AliLog::SetClassDebugLevel("AliESDEvent",0);

  // config GRP
  printf("runNumber from runCalibTrain = %d\n",runNumber);
  ConfigCalibTrain(runNumber, defaultOCDBstorage.Data());
  //
  //
  // check the presence of the detectors
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject()); 
  if (!grpData) {printf("Failed to get GRP data for run",runNumber); return;}
  Int_t activeDetectors = grpData->GetDetectorMask(); 
  TString detStr = AliDAQ::ListOfTriggeredDetectors(activeDetectors);
  printf("Detectors in the data:\n%s\n",detStr.Data());
  printf("Monalisa: Detectors in the data:\t%s\n",detStr.Data());


  // Steering Tasks - set output storage
  // DefaultStorage set already before - in ConfigCalibTrain.C
//ocdbStorage+="?se=ALICE::CERN::SE";

  AliCDBManager::Instance()->SetSpecificStorage("*/*/*",ocdbStorage.Data());
  if (gSystem->AccessPathName("TPC", kFileExists)==0) {  
    AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Correction","local://");
  }
  // set OCDB storage
  if (ocdbStorage.Length()==0) ocdbStorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";

  // TPC part
  TFile fcalib("CalibObjects.root");
  AliTPCPreprocessorOffline *procesTPC = 0;
  if  (detStr.Contains("TPC")){ 
    //if  (0){ 
    procesTPC = new AliTPCPreprocessorOffline; 
    // switch on parameter validation
    procesTPC->SetTimeGainRange(0.5,3.0);
    procesTPC->SwitchOnValidation();
    
    // Make timegain calibration
    //proces.CalibTimeGain("CalibObjects.root", runNumber,AliCDBRunRange::Infinity(),ocdbStorage);
    procesTPC->CalibTimeGain("CalibObjects.root", runNumber,runNumber,ocdbStorage);
    
    // Make vdrift calibration
    //proces.CalibTimeVdrift("CalibObjects.root",runNumber,AliCDBRunRange::Infinity(),ocdbStorage);
    procesTPC->CalibTimeVdrift("CalibObjects.root",runNumber,runNumber,ocdbStorage);
  }
  //
  // TOF part
  //
  AliTOFAnalysisTaskCalibPass0 *procesTOF=0;
  if ( detStr.Contains("TOF")){    
    procesTOF = new AliTOFAnalysisTaskCalibPass0;
    Printf("Calibrating TOF");
    procesTOF->ProcessOutput("CalibObjects.root", ocdbStorage);
  }
  //
  //  
  // T0 part
  AliT0PreprocessorOffline *procesT0= 0;
  if ( detStr.Contains("T0")) {
    // Make  calibration of channels offset
    // procesT0.setDArun(179000);
    procesT0= new  AliT0PreprocessorOffline;
    procesT0->Process("CalibObjects.root",runNumber, runNumber, ocdbStorage);
  }
  //
  //TRD part
  //
  AliTRDPreprocessorOffline *procesTRD = 0;
  if ( detStr.Contains("TRD")){
    //if  (0){ 
    procesTRD = new  AliTRDPreprocessorOffline;
    procesTRD->SetLinearFitForVdrift(kTRUE);
    procesTRD->SetMinStatsVdriftT0PH(600*10);
    procesTRD->SetMinStatsVdriftLinear(50);
    procesTRD->SetMinStatsGain(600);
    procesTRD->SetLimitValidateNoData(60);
    procesTRD->SetLimitValidateBadCalib(60);
    procesTRD->SetAlternativeDriftVelocityFit(kTRUE);
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
    procesTRD->Process("CalibObjects.root",runNumber,runNumber,ocdbStorage);
  }
  
  //Mean Vertex
  AliMeanVertexPreprocessorOffline * procesMeanVtx=0;
  if ( detStr.Contains("ITSSPD")) {
    procesMeanVtx =  new AliMeanVertexPreprocessorOffline;
    procesMeanVtx->ProcessOutput("CalibObjects.root", ocdbStorage, runNumber);
  }

  //
  // Print calibration status into the stdout
  //
  Int_t trdStatus = (procesTRD) ?  procesTRD->GetStatus():0;
  Int_t tofStatus = (procesTOF) ?  procesTOF->GetStatus():0;
  Int_t t0Status  = (procesT0)  ?  procesT0->GetStatus():0;
  Int_t tpcStatus = (procesTPC) ? ((procesTPC->ValidateTimeDrift() || procesTPC->ValidateTimeGain())==kFALSE):0;
  //
  printf("\n\n\n\n");
  printf("CPass0 calibration status\n");
  printf("TRD calibration status=%d\n",trdStatus);
  printf("TOF calibration status=%d\n",tofStatus);
  printf("TPC calibration status=%d\n",tpcStatus);
  printf("T0  calibration status=%d\n",t0Status);
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("cpassStat.root","recreate");
  printCalibStat(runNumber, "CalibObjects.root",pcstream);
  delete pcstream;
  return;
}







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
  //
  TList * keyList = fin->GetListOfKeys();
  Int_t nkeys=keyList->GetEntries();
  Double_t *hisEntries = new Double_t[kMaxHis];
  Double_t *hisMean = new Double_t[kMaxHis];
  Double_t *hisMeanError = new Double_t[kMaxHis];
  Double_t *hisRMS = new Double_t[kMaxHis];
  Double_t *hisMaxBin = new Double_t[kMaxHis];
  Int_t counter=0;
  
  if (pcstream) (*pcstream)<<"calibStatAll"<<"run="<<run;
  {for (Int_t ikey=0; ikey<nkeys; ikey++){
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
  }
  //
  //
  // Expert dump example (MI first itteration):
  //
  // 0.)  TOF dump
  //
  Int_t tofEvents=0;
  Int_t tofTracks=0;
  TList * TOFCalib = (TList*)fin->Get("TOFHistos");      
  {if (TOFCalib) {
      TH1  *histoEvents = (TH1*) TOFCalib->FindObject("hHistoVertexTimestamp");
      TH1  *histoTracks = (TH1*)TOFCalib->FindObject("hHistoDeltatTimestamp");
      if (histoEvents && histoTracks){
	tofEvents= histoEvents->GetEntries();
	tofTracks= histoTracks->GetEntries();
      }
      delete TOFCalib;
    }}
  printf("Monalisa TOFevents\t%f\n",tofEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TOFevents="<<tofEvents;
  printf("Monalisa TOFtracks\t%f\n",tofTracks);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TOFtracks="<<tofTracks;
  //
  // 1.)  TPC  dump - usefull events/tracks  for the calibration
  //
  Int_t tpcEvents=0;
  Int_t tpcTracks=0;
  TList * TPCCalib = (TList*)fin->Get("TPCCalib");      
  {if (TPCCalib) {
    AliTPCcalibTime  * calibTime = (AliTPCcalibTime  *)TPCCalib->FindObject("calibTime");
    if (calibTime){
      tpcEvents = calibTime->GetTPCVertexHisto(0)->GetEntries();
      tpcTracks = calibTime->GetResHistoTPCITS(0)->GetEntries();
    }
    }}
  printf("Monalisa TPCevents\t%f\n",tpcEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TPCevents="<<tpcEvents;
  printf("Monalisa TPCtracks\t%f\n",tpcTracks);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TPCtracks="<<tpcTracks;
  //
  // 2. TRD dump 
  //
  Int_t trdEvents=0;
  Int_t trdTracks=0;
  TList * TRDCalib = (TList*)fin->Get("TRDCalib");      
  {if (TRDCalib) {
      TH1  *histoEvents = (TH1*) TRDCalib->FindObject("NEventsInput_AliTRDCalibTask");
      TH1  *histoTracks = (TH1*)TRDCalib->FindObject("AbsoluteGain_AliTRDCalibTask");
      if (histoEvents && histoTracks){
	trdEvents= histoEvents->GetEntries();
	trdTracks= histoTracks->GetEntries();
      }
      delete TRDCalib;
    }}
  printf("Monalisa TRDevents\t%f\n",trdEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TRDevents="<<trdEvents;
  printf("Monalisa TRDtracks\t%f\n",trdTracks);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"TRDtracks="<<trdTracks;
  //
  // 3. T0 dump 
  //
  Int_t T0Events=0;
  TList * T0Calib = (TList*)fin->Get("T0Calib");      
  {if (T0Calib) {
      TH1  *histoEvents = (TH1*) T0Calib->FindObject("fTzeroORAplusORC");
      if (histoEvents && histoTracks){
	T0Events= histoEvents->GetEntries();
      }
      delete T0Calib;
    }}
  printf("Monalisa T0events\t%f\n",T0Events);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"T0events="<<T0Events;
  //
  // 3. Mean vertex -   dump 
  //
  Int_t MeanVertexEvents=0;
  TList * MeanVertexCalib = (TList*)fin->Get("MeanVertex");      
  {if (MeanVertexCalib) {
      TH1  *histoEvents = (TH1*) MeanVertexCalib->FindObject("hSPDVertexX");
      if (histoEvents && histoTracks){
	MeanVertexEvents= histoEvents->GetEntries();
      }
      delete MeanVertexCalib;
    }}
  printf("Monalisa MeanVertexevents\t%f\n",MeanVertexEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"MeanVertexevents="<<MeanVertexEvents;
  //
  // 4. SDD dump 
  //
  Int_t sddEvents=0;
  Int_t sddTracks=0;
  TList * SDDCalib = (TList*)fin->Get("clistSDDCalib");      
  {if (SDDCalib) {
      TH1  *histoEvents = (TH1*) SDDCalib->FindObject("hNEvents");
      if (histoEvents ){
	sddEvents= histoEvents->GetBinContent(4);
	sddTracks= histoEvents->GetBinContent(5);
      }
      delete SDDCalib;
    }}
  printf("Monalisa SDDevents\t%f\n",sddEvents);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"SDDevents="<<sddEvents;
  printf("Monalisa SDDtracks\t%f\n",sddTracks);
  if (pcstream) (*pcstream)<<"calibStatAll"<<"SDDtracks="<<sddTracks;
  //
  //
  if (pcstream) (*pcstream)<<"calibStatAll"<<"\n";
  delete fin;

}

