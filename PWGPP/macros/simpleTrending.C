void processContainer(TObject* object, TTreeSRedirector* debugStreamer, TString name);
void loadLibraries();
void simpleTrending(TString inputFileName, Int_t run, TString filterExpr=".*", TString trendingFileName="trending.root", TString treeName="trending", TString fileOpenMode="update" );

TString treeName;

void simpleTrending(TString inputFileName, Int_t run, TString filterExpr, TString trendingFileName, TString debugTreeName, TString fileOpenMode )
{

  // Dump the statistical information about all histograms in the file
  //    into a tree
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
  //                                      - updated for QA by MK
  //  
  // 

  loadLibraries();

  treeName=debugTreeName;

  TFile *inputFile = TFile::Open(inputFileName.Data());
  if (!inputFile) return;
  
  TList * keyList = inputFile->GetListOfKeys();
  Int_t nkeys=keyList->GetEntries();
 
  TRegexp filterRegexp=filterExpr.Data();

  //check if we have a matching container, only then create the output file
  TList* keyList=inputFile->GetListOfKeys();
  Int_t nkeys=keyList->GetEntries();
  Bool_t containerExists=kFALSE;
  for (Int_t i=0; i<nkeys; i++)
  {
    TObject* object=keyList->At(i);
    if (!object) continue;
    TString name=object->GetName();
    if (name.Contains(filterRegexp))
    {
      containerExists=kTRUE;
      break;
    }
  }
  if (!containerExists) 
  {
    printf("container %s does not exist in %s\n",filterExpr.Data(),inputFileName.Data());
    return;
  }

  TTreeSRedirector *pcstream = new TTreeSRedirector(trendingFileName,fileOpenMode);
  (*pcstream)<<treeName.Data()<<"run="<<run;
  
  //main loop over the top level objects, filtering is done here
  for (Int_t i=0; i<nkeys; i++)
  {
    TObject * object = inputFile->Get(keyList->At(i)->GetName());
    if (!object) continue;
    TString name=object->GetName();
    if (!name.Contains(filterRegexp)) continue;
    processContainer(object,pcstream,name);
  }
  
  //
  // Expert dump example (MI first iteration):
  //
  // 0.)  TOF dump
  //

  Int_t tofEvents=0;
  Int_t tofTracks=0;
  TList * TOFCalib = (TList*)inputFile->Get("TOFHistos");      
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
  if (pcstream) (*pcstream)<<treeName.Data()<<"TOFevents="<<tofEvents;
  printf("Monalisa TOFtracks\t%d\n",tofTracks);
  if (pcstream) (*pcstream)<<treeName.Data()<<"TOFtracks="<<tofTracks;

  //
  // 1.)  TPC  dump - usefull events/tracks  for the calibration
  //
  Int_t tpcEvents=0;
  Int_t tpcTracks=0;
  TObject* obj = dynamic_cast<TObject*>(inputFile->Get("TPCCalib"));
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
  if (pcstream) (*pcstream)<<treeName.Data()<<"TPCevents="<<tpcEvents;
  printf("Monalisa TPCtracks\t%d\n",tpcTracks);
  if (pcstream) (*pcstream)<<treeName.Data()<<"TPCtracks="<<tpcTracks;

  //
  // 2. TRD dump 
  //
  Int_t trdEvents=0;
  Int_t trdTracks=0;
  TList * TRDCalib = (TList*)inputFile->Get("TRDCalib");      
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
  if (pcstream) (*pcstream)<<treeName.Data()<<"TRDevents="<<trdEvents;
  printf("Monalisa TRDtracks\t%d\n",trdTracks);
  if (pcstream) (*pcstream)<<treeName.Data()<<"TRDtracks="<<trdTracks;

  //
  // 3. T0 dump 
  //
  Int_t T0Events=0;
  TList * T0Calib = (TList*)inputFile->Get("T0Calib");      
  if (T0Calib) {
    TH1  *histoEvents = (TH1*) T0Calib->FindObject("fTzeroORAplusORC");
    if (histoEvents){
      T0Events= TMath::Nint(histoEvents->GetEntries());
    }
    delete T0Calib;
  }
  printf("Monalisa T0events\t%d\n",T0Events);
  if (pcstream) (*pcstream)<<treeName.Data()<<"T0events="<<T0Events;

  //
  // 4. Mean vertex -   dump 
  // Not present in CPass1
  /*
    Int_t meanVertexEvents=0;
  TList * meanVertexCalib = (TList*)inputFile->Get("MeanVertex");      
  if (meanVertexCalib) {
    TH1  *histoEvents = (TH1*) meanVertexCalib->FindObject("hTRKVertexX");
    if (histoEvents){
      meanVertexEvents = TMath::Nint(histoEvents->GetEntries());
    }
    delete meanVertexCalib;
  }
  printf("Monalisa MeanVertexevents\t%d\n",meanVertexEvents);
  if (pcstream) (*pcstream)<<treeName.Data()<<"MeanVertexevents="<<meanVertexEvents;
  */

  //
  // 5. SDD dump 
  //
  Int_t sddEvents=0;
  Int_t sddTracks=0;
  TList * SDDCalib = (TList*)inputFile->Get("clistSDDCalib");      
  if (SDDCalib) {
    TH1  *histoEvents = (TH1*) SDDCalib->FindObject("hNEvents");
    if (histoEvents ){
      sddEvents = TMath::Nint(histoEvents->GetBinContent(4));
      sddTracks = TMath::Nint(histoEvents->GetBinContent(5));
    }
    delete SDDCalib;
  }
  printf("Monalisa SDDevents\t%d\n",sddEvents);
  if (pcstream) (*pcstream)<<treeName.Data()<<"SDDevents="<<sddEvents;
  printf("Monalisa SDDtracks\t%d\n",sddTracks);
  if (pcstream) (*pcstream)<<treeName.Data()<<"SDDtracks="<<sddTracks;

  //
  if (pcstream) (*pcstream)<<treeName.Data()<<"\n";
  delete pcstream;
  delete inputFile;

}

void processContainer(TObject* inputObject, TTreeSRedirector* pcstream, TString parentname)
{
  //recursively process the contents of an object:
  //might be a TDirectory
  //or a TCollection
  //pus information about the contained histograms in the debugStreamer

  TDirectory* inputDir=NULL;
  TSeqCollection* inputCollection=NULL;
  TH1* inputHistogram=NULL;
  
  TString inputObjectName=inputObject->GetName();

  //TDirectory* inputDir=dynamic_cast<TDirectory*>(inputObject);
  //TSeqCollection* inputCollection=dynamic_cast<TSeqCollection*>(inputObject);
  //TH1* inputHistogram=dynamic_cast<TH1*>(inputObject);
  if (inputObject->InheritsFrom("TDirectory")) inputDir=dynamic_cast<TDirectory*>(inputObject);
  if (inputObject->InheritsFrom("TSeqCollection")) inputCollection=dynamic_cast<TSeqCollection*>(inputObject);
  if (inputObject->InheritsFrom("TH1")) inputHistogram=dynamic_cast<TH1*>(inputObject);

  if (inputCollection)
  {
    printf("processing collection: %s\n",inputCollection->GetName());
    Int_t nentries= inputCollection->GetEntries();
    for (Int_t i=0; i<nentries; i++)
    {
      TObject * object = inputCollection->At(i);
      if (!object) continue;
      TString name=parentname+"/"+object->GetName();
      processContainer(object,pcstream,name);
    }
  } 
  else if (inputDir)
  {
    printf("processing directory: %s\n",inputDir->GetName());
    TList* keyList=inputDir->GetListOfKeys();
    Int_t nkeys=keyList->GetEntries();
    for (Int_t i=0; i<nkeys; i++)
    {
      TObject * object = inputDir->Get(keyList->At(i)->GetName());
      if (!object) continue;
      TString name=parentname+"/"+object->GetName();
      processContainer(object,pcstream,name);
    }

  }
  else if (inputHistogram)
  {
    Double_t hisEntries;
    Double_t hisMean;
    Double_t hisMeanError;
    Double_t hisRMS;
    Double_t hisMaxBin;
    if (inputHistogram->InheritsFrom("TH1")==0) continue;
    TH1* phis = (TH1*)inputHistogram;
    hisEntries=phis->GetEntries();	
    Int_t idim=1;
    if (inputHistogram->InheritsFrom("TH2")) idim=2;
    if (inputHistogram->InheritsFrom("TH3")) idim=3;        
    hisMean=phis->GetMean(idim);	
    hisMeanError=phis->GetMeanError(idim);	
    hisRMS=phis->GetRMS(idim);	
    hisMaxBin=phis->GetBinCenter(phis->GetMaximumBin());
    TString name=parentname;
    printf("histogram: %s\n",name.Data());
    if (pcstream) (*pcstream)<<treeName.Data()<<
        Form("%s_Entries=",name.Data())<<hisEntries<<	
        Form("%s_Mean=",name.Data())<<hisMean<<	
        Form("%s_MeanError=",name.Data())<<hisMeanError<<	
        Form("%s_RMS=",name.Data())<<hisRMS<<	
        Form("%s_MaxBin=",name)<<hisMaxBin;	
  }

}

void loadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWGPP -     I$ALICE_ROOT/PWGPP/TRD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libESDfilter");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libAliHLTTrigger.so");
  gSystem->Load("libPWGTools");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libPWGGACaloTrackCorrelations");
  gSystem->Load("libPWGGACaloTasks");
  gSystem->Load("libPWGGAPHOSTasks");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");
  gSystem->Load("libPWGGAEMCALTasks");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libPWGPPMUONlite");
  gSystem->Load("libPWGmuondep");
  gSystem->Load("libPWGLFforward2");
}

