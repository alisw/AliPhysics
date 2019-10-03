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

  if (inputCollection){
    printf("processing collection: %s\n",inputCollection->GetName());
    Int_t nentries= inputCollection->GetEntries();
    for (Int_t i=0; i<nentries; i++)
      {
	TObject * object = inputCollection->At(i);
	if (!object) continue;
	TString name=parentname+"."+object->GetName();
	processContainer(object,pcstream,name);
      }
  } 
  else if (inputDir){
    printf("processing directory: %s\n",inputDir->GetName());
    TList* keyList=inputDir->GetListOfKeys();
    Int_t nkeys=keyList->GetEntries();
    for (Int_t i=0; i<nkeys; i++){
      TObject * object = inputDir->Get(keyList->At(i)->GetName());
      if (!object) continue;
      TString name=parentname+"."+object->GetName();
      processContainer(object,pcstream,name);
    }
  }
  else if (inputHistogram){
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
      Form("%s.Entries=",name.Data())<<hisEntries<<	
      Form("%s.Mean=",name.Data())<<hisMean<<	
      Form("%s.MeanError=",name.Data())<<hisMeanError<<	
      Form("%s.RMS=",name.Data())<<hisRMS<<	
      Form("%s.MaxBin=",name)<<hisMaxBin;	
  }
}

void loadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/PWGPP -     I$ALICE_PHYSICS/PWGPP/TRD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libESDfilter");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  gSystem->Load("libAliHLTTrigger");
  gSystem->Load("libPWGTools");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libPWGGACaloTrackCorrelations");
  gSystem->Load("libPWGGAPHOSTasks");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCALbase");
  gSystem->Load("libPWGGAEMCALTasks");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libPWGPPMUONlite");
  gSystem->Load("libPWGmuondep");
  gSystem->Load("libPWGLFforward2");
}

