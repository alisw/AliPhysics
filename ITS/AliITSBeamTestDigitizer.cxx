////////////////////////////////////////////////////
//  Class to manage the                           //
//  ITS beam test conversion from rawdata         //
//  to digits. It executes the digitization for   //
//  SPD, SDD and SSD.                             //
//  Origin:  E. Crescio crescio@to.infn.it        //
//           J. Conrad  Jan.Conrad@cern.ch        //
////////////////////////////////////////////////////
#include "AliHeader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliITSEventHeader.h"
#include "AliITSLoader.h"
#include "AliITSBeamTest.h"
#include "AliITSBeamTestDigSDD.h"
#include "AliITSBeamTestDigSPD.h"
#include "AliITSBeamTestDigSSD.h"
#include "AliITSBeamTestDigitizer.h"
#include "AliRawReaderDate.h"


const TString AliITSBeamTestDigitizer::fgkDefaultDigitsFileName="ITSbt.Digits.root";  

ClassImp(AliITSBeamTestDigitizer)


//_____________________________________________________________
AliITSBeamTestDigitizer::AliITSBeamTestDigitizer():TTask() 
{  
  //
  // Default constructor
  //
  fRunLoader = 0;
  fLoader =0;
  fEvIn=0;
  fEvFin=0;
  fFlagHeader=kTRUE;
  fDATEEvType=7;
  fRunNumber=-1;
  SetFlagInit();
  fBt=0;
  SetBeamTestPeriod();
} 

//_____________________________________________________________
  AliITSBeamTestDigitizer::AliITSBeamTestDigitizer(const Text_t* name, const Text_t* title):TTask(name,title) 
{  
  //
  // Standard constructor 
  //
  Init();
  fEvIn=0;
  fEvFin=0;
  fDATEEvType=7;
  fFlagHeader=kTRUE;
  fRunNumber=-1;
  SetBeamTestPeriod();
 } 

//_____________________________________________________________
  AliITSBeamTestDigitizer::AliITSBeamTestDigitizer(const Text_t* name, const Text_t* title, Int_t run):TTask(name,title) 

{  
  //
  // Constructor 
  //
  Init();
  fEvIn=0;
  fEvFin=0;
  fDATEEvType=7;
  fFlagHeader=kTRUE;
  fRunNumber=run;
  SetBeamTestPeriod();
 } 

//___________________________________________________________
void AliITSBeamTestDigitizer::Init(){

  //
  //Initialization of run loader and its loader 
  //creation of galice.root
  //
  fRunLoader = AliRunLoader::Open("galice.root",
				  AliConfig::GetDefaultEventFolderName(),"recreate");
  
  gAlice->SetRunLoader(fRunLoader);    
  fRunLoader->SetEventFolderName();
  fBt = new AliITSBeamTest("ITS","ITS beam test");
  fBt->SetDefaults();
  gAlice->AddModule(fBt);
  fRunLoader->AddLoader(fBt);
  fLoader = (AliITSLoader*)fRunLoader->GetLoader("ITSLoader");
  fRunLoader->MakeTree("E");  

  fRunLoader->WriteRunLoader("OVERWRITE");
  fRunLoader->WriteAliRun("OVERWRITE");
  fDigitsFileName=fgkDefaultDigitsFileName;
  this->Add(new AliITSBeamTestDigSPD("DigSPD","SPD Digitization")); 
  this->Add(new AliITSBeamTestDigSDD("DigSDD","SDD Digitization")); 
  this->Add(new AliITSBeamTestDigSSD("DigSSD","SSD Digitization"));

  SetFlagInit(kTRUE);
}

//_____________________________________________________________
  AliITSBeamTestDigitizer::AliITSBeamTestDigitizer(const char* filename)
{
  //
  // Constructor for reading (reads galice.root)
  //

  fRunLoader = AliRunLoader::Open(filename);
  if (fRunLoader == 0x0)
    {
      Error("AliITSBeamTestDigitizer","Can not load the session",filename);
      return;
    }
  fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();
  
  if(!gAlice) {
    Error("AliITSBeamTestDigitizer","gAlice not found on file. Aborting.");
    return;
  } 
  
  fLoader = (AliITSLoader*)fRunLoader->GetLoader("ITSLoader");
  
  fBt = (AliITSBeamTest*)gAlice->GetDetector("ITS");
  fBt->SetDefaults();

  fDigitsFileName=fgkDefaultDigitsFileName;

  fEvIn=0;
  fEvFin=0;
  
}

//______________________________________________________________________
AliITSBeamTestDigitizer::AliITSBeamTestDigitizer(const AliITSBeamTestDigitizer &bt):TTask(bt){
  // Copy constructor. 
  //not allowed
  if(this==&bt) return;
  Error("Copy constructor",
	"You are not allowed to make a copy of the AliITSBeamTestDigitizer");
  exit(1);

}
//______________________________________________________________________
AliITSBeamTestDigitizer& AliITSBeamTestDigitizer::operator=(AliITSBeamTestDigitizer &bt){
    // Assignment operator. This is a function which is not allowed to be
    // done to the ITS beam test digitizer. It exits with an error.
    // Inputs:
    if(this==&bt) return *this;
    Error("operator=","You are not allowed to make a copy of the AliITSBeamTestDigitizer");
    exit(1);
    return *this; //fake return
}


//______________________________________________________________
AliITSBeamTestDigitizer::~AliITSBeamTestDigitizer(){

  //Destructor
  if(fBt) delete fBt;
  if(fLoader) delete fLoader;
  if(fRunLoader) delete fRunLoader;
} 


//_____________________________________________________________
void AliITSBeamTestDigitizer::SetNumberOfEventsPerFile(Int_t nev)
{
  //Sets number of events per file

  if(fRunLoader) fRunLoader->SetNumberOfEventsPerFile(nev);
  else Warning("SetNumberOfEventsPerFile","fRunLoader is 0");
}


//____________________________________________________
void AliITSBeamTestDigitizer::ExecDigitization(){

  // Execution of digitisation for SPD,SDD and SSD

  if(!GetFlagInit()){
    Warning("ExecDigitization()","Run Init() please..");
    return;
  }
  fLoader->SetDigitsFileName(fDigitsFileName);
  fLoader->LoadDigits("recreate");
 
  AliRawReaderDate rd(fRawdataFileName,fEvIn);
  AliHeader* header = fRunLoader->GetHeader();
  
  Int_t iev=fEvIn-1;

  
  AliITSBeamTestDigSDD* digSDD = (AliITSBeamTestDigSDD*)fTasks->FindObject("DigSDD");
  AliITSBeamTestDigSPD* digSPD = (AliITSBeamTestDigSPD*)fTasks->FindObject("DigSPD");
  AliITSBeamTestDigSSD* digSSD = (AliITSBeamTestDigSSD*)fTasks->FindObject("DigSSD");


  do{
    iev++;
    if(fEvFin!=0){
      if(iev>fEvFin) break;
    } 
    AliITSEventHeader* itsh = new AliITSEventHeader("ITSHeader");
    fRunLoader->SetEventNumber(iev);
   
    rd.RequireHeader(fFlagHeader);
    rd.SelectEvents(fDATEEvType);
 
    digSDD->SetRawReaderDate(&rd);
    digSPD->SetRawReaderDate(&rd);
    digSSD->SetRawReaderDate(&rd);
    
    if(fLoader->TreeD() == 0x0) fLoader->MakeTree("D");

    TTree* treeD = (TTree*)fLoader->TreeD();
   
    // Make branches outside the dig-classes

    TClonesArray* digitsSPD = new TClonesArray("AliITSdigitSPD",1000);
    treeD->Branch("ITSDigitSPD",&digitsSPD);
 
    TClonesArray* digitsSDD = new TClonesArray("AliITSdigitSDD",1000);
    treeD->Branch("ITSDigitSDD",&digitsSDD);
   
    TClonesArray* digitsSSD = new TClonesArray("AliITSdigitSSD",1000);
    treeD->Branch("ITSDigitSSD",&digitsSSD);


    digSSD->SetTree(treeD);
    digSDD->SetTree(treeD);
    digSPD->SetTree(treeD);

    digSSD->SetBeamTest(fBt);
    digSDD->SetBeamTest(fBt);
    digSPD->SetBeamTest(fBt);

    digSSD->SetITSEventHeader(itsh);
    digSDD->SetITSEventHeader(itsh);
    digSPD->SetITSEventHeader(itsh);

    digSDD->SetBtPeriod(GetBeamTestPeriod());
    digSDD->SetThreshold(16);

    ExecuteTask(0);  

    header->SetEventNrInRun(iev);
    header->SetEvent(iev);
    header->SetRun(fRunNumber);
    fRunLoader->GetHeader()->AddDetectorEventHeader(itsh);
    fRunLoader->TreeE()->Fill();
    header->Reset(fRunNumber,iev);
    
    delete digitsSPD;
    delete digitsSDD;
    delete digitsSSD;

   }while(rd.NextEvent());

  fRunLoader->WriteHeader("OVERWRITE");
  fRunLoader->WriteRunLoader("OVERWRITE");
  fLoader->UnloadDigits();
  fLoader->UnloadRawClusters();
  fRunLoader->UnloadHeader();
  
}



//_______________________________________________
void AliITSBeamTestDigitizer:: SetActive(const TString& subdet,Bool_t value){

  //Sets active sub-tasks (detectors)
  
  Bool_t sdd = subdet.Contains("SDD");
  Bool_t spd = subdet.Contains("SPD");
  Bool_t ssd = subdet.Contains("SSD");

  if(sdd){
  AliITSBeamTestDigSDD* digSDD = (AliITSBeamTestDigSDD*)fTasks->FindObject("DigSDD");
  digSDD->SetActive(value);
  }
  
 if(spd){
  AliITSBeamTestDigSPD* digSPD = (AliITSBeamTestDigSPD*)fTasks->FindObject("DigSPD");
  digSPD->SetActive(value);
  
  }
 
  if(ssd){
  AliITSBeamTestDigSSD* digSSD = (AliITSBeamTestDigSSD*)fTasks->FindObject("DigSSD");
  digSSD->SetActive(value);
  
  }



}

