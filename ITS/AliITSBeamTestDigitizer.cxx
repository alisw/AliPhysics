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
#include "AliITSvBeamTestITS04.h"
#include "AliITSvSDD03.h"
#include "AliITSBeamTestDigSDD.h"
#include "AliITSBeamTestDigSPD.h"
#include "AliITSBeamTestDigSSD.h"
#include "AliITSBeamTestDigitizer.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include <TVirtualMC.h>
//#include <TGeant3.h>

const TString AliITSBeamTestDigitizer::fgkDefaultDigitsFileName="ITS.Digits.root";  

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
  SetOptDate();
  fBt=0;
  fPeriod=kNov04;
} 

//_____________________________________________________________
  AliITSBeamTestDigitizer::AliITSBeamTestDigitizer(const Text_t* name, const Text_t* title, Char_t* opt):TTask(name,title) 
{  
  //
  // Standard constructor 
  //
  fRunLoader=0;
  fLoader=0;
  fBt=0;
  fEvIn=0;
  fEvFin=0;
  fDATEEvType=7;
  fFlagHeader=kTRUE;
  fRunNumber=-1;
  SetOptDate();

  TString choice(opt);
  Bool_t aug04 = choice.Contains("Aug04");
  Bool_t nov04 = choice.Contains("Nov04");
  if(aug04) fPeriod=kAug04;
  if(nov04) fPeriod=kNov04;
  Init();
 } 

//_____________________________________________________________
  AliITSBeamTestDigitizer::AliITSBeamTestDigitizer(const Text_t* name, const Text_t* title, Int_t run, Char_t* opt):TTask(name,title) 

{  
  //
  // Constructor 
  //
  fRunLoader=0;
  fLoader=0;
  fBt=0;
  fEvIn=0;
  fEvFin=0;
  fDATEEvType=7;
  fFlagHeader=kTRUE;
  fRunNumber=run;
  SetOptDate();
  TString choice(opt);
  Bool_t aug04 = choice.Contains("Aug04");
  Bool_t nov04 = choice.Contains("Nov04");
  if(aug04) fPeriod=kAug04;
  if(nov04) fPeriod=kNov04;

  Init();
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
  if(GetBeamTestPeriod()==kNov04){
    if(gGeoManager) delete gGeoManager;
    gGeoManager = new TGeoManager("ITSGeometry","ITS Simulation Geometry Manager");
    TGeoMaterial *vacmat = new TGeoMaterial("Vacuum",0,0,0);
    TGeoMedium   *vacmed = new TGeoMedium("Vacuum_med",1,vacmat);
    TGeoVolume *aLICE = gGeoManager->MakeBox("ALICE",vacmed,100.,100.,200.);
    gGeoManager->SetTopVolume(aLICE);  
    fBt = new AliITSvBeamTestITS04("ITS","ITS beam test");
    fBt->CreateGeometry();
    fBt->Init();

  }
  if(GetBeamTestPeriod()==kAug04){
    fBt = new AliITSvSDD03("ITS",2004);
    gSystem->Load("libgeant321");
    //    new TGeant3("C++ Interface to Geant3");
    if(strcmp(gMC->GetName(),"TGeant3")) {
	Fatal("Init","TGeant3 should be instantiated in advance");
	return;
    } 
    fBt->CreateMaterials();
    fBt->CreateGeometry();
    fBt->Init();
  }
  gAlice->AddModule(fBt);
  fBt->SetDefaults();   
  fRunLoader->AddLoader(fBt);
  fLoader = (AliITSLoader*)fRunLoader->GetLoader("ITSLoader");
  fRunLoader->MakeTree("E");  
  fRunLoader->WriteHeader("OVERWRITE");
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
  
  //fBt = (AliITSBeamTest*)gAlice->GetDetector("ITS");
  fBt->SetDefaults();

  fDigitsFileName=fgkDefaultDigitsFileName;

  fEvIn=0;
  fEvFin=0;
  SetOptDate();
  
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
AliITSBeamTestDigitizer& AliITSBeamTestDigitizer::operator=(const AliITSBeamTestDigitizer &source){
    // Assignment operator. This is a function which is not allowed to be
    // done to the ITS beam test digitizer. It exits with an error.
    // Inputs:
    if(this==&source) return *this;
    Error("operator=","You are not allowed to make a copy of the AliITSBeamTestDigitizer");
    exit(1);
    return *this; //fake return
}


//______________________________________________________________
AliITSBeamTestDigitizer::~AliITSBeamTestDigitizer(){

  //Destructor
  if(fBt) delete fBt;
  if(fLoader) delete fLoader;
  if(fHeader) delete fHeader;
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
 
  AliRawReader* rd;

  if(GetOptDate()) rd = new AliRawReaderDate(fRawdataFileName,fEvIn);
  else rd = new AliRawReaderRoot(fRawdataFileName,fEvIn);

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
   
    rd->RequireHeader(fFlagHeader);
    rd->SelectEvents(fDATEEvType);
 
    digSDD->SetRawReader(rd);
    digSPD->SetRawReader(rd);
    digSSD->SetRawReader(rd);
    
    if(fLoader->TreeD() == 0x0) fLoader->MakeTree("D");

    TTree* treeD = (TTree*)fLoader->TreeD();
   
    // Make branches outside the dig-classes

    TClonesArray* digitsSPD = new TClonesArray("AliITSdigitSPD",1000);
    treeD->Branch("ITSDigitsSPD",&digitsSPD);
 
    TClonesArray* digitsSDD = new TClonesArray("AliITSdigitSDD",1000);
    treeD->Branch("ITSDigitsSDD",&digitsSDD);
   
    TClonesArray* digitsSSD = new TClonesArray("AliITSdigitSSD",1000);
    treeD->Branch("ITSDigitsSSD",&digitsSSD);


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
    if(GetBeamTestPeriod()==1)digSDD->SetThreshold(16);
    else digSDD->SetThreshold(0);
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

   }while(rd->NextEvent());

  
  fRunLoader->WriteHeader("OVERWRITE");
  fRunLoader->WriteRunLoader("OVERWRITE");

  delete rd;
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

